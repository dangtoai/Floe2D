from math import acos, sin
from scipy.sparse import coo_matrix
from scipy.integrate import solve_ivp
from scipy.spatial import Delaunay
import time
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm
# from graph import *

# import matplotlib.animation as animation


# 4 classes Node-> Spring-> Floe-> Percussion/Percussion_wall to simulate the collision between a floe and a wall

N_T = 1600  # time's discretization


class Node:
    """A class representing one node of an ice floe"""

    def __init__(self, position, velocity=np.array([0, 0]), id_number=None):
        self.x, self.y = position
        # self.x0, self.y0 = position  # Initial position needed for plots
        self.vx, self.vy = velocity
        self.id = id_number
        self.parent_ice_floe = None  # Ice floe to which this node belongs

    def position(self):
        """position of the node"""
        return np.array([self.x, self.y])

    def velocity(self):
        """ velocity of the node """
        return np.array([self.vx, self.vy])

    def get_details(self):
        """ information of the node """
        return self.position(), self.velocity()

    def change_velocity(self, new_velocity: np.array):
        self.velocity = new_velocity


class Spring:
    """
    A class representing one spring of an ice floe
    A spring containts 2 nodes

    """

    def __init__(self, node1: Node, node2: Node, id_number):
        self.node1 = node1
        self.node2 = node2
        self.L0 = norm(node1.position() - node2.position())
        self.id = id_number
        self.parent_ice_floe = None  # Ice floe to which this spring belongs

    def get_nodes(self):
        """ Positions of 2 nodes in this spring """
        return np.array([self.node1.position(), self.node2.position()])

    def get_details(self):
        """ information of the node: positions and length """
        return [self.node1.position(), self.node2.position(), self.L0]


class Floe:
    """
    A class representing an ice floe
    """

    def __init__(self, nodes=None, springs=None, mass=28274333.88230814, stiffness=372350.3982346764, torsion_stiff=32193.692227271033, viscosity=0, tenacity=0,
                 id_number=None, impact_node=False):
        if nodes:
            self.nodes = nodes
            self.n = len(nodes)
        else:
            print("Ice floe "+str(id_number)+" created but contains no nodes")

        if springs:
            self.springs = springs
        else:
            print("Ice floe "+str(id_number) +
                  " created but contains no springs")

        self.E = 8.95*10**9  # Young modulus
        self.nu = 0.295  # Poisson ratio

        # Lame parameters
        self.lamb = self.E*self.nu / ((1+self.nu)*(1-2*self.nu))
        self.g = self.E/(2*(1+self.nu))

        self.m = mass   # mass network
        self.mass_node = self.m/self.n  # mass of each node
        self.k = stiffness
        self.G = torsion_stiff
        self.mu = viscosity
        # self.L = tenacity
        self.id = id_number
        if springs is None:
            possible = []
            for triangle in self.simplices():
                for index1 in triangle:
                    for index2 in triangle:
                        if index1 != index2:
                            possible.append(
                                (min(index1, index2), max(index1, index2)))
            self.springs = set(possible)

        # self.impact_mass = 200  # mass of collided object.

        self.impact_mass = 1e5
        # necessary when collision
        self.mass_nodes = np.full(self.n, self.mass_node)

        if impact_node == True:
            self.mass_nodes[self.n - 1] += self.impact_mass
            # print(self.mass_nodes)

    def degree(self, i):
        count = 0
        for s in self.springs:
            if i in s:
                count += 1
        return count

    def deg(self):
        return [self.degree(i) for i in range(self.n)]

    def maxdegree(self):
        L = [self.degree(i) for i in range(self.n)]
        return max(L), L.index(max(L))

    def meandegree(self):
        L = [self.degree(i) for i in range(self.n)]
        return np.mean(L)

    # def volume(self):
    #     BV = self.border_nodes_index()
    #     return 0.

    def generate_springs(self):
        """ generates springs if the set of springs is already constructed """
        l = []
        for s in self.springs:
            l.append(Spring(self.nodes[s[0]], self.nodes[s[1]], None))
        return l

    def center(self):
        """ center of gravitation of the floe """
        center = sum(self.get_nodes())/self.n
        return np.array([center[0], center[1]])

    def border_edges_index(self):
        All_vertices = []
        for triangles in self.simplices():
            for i in range(len(triangles)):
                for j in range(i+1, len(triangles)):
                    All_vertices.append(
                        (min(triangles[i], triangles[j]), max(triangles[i], triangles[j])))

        S = set(All_vertices)
        Border = ([e for e in S if All_vertices.count(e) == 1])
        return Border

    def border_nodes_index(self):
        BV = self.border_edges_index()

        def find_cycle(edges):
            # Create a dictionary to represent the graph
            graph = {}
            for edge in edges:
                u, v = edge
                if u not in graph:
                    graph[u] = []
                if v not in graph:
                    graph[v] = []
                graph[u].append(v)
                graph[v].append(u)

            def dfs(node, visited, path):
                visited[node] = True
                path.append(node)
                for neighbor in graph[node]:
                    if not visited[neighbor]:
                        dfs(neighbor, visited, path)

            start_node = edges[0][0]
            visited = {node: False for node in graph}
            cycle = []

            dfs(start_node, visited, cycle)

            # To close the cycle, add the start_node at the end
            cycle.append(start_node)

            return cycle

        return find_cycle(BV)
    
    def build_neighborhood(self, k, order):
        """
        springs: list of (i, j) pairs
        k: initial node
        order: neighborhood depth

        returns: set of node IDs in N^order(k)
        """
        neighborhood = set([k])
        frontier = set([k])
        springs = self.springs
        for _ in range(order):
            new_frontier = set()
            for i, j in springs:
                if i in frontier and j not in neighborhood:
                    new_frontier.add(j)
                elif j in frontier and i not in neighborhood:
                    new_frontier.add(i)
            neighborhood.update(new_frontier)
            frontier = new_frontier

        return neighborhood
    
    def triangle_neighbor(self, i):
        return [s for s in self.simplices() if i in s]
    
    def effective(self, traction_mat, torsion_mat, order):
        """
        Compute M_eff^(order) and K_eff^(order) around node k.

        Parameters:
        - springs: list or set of (i, j) node ID pairs
        - k: contact node ID (int)
        - order: neighborhood depth (int)
        - mass_node: mass of each node (assumed constant)
        - spring_k: stiffness of each spring (assumed constant)

        Returns:
        - M_eff: total effective mass
        - K_eff: total effective stiffness
        """
        # Step 1: build the neighborhood
        k = self.n-1
        springs = self.springs
        # traction_mat = self.traction_mat()
        mass_node = self.mass_node
        neighborhood = self.build_neighborhood(k, order)

        # Step 2: count mass
        M_eff = len(neighborhood) * mass_node + self.impact_mass
        # print(M_eff)
        # Step 3: sum stiffness of springs touching neighborhood
        neighborhood = self.build_neighborhood(k, order-1)
        K_eff = 0.0
        for i, j in springs:
            if i in neighborhood or j in neighborhood:
                # print(i,j)
                K_eff += traction_mat[i,j]
        
        # triangle_contain_x0 = [s for s in self.simplices() if k in s]
        # G_eff = 0.0
        # for triangle in triangle_contain_x0:
        #     i,j,k = triangle
        #     G_eff += torsion_mat[i, j, k]   # angle at j
        #     G_eff += torsion_mat[j, k, i]   # angle at k
        #     G_eff += torsion_mat[k, i, j]   # angle at i
                    
                    
        return M_eff, K_eff #+ G_eff
    
    
    # def fractures_admissible(self):
    #     G = nx.Graph()
    #     Border_edges = self.border_edges_index()
    #     Border_nodes = self.border_nodes_index()
    #     l = list( range(self.n) )
    #     triangle_list = [list(e) for e in self.simplices()]
    #     G.add_nodes_from(l)
    #     G.add_edges_from(self.springs)

    #     Fractures_admissible = []
    #     for i, j in combinations(self.border_nodes_index(), 2):
    #         for path in nx.all_simple_edge_paths(G, i, j, cutoff=(self.n-1)):
    #             if len(path) in range(3, int(self.n)):
    #                 if (tuple(sorted((path[0]))) in Border_edges) and (tuple(sorted(path[-1])) in Border_edges):
    #                     l = path[0]
    #                     for k in range(1, len(path)):
    #                         l += path[k]
    #                     l = list(l)
    #                     # print(l)
    #                     if all(element not in Border_nodes for element in l[3:-3]):
    #                         l = [sorted((l[i], l[i+2], l[i+3]))
    #                              for i in range(0, len(l)-2, 2)]
    #                         # print(l)
    #                         if all([l[i] in self.simplices() for i in range(len(l))]):
    #                             # print(l)
    #                             l = [
    #                                 l[i] in triangle_list for i in range(len(l))]
    #                             # print(l, all(l))
    #                             if all(l):
    #                                 Fractures_admissible.append(path)

        # def Thales_edge(t1, t2):
        #     """
        #     Parameters
        #     ----------
        #     t1, t2: first and second edge of a triangle
        #     Returns
        #     -------
        #     the third edge of triangle
        #     """
        #     return (t1[0], t2[1])

        # def fractures_length(Fractures_admissible):
        #     length = np.zeros(len(Fractures_admissible))
        #     for j in range(len(Fractures_admissible)):
        #         l = [Thales_edge(Fractures_admissible[j][i], Fractures_admissible[j][i+1])
        #              for i in range(len(Fractures_admissible[j])-1)]
        #         l = [0.5 * Spring(self.nodes[l[i][0]], self.nodes[l[i]
        #                           [1]], None).L0 for i in range(len(l))]
        #         length[j] += sum(l)
        #     return length

        # return Fractures_admissible, fractures_length(Fractures_admissible)

    # def Route(self):
    #     g = UndirectedGraph(self.n)
    #     for v1, v2 in self.springs:
    #         g.AddEdge(v1, v2)
    #     while True:
    #         try:
    #             Route = g.RouteInspection()
    #             break
    #         except:
    #             continue
    #     return Route

    def get_nodes(self):
        """return position of all nodes"""
        return np.array([node.position() for node in self.nodes])

    def simplices(self):
        """
        Returns all triangles (simplices) in the mass-spring system
        -------
        TYPE
            DESCRIPTION.
        """
        Points = self.get_nodes()
        tri = Delaunay(Points)
        return np.sort(tri.simplices)

    def get_velocity(self):
        """
        informations of all velocites
        """
        return np.array([node.velocity() for node in self.nodes])

    def connexe_mat(self):
        """ Connectivity matrix """

        mat = np.zeros((self.n, self.n))
        for (i, j) in self.springs:
            mat[i, j] = 1
            mat[j, i] = mat[i, j]
        mat = coo_matrix(mat)
        return mat.todok()

    def length_mat(self):
        """ initial length of traction spring """
        mat = np.zeros((self.n, self.n))
        for (i, j) in self.springs:
            mat[i, j] = Spring(self.nodes[i], self.nodes[j], None).L0
            mat[j, i] = mat[i, j]
        mat = coo_matrix(mat)
        return mat.todok()

    def angle_init(self):
        """ Initial angle of torsion spring """

        mat = np.zeros((self.n, self.n, self.n))
        Nodes_positions = self.get_nodes()
        for i, j, k in self.simplices():
            mat[i, j, k] = Angle(Nodes_positions[i],
                                 Nodes_positions[j], Nodes_positions[k])
            mat[k, j, i] = mat[i, j, k]

            mat[i, k, j] = Angle(Nodes_positions[i],
                                 Nodes_positions[k], Nodes_positions[j])
            mat[j, k, i] = mat[i, k, j]

            mat[j, i, k] = Angle(Nodes_positions[j],
                                 Nodes_positions[i], Nodes_positions[k])
            mat[k, i, j] = mat[j, i, k]
        return mat

    def torsion_mat(self):
        """ stiffness constant of every torsion spring """
        G = self.G
        mat = np.zeros((self.n, self.n, self.n))
        for i, j, k in self.simplices():
            mat[i, j, k] = G * self.length_mat()[i, j] * self.length_mat()[j,
                                                                           k] / sin(self.angle_init()[i, j, k])
            mat[k, j, i] = mat[i, j, k]

            mat[i, k, j] = G * self.length_mat()[i, k] * self.length_mat()[j,
                                                                           k] / sin(self.angle_init()[i, k, j])
            mat[j, k, i] = mat[i, k, j]

            mat[j, i, k] = G * self.length_mat()[i, j] * self.length_mat()[i,
                                                                           k] / sin(self.angle_init()[j, i, k])
            mat[k, i, j] = mat[j, i, k]

        return mat

    def traction_mat(self):
        """ stiffness constant of traction spring """
        mat = np.zeros((self.n, self.n))
        for i, j, k in self.simplices():
            # start = time.time()

            mat[i, j] += self.k / sin(self.angle_init()[i, k, j])
            mat[j, i] = mat[i, j]
            # end = time.time()

            mat[i, k] += self.k / sin(self.angle_init()[i, j, k])
            mat[k, i] = mat[i, k]

            mat[j, k] += self.k / sin(self.angle_init()[k, i, j])
            mat[k, j] = mat[j, k]
        mat = coo_matrix(mat)
        return mat.todok()


    def Neighbors(self):
        edge_set = self.springs
        edge_array = np.array(list(edge_set), dtype=int)
        # print(edge_array)
        reverse_edges = edge_array[:, [1, 0]]
        all_edges = np.vstack((edge_array, reverse_edges))

        # Number of nodes (infer from max index)
        num_nodes = all_edges.max() + 1

        # Sort edges by the source node
        sorted_idx = np.argsort(all_edges[:, 0])
        all_edges = all_edges[sorted_idx]

        # Group neighbors by node
        unique_nodes, start_idx = np.unique(all_edges[:, 0], return_index=True)
        neighbors = []
        for i in range(num_nodes):
            if i in unique_nodes:
                idx = np.where(unique_nodes == i)[0][0]
                start = start_idx[idx]
                end = start_idx[idx + 1] if idx + \
                    1 < len(start_idx) else len(all_edges)
                neighbors.append(np.unique(all_edges[start:end, 1]))
            else:
                neighbors.append(np.array([], dtype=int))

        # print(neighbors)
        return neighbors

    def neighbors(self, i):
        return self.Neighbors()[i]

    def Move(self, time_end: float, Traction_mat, Length_mat, Torsion_mat, Angle_mat):
        t = np.linspace(0, time_end, N_T)
        mass_nodes = self.mass_nodes
        simplices = self.simplices()
        CM = self.connexe_mat()
        All_pos = self.get_nodes()
        All_vel = self.get_velocity()

        Y0_ = np.empty((2 * self.n, All_pos.shape[1]))
        Y0_[0::2] = All_pos
        Y0_[1::2] = All_vel
        Y0_ = Y0_.flatten()

        # Sol = explicit_euler(System(t, Y0_, Y0_, self.n, CM, Length_mat, self.m, self.mu, Traction_mat, Torsion_mat, Angle_mat, self.simplices()), Y0_, 0, time_end, 1./N_T)
        Sol = solve_ivp(System, [0, time_end], Y0_, t_eval=t,
                        args=(self.n, CM, Length_mat, mass_nodes,
                              self.mu, Traction_mat, Torsion_mat, Angle_mat, simplices), method='BDF')

        return Sol

    def Move_stable_1(self, time_end: float, Traction_mat, Length_mat, Torsion_mat, Angle_mat, contact_node):
        t = np.linspace(0, time_end, N_T)
        CM = self.connexe_mat()
        All_pos = self.get_nodes()
        All_vel = self.get_velocity()
        Y0_ = np.array([])
        for i in range(self.n):
            Y0_ = np.append(Y0_, All_pos[i])
            Y0_ = np.append(Y0_, All_vel[i])

        Sol = solve_ivp(System_stable_1, [0, time_end], Y0_, t_eval=t,
                        args=(self.n, CM, Length_mat, self.m,
                              self.mu, Traction_mat, Torsion_mat, Angle_mat, self.simplices(), contact_node))
        return Sol

    def Move_stable_neighbor(self, time_end: float, Traction_mat, Length_mat, Torsion_mat, Angle_mat, contact_node):
        t = np.linspace(0, time_end, N_T)
        CM = self.connexe_mat()
        All_pos = self.get_nodes()
        All_vel = self.get_velocity()
        Y0_ = np.array([])
        for i in range(self.n):
            Y0_ = np.append(Y0_, All_pos[i])
            Y0_ = np.append(Y0_, All_vel[i])

        Sol = solve_ivp(System_stable_neighbor, [0, time_end], Y0_, t_eval=t,
                        args=(Y0_, self.n, CM, Length_mat, self.m,
                              self.mu, Traction_mat, Torsion_mat, Angle_mat, self.simplices(), contact_node))
        return Sol

    def plot_init(self, figax=None):
        plt.figure()
        for (i, j) in self.springs:

            plt.plot([self.nodes[i].position()[0], self.nodes[j].position()[0]],
                     [self.nodes[i].position()[1], self.nodes[j].position()[1]],  color='blue')
            plt.text(self.nodes[i].position()[0],
                     self.nodes[i].position()[1], self.nodes[i].id)
            plt.text(self.nodes[j].position()[0],
                     self.nodes[j].position()[1], self.nodes[j].id)
        # circle = plt.Circle((0., 0.) , 100. , color = 'black', fill=False)
        # plt.gca().add_patch(circle)
        # plt.tight_layout()
        # plt.gca().set_aspect('equal')
        plt.show()

    def plot_border(self):
        plt.figure()
        l = self.border_edges_index()
        for (i, j) in l:
            plt.plot([self.nodes[i].position()[0], self.nodes[j].position()[0]],
                     [self.nodes[i].position()[1], self.nodes[j].position()[1]],  color='blue')
            plt.text(self.nodes[i].position()[0],
                     self.nodes[i].position()[1], self.nodes[i].id)
            plt.text(self.nodes[j].position()[0],
                     self.nodes[j].position()[1], self.nodes[j].id)

    def plot_displacements(self, time_end):
        """
        Parameters
        ----------
        time_end : .

        Plot the position's evolution of all nodes.

        """
        time = self.Move(time_end, self.traction_mat(), self.length_mat(),
                         self.torsion_mat(), self.angle_init()).t
        solution = self.Move(time_end, self.traction_mat(), self.length_mat(),
                             self.torsion_mat(), self.angle_init()).y
        Index_x = np.arange(0, 4*self.n, 4)
        Index_y = Index_x + 1
        plt.figure()
        for i in Index_x:
            plt.plot(time, solution[i])
        for i in Index_y:
            plt.plot(time, solution[i])
        plt.xlabel("time(s)")

        # plt.legend()

    def plot_displacements_1free(self, time_end: float, Traction_mat, Length_mat, Torsion_mat, Angle_mat, contact_node):
        """
        Parameters
        ----------
        time_end : .

        Plot the position's evolution of all nodes.

        """
        time = self.Move_stable_1(
            time_end, Traction_mat, Length_mat, Torsion_mat, Angle_mat, contact_node).t
        solution = self.Move_stable_1(
            time_end, Traction_mat, Length_mat, Torsion_mat, Angle_mat, contact_node).y
        Index_x = np.arange(0, 4*self.n, 4)
        Index_y = Index_x + 1
        plt.figure()
        for i in Index_x:
            plt.plot(time, solution[i])
        for i in Index_y:
            plt.plot(time, solution[i])
        plt.xlabel("time(s)")

    def plot_displacements_Neighborfree(self, time_end: float, Traction_mat, Length_mat, Torsion_mat, Angle_mat, contact_node):
        """
        Parameters
        ----------
        time_end : .

        Plot the position's evolution of all nodes.

        """
        time = self.Move_stable_neighbor(
            time_end, Traction_mat, Length_mat, Torsion_mat, Angle_mat, contact_node).t
        solution = self.Move_stable_neighbor(
            time_end, Traction_mat, Length_mat, Torsion_mat, Angle_mat, contact_node).y
        Index_x = np.arange(0, 4*self.n, 4)
        Index_y = Index_x + 1
        plt.figure()
        for i in Index_x:
            plt.plot(time, solution[i])
        for i in Index_y:
            plt.plot(time, solution[i])
        plt.xlabel("time(s)")

    def evolution(self, time_end, time, Traction_mat, Length_mat, Torsion_mat, Angle_mat):
        """
        Parameters
        ----------
        time_end.
        At each time's step, save the position of all node as a floe.
        Returns
        -------
        Positions of all nodes at time t .

        """
        assert 0 <= time <= time_end, "time must be inside the time discretisation"
        time = int(time*20/time_end)-1
        Res = self.Move(time_end, Traction_mat, Length_mat,
                        Torsion_mat, Angle_mat).y
        New_Nodes_positions = []
        New_Nodes_velocity = []
        Nodes = []
        for i in range(0, self.n*4, 4):
            # print(i)
            New_Nodes_positions.append(
                np.array([Res[i][time], Res[i+1][time]]))
            New_Nodes_velocity.append(
                np.array([Res[i+2][time], Res[i+3][time]]))
        for i in range(self.n):
            Nodes.append(
                Node(New_Nodes_positions[i], New_Nodes_velocity[i], i))
        New_floe = Floe(nodes=Nodes, springs=self.springs, stiffness=self.k)
        return New_floe

    def energy_evolution(self, time_end):
        Length_mat = self.length_mat()
        Traction_mat = self.traction_mat()
        Torsion_mat = self.torsion_mat()
        Angle_mat = self.angle_init()
        Sol = self.Move(time_end, Traction_mat,
                        Length_mat, Torsion_mat, Angle_mat)

        Traction_energy = np.zeros(len(Sol.t))
        Torsion_energy = np.zeros(len(Sol.t))

        for index in range(len(Sol.t)):
            Sum1 = 0.
            Sum2 = 0.
            for i, j, k in self.simplices():
                Qi = np.array([Sol.y[4*i][index],
                              Sol.y[4*i+1][index]])
                Qj = np.array([Sol.y[4*j][index],
                              Sol.y[4*j+1][index]])
                Qk = np.array([Sol.y[4*k][index],
                              Sol.y[4*k+1][index]])

                Sum1 += 0.5 * (Traction_mat[i, j] * (norm(Qi-Qj) - Length_mat[i, j])**2
                               + Traction_mat[i, k] *
                               (norm(Qi-Qk) - Length_mat[i, k])**2
                               + Traction_mat[j, k] * (norm(Qj-Qk) - Length_mat[j, k])**2)

                Sum2 += 0.5 * (Torsion_mat[i, j, k] * (Angle(Qi, Qj, Qk) - Angle_mat[i, j, k])**2
                               + Torsion_mat[i, k, j] *
                               (Angle(Qi, Qk, Qj) - Angle_mat[i, k, j])**2
                               + Torsion_mat[j, i, k] * (Angle(Qj, Qi, Qk) - Angle_mat[j, i, k])**2)

            Traction_energy[index] = Sum1
            Torsion_energy[index] = Sum2
            Total_energy = Traction_energy + Torsion_energy
        return Traction_energy, Torsion_energy, Total_energy

    def energy_evolution_stable(self, time_end, contact_node):
        Length_mat = self.length_mat()
        Traction_mat = self.traction_mat()
        Torsion_mat = self.torsion_mat()
        Angle_mat = self.angle_init()
        Sol = self.Move_stable_1(time_end, Traction_mat,
                                 Length_mat, Torsion_mat, Angle_mat, contact_node)

        Traction_energy = np.zeros(len(Sol.t))
        Torsion_energy = np.zeros(len(Sol.t))

        for index in range(len(Sol.t)):
            Sum1 = 0.
            Sum2 = 0.
            for i, j, k in self.simplices():
                Qi = np.array([Sol.y[4*i][index],
                              Sol.y[4*i+1][index]])
                Qj = np.array([Sol.y[4*j][index],
                              Sol.y[4*j+1][index]])
                Qk = np.array([Sol.y[4*k][index],
                              Sol.y[4*k+1][index]])

                Sum1 += 0.5 * (Traction_mat[i, j] * (norm(Qi-Qj) - Length_mat[i, j])**2
                               + Traction_mat[i, k] *
                               (norm(Qi-Qk) - Length_mat[i, k])**2
                               + Traction_mat[j, k] * (norm(Qj-Qk) - Length_mat[j, k])**2)

                Sum2 += 0.5 * (Torsion_mat[i, j, k] * (Angle(Qi, Qj, Qk) - Angle_mat[i, j, k])**2
                               + Torsion_mat[i, k, j] *
                               (Angle(Qi, Qk, Qj) - Angle_mat[i, k, j])**2
                               + Torsion_mat[j, i, k] * (Angle(Qj, Qi, Qk) - Angle_mat[j, i, k])**2)

            Traction_energy[index] = Sum1
            Torsion_energy[index] = Sum2
            Total_energy = Traction_energy + Torsion_energy
        return Traction_energy, Torsion_energy, Total_energy

    def update_velocity_node(self, i, new_velocity):
        self.nodes[i] = Node(self.nodes[i].position(),
                             new_velocity, self.nodes[i].id)


class Percussion_Wall:
    def __init__(self, floe: Floe, Wall=2., restitution_coef=0.9, time_end=4., eta=0.01):
        self.t_end = time_end
        self.t = np.linspace(0, time_end, 800)
        self.floe = floe
        self.Wall = Wall
        self.eps = restitution_coef
        self.eta = eta  # collision distance

    def simulation(self):
        Length_mat = self.floe.length_mat()
        Traction_mat = self.floe.traction_mat()
        Torsion_mat = self.floe.torsion_mat()
        Angle_mat = self.floe.angle_init()

        floe = self.floe
        Sol = self.floe.Move(self.t_end, Traction_mat,
                             Length_mat, Torsion_mat, Angle_mat)
        Positions = Sol.y

        Ix = list(range(0, self.floe.n*4, 4))
        All_x_positions = []
        # save all positions of all nodes
        All_positions_velocities = [[]]*self.floe.n*4
        for i in range(self.floe.n*4):
            All_positions_velocities[i] = Positions[i]

        # save positions and velocity of the first simulation
        for i in Ix:
            All_x_positions = np.append(All_x_positions, Positions[i])

        j = 0
        while np.any(All_x_positions > self.Wall - self.eta) and j <= 8:
            m = len(All_positions_velocities[0])
            liste = []
            for i in Ix:
                liste.append(
                    np.where(All_positions_velocities[i] >= self.Wall - self.eta)[0])
            for i in range(len(liste)):
                if len(liste[i]) == 0:
                    liste[i] = [1000]
            for i in range(len(liste)):
                liste[i] = min(liste[i])
            k = min(liste)
            index_nodes_contact = np.where(liste == min(liste))[0][0]
            for i in range(self.floe.n*4):
                All_positions_velocities[i] = All_positions_velocities[i][:k]

            After_shock_floe = floe.evolution(
                self.t_end, self.t_end * ((k-m) % 800)/800., Traction_mat, Length_mat, Torsion_mat, Angle_mat)
            After_shock_floe = Floe(
                After_shock_floe.nodes, After_shock_floe.springs)
            # update velocity of nodes when collision!!!
            velocity_after_shock = -self.eps * \
                After_shock_floe.get_velocity()[index_nodes_contact]
            After_shock_floe.update_velocity_node(
                index_nodes_contact, velocity_after_shock)
            # run simulation again with the floe after shock
            floe = After_shock_floe
            New_position_velocities = After_shock_floe.Move(
                self.t_end, Traction_mat, Length_mat, Torsion_mat, Angle_mat).y
            for i in range(self.floe.n*4):
                All_positions_velocities[i] = np.append(
                    All_positions_velocities[i], New_position_velocities[i])

            All_x_positions = []
            for i in Ix:
                All_x_positions = np.append(
                    All_x_positions, All_positions_velocities[i])
            j += 1

        return All_positions_velocities

    # def simulation_with_fracture(self):
    #     All_positions_velocities = self.simulation()
    #     E_nofrac = Energy_studies(All_positions_velocities, self.floe)[-1]
    #     E_allfrac = Energy_studies_fr(All_positions_velocities, self.floe)[-1]
    #     frac_ind, last_step_bef_frac = Find_frac_index(
    #         E_allfrac, E_nofrac)[-2:]

    #     for i in range(self.floe.n*4):
    #         All_positions_velocities[i] = All_positions_velocities[i][:last_step_bef_frac+1]

    #     frac = self.floe.fractures_admissible()[0][frac_ind]
    #     frac = [tuple(set(frac[i])) for i in range(len(frac))]
    #     Springs = [el for el in self.floe.springs if el not in frac]

    #     G = nx.Graph()
    #     G.add_edges_from(Springs)
    #     Springs_new1, Springs_new2 = nx.biconnected_component_edges(G)
    #     Springs_new1, Springs_new2 = set(Springs_new1), set(Springs_new2)

    #     # create 2 new floes in order to respect the fracture
    #     IndexNewfloe1 = list([el for el in nx.connected_components(G)][0])
    #     IndexNewfloe2 = list([el for el in nx.connected_components(G)][1])

    #     Points_new = [np.array([All_positions_velocities[4*i][-1],
    #                             All_positions_velocities[4*i+1][-1]]) for i in range(self.floe.n)]
    #     Vel_new = [np.array([All_positions_velocities[4*i+2][-1],
    #                          All_positions_velocities[4*i+3][-1]]) for i in range(self.floe.n)]
    #     nodes = []
    #     Springs1 = set()
    #     for i, j in Springs_new1:
    #         i = IndexNewfloe1.index(i)
    #         j = IndexNewfloe1.index(j)
    #         Springs1.add((i, j))
    #     for i in range(len(IndexNewfloe1)):
    #         nodes.append(
    #             Node(Points_new[IndexNewfloe1[i]], Vel_new[IndexNewfloe1[i]], IndexNewfloe1[i]))
    #     New_floe1 = Floe(nodes, Springs1, stiffness=self.floe.k,
    #                      viscosity=self.floe.mu, id_number=1)

    #     nodes = []
    #     Springs2 = set()
    #     for i, j in Springs_new2:
    #         i = IndexNewfloe2.index(i)
    #         j = IndexNewfloe2.index(j)
    #         Springs2.add((i, j))
    #     for i in range(len(IndexNewfloe2)):
    #         nodes.append(
    #             Node(Points_new[IndexNewfloe2[i]], Vel_new[IndexNewfloe2[i]], IndexNewfloe2[i]))
    #     New_floe2 = Floe(nodes, Springs2, stiffness=self.floe.k,
    #                      viscosity=self.floe.mu, id_number=2)

    #     Solution1 = Percussion_Wall(New_floe1).simulation()
    #     Solution2 = Percussion_Wall(New_floe2).simulation()

    #     return All_positions_velocities, New_floe1, New_floe2, Solution1, Solution2, last_step_bef_frac

    def position_at_time(self, time_step):
        """
        Parameters
        ----------
        time_step :.
        Returns
        -------
        Positions of each nodes at time_step of simulation
        """
        All_positions_velocities = self.simulation()
        I = list(range(0, self.floe.n*4, 4))

        Pos = [np.array([All_positions_velocities[I[i]][time_step],
                         All_positions_velocities[I[i]+1][time_step]]) for i in range(self.floe.n)]

        Vel = [np.array([All_positions_velocities[I[i]+2][time_step],
                         All_positions_velocities[I[i]+3][time_step]]) for i in range(self.floe.n)]
        return Pos, Vel


def Energy_studies(All_positions_velocities, floe, Length_mat, Angle_mat, Torsion_mat, T_end=1.5):
    """
    Optimized version of Energy_studies function.
    """
    start_time = time.time()
    K = floe.k
    Connect = floe.connexe_mat()
    N_T = len(All_positions_velocities[0])
    
    # Precompute indices for vectorization
    simplices = floe.simplices()
    # num_simplices = len(simplices)
    
    # Preallocate arrays
    Traction_energy = np.zeros(N_T)
    Torsion_energy = np.zeros(N_T)
    
    # Precompute all positions for all time steps
    positions = np.zeros((N_T, len(All_positions_velocities)//4, 2))
    for i in range(len(All_positions_velocities)//4):
        positions[:, i, 0] = All_positions_velocities[4*i]
        positions[:, i, 1] = All_positions_velocities[4*i+1]
    
    for index in range(N_T):
        Traction_en = 0.
        Torsion_en = 0.
        
        for i, j, k in simplices:
            Qi = positions[index, i]
            Qj = positions[index, j]
            Qk = positions[index, k]
            
            # Vector differences
            dij = Qi - Qj
            dik = Qi - Qk
            djk = Qj - Qk
            
            # Lengths
            l_ij = norm(dij)
            l_ik = norm(dik)
            l_jk = norm(djk)
            
            # Traction energy terms
            if Connect[i, j]:
                Traction_en += 0.5 * (K / np.sin(Angle_mat[i, k, j])) * (l_ij - Length_mat[i, j])**2
            if Connect[i, k]:
                Traction_en += 0.5 * (K / np.sin(Angle_mat[i, j, k])) * (l_ik - Length_mat[i, k])**2
            if Connect[j, k]:
                Traction_en += 0.5 * (K / np.sin(Angle_mat[j, i, k])) * (l_jk - Length_mat[j, k])**2
            
            # Angle calculations
            angle_ijk = Angle(Qi, Qj, Qk)
            angle_ikj = Angle(Qi, Qk, Qj)
            angle_jik = Angle(Qj, Qi, Qk)
            
            # Torsion energy terms
            Torsion_en += 0.5 * (Torsion_mat[i, j, k] * (angle_ijk - Angle_mat[i, j, k])**2 +
                                 Torsion_mat[i, k, j] * (angle_ikj - Angle_mat[i, k, j])**2 +
                                 Torsion_mat[j, i, k] * (angle_jik - Angle_mat[j, i, k])**2)
        
        Traction_energy[index] = Traction_en
        Torsion_energy[index] = Torsion_en
    
    E_tot = Traction_energy + Torsion_energy
    end_time = time.time()
    print("time of computation = ", start_time-end_time)
    return Traction_energy, Torsion_energy, E_tot

def Angular_speed_list(All_positions_velocities, floe, T_end=1.2, N_T=N_T):
    n = len(All_positions_velocities[0])
    # wi_list = []
    wj_list = np.zeros(n-1)
    # wk_list = []
    angle = floe.angle_init()
    torsion = Torsion_mat(floe)
    length = Length_mat(floe)
    h = T_end/N_T
    print("h=", h)
    theta = np.zeros(n-1)
    for index in range(n-1):
        # wi, wj, wk = 0., 0., 0.
        wj = 0.
        # v[index] = 0.5*norm(np.array([All_positions_velocities[10][index],
        #               All_positions_velocities[11][index]]))**2
        for i, j, k in floe.simplices():
            print(i, j, k, index)
            Qi_old = np.array([All_positions_velocities[4*i][index],
                               All_positions_velocities[4*i+1][index]])

            Qj_old = np.array([All_positions_velocities[4*j][index],
                               All_positions_velocities[4*j+1][index]])

            Qk_old = np.array([All_positions_velocities[4*k][index],
                               All_positions_velocities[4*k+1][index]])

            Qi_new = np.array([All_positions_velocities[4*i][index+1],
                               All_positions_velocities[4*i+1][index+1]])

            Qj_new = np.array([All_positions_velocities[4*j][index+1],
                               All_positions_velocities[4*j+1][index+1]])

            Qk_new = np.array([All_positions_velocities[4*k][index+1],
                               All_positions_velocities[4*k+1][index+1]])

            # r_ij = norm(Qi_old-Qj_old)
            r_ik = norm(Qi_old-Qk_old)
            # r_jk = norm(Qj_old-Qk_old)

            print(Qi_new, Qj_new, Qk_new, Qi_old, Qj_old, Qk_old)
            print(Angle(Qk_new, Qi_new, Qj_new), Angle(Qk_old, Qi_old, Qj_old))
            print(Angle(Qk_new, Qi_new, Qj_new) -
                  Angle(Qk_old, Qi_old, Qj_old))
            # wi += (Angle(Qj_new, Qi_new, Qk_new) - Angle(Qj_old, Qi_old, Qk_old))/h
            # wj += (Angle(Qi_new, Qj_new, Qk_new) - Angle(Qi_old, Qj_old, Qk_old))/h
            wj = (Angle(Qk_new, Qi_new, Qj_new) -
                  Angle(Qk_old, Qi_old, Qj_old))/(h*r_ik)
            print(wj)
            # wk += (Angle(Qj_new, Qi_new, Qk_new) - Angle(Qi_old, Qj_old, Qk_old))/h
            # wi_list.append( wi**2 * (r_ij**2 + r_ik**2 ))
            # wj_list.append( 0.5*wj**2 * (r_ij**2))
            wj_list[index] = np.rad2deg(wj)
            theta[index] = np.rad2deg(
                Angle(Qk_new, Qi_new, Qj_new) - angle[k, i, j])
            # R[index]     = 0.5 * K[i,k]*(r_ik - 1.) ** 2
            # wk_list.append(wk)

    plt.figure()
    t = np.linspace(0, T_end, N_T-1)
    plt.plot(t, wj_list, label='rotation speed')
    plt.plot(t, theta, label='angle difference')
    plt.legend()

    kinematic = 0.5*length[0, 2]*np.deg2rad(wj_list)**2
    elastic = 0.5*torsion[k, i, j] * np.deg2rad(theta) ** 2

    total = kinematic + elastic  # + R + v

    plt.figure()
    plt.plot(t, elastic, label='elastic torsion')
    plt.plot(t, kinematic, label='kinematic torsion')
    plt.plot(t, total, label='total torsion')
    plt.legend()
    plt.show()
    return theta, wj_list

def compute_kinetic_energy(All_positions_velocities, floe):
    """
    Compute kinetic energy for the mass-spring network.
    
    Kinetic energy: E_kin = 0.5 * sum_i m_i * ||v_i||^2
    
    Parameters:
    - All_positions_velocities: Array of positions and velocities [x0, y0, vx0, vy0, x1, y1, ...]
    - floe: Network object with node masses
    
    Returns:
    - kinetic_energy: Array of kinetic energy values at each time step
    """
    
    # N_T = len(All_positions_velocities[0])  # Number of time steps
    n_nodes = floe.n  # Number of nodes
    
    # Preallocate kinetic energy array
    kinetic_energy = np.zeros(N_T)
    
    # Extract velocities for all time steps [time, node, velocity_component]
    velocities = np.zeros((N_T, n_nodes, 2))
    for i in range(n_nodes):
        velocities[:, i, 0] = All_positions_velocities[4*i + 2]  # vx
        velocities[:, i, 1] = All_positions_velocities[4*i + 3]  # vy
    
    # Compute kinetic energy at each time step
    for time_index in range(N_T):
        total_ke = 0.0
        
        for i in range(n_nodes):
            # Get velocity vector for node i at this time
            vx = velocities[time_index, i, 0]
            vy = velocities[time_index, i, 1]
            
            # Kinetic energy for this node: 0.5 * m * v^2
            node_ke = 0.5 * floe.mass_nodes[i] * (vx**2 + vy**2)
            total_ke += node_ke
        
        kinetic_energy[time_index] = total_ke
    
    return kinetic_energy



def Energy_studies_fr(All_positions_velocities, floe):
    """
    Parameters
    ----------
    All_positions_velocities : Evolution of all nodes and velocities
        DESCRIPTION.
    Returns
    -------
    Traction_energy, Torsion_energy, E_tot when each fracture situation happen
    """
    all_frac, length_frac = floe.fractures_admissible()
    alpha = 1.                                              # ductibility coef
    Length_mat = floe.length_mat()
    Torsion_mat = floe.torsion_mat()
    Angle_mat = floe.angle_init()
    All_Traction_energy = np.zeros(
        (len(all_frac), len(All_positions_velocities[0])))
    All_Torsion_energy = np.zeros_like(All_Traction_energy)
    triangle_list = [list(triangle) for triangle in floe.simplices()]

    for l in range(len(all_frac)):
        Triangle_after_eliminate = [
            el for el in triangle_list if el not in frac_triangle(all_frac[l])]

        for index in range(len(All_positions_velocities[0])):
            Sum1 = 0.
            Sum2 = 0.
            for i, j, k in Triangle_after_eliminate:
                Qi = np.array([All_positions_velocities[4*i][index],
                               All_positions_velocities[4*i+1][index]])
                Qj = np.array([All_positions_velocities[4*j][index],
                               All_positions_velocities[4*j+1][index]])
                Qk = np.array([All_positions_velocities[4*k][index],
                               All_positions_velocities[4*k+1][index]])

                Sum1 += 0.5 * ((floe.k/sin(Angle_mat[i, k, j])) * (norm(Qi-Qj) - Length_mat[i, j])**2
                               + (floe.k/sin(Angle_mat[i, j, k])) *
                               (norm(Qi-Qk) - Length_mat[i, k])**2
                               + (floe.k/sin(Angle_mat[j, i, k])) * (norm(Qj-Qk) - Length_mat[j, k])**2)

                Sum2 += 0.5 * (Torsion_mat[i, j, k] * (Angle(Qi, Qj, Qk) - Angle_mat[i, j, k])**2
                               + Torsion_mat[i, k, j] *
                               (Angle(Qi, Qk, Qj) - Angle_mat[i, k, j])**2
                               + Torsion_mat[j, i, k] * (Angle(Qj, Qi, Qk) - Angle_mat[j, i, k])**2)

            All_Traction_energy[l][index] = Sum1 + alpha * length_frac[l]
            All_Torsion_energy[l][index] = Sum2

    All_E_tot = All_Traction_energy + All_Torsion_energy
    return All_Traction_energy, All_Torsion_energy, All_E_tot


def Find_frac_index(E_tot_frac, E_tot):
    """
    Parameters
    ----------
    E_tot_frac : array (n* nb of time step) of all possible energy when fracture,
    E_tot_frac[i] correspond to total energy if fracture number i happen .
    E_tot : total energy with no fracture.

    Returns
    -------
    The index of fracture that happen, i.e i0 in {0,1,..,n-1}

    """

    ar = np.zeros(len(E_tot_frac)) + len(E_tot)
    for i in range(len(E_tot_frac)):
        if len(np.where(E_tot_frac[i] < E_tot)[0]) != 0:
            ar[i] = np.where(E_tot_frac[i] < E_tot)[0][0]
    time_step_frac = min(ar)
    frac_number = np.where(ar == time_step_frac)[0][0]

    return ar, frac_number, int(time_step_frac)


def frac_triangle(l):
    """
    Parameters
    ----------
    l : fracture road.

    Returns
    -------
    simplices that have to be remove when compute energy

    """
    triangle_to_del = []
    for i in range(len(l)-1):
        triangle_to_del.append(list(set(l[i]+l[i+1])))

    return triangle_to_del


def node_to_node(node1: Node, node2: Node):
    position1 = node1.position()
    position2 = node2.position()
    return norm(position1-position2)


def node_to_floe(node: Node, floe: Floe):
    """
    Parameters
    ----------
    node : Node
    floe : Floe

    Returns
    -------
    float
        Distance between a node and a floe.

    """
    dist = [norm(node.position() - floe.nodes[i].position())
            for i in range(floe.n)]
    return min(dist)


def Unit_vect(vect1, vect2):
    if (vect1[0] == vect2[0] and vect1[1] == vect2[1]):
        return 0.
    return (vect2-vect1)/norm(vect2-vect1)


def Orthogonal_vect(v):
    return np.array([-v[1], v[0]])


def Connexe_mat(floe: Floe):
    """ Connectivity matrix """

    mat = np.zeros((floe.n, floe.n))
    for (i, j) in floe.springs:
        mat[i, j] = 1
        mat[j, i] = mat[i, j]
    mat = coo_matrix(mat)
    return mat.todok()

def Angle(A,B,C):
    """
    find angle by dot product: 
        theta = acos(u.v/(normu.normv))
    """
    u = A - B
    v = C - B
    cos_theta = np.dot(u, v) / (norm(u) * norm(v))
    return acos(cos_theta)

def Torsion_mat(floe: Floe):
    """ stiffness constant of every torsion spring """
    mat = np.zeros((floe.n, floe.n, floe.n))
    G = floe.G
    simplices = floe.simplices()
    lengths = floe.length_mat()
    angle_init = floe.angle_init()

    for i, j, k in simplices:
        mat[i, j, k] = G * lengths[i, j] * lengths[j,
                                                   k] / sin(angle_init[i, j, k])
        mat[k, j, i] = mat[i, j, k]

        mat[i, k, j] = G * lengths[i, k] * lengths[j,
                                                   k] / sin(angle_init[i, k, j])
        mat[j, k, i] = mat[i, k, j]

        mat[j, i, k] = G * lengths[i, j] * lengths[i,
                                                   k] / sin(angle_init[j, i, k])
        mat[k, i, j] = mat[j, i, k]
    return mat


def Traction_mat(floe: Floe, angle_mat):
    mat = np.zeros((floe.n, floe.n))
    for i, j, k in floe.simplices():
        mat[i, j] += floe.k / sin(angle_mat[i, k, j])
        mat[j, i] = mat[i, j]

        mat[i, k] += floe.k / sin(angle_mat[i, j, k])
        mat[k, i] = mat[i, k]

        mat[j, k] += floe.k / sin(angle_mat[k, i, j])
        mat[k, j] = mat[j, k]
    mat = coo_matrix(mat)
    return mat.todok()


def Angle_mat(floe: Floe):
    mat = np.zeros((floe.n, floe.n, floe.n))
    Nodes_positions = floe.get_nodes()
    for i, j, k in floe.simplices():
        mat[i, j, k] = Angle(Nodes_positions[i],
                             Nodes_positions[j], Nodes_positions[k])
        mat[k, j, i] = mat[i, j, k]

        mat[i, k, j] = Angle(Nodes_positions[i],
                             Nodes_positions[k], Nodes_positions[j])
        mat[j, k, i] = mat[i, k, j]

        mat[j, i, k] = Angle(Nodes_positions[j],
                             Nodes_positions[i], Nodes_positions[k])
        mat[k, i, j] = mat[j, i, k]
    return mat


def Length_mat(floe: Floe):
    mat = np.zeros((floe.n, floe.n))
    for (i, j) in floe.springs:
        mat[i, j] = Spring(floe.nodes[i], floe.nodes[j], None).L0
        mat[j, i] = mat[i, j]
    mat = coo_matrix(mat)
    return mat.todok()


"""
Main system describes all nodes evolutions. 
Each node i depends on TRACTION's spring (i.e spring between node i and node j neighbor)
and TORSION's spring ( formed by the triangle it belongs to)
"""


def mass_nodes(floe: Floe):
    mass_nodes = floe.mass_nodes
    return mass_nodes

# def grad_theta_i_all(qi, qj, qk):
#     """
#     Compute the gradient of angle at qi with respect to qi, qj, qk
#     Returns: (∇_qi θi, ∇_qj θi, ∇_qk θi)
#     """
#     u = qj - qi
#     v = qk - qi
#     ru = norm(u)
#     rv = norm(v)

#     # if ru < 1e-12 or rv < 1e-12:
#     #     zero = np.zeros_like(qi)
#     #     return zero, zero, zero
    
#     dot_uv = np.dot(u, v)
#     costheta = dot_uv / (ru * rv)
#     # costheta = np.clip(costheta, -1.0, 1.0)
#     sintheta = np.sqrt(1.0 - costheta**2)


#     # ∇_qj
#     # grad_qj = (v * ru * rv - dot_uv * (rv / ru) * u) / (ru**2 * rv**2) / -sintheta
#     grad_qj = (dot_uv*u/ru**3*rv - v/(ru*rv))/sintheta
    
#     # ∇_qk
#     # grad_qk = (u * ru * rv - dot_uv * (ru / rv) * v) / (ru**2 * rv**2) / -sintheta
#     grad_qk = (dot_uv*v/rv**3*ru - u/(ru*rv))/sintheta
    
#     # ∇_qi
#     # term1 = -(u + v) * (ru * rv)
#     # term2 = dot_uv * ((rv / ru) * u + (ru / rv) * v)
#     grad_qi = - grad_qj - grad_qk
#     return grad_qi, grad_qj, grad_qk

def grad_theta_i_all(qi, qj, qk):
    """
    Gradient of angle θ_i = angle(qj - qi, qk - qi) using orthogonal unit vector
    """
    u = qj - qi
    v = qk - qi

    ru = norm(u)
    rv = norm(v)

    u_hat = u / ru
    v_hat = v / rv

    # cos_theta = np.dot(u_hat, v_hat)
    # cos_theta = np.clip(cos_theta, -1.0, 1.0)
    # sin_theta = np.sqrt(1 - cos_theta**2)

    # Direction vector orthogonal to u_hat
    u_perp = Orthogonal_vect(u_hat)
    v_perp = Orthogonal_vect(v_hat)

    # Use geometric form:
    grad_qj = - (1 / ru) * u_perp  # Or more precisely: -1 / (ru * sinθ) * component in u^
    grad_qk =  (1 / rv) * v_perp
    grad_qi = - grad_qj - grad_qk

    return grad_qi, grad_qj, grad_qk

def System(t, Y, nb_nodes, Connex_mat, Length_mat, Mass_mat, mu,
           Traction_mat, Torsion_mat, Angle_init, Triangle_list):
    """
    Parameters
    ----------
    t : time discretisation.
    Y : TYPE
        DESCRIPTION.
    Y0 : init condition of system (init position and velocity).
    nb_nodes : number of nodes.
    Connexe_mat : connectivity between nodes's matrix.
    Length_mat : length matrix.
    m : matrix contains mass of each node.
    mu : viscosity const.
    k : stiffness const of traction spring.
    Torsion_mat: stiffness constant of torsion spring at (i,j,k)

    Returns
    -------
    (evolution of node_i, velo_i
      for i in range nb_nodes as a dynamical system).

    """
    # u = np.zeros((nb_nodes, nb_nodes, 2))
    Q = np.reshape(Y, (nb_nodes*2, 2))
    Y_ = np.zeros_like(Q)
    K = Traction_mat
    G = Torsion_mat
    Theta0 = Angle_init
    inv_m = 1./Mass_mat

    # Compute elastic force between q_i and q_j
    for i in range(0, nb_nodes):
        Y_[2*i] = Q[2*i+1]
        for j in range(i+1, nb_nodes):
            u_ij = Unit_vect(Q[2*i], Q[2*j])
            force = Connex_mat[i, j] * (K[i, j] * (norm(Q[2*j]-Q[2*i]) - Length_mat[i, j]) * u_ij
                                        + mu * (Q[2*j+1] - Q[2*i+1]) @ u_ij * u_ij)
            Y_[2*i+1] += inv_m[i] * force
            Y_[2*j+1] -= inv_m[i] * force

    # to verify this part of calculation,
    # using 1 floe of 3 nodes
    # try different initial speed of each node

    for i, j, k in Triangle_list:
        if orientation(Q[2*i], Q[2*j], Q[2*k]) != 1:
            (j, k) = k, j

        Theta_i = Angle(Q[2*j], Q[2*i], Q[2*k])
        Theta_j = Angle(Q[2*i], Q[2*j], Q[2*k])
        Theta_k = Angle(Q[2*i], Q[2*k], Q[2*j])

        G_i, G_j, G_k = G[j, i, k], G[i, j, k], G[i, k, j]
        
        ### Torsional force by gradient of angle
        # Gradients of angle_i
        grad_i_i, grad_i_j, grad_i_k = grad_theta_i_all(Q[2*i], Q[2*j], Q[2*k])
        F_i = - G_i * (Theta_i - Theta0[j, i, k])
        
        # Gradients of angle_j
        grad_j_j, grad_j_k, grad_j_i = grad_theta_i_all(Q[2*j], Q[2*k], Q[2*i])
        F_j = - G_j * (Theta_j - Theta0[i, j, k])
        
        # Gradients of angle_k
        grad_k_k, grad_k_i, grad_k_j = grad_theta_i_all(Q[2*k], Q[2*i], Q[2*j])
        F_k = - G_k * (Theta_k - Theta0[i, k, j])
        
        #update forces
        Y_[2*i+1] += inv_m[i] * (F_i * grad_i_i + F_j * grad_j_i + F_k * grad_k_i)
        Y_[2*j+1] += inv_m[j] * (F_i * grad_i_j + F_j * grad_j_j + F_k * grad_k_j)
        Y_[2*k+1] += inv_m[k] * (F_i * grad_i_k + F_j * grad_j_k + F_k * grad_k_k)
        
        #### Force independant of traction's length
        # unit vector
        # u_ij = Unit_vect(Q[2*i], Q[2*j])
        # u_ji = -u_ij
        # u_ik = Unit_vect(Q[2*i], Q[2*k])
        # u_ki = -u_ik
        # u_jk = Unit_vect(Q[2*j], Q[2*k])
        # u_kj = -u_jk
        # # print(u_ij, u_ik, u_jk)
        # # length
        # l_ij = norm(Q[2*j]-Q[2*i])
        # l_ik = norm(Q[2*k]-Q[2*i])
        # l_kj = norm(Q[2*j]-Q[2*k])

        # Theta_i = Angle(Q[2*j], Q[2*i], Q[2*k])
        # Theta_j = Angle(Q[2*i], Q[2*j], Q[2*k])
        # Theta_k = Angle(Q[2*i], Q[2*k], Q[2*j])
        # force_i = G_j * ((Theta_j - Theta0[i, j, k]) * Orthogonal_vect(u_ij) / l_ij) + G_k * (
        #     (Theta_k - Theta0[i, k, j]) * Orthogonal_vect(u_ki) / l_ik)

        # force_j = G_i * ((Theta_i - Theta0[j, i, k]) * Orthogonal_vect(u_ij) / l_ij) + G_k * (
        #     (Theta_k - Theta0[i, k, j]) * Orthogonal_vect(u_jk) / l_kj)

        # force_k = G_i * ((Theta_i - Theta0[j, i, k]) * Orthogonal_vect(
        #     u_ki) / l_ik) + G_j * ((Theta_j - Theta0[i, j, k]) * Orthogonal_vect(u_jk)/l_kj)

        # Y_[2*i+1] += inv_m[i] * force_i  # Force on node i
        # Y_[2*j+1] += inv_m[j] * force_j  # Force on node j
        # Y_[2*k+1] += inv_m[k] * force_k  # Force on node k

    return np.reshape(Y_, (nb_nodes * 4))


def System_stable_1(t, Y, Y0, nb_nodes, Connex_mat, Length_mat, m, mu, Traction_mat, Torsion_mat, Angle_init, Triangle_list, contact_node):
    """
    Parameters
    ----------
    t : time discretisation.
    Y : TYPE
        DESCRIPTION.
    Y0 : init condition of system (init position and velocity).
    nb_nodes : number of nodes.
    Connexe_mat : connectivity between nodes's matrix.
    Length_mat : length matrix.
    m : mass of each node.
    mu : viscosity const.
    k : stiffness const of traction spring.
    Torsion_mat: stiffness constant of torsion spring at (i,j,k)
    contac_nodes: Nodes in contact with objects

    Returns
    -------
    (evolution of node_i, velo_i
      for i in range nb_nodes as a dynamical system).

    """
    u = np.zeros((nb_nodes, nb_nodes, 2))
    Q = np.reshape(Y, (nb_nodes*2, 2))
    Y_ = np.zeros_like(Q)
    k = Traction_mat
    G = Torsion_mat
    Theta0 = Angle_init
    inv_m = 1./m

    # find the neigborhood of contact node
    Neighbor_contact = [np.any(Triangle_list[i] == contact_node)
                        for i in range(len(Triangle_list))]
    Neighbor_contact = [i for i, val in enumerate(
        Neighbor_contact) if val is True]
    Neighbor_contact = Triangle_list[Neighbor_contact]

    for i in range(nb_nodes):
        if i == contact_node:
            Y_[2*i] = Q[2*i+1]
            for j in range(i+1, i+nb_nodes):
                j = j % nb_nodes
                u_ij = Unit_vect(Q[2*i], Q[2*j])
                Y_[2*i+1] += inv_m * Connex_mat[i, j] * (k[i, j] * (norm(Q[2*j]-Q[2*i]) - Length_mat[i, j]) * u_ij
                                                         + mu * (Q[2*j+1] - Q[2*i+1]) @ u_ij * u_ij)

    for i, j, k in Triangle_list:
        u_ij = Unit_vect(Q[2*i], Q[2*j])
        u_ji = -u_ij
        u_ik = Unit_vect(Q[2*i], Q[2*k])
        u_ki = -u_ik
        u_jk = Unit_vect(Q[2*j], Q[2*k])
        u_kj = -u_jk
    # all nodes stable, only the contact node can move
        if i == contact_node or j == contact_node or k == contact_node:
            Y_[2*contact_node+1] += inv_m * (G[i, j, k] * (Angle(Q[2*i], Q[2*j], Q[2*k]) - Theta0[i, j, k]) /
                                             (norm(Q[2*i] - Q[2*j])) * u_ik
                                             + G[i, k, j] * (Angle(Q[2*i], Q[2*k], Q[2*j]) - Theta0[i, k, j]) /
                                             norm(Q[2*i] - Q[2*k]) * u_ij)

    return np.reshape(Y_, (nb_nodes*4))


def orientation(A, B, C):
    """
    Returns the orientation of the triplet (A, B, C).
    - +1: Counterclockwise (CCW)
    - -1: Clockwise (CW)
    -  0: Collinear
    """
    cross = (B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0])
    if cross > 0:
        return 1  # CCW
    elif cross < 0:
        return -1  # CW
    else:
        return 0  # Collinear


def System_stable_neighbor(t, Y, Y0, nb_nodes, Connex_mat, Length_mat, m, mu, Traction_mat, Torsion_mat, Angle_init, Triangle_list, contact_node):
    """
    Parameters
    ----------
    t : time discretisation.
    Y : TYPE
        DESCRIPTION.
    Y0 : init condition of system (init position and velocity).
    nb_nodes : number of nodes.
    Connexe_mat : connectivity between nodes's matrix.
    Length_mat : length matrix.
    m : mass of each node.
    mu : viscosity const.
    k : stiffness const of traction spring.
    Torsion_mat: stiffness constant of torsion spring at (i,j,k)
    contac_nodes: Nodes in contact with objects

    Returns
    -------
    (evolution of node_i, velo_i
      for i in range nb_nodes as a dynamical system for a neighborhood of contact node).
    """
    # not done yet!!!

    Q = np.reshape(Y, (nb_nodes*2, 2))
    Y_ = np.zeros_like(Q)
    k = Traction_mat
    inv_m = 1./m
    # find the neigborhood of contact node
    Neighbor_contact = [np.any(Triangle_list[i] == contact_node)
                        for i in range(len(Triangle_list))]
    Neighbor_contact = [i for i, val in enumerate(
        Neighbor_contact) if val is True]
    Neighbor_contact = Triangle_list[Neighbor_contact].tolist()
    list_contact = []
    for i in range(len(Neighbor_contact)):
        list_contact = list_contact + Neighbor_contact[i]

    for i in list_contact:
        # if i == contact_node:
        Y_[2*i] = Q[2*i+1]
        for j in range(i+1, i+nb_nodes):

            j = j % nb_nodes
            u_ij = Unit_vect(Q[2*i], Q[2*j])
            Y_[2*i+1] += inv_m * Connex_mat[i, j] * (k[i, j] * (norm(Q[2*j]-Q[2*i]) - Length_mat[i, j]) * u_ij
                                                     + mu * (Q[2*j+1] - Q[2*i+1]) @ u_ij * u_ij)

    return np.reshape(Y_, (nb_nodes*4))


T_LIMIT = 1850  # (s) limit time of simulation

def compute_dissipated_energy(All_positions_velocities, floe, T_end, mu):
    """
    Compute dissipated energy over time from position data
    """
    n_timesteps = N_T
    n_nodes = floe.n
    dissipated_energy = np.zeros(n_timesteps)
    C = floe.connexe_mat()
    # Precompute all spring length derivatives
    for t in range(1, n_timesteps):
        current_Ed = 0
        
        # Loop over all edges in the network
        for i in range(n_nodes):
            for j in range(i+1, n_nodes):  # Avoid double counting
                if C[i, j] == 1:  # If nodes are connected
                    # Current spring length
                    pos_i = All_positions_velocities[4*i:4*i+2, t]
                    pos_j = All_positions_velocities[4*j:4*j+2, t]
                    current_length = np.linalg.norm(pos_i - pos_j)
                    
                    # Previous spring length  
                    pos_i_prev = All_positions_velocities[4*i:4*i+2, t-1]
                    pos_j_prev = All_positions_velocities[4*j:4*j+2, t-1]
                    prev_length = np.linalg.norm(pos_i_prev - pos_j_prev)
                    
                    # Time derivative (finite difference)
                    dt = T_end/N_T
                    length_derivative = (current_length - prev_length) / dt
                    
                    # Add to dissipated energy (trapezoidal rule integration)
                    current_Ed += mu * (length_derivative ** 2) * dt
        
        dissipated_energy[t] = dissipated_energy[t-1] + current_Ed
    
    return dissipated_energy
def timeout_handler(signum, frame):
    raise TimeoutError("Simulation exceeded time limit and was stopped!")
