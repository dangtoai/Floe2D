import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm
from scipy.integrate import solve_ivp
from scipy.spatial import Voronoi, voronoi_plot_2d, Delaunay
from graph import *
from math import acos, sin
from scipy.sparse import coo_matrix
import networkx as nx
from itertools import combinations
# from collections import Counter
import matplotlib.animation as animation


# 4 classes Node-> Spring-> Floe-> Percussion/Percussion_wall to simulate the collision between a floe and a wall


class Node:
    """A class representing one node of an ice floe"""

    def __init__(self, position, velocity: np.array([0, 0]), id_number=None):
        self.x, self.y = position
        self.x0, self.y0 = position  # Initial position needed for plots
        self.vx, self.vy = velocity
        self.id = id_number
        self.parent_ice_floe = None  # Ice floe to which this node belongs

    def position(self):
        return np.array([self.x, self.y])

    def velocity(self):
        return np.array([self.vx, self.vy])

    def get_details(self):
        return self.position(), self.velocity()


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
        return [self.node1.position(), self.node2.position()]

    def get_details(self):
        return [self.node1.position(), self.node2.position(), self.L0]


class Floe:
    """
    A class representing an ice floe
    """

    def __init__(self, nodes=None, springs=None, mass=1., stiffness=15, viscosity=2.0, tenacity=1.0, id_number=None):
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

        self.m = mass  # mass of each node
        self.k = stiffness
        self.mu = viscosity
        self.L = tenacity
        self.id = id_number

    def generate_springs(self):
        l = []
        for s in self.springs:
            l.append(Spring(self.nodes[s[0]], self.nodes[s[1]], None))
        return l

    def center(self):
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
        l = []
        for i, j in BV:
            l.append(i)
            l.append(j)
        return list(set(l))

    def fractures_admissible(self):
        G = nx.Graph()
        Border_edges = self.border_edges_index()
        Border_nodes = self.border_nodes_index()
        l = [i for i in range(self.n)]
        triangle_list = [list(e) for e in self.simplices()]
        G.add_nodes_from(l)
        G.add_edges_from(self.springs)

        Fractures_admissible = []
        for i, j in combinations(self.border_nodes_index(), 2):
            for path in nx.all_simple_edge_paths(G, i, j, cutoff=(self.n-1)):
                if len(path) in range(3, int(self.n)):
                    if (tuple(sorted((path[0]))) in Border_edges) and (tuple(sorted(path[-1])) in Border_edges):
                        l = path[0]
                        for k in range(1, len(path)):
                            l += path[k]
                        l = list(l)
                        # print(l)
                        if all(element not in Border_nodes for element in l[3:-3]):
                            l = [sorted((l[i], l[i+2], l[i+3]))
                                 for i in range(0, len(l)-2, 2)]
                            # print(l)
                            if all([l[i] in self.simplices() for i in range(len(l))]):
                                # print(l)
                                l = [
                                    l[i] in triangle_list for i in range(len(l))]
                                # print(l, all(l))
                                if all(l):
                                    Fractures_admissible.append(path)

        def Thales_edge(t1, t2):
            """
            Parameters
            ----------
            t1, t2: first and second edge of a triangle 
            Returns
            -------
            the third edge of triangle
            """
            return (t1[0], t2[1])

        def fractures_length(Fractures_admissible):
            length = np.zeros(len(Fractures_admissible))
            for j in range(len(Fractures_admissible)):
                l = [Thales_edge(Fractures_admissible[j][i], Fractures_admissible[j][i+1])
                     for i in range(len(Fractures_admissible[j])-1)]
                l = [0.5 * Spring(self.nodes[l[i][0]], self.nodes[l[i]
                                  [1]], None).L0 for i in range(len(l))]
                length[j] += sum(l)
            return length

        return Fractures_admissible, fractures_length(Fractures_admissible)

    def First_radius(self):
        """
        Returns the radius of the first circle circumscribed of the floe.

        """
        center = np.array(self.center())
        all_radius = [norm(center-node.position()) for node in self.nodes]
        return max(all_radius)

    def radius_n_circle(self, n):
        return self.First_radius()/n

    def center_velocity(self):
        all_velocity = self.get_velocity()
        return sum(all_velocity/self.n)

    def Route(self):
        g = UndirectedGraph(self.n)
        for v1, v2 in self.springs:
            g.AddEdge(v1, v2)
        while True:
            try:
                Route = g.RouteInspection()
                break
            except:
                continue
        return Route

    def get_nodes(self):
        # return position of all nodes
        return [node.position() for node in self.nodes]

    def simplices(self):
        """
        Returns all triangle (simplices) in the mass-spring system
        -------
        TYPE
            DESCRIPTION.
        """
        Points = self.get_nodes()
        tri = Delaunay(Points)
        return np.sort(tri.simplices)

    def get_velocity(self):
        return np.array([node.velocity() for node in self.nodes])

    def get_springs(self):  # not necessary yet ?
        return [(self.get_nodes()[i], self.get_nodes()[j]) for (i, j) in self.springs]

    def connexe_mat(self):
        k = max(max(self.springs))+1
        Mat = np.zeros((k, k))
        for (i, j) in self.springs:
            Mat[i, j] = 1
            Mat[j, i] = Mat[i, j]
        Mat = coo_matrix(Mat)
        return Mat.todok()

    def length_mat(self):
        # k = max(max(self.springs))+1
        # Mat = np.zeros((k, k))
        Mat = np.zeros((self.n, self.n))
        for (i, j) in self.springs:
            Mat[i, j] = Spring(self.nodes[i], self.nodes[j], None).L0
            Mat[j, i] = Mat[i, j]
        Mat = coo_matrix(Mat)
        return Mat.todok()

    def angle_init(self):
        # k = max(max(self.springs))+1
        # Mat = np.zeros((k, k, k))
        Mat = np.zeros((self.n, self.n, self.n))
        Nodes_positions = self.get_nodes()
        for i, j, k in self.simplices():
            Mat[i, j, k] = Angle(Nodes_positions[i],
                                 Nodes_positions[j], Nodes_positions[k])
            Mat[k, j, i] = Mat[i, j, k]

            Mat[i, k, j] = Angle(Nodes_positions[i],
                                 Nodes_positions[k], Nodes_positions[j])
            Mat[j, k, i] = Mat[i, k, j]

            Mat[j, i, k] = Angle(Nodes_positions[j],
                                 Nodes_positions[i], Nodes_positions[k])
            Mat[k, i, j] = Mat[j, i, k]

        return Mat

    def torsion_mat(self):
        """ stiffness constant of every torsion spring """
        # k = max(max(self.springs))+1
        # Mat = np.zeros((k, k, k))
        Mat = np.zeros((self.n, self.n, self.n))
        for i, j, k in self.simplices():
            Mat[i, j, k] = G * self.length_mat()[i, j] * self.length_mat()[j,
                                                                           k] / sin(self.angle_init()[i, j, k])
            Mat[k, j, i] = Mat[i, j, k]

            Mat[i, k, j] = G * self.length_mat()[i, k] * self.length_mat()[j,
                                                                           k] / sin(self.angle_init()[i, k, j])
            Mat[j, k, i] = Mat[i, k, j]

            Mat[j, i, k] = G * self.length_mat()[i, j] * self.length_mat()[i,
                                                                           k] / sin(self.angle_init()[j, i, k])
            Mat[k, i, j] = Mat[j, i, k]
        return Mat

    def traction_mat(self):
        # k = max(max(self.springs))+1
        # Mat = np.zeros((k, k))
        Mat = np.zeros((self.n, self.n))
        for i, j, k in self.simplices():
            Mat[i, j] += self.k / sin(self.angle_init()[i, k, j])
            Mat[j, i] = Mat[i, j]

            Mat[i, k] += self.k / sin(self.angle_init()[i, j, k])
            Mat[k, i] = Mat[i, k]

            Mat[j, k] += self.k / sin(self.angle_init()[k, i, j])
            Mat[k, j] = Mat[j, k]
        Mat = coo_matrix(Mat)
        return Mat.todok()

    def Neighbor_contact(self, contact_node):
        Neighbor_contact = [np.any(self.simplices()[i] == contact_node)
                            for i in range(len(self.simplices()))]
        Neighbor_contact = [i for i, val in enumerate(
            Neighbor_contact) if val == True]
        Neighbor_contact = self.simplices()[Neighbor_contact]
        return Neighbor_contact

    def Move(self, time_end: float, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat):
        N = 800
        t = np.linspace(0, time_end, N)
        CM = self.connexe_mat()
        All_pos = self.get_nodes()
        All_vel = self.get_velocity()
        Y0_ = np.array([])
        for i in range(self.n):
            Y0_ = np.append(Y0_, All_pos[i])
            Y0_ = np.append(Y0_, All_vel[i])

        Sol = solve_ivp(System, [0, time_end], Y0_, t_eval=t,
                        args=(Y0_, self.n, CM, Length_Mat, self.m,
                              self.mu, Traction_Mat, Torsion_Mat, Angle_Mat, self.simplices()))

        return Sol

    def Move_stable_1(self, time_end: float, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node):
        N = 800
        t = np.linspace(0, time_end, N)
        CM = self.connexe_mat()
        All_pos = self.get_nodes()
        All_vel = self.get_velocity()
        Y0_ = np.array([])
        for i in range(self.n):
            Y0_ = np.append(Y0_, All_pos[i])
            Y0_ = np.append(Y0_, All_vel[i])

        Sol = solve_ivp(System_stable_1, [0, time_end], Y0_, t_eval=t,
                        args=(Y0_, self.n, CM, Length_Mat, self.m,
                              self.mu, Traction_Mat, Torsion_Mat, Angle_Mat, self.simplices(), contact_node))
        return Sol

    def Move_stable_neighbor(self, time_end: float, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node):
        N = 800
        t = np.linspace(0, time_end, N)
        CM = self.connexe_mat()
        All_pos = self.get_nodes()
        All_vel = self.get_velocity()
        Y0_ = np.array([])
        for i in range(self.n):
            Y0_ = np.append(Y0_, All_pos[i])
            Y0_ = np.append(Y0_, All_vel[i])

        Sol = solve_ivp(System_stable_neighbor, [0, time_end], Y0_, t_eval=t,
                        args=(Y0_, self.n, CM, Length_Mat, self.m,
                              self.mu, Traction_Mat, Torsion_Mat, Angle_Mat, self.simplices(), contact_node))
        return Sol

    def plot_init(self):
        plt.figure()
        # l = np.array([node.id for node in self.nodes])
        for (i, j) in self.springs:

            plt.plot([self.nodes[i].position()[0], self.nodes[j].position()[0]],
                     [self.nodes[i].position()[1], self.nodes[j].position()[1]],  color='blue')
            plt.text(self.nodes[i].position()[0],
                     self.nodes[i].position()[1], self.nodes[i].id)
            plt.text(self.nodes[j].position()[0],
                     self.nodes[j].position()[1], self.nodes[j].id)

        plt.plot(self.nodes[3].position()[0], self.nodes[3].position()[
                 1], marker="o", markersize=10, markeredgecolor="red", markerfacecolor="green")
        plt.plot(self.nodes[6].position()[0], self.nodes[6].position()[
                 1], marker="o", markersize=10, markeredgecolor="red", markerfacecolor="green")
        plt.plot(self.nodes[1].position()[0], self.nodes[1].position()[
                 1], marker="o", markersize=10, markeredgecolor="red", markerfacecolor="green")
        plt.plot(self.nodes[0].position()[0], self.nodes[0].position()[
                 1], marker="o", markersize=10, markeredgecolor="red", markerfacecolor="green")

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

    def plot_displacements_1free(self, time_end: float, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node):
        """
        Parameters
        ----------
        time_end : .

        Plot the position's evolution of all nodes.

        """
        time = self.Move_stable_1(
            time_end, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node).t
        solution = self.Move_stable_1(
            time_end, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node).y
        Index_x = np.arange(0, 4*self.n, 4)
        Index_y = Index_x + 1
        plt.figure()
        for i in Index_x:
            plt.plot(time, solution[i])
        for i in Index_y:
            plt.plot(time, solution[i])
        plt.xlabel("time(s)")

    def plot_displacements_Neighborfree(self, time_end: float, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node):
        """
        Parameters
        ----------
        time_end : .

        Plot the position's evolution of all nodes.

        """
        time = self.Move_stable_neighbor(
            time_end, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node).t
        solution = self.Move_stable_neighbor(
            time_end, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node).y
        Index_x = np.arange(0, 4*self.n, 4)
        Index_y = Index_x + 1
        plt.figure()
        for i in Index_x:
            plt.plot(time, solution[i])
        for i in Index_y:
            plt.plot(time, solution[i])
        plt.xlabel("time(s)")

    def evolution(self, time_end, time, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat):
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
        time = int(time*800/time_end)-1
        Res = self.Move(time_end, Traction_Mat, Length_Mat,
                        Torsion_Mat, Angle_Mat).y
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
        Length_Mat = self.length_mat()
        Traction_Mat = self.traction_mat()
        Torsion_Mat = self.torsion_mat()
        Angle_Mat = self.angle_init()
        Sol = self.Move(time_end, Traction_Mat,
                        Length_Mat, Torsion_Mat, Angle_Mat)

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

                Sum1 += 0.5 * (Traction_Mat[i, j] * (norm(Qi-Qj) - Length_Mat[i, j])**2
                               + Traction_Mat[i, k] *
                               (norm(Qi-Qk) - Length_Mat[i, k])**2
                               + Traction_Mat[j, k] * (norm(Qj-Qk) - Length_Mat[j, k])**2)

                Sum2 += 0.5 * (Torsion_Mat[i, j, k] * (Angle(Qi, Qj, Qk) - Angle_Mat[i, j, k])**2
                               + Torsion_Mat[i, k, j] *
                               (Angle(Qi, Qk, Qj) - Angle_Mat[i, k, j])**2
                               + Torsion_Mat[j, i, k] * (Angle(Qj, Qi, Qk) - Angle_Mat[j, i, k])**2)

            Traction_energy[index] = Sum1
            Torsion_energy[index] = Sum2
            Total_energy = Traction_energy + Torsion_energy
        return Traction_energy, Torsion_energy, Total_energy

    def energy_evolution_stable(self, time_end, contact_node):
        Length_Mat = self.length_mat()
        Traction_Mat = self.traction_mat()
        Torsion_Mat = self.torsion_mat()
        Angle_Mat = self.angle_init()
        Sol = self.Move_stable_1(time_end, Traction_Mat,
                                 Length_Mat, Torsion_Mat, Angle_Mat, contact_node)

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

                Sum1 += 0.5 * (Traction_Mat[i, j] * (norm(Qi-Qj) - Length_Mat[i, j])**2
                               + Traction_Mat[i, k] *
                               (norm(Qi-Qk) - Length_Mat[i, k])**2
                               + Traction_Mat[j, k] * (norm(Qj-Qk) - Length_Mat[j, k])**2)

                Sum2 += 0.5 * (Torsion_Mat[i, j, k] * (Angle(Qi, Qj, Qk) - Angle_Mat[i, j, k])**2
                               + Torsion_Mat[i, k, j] *
                               (Angle(Qi, Qk, Qj) - Angle_Mat[i, k, j])**2
                               + Torsion_Mat[j, i, k] * (Angle(Qj, Qi, Qk) - Angle_Mat[j, i, k])**2)

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
        Length_Mat = self.floe.length_mat()
        Traction_Mat = self.floe.traction_mat()
        Torsion_Mat = self.floe.torsion_mat()
        Angle_Mat = self.floe.angle_init()

        floe = self.floe
        Sol = self.floe.Move(self.t_end, Traction_Mat,
                             Length_Mat, Torsion_Mat, Angle_Mat)
        Positions = Sol.y

        Ix = [j for j in range(0, self.floe.n*4, 4)]
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
            print(k)
            for i in range(self.floe.n*4):
                All_positions_velocities[i] = All_positions_velocities[i][:k]

            After_shock_floe = floe.evolution(
                self.t_end, self.t_end * ((k-m) % 800)/800., Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat)
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
                self.t_end, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat).y
            for i in range(self.floe.n*4):
                All_positions_velocities[i] = np.append(
                    All_positions_velocities[i], New_position_velocities[i])

            All_x_positions = []
            for i in Ix:
                All_x_positions = np.append(
                    All_x_positions, All_positions_velocities[i])
            j += 1

        return All_positions_velocities

    def simulation_with_fracture(self):
        All_positions_velocities = self.simulation()
        E_nofrac = Energy_studies(All_positions_velocities, self.floe)[-1]
        E_allfrac = Energy_studies_fr(All_positions_velocities, self.floe)[-1]
        frac_ind, last_step_bef_frac = Find_frac_index(
            E_allfrac, E_nofrac)[-2:]

        for i in range(self.floe.n*4):
            All_positions_velocities[i] = All_positions_velocities[i][:last_step_bef_frac+1]

        frac = self.floe.fractures_admissible()[0][frac_ind]
        frac = [tuple(set(frac[i])) for i in range(len(frac))]
        Springs = [el for el in self.floe.springs if el not in frac]

        G = nx.Graph()
        G.add_edges_from(Springs)
        Springs_new1, Springs_new2 = nx.biconnected_component_edges(G)
        Springs_new1, Springs_new2 = set(Springs_new1), set(Springs_new2)

        # create 2 new floes in order to respect the fracture
        IndexNewfloe1 = list([el for el in nx.connected_components(G)][0])
        IndexNewfloe2 = list([el for el in nx.connected_components(G)][1])

        Points_new = [np.array([All_positions_velocities[4*i][-1],
                                All_positions_velocities[4*i+1][-1]]) for i in range(self.floe.n)]
        Vel_new = [np.array([All_positions_velocities[4*i+2][-1],
                             All_positions_velocities[4*i+3][-1]]) for i in range(self.floe.n)]
        nodes = []
        Springs1 = set()
        for i, j in Springs_new1:
            i = IndexNewfloe1.index(i)
            j = IndexNewfloe1.index(j)
            Springs1.add((i, j))
        for i in range(len(IndexNewfloe1)):
            nodes.append(
                Node(Points_new[IndexNewfloe1[i]], Vel_new[IndexNewfloe1[i]], IndexNewfloe1[i]))
        New_floe1 = Floe(nodes, Springs1, stiffness=self.floe.k,
                         viscosity=self.floe.mu, id_number=1)

        nodes = []
        Springs2 = set()
        for i, j in Springs_new2:
            i = IndexNewfloe2.index(i)
            j = IndexNewfloe2.index(j)
            Springs2.add((i, j))
        for i in range(len(IndexNewfloe2)):
            nodes.append(
                Node(Points_new[IndexNewfloe2[i]], Vel_new[IndexNewfloe2[i]], IndexNewfloe2[i]))
        New_floe2 = Floe(nodes, Springs2, stiffness=self.floe.k,
                         viscosity=self.floe.mu, id_number=2)

        Solution1 = Percussion_Wall(New_floe1).simulation()
        Solution2 = Percussion_Wall(New_floe2).simulation()

        return All_positions_velocities, New_floe1, New_floe2, Solution1, Solution2, last_step_bef_frac

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
        I = [j for j in range(0, self.floe.n*4, 4)]

        Pos = [np.array([All_positions_velocities[I[i]][time_step],
                         All_positions_velocities[I[i]+1][time_step]]) for i in range(self.floe.n)]

        Vel = [np.array([All_positions_velocities[I[i]+2][time_step],
                         All_positions_velocities[I[i]+3][time_step]]) for i in range(self.floe.n)]

        return Pos, Vel


def Energy_studies(All_positions_velocities, floe):
    """
    Parameters
    ----------
    All_positions_velocities : Evolution of all nodes and velocities
    floe : TYPE
        DESCRIPTION.

    Returns
    -------
    Traction_energy 
    Torsion_energy 
    E_tot

    """
    Length_Mat = floe.length_mat()
    Torsion_Mat = floe.torsion_mat()
    Angle_Mat = floe.angle_init()
    Traction_energy = np.zeros(len(All_positions_velocities[0]))
    Torsion_energy = np.zeros(len(All_positions_velocities[0]))

    for index in range(len(All_positions_velocities[0])):
        Sum1 = 0.
        Sum2 = 0.
        for i, j, k in floe.simplices():
            Qi = np.array([All_positions_velocities[4*i][index],
                          All_positions_velocities[4*i+1][index]])
            Qj = np.array([All_positions_velocities[4*j][index],
                          All_positions_velocities[4*j+1][index]])
            Qk = np.array([All_positions_velocities[4*k][index],
                          All_positions_velocities[4*k+1][index]])

            Sum1 += 0.5 * ((floe.k/sin(Angle_Mat[i, k, j])) * (norm(Qi-Qj) - Length_Mat[i, j])**2
                           + (floe.k/sin(Angle_Mat[i, j, k])) *
                           (norm(Qi-Qk) - Length_Mat[i, k])**2
                           + (floe.k/sin(Angle_Mat[j, i, k])) * (norm(Qj-Qk) - Length_Mat[j, k])**2)

            Sum2 += 0.5 * (Torsion_Mat[i, j, k] * (Angle(Qi, Qj, Qk) - Angle_Mat[i, j, k])**2
                           + Torsion_Mat[i, k, j] *
                           (Angle(Qi, Qk, Qj) - Angle_Mat[i, k, j])**2
                           + Torsion_Mat[j, i, k] * (Angle(Qj, Qi, Qk) - Angle_Mat[j, i, k])**2)

        Traction_energy[index] = Sum1
        Torsion_energy[index] = Sum2

    E_tot = Traction_energy + Torsion_energy
    return Traction_energy, Torsion_energy, E_tot


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
    alpha = 1.  # ductibility coef
    Length_Mat = floe.length_mat()
    Torsion_Mat = floe.torsion_mat()
    Angle_Mat = floe.angle_init()
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

                Sum1 += 0.5 * ((floe.k/sin(Angle_Mat[i, k, j])) * (norm(Qi-Qj) - Length_Mat[i, j])**2
                               + (floe.k/sin(Angle_Mat[i, j, k])) *
                               (norm(Qi-Qk) - Length_Mat[i, k])**2
                               + (floe.k/sin(Angle_Mat[j, i, k])) * (norm(Qj-Qk) - Length_Mat[j, k])**2)

                Sum2 += 0.5 * (Torsion_Mat[i, j, k] * (Angle(Qi, Qj, Qk) - Angle_Mat[i, j, k])**2
                               + Torsion_Mat[i, k, j] *
                               (Angle(Qi, Qk, Qj) - Angle_Mat[i, k, j])**2
                               + Torsion_Mat[j, i, k] * (Angle(Qj, Qi, Qk) - Angle_Mat[j, i, k])**2)

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

    # to complete

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
    else:
        return (vect2-vect1)/norm(vect2-vect1)


def Orthogonal_vect(vect):
    return np.array([-vect[1], vect[0]])


def find_simplice(v1, v2):
    """
    Parameters
    ----------
    v1,v2 : 1st and 2nd springs.
    Returns
    the simplice that its belong to
    """
    l = ((v1[0], v1[1], v2[1]))
    return l


def Angle(A, B, C):
    """
    Input: coordinates of 3 nodes (A,B,C)
    Returns
    Angle ABC by al-kashi's formula 
    """
    # X,Y,Z = A.position(), B.position(), C.position()
    a = norm(C-A)
    b = norm(B-A)
    c = norm(C-B)
    return acos((-a**2+c**2+b**2)/(2*b*c))


def Torsion_Mat(floe: Floe):
    """ stiffness constant of every torsion spring """
    Mat = np.zeros((floe.n, floe.n, floe.n))
    for i, j, k in floe.simplices():
        Mat[i, j, k] = G * floe.length_mat()[i, j] * floe.length_mat()[j,
                                                                       k] / sin(floe.angle_init()[i, j, k])
        Mat[k, j, i] = Mat[i, j, k]

        Mat[i, k, j] = G * floe.length_mat()[i, k] * floe.length_mat()[j,
                                                                       k] / sin(floe.angle_init()[i, k, j])
        Mat[j, k, i] = Mat[i, k, j]

        Mat[j, i, k] = G * floe.length_mat()[i, j] * floe.length_mat()[i,
                                                                       k] / sin(floe.angle_init()[j, i, k])
        Mat[k, i, j] = Mat[j, i, k]
    return Mat


def Angle_Mat(floe: Floe):
    Mat = np.zeros((floe.n, floe.n, floe.n))
    Nodes_positions = floe.get_nodes()
    for i, j, k in floe.simplices():
        Mat[i, j, k] = Angle(Nodes_positions[i],
                             Nodes_positions[j], Nodes_positions[k])
        Mat[k, j, i] = Mat[i, j, k]

        Mat[i, k, j] = Angle(Nodes_positions[i],
                             Nodes_positions[k], Nodes_positions[j])
        Mat[j, k, i] = Mat[i, k, j]

        Mat[j, i, k] = Angle(Nodes_positions[j],
                             Nodes_positions[i], Nodes_positions[k])
        Mat[k, i, j] = Mat[j, i, k]
    return Mat


def Length_Mat(floe: Floe):
    Mat = np.zeros((floe.n, floe.n))
    for (i, j) in floe.springs:
        Mat[i, j] = Spring(floe.nodes[i], floe.nodes[j], None).L0
        Mat[j, i] = Mat[i, j]
    return Mat


""" 
Main system describes all nodes evolutions. 
Each node i depends on TRACTION's spring (i.e spring between node i and node j neighbor)
and TORSION's spring ( formed by the triangle it belongs to)
"""


def System(t, Y, Y0, nb_nodes, Connex_Mat, Length_Mat, m, mu, Traction_mat, Torsion_mat, Angle_init, Triangle_list):
    """
    Parameters
    ----------
    t : time discretisation.
    Y : TYPE
        DESCRIPTION.
    Y0 : init condition of system (init position and velocity).
    nb_nodes : number of nodes.
    Connexe_Mat : connectivity between nodes's matrix.
    Length_Mat : length matrix.
    m : mass of each node.
    mu : viscosity const.
    k : stiffness const of traction spring.
    Torsion_mat: stiffness constant of torsion spring at (i,j,k)

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

    for i in range(0, nb_nodes):
        Y_[2*i] = Q[2*i+1]
        for j in range(i+1, i+nb_nodes):
            j = j % nb_nodes
            u[i, j] = Unit_vect(Q[2*i], Q[2*j])
            Y_[2*i+1] += (1./m) * Connex_Mat[i, j] * (k[i, j] * (norm(Q[2*j]-Q[2*i]) - Length_Mat[i, j]) * u[i, j]
                                                      + mu * (Q[2*j+1] - Q[2*i+1]) @ u[i, j] * u[i, j])

    for i, j, k in Triangle_list:
        u[i, j] = Unit_vect(Q[2*i], Q[2*j])
        u[j, i] = -u[i, j]
        u[i, k] = Unit_vect(Q[2*i], Q[2*k])
        u[k, i] = -u[i, k]
        u[j, k] = Unit_vect(Q[2*j], Q[2*k])
        u[k, j] = -u[j, k]
        Y_[2*i+1] += (1./m) * (G[i, j, k] * (Angle(Q[2*i], Q[2*j], Q[2*k]) - Theta0[i, j, k])/(norm(Q[2*i] - Q[2*j])) * u[i, k]
                               + G[i, k, j] * (Angle(Q[2*i], Q[2*k], Q[2*j]) - Theta0[i, k, j])/norm(Q[2*i] - Q[2*k]) * u[i, j])

        Y_[2*j+1] += (1./m) * (G[j, i, k] * (Angle(Q[2*j], Q[2*i], Q[2*k]) - Theta0[j, i, k])/norm(Q[2*i] - Q[2*j]) * u[j, k]
                               + G[i, k, j] * (Angle(Q[2*i], Q[2*k], Q[2*j]) - Theta0[j, k, i])/norm(Q[2*i] - Q[2*k]) * u[j, i])

        Y_[2*k+1] += (1./m) * (G[j, i, k] * (Angle(Q[2*j], Q[2*i], Q[2*k]) - Theta0[j, i, k])/norm(Q[2*i] - Q[2*k]) * u[k, j]
                               + G[i, j, k] * (Angle(Q[2*i], Q[2*j], Q[2*k]) - Theta0[i, j, k])/norm(Q[2*i] - Q[2*j]) * u[k, i])

    return np.reshape(Y_, (nb_nodes*4))


def System_stable_1(t, Y, Y0, nb_nodes, Connex_Mat, Length_Mat, m, mu, Traction_mat, Torsion_mat, Angle_init, Triangle_list, contact_node):
    """
    Parameters
    ----------
    t : time discretisation.
    Y : TYPE
        DESCRIPTION.
    Y0 : init condition of system (init position and velocity).
    nb_nodes : number of nodes.
    Connexe_Mat : connectivity between nodes's matrix.
    Length_Mat : length matrix.
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
    # find the neigborhood of contact node
    Neighbor_contact = [np.any(Triangle_list[i] == contact_node)
                        for i in range(len(Triangle_list))]
    Neighbor_contact = [i for i, val in enumerate(
        Neighbor_contact) if val == True]
    Neighbor_contact = Triangle_list[Neighbor_contact]

    for i in range(nb_nodes):
        if i == contact_node:
            Y_[2*i] = Q[2*i+1]
            for j in range(i+1, i+nb_nodes):
                j = j % nb_nodes
                u[i, j] = Unit_vect(Q[2*i], Q[2*j])
                Y_[2*i+1] += (1./m) * Connex_Mat[i, j] * (k[i, j] * (norm(Q[2*j]-Q[2*i]) - Length_Mat[i, j]) * u[i, j]
                                                          + mu * (Q[2*j+1] - Q[2*i+1]) @ u[i, j] * u[i, j])

    for i, j, k in Triangle_list:
        u[i, j] = Unit_vect(Q[2*i], Q[2*j])
        u[j, i] = -u[i, j]
        u[i, k] = Unit_vect(Q[2*i], Q[2*k])
        u[k, i] = -u[i, k]
        u[j, k] = Unit_vect(Q[2*j], Q[2*k])
        u[k, j] = -u[j, k]
    # all nodes stable, only the contact node can move
        if i == contact_node or j == contact_node or k == contact_node:
            Y_[2*contact_node+1] += (1./m) * (G[i, j, k] * (Angle(Q[2*i], Q[2*j], Q[2*k]) - Theta0[i, j, k]) /
                                              (norm(Q[2*i] - Q[2*j])) * u[i, k]
                                              + G[i, k, j] * (Angle(Q[2*i], Q[2*k], Q[2*j]) - Theta0[i, k, j]) /
                                              norm(Q[2*i] - Q[2*k]) * u[i, j])

    return np.reshape(Y_, (nb_nodes*4))


def System_stable_neighbor(t, Y, Y0, nb_nodes, Connex_Mat, Length_Mat, m, mu, Traction_mat, Torsion_mat, Angle_init, Triangle_list, contact_node):
    """
    Parameters
    ----------
    t : time discretisation.
    Y : TYPE
        DESCRIPTION.
    Y0 : init condition of system (init position and velocity).
    nb_nodes : number of nodes.
    Connexe_Mat : connectivity between nodes's matrix.
    Length_Mat : length matrix.
    m : mass of each node.
    mu : viscosity const.
    k : stiffness const of traction spring.
    Torsion_mat: stiffness constant of torsion spring at (i,j,k)
    contac_nodes: Nodes in contact with objects

    Returns
    -------
    (evolution of node_i, velo_i
      for i in range nb_nodes as a dynamical system for a neighborhood of contact node) .

    """
    u = np.zeros((nb_nodes, nb_nodes, 2))
    Q = np.reshape(Y, (nb_nodes*2, 2))
    Y_ = np.zeros_like(Q)
    k = Traction_mat
    G = Torsion_mat
    Theta0 = Angle_init
    # find the neigborhood of contact node
    Neighbor_contact = [np.any(Triangle_list[i] == contact_node)
                        for i in range(len(Triangle_list))]
    Neighbor_contact = [i for i, val in enumerate(
        Neighbor_contact) if val == True]
    Neighbor_contact = Triangle_list[Neighbor_contact].tolist()
    list_contact = []
    for i in range(len(Neighbor_contact)):
        list_contact = list_contact + Neighbor_contact[i]

    for i in list_contact:
        # if i == contact_node:
        Y_[2*i] = Q[2*i+1]
        for j in range(i+1, i+nb_nodes):

            j = j % nb_nodes
            u[i, j] = Unit_vect(Q[2*i], Q[2*j])
            Y_[2*i+1] += (1./m) * Connex_Mat[i, j] * (k[i, j] * (norm(Q[2*j]-Q[2*i]) - Length_Mat[i, j]) * u[i, j]
                                                      + mu * (Q[2*j+1] - Q[2*i+1]) @ u[i, j] * u[i, j])

    return np.reshape(Y_, (nb_nodes*4))


dt = 0.00125
G = 1000.  # stiffness of torsion's spring
