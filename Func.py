
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm
from scipy.integrate import solve_ivp
from scipy.spatial import Voronoi, voronoi_plot_2d, Delaunay
from graph import *
from math import acos, sin
from collections import Counter
import matplotlib.animation as animation

# 4 classes Node-> Spring-> Ice-Floe-> Percussion to generate the Percussion of 2 floes


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

    def border_vertices_index(self):

        All_vertices = []
        for triangles in self.simplices():
            for i in range(len(triangles)):
                for j in range(i+1, len(triangles)): 
                    All_vertices.append((min(triangles[i], triangles[j]),max(triangles[i], triangles[j])))
                    
        S = set(All_vertices)
        Border = [e for e in S if All_vertices.count(e) == 1]
        
        return Border

    def border_nodes_index(self):
        BV = self.border_vertices_index()
        l = []
        for i,j in BV: 
            l.append(i)
            l.append(j)
        return set(l)

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

    def E0_traction(self):
        Sum = 0
        L = self.length_mat()
        for i, j, k in self.simplices():
            Sum += self.k * (L[i, j]**2 + L[i, k]**2 + L[j, k]**2)
        return Sum/2.

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
        return tri.simplices

    def get_velocity(self):
        return np.array([node.velocity() for node in self.nodes])

    def get_springs(self):  # not necessary yet ?
        return [(self.get_nodes()[i], self.get_nodes()[j]) for (i, j) in self.springs]

    def connexe_mat(self):
        Mat = np.zeros((self.n, self.n))
        for (i, j) in self.springs:
            Mat[i, j] = 1
            Mat[j, i] = Mat[i, j]
        return Mat

    def length_mat(self):
        Mat = np.zeros((self.n, self.n))
        for (i, j) in self.springs:
            Mat[i, j] = Spring(self.nodes[i], self.nodes[j], None).L0
            Mat[j, i] = Mat[i, j]
        return Mat

    def angle_init(self):
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
        Mat = np.zeros((self.n, self.n))
        for i, j, k in self.simplices():
            Mat[i, j] += self.k / sin(self.angle_init()[i, k, j])
            Mat[j, i] = Mat[i, j]

            Mat[i, k] += self.k / sin(self.angle_init()[i, j, k])
            Mat[k, i] = Mat[i, k]

            Mat[j, k] += self.k / sin(self.angle_init()[k, i, j])
            Mat[k, j] = Mat[j, k]

        return Mat

    def Move(self, time_end: float, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat):
        N = 800
        t = np.linspace(0, time_end, N)
        All_pos = self.get_nodes()
        All_vel = self.get_velocity()
        Y0_ = np.array([])
        for i in range(self.n):
            Y0_ = np.append(Y0_, All_pos[i])
            Y0_ = np.append(Y0_, All_vel[i])

        Sol = solve_ivp(System, [0, time_end], Y0_, t_eval=t,
                        args=(Y0_, self.n, self.connexe_mat(), Length_Mat, self.m,
                              self.mu, Traction_Mat, Torsion_Mat, Angle_Mat, self.simplices()))
        return Sol

    def plot_init(self):
        plt.figure()
        for (i, j) in self.springs:
            plt.plot([self.nodes[i].position()[0], self.nodes[j].position()[0]],
                     [self.nodes[i].position()[1], self.nodes[j].position()[1]])
            plt.text(self.nodes[i].position()[0],
                     self.nodes[i].position()[1], self.nodes[i].id)
            plt.text(self.nodes[j].position()[0],
                     self.nodes[j].position()[1], self.nodes[j].id)
        plt.plot(self.center()[0], self.center()[1], 'o', color='red')
        theta = np.linspace(0, 2*np.pi, 100)
        x = self.First_radius()*np.cos(theta) + self.center()[0]
        y = self.First_radius()*np.sin(theta) + self.center()[1]
        plt.plot(x, y, "--")

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

                # change self.k by self.traction_mat()[i]
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


class Percussion:
    def __init__(self, floe1: Floe, floe2: Floe, restitution_coef=0.4, time_end=4., eta=0.1):
        self.t = np.linspace(0, time_end, 800)
        self.floe1 = floe1
        self.floe2 = floe2
        self.eps = restitution_coef
        self.eta = eta

    def compute_before_contact(self):
        r1 = self.floe1.First_radius()
        r2 = self.floe1.First_radius()
        contact_distance = r1 + r2 + self.eta
        center1_position = [self.floe1.center(
        ) + self.t[i]*self.floe1.center_velocity() for i in range(self.t.size)]
        center2_position = [self.floe2.center(
        ) + self.t[i]*self.floe2.center_velocity() for i in range(self.t.size)]
        distance_evolution = [
            norm(center1_position[i]-center2_position[i]) for i in range(self.t.size)]
        return distance_evolution < contact_distance


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


dt = 0.005
G = 100.  # stiffness of torsion's spring
