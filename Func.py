
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm
from scipy.integrate import solve_ivp
import matplotlib.animation as animation
from graph import *

# 4 classes Node-> Spring-> Ice-Floe-> Percussion to generate the Percussion of 2 floes


class Node:
    """A class representing one node of an ice floe"""

    def __init__(self, position, velocity: np.array([0,0]), id_number = None):
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
    def __init__(self, nodes=None, springs=None, mass=1.0, stiffness=15, viscosity=2.0, tenacity=1.0, id_number=None):
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

        self.m = mass
        self.k = stiffness
        self.mu = viscosity
        self.L = tenacity
        # self.v0 = rigid_velocity  # One velocity for all nodes
        self.id = id_number
    
    def center(self):
        return 0
    
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
        # return np.array([node.position() for node in self.nodes])
        return [node.position() for node in self.nodes]
    
    def get_velocity(self):
        return np.array([node.velocity() for node in self.nodes])
    
    def get_springs(self):          #not necessary yet ?
        return [(self.get_nodes()[i], self.get_nodes()[j]) for (i,j) in self.springs]
        
    def connexe_mat(self):
        Mat = np.zeros((self.n, self.n))
        for (i,j) in self.springs:
            Mat[i,j] = 1
            Mat[j,i] = Mat[i,j]
        return Mat
    
    def length_mat(self):
        # v = np.array([0, 0])
        Mat = np.zeros((self.n, self.n))
        for (i,j) in self.springs:
            Mat[i,j] = Spring(self.nodes[i], self.nodes[j], None).L0
            Mat[j,i] = Mat[i,j]
        return Mat

    def Move(self, time_end: float()):
        N = 800
        t = np.linspace(0, time_end, N)
        All_pos = self.get_nodes()
        All_vel = self.get_velocity()
        Y0_ = np.array([])
        for i in range(self.n):
            Y0_ = np.append(Y0_, All_pos[i])
            Y0_ = np.append(Y0_, All_vel[i])
        
        Sol = solve_ivp(System, [0, time_end], Y0_, t_eval=t, 
                        args=( Y0_, self.n, self.connexe_mat(), self.length_mat(), self.m, self.mu, self.k ))
        return Sol
    
    def plot_init(self):
        plt.figure()
        for (i,j) in self.springs:
            plt.plot([self.nodes[i].position()[0], self.nodes[j].position()[0]],
                     [self.nodes[i].position()[1], self.nodes[j].position()[1]])
            plt.text(self.nodes[i].position()[0], self.nodes[i].position()[1], self.nodes[i].id)
            plt.text(self.nodes[j].position()[0], self.nodes[j].position()[1], self.nodes[j].id)
    
    def plot_displacements(self, time_end):
        """
        Parameters
        ----------
        time_end : .

        Plot the position's evolution of all nodes.

        """
        time = self.Move(time_end).t
        solution = self.Move(time_end).y
        Index_x = np.arange(0, 4*self.n, 4)
        Index_y = Index_x + 1
        plt.figure()
        for i in Index_x:
            plt.plot(time, solution[i])
        for i in Index_y:
            plt.plot(time, solution[i])
        plt.xlabel("time(s)")
        plt.legend()
    
    def evolution(self, time_end, time):
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
        Res = self.Move(time_end).y
        
        New_Nodes_positions = []
        New_Nodes_velocity = []
        Nodes = []
        for i in range(0, self.n*4, 4):
            # print(i)
            New_Nodes_positions.append(np.array([Res[i ][time], Res[i+1][time]]))
            New_Nodes_velocity.append(np.array([Res[i+2][time], Res[i+3][time]]))
        
        for i in range(self.n):
            Nodes.append(Node(New_Nodes_positions[i], New_Nodes_velocity[i], i))
        
        New_floe = Floe(nodes = Nodes, springs=self.springs )
        return New_floe
    
    def position_at_time(self, time_end, time):
        """
        Parameters
        ----------
        time_end.
        
        At each time's step, save the position of all node as a floe.
        
        Returns
        -------
        Position of all nodes at time t .

        """
        assert 0 <= time <= time_end, "time must be inside simulation interval"
        return self.New_floe(time_end, time).get_nodes()
    
    def velocity_at_time(self, time_end, time):
        """
        Parameters
        ----------
        time_end.
        
        At each time's step, save the position of all node as a floe.
        
        Returns
        -------
        velocity of all nodes at time t .

        """
        assert 0 <= time <= time_end, "time must be inside the simulation interval"
        return self.New_floe(time_end, time).get_velocity()
    
class Percussion:
    def __init__(self, floe1:Floe, floe2:Floe, restitution_coef=0.4, time_end = 4., eta = 0.0001):
        self.t = np.linspace(0, time_end, 1000)
        self.floe1 = floe1
        self.floe2 = floe2
        self.eps = restitution_coef
    
    def check_collision(self):
        collide = False
        for i in range(len(self.t)):
            0
        return collide
        

def node_to_floe(node: Node , floe: Floe):
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
    dist = [norm(node.position() - floe.nodes[i].position()) for i in range(floe.n)]
    return min(dist)

def Unit_vect(vect1, vect2):
    if (vect1[0] == vect2[0] and vect1[1] == vect2[1]):
        return 0.
    else:
        return (vect2-vect1)/norm(vect2-vect1)


def System(t, Y, Y0, nb_nodes, Connex_Mat, Length_Mat, m, mu, k):
    """
    
    Parameters
    ----------
    t : time discretisation.
    Y : TYPE
        DESCRIPTION.
    Y0 : init condition.
    nb_nodes : number of nodes.
    Connexe_Mat : contact matrix.
    Length_Mat : length matrix.
    m : masse of each node.
    mu : viscosity const.
    k : stiffness const.

    Returns
    -------
    (evolution of node_i, velo_i
     for i in range nb_nodes as a dynamical system).
    
    """
    u = np.zeros((nb_nodes, nb_nodes, 2))
    Q = np.reshape(Y, (nb_nodes*2, 2))
    Y_ = np.zeros_like(Q)

    # if node0 is stable, its velocity and acceleration = 0
    # begin at 1
    for i in range(0, nb_nodes):
        Y_[2*i] = Q[2*i+1]
        for j in range(i+1, i+nb_nodes):  # for j in range(i+1, i+nb_nodes)
            j = j % nb_nodes
            u[i, j] = Unit_vect(Q[2*i], Q[2*j])
            Y_[2*i+1] += (1./m)*Connex_Mat[i, j]*(k*(norm(Q[2*j]-Q[2*i]) - Length_Mat[i, j])*u[i, j]
                                                   + mu*(Q[2*j+1] - Q[2*i+1])@u[i, j]*u[i, j])
    return np.reshape(Y_, (nb_nodes*4))

dt = 0.005
