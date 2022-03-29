# -*- coding: utf-8 -*-
from Func import *
# from graph import *
import numpy as np
# from scipy.spatial import Voronoi, voronoi_plot_2d, Delaunay 

if __name__ == '__main__':
    ######################
    ##### An ice floe ####
    ######################
    
    #first floe
    Points = np.array([[0.44080984, 0.55885409],
                       [0.02987621, 0.25925245],
                       [0.45683322, 0.4151012 ],
                       [0.64914405, 0.28352508],
                       [0.27848728, 0.69313792],
                       [0.6762549, 0.44045372 ],
                       [0.59086282, 0.15686774],
                       [0.02398188, 0.54464902]])

    Nodes = []
    t_end = 4.
    V0 = np.array([1, 0.])
    # V1 = np.array([0.5, 0.])

    for i in range(len(Points)):
        Nodes.append(Node(Points[i], V0, i))
    # Nodes[1] = Node(Points[1], V1, 1)
    # Nodes[5] = Node(Points[5], -V1, 5)
    
    Springs = {(0, 7), (1, 2), (0, 4), (2, 7), (2, 3), (1, 7), (0, 2), (2, 6),
               (4, 5), (0, 5), (3, 6), (1, 6), (2, 5), (4, 7), (3, 5)}
    k = 100.
    First_floe = Floe(nodes=Nodes, springs=Springs, stiffness=k, viscosity=k/10.  ,id_number = 1 )

    
    # First_floe.plot_init()
    # print("all triangle in system", First_floe.simplices())
    
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=True, xlim=(-.1, 5.), ylim=(-.0, 1.1))
    ax.set_aspect('equal')
    ax.grid()
    line1, = ax.plot([], [], '.-', lw= .5)
    time_template = 'time = % 10fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    Route = First_floe.Route()
    
    Sol = First_floe.Move(t_end)
    
    ### Perturbation avec une masse ponctuelle
    collision_floe = First_floe.evolution(t_end, t_end)
    # print(collision_floe.get_velocity())
    collision_floe.update_velocity_node(5, np.array([-.125,  0]))
    collision_floe.update_velocity_node(3, np.array([-.125,  0]))
    collision_floe.update_velocity_node(6, np.array([-.125,  0]))
    # print(collision_floe.get_velocity())
    collision_floe.plot_init()
    # collision_floe.plot_displacements(4)
    Sol2 = collision_floe.Move(t_end)
    
    def init():
        line1.set_data([], [])
        time_text.set_text('')
        return line1, time_text
    
    def animate_spring(i):
        Ix = [j for j in range(0, First_floe.n*4, 4)]
        Iy = [j for j in range(1, First_floe.n*4, 4)]
        thisx = []
        thisy = []
        for j in Ix:
            thisx = np.append(thisx, np.append(Sol.y[j], Sol2.y[j])[i])
        for j in Iy:
            thisy = np.append(thisy, np.append(Sol.y[j], Sol2.y[j])[i])
        for k in Route:
            thisx = np.append(thisx,thisx[k])
            thisy = np.append(thisy,thisy[k])
        line1.set_data(thisx[First_floe.n:], thisy[First_floe.n:])
        time_text.set_text(time_template % (i*dt))
        return line1, time_text
    
    ani = animation.FuncAnimation(fig, animate_spring, 
                                    np.arange(600, 1000), interval=25, blit=False)

            