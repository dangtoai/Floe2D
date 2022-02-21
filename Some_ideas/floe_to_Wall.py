# -*- coding: utf-8 -*-

from Func import *
from graph import *
import numpy as np

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
    V0 = np.array([0.25, 0.])
    # V1 = np.array([0.5, 0.])

    for i in range(len(Points)):
        Nodes.append(Node(Points[i], V0, i))
    # Nodes[1] = Node(Points[1], V1, 1)
    # Nodes[5] = Node(Points[5], -V1, 5)
    
    Springs = {(0, 7), (1, 2), (0, 4), (2, 7), (2, 3), (1, 7), (0, 2), (2, 6),
               (4, 5), (0, 5), (3, 6), (1, 6), (2, 5), (4, 7), (3, 5)}
    k = 100.
    First_floe = Floe(nodes=Nodes, springs=Springs, stiffness=k, viscosity=k/10.  ,id_number = 1 )

    t_end = 6.
    collision_dist   = 0.01
    coef_restitution = 0.99
    
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=True, xlim=(-2.1, 8.5), ylim=(-1.5, 1.5))
    ax.set_aspect('equal')
    ax.grid()
    plt.axvline(x=5., color = "red")
    line1, = ax.plot([], [], '.-', lw= .95)
    time_template = 'time = % 10fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    Route = floe.Route()
    
    Sol = floe.Move(t_end)
    
    Positions = Sol.y
    
    # tant que: la position en x de n'importe quel noeud depasse Wall:
    #   Faire:  tronquer toutes les vitesses et positions jusqu'a cette indice
    #           update velocity of collision'nodes
    #           lancer la simulation et concatener avec les positions et noeuds.
    
    
    
    
    
    def init():
        line1.set_data([], [])
        time_text.set_text('')
        return line1, time_text
    
    def animate_spring(i):
        Ix = [j for j in range(0, floe.n*4, 4)]
        Iy = [j for j in range(1, floe.n*4, 4)]
        thisx = []
        thisy = []
        for j in Ix:
            thisx = np.append(thisx, [Node1x[i], Node2x[i], Node3x[i]])
        for j in Iy:
            thisy = np.append(thisy, [Node1y[i], Node2y[i], Node3y[i]])
        for k in Route:
            thisx = np.append(thisx,thisx[k])
            thisy = np.append(thisy,thisy[k])
        line1.set_data(thisx[floe.n:], thisy[floe.n:])
        time_text.set_text(time_template % (i*dt))
        # return thisx, thisy
        return line1, time_text
    
    ani = animation.FuncAnimation(fig, animate_spring, 
                                    np.arange(400,len(Node1x)), interval=25, blit=False)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

