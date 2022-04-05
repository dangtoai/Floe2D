# -*- coding: utf-8 -*-
from Func import *
# from graph import *
import numpy as np
from matplotlib.animation import PillowWriter
# from math import acos, degrees

if __name__ == '__main__':
    ######################
    ##### An ice floe ####
    ######################
    
    #global constant:   
 
    Points = np.array([[0, 0.], [0., 1.], [1., 0.]])
    V0     = np.array([0., 0.])
    V1     = np.array([0.75, 0.0])
    Nodes  = []
    for i in range(len(Points)):
        Nodes.append(Node(Points[i], V0, i))
    Nodes[1] = Node(Points[1], V1, 1)
    Springs = {(0,1),(0,2)}
    k = 1000.
    floe = Floe(nodes=Nodes, springs=Springs, stiffness=k, viscosity=k/10.,  id_number = 1 )
    
    t_end = 4.
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=True, xlim=(-.1, 2.5), ylim=(-1.5, 1.5))
    ax.set_aspect('equal')
    ax.grid()
    # plt.axvline(x=5., color = "red")
    line1, = ax.plot([], [], '.-', lw= .95)
    time_template = 'time = % 10fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    Route = floe.Route()

    #compute after contact
    Sol = floe.Move(t_end)
    # floe.plot_displacements(t_end)
    Node2x = Sol.y[4]
    Node2y = Sol.y[5]
    Node1x = Sol.y[0]
    Node1y = Sol.y[1]
    Node3x = Sol.y[8]
    Node3y = Sol.y[9]
    
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
        line1.set_data(thisx[-floe.n:], thisy[-floe.n:])
        time_text.set_text(time_template % (i*dt))
        # return thisx, thisy
        return line1, time_text
    ani = animation.FuncAnimation(fig, animate_spring, 
                                    np.arange(0,len(Node1x)), interval=25, blit=False)
    
    # ani.save("0torsion.gif", writer=PillowWriter(fps=25))

