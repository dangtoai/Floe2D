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
                       [0.45683322, 0.4151012],
                       [0.64914405, 0.28352508],
                       [0.27848728, 0.69313792],
                       [0.6762549, 0.44045372],
                       [0.59086282, 0.15686774],
                       [0.02398188, 0.54464902]])

    Nodes = []
    
    t_end = 4.
    # dt = 0.01
    V0 = np.array([0.5, 0.])
    # V1 = np.array([0.75, 0.])

    for i in range(len(Points)):
        Nodes.append(Node(Points[i], V0, i))
    # Nodes[1] = Node(Points[1], V1, 1)
    # Nodes[5] = Node(Points[5], -V1, 5)
    
    Springs = {(0, 7), (1, 2), (0, 4), (2, 7), (2, 3), (1, 7), (0, 2), (2, 6),
               (4, 5), (0, 5), (3, 6), (1, 6), (2, 5), (4, 7), (3, 5)}
    
    k = 1000
    New_floe = Floe(nodes=Nodes, springs=Springs, stiffness=k)
    
    #second floe
    Points = np.array([[3.52981736, 0.73588211],
           [3.41880743, 0.51803641],
           [3.33540785, 0.5788586 ],
           [3.62251943, 0.6453551 ],
           [3.43814143, 0.99022427]])
    Nodes = []
    V0_ = -V0
    for i in range(len(Points)):
        Nodes.append(Node(Points[i], V0_, i))
    Springs = {(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (2, 4), (3, 4)}
    Second_floe = Floe(nodes=Nodes, springs=Springs, stiffness=k/2.)
    
    # print("connexe matrix = \n", New_floe.connexe_mat(), "\n")
    # print("length matrix = \n",New_floe.length_mat(), "\n")

    # print(" Initial position ", New_floe.get_nodes())
    # print(" Initial velocity ",New_floe.get_velocity())

    New_floe.plot_init()
    Second_floe.plot_init()
    # New_floe.plot_displacements(t_end)
    # Route = New_floe.Route()
    Sol = New_floe.Move(t_end)
    Sol2= Second_floe.Move(t_end)
    
    ### re-implemented in Func.py!!! ###
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-.5, 4.), ylim=(-.5, 1.5))
    ax.set_aspect('equal')
    ax.grid()
    line1, = ax.plot([], [], 'o-', lw=1.)
    line2, = ax.plot([], [], 'o-', lw=1.)
    time_template = 'time = %.9fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    Route_ = Second_floe.Route()
    Route = New_floe.Route()
    
    def init():
        line1.set_data([], [])
        # time_text.set_text('')
        return line1, 
    
    def animate_spring(i):
        Ix = [j for j in range(0, New_floe.n*4, 4)]
        Iy = [j for j in range(1, New_floe.n*4, 4)]
        thisx = []
        thisy = []
        for j in Ix:
            thisx = np.append(thisx, Sol.y[j][i])
        for j in Iy:
            thisy = np.append(thisy, Sol.y[j][i])
        for k in Route:
            thisx = np.append(thisx,thisx[k])
            thisy = np.append(thisy,thisy[k])
            
        Ix_ = [j for j in range(0, Second_floe.n*4, 4)]
        Iy_ = [j for j in range(1, Second_floe.n*4, 4)]
        thisx_ = []
        thisy_ = []
        for j in Ix_:
            thisx_ = np.append(thisx_, Sol2.y[j][i])
        for j in Iy_:
            thisy_ = np.append(thisy_, Sol2.y[j][i])
        for k in Route_:
            thisx_ = np.append(thisx_,thisx_[k])
            thisy_ = np.append(thisy_,thisy_[k])
        line2.set_data(thisx_[Second_floe.n:], thisy_[Second_floe.n:])
        line1.set_data(thisx[New_floe.n:], thisy[New_floe.n:])
        
        return line1, line2, time_text.set_text(time_template % (i*dt))

    ani1 = animation.FuncAnimation(fig, animate_spring, 
                                    np.arange(0, len(Sol.y[0])), interval=25, blit=False)

    
