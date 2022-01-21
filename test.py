# -*- coding: utf-8 -*-
from Func import *
from graph import *

if __name__ == '__main__':
    ######################
    ##### An ice floe ####
    ######################
    Points = np.array([[0.44080984, 0.55885409],
                       [0.02987621, 0.25925245],
                       [0.45683322, 0.4151012],
                       [0.64914405, 0.28352508],
                       [0.27848728, 0.69313792],
                       [0.6762549, 0.44045372],
                       [0.59086282, 0.15686774],
                       [0.02398188, 0.54464902]])

    Nodes = []
    
    t_end = 5
    # dt = 0.01
    V0 = np.array([0., 0.])
    V1 = np.array([0.5, 0.5])

    for i in range(len(Points)):
        Nodes.append(Node(Points[i], V0, i))
    Nodes[1] = Node(Points[1], V1, 1)
    
    Springs = {(0, 7), (1, 2), (0, 4), (2, 7), (2, 3), (1, 7), (0, 2), (2, 6),
               (4, 5), (0, 5), (3, 6), (1, 6), (2, 5), (4, 7), (3, 5)}
    
    k = 1000
    New_floe = Floe(nodes=Nodes, springs=Springs, stiffness=k)
    

    # print("connexe matrix = \n", New_floe.connexe_mat(), "\n")
    # print("length matrix = \n",New_floe.length_mat(), "\n")

    print(" Initial position ", New_floe.get_nodes())
    print(" Initial velocity ",New_floe.get_velocity())

    New_floe.plot_init()
    Route = New_floe.Route()
    
    Sol = New_floe.Move(t_end)
    # New_floe.animation_move(t_end)
    
    ### re-implemented in Func.py!!! ###
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-.5, 1.5), ylim=(-.5, 1.5))
    ax.set_aspect('equal')
    ax.grid()
    line, = ax.plot([], [], 'o-', lw=2)
    center,= ax.plot([], [], 'o', color='r')
    time_template = 'time = %.9fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    def animate_spring(i):
        Ix = [i for i in range(0, New_floe.n*4, 4)]
        Iy = [i for i in range(1, New_floe.n*4, 4)]
        thisx = []
        thisy = []
        for j in Ix:
            thisx = np.append(thisx, Sol.y[j][i])
        for j in Iy:
            thisy = np.append(thisy, Sol.y[j][i])
        for i in Route:
            thisx = np.append(thisx,thisx[i])
            thisy = np.append(thisy,thisy[i])
        line.set_data(thisx[New_floe.n:], thisy[New_floe.n:])
        # time_text.set_text(time_template % (i*dt))
        return line, 
    ani = animation.FuncAnimation(fig, animate_spring, np.arange(1, len(Sol.y[0])),
                                interval=2.5, blit=False)
