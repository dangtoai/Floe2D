# -*- coding: utf-8 -*-
from Func import *
from graph import *
import numpy as np

if __name__ == '__main__':
    ######################
    ##### An ice floe ####
    ######################
    
    Points = np.array([[1, 1.], [0., -1.], [-0.5, 0.]])
    # Points = np.array([[1, 0.], [0., -1.], [-0., 1.]])
    V0     = np.array([1., 0.])
    # V1     = np.array([-1.0, 0.])
    Nodes  = []
    for i in range(len(Points)):
        Nodes.append(Node(Points[i], V0, i))

    Springs = {(0,1),(1,2),(0,2)}
    k = 1.
    
    floe = Floe(nodes=Nodes, springs=Springs, stiffness=k, viscosity=k/10.,  id_number = 1 )
    
    Problem = Percussion_Wall(floe, Wall = 5., restitution_coef=0.9, time_end = 4., eta=0.01)
    

    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=True, xlim=(-.1, 5.5), ylim=(-1.9, 1.5))
    ax.set_aspect('equal')
    ax.grid()
    plt.axvline(x=5., color = "red")
    line1, = ax.plot([], [], '.-', lw= .95)
    time_template = 'time = % 10fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    Route = floe.Route()
    
    All_positions_velocities = Problem.simulation()
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
            thisx = np.append(thisx, All_positions_velocities[j][i])
        for j in Iy:
            thisy = np.append(thisy, All_positions_velocities[j][i])
        for k in Route:
            thisx = np.append(thisx, thisx[k])
            thisy = np.append(thisy, thisy[k])
        line1.set_data(thisx[floe.n:], thisy[floe.n:])
        time_text.set_text(time_template % (i*dt))
        return line1, time_text

    ani = animation.FuncAnimation(fig, animate_spring,
                                  np.arange(0, len(All_positions_velocities[0])-200), interval=25, blit=False)

    
    