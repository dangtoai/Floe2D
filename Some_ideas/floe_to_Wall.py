# -*- coding: utf-8 -*-
from Func import *
from graph import *
import numpy as np
import matplotlib.animation as animation
from matplotlib.animation import FFMpegFileWriter
from IPython import display
from datetime import datetime

if __name__ == '__main__':
    # start = datetime.now()
    ######################
    ##### An ice floe ####
    ######################

    # first floe
    Points = np.array([[0.25298236, 0.59733394],
      [0.43479153, 0.0089861 ],
      [0.77938292, 0.38657128],
      [0.19768507, 0.04416006],
      [0.86299324, 0.95665297],
      [0.98340068, 0.43614665],
      [0.16384224, 0.94897731]])

    # Points = np.array([[0.33033482, 0.26682728],
    #                     [0.20464863, 0.62113383],
    #                     [0.61927097, 0.52914209],
    #                     [0.29965467, 0.13457995]])

    Nodes = []
    V0 = np.array([0.85, 0.])

    for i in range(len(Points)):
        Nodes.append(Node(Points[i], V0, i))

    Springs = {(0, 1), (2, 4), (1, 2), (0, 4), (1, 5), (4, 6), (0, 3), (0, 6), (4, 5), (0, 2), (3, 6), (2, 5), (1, 3)}
    # Springs = {(0, 7), (1, 2), (0, 4), (2, 7), (2, 3), (1, 7), (0, 2), (2, 6), (4, 5), (0, 5), (3, 6), (1, 6), (2, 5), (4, 7), (3, 5)}
    # Springs = {(0, 1), (1, 2), (0, 3), (2, 3), (0, 2), (1, 3)}
    k = 100.
    floe = Floe(nodes=Nodes, springs=Springs,
                stiffness=k, viscosity=k/5., id_number=1)

    Problem = Percussion_Wall(floe, Wall = 2.) 
    # F1,F2 = Problem.simulation_with_fracture()
    # print(S,F1,F2) 
    
    # All_positions_velocities = Problem.simulation()
    All_positions_velocities, F1, F2, AP1, AP2, last_step = Problem.simulation_with_fracture()
    
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=True,
                          xlim=(.5, 2.1), ylim=(-1.1, 1.11))
    ax.set_aspect('equal')
    ax.grid()
    plt.axvline(x=2., color="red")
    line1, = ax.plot([], [], '.-', lw=1.95)
    line2, = ax.plot([], [], '.-', lw=1.95)
    time_template = 'time = % 10fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    Route = floe.Route()
    Route1 = F1.Route()
    Route2 = F2.Route()
    
    # animation
    def init():
        line1.set_data([], [])
        time_text.set_text('')
        return line1, time_text

    def animate_spring(i):
        if i <= last_step:
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
            line2.set_data(None,None)
            time_text.set_text(time_template % (i*dt))
        
        if i<= 1024 and i> last_step: 
            Ix = [j for j in range(0, F1.n*4, 4)]
            Iy = [j for j in range(1, F1.n*4, 4)]
            thisx = []
            thisy = []
            for j in Ix:
                thisx = np.append(thisx, AP1[j][i-last_step])
            for j in Iy:
                thisy = np.append(thisy, AP1[j][i-last_step])
            for k in Route1:
                thisx = np.append(thisx, thisx[k])
                thisy = np.append(thisy, thisy[k])
            line1.set_data(thisx[F1.n:], thisy[F1.n:])
            
            Ix = [j for j in range(0, F2.n*4, 4)]
            Iy = [j for j in range(1, F2.n*4, 4)]
            thisx = []
            thisy = []
            for j in Ix:
                thisx = np.append(thisx, AP2[j][i-last_step])
            for j in Iy:
                thisy = np.append(thisy, AP2[j][i-last_step])
            for k in Route2:
                thisx = np.append(thisx, thisx[k])
                thisy = np.append(thisy, thisy[k])
            line2.set_data(thisx[F2.n:], thisy[F2.n:])
            time_text.set_text(time_template % (i*dt))
            return line1, line2, time_text
        

    ani = animation.FuncAnimation(fig, animate_spring,
                                  np.arange(0, len(AP2[0])+len(All_positions_velocities[0])-200), interval=25, blit=False)
    
