#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 11:07:24 2022

@author: phandangtoai
"""
from Func import *
from graph import *
import numpy as np
import matplotlib.animation as animation
from matplotlib.animation import FFMpegFileWriter
from IPython import display
from datetime import datetime

if __name__ == '__main__':
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

    # Points = np.array([[0.44080984, 0.55885409],
    #         [0.02987621, 0.25925245],
    #         [0.45683322, 0.4151012 ],
    #         [0.64914405, 0.28352508],
    #         [0.27848728, 0.69313792],
    #         [0.6762549 , 0.44045372],
    #         [0.59086282, 0.15686774],
    #         [0.02398188, 0.54464902]])
    
    Nodes = []
    V0 = np.array([0., 0.])
    V1 = np.array([-35.5, 0.])
    for i in range(len(Points)):
        Nodes.append(Node(Points[i], V0, i))
    
    Nodes[5] = Node(Points[5], V1, 5)
    Springs = {(0, 1), (2, 4), (1, 2), (0, 4), (1, 5), (4, 6), (0, 3), (0, 6), (4, 5), (0, 2), (3, 6), (2, 5), (1, 3)}
    # Springs = {(0, 7), (1, 2), (0, 4), (2, 7), (2, 3), (1, 7), (0, 2), (2, 6), (4, 5), (0, 5), (3, 6), (1, 6), (2, 5), (4, 7), (3, 5)}
    # Springs = {(0, 1), (1, 2), (0, 3), (2, 3), (0, 2), (1, 3)}
    k = 1000.
    floe = Floe(nodes=Nodes, springs=Springs,
                stiffness=k, viscosity=k/5., id_number=1)
    
    Traction_Mat = floe.traction_mat() 
    Length_Mat   = floe.length_mat()
    Torsion_Mat  = floe.torsion_mat()
    Angle_Mat    = floe.angle_init()

    All_positions_velocities = floe.Move(1., Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat).y
    
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=True,
                          xlim=(-1.5, 1.1), ylim=(-.75, 1.51))
    ax.set_aspect('equal')
    ax.grid()
    # plt.axvline(x=2., color="red")
    line1, = ax.plot([], [], '.-', lw=1.95)
    line2, = ax.plot([], [], '.-', lw=1.95)
    time_template = 'time = % 10fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    Route = floe.Route()
    def init():
        line1.set_data([], [])
        time_text.set_text(' ')
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
            
        line1.set_data(thisx[floe.n: ], 
                       thisy[floe.n: ])
        
        
        time_text.set_text(time_template % (i*dt))
        return line1, time_text

    ani = animation.FuncAnimation(fig, animate_spring,
                                  np.arange(0, len(All_positions_velocities[0])), interval=2, blit=False)















