#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 16:36:29 2022

@author: phandangtoai
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 11:07:24 2022

@author: phandangtoai
"""
# from matplotlib.animation import FFMpegFileWriter
# from IPython import display
# from datetime import datetime

from Func import *
from graph import *
import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib import transforms
if __name__ == '__main__':
    ######################
    ##### An ice floe ####
    ######################

    # first floe
    Points = np.array([[0.29090474, 0.78031476],
           [0.51082761, 0.30636353],
           [0.89294695, 0.22195788],
           [0.89629309, 0.38797126],
           [0.12558531, 0.93638365],
           [0.20724288, 0.97599542],
           [0.0514672 , 0.67238368],
           [0.44080984, 0.90283411],
           [0.02987621, 0.84575087],
           [0.45683322, 0.37799404],
           [0.64914405, 0.09221701],
           [0.27848728, 0.6534109 ],
           [0.6762549 , 0.55784076],
           [0.59086282, 0.36156476],
           [0.02398188, 0.2250545 ],
           [0.55885409, 0.40651992],
           [0.25925245, 0.46894025],
           [0.4151012 , 0.26923558],
           [0.28352508, 0.29179277],
           [0.69313792, 0.4576864 ],
           [0.44045372, 0.86053391],
           [0.15686774, 0.5862529 ],
           [0.54464902, 0.28348786]])
    
    # Points = np.array([[0.76590786, 0.1441643 ],
    #         [0.51841799, 0.16561286],
    #         [0.2968005 , 0.96393053],
    #         [0.18772123, 0.96022672],
    #         [0.08074127, 0.18841466],
    #         [0.7384403 , 0.02430656],
    #         [0.44130922, 0.20455555],
    #         [0.15830987, 0.69984361],
    #         [0.87993703, 0.77951459],
    #         [0.27408646, 0.02293309],
    #         [0.41423502, 0.57766286],
    #         [0.29607993, 0.00164217],
    #         [0.62878791, 0.51547261],
    #         [0.57983781, 0.63979518],
    #         [0.5999292 , 0.9856244 ],
    #         [0.26581912, 0.2590976 ],
    #         [0.28468588, 0.80249689],
    #         [0.25358821, 0.87048309],
    #         [0.32756395, 0.92274961]])

    # Points = np.array([[0.97627445, 0.8662893],
    #                     [0.00623026, 0.17316542],
    #                     [0.25298236, 0.07494859],
    #                     [0.43479153, 0.60074272],
    #                     [0.77938292, 0.16797218],
    #                     [0.19768507, 0.73338017],
    #                     [0.86299324, 0.40844386],
    #                     [0.98340068, 0.52790882],
    #                     [0.16384224, 0.93757158],
    #                     [0.59733394, 0.52169612],
    #                     [0.0089861, 0.10819338],
    #                     [0.38657128, 0.15822341],
    #                     [0.04416006, 0.54520265],
    #                     [0.95665297, 0.52440408],
    #                     [0.43614665, 0.63761024],
    #                     [0.94897731, 0.40149544],
    #                     [0.78630599, 0.64980511]])

    Nodes = []
    V0 = np.array([0., 0.])
    V1 = np.array([5.5, .5])
    # V1 = np.array([-30.5, 0.])
    for i in range(len(Points)):
        Nodes.append(Node(Points[i], V0, i))

    contact_node =  14 #4 #12

    Nodes[contact_node] = Node(Points[contact_node], V1, contact_node)
    
    Springs = {(0, 4),
     (0, 5),
     (0, 6),
     (0, 7),
     (0, 8),
     (0, 11),
     (0, 20),
     (1, 9),
     (1, 13),
     (1, 15),
     (1, 17),
     (1, 22),
     (2, 3),
     (2, 10),
     (2, 13),
     (2, 19),
     (3, 7),
     (3, 12),
     (3, 19),
     (4, 5),
     (4, 8),
     (5, 7),
     (6, 8),
     (6, 11),
     (6, 14),
     (6, 21),
     (7, 12),
     (7, 20),
     (8, 14),
     (9, 11),
     (9, 15),
     (9, 16),
     (9, 17),
     (9, 18),
     (10, 13),
     (10, 14),
     (10, 17),
     (10, 18),
     (10, 22),
     (11, 12),
     (11, 15),
     (11, 16),
     (11, 20),
     (11, 21),
     (12, 15),
     (12, 19),
     (12, 20),
     (13, 15),
     (13, 19),
     (13, 22),
     (14, 16),
     (14, 18),
     (14, 21),
     (15, 19),
     (16, 18),
     (16, 21),
     (17, 18),
     (17, 22)}
    
    # Springs = {(0, 1),
    #   (0, 5),
    #   (0, 8),
    #   (0, 12),
    #   (1, 5),
    #   (1, 6),
    #   (1, 11),
    #   (1, 12),
    #   (2, 3),
    #   (2, 14),
    #   (2, 17),
    #   (2, 18),
    #   (3, 4),
    #   (3, 7),
    #   (3, 14),
    #   (3, 17),
    #   (4, 7),
    #   (4, 9),
    #   (4, 11),
    #   (4, 15),
    #   (5, 11),
    #   (6, 9),
    #   (6, 10),
    #   (6, 11),
    #   (6, 12),
    #   (6, 15),
    #   (7, 10),
    #   (7, 15),
    #   (7, 16),
    #   (7, 17),
    #   (8, 12),
    #   (8, 13),
    #   (8, 14),
    #   (9, 11),
    #   (9, 15),
    #   (10, 12),
    #   (10, 13),
    #   (10, 15),
    #   (10, 16),
    #   (12, 13),
    #   (13, 14),
    #   (13, 16),
    #   (13, 18),
    #   (14, 18),
    #   (16, 17),
    #   (16, 18),
    #   (17, 18)}
    
    # Springs = {(0, 7),
    #   (0, 8),
    #   (0, 14),
    #   (0, 16),
    #   (1, 2),
    #   (1, 10),
    #   (1, 11),
    #   (1, 12),
    #   (2, 4),
    #   (2, 10),
    #   (2, 11),
    #   (3, 5),
    #   (3, 9),
    #   (3, 11),
    #   (3, 12),
    #   (3, 14),
    #   (4, 6),
    #   (4, 9),
    #   (4, 11),
    #   (4, 15),
    #   (5, 8),
    #   (5, 12),
    #   (5, 14),
    #   (6, 9),
    #   (6, 13),
    #   (6, 15),
    #   (6, 16),
    #   (7, 13),
    #   (7, 15),
    #   (7, 16),
    #   (8, 12),
    #   (8, 14),
    #   (9, 11),
    #   (9, 14),
    #   (9, 16),
    #   (11, 12),
    #   (13, 15),
    #   (13, 16),
    #   (14, 16)}
    
    

    k = 1000.
    floe = Floe(nodes=Nodes, springs=Springs,
                stiffness=k, viscosity=k/5, id_number=1)
    
    Traction_Mat = floe.traction_mat()
    Length_Mat = floe.length_mat()
    Torsion_Mat = floe.torsion_mat()
    Angle_Mat = floe.angle_init()
    
    T_end = 1.
    # All_positions_velocities = floe.Move_stable_1(1., Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node).y
    All_positions_velocities = floe.Move(
        T_end, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat).y
    # All_positions_velocities = floe.Move_stable_neighbor(T_end, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node).y

    Energy = Energy_studies(All_positions_velocities, floe)
    #finding time T* when energy elastic is max
    M = np.where(Energy[-1] == max(Energy[-1]))[0][0]
    
    t = np.linspace(0, T_end, 100)
    
    # plt.figure()
    # plt.plot(t, Energy[0], label = "traction")
    # plt.plot(t, Energy[1], label = "torsion")
    # plt.plot(t, Energy[2], label = "elastic")
    # plt.legend()
    # plt.show()
    
    # plt.figure()
    # plt.plot(t[:M], Energy[0][:M])
    # plt.plot(t[:M], Energy[1][:M])
    # plt.plot(t[:M], Energy[2][:M])

    All_pos_vel = All_positions_velocities[:, :M]
    
    deformation_field = np.zeros((floe.n * 2, M))
    
    for i in range(0, floe.n):
        deformation_field[2*i] = All_pos_vel[4*i]
        deformation_field[2*i+1] = All_pos_vel[4*i+1]

    for i in range(floe.n*2):
        deformation_field[i] = deformation_field[i]-deformation_field[i][0]  # compute displacement field

    plt.figure()
    for i in range(len(Points)):
        plt.quiver(Points[i][0], Points[i][1],
                   deformation_field[2*i][M-1], deformation_field[2*i+1][M-1])
    
    #compute the norm of deformation field at each 

    data_deformation = deformation_field[:,-1].reshape(floe.n,2)
    deformation_norm = np.zeros(floe.n)
    for i in range(floe.n):
        deformation_norm[i] = norm(data_deformation[i])
    deformation_norm = np.array([deformation_norm])
    
    x,y = np.array([Points[:,0]]), np.array([Points[:,1]])
    
    #plot norm of deformation at each point 

    
    #plot heat map of the norm of deformation 
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, deformation_norm, c=deformation_norm, alpha=0.5)
    # plt.colorbar()
    plt.show()
    
    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_trisurf(x[0], y[0], deformation_norm[0],
                    cmap='Blues', edgecolor='none');
    plt.show()
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111, autoscale_on=True,
    #                       xlim=(-.2, 1.2), ylim=(-.75, 1.51))
    # ax.set_aspect('equal')
    # ax.grid()
    # line1, = ax.plot([], [], '.-', lw=1.95)
    # time_template = 'time = % 10fs'
    # time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    # Route = floe.Route()

    # def init():
    #     line1.set_data([], [])
    #     time_text.set_text(' ')
    #     return line1, time_text

    # def animate_spring(i):
    #     Ix = [j for j in range(0, floe.n*4, 4)]
    #     Iy = [j for j in range(1, floe.n*4, 4)]
    #     thisx = []
    #     thisy = []
    #     for j in Ix:
    #         thisx = np.append(thisx, All_positions_velocities[j][i])
    #     for j in Iy:
    #         thisy = np.append(thisy, All_positions_velocities[j][i])
    #     for k in Route:
    #         thisx = np.append(thisx, thisx[k])
    #         thisy = np.append(thisy, thisy[k])

    #     line1.set_data(thisx[floe.n:],
    #                     thisy[floe.n:])

    #     time_text.set_text(time_template % (i*dt))
    #     return line1, time_text

    # ani = animation.FuncAnimation(fig, animate_spring,
    #                               np.arange(0, M+1), interval=2, blit=False)



    # floe.plot_displacements_1free(1.,Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node)
    # floe.plot_displacements_Neighborfree(1.,Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node)
