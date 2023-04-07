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
    # Points = np.array([[0.25298236, 0.59733394],
    #   [0.43479153, 0.0089861 ],
    #   [0.77938292, 0.38657128],
    #   [0.19768507, 0.04416006],
    #   [0.86299324, 0.95665297],
    #   [0.98340068, 0.43614665],
    #   [0.16384224, 0.94897731]])

    Points = np.array([[0.97627445, 0.8662893],
                       [0.00623026, 0.17316542],
                       [0.25298236, 0.07494859],
                       [0.43479153, 0.60074272],
                       [0.77938292, 0.16797218],
                       [0.19768507, 0.73338017],
                       [0.86299324, 0.40844386],
                       [0.98340068, 0.52790882],
                       [0.16384224, 0.93757158],
                       [0.59733394, 0.52169612],
                       [0.0089861, 0.10819338],
                       [0.38657128, 0.15822341],
                       [0.04416006, 0.54520265],
                       [0.95665297, 0.52440408],
                       [0.43614665, 0.63761024],
                       [0.94897731, 0.40149544],
                       [0.78630599, 0.64980511]])

    Nodes = []
    V0 = np.array([0., 0.])
    V1 = np.array([3.5, 3.5])
    # V1 = np.array([-30.5, 0.])
    for i in range(len(Points)):
        Nodes.append(Node(Points[i], V0, i))

    contact_node = 10

    Nodes[contact_node] = Node(Points[contact_node], V1, contact_node)
    # Springs = {(0, 1),
    #   (0, 2),
    #   (0, 3),
    #   (0, 4),
    #   (0, 6),
    #   (1, 2),
    #   (1, 3),
    #   (1, 5),
    #   (2, 4),
    #   (2, 5),
    #   (3, 6),
    #   (4, 5),
    #   (4, 6)}

    
    Springs = {(0, 7),
     (0, 8),
     (0, 14),
     (0, 16),
     (1, 2),
     (1, 10),
     (1, 11),
     (1, 12),
     (2, 4),
     (2, 10),
     (2, 11),
     (3, 5),
     (3, 9),
     (3, 11),
     (3, 12),
     (3, 14),
     (4, 6),
     (4, 9),
     (4, 11),
     (4, 15),
     (5, 8),
     (5, 12),
     (5, 14),
     (6, 9),
     (6, 13),
     (6, 15),
     (6, 16),
     (7, 13),
     (7, 15),
     (7, 16),
     (8, 12),
     (8, 14),
     (9, 11),
     (9, 14),
     (9, 16),
     (11, 12),
     (13, 15),
     (13, 16),
     (14, 16)}

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
    M = np.where(Energy[-1] == max(Energy[-1]))[0][0]

    t = np.linspace(0, T_end, 500)

    # plt.figure()
    # plt.plot(t, Energy[0])
    # plt.plot(t, Energy[1])
    # plt.plot(t, Energy[2])

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
    plt.colorbar()
    plt.show()
    
    # plt.figure()
    # ax = plt.axes(projection='3d')
    # ax.plot_trisurf(x[0], y[0], deformation_norm[0],
    #                 cmap='reds', edgecolor='none');
    # plt.show()
    # fig = plt.figure()
    # ax = fig.add_subplot(111, autoscale_on=True,
    #                      xlim=(-.2, 1.2), ylim=(-.75, 1.51))
    # ax.set_aspect('equal')
    # ax.grid()
    # # plt.axvline(x=2., color="red")
    # line1, = ax.plot([], [], '.-', lw=1.95)
    # # line2, = ax.plot([], [], '.-', lw=1.95)
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
    #                    thisy[floe.n:])

    #     time_text.set_text(time_template % (i*dt))
    #     return line1, time_text

    # ani = animation.FuncAnimation(fig, animate_spring,
    #                               np.arange(0, 500), interval=2, blit=False)

    # plt.figure()
    # Traction_energy, Torsion_energy, Total_energy = Energy_studies(All_positions_velocities, floe)
    # # Traction_energy, Torsion_energy, Total_energy = floe.energy_evolution_stable(1., contact_node)
    # # Traction_energy, Torsion_energy, Total_energy = floe.energy_evolution(1.)
    # t = np.linspace(0,1,800)
    # plt.plot(t, Traction_energy, label = "Traction energy")
    # plt.plot(t, Torsion_energy, label = "Torsion energy")
    # plt.plot(t, Total_energy, label = "Total energy")
    # plt.legend()

    # plt.show()

    # floe.plot_displacements_1free(1.,Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node)
    # floe.plot_displacements_Neighborfree(1.,Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node)
