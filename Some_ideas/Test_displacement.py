#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:15:23 2023

@author: phandangtoai
"""

import os, sys
sys.path.append("/Users/phandangtoai/Documents/Floe2D/")

from Func import *
from graph import *
import numpy as np
from numpy import array
# import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib import transforms
import time
from matplotlib import animation

if __name__ == '__main__':
    start_time = time.time()
    
    
    ######################
    ##### An ice floe ####
    ######################

    # first floe
    
    # Points =     np.array([[0.20671916, 0.84355497],
    #        [0.91861091, 0.34602815],
    #        [0.48841119, 0.10082727],
    #        [0.61174386, 0.38340907],
    #        [0.76590786, 0.5103548 ],
    #        [0.51841799, 0.96110308],
    #        [0.2968005 , 0.37151262],
    #        [0.18772123, 0.01236941],
    #        [0.08074127, 0.85970689],
    #        [0.7384403 , 0.11111075],
    #        [0.44130922, 0.47833904],
    #        [0.15830987, 0.84998003],
    #        [0.87993703, 0.51473797],
    #        [0.27408646, 0.44660783],
    #        [0.41423502, 0.80047642],
    #        [0.29607993, 0.02039138],
    #        [0.62878791, 0.57261865],
    #        [0.57983781, 0.41138362],
    #        [0.5999292 , 0.9851368 ],
    #        [0.26581912, 0.80140153],
    #        [0.28468588, 0.0539621 ],
    #        [0.25358821, 0.19047777],
    #        [0.32756395, 0.45241885],
    #        [0.1441643 , 0.70294208],
    #        [0.16561286, 0.33204815],
    #        [0.96393053, 0.3599832 ],
    #        [0.96022672, 0.92147057],
    #        [0.18841466, 0.95363051],
    #        [0.02430656, 0.40768573],
    #        [0.20455555, 0.89857115],
    #        [0.69984361, 0.33025325],
    #        [0.77951459, 0.08273857],
    #        [0.02293309, 0.52671757],
    #        [0.57766286, 0.66084439],
    #        [0.00164217, 0.89298429],
    #        [0.51547261, 0.96515755],
    #        [0.63979518, 0.76993268],
    #        [0.9856244 , 0.75909912],
    #        [0.2590976 , 0.71004905],
    #        [0.80249689, 0.70155283],
    #        [0.87048309, 0.76733138],
    #        [0.92274961, 0.97434948],
    #        [0.00221421, 0.37371452],
    #        [0.46948837, 0.08305444],
    #        [0.98146874, 0.23964044],
    #        [0.3989448 , 0.22148276],
    #        [0.81373248, 0.3635998 ],
    #        [0.5464565 , 0.81031423],
    #        [0.77085409, 0.06009617],
    #        [0.48493107, 0.44974356],
    #        [0.02911156, 0.81317053],
    #        [0.08652569, 0.26423837],
    #        [0.11145381, 0.06340098],
    #        [0.25124511, 0.24210767],
    #        [0.96491529, 0.08507046],
    #        [0.63176605, 0.80777744],
    #        [0.8166602 , 0.17025008],
    #        [0.566082  , 0.19534463],
    #        [0.63535621, 0.81464201],
    #        [0.81190239, 0.81028553],
    #        [0.92668262, 0.58937388],
    #        [0.91262676, 0.91473434],
    #        [0.82481072, 0.05982164],
    #        [0.09420273, 0.96499664],
    #        [0.36104842, 0.57097522],
    #        [0.03550903, 0.30251811],
    #        [0.54635835, 0.82570583],
    #        [0.79614272, 0.65941728],
    #        [0.0511428 , 0.98650144],
    #        [0.18866774, 0.10748953],
    #        [0.36547777, 0.58091853],
    #        [0.24429087, 0.47282816],
    #        [0.79508747, 0.65227079],
    #        [0.35209494, 0.24185905],
    #        [0.63887768, 0.03113006],
    #        [0.49341505, 0.54423235],
    #        [0.58349974, 0.36471018],
    #        [0.93929935, 0.89233328],
    #        [0.94354008, 0.45934559],
    #        [0.11169243, 0.4186541 ]])*100
    
    Points =  [array([1027.51049423,  115.02981638]),
     array([1008.83061726,  120.51666235]),
     array([1037.78820906,   97.13480391]),
     array([1035.02215601,  120.76515947]),
     array([1015.68515921,   99.63051419]),
     array([1008.20004852,  121.6164766 ]),
     array([1022.86220681,  116.02585528]),
     array([1020.06117011,  111.11561771]),
     array([1019.39223425,  119.73798945]),
     array([1034.5752819 ,  114.36058709]),
     array([1024.9244897 ,  119.61649976]),
     array([1013.81962978,  103.12252245]),
     array([1014.55549189,  100.70315335]),
     array([1024.7685969 ,  100.23957787]),
     array([1024.37984484,  105.15353873]),
     array([1028.76624654,  100.92891029]),
     array([1016.59149878,  122.00517625]),
     array([1031.59989658,  122.63117248]),
     array([1021.95447836,  113.36603073]),
     array([1033.02148787,  112.40002039]),
     array([1035.02438626,  101.97046968]),
     array([1010.84974601,  114.47814689]),
     array([1035.00713767,  101.46599471]),
     array([1015.46557085,  102.99019552]),
     array([1028.24151437,  119.54719156]),
     array([1006.35386329,  119.65663572]),
     array([1023.9693291 ,   97.14951044]),
     array([1029.39535251,  109.69193681]),
     array([1036.30495577,  105.09943046]),
     array([1035.8004482 ,  119.02325002]),
     array([1033.06029494,  121.94473199]),
     array([1039.95380338,  111.17799553]),
     array([1039.59730348,  115.40706503]),
     array([1018.07660579,  117.01292055]),
     array([1008.60974948,  113.92335679]),
     array([1013.46135659,  102.28110259]),
     array([1032.99694059,  101.18457632]),
     array([1034.26970065,   94.17329597]),
     array([1029.23708876,  105.10960151]),
     array([1016.25233069,   98.98781346]),
     array([1029.09644768,  104.41298897]),
     array([1010.60675063,  109.76496401]),
     array([1020.33058072,  103.46937577]),
     array([1033.29518656,  100.42834527]),
     array([1033.85347247,   95.44271078]),
     array([1033.24827229,  106.13327542]),
     array([1017.18898229,  104.29408332]),
     array([1014.73128549,  106.76487661]),
     array([1014.79351745,  104.41239607]),
     array([1011.5676236 ,  122.92165197]),
     array([1023.20116724,  107.80270782]),
     array([1028.35979043,  115.88399794]),
     array([1024.54882415,  126.88045697]),
     array([1032.80378498,  119.18948067]),
     array([1018.73962309,   98.92552814]),
     array([1029.29404047,   94.36506344]),
     array([1029.53367119,  106.76741578]),
     array([1011.83690027,  122.69025661]),
     array([1022.51156849,  107.5503004 ]),
     array([1023.48641538,  123.4936645 ]),
     array([1008.88565605,  109.21272937]),
     array([1025.0246706 ,   97.21917478]),
     array([1023.91349025,  125.71626308]),
     array([1013.05724619,  117.40141402]),
     array([1013.76645421,  106.87863897]),
     array([1022.07369571,   98.03843031]),
     array([1036.41226329,  108.36667996]),
     array([1010.84795912,  110.93606847])]
    
    Points = np.array([list(el) for el in Points])
    
    Nodes = []
    V0 = np.array([0., 0.])
    V1 = np.array([10., 10])

    for i in range(len(Points)):
        Nodes.append(Node(Points[i], V0, i))

    contact_node = 44 #8 #14 #4 #12
    contact_node = 25

    Nodes[contact_node] = Node(Points[contact_node], V1, contact_node)
    

    # Springs = {(0, 11),
    #  (0, 19),
    #  (0, 23),
    #  (0, 29),
    #  (1, 25),
    #  (1, 44),
    #  (1, 46),
    #  (1, 56),
    #  (1, 78),
    #  (2, 43),
    #  (2, 45),
    #  (2, 57),
    #  (2, 74),
    #  (3, 4),
    #  (3, 16),
    #  (3, 17),
    #  (3, 30),
    #  (3, 76),
    #  (4, 12),
    #  (4, 16),
    #  (4, 30),
    #  (4, 46),
    #  (4, 72),
    #  (5, 14),
    #  (5, 18),
    #  (5, 35),
    #  (5, 66),
    #  (6, 10),
    #  (6, 13),
    #  (6, 22),
    #  (6, 24),
    #  (6, 49),
    #  (6, 53),
    #  (6, 73),
    #  (7, 15),
    #  (7, 20),
    #  (7, 52),
    #  (7, 69),
    #  (7, 74),
    #  (8, 11),
    #  (8, 23),
    #  (8, 34),
    #  (8, 50),
    #  (8, 63),
    #  (9, 30),
    #  (9, 31),
    #  (9, 48),
    #  (9, 56),
    #  (9, 57),
    #  (9, 74),
    #  (10, 22),
    #  (10, 49),
    #  (10, 64),
    #  (10, 75),
    #  (11, 23),
    #  (11, 27),
    #  (11, 29),
    #  (11, 63),
    #  (12, 46),
    #  (12, 60),
    #  (12, 72),
    #  (12, 78),
    #  (13, 22),
    #  (13, 24),
    #  (13, 71),
    #  (14, 19),
    #  (14, 29),
    #  (14, 33),
    #  (14, 35),
    #  (14, 38),
    #  (14, 47),
    #  (14, 66),
    #  (14, 70),
    #  (15, 20),
    #  (15, 43),
    #  (15, 74),
    #  (16, 17),
    #  (16, 33),
    #  (16, 72),
    #  (16, 75),
    #  (17, 49),
    #  (17, 75),
    #  (17, 76),
    #  (18, 35),
    #  (18, 41),
    #  (18, 58),
    #  (18, 59),
    #  (18, 61),
    #  (18, 66),
    #  (18, 68),
    #  (19, 23),
    #  (19, 29),
    #  (19, 38),
    #  (20, 21),
    #  (20, 43),
    #  (20, 45),
    #  (20, 69),
    #  (21, 45),
    #  (21, 51),
    #  (21, 53),
    #  (21, 69),
    #  (21, 73),
    #  (22, 64),
    #  (22, 71),
    #  (23, 32),
    #  (23, 38),
    #  (23, 50),
    #  (23, 71),
    #  (24, 51),
    #  (24, 53),
    #  (24, 65),
    #  (24, 71),
    #  (24, 79),
    #  (25, 37),
    #  (25, 44),
    #  (25, 78),
    #  (26, 37),
    #  (26, 41),
    #  (26, 61),
    #  (26, 77),
    #  (27, 29),
    #  (27, 35),
    #  (27, 63),
    #  (27, 68),
    #  (28, 32),
    #  (28, 42),
    #  (28, 65),
    #  (28, 79),
    #  (29, 35),
    #  (30, 46),
    #  (30, 56),
    #  (30, 57),
    #  (30, 76),
    #  (31, 48),
    #  (31, 56),
    #  (31, 62),
    #  (32, 34),
    #  (32, 42),
    #  (32, 50),
    #  (32, 71),
    #  (32, 79),
    #  (33, 36),
    #  (33, 47),
    #  (33, 70),
    #  (33, 72),
    #  (33, 75),
    #  (34, 42),
    #  (34, 50),
    #  (34, 63),
    #  (34, 68),
    #  (35, 68),
    #  (36, 39),
    #  (36, 47),
    #  (36, 55),
    #  (36, 58),
    #  (36, 59),
    #  (36, 67),
    #  (36, 72),
    #  (37, 40),
    #  (37, 44),
    #  (37, 60),
    #  (37, 77),
    #  (37, 78),
    #  (38, 64),
    #  (38, 70),
    #  (38, 71),
    #  (39, 40),
    #  (39, 59),
    #  (39, 60),
    #  (39, 67),
    #  (40, 59),
    #  (40, 60),
    #  (40, 77),
    #  (41, 61),
    #  (42, 52),
    #  (42, 65),
    #  (43, 45),
    #  (43, 74),
    #  (44, 54),
    #  (44, 56),
    #  (45, 49),
    #  (45, 57),
    #  (45, 73),
    #  (45, 76),
    #  (46, 56),
    #  (46, 78),
    #  (47, 55),
    #  (47, 66),
    #  (48, 62),
    #  (48, 74),
    #  (49, 73),
    #  (49, 75),
    #  (49, 76),
    #  (51, 52),
    #  (51, 53),
    #  (51, 65),
    #  (51, 69),
    #  (52, 65),
    #  (52, 69),
    #  (53, 73),
    #  (54, 56),
    #  (54, 62),
    #  (55, 58),
    #  (55, 66),
    #  (56, 62),
    #  (57, 74),
    #  (57, 76),
    #  (58, 59),
    #  (58, 66),
    #  (59, 61),
    #  (59, 77),
    #  (60, 67),
    #  (60, 72),
    #  (60, 78),
    #  (61, 77),
    #  (62, 74),
    #  (63, 68),
    #  (64, 70),
    #  (64, 71),
    #  (64, 75),
    #  (65, 79),
    #  (67, 72),
    #  (70, 75),
    #  (71, 79)}
    

    k = 1000.
    floe = Floe(nodes=Nodes, springs=None,
                stiffness=k, viscosity=k/5, id_number=1)
    
    # floe.plot_border()
    
    Traction_Mat = floe.traction_mat()
    Length_Mat = floe.length_mat()
    Torsion_Mat = floe.torsion_mat()
    Angle_Mat = floe.angle_init()
    
    T_end = 0.1
    dt = 1./N_T
    # All_positions_velocities = floe.Move_stable_1(1., Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node).y
    All_positions_velocities = floe.Move(
        T_end, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat).y
    # All_positions_velocities = floe.Move_stable_neighbor(T_end, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node).y

    Energy = Energy_studies(All_positions_velocities, floe)
    #finding time T* when energy elastic is max
    M = np.where(Energy[-1] == max(Energy[-1]))[0][0]
    
    t = np.linspace(0, T_end, N_T)
    
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
    end_time = time.time()
    elapsed_time = end_time - start_time
    
    print("time of computation for this simulation:", elapsed_time)
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
    

    
    # DD = [deformation_field[2*i][M-1] for i in range(len(Points))]
    # DDD = [deformation_field[2*i+1][M-1] for i in range(len(Points))]
    # plt.figure()
    # plt.quiver(Points[:,0], Points[:,1], DD, DDD)
    #compute the norm of deformation field at each 

    data_deformation = deformation_field[:,-1].reshape(floe.n,2) #
    x,y = np.array([Points[:,0]]), np.array([Points[:,1]])
    #plot displacement field in the x-direction
    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_trisurf(x[0], y[0], data_deformation[:,0], 
                    cmap='RdYlBu', edgecolor='none')
    plt.title("$u_1$")
    plt.show()
    
    
    #plot displacement field in the y-direction
    
    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_trisurf(x[0], y[0], data_deformation[:,1],
                    cmap='bwr', edgecolor='none');
    plt.title("$u_2$")
    plt.show()
    
    
    deformation_norm = np.zeros(floe.n)
    for i in range(floe.n):
        deformation_norm[i] = norm(data_deformation[i])
    deformation_norm = np.array([deformation_norm])
    
    DD = np.array([deformation_field[2*i][M-1] for i in range(len(Points))])
    DDD = np.array([deformation_field[2*i+1][M-1] for i in range(len(Points))])
    plt.figure()
    plt.quiver(Points[:,0], Points[:,1], DD, DDD)
    plt.tight_layout()
    
    #only plot displacement at border : 
    plt.figure()
    plt.quiver(Points[:,0][floe.border_nodes_index()], Points[:,1][floe.border_nodes_index()], 
               DD[floe.border_nodes_index()], DDD[floe.border_nodes_index()])
    plt.tight_layout()
    # x,y = np.array([Points[:,0]]), np.array([Points[:,1]])
    
    #plot norm of deformation at each point 
    #plot heat map of the norm of deformation 
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, deformation_norm, c=deformation_norm, alpha=0.5)
    plt.show()
    
    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_trisurf(x[0], y[0], deformation_norm[0],
                    cmap='Reds', edgecolor='none')
    plt.title("Displacement's norm")
    plt.show()
    
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=True,
                          # xlim=(-0.2, 101.2), ylim=(-0.75, 101.51))
                          xlim=(1000.2, 1040.2), ylim=(90, 130.51))
    ax.set_aspect('equal')
    ax.grid()
    line1, = ax.plot([], [], '.-', lw=1.95)
    time_template = 'time = % 10fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    Route = floe.border_nodes_index()

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

        line1.set_data(thisx[floe.n:],
                        thisy[floe.n:])

        time_text.set_text(time_template % (i*dt))
        return line1, time_text

    ani = animation.FuncAnimation(fig, animate_spring,
                                  np.arange(0, len(t)), interval=200, blit=False)
    plt.show()
    

    # floe.plot_displacements_1free(1.,Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node)
    # floe.plot_displacements_Neighborfree(1.,Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat, contact_node)
