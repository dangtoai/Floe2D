#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 11:23:28 2025

@author: phandangtoai
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 17:16:06 2023

@author: phandangtoai
"""

import csv
import time
import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm
from matplotlib import animation
sys.path.append("/Users/phandangtoai/Documents/Floe2D/Griffith-master_Dimitri")
from griffith.geometry import Point, Polygon, dist
from Func import Node, Floe, Energy_studies, N_T

if __name__ == '__main__':
    Nodes = []
    Points = []

    with open('masses-springs_circle.csv', mode='r') as csv_file:
        csv_reader = csv.reader(csv_file)
        i = 0
        for row in csv_reader:
            row = np.array(row, dtype=float)
            Nodes.append(Node(row, id_number = i))
            i += 1
            Points.append(Point(row[0], row[1]))

    floe = Floe(nodes=Nodes, springs=None, id_number=1)

    # floe.plot_init()
    point_input = Point(1., 0)
    V0 = np.array([-0.1 , 0.])
    polygon = Polygon(Points)
    distances = [dist(point_input, Points[i])
                  for i in floe.border_nodes_index()]
    index_contact = floe.border_nodes_index()[distances.index(min(distances))]
    print("the nearest collision at point ",
          Points[index_contact], " of the floe")
    Nodes[index_contact] = Node(Nodes[index_contact].position(), V0, id_number = index_contact)
    
    

    # Nodes = [Node(np.array([0,-1]),id_number = 0), 
    #           Node(np.array([0,1]), id_number = 1),
    #           Node(np.array([-1,0]),velocity= V0, id_number = 2)]

    # Nodes = [Node(np.array([0,-1]), id_number = 0), 
    #           Node(np.array([0,1]), id_number = 1),
    #           Node(np.array([1, 0.]), id_number = 2),
    #           Node(np.array([-1, 0.]), id_number = 3, velocity= V0)]
    
    
    floe = Floe(nodes=Nodes, springs=None,
                stiffness= 10., viscosity=1., id_number=1, impact_node= True)
    
    # floe.plot_init()
    
    # Simulation of masses-springs network
    Traction_Mat = floe.traction_mat()
    Length_Mat = floe.length_mat()
    Torsion_Mat = floe.torsion_mat()
    Angle_Mat = floe.angle_init()
    # Mass_mat = floe.mass_nodes

    T_END = 1.  # time end

    dt = T_END/N_T  # time's step

    start_time = time.time()
    All_positions_velocities = floe.Move(
        T_END, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat).y
    end_time = time.time()

    print("time of simulation = ", end_time - start_time)

    Energy = Energy_studies(All_positions_velocities, floe)
    # # finding time T* := argmax of the energy elastic
    M = np.where(Energy[-1] == max(Energy[-1]))[0][0]
    # M = N_T

    t = np.linspace(0, T_END, N_T)
    
    # cutting the data after T*
    All_pos_vel = All_positions_velocities[:, :M]
    
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111, autoscale_on=True,
    #                       xlim=(-1.1, 1.1), ylim=(-1.2, 1.1))
    # ax.set_aspect('equal')
    # # ax.grid()
    # line1, = ax.plot([], [], '.-', lw=1.95)
    # time_template = 'time = % 10fs'
    # time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    # # plot the boundary before and after applying the displacement

    Route = floe.border_nodes_index()

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
    #                               np.arange(0, len(t)), interval=20, blit=False)
    # plt.show()
    
    # Route = floe.border_nodes_index()

    ### plot the boundary before and after applying the displacement
    plt.figure()
    pos_init = All_pos_vel[:,0].reshape(floe.n, 4)
    pos_after = All_pos_vel[:,-1].reshape(floe.n, 4)
    plt.plot(pos_init[:,0][Route], pos_init[:,1][Route])
    plt.plot(pos_after[:,0][Route], pos_after[:,1][Route])
    
    # computing displacement field
    deformation_field = np.zeros((floe.n * 2, M))

    for i in range(0, floe.n):
        deformation_field[2*i] = All_pos_vel[4*i]
        deformation_field[2*i+1] = All_pos_vel[4*i+1]

    for i in range(floe.n*2):
        deformation_field[i] = deformation_field[i]-deformation_field[i][0]

    X = np.array([p.x for p in Nodes])
    Y = np.array([p.y for p in Nodes])
    U = [deformation_field[2*i][M-1] for i in range(len(Nodes))]  # List of x-components of deformation
    V = [deformation_field[2*i+1][M-1] for i in range(len(Nodes))]  # List of y-components of deformation
    # plotting the displacement field
    plt.figure()
    plt.quiver(X,Y,U,V)
    plt.show()
    
    data_deformation = deformation_field[:, -1].reshape(floe.n, 2)

    # plot displacement field in the 1st-direction
    # plt.figure()
    # ax = plt.axes(projection='3d')
    # ax.plot_trisurf(X, Y, data_deformation[:, 0],
    #                 cmap='cool', edgecolor='none')
    # plt.title("$u_1$")

#     # plotting displacement field in the 2nd-direction
# #bwr
    # plt.figure()
    # ax = plt.axes(projection='3d')
    # ax.plot_trisurf(X, Y, data_deformation[:, 1], cmap='cool', edgecolor='none')
    # plt.title("$u_2$")

    # Localisation of the deformation
    # deformation_norm = norm(data_deformation, axis=1)

    # plt.figure()
    # ax = plt.axes(projection='3d')
    # ax.plot_trisurf(X, Y, deformation_norm, cmap='Reds', edgecolor='none')
    # plt.title("$||u||$")

    # write csv file to save the displacement field
    # only data on the boundary is needed
    index_boundary = floe.border_nodes_index()[:-1]
    Xboundary, Yboundary = X[index_boundary], Y[index_boundary]
    data_x, data_y = data_deformation[index_boundary][:,0], data_deformation[index_boundary][:,1]

    PRECISION = 5 # precision of the data: 5 decimals.
    Xboundary, Yboundary, data_x, data_y = np.around((Xboundary, Yboundary, data_x, data_y), decimals= PRECISION)

    data = zip(Xboundary, Yboundary, data_x, data_y)

    New_Nodes = []

    with open('boundary_data.csv', mode = 'w', newline = '', encoding='utf-8') as csv_file:
        # write the contact node at the boundary of the network
        csv_writer = csv.writer(csv_file, delimiter=' ')
        csv_writer.writerow( "Contact region in coninuum domain:" )
        csv_writer.writerow((Points[index_contact].x, Points[index_contact].y))
        csv_writer.writerow(['X', 'Y', 'Data_X', 'Data_Y'])
        csv_writer.writerows(data)

    #new ice floe without the boundary
    for i in sorted(index_boundary, reverse = True):
        Nodes.remove(Nodes[i])
        X = np.delete(X, i)
        Y = np.delete(Y, i)
        data_deformation = np.delete(data_deformation, i, axis=0)

    for i in range(len(Nodes)):
        New_Nodes.append(Node(Nodes[i].position(), id_number = i ))
        i += 1

    floe1 = Floe(nodes=New_Nodes, springs=None, id_number=0)

    index_boundary = floe1.border_nodes_index()[:-1]
    Xboundary, Yboundary = X[index_boundary], Y[index_boundary]
    data_x, data_y = data_deformation[index_boundary][:,0], data_deformation[index_boundary][:,1]
    PRECISION = 5 # data's precision: 5 decimals.
    Xboundary, Yboundary, data_x, data_y = np.around((Xboundary, Yboundary, data_x, data_y), decimals= PRECISION)
    data = zip(Xboundary, Yboundary, data_x, data_y)

    with open('boundary_data.csv', mode = 'a', newline = '', encoding='utf-8') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=' ')
        csv_writer.writerows(data)

