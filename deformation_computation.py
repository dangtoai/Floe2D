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
sys.path.append("/Users/phandangtoai/Documents/Floe2D/Griffith-master_Dimitri")
from griffith.geometry import Point, Polygon, dist
from Func import Node, Floe, Energy_studies

if __name__ == '__main__':
    Nodes = []
    Points = []

    with open('masses-springs.csv', mode='r') as csv_file:
        csv_reader = csv.reader(csv_file)
        i = 0
        for row in csv_reader:
            row = np.array(row, dtype=float)
            Nodes.append(Node(row,id_number = i))
            i += 1
            Points.append(Point(row[0], row[1]))

    floe = Floe(nodes=Nodes, springs=None,
                stiffness=1, viscosity=1., id_number=1)

    Input_pos = input(
        "Enter a point separated by a space (e.g., 'a b') to run the simulation: ")
    Input_vel = input(
        "Enter a velocity separated by a space (e.g., 'a b') to run the simulation: ")
    point_input = Input_pos.split()
    vel_input = Input_vel.split()
    if len(point_input) == len(vel_input) == 2:
        a, b = float(point_input[0]), float(point_input[1])
        point_input = Point(a, b)
        a, b = float(vel_input[0]), float(vel_input[1])
        print("the point of contact is: ")
        print(point_input)
        V0 = np.array([a, b])
    else:
        print("Invalid format.")

    polygon = Polygon(Points)
    distances = [dist(point_input, Points[i])
                 for i in floe.border_nodes_index()]
    index_contact = distances.index(min(distances))
    print("the nearest collision at point ",
          Points[index_contact], " of the floe")
    Nodes[index_contact] = Node(Nodes[index_contact].position(), V0)

    # Simulation of masses-springs network
    Traction_Mat = floe.traction_mat()
    Length_Mat = floe.length_mat()
    Torsion_Mat = floe.torsion_mat()
    Angle_Mat = floe.angle_init()

    T_END = 0.1  # time end
    N_T = 20        #
    dt = 1./N_T  # time's step
    start_time = time.time()
    All_positions_velocities = floe.Move(
        T_END, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat).y
    end_time = time.time()

    print("time of simulation = ", end_time - start_time)

    Energy = Energy_studies(All_positions_velocities, floe)
    # finding time T* when energy elastic is max
    M = np.where(Energy[-1] == max(Energy[-1]))[0][0]

    t = np.linspace(0, T_END, N_T)

    All_pos_vel = All_positions_velocities[:, :M]

    # computing displacement field
    deformation_field = np.zeros((floe.n * 2, M))

    for i in range(0, floe.n):
        deformation_field[2*i] = All_pos_vel[4*i]
        deformation_field[2*i+1] = All_pos_vel[4*i+1]

    for i in range(floe.n*2):
        deformation_field[i] = deformation_field[i]-deformation_field[i][0]

    X = np.array([p.x for p in Points])
    Y = np.array([p.y for p in Points])

    # plotting the displacement field
    plt.figure()
    for i in range(len(Points)):
        plt.quiver(X[i], Y[i], deformation_field[2*i]
                   [M-1], deformation_field[2*i+1][M-1])

    data_deformation = deformation_field[:, -1].reshape(floe.n, 2)

    # plot displacement field in the 1st-direction
    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_trisurf(X, Y, data_deformation[:, 0],
                    cmap='RdYlBu', edgecolor='none')
    plt.title("$u_1$")

    # plotting displacement field in the 2nd-direction

    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_trisurf(X, Y, data_deformation[:, 1], cmap='bwr', edgecolor='none')
    plt.title("$u_2$")

    # Localisation of the deformation
    deformation_norm = norm(data_deformation, axis=1)

    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_trisurf(X, Y, deformation_norm, cmap='Reds', edgecolor='none')
    plt.title("$||u||$")

    # write csv file to save the displacement field
    # only data on the boundary is needed
    index_boundary = floe.border_nodes_index()
    Xboundary, Yboundary = X[index_boundary], Y[index_boundary]
    data_x, data_y = data_deformation[index_boundary][:,0], data_deformation[index_boundary][:,1]

    PRECISION = 5 # precision of the data: 5 decimals.
    Xboundary, Yboundary, data_x, data_y = np.around((Xboundary, Yboundary, data_x, data_y), decimals= PRECISION)

    data = zip(Xboundary, Yboundary, data_x, data_y)

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
        np.delete(X, i)
        np.delete(Y, i)
        np.delete(data_deformation, i)

    floe1 = Floe(nodes=Nodes, springs=None,
                stiffness=1, viscosity=1., id_number=0)

    index_boundary = floe1.border_nodes_index()
    Xboundary, Yboundary = X[index_boundary], Y[index_boundary]
    data_x, data_y = data_deformation[index_boundary][:,0], data_deformation[index_boundary][:,1]
    PRECISION = 5 # precision of the data: 5 decimals.
    Xboundary, Yboundary, data_x, data_y = np.around((Xboundary, Yboundary, data_x, data_y), decimals= PRECISION)
    data = zip(Xboundary, Yboundary, data_x, data_y)

    with open('boundary_data.csv', mode = 'a', newline = '', encoding='utf-8') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=' ')
        csv_writer.writerows(data)
