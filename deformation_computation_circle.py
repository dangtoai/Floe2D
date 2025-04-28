#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 11:23:28 2025

@author: phandangtoai
"""

from mpi4py import MPI
import csv
import signal
import time
import sys
import numpy as np
import matplotlib.pyplot as plt
# from numpy.linalg import norm
from matplotlib import animation
sys.path.append("/Users/phandangtoai/Documents/Floe2D/Griffith-master_Dimitri")
from griffith.geometry import Point, Polygon, dist
from Func import Angle_mat, Traction_mat, Length_mat, Torsion_mat
from Func import Node, Floe, Energy_studies, N_T, timeout_handler, T_LIMIT
# from Circle_geometry_simulation import radius


if __name__ == '__main__':

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()  # Process ID
    size = comm.Get_size()  # Total number of processes

    radius = 100.
    Nodes = []
    Points = []

    rho = 900. #kg/m^3 density of sea ice
    Volume = np.pi * radius**2

    Mass = Volume * rho

    matrix = np.array([[32*Volume/(9*np.pi**2), 3*Volume/4],
                        [32*Volume/(9*np.pi**2), -3*Volume/4]])
    
    ### Young modulus E and Poisson ration nu
    
    E, nu = 8.95*10**9, 0.295
    
    ### Lame parameters
    lamb, mu = E*nu/ ((1+nu)*(1-2*nu)), E/(2*(1+nu))
        
    ### K,G springs stiffness
    traction_stiff, torsion_stiff = np.linalg.solve(matrix, np.array([lamb, mu])) 
    # torsion_stiff = 0
    

    filename = f"masses-springs_limit_circle_{rank+1}.csv"
    
    print(f" reading network from Process {rank} in {filename}")
    
    with open(filename, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file)
        i = 0
        for row in csv_reader:
            row = np.array(row, dtype=float)
            Nodes.append(Node(row, id_number = i))
            i += 1
            Points.append(Point(row[0], row[1]))

    floe = Floe(nodes=Nodes, springs=None, id_number=1)

    # # floe.plot_init()
    point_input = Point(radius, 0)
    V0 = np.array([-1. , 0.])
    # polygon = Polygon(Points)
    distances = [dist(point_input, Points[i])
                  for i in floe.border_nodes_index()]
    index_contact = floe.border_nodes_index()[distances.index(min(distances))]
    # # print("the nearest collision at point ",
    #       # Points[index_contact], " of the floe")
    Nodes[index_contact] = Node(Nodes[index_contact].position(), V0, id_number = index_contact)
    floe = Floe(nodes=Nodes, springs=None, mass= Mass,
                stiffness= traction_stiff, torsion_stiff= torsion_stiff ,viscosity=1000, id_number=1, impact_node= True)
    
    
    
    # V0 = np.array([10. , 0.])

    # Nodes = [Node(np.array([0,-1]), id_number = 0), 
    #           Node(np.array([0,1]), id_number = 1),
    #           Node(np.array([1,0]),velocity= V0, id_number = 2)]
    # index_contact = 2

    # Nodes = [Node(np.array([0,-1]), id_number = 0), 
    #           Node(np.array([0,1]), id_number = 1),
    #           Node(np.array([1, 0.]), V0, id_number = 2),
    #           Node(np.array([-1, 0.]), id_number = 3)]
    
    # floe = Floe(nodes=Nodes, springs=None, mass= 1000,
    #             stiffness= traction_stiff, torsion_stiff= torsion_stiff ,viscosity=1000, id_number=1, impact_node= True)

    
    # floe.plot_init()
    # Simulation of masses-springs network
    
    Angle_Mat = Angle_mat(floe)
    Traction_Mat = Traction_mat(floe, Angle_Mat)
    Length_Mat = Length_mat(floe)
    Torsion_Mat = Torsion_mat(floe)
    
    
    # Traction_Mat = floe.traction_mat()
    # Length_Mat = floe.length_mat()
    # Torsion_Mat = floe.torsion_mat()
    # Angle_Mat = floe.angle_init()
    
    
    
    # Mass_mat = floe.mass_nodes
    
    #### compute effective stiffness
    Neighbors = floe.Neighbors()[floe.n-1]
    K_neighbors = [Traction_Mat[floe.n-1, i] for i in Neighbors]
    K_eff = sum(K_neighbors)
    
    print("max displacement of q_0:", -np.sqrt(floe.mass_nodes[-1]/K_eff))
    
    T_END = 1.5  # time end

    dt = T_END/N_T  # time's step
    
    ### control execution time. 
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(T_LIMIT)  # Trigger timeout after T_LIMIT seconds
    
    try:
        start_time = time.time()

        # Run the simulation (this is now interruptible)
        All_positions_velocities = floe.Move(
            T_END, Traction_Mat, Length_Mat, Torsion_Mat, Angle_Mat).y

        end_time = time.time()
        print("time of simulation = ", end_time - start_time)

        # Disable alarm after successful execution
        signal.alarm(0)
        
    except TimeoutError as e: print(e)

    Energy = Energy_studies(All_positions_velocities, floe)
    
    ##### finding time T* := argmax of the energy elastic
    M = np.where(Energy[-1] == max(Energy[-1]))[0][0]
    print("characteristic collision time = ", M*T_END/N_T)
    # M = N_T
    if M == N_T: 
        print(f"Stopping execution because M == {N_T}")
        sys.exit()
    
    t = np.linspace(0, T_END, N_T)
    
    # cutting the data after T*
    All_pos_vel = All_positions_velocities[:, :M]
    
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111, autoscale_on=True,
    #                       xlim=(-radius-1., radius + 2.), ylim=(-radius-1.2,radius+ 1.1))
    
    # ax.set_aspect('equal')
    # # ### ax.grid()
    # line1, = ax.plot([], [], '.-', lw=1.95)
    # time_template = 'time = % 10fs'
    # time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    # ### plot the boundary before and after applying the displacement

    # Route = floe.border_nodes_index()

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
    

    ### plot the boundary before and after applying the displacement
    
    # plt.figure()
    # pos_init = All_pos_vel[:,0].reshape(floe.n, 4)
    # pos_after = All_pos_vel[:,-1].reshape(floe.n, 4)
    # plt.plot(pos_init[:,0][Route], pos_init[:,1][Route])
    # plt.plot(pos_after[:,0][Route], pos_after[:,1][Route])
    
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
    
    print("real displacement of q_0:",U[-1])
    ### plotting the displacement field
    # plt.figure()
    # plt.quiver(X,Y,U,V,angles='xy', scale_units='xy', scale=0.03)
    # plt.show()
    
    data_deformation = deformation_field[:, -1].reshape(floe.n, 2)

    ### plot displacement field in the 1st-direction
    # min_val = min(np.min(data_deformation[:, 0]), np.min(data_deformation[:, 1]))
    # max_val = max(np.max(data_deformation[:, 0]), np.max(data_deformation[:, 1]))

    # fig = plt.figure(figsize=(12, 6))

    # First subplot (for u1)
    # ax1 = fig.add_subplot(121, projection='3d')
    # trisurf1 = ax1.plot_trisurf(X, Y, data_deformation[:, 0], 
    #                         cmap='Reds', edgecolor='none', vmin=min_val, vmax=max_val)
    # ax1.set_title("$u_1$")

    # # Second subplot (for u2)
    # ax2 = fig.add_subplot(122, projection='3d')
    # trisurf2 = ax2.plot_trisurf(X, Y, data_deformation[:, 1], 
    #                         cmap='Reds', edgecolor='none', vmin=min_val, vmax=max_val)
    # ax2.set_title("$u_2$")

    # # Add a colorbar that is shared by both plots
    # fig.colorbar(trisurf1, ax=[ax1, ax2], location = "right", shrink=0.6, aspect=40)
    # plt.show()

    #### Localisation of the deformation
    # deformation_norm = norm(data_deformation, axis=1)

    # plt.figure()
    # ax = plt.axes(projection='3d')
    # ax.plot_trisurf(X, Y, deformation_norm, cmap='Reds', edgecolor='none')
    # plt.title("$||u||$")

    ### write csv file to save the displacement field
    ### only data on the boundary is needed
    # index_boundary = floe.border_nodes_index()[:-1]
    
    index_boundary = np.arange(floe.n)
    Xboundary, Yboundary = X[index_boundary], Y[index_boundary]
    data_x, data_y = data_deformation[index_boundary][:,0], data_deformation[index_boundary][:,1]

    PRECISION = 5 # precision of the data: 5 decimals.
    Xboundary, Yboundary, data_x, data_y = np.around((Xboundary, Yboundary, data_x, data_y), decimals= PRECISION)

    data = zip(Xboundary, Yboundary, data_x, data_y)
    
    filename = f"data_spring_{rank+1}.txt"
    with open(filename, mode = 'a', newline = '', encoding='utf-8') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=' ')
        csv_writer.writerow(['T* ', M*T_END/N_T])
        csv_writer.writerows(data)
        csv_writer.writerow("-----------------------------------------")
        
        
        
    # print(list(data))
    Xboundary, Yboundary, data_x, data_y = np.around((Xboundary, Yboundary, data_x, data_y), decimals= PRECISION)
    data = zip(Xboundary, Yboundary, data_x, data_y)
    
    filename = f"boundary_data_limit_circle_{rank+1}.csv"
    
    with open(filename, mode = 'w', newline = '', encoding='utf-8') as csv_file:
        # write the contact node at the boundary of the network
        csv_writer = csv.writer(csv_file, delimiter=' ')
        csv_writer.writerow( ["Contact region in coninuum domain and characteristic time" ])
        csv_writer.writerow((Points[index_contact].x, Points[index_contact].y, M*T_END/N_T))
        csv_writer.writerow(['X', 'Y', 'Data_X', 'Data_Y'])
        csv_writer.writerows(data)



    with open(filename, mode = 'a', newline = '', encoding='utf-8') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=' ')
        csv_writer.writerows(data)

    print(f"Process {rank} saved data to {filename}")
    
           
    
    