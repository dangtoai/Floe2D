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
from numpy.linalg import norm
import matplotlib.tri as tri
# from matplotlib import animation
sys.path.append("/Users/phandangtoai/Documents/Floe2D/Griffith-master_Dimitri")
from griffith.geometry import Point, dist
from Func import Angle_mat, Traction_mat, Length_mat, Torsion_mat
from Func import Node, Floe, Energy_studies, N_T, timeout_handler, T_LIMIT
from Func import compute_dissipated_energy, compute_kinetic_energy

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
    
    # filename = f"masses-springs_limit_circle_{rank+1}.csv"
    filename = f"masses-springs_limit_circle_1.csv"

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
    V0 = np.array([-0.5 , 0.])
    # polygon = Polygon(Points)
    distances = [dist(point_input, Points[i])
                  for i in floe.border_nodes_index()]
    index_contact = floe.border_nodes_index()[distances.index(min(distances))]
    # # print("the nearest collision at point ",
    #       # Points[index_contact], " of the floe")
    Nodes[index_contact] = Node(Nodes[index_contact].position(), V0, id_number = index_contact)
    floe = Floe(nodes=Nodes, springs=None, mass= Mass,
                stiffness= traction_stiff, torsion_stiff= 0*torsion_stiff ,viscosity=100, id_number=1, impact_node= True)
    
    
    Angle_Mat = Angle_mat(floe)
    Traction_Mat = Traction_mat(floe, Angle_Mat)
    Length_Mat = Length_mat(floe)
    Torsion_Mat = Torsion_mat(floe)
    
    T_END = 2.4  # time end

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

    Energy = Energy_studies(All_positions_velocities, floe, Length_Mat, Angle_Mat, Torsion_Mat)
    
    ##### finding time T* := argmax of the energy elastic
    M = np.where(Energy[-1] == max(Energy[-1]))[0][0]
    print("characteristic collision time = ", M*T_END/N_T)

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
    
    # min_val= -0.16
    # max_val = 0.16
    # fig = plt.figure(figsize=(12, 6))

    #####First subplot (for u1)
    # ax1 = fig.add_subplot(121, projection='3d')
    # trisurf1 = ax1.plot_trisurf(X, Y, data_deformation[:, 0], 
    #                         cmap='bwr', edgecolor='none', vmin=min_val, vmax=max_val)
    # ax1.set_title("$u_1$")

    # #### Second subplot (for u2)
    # ax2 = fig.add_subplot(122, projection='3d')
    # trisurf2 = ax2.plot_trisurf(X, Y, data_deformation[:, 1], 
    #                         cmap='bwr', edgecolor='none', vmin=min_val, vmax=max_val)
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
    times_to_plot = [0.25, 0.8, 1.5, 2.3]
    steps = [int(t * N_T / T_END) for t in times_to_plot]

    # Compute all displacements
    displacements_u1 = []
    displacements_u2 = []

    for step in steps:
        data_deformation = np.zeros((floe.n, 2))
        for i in range(floe.n):
            x0 = All_positions_velocities[4*i][0]
            y0 = All_positions_velocities[4*i+1][0]
            x_now = All_positions_velocities[4*i][step]
            y_now = All_positions_velocities[4*i+1][step]
            data_deformation[i, 0] = x_now - x0  # u1
            data_deformation[i, 1] = y_now - y0  # u2

        displacements_u1.append(data_deformation[:, 0])
        displacements_u2.append(data_deformation[:, 1])

    # Calculate separate color limits
    u1_all = np.concatenate(displacements_u1)
    u2_all = np.concatenate(displacements_u2)

    max_abs_u1 = np.max(np.abs(u1_all))
    max_abs_u2 = np.max(np.abs(u2_all))

    vmin_u1, vmax_u1 = -max_abs_u1, max_abs_u1
    vmin_u2, vmax_u2 = -max_abs_u2, max_abs_u2

    # Create triangulation
    triang = tri.Triangulation(X, Y)

    # Create 4x2 figure (4 rows, 2 columns)
    fig, axes = plt.subplots(4, 2, figsize=(12, 16))

    for idx, (step, u1, u2) in enumerate(zip(steps, displacements_u1, displacements_u2)):
        # Left column: u1
        ax_left = axes[idx, 0]
        tcf1 = ax_left.tricontourf(triang, u1, cmap='RdBu_r', levels=50, 
                                  vmin=vmin_u1, vmax=vmax_u1)
        ax_left.tricontour(triang, u1, colors='k', linewidths=0.3, levels=10, alpha=0.3)
        ax_left.triplot(triang, 'k-', linewidth=0.2, alpha=0.2)
        ax_left.plot(X[-1], Y[-1], 'ro', markersize=6, markeredgecolor='white', markeredgewidth=1)
        ax_left.set_xlabel('X (m)')
        ax_left.set_ylabel('Y (m)')
        ax_left.set_title(f't = {times_to_plot[idx]:.2f} s')
        ax_left.set_aspect('equal')

        # Right column: u2
        ax_right = axes[idx, 1]
        tcf2 = ax_right.tricontourf(triang, u2, cmap='RdBu_r', levels=50, 
                                   vmin=vmin_u2, vmax=vmax_u2)
        ax_right.tricontour(triang, u2, colors='k', linewidths=0.3, levels=10, alpha=0.3)
        ax_right.triplot(triang, 'k-', linewidth=0.2, alpha=0.2)
        ax_right.plot(X[-1], Y[-1], 'ro', markersize=6, markeredgecolor='white', markeredgewidth=1)
        ax_right.set_xlabel('X (m)')
        ax_right.set_ylabel('Y (m)')
        ax_right.set_title(f't = {times_to_plot[idx]:.2f} s')
        ax_right.set_aspect('equal')

    # Add column labels
    axes[0, 0].text(0.5, 1.1, r'$u_1$ ', transform=axes[0, 0].transAxes, 
                    fontsize=14, ha='center', va='bottom', fontweight='bold')
    axes[0, 1].text(0.5, 1.1, r'$u_2$ ', transform=axes[0, 1].transAxes, 
                    fontsize=14, ha='center', va='bottom', fontweight='bold')

    # Add separate colorbars - LEFT for u1, RIGHT for u2
    fig.subplots_adjust(left=0.08, right=0.9, wspace=0.3)

    # Colorbar for u1 (LEFT side)
    cbar_ax1 = fig.add_axes([0.05, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
    cbar1 = fig.colorbar(tcf1, cax=cbar_ax1)
    cbar1.set_label(r'$u_1$ (m)', fontsize=12)

    # Colorbar for u2 (RIGHT side)  
    cbar_ax2 = fig.add_axes([0.95, 0.15, 0.02, 0.7])
    cbar2 = fig.colorbar(tcf2, cax=cbar_ax2)
    cbar2.set_label(r'$u_2$ (m)', fontsize=12)

    plt.tight_layout()
    fig.savefig("displacement_comparison_4x2.png", dpi=300, bbox_inches='tight')
    # plt.show()
    
    
    D = compute_dissipated_energy(All_positions_velocities, floe, T_END, floe.mu)
    K = compute_kinetic_energy(All_positions_velocities, floe)
    plt.figure()
    plt.plot(D)
    plt.plot(K)
    plt.plot(Energy[-1])
    plt.plot(D+K+Energy[-1])

