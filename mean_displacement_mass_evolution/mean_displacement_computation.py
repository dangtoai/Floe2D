#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 10:54:02 2025

@author: phandangtoai
"""

import numpy as np
from numpy import pi,cos,sin
import matplotlib.pyplot as plt
# import sys
import os

directories = [
    "data circle network v0 = (-1,0) m=1e2",
    # "data circle network v0 = (-1,0) m=150",
    "data circle network v0 = (-1,0) m=1e3",
    "data circle network v0 = (-1,0) m=1e4",
    "data circle network v0 = (-1,0) m=1e5",
    "data circle network v0 = (-1,0) m=1e6",
    "data circle network v0 = (-1,0) m=1e7",
    # "data circle network v0 = (-1,0) m=1e8",
]
num_points = 41

radius = 100.
x_data = np.array([radius * cos(theta) for theta in np.linspace(0, pi/2, 21) ])
y_data = np.array([radius * sin(theta) for theta in np.linspace(0, pi/2, 21) ])
x0, y0 = 100, 0  # Adjust based on your setup
v0_x, v0_y = -0.5, 0.  # Initial velocity
theta_0 = np.arctan2(v0_y, v0_x)  # This aligns with impact direction
r_data = np.sqrt((x_data - x0) ** 2 + (y_data - y0) ** 2)
theta_data = np.arctan2(y_data - y0, x_data - x0)


# Filename to read in each directory
filename = "data_circle_limit_RBF_1.txt"

all_data = []

for directory in directories:
    file_path = os.path.join(directory, filename)
    # file_path = directory
    # print(file_path)
    if not os.path.exists(file_path):
        print(f"Warning: File not found in {directory}")
        # break
        # continue  # Skip missing files

    # # Read the file and convert it into a NumPy array
    data_list = []
    with open(file_path, "r") as file:
        for line in file:
            values = list(map(float, line.split()))  # Convert strings to floats
            data_list.append(values)

    data_array = np.array(data_list)  # Convert list to NumPy array
    all_data.append(data_array)
    
    #### if plot every simulation results
    # n = int(len(data_array)/num_points)
    # print("total simulation = " , n)
    # plt.figure()
    # for i in range(n):
    #     plt.plot(np.arange(0, num_points), data_array[i*num_points: i*num_points+41][:,0],'-')
    # # plt.plot(np.arange(0, num_points), mean_values[:, 0])


    # plt.figure()
    # for i in range(n):
    #     plt.plot(np.arange(0, num_points), data_array[i*num_points: i*num_points+41][:,1],'x')
    # # plt.plot(np.arange(0, num_points), mean_values[:, 1])


num_points = 41

# n = int(len(data_array)/num_points)
# print("total simulation = " , n)
# plt.figure()
# # for i in range(n):
# #     plt.plot(np.arange(0, num_points), data_array[i*num_points: i*num_points+41][:,0],'x')
# plt.plot(np.arange(0, num_points), mean_values[:, 0])


# plt.figure()
# # for i in range(n):
# #     plt.plot(np.arange(0, num_points), data_array[i*num_points: i*num_points+41][:,1],'x')
# plt.plot(np.arange(0, num_points), mean_values[:, 1])


print("number of simulation", [len(all_data[i])/num_points for i in range(len(directories))])

mean_all_data = []
for data in all_data:
    num_rows = data.shape[0]
    
    if num_rows % num_points != 0:
        print(f"Warning: Data size {num_rows} is not a multiple of {num_points}, some data may be ignored.")
    
    # Compute mean every 41 rows
    mean_values = np.array([np.mean(data[i::num_points], axis=0) for i in range(num_points)])
    
    # Store the mean results for each dataset
    mean_all_data.append(mean_values)

# Convert list to NumPy array for easier handling
mean_all_data = np.array(mean_all_data)  # Shape: (5, 41, num_features)


# plt.figure()
# for i in range(len(directories)):
#     # plt.plot(np.arange(41), mean_all_data[i][:,0], label = f'$P=10^{i+2}(kg)$')

#     plt.plot(r_data, mean_all_data[i][:,0][20:], label = f'$P=10^{i+2}(kg)$')
# plt.ylabel('$u_x$')
# plt.xlabel("Distance to $x_0$ (m)")
# plt.legend()

# print([mean_all_data[i,20] for i in range(len(directories))])

# plt.figure()
# for i in range(len(directories)):
#     # plt.plot(np.arange(41), mean_all_data[i][:,1], label = f'$P=10^{i+2}(kg)$')
#     plt.plot(r_data, mean_all_data[i][:,1][20:], label = f'$P=10^{i+2}(kg)$')
# plt.ylabel('$u_x$')
# plt.xlabel("Distance to $x_0$ (m)")
# plt.legend()


fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

# Plot u_x
for i in range(len(directories)):
    axs[0].plot(r_data, mean_all_data[i][:,0][20:], label=fr'$P = 10^{i+2}~\mathrm{{kg}}$')
axs[0].set_ylabel(r'$u_x$')
axs[0].legend()
axs[0].grid(True, linestyle='--', alpha=0.3)

# Plot u_y
for i in range(len(directories)):
    axs[1].plot(r_data, mean_all_data[i][:,1][20:], label=fr'$P = 10^{i+2}~\mathrm{{kg}}$')
axs[1].set_ylabel(r'$u_y$')
axs[1].set_xlabel(r'Distance to $x_0$ (m)')
axs[1].legend()
axs[1].grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()

max_displacement = np.array([data[20] for data in mean_all_data])


def u1_model(r, theta, A, B):
    return A * np.exp(-B * r) #* np.cos(theta + theta_0)


from scipy.optimize import curve_fit

for i in range(len(directories)):
    
    z1_data, z2_data = mean_all_data[i][:,0][20:], mean_all_data[i][:,1][20:]
    theta_flat = theta_data.flatten()
    r_flat = r_data.flatten()
    z1_flat = z1_data.flatten()
    z2_flat = z2_data.flatten()
    
    # Fit u1
    params_u1, _ = curve_fit(lambda r, A, B: u1_model(r, theta_flat, A, B), r_flat, z1_flat, p0=[-1, 1])
    # params_u2, _ = curve_fit(lambda r, A, B: u1_model(r, theta_flat, A, B), r_flat, z2_flat, p0=[-1, 1])

    
    A_opt, B_opt = params_u1
    # AA_opt, BB_opt = params_u2
    print(f"Optimal A: {A_opt}, B: {B_opt}")
    # print(f"Optimal A: {AA_opt}, B: {BB_opt}")
    
    plt.figure()
    plt.plot(r_flat, z1_flat, label= "Mean displacement field")
    plt.plot(r_flat, u1_model(r_flat, theta_flat, A_opt, B_opt),'--', label = "Curve fit")
    plt.legend()
    break 


# plt.figure()
# for i in range(len(all_data)):
#     plt.plot(t, max_displacement[:,0], 'x')
#     # plt.plot(t, S)
# # plt.ylabel('$-|u|_\infty (m)$')
# plt.xlabel('$P(kg)$')

# plt.figure()
# for i in range(5):
#     plt.plot(np.arange(5), max_displacement[:,0], 'x')
#     plt.plot(np.arange(5), S)



