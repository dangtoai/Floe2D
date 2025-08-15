#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 11:35:29 2025

@author: phandangtoai
"""

import numpy as np
from numpy import pi,cos,sin
import matplotlib.pyplot as plt
# import sys
import os
from scipy.optimize import curve_fit


directories = [
    "data circle network v0 = (-0.5,0) m=1e5",
    "data circle network v0 = (-1,0) m=1e5",
    "data circle network v0 = (-1.5,0) m=1e5",
    "data circle network v0 = (-2,0) m=1e5",
    "data circle network v0 = (-2.5,0) m=1e5",
]

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
        continue  # Skip missing files

    # # Read the file and convert it into a NumPy array
    data_list = []
    with open(file_path, "r") as file:
        for line in file:
            values = list(map(float, line.split()))  # Convert strings to floats
            data_list.append(values)

    data_array = np.array(data_list)  # Convert list to NumPy array
    all_data.append(data_array)

num_points = 41


radius = 100.
x_data = np.array([radius * cos(theta) for theta in np.linspace(0, pi/2, 21) ])
y_data = np.array([radius * sin(theta) for theta in np.linspace(0, pi/2, 21) ])
x0, y0 = 100, 0  # Adjust based on your setup
v0_x, v0_y = -0.5, 0.  # Initial velocity
theta_0 = np.arctan2(v0_y, v0_x)  # This aligns with impact direction
r_data = np.sqrt((x_data - x0) ** 2 + (y_data - y0) ** 2)
theta_data = np.arctan2(y_data - y0, x_data - x0)



print("number of simulations, ", [len(all_data[i])/num_points for i in range(len(all_data))])
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
# for i in range(len(all_data)):
#     plt.plot(r_data, mean_all_data[i][:,0][20:], label = f'$v_0=(-{0.5*(i+1)},0)$')
# plt.ylabel('$u_x$')
# plt.xlabel("Distance to $x_0$ (m)")
# plt.legend()

# plt.figure()
# for i in range(len(all_data)):
#     plt.plot(r_data, mean_all_data[i][:,1][20:], label = f'$v_0=(-{0.5*(i+1)},0)$')
# plt.ylabel('$u_y$')
# plt.xlabel("Distance to $x_0$ (m)")
# plt.legend()

fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)  # 2 rows, 1 column

# Plot u_x (first component)
for i in range(len(all_data)):
    axs[0].plot(r_data, mean_all_data[i][:, 0][20:], label=fr'$v_0 = (-{0.5*(i+1)}, 0)$')
axs[0].set_ylabel(r'$u_x$')
axs[0].legend()
axs[0].grid(True, linestyle='--', alpha=0.3)

# Plot u_y (second component)
for i in range(len(all_data)):
    axs[1].plot(r_data, mean_all_data[i][:, 1][20:], label=fr'$v_0 = (-{0.5*(i+1)}, 0)$')
axs[1].set_ylabel(r'$u_y$')
axs[1].set_xlabel(r'Distance to $x_0$ (m)')
axs[1].legend()
axs[1].grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()


max_displacement = np.array([data[20] for data in mean_all_data])


def u1_model(r, A, B):
    return A * np.exp(-B * r) #* np.cos(theta + theta_0)


for i in range(len(directories)):
    
    z1_data, z2_data = mean_all_data[i][:,0][20:], mean_all_data[i][:,1][20:]
    theta_flat = theta_data.flatten()
    r_flat = r_data.flatten()
    z1_flat = z1_data.flatten()
    z2_flat = z2_data.flatten()
    
    # Fit u1
    params_u1, _ = curve_fit(lambda r, A, B: u1_model(r, A, B), r_flat, z1_flat, p0=[-1, 1])
    
    A_opt, B_opt = params_u1
    print(f"Optimal A: {A_opt}, B: {B_opt}")



# r_fit = np.linspace(min(r_flat), max(r_flat), 100)
# u1_fit = u1_model(r_fit, 0, A_opt, B_opt)

# plt.figure()
# plt.scatter(r_flat, z1_flat, label="Data u1", s=10)
# plt.plot(r_fit, u1_fit,'--' , label="Fit u1", color="r")
# plt.legend()
# plt.xlabel("r")
# plt.ylabel("u1")
# plt.title("Fit for u1")
# plt.show()

fig, axs = plt.subplots(2, 1, figsize=(8, 6))
fitted_params = []

def null_func(t):
    return 0*t

for i in range(len(mean_all_data)):
    vx_label = fr'$v_0 = (-{0.5*(i+1)}, 0)$'

    r_fit = r_data#[20:]
    u_x_fit = mean_all_data[i][:, 0][20:]
    u_y_fit = mean_all_data[i][:, 1][20:]

    # Fit u_x
    params_u1, _ = curve_fit(u1_model, r_fit, u_x_fit, p0=[-1, 1])
    A_opt, B_opt = params_u1
    fitted_params.append((0.5 * (i+1), A_opt, B_opt))

    # Plot u_x data and fit
    axs[0].plot(r_fit, u_x_fit, label=rf'{vx_label} (data)')
    axs[0].plot(r_fit, u1_model(r_fit, A_opt, B_opt), '--', label=rf'{vx_label} (fit)')
    
    # Plot u_y only (no fit here)
    axs[1].plot(r_fit, u_y_fit, label=rf'{vx_label}')
    axs[1].plot(r_fit, null_func(r_fit),'--', label=rf'{vx_label}')

axs[0].set_ylabel(r'$u_x$')
axs[0].legend()

axs[0].grid(True, linestyle='--', alpha=0.3)

axs[1].set_ylabel(r'$u_y$')
axs[1].set_xlabel(r'Distance to $x_0$ (m)')
axs[1].legend()
axs[1].grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()








