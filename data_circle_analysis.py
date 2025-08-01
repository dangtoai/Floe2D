#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 11:06:57 2025

@author: phandangtoai
"""

import numpy as np
from numpy import pi,cos,sin
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit
from numpy.linalg import norm
sys.path.append("/Users/phandangtoai/Documents/Floe2D/mean displacement mass evolution/data circle network v0 = (-1,0) m=150/")


data_list = []
# with open("data_circle_2.txt", "r") as file:
with open("data_circle_limit_RBF_1.txt", "r") as file:
    
# with open("data_circle_geometry3_1.txt", "r") as file:
    for line in file:
        values = list(map(float, line.split()))  # Convert strings to floats
        data_list.append(np.array(values))  # Store as NumPy array

# Convert list of arrays into a NumPy array if needed
data_array = np.array(data_list)

# Print the data
radius = 100
# Number of points (41 points)
num_points = 41

# Compute mean for each point
mean_values = np.zeros((num_points, data_array.shape[1]))  # Initialize mean array

for i in range(num_points):
    mean_values[i] = np.mean(data_array[i::num_points], axis=0)  # Take mean of every num_points row

# Print result
# print("Mean values for each point:")
# print(mean_values)

n = int(len(data_array)/num_points)
print("total simulation = " , n)

#convergence ponctuel 
val_x0 = data_array[20::num_points]
mean_x0 = np.zeros_like(val_x0)
for i in range(1,n): 
    mean_x0[i] = abs(np.mean(val_x0[:i], axis=0) - mean_values[20])
    
print(" mean Ud(x0)= ",mean_values[20])
# plt.figure()
# plt.plot(np.arange(1, n), mean_x0[:,0][1:])
    


# convergence of sequenc Sn = 1/n sum_1^n 41 valeurs 
# for i in range(num_points):
#     mean_values[i] = np.mean(data_array[i::num_points], axis=0)  # Take mean of every num_points row



# # mean value of boundary condition simulation

plt.figure()
plt.plot(np.arange(0, num_points), mean_values[:, 0], '-')


plt.figure()
plt.plot(np.arange(0, num_points), mean_values[:, 1],'-')


# mean and each simulation


plt.figure()
for i in range(n):
    plt.plot(np.arange(0, num_points), data_array[i*num_points: i*num_points+41][:,0],'x')
# plt.plot(np.arange(0, num_points), mean_values[:, 0])


plt.figure()
for i in range(n):
    plt.plot(np.arange(0, num_points), data_array[i*num_points: i*num_points+41][:,1],'x')
# plt.plot(np.arange(0, num_points), mean_values[:, 1])

###plot max displacement

# K_eff = 6172539.645169568

# Masses = [150, 1e3, 1e4, 1e5, 1e6]
# Sol = [-np.sqrt((288513.6110439606+m)/K_eff) for m in Masses]

# data = np.append(data_array[20::41][:,0], [-0.34134081])

# plt.figure()
# plt.plot(np.arange(0,5), data, label =  "valeur numerique")
# plt.plot(np.arange(0,5), Sol, label = "estimation")
# plt.ylim(-0.5,-0)
# plt.legend()

#plot on demi circle

new_xdata = [radius * cos(theta) for theta in np.linspace(-pi/2, pi/2, 41) ]
new_ydata = [radius * sin(theta) for theta in np.linspace(-pi/2, pi/2, 41) ]


plt.figure()
plt.quiver(new_xdata, new_ydata, mean_values[:,0], mean_values[:,1], angles='xy', scale_units='xy', scale=0.01)



# initial_guess = (-1, 0)
# # # bounds = ([0, -np.inf], [np.inf, np.inf])

# popt, pcov = curve_fit(Model, distances, z_data, p0=initial_guess)

# A_opt, B_opt = popt
# print(f"Estimated parameters: A = {A_opt:.4f}, B = {B_opt:.4f}")

# plt.figure()
# plt.plot(distances, z_data, label = 'data')
# plt.plot(distances, Model(distances, A_opt, B_opt), label = 'fit')
# plt.legend()
# plt.show()


# z_model = Model((x_data, y_data), A_opt, B_opt)

# plt.figure()
# plt.plot(np.arange(0, num_points), mean_values[:, 0], '-')
# plt.plot(np.arange(0, num_points), )


# plt.figure()
# plt.plot()


# x_fit = np.linspace(min(x_data), max(x_data), 100)
# y_fit = np.linspace(min(y_data), max(y_data), 100)
# X_fit, Y_fit = np.meshgrid(x_fit, y_fit)
# Z_fit = Model((X_fit, Y_fit), A_opt, B_opt)

# # Plot the data and fitted surface
# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(x_data, y_data, z_data, color='red', label='Data')
# # ax.plot_surface(X_fit, Y_fit, Z_fit, cmap='viridis', alpha=0.5)
# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("z")
# plt.legend()
# plt.show()





