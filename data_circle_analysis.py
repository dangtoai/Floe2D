#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 11:06:57 2025

@author: phandangtoai
"""

import numpy as np
from numpy import pi,cos,sin
import matplotlib.pyplot as plt


data_list = []
# with open("data_circle_2.txt", "r") as file:
with open("data_circle_RBF.txt", "r") as file:
    for line in file:
        values = list(map(float, line.split()))  # Convert strings to floats
        data_list.append(np.array(values))  # Store as NumPy array

# Convert list of arrays into a NumPy array if needed
data_array = np.array(data_list)

# Print the data
# print(data_array)
radius = 1
# Number of points (41 points)
num_points = 41

# Compute mean for each point
mean_values = np.zeros((num_points, data_array.shape[1]))  # Initialize mean array

for i in range(num_points):
    mean_values[i] = np.mean(data_array[i::num_points], axis=0)  # Take mean of every num_points row

# Print result
print("Mean values for each point:")
print(mean_values)


# mean value of boundary condition simulation

plt.figure()
plt.plot(np.arange(0, num_points), mean_values[:, 0], '-')


plt.figure()
plt.plot(np.arange(0, num_points), mean_values[:, 1],'-')


# mean and each simulation


plt.figure()
for i in range(100):
    plt.plot(np.arange(0, num_points), data_array[i*num_points: i*num_points+41][:,0],'--')
plt.plot(np.arange(0, num_points), mean_values[:, 0])


plt.figure()
for i in range(100):
    plt.plot(np.arange(0, num_points), data_array[i*num_points: i*num_points+41][:,1],'--')
plt.plot(np.arange(0, num_points), mean_values[:, 1])

#plot on demi circle

new_xdata = [radius * cos(theta) for theta in np.linspace(-pi/2, pi/2, 41) ]
new_ydata = [radius * sin(theta) for theta in np.linspace(-pi/2, pi/2, 41) ]


plt.figure()
plt.quiver(new_xdata, new_ydata, mean_values[:,0], mean_values[:,1], angles='xy', scale_units='xy', scale=0.5)





