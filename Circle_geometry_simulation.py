#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 10:19:36 2025

@author: phandangtoai
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 15:51:18 2025

@author: phandangtoai
"""

from mpi4py import MPI
import numpy as np
import csv
# import matplotlib.pyplot as plt
from scipy.spatial._qhull import Delaunay
import sys
sys.path.append("/Users/phandangtoai/Documents/Floe2D/Griffith-master_Dimitri/")
# from griffith.geometry import Polygon, Point
# from scipy.spatial import Voronoi, voronoi_plot_2d

# np.random.seed(1212)

lambda_intensity= 0.0025 
radius = 100.

area = np.pi * radius**2  # Area of the disk
num_points = np.random.poisson(lambda_intensity * area)  # Poisson-distributed number of points

print(num_points)
# Generate uniform random points in polar coordinates
r = np.sqrt(np.random.uniform(0, radius**2, num_points))  # Radius sampled with sqrt to maintain uniformity
theta = np.random.uniform(0, 2 * np.pi, num_points)  # Angle uniformly distributed

# Convert polar to Cartesian coordinates
x = r * np.cos(theta)
y = r * np.sin(theta)
x = np.append(x, radius)
y = np.append(y, 0.)

Points = np.column_stack((x, y))

tri = Delaunay(Points)  # triangulation step
Springs = set()
for triangle in tri.simplices:
    for index1 in triangle:
        for index2 in triangle:
            if index1 != index2:
                Springs.add((min(index1, index2), max(index1, index2)))



# poly.plot()
# plt.figure()
# circle = plt.Circle((0., 0.) , 1. , color = 'blue', fill=False)
# plt.plot(Points[:,0], Points[:,1],'o', color = 'orange')
# plt.gca().add_patch(circle)

# # plt.figure()
# V = Voronoi(Points)
# voronoi_plot_2d(V)


# plt.figure()
# plt.plot(Points[:,0], Points[:,1],'o', color = 'orange')
# plt.triplot(Points[:,0], Points[:,1], tri.simplices, color = 'b')
# plt.show()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()  # Process ID
size = comm.Get_size()  # Total number of processes

# Unique output filename per MPI process
filename = f"masses-springs_circle_{rank+1}.csv"

with open(filename, mode = 'w', newline = '', encoding='utf-8') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter=',')
    csv_writer.writerows(Points)

print(f"Process {rank} saved data to {filename}")








