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

import numpy as np
# import scipy as sp
import matplotlib.pyplot as plt
from scipy import stats
from scipy.spatial._qhull import Delaunay
import sys
sys.path.append("/Users/phandangtoai/Documents/Floe2D/Griffith-master_Dimitri/")
from griffith.geometry import Polygon, Point
from scipy.spatial import Voronoi, voronoi_plot_2d

# np.random.seed(1212)

lambda_intensity= 15. 
radius = 1.

area = np.pi * radius**2  # Area of the disk
num_points = np.random.poisson(lambda_intensity * area)  # Poisson-distributed number of points

# Generate uniform random points in polar coordinates
r = np.sqrt(np.random.uniform(0, radius**2, num_points))  # Radius sampled with sqrt to maintain uniformity
theta = np.random.uniform(0, 2 * np.pi, num_points)  # Angle uniformly distributed

# Convert polar to Cartesian coordinates
x = r * np.cos(theta)
y = r * np.sin(theta)
x = np.append(x, 1.)
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

import csv

with open('masses-springs_circle.csv', mode = 'w', newline = '', encoding='utf-8') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter=',')
    csv_writer.writerows(Points)










