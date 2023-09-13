#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:15:23 2023

@author: phandangtoai
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
# from scipy.spatial import Voronoi, voronoi_plot_2d, Delaunay
import scipy
from Func import *


sys.path.append("/Users/phandangtoai/Documents/Floe2D/")

# from scipy import stats


np.random.seed(5)


# # unit square poisson point process
XMIN = 0
XMAX = 1
YMIN = 0
YMAX = 1
XDELTA = XMAX-XMIN
YDELTA = YMAX-YMIN  # rectangle dimensions
TOTAL_AREA = XDELTA*YDELTA

# Point process parameters
LAMBDA = 88  # intensity

# Simulate Poisson point process:
numbPoints = scipy.stats.poisson( LAMBDA*TOTAL_AREA).rvs()  # Poisson number of points
# x coordinates of Poisson points
xx = XDELTA*scipy.stats.uniform.rvs(0, 1, ((numbPoints, 1))) + XMIN
# y coordinates of Poisson points
yy = YDELTA*scipy.stats.uniform.rvs(0, 1, ((numbPoints, 1))) + YMIN

Points = []  # points coordinates

for i in range(numbPoints):
    Points.append([xx[i][0], yy[i][0]])
np.array(Points)
# vor = Voronoi(Points)
tri = Delaunay(Points)
# print("number of nodes = ", numbPoints)
nb_nodes = len(tri.points)

possible = []
for triangle in tri.simplices:
    for index1 in triangle:
        for index2 in triangle:
            if index1 != index2:
                possible.append((min(index1, index2), max(index1, index2)))


# Points = np.array(Points)
Springs = set(possible)

Points = tri.points

Nodes = []
V0 = np.array([0.55, 0.])
for i in range(numbPoints):
    Nodes.append(Node(Points[i], V0, i))

k = 100.
floe = Floe(nodes=Nodes, springs=Springs,
            stiffness=k, viscosity=k/10., id_number=1)
floe.plot_init()
plt.show()
# print(floe.fractures_admissible())
# nx.draw(G, with_labels=True, font_weight='bold')
