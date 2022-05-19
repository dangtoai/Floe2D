# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.animation as animation
from scipy.spatial import Voronoi, voronoi_plot_2d, Delaunay
from Func import *
from graph import *
import scipy
from scipy import stats


np.random.seed(4)
xMin = 0
xMax = 1
yMin = 0
yMax = 1
xDelta = xMax-xMin
yDelta = yMax-yMin  # rectangle dimensions
areaTotal = xDelta*yDelta

# Point process parameters
lambda0 = 7  # intensity

# Simulate Poisson point process
numbPoints = scipy.stats.poisson(
    lambda0*areaTotal).rvs()  # Poisson number of points
# x coordinates of Poisson points
xx = xDelta*scipy.stats.uniform.rvs(0, 1, ((numbPoints, 1))) + xMin
# y coordinates of Poisson points
yy = yDelta*scipy.stats.uniform.rvs(0, 1, ((numbPoints, 1))) + yMin
Points = []  # points coordinates
for i in range(len(xx)):
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
for i in range(len(Points)):
    Nodes.append(Node(Points[i], V0, i))

k = 100.
floe = Floe(nodes=Nodes, springs=Springs,
            stiffness=k, viscosity=k/10., id_number=1)
floe.plot_init()

print(floe.fractures_admissible())

# nx.draw(G, with_labels=True, font_weight='bold')


                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
