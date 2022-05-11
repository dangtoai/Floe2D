# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.animation as animation
from scipy.spatial import Voronoi, voronoi_plot_2d, Delaunay 
from Func import *
from graph import *
import scipy
from scipy import stats


np.random.seed(7)
xMin=0;xMax=1
yMin=0;yMax=1
xDelta=xMax-xMin
yDelta=yMax-yMin; #rectangle dimensions
areaTotal=xDelta*yDelta;

#Point process parameters
lambda0 = 5                                                      #intensity 

#Simulate Poisson point process
numbPoints = scipy.stats.poisson( lambda0*areaTotal ).rvs()       #Poisson number of points
xx = xDelta*scipy.stats.uniform.rvs(0,1,((numbPoints,1))) + xMin  #x coordinates of Poisson points
yy = yDelta*scipy.stats.uniform.rvs(0,1,((numbPoints,1))) + yMin  #y coordinates of Poisson points
Points = []                                                       #points coordinates
for i in range(len(xx)):
    Points.append([xx[i][0],yy[i][0]])
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
            stiffness= k, viscosity=k/10., id_number=1)
floe.plot_init()

# print(floe.fractures_admissible())

# # floe.fracture_admissible()

G = nx.Graph()
BE = floe.border_edges_index()
l = [i for i in range(floe.n)]
l2 = [list(e) for e in floe.simplices()]
G.add_nodes_from(l)
G.add_edges_from(floe.springs)

# nx.draw(G, with_labels=True, font_weight='bold')

# floe.plot_init()
FA = []
for i,j in combinations(floe.border_nodes_index(), 2):
    for path in nx.all_simple_edge_paths(G, i, j, cutoff=(floe.n-2)): 
        path = np.sort(np.array(path))
        if len(path) in range(3, int(floe.n - 1)):
            # or len(path) == floe.n-2 
            if ((tuple(path[0]) in BE) == True and (tuple(path[-1]) in BE) == True ):
                # FA.append(path)
                if np.all([list(set(np.append(path[i], path[i+1]))) in l2 for i in range(len(path)-1)]):
                    FA.append(path)
                
# # print([list(set(np.append(FA[1][i], FA[1][i+1]))) in l2 for i in range(len(FA[1])-1)])








