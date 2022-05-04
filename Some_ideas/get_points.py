# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.animation as animation
from scipy.spatial import Voronoi, voronoi_plot_2d, Delaunay 
from Func import *
from graph import *
import scipy
from scipy import stats


np.random.seed(2)
xMin=0;xMax=1
yMin=0;yMax=1
xDelta=xMax-xMin
yDelta=yMax-yMin; #rectangle dimensions
areaTotal=xDelta*yDelta;

#Point process parameters
lambda0 = 8                                                      #intensity 

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
print("number of nodes = ", numbPoints)
nb_nodes = len(tri.points)


#plt.figure()
# for triangle in tri.simplices:
#     for index1 in triangle:
#         for index2 in triangle:
#             if index1 != index2 :
#                 plt.plot([tri.points[index1][0], tri.points[index2][0]], 
#                          [tri.points[index1][1], tri.points[index2][1]])
#                 # plt.text(tri.points[index1][0], tri.points[index1][1], str(index1), color = "red")

possible = []
for triangle in tri.simplices:
    for index1 in triangle:
        for index2 in triangle:
            if index1 != index2:
                possible.append((min(index1, index2), max(index1, index2)))
#print(possible)
# print("set of all edge=", set(possible))


# print("all points", tri.points)
print("all edges", set(possible))

# Points = np.array(Points)
Springs = set(possible)

#calculer les angles initiale
# idee pour chercher la fracture: chercher les simplexes au bord tel que leurs degres egale au moins a 3
# commencer uniquement par les arretes au bords 
# chercher les arretes dans les simplexes voisins, 
# chercher les arretes voisins 

Points = tri.points
# Springs = {(0,1),(1,2),(2,3),(3,4),(4,5),(0,5),(1,5),(0,5),(2,5),(2,4)}


Nodes = []
V0 = np.array([0.55, 0.])
for i in range(len(Points)):
    Nodes.append(Node(Points[i], V0, i))

k = 100.
floe = Floe(nodes=Nodes, springs=Springs,
            stiffness= k, viscosity=k/10., id_number=1)



c = floe.border()









