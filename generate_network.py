"""
Created on Thu May  4 16:15:23 2023

@author: phandangtoai
"""
import csv
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.spatial._qhull import Delaunay
# from scipy.spatial import Voronoi, voronoi_plot_2d
from mpi4py import MPI


from scipy import stats
# sys.path.append("/Users/phandangtoai/Documents/Floe2D/Griffith-master_Dimitri/")
from griffith.geometry import Polygon, Point

comm = MPI.COMM_WORLD
rank = comm.Get_rank()  # Process ID
size = comm.Get_size()  # Total number of processes
# np.random.seed(1212)


### set of points is in the counter clock sens

# print("alles gut")
# Input = input("Select a geometry to run the random process: ")
# Input = int(Input)
# print(f"the geometry number {Input}")

Input = 45
mat = sp.io.loadmat("Biblio_Floes.mat")['G'][Input][0] # 100 geometries of ice floes

poly = Polygon( [ Point(mat[i][0], mat[i][1]) for i in range(len(mat) - 1) ])
# print("area = ", poly.area())

# random points process:

xMin = min(mat[:,0])
xMax = max(mat[:,0])
yMin = min(mat[:,1])
yMax = max(mat[:,1])

xDelta = xMax - xMin
yDelta = yMax - yMin
# print(xMin, xMax, yMin, yMax)

LAMBDA = 0.002 # intensity of random process
num_points = int(LAMBDA * (poly.area()  )) # number of random points inside of the polygon
print(num_points)
# x coordinates of Poisson points
xx = xDelta * stats.uniform.rvs(0, 1, ((num_points, 1))) + xMin
# y coordinates of Poisson points
yy = yDelta * stats.uniform.rvs(0, 1, ((num_points, 1))) + yMin




Points = []
OutPoints = []
for i in range(num_points):
    # select interior point 
    if poly.has_point(Point(xx[i][0], yy[i][0])):
        Points.append( np.array([xx[i][0], yy[i][0]]))
    else: OutPoints.append(np.array([xx[i][0], yy[i][0]]))

# P = np.array([1010.5, 106.0])
# Points.append(P)

Points = np.array(Points)



# Points = np.append(Points, P)

# print(Points[-1])
OutPoints = np.array(OutPoints)

tri = Delaunay(Points)  # triangulation step
Springs = set()
for triangle in tri.simplices:
    for index1 in triangle:
        for index2 in triangle:
            if index1 != index2:
                Springs.add((min(index1, index2), max(index1, index2)))


filename = f"masses-springs_geometry3_{rank+1}.csv"

with open(filename, mode = 'w', newline = '', encoding='utf-8') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter=',')
    csv_writer.writerows(Points)

print(f"Process {rank} saved data to {filename}")

# with open('masses-springs.csv', mode = 'w', newline = '', encoding='utf-8') as csv_file:
#     csv_writer = csv.writer(csv_file, delimiter=',')
#     csv_writer.writerows(Points)

# poly.plot()
plt.figure()
# poly.plot()
plt.plot(Points[:,0], Points[:,1],'o', color = 'orange')
plt.plot(OutPoints[:,0], OutPoints[:,1], 'o', color = 'green')

# fig, ax = plt.subplots()
# V = Voronoi(Points)
# fig, ax = poly.plot()
# voronoi_plot_2d(V, ax = ax)
# ax.set_xlim(xMin-1, xMax+1)
# ax.set_ylim(yMin-1, yMax+1)

plt.figure()
poly.plot()
plt.triplot(Points[:,0], Points[:,1], tri.simplices, color = 'b')
plt.show()

# sys.exit()
