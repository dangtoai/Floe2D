# -*- coding: utf-8 -*-

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from numpy.linalg import norm
from scipy.integrate import solve_ivp
import matplotlib.animation as animation
from scipy.spatial import Voronoi, voronoi_plot_2d, Delaunay 
import os
import imageio
# from Func import *
from graph import *
import random
# from test import *

# Constant parameters
T_f      = 4.       #final time simulation
dt = 0.005          #time's step
N = 500

m  = 6.2
mu = 10.
k  = 100.3

#Simulation window parameters
np.random.seed(6)
xMin=3;xMax=4
yMin=0;yMax=1
xDelta=xMax-xMin
yDelta=yMax-yMin; #rectangle dimensions
areaTotal=xDelta*yDelta;

#Point process parameters
lambda0 = 9                                                      #intensity 

#Simulate Poisson point process
numbPoints = scipy.stats.poisson( lambda0*areaTotal ).rvs()       #Poisson number of points
xx = xDelta*scipy.stats.uniform.rvs(0,1,((numbPoints,1))) + xMin  #x coordinates of Poisson points
yy = yDelta*scipy.stats.uniform.rvs(0,1,((numbPoints,1))) + yMin  #y coordinates of Poisson points
Point_2d = []                                                     #points coordinates
for i in range(len(xx)):
    Point_2d.append([xx[i][0],yy[i][0]])
np.array(Point_2d)
# vor = Voronoi(Point_2d)
tri = Delaunay(Point_2d)
print("number of nodes = ", numbPoints)
nb_nodes = len(tri.points)
# voronoi_plot_2d(vor)  #Voronoi
# plt.scatter(xx,yy, edgecolor='b',  alpha=0.5 )
# plt.xlim(-0,1.)
# plt.ylim(-0.,1.)

# plt.triplot(xx[:,0], yy[:,0]) #Delaunay
# plt.xlabel("x")
# plt.ylabel("y")
# plt.tight_layout()
# plt.show()
# plt.savefig("Carre")

#plt.figure()
for triangle in tri.simplices:
    for index1 in triangle:
        for index2 in triangle:
            if index1 != index2 :
                plt.plot([tri.points[index1][0], tri.points[index2][0]], 
                         [tri.points[index1][1], tri.points[index2][1]])
                plt.text(tri.points[index1][0], tri.points[index1][1], str(index1), color = "red")

#contact matrix
Contact_Mat = np.zeros(nb_nodes*nb_nodes)
# I = np.eye(nb_nodes)
Contact_Mat = Contact_Mat.reshape(nb_nodes,nb_nodes)
# print(Contact_Mat)
for triangle in tri.simplices:
    for index1 in triangle:
        for index2 in triangle:
            if index1 != index2 :
                Contact_Mat[index1, index2] = 1
# print(Contact_Mat)

#initial position
Q0 =  tri.points
# length matrix
L = np.zeros((nb_nodes, nb_nodes))
#print(L)
for i in range(nb_nodes):
    for j in range(nb_nodes):
        L[i,j] = norm(Q0[j]-Q0[i])
        L[j,i] = L[i,j]

#initial velocity
# V0 = np.zeros((nb_nodes, 2))+0.5
V0 = np.zeros((nb_nodes, 2))
V0[3]=np.array([.2, .2]) #seed(3), lambda = 5
# V0[5]=np.array([1., 1.]) #seed(3), lambda = 5

# V0[1]=np.array([.5, .5]) #seed(4), lambda = 9
# V0[7]=np.array([0.5, .5])
# V0[5]=np.array([.5, .5]) #seed(3), lambda = 9

def Unit_vect(vect1, vect2):
    if (vect1[0] == vect2[0] and vect1[1] == vect2[1] ): return 0.
    else : return (vect2-vect1)/norm(vect2-vect1)

def System(t, Y):
    u = np.zeros((nb_nodes, nb_nodes, 2))
    Q = np.reshape(Y, (nb_nodes*2, 2))
    Y_ = np.zeros_like(Q)
    # if node0 is stable, its velocity and acceleration = 0
    # begin at 1
    for i in range(0, nb_nodes):
        Y_[2*i] = Q[2*i+1] 
        for j in range(i+1, i+nb_nodes):
            j = j % nb_nodes
            u[i,j] = Unit_vect(Q[2*i], Q[2*j])
            Y_[2*i+1] += (1./m)*Contact_Mat[i,j]*( k*(norm(Q[2*j]-Q[2*i]) - L[i,j])*u[i,j]
                                            +  mu*(Q[2*j+1] - Q[2*i+1])@u[i,j]*u[i,j] )                         
    return np.reshape(Y_, (nb_nodes*4))  

t = np.linspace(0, T_f, N)
#initial condition for system
Y0_ = np.array([])
for i in range(nb_nodes):
    Y0_ = np.append(Y0_, tri.points[i])
    Y0_ = np.append(Y0_, V0[i])
print("initial condition :",Y0_)
Sol = solve_ivp(System, [0, T_f], Y0_, t_eval=t)

#animation ice floe:
fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(2.5, 4.5), ylim=(-.5, 1.5))
ax.set_aspect('equal')
ax.grid()
line, = ax.plot([], [], 'o-', lw=2)
center,= ax.plot([], [], 'o', color='r')
time_template = 'time = %.9fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def init():
    line.set_data([], [])
    time_text.set_text(' ')
    return line, time_text

possible = []
for triangle in tri.simplices:
    for index1 in triangle:
        for index2 in triangle:
            if index1 != index2:
                possible.append((min(index1, index2), max(index1, index2)))
#print(possible)
print("set of all edge=", set(possible))
#define a graph
g = UndirectedGraph(nb_nodes)
for v1,v2 in set(possible):
    g.AddEdge(v1, v2)
print("set graph=",g.AsSet)

while True:
    try:
        Route = g.RouteInspection()
        break
    except:
        continue
print("chinese postman solution ",Route) #optimized route to use in line.setdata 

# Route = g.RouteInspection()
# Route = [0, 3, 1, 5, 2, 1, 0, 2, 4, 0, 3, 6, 0, 6, 4, 5]
def animate_spring(i):
    Ix = [i for i in range(0, nb_nodes*4, 4)]
    Iy = [i for i in range(1, nb_nodes*4, 4)]
    thisx = []
    thisy = []
    for j in Ix:
        thisx = np.append(thisx, Sol.y[j][i])
    for j in Iy:
        thisy = np.append(thisy, Sol.y[j][i])
    for k in Route:
        thisx = np.append(thisx,thisx[k])
        thisy = np.append(thisy,thisy[k])
    line.set_data(thisx[nb_nodes:], thisy[nb_nodes:])
    time_text.set_text(time_template % (i*dt))
    return line, time_text

ani = animation.FuncAnimation(fig, animate_spring, np.arange(1, len(Sol.y[0])),
                                interval=2.5, blit=False, init_func=init)
# ani.save('intensity=5_2.mp4', fps=100)
plt.show()








