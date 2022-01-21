# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 14:06:04 2021

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
from Func import *

#constant in this problem
nb_nodes = 3        #number of nodes 
# T_f      = 4.       #final time simulation
T_f = 4.
N = 500
dt = T_f/N          #time's step

m  = 6.2
mu = 3.
k  = 23.3

#Contact matrix
Contact_Mat = np.ones(nb_nodes*nb_nodes)
I = np.eye(nb_nodes)
Contact_Mat = Contact_Mat.reshape(nb_nodes,nb_nodes) - I
#print(Contact_Mat)

#Initial velocity 
# v1 = 1.5
# v2 = 1.3
# v3 = 1.2

v1 = 0.3
v2 = 0.1
v3 = 0.1

theta1 = 180
theta2 = 270
theta3 = 240
theta1, theta2, theta3 = np.deg2rad(theta1), np.deg2rad(theta2), np.deg2rad(theta3)

#intial position
q1_init = np.array([0.8, 0.])
q2_init = np.array([0.5, 0.5])
q3_init = np.array([-0.5, 0.])

#initial velocity
dq1_init = v1*np.array([np.cos(theta1), np.sin(theta1)])
dq2_init = v2*np.array([np.cos(theta2), np.sin(theta2)])
dq3_init = v3*np.array([np.cos(theta3), np.sin(theta3)])

# initial condition 
Y0 = np.stack([q1_init,dq1_init,q2_init,dq2_init,q3_init,dq3_init])
Y0_ = Y0.reshape((nb_nodes*4))

#Length matrix init
L = np.zeros((nb_nodes, nb_nodes))
for i in range(nb_nodes):
    for j in range(nb_nodes):
        L[i,j] = norm(Y0[2*j]-Y0[2*i])

t = np.linspace(0, T_f, N)
Sol = solve_ivp(System, [0, T_f], Y0_, t_eval=t, args=( Y0, nb_nodes, Contact_Mat, L, m, mu, k ))

Label = ["$q_{0,x}$","$q_{0,y}$",
         "$\dot{q_{0,x}}$","$\dot{q_{0,y}}$",
         "$q_{1,x}$","$q_{1,y}$",
         "$\dot{q_{1,x}}$","$\dot{q_{1,y}}$",
         "$q_{2,x}$","$q_{2,y}$",
         "$\dot{q_{2,x}}$","$\dot{q_{2,y}}$"]

I = [0,1,4,5,8,9] #indice for plot
for i in I:
    plt.plot(t, Sol.y[i], label = Label[i])
#plt.ylim([-1,1])
plt.xlabel("time(s)")
plt.ylabel("position")
plt.legend()
plt.tight_layout()
# plt.savefig("3nodes_1fixed")

#animation ice floe:
fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=True, xlim=(-1, 1), ylim=(-1, 1))
# ax = fig.add_subplot(111, autoscale_on=False)
ax.set_aspect('equal')
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
center,= ax.plot([], [], 'o', color='r')
time_template = 'time = %.9fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

Ix = [0,4,8,0]
Iy = [1,5,9,1]
def animate_spring(i):
    thisx = []
    thisy = []
    for j in Ix:
        thisx = np.append(thisx, Sol.y[j][i])
    for j in Iy:
        thisy = np.append(thisy, Sol.y[j][i])
    thisx = thisx.tolist()
    thisy = thisy.tolist()
    line.set_data(thisx, thisy)
    time_text.set_text(time_template % (i*dt))
    return line, time_text

def animate_center(i):
    thisx = []
    thisy = []
    for j in Ix:
        thisx = np.append(thisx, Sol.y[j][i])
    for j in Iy:
        thisy = np.append(thisy, Sol.y[j][i])
    thisx = thisx.tolist()
    thisx = [sum(thisx)/len(thisx)]
    thisy = thisy.tolist()
    thisy = [sum(thisy)/len(thisy)]
    center.set_data(thisx[:], thisy[:])
    return center, 

ani = animation.FuncAnimation(fig, animate_spring, np.arange(1, len(Sol.y[0])),
                                interval=.1, blit=False, init_func=init)
ani2 = animation.FuncAnimation(fig, animate_center, np.arange(1, len(Sol.y[0])),
                                interval=.1, blit=False, init_func=init)
# ani.save('3nodes1stable_bigk.mp4', fps=1000)
# plt.show()


