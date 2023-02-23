#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 14:32:42 2023

@author: phandangtoai
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# import matplotlib
# from collections import deque

N = 20       #nombre de noeuds 
L = 1.      #longueur du reseau
k = 10000.     #raideur
m = 1.      #masse
mu = 400.     #viscosity
x0 = np.array(np.linspace(0,L,N))       #position initiale des noeuds
v0 = np.zeros(N-1)                  
v  = np.array([-0.25*10])               #vitesse init noeud libre
v0 = np.concatenate([v0, v])            #vitesse initiale des noeuds
L0 = L/(N-1)*np.ones(N-1)               #longueurs a vides de chaque ressorts
t_simu = 2.
n  = 500
dt = t_simu/n
history_len = 500

#on utilise la meme notation que Desmond dans son rapport

def Model(L, N, x0, v0):
    diagB = -2.0 * k * np.ones((N))
    diagB[0] = -k
    diagB[-1] = -k
    B = np.diag(diagB / m) + np.diag(k * np.ones((N - 1)) / m, 1) + np.diag(k * np.ones((N - 1)) / m, -1)

    diagC = -2.0 * mu * np.ones((N))
    diagC[0] = -mu
    diagC[-1] = -mu
    C = np.diag(diagC / m) + np.diag(mu * np.ones((N - 1)) / m, 1) + np.diag(mu * np.ones((N - 1)) / m, -1)
    
    E = np.zeros((2 * N, 2 * N))
    E[:N, N:] = np.identity(N)
    E[N:, :N] = B
    E[N:, N:] = C
    E[0, N]=0.
    E[N,:] = E[N,:]*0
    
    F = np.zeros((2*N, 2*N-2))
    F[N:, :N-1] = (np.diag(-k * np.ones((N)) / m) + np.diag(k * np.ones((N-1)) / m, -1))[:, :N-1]
    F[N,0] = 0. 
    
    Y0 = np.concatenate([x0, v0])
    t = np.linspace(0, t_simu, n + 1)
    full_L0 = np.concatenate([L0, np.zeros((N-1))])
    # return E

    def model(Y, t):
        return E @ Y + F @ full_L0

    return t, odeint(model, Y0, t)
    # return E, F, Y0, full_L0

t,y  = Model(L, N, x0, v0)
# E, F, Y0, LL = Model(L, N, x0, v0)
# print(E)
# print(F)
# print(Y0)
# print(LL)

# plt.figure()
# for i in range(N):
#     plt.plot(t,y[:,i],label = i)
# plt.xlabel("times")
# plt.ylabel("positions")    
# plt.legend()

# plt.figure()
# for i in range(N,2*N):
#     plt.plot(t,y[:,i],label = i-N)
# plt.xlabel("times")
# plt.ylabel("velocities")   
# plt.legend()

# Variation = []          #variation de l'atome i par rapport a sa pos init
# for i in range(N):
#     Variation.append(max(abs(y[:,i]-x0[i])))

# plt.figure()
# plt.bar(np.arange(N), Variation, width = 0.1, label='123')
# plt.xlabel("i-th node")
# plt.ylabel("$\epsilon$")
# plt.title("Stable zone of each node")

##mesure la deformation. 

# print(y[500,:N]) # la position courant des noeuds à l'instant 500 
plt.figure()
for i in range(0,50,5):
    plt.plot(x0, y[i,:N], "--")
plt.xlabel("x") 
plt.ylabel("$\phi$- Deformation field")

## mesure le champs de deplacement 1D
plt.figure()
for i in range(0,50,5):
    plt.plot(x0, y[i,:N]-y[0,:N], "--")
plt.xlabel("x") 
plt.ylabel("$u=\phi-Id$- Displacement field")



#animation positions
# fig = plt.figure(figsize=(5, 4))
# ax = fig.add_subplot(autoscale_on=False, xlim=(-0.1, L+0.2), ylim=(-.2, .2))
# ax.set_aspect('equal')
# ax.grid()
# line, = ax.plot([], [], 'o-', lw=2)
# # trace, = ax.plot([], [], '.-', lw=1, ms=2)
# time_template = 'time = %.1fs'
# time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
# history_x, history_y = deque(maxlen=history_len), deque(maxlen=history_len)
# plt.title("Floe's evolution after choc")



# def animate(i):
#     thisx = [y[:,j][i] for j in range(N)]
#     thisx.insert(0,0)
#     thisy = np.zeros_like(thisx)
#     if i == 0:
#         history_x.clear()
#         history_y.clear()

#     history_x.appendleft(thisx[2])
#     history_y.appendleft(thisy[2])

#     line.set_data(thisx, thisy)
#     # trace.set_data(history_x, history_y)
#     time_text.set_text(time_template % (i*dt))
#     return line, time_text

# ani = animation.FuncAnimation(fig, animate, len(y), interval=dt*1000, blit=True)
# writergif = animation.PillowWriter(fps=30)
# ani.save('filename.gif',writer=writergif)
# plt.show()



#animation deformation

# fig1, ax1 = plt.subplots()
# deformation, = ax1.plot(x0, y[0, :N], )
# time_template = 'time = %.1fs'
# time_text = ax1.text(0.05, 0.9, '', transform=ax1.transAxes)
# plt.title("$\phi$ - Deformation function")

# def animate1(i):
#     deformation.set_ydata(y[i,:N])  # update the data
#     time_text.set_text(time_template % (i*dt))
#     return deformation, time_text

# # Init only required for blitting to give a clean slate.
# def init():
#     deformation.set_ydata(np.ma.array(x0, mask=True))
#     return deformation, 

# ani1 = animation.FuncAnimation(fig1, animate1, len(y), init_func=init,
#                                interval=dt*1000, blit=True)
# writergif = animation.PillowWriter(fps=30)
# ani1.save('filename.gif',writer=writergif)
# # plt.show()

# animation displacement field
# fig2, ax2 = plt.subplots()
# ax2.set_ylim([-0.075, 0.075])

# displacement, = ax2.plot(x0, y[0, :N] - y[0, :N])
# time_template = 'time = %.1fs'
# time_text = ax2.text(0.05, 0.9, '', transform=ax2.transAxes)
# plt.title("$u=\phi-Id$ - Displacement function")

# def animate2(i):
#     displacement.set_ydata(y[i, :N] - y[0, :N])  # update the data
#     time_text.set_text(time_template % (i*dt))
#     return displacement, time_text


# # # Init only required for blitting to give a clean slate.
# def init2():
#     displacement.set_ydata(np.ma.array(x0, mask=True))
#     return displacement,

# ani2 = animation.FuncAnimation(fig2, animate2, len(y), init_func=init2,
#                                interval=dt*1000, blit=True)
# ani2.save('filename.gif',writer=writergif)

# plt.show()
