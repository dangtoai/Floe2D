#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 14:32:42 2023

@author: phandangtoai
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


N = 10       #nombre de noeuds 
L = 1.      #longueur du reseau
k = 160.     #raideur
m = 1.      #masse
mu = 40.     #viscosity
x0 = np.array(np.linspace(0,L,N))       #position initiale des noeuds
v0 = np.zeros(N-1)                  
v  = np.array([-0.25])
v0 = np.concatenate([v0, v])            #vitesse initiale des noeuds
L0 = L/(N-1)*np.ones(N-1)               #longueurs vide de chaque ressorts
t_simu = 2.
n  = 500

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

plt.figure()
for i in range(N):
    plt.plot(t,y[:,i])
plt.show()

Variation = []          #variation de l'atome i par rapport a sa pos init
for i in range(N):
    Variation.append(max(abs(y[:,i]-x0[i])))

plt.figure()
plt.bar(np.arange(N), Variation, width = 0.05)

###mesure la deformation. 

# print(y[500,:N]) # la position de chaque noeud à l'instant 500 du discretisation de temps






