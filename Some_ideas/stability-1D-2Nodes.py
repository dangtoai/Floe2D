#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 16:53:41 2022

@author: phandangtoai
"""

from Func import *
# from graph import *
import numpy as np
from scipy.integrate import odeint

T = 3.      #time end
t = np.linspace(0,T,1000)

l0 = 1.
k  = 10.
m  = 1.
v0 = 3.15



def position(t):
    return v0* np.sqrt(m/k) * np.sin(t*np.sqrt(k/m)) + l0

def velocity(t):
    return v0* np.cos(t*np.sqrt(k/m))

def E_c(t):
    return 0.5 * m * velocity(t)**2

def E_el(t):
    return 0.5 * k * (position(t)-l0)**2

def E_tot(t):
    return E_c(t) + E_el(t)

def E_max(l_ ,k_):
    return 0.5 * l_ * k_

# plt.figure()
# plt.plot(t, position(t), label = "position")
# plt.plot(t, velocity(t), label = "velocity")
# plt.xlabel("time(s)")
# plt.legend()
# plt.show()

# plt.figure()
# plt.plot(t, E_c(t), label = "$E_c$")
# plt.plot(t, E_el(t), label = "$E_{el}$")
# plt.plot(t, E_tot(t), label = "$E_{tot}$")
# plt.hlines(E_max(l0, k), 0, T,linestyles='dashed', label = "$E_{max}$")
# plt.xlabel("time(s)")
# plt.legend()
# plt.show()

# plt.figure()

#resolution l equation de conservation de lenergie

def SecondMembre(Y, t):
    return [ Y[1], -(k/m) * Y[0] + k*l0/m ]

sol = odeint(SecondMembre, [l0, v0], t)

# plt.figure()
# plt.plot(t, sol[:,0], label= "position")
# plt.plot(t, sol[:,1], label= "velocity")
# plt.xlabel("time(s)")
# plt.legend()
# plt.show()


#resolution du systeme dynamique de masse-ressort sans viscosite

def System_(Y, t):
    dydt = [Y[2], 0, k/m *(abs(Y[0]-Y[1]) - l0)*-np.sign(Y[0]-Y[1]), 0]
    
    return dydt

sol_ = odeint(System_, [l0, 0, v0+1, 0], t)



plt.figure()
plt.plot(t, sol_[:,0], label= "position")
# plt.plot(t, sol_[:,2], label= "velocity")
plt.xlabel("time(s)")
plt.legend()
plt.show()










