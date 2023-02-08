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

T = 6.  # time end
t = np.linspace(0, T, 1000)

l0 = 1.
k = 100.
m = 1.
v0 = 31.15
mu = 10.


def position(t):
    return v0 * np.sqrt(m/k) * np.sin(t*np.sqrt(k/m)) + l0


def velocity(t):
    return v0 * np.cos(t*np.sqrt(k/m))


def E_c(t):
    return 0.5 * m * velocity(t)**2


def E_el(t):
    return 0.5 * k * (position(t)-l0)**2


def E_tot(t):
    return E_c(t) + E_el(t)


def E_max(l_, k_):
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


# resolution l equation de conservation de lenergie

def SecondMembre(Y, t):
    return [Y[1], -(k/m) * Y[0] + k*l0/m]

# sol = odeint(SecondMembre, [l0, v0], t)

# plt.figure()
# plt.plot(t, sol[:,0], label= "position")
# plt.plot(t, sol[:,1], label= "velocity")
# plt.xlabel("time(s)")
# plt.legend()
# plt.show()


# resolution du systeme dynamique de masse-ressort sans viscosite

def System_(Y, t):
    dydt = [Y[2], 0, k/m * (abs(Y[0]-Y[1]) - l0)*-
            np.sign(Y[0]-Y[1]) + mu/m * (Y[3]-Y[2]), 0]

    return dydt


sol = odeint(System_, [l0, 0, -v0-1, 0], t)


plt.figure()
plt.plot(t, sol[:, 0], label="position")
plt.plot(t, sol[:, 2], label="velocity")
plt.xlabel("time(s)")
plt.legend()
plt.show()


def E_c_num(sol_):
    return 0.5 * m * sol_**2


def E_el_num(sol_):
    return 0.5 * k * (sol_-l0)**2


def E_tot_num(sol_, sol__):
    return E_c_num(sol__) + E_el_num(sol_)

# plt.figure()
# plt.plot(t, position(t), label = "position")
# plt.plot(t, velocity(t), label = "velocity")
# plt.xlabel("time(s)")
# plt.legend()
# plt.show()


plt.figure()
plt.plot(t, E_c_num(sol[:, 2]), label="$E_c$")
plt.plot(t, E_el_num(sol[:, 0]), label="$E_{el}$")
plt.plot(t, E_tot_num(sol[:, 0], sol[:, 2]), label="$E_{tot}$")
plt.hlines(E_max(l0, k), 0, T, linestyles='dashed', label="$E_{max}$")
plt.xlabel("time(s)")
plt.legend()
plt.show()


# fig = plt.figure()
# ax = fig.add_subplot(111, autoscale_on=True,)
# line1, = ax.plot([], [], '.-', lw=1.95)
# line2, = ax.plot([], [], '.-', lw=1.95)
# time_template = 'time = % 10fs'
# time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
# Route = floe.Route()
# def init():
#     line1.set_data([], [])
#     time_text.set_text(' ')
#     return line1, time_text

# def animate_spring(i):
#     Ix = [j for j in range(0, floe.n*4, 4)]
#     Iy = [j for j in range(1, floe.n*4, 4)]
#     thisx = []
#     thisy = []
#     for j in Ix:
#         thisx = np.append(thisx, All_positions_velocities[j][i])
#     for j in Iy:
#         thisy = np.append(thisy, All_positions_velocities[j][i])
#     for k in Route:
#         thisx = np.append(thisx, thisx[k])
#         thisy = np.append(thisy, thisy[k])

#     line1.set_data(thisx[floe.n:],
#                    thisy[floe.n:])

#     time_text.set_text(time_template % (i*dt))
#     return line1, time_text

# ani = animation.FuncAnimation(fig, animate_spring,
#                                   np.arange(0, len(All_positions_velocities[0])), interval=2, blit=False)
