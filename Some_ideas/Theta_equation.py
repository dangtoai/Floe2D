# -*- coding: utf-8 -*-
import numpy as np
from Func import *
from scipy.integrate import odeint
# import math
import matplotlib.pyplot as plt

J = 1.
C = 100.     #stiffness of torsion spring
f = 1.
T_end = 7.

R = 1. # traction spring's length
theta0 = np.pi/3.

def deriv(y,t):
    dxdt = y[1]
    dydt = -f/J *y[1] - C/J* y[0]
    dzdt = [dxdt, dydt]
    return dzdt

y0 = [30, 0.]
t = np.linspace(0,T_end,1000)

Sol = odeint(deriv, y0, t)

plt.figure()
plt.plot(t,Sol[:,0], label='$\Theta$(rad)')
plt.plot(t,Sol[:,1], label="$\dot{\Theta}$(rad/s)")
plt.tight_layout()
plt.xlabel("temps(s)")
plt.ylabel("rad")
plt.legend()

#calculer les coordonnés des noeuds en mouvement:

Node0x = np.zeros_like(t)
Node0y = np.zeros_like(t)

Node1x = np.zeros_like(t) + R
Node1y = np.zeros_like(t) 

Node2x = R*np.cos(theta0 + np.deg2rad(Sol[:,0]))
Node2y = R*np.sin(theta0 + np.deg2rad(Sol[:,0]))

fig = plt.figure()

ax = fig.add_subplot(111, autoscale_on=True, xlim=(-.6, 1.1), ylim=(-.5, 1.5))
ax.set_aspect('equal')
plt.plot([0, R*np.cos(theta0)],[0,R*np.sin(theta0)], marker = "o")
ax.grid()

line1, = ax.plot([], [], '.-', lw= .95)
time_template = 'time = % 10fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

dt = T_end/1000.

def init():
    line1.set_data([], [])
    time_text.set_text('')
    return line1, time_text

def animate_spring(i):
    thisx = []
    thisy = []
    I = [j for j in range(3)]
    for j in I:
        thisx = np.append(thisx, [Node1x[i], Node0x[i], Node2x[i]])
        thisy = np.append(thisy, [Node1y[i], Node0y[i], Node2y[i]])
    line1.set_data(thisx[-3:], thisy[-3:])
    time_text.set_text(time_template % (i*dt))
    # return thisx, thisy
    return line1, time_text


ani = animation.FuncAnimation(fig, animate_spring, 
                                np.arange(0,len(Node1x)), interval=25, blit=False)

# ani.save("1Torsion.gif", writer=PillowWriter(fps=50))



















