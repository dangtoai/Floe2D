# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 14:31:44 2022

@author: user
"""

from Func import *
from graph import *
import numpy as np
import matplotlib.animation as animation
from matplotlib.animation import FFMpegFileWriter
from IPython import display

if __name__ == '__main__':
    ######################
    ##### An ice floe ####
    ######################

    # first floe
    Points = np.array([[0., 0.], [1.1, 0], [1, 1], [0., 1]])
    Springs = {(0,1),(0,2),(0,3),(2,3),(1,2)}
    
    Nodes = []
    V0 = np.array([0.75, 0.])
    for i in range(len(Points)):
        Nodes.append(Node(Points[i], V0, i))
    
    k = 100.
    floe = Floe(nodes=Nodes, springs=Springs,
                stiffness= k, viscosity=k/10., id_number=1)
    
    Problem = Percussion_Wall(floe)


    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=True,
                          xlim=(0., 2.1), ylim=(-.5, 1.5))
    ax.set_aspect('equal')
    ax.grid()
    plt.axvline(x=2., color="red")
    line1, = ax.plot([], [], '.-', lw=.95)
    time_template = 'time = % 10fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    Route = floe.Route()
    
    All_positions_velocities = Problem.simulation()
    
    # Energy study:
 
    Traction_energy = np.zeros(len(All_positions_velocities[0]))
    Torsion_energy = np.zeros(len(All_positions_velocities[0]))

    for index in range(len(All_positions_velocities[0])):
        Sum1 = 0.
        Sum2 = 0.
        for i, j, k in floe.simplices():
            Qi = np.array([All_positions_velocities[4*i] [index],
                          All_positions_velocities[4*i+1][index]])
            Qj = np.array([All_positions_velocities[4*j] [index],
                          All_positions_velocities[4*j+1][index]])
            Qk = np.array([All_positions_velocities[4*k] [index],
                          All_positions_velocities[4*k+1][index]])

            Sum1 += 0.5 * ( (floe.k/sin(Angle_Mat[i, k, j])) * (norm(Qi-Qj) - Length_Mat[i, j])**2
                                    + (floe.k/sin(Angle_Mat[i, j, k])) * (norm(Qi-Qk) - Length_Mat[i, k])**2
                                    + (floe.k/sin(Angle_Mat[j, i, k])) * (norm(Qj-Qk) - Length_Mat[j, k])**2)

            Sum2 += 0.5 * (Torsion_Mat[i, j, k] * (Angle(Qi, Qj, Qk) - Angle_Mat[i, j, k])**2
                            + Torsion_Mat[i, k, j] * (Angle(Qi, Qk, Qj) - Angle_Mat[i, k, j])**2
                            + Torsion_Mat[j, i, k] * (Angle(Qj, Qi, Qk) - Angle_Mat[j, i, k])**2)

        Traction_energy[index] = Sum1
        Torsion_energy[index] = Sum2

    E_tot = Traction_energy + Torsion_energy

    plt.figure()
    plt.plot(Traction_energy, "-", label="Traction")
    plt.plot(Torsion_energy, "-", label="Torsion")
    plt.plot(E_tot, "-", label="Total")
    plt.legend()

    #animation

    def init():
        line1.set_data([], [])
        time_text.set_text('')
        return line1, time_text

    def animate_spring(i):
        Ix = [j for j in range(0, floe.n*4, 4)]
        Iy = [j for j in range(1, floe.n*4, 4)]
        thisx = []
        thisy = []
        for j in Ix:
            thisx = np.append(thisx, All_positions_velocities[j][i])
        for j in Iy:
            thisy = np.append(thisy, All_positions_velocities[j][i])
        for k in Route:
            thisx = np.append(thisx, thisx[k])
            thisy = np.append(thisy, thisy[k])
        line1.set_data(thisx[floe.n:], thisy[floe.n:])
        time_text.set_text(time_template % (i*dt))
        return line1, time_text

    ani = animation.FuncAnimation(fig, animate_spring,
                                  np.arange(0, len(All_positions_velocities[0])-200), interval=25, blit=False)

    
    
    
    
    
    
    
    
    
    