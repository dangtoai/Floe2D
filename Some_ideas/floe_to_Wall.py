# -*- coding: utf-8 -*-
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
    # Points = np.array([[0.25298236, 0.59733394],
    #   [0.43479153, 0.0089861 ],
    #   [0.77938292, 0.38657128],
    #   [0.19768507, 0.04416006],
    #   [0.86299324, 0.95665297],
    #   [0.98340068, 0.43614665],
    #   [0.16384224, 0.94897731]])

    Points = np.array([[0.44080984, 0.55885409],
            [0.02987621, 0.25925245],
            [0.45683322, 0.4151012 ],
            [0.64914405, 0.28352508],
            [0.27848728, 0.69313792],
            [0.6762549 , 0.44045372],
            [0.59086282, 0.15686774],
            [0.02398188, 0.54464902]])

    # Points = np.array([[0.33033482, 0.26682728],
    #                     [0.20464863, 0.62113383],
    #                     [0.61927097, 0.52914209],
    #                     [0.29965467, 0.13457995]])

    Nodes = []
    V0 = np.array([0.75, 0.])

    for i in range(len(Points)):
        Nodes.append(Node(Points[i], V0, i))

    # Springs = {(0, 1), (2, 4), (1, 2), (0, 4), (1, 5), (4, 6), (0, 3), (0, 6), (4, 5), (0, 2), (3, 6), (2, 5), (1, 3)}
    Springs = {(0, 7), (1, 2), (0, 4), (2, 7), (2, 3), (1, 7), (0, 2), (2, 6), (4, 5), (0, 5), (3, 6), (1, 6), (2, 5), (4, 7), (3, 5)}
    # Springs = {(0, 1), (1, 2), (0, 3), (2, 3), (0, 2), (1, 3)}
    k = 1000.
    floe = Floe(nodes=Nodes, springs=Springs,
                stiffness=k, viscosity=k/10., id_number=1)

    # t_end = 4.
    # collision_dist = 0.01
    # coef_restitution = 0.9
    
    Problem = Percussion_Wall(floe, Wall = 2.) 
    
    All_positions_velocities = Problem.simulation()
    

    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=True,
                          xlim=(.5, 2.1), ylim=(-.5, 1.1))
    ax.set_aspect('equal')
    ax.grid()
    plt.axvline(x=2., color="red")
    line1, = ax.plot([], [], '.-', lw=.95)
    time_template = 'time = % 10fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    Route = floe.Route()


    # Energy study:
    # Traction's energy, same as Dimitri
    Traction_energy = np.zeros(len(All_positions_velocities[0]))
    Torsion_energy = np.zeros(len(All_positions_velocities[0]))
    
    
    Traction_Mat = floe.traction_mat()
    Length_Mat  = floe.length_mat()
    Torsion_Mat = floe.torsion_mat()
    Angle_Mat = floe.angle_init()
    
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

            Sum1 += 0.5 * ((floe.k/sin(Angle_Mat[i, k, j])) * (norm(Qi-Qj) - Length_Mat[i, j])**2
                                    + (floe.k/sin(Angle_Mat[i, j, k])) * (norm(Qi-Qk) - Length_Mat[i, k])**2
                                    + (floe.k/sin(Angle_Mat[j, i, k])) * (norm(Qj-Qk) - Length_Mat[j, k])**2)

            Sum2 += 0.5 * (Torsion_Mat[i, j, k] * (Angle(Qi, Qj, Qk) - Angle_Mat[i, j, k])**2
                            + Torsion_Mat[i, k, j] * (Angle(Qi, Qk, Qj) - Angle_Mat[i, k, j])**2
                            + Torsion_Mat[j, i, k] * (Angle(Qj, Qi, Qk) - Angle_Mat[j, i, k])**2)

        Traction_energy[index] = Sum1
        Torsion_energy[index] = Sum2

    E_tot = Traction_energy + Torsion_energy

    plt.figure()
    plt.plot(Traction_energy [200:400], "-",color="green", label="Traction")
    plt.plot(Torsion_energy[200:400], "-", color="blue",  label="Torsion")
    plt.plot(E_tot[200:400], "--", color= "r" ,label="Total")
    plt.xlabel("time's step")
    plt.ylabel("Energy")
    plt.tight_layout()
    plt.legend()
    # animation

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
                                  np.arange(200, len(All_positions_velocities[0])-200), interval=25, blit=False)
    # ani.save("floe_to_wall.gif", writer='pillow')
