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
    
    Points = np.array([[0.25298236, 0.59733394],
      [0.43479153, 0.0089861 ],
      [0.77938292, 0.38657128],
      [0.19768507, 0.04416006],
      [0.86299324, 0.95665297],
      [0.98340068, 0.43614665],
      [0.16384224, 0.94897731]])
    
    Nodes = []
    V0 = np.array([0.25, 0.])
    # V1 = np.array([0.5, 0.])

    for i in range(len(Points)):
        Nodes.append(Node(Points[i], V0, i))
    # Nodes[1] = Node(Points[1], V1, 1)
    # Nodes[5] = Node(Points[5], -V1, 5)

    Springs = {(0, 1), (2, 4), (1, 2), (0, 4), (1, 5), (4, 6), (0, 3), (0, 6), (4, 5), (0, 2), (3, 6), (2, 5), (1, 3)}
    
    k = 10.
    floe = Floe(nodes=Nodes, springs=Springs,
                stiffness=k, viscosity=k/2., id_number=1)
    # floe.plot_init()

    t_end = 6.
    collision_dist = 0.01
    coef_restitution = 0.9

    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=True, xlim=(-.1, 2.1), ylim=(-0.5, 1.1))
    ax.set_aspect('equal')
    ax.grid()
    plt.axvline(x=2., color="red")
    line1, = ax.plot([], [], '.-', lw=.95)
    time_template = 'time = % 10fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    Route = floe.Route()

    Sol = floe.Move(t_end)

    Positions = Sol.y

    # tant que: la position en x de n'importe quel noeud depasse Wall:
    #   Faire:  tronquer toutes les vitesses et positions jusqu'a cette indice
    #           update velocity of collision'nodes
    #           lancer la simulation et concatener avec les positions et noeuds.

    Ix = [j for j in range(0, floe.n*4, 4)]

    All_x_positions = []
    All_positions_velocities = [[]]*floe.n*4 #save all positions of all nodes
    for i in range(floe.n*4):
        All_positions_velocities[i] = Positions[i]
    
    # save positions and velocity of the first simulation
    for i in Ix:
        All_x_positions = np.append(All_x_positions, Positions[i])

    j = 0
    while np.any(All_x_positions > 2-collision_dist) and j <=10: 
        m = len(All_positions_velocities[0])
        # print([np.where(Positions[i] >= 2-collision_dist)[0][0] for i in Ix])
        liste = []
        for i in Ix:
            liste.append(np.where(All_positions_velocities[i] >= 2-collision_dist)[0])
        for i in range(len(liste)):
            if len(liste[i])==0:liste[i] = [10000]
        # liste = [np.where(All_positions_velocities[i] >= 2-collision_dist)[0][0] for i in Ix]
        for i in range(len(liste)): liste[i] = min(liste[i])
        k = min(liste)
        index_nodes_contact = np.where(liste == min(liste))[0][0]
        print(k)
        for i in range(floe.n*4):
            All_positions_velocities[i] = All_positions_velocities[i][:k]
        
        After_shock_floe = floe.evolution(t_end, t_end *((k-m)%800)/800.)
        print("last position of floe before collision \n", After_shock_floe.get_nodes())
        
        #update velocity of nodes at collision!!!
        velocity_after_shock =  -coef_restitution*After_shock_floe.get_velocity()[index_nodes_contact]
        After_shock_floe.update_velocity_node(index_nodes_contact, velocity_after_shock)
        #
        floe = After_shock_floe
        New_position_velocities = After_shock_floe.Move(t_end).y
        for i in range(floe.n*4):
            All_positions_velocities[i] = np.append(All_positions_velocities[i], New_position_velocities[i])
        # break
        All_x_positions = []
        for i in Ix:
            All_x_positions = np.append(All_x_positions, All_positions_velocities[i])
        j = j+1
        
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
                                  np.arange(0, len(All_positions_velocities[0])), interval=25, blit=False)

    # ani.save("floe_to_wall.gif", writer=PillowWriter(fps=50))
    
    
    
    
    
    
    
    
    