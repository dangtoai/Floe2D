# -*- coding: utf-8 -*-
from Func import *
from graph import *
import numpy as np

if __name__ == '__main__':
    ######################
    ##### An ice floe ####
    ######################
    
    Points = np.array([[1, 1.], [0., -1.], [-0.5, 0.]])
    # Points = np.array([[1, 0.], [0., -1.], [-0., 1.]])
    V0     = np.array([1., 0.])
    # V1     = np.array([-1.0, 0.])
    Nodes  = []
    for i in range(len(Points)):
        Nodes.append(Node(Points[i], V0, i))

    Springs = {(0,1),(1,2),(0,2)}
    k = 100000.
    
    floe = Floe(nodes=Nodes, springs=Springs, stiffness=k, viscosity=k/100.,  id_number = 1 )
    # floe.plot_init()
    
    t_end = 6.
    collision_dist   = 0.01
    coef_restitution = 0.7
    
    # After_shock_floe = floe.evolution(t_end, t_end)
    # After_shock_floe.update_velocity_node(0, np.array([-1.5,  0]))
    # print(After_shock_floe.get_velocity())

    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=True, xlim=(-2.1, 6.5), ylim=(-1., 3.5))
    ax.set_aspect('equal')
    ax.grid()
    plt.axvline(x=5., color = "red")
    line1, = ax.plot([], [], '.-', lw= .95)
    time_template = 'time = % 10fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    Route = floe.Route()
    
    # compute at contact
    # update the new velocity
    
    #compute after contact
    Sol = floe.Move(t_end)
    # floe.plot_displacements(t_end)
    Node2x = Sol.y[4]
    Node2y = Sol.y[5]
    Node1x = Sol.y[0]
    Node1y = Sol.y[1]
    Node3x = Sol.y[8]
    Node3y = Sol.y[9]
    
    j = 0
    # tant que (node1x ou node2x depasse (5-collision_dist)): 
    #   do: enregistrer positions noeuds 
    #       calculer nouveau vitesse du noeud en contact (node1x ou node2x)
    #       concatener avec ancient position.     
    
    while np.any(np.append(Node1x,Node2x) >5. - collision_dist) and j<= 1 :
        
        k = np.where(Node2x >= 5. - collision_dist)[0][0] 
        l = np.where(Node1x >= 5. - collision_dist)[0][0]
        
        print(k,l, min(k,l))
        
        print("truncate after " , min(k,l))
        m = len(Node1x)
        Node1x = Node1x[:min(k,l)]
        Node1y = Node1y[:min(k,l)]
        Node2x = Node2x[:min(k,l)]
        Node2y = Node2y[:min(k,l)]
        Node3x = Node3x[:min(k,l)]
        Node3y = Node3y[:min(k,l)]
        
        After_shock_floe = floe.evolution(t_end, t_end *((min(k,l)-m)%800/800.)) 
        print("last position of floe before collision \n", After_shock_floe.get_nodes())
        # print("last velocity of floe before collision \n", After_shock_floe.get_velocity())
        
        if min(k,l) == l:
            velocity_after_shock =  -coef_restitution*After_shock_floe.get_velocity()[0]
            After_shock_floe.update_velocity_node(0, velocity_after_shock)
            
        if min(k,l) == k:
            velocity_after_shock =  -coef_restitution*After_shock_floe.get_velocity()[1]
            After_shock_floe.update_velocity_node(1, velocity_after_shock)
        # print("new velocity of floe \n" , After_shock_floe.get_velocity())
        
        floe = After_shock_floe
        Node1x = np.append(Node1x, After_shock_floe.Move(t_end).y[0])
        Node1y = np.append(Node1y, After_shock_floe.Move(t_end).y[1])
        Node2x = np.append(Node2x, After_shock_floe.Move(t_end).y[4])
        Node2y = np.append(Node2y, After_shock_floe.Move(t_end).y[5])
        Node3x = np.append(Node3x, After_shock_floe.Move(t_end).y[8])
        Node3y = np.append(Node3y, After_shock_floe.Move(t_end).y[9])
        # break
        j = j+1
    print("there are", j, "collisions")
        
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
            thisx = np.append(thisx, [Node1x[i], Node2x[i], Node3x[i]])
        for j in Iy:
            thisy = np.append(thisy, [Node1y[i], Node2y[i], Node3y[i]])
        for k in Route:
            thisx = np.append(thisx,thisx[k])
            thisy = np.append(thisy,thisy[k])
        line1.set_data(thisx[floe.n:], thisy[floe.n:])
        time_text.set_text(time_template % (i*dt))
        return line1, time_text
    
    ani = animation.FuncAnimation(fig, animate_spring, 
                                    np.arange(400,1000), interval=25, blit=False)
    # ani.save("3nodes_to_Wall_2.mp4", fps=100)
