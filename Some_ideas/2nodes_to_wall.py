# -*- coding: utf-8 -*-
from Func import *
from graph import *
import numpy as np

if __name__ == '__main__':
    ######################
    ##### An ice floe ####
    ######################
    
    #global constant
    coef_restitution = .6
    collision_dist   = 0.01
    
    Points = np.array([[0., 0.],
                       [1.,0.]])
    Nodes = []
    t_end = 4.
    V0 = [np.array([0.25, 0.]), np.array([0.25, 0.])]
    # V0 = [np.array([0.5, 0.]), np.array([0.5, 0.])]
    for i in range(len(Points)):
        Nodes.append(Node(Points[i], V0[i], i))
    Springs = {(0, 1)}
    floe = Floe(nodes=Nodes, springs=Springs, stiffness=200., viscosity=20.0, id_number = 1 )
    
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=True, xlim=(-.1, 3.1), ylim=(-.5, .5))
    plt.plot(2,0, 'o', color='red')
    ax.set_aspect('equal')
    ax.grid()
    line1, = ax.plot([], [], '.-', lw= 2.5)
    center, = ax.plot([], [], 'o', lw= 1.5)
    time_template = 'time = % 10fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    
    Sol = floe.Move(t_end)
    # After_shock_floe = floe.evolution(t_end, t_end)
    
    #compute new velocity!!! 
    # First model using only coef of restitution because m' = \infty, v_0' and V_0' = (0,0) 
    
    # tant que Collision ( il existe un i telque Sol.y[4][i]>2-0.001)
    # Faire: delete tous les positions et vitesses a partir de i+1
    # recalculer la vitesse et position à partir de cet indice
    
    Node2x = Sol.y[4]
    Node2y = Sol.y[5]
    Node1x = Sol.y[0]
    Node1y = Sol.y[1]
    
    velo1x = Sol.y[2]
    velo1y = Sol.y[3]
    velo2x = Sol.y[6]
    velo2y = Sol.y[7]
    
    # print(Node1y)
    
    # print(Node2x.size)
    j = 0
    m = 0
    while Node2x.all() <2. - collision_dist and j < 2:
        k = np.where(Node2x >= 2. - collision_dist)[0][0] 
        print("tronquer apres" ,k)
        m += k
        Node2x = Node2x[:k]
        Node2y = Node2y[:k]
        Node1x = Node1x[:k]
        Node1y = Node1y[:k]
        
        velo1x = velo1x[:k]
        velo1y = velo1y[:k]
        velo2x = velo2x[:k]
        velo2y = velo2y[:k]
        # print(Node2x[-1])
        
        After_shock_floe = floe.evolution(t_end, t_end *(k%800)/800.)
        ###
        
        print("last position of floe before collision \n", After_shock_floe.get_nodes())
        print("last velocity of floe before collision \n", After_shock_floe.get_velocity())
        floe = After_shock_floe
        
        # print(Node2x[-20:])
        # print("last position of floe before collision \n", After_shock_floe.get_nodes())

        velocity_after_shock =  -coef_restitution*After_shock_floe.get_velocity()[-1]
        After_shock_floe.update_velocity_node(1, velocity_after_shock)
        print("new velocity of floe \n" , After_shock_floe.get_velocity())
        
        # print(" ------ " , After_shock_floe.Move(t_end).y[6][:5])
        Node2x = np.append(Node2x, After_shock_floe.Move(t_end).y[4])
        Node2y = np.append(Node2y, After_shock_floe.Move(t_end).y[5])
        Node1x = np.append(Node1x, After_shock_floe.Move(t_end).y[0])
        Node1y = np.append(Node1y, After_shock_floe.Move(t_end).y[1])
        
        velo1x = np.append(velo1x, floe.Move(t_end).y[2])
        velo1y = np.append(velo1y, floe.Move(t_end).y[3])
        velo2x = np.append(velo2x, After_shock_floe.Move(t_end).y[6])
        velo2y = np.append(velo2y, After_shock_floe.Move(t_end).y[7])
        
        # print(velo2x[792: 797])
        j += 1
        

    
    plt.figure()
    plt.plot(np.linspace(0,8,len(Node2x)), Node2x, label = "$q_{1x}$")
    plt.plot(np.linspace(0,8,len(Node1x)), Node1x, label = "$q_{2x}$")
    plt.xlabel("time (s)")
    plt.ylabel("x")
    plt.legend()
    
    
    
    plt.figure()
    plt.plot(np.linspace(0,8,len(Node2x)), velo1x, label = "$v_{1x}$")
    plt.plot(np.linspace(0,8,len(velo2x)), velo2x, label = "$v_{2x}$")
    plt.plot(np.linspace(0,8,len(Node2x)), velo1y, label = "$v_{1y}$")
    plt.plot(np.linspace(0,8,len(velo2y)), velo2y, label = "$v_{2y}$")
    plt.xlabel("time (s)")
    plt.ylabel("velocity")
    plt.legend()
    
    
    #animation ice floe:    
    def init():
        line1.set_data([], [])
        time_text.set_text('')
        return line1, time_text

    Ix = [0,4]
    Iy = [1,5]
    
    def animate2(i):
        thisx = []
        thisy = []
        # for j in Ix:
        #     thisx = np.append(thisx, np.append(Sol.y[j], Sol2.y[j])[i])
        # for j in Iy:
        #     thisy = np.append(thisy, np.append(Sol.y[j], Sol2.y[j])[i])
        thisx = np.append(thisx, np.append(Node1x[i], Node2x[i]))
        thisy = np.append(thisy, np.append(Node1y[i], Node2y[i]))
        thisx = thisx.tolist()
        thisy = thisy.tolist()
        line1.set_data(thisx[:], thisy[:])
        time_text.set_text(time_template % (i*dt))
        center.set_data([sum(thisx[:])/2.], [sum(thisy[:])/2.])
        # print(i)
        # print(thisx[:])
        return line1, center, time_text
        # return thisx, thisy
    
    ani = animation.FuncAnimation(fig, animate2, np.arange(600, 1000),
                                    interval=.1, blit=False, init_func=init)
    
    
    
    
    
    
    
    
    
    
    