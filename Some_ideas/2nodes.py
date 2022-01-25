# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 14:03:32 2022

@author: user
"""

from Func import *

if __name__=="__main__":
    
    #global constant
    nb_nodes = 2        #number of nodes 
    # T_f      = 4.       #final time simulation
    T_f = 4.
    N = 1000
    dt = T_f/N          #time's step
    m  = 6.2
    mu = 3.
    k  = 230.3
    
    #Desmond's parameters
    # m  = 1
    # mu = 1.3
    # k  = 18
    
    #Contact matrix
    Contact_Mat = np.ones(nb_nodes*nb_nodes)
    I = np.eye(nb_nodes)
    Contact_Mat = Contact_Mat.reshape(nb_nodes,nb_nodes) - I
    
    #Initial condition
    v1 = 0.
    v2 = 3.
    # v1 = 0.8
    # v2 = -0.8
    
    theta1 = 180
    theta2 = 90
    # theta1 = 90
    # theta2 = 90
    theta1, theta2 = np.deg2rad(theta1), np.deg2rad(theta2)

    #intial position
    #1nodefree
    # q1_init = np.array([0., 1.])
    # q2_init = np.array([0., 0.])
    
    #2nodesfree
    q1_init = np.array([0., 1.])
    q2_init = np.array([0., 0.])
    
    # q1_init = np.array([0., 0.])
    # q2_init = np.array([0., 0.])

    #initial velocity
    dq1_init = v1*np.array([np.cos(theta1), np.sin(theta1)])
    dq2_init = v2*np.array([np.cos(theta2), np.sin(theta2)])
    Y0 = np.stack([q1_init,dq1_init,q2_init,dq2_init])
    Y0_ = Y0.reshape((nb_nodes*4))
    
    #Length matrix init
    L = np.zeros((nb_nodes, nb_nodes))
    for i in range(nb_nodes):
        for j in range(nb_nodes):
            L[i,j] = norm(Y0[2*j]-Y0[2*i])
    
    
    t = np.linspace(0, T_f, N)
    Sol = solve_ivp(System, [0, T_f], Y0_, method='RK45', t_eval=t, args=( Y0, nb_nodes, Contact_Mat, L, m, mu, k ) )
    Label = ["$q_{0,x}$","$q_{0,y}$",
             "$\dot{q_{0,x}}$","$\dot{q_{0,y}}$",
             "$q_{1,x}$","$q_{1,y}$",
             "$\dot{q_{1,x}}$","$\dot{q_{1,y}}$",
             "$q_{2,x}$","$q_{2,y}$",
             "$\dot{q_{2,x}}$","$\dot{q_{2,y}}$"]


    I = [0,1,4,5] #index for plot
    I_ = [2,3,6,7] #index for plot
    Cx = (Sol.y[0]-Sol.y[4])**2
    Cy = (Sol.y[1]-Sol.y[5])**2
    spring_length = np.sqrt(Cx+Cy)
    
    for i in I:
        plt.plot(t, Sol.y[i], label = Label[i])
    plt.xlabel("time(s)")
    plt.ylabel("position")
    plt.legend()
    plt.tight_layout()
    # plt.savefig("3nodes_1fixed")
    
    plt.figure()
    for i in I_:
        plt.plot(t, Sol.y[i], label = Label[i])
    # plt.hlines(3, 0, T_f)
    plt.xlabel("time(s)")
    plt.ylabel("velocity")
    plt.legend()
    plt.tight_layout()
    
    plt.figure()
    plt.plot(t, spring_length, label= "spring length")
    plt.hlines(norm(q1_init-q2_init), 0, T_f, label="initial length")
    plt.xlabel("time(s)")
    plt.legend()
    plt.tight_layout()
    
    #animation ice floe:    
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-8, 8), ylim=(-1, 7))
    ax.set_aspect('equal')
    ax.grid()

    line,  = ax.plot([], [], 'o-', color='b', lw=2)
    center,= ax.plot([], [], 'o', color='r')
    time_template = 'time = %.9fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    def init():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text

    Ix = [0,4]
    Iy = [1,5]
    
    def animate2(i):
        thisx = []
        thisy = []
        for j in Ix:
            thisx = np.append(thisx, Sol.y[j][i])
        for j in Iy:
            thisy = np.append(thisy, Sol.y[j][i])
        thisx = thisx.tolist()
        thisy = thisy.tolist()
        line.set_data(thisx[:], thisy[:])
        time_text.set_text(time_template % (i*dt))
        # print(i)
        return line, time_text
    
    def animate_center(i):
        thisx = []
        thisy = []
        for j in Ix:
            thisx = np.append(thisx, Sol.y[j][i])
        for j in Iy:
            thisy = np.append(thisy, Sol.y[j][i])
        thisx = thisx.tolist()
        thisx = [sum(thisx)/len(thisx)]
        thisy = thisy.tolist()
        thisy = [sum(thisy)/len(thisy)]
        center.set_data(thisx[:], thisy[:])
        return center, 
    
    ani = animation.FuncAnimation(fig, animate2, np.arange(1, len(Sol.y[0])),
                                    interval=.1, blit=False, init_func=init)
    ani2 = animation.FuncAnimation(fig, animate_center, np.arange(1, len(Sol.y[0])),
                                    interval=.1, blit=False, init_func=init)
    
    # ani.save('v=10.mp4', fps=1000)
    # plt.show()



