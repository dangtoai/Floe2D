#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 10:52:55 2023

@author: phandangtoai
"""
import time

import matplotlib.pyplot as plt
from math import pi
import plotrc
from griffith import solver, log, problem_data
from griffith.mesh import Mesh, Boundary_Mesh, Triangle

from griffith.mesh import *
from griffith.geometry import Point, Segment
from scipy.interpolate import Rbf
from scipy.interpolate import LinearNDInterpolator
import numpy as np
# from griffith.mesh import 
# from Func import * 


if __name__ == '__main__':
    
    # define Lame tensor. 
    # T = problem_data.lame_tensor_ice
    # print(" Lame constant are: ", T._lambda, T._mu)
    
    #boudary condition, 
    # boundary_displacement = problem_data.Constant_Displacement_On_Y( traction_coefficient= 3.)
    
    #utiliser mesh obtenu par GMSH 
    mesh = Mesh("square.msh")
    
    # time_discretization = "smart time"
    # fracture_discretization = problem_data.Fracture_Discretization(args.angular_step, args.lengh_step, boundary_point=args.boundary_point, boundary_step=args.boundary_step,
                                                                   # interior_step=args.interior_step, interior_fast_init=args.interior_fast_init, interior_fast_step=args.interior_fast_step, interior_start_angle=args.interior_start_angle)
    # discretization_data =  problem_data.Discretization_Data(mesh, time_discretization, fracture_discretization, tip_enrichement=False)

    # test file geometry.py and mesh.py
    # print(mesh.nodes[0])
    # print(type(mesh.nodes[0].array))
    # print(mesh.boundary_edges_list())
    # mesh.plot()
    # print(len(mesh.triangles))
    # print(type(mesh.triangles[0]))
    
    # figax = plt.subplots()
    # fig, ax = figax
    # ax.set_aspect('equal')
    # for i in range(50):
    #     mesh.triangles[i].plot(figax)
    # mesh.triangles[2].plot(figax)
    
    # print(mesh.boundary_mesh.edges) 
    
    # print(mesh.boundary_mesh.edges[0]) # an edge on boundary
    # print(mesh.boundary_mesh.dirichlet_parts) # Dirichlet parts 
    # print(mesh.boundary_mesh.points[0]) #list of nodes on boundary
    
    # print(mesh.boundary_mesh.nodes[0].of_triangles) #set of triangle that this node belongs to
    # print(list(mesh.boundary_mesh.nodes[0].of_triangles)[0])
    # m = list(mesh.boundary_mesh.nodes[0].of_triangles)[0]
    
    # print(m.has_point(mesh.boundary_mesh.points[0])) #test if a triangle contains a point. 
    
    


    ##############################################
    
    # triangles = mesh.triangles
    # point = mesh.nodes[75]
    # triangle = mesh.triangles[0]
    
    # print(triangle)
    # print(mesh.adjacent_element(point, triangle))
    
    # point_test =  Point(100*np.random.rand(), 100*np.random.rand())
    
    # print(mesh.adjacent_element(point, triangle))
    # for i in range(4):
    #     point_test =  Point(100*np.random.rand(), 100*np.random.rand())
    #     print(point_test)
    #     l = mesh.find_triangle(point_test)
    #     figax = plt.subplots()
    #     fig, ax = figax
    #     ax.set_aspect('equal')
    #     ax.plot(point_test.x ,point_test.y ,'o')
    #     for e in l[:]:
    #         e.plot(figax)
        # print(len(l))
    # print(l)
    
    # print(triangle.nodes[75].array)
    # print(triangle.nodes[75]._of_triangles)
    # print(triangles[0])
    
    # figax = plt.subplots()
    # fig, ax = figax
    
    # ax.set_aspect('equal')
    # mesh.plot(figax)
    # for i in point._of_triangles:
    #     mesh.triangles[i].plot(figax)
    
    # ax.plot(point_test.x ,point_test.y ,'o')
    # l[-1].plot(figax)
    # triangle.plot(figax)
    # # mesh.adjacent_element(point, triangle).plot(figax)
    
    # for e in l[:]:
    #     e.plot(figax)
    # print(len(l))
    # mesh.triangles[2].plot(figax)
    
    # print(mesh.boundary_mesh.nodes[0].of_triangles[:])
    # print(mesh.boundary_mesh.nodes)
    # mesh.boundary_mesh.plot()
    
    Points =     np.array([[0.20671916, 0.84355497],
           [0.91861091, 0.34602815],
           [0.48841119, 0.10082727],
           [0.61174386, 0.38340907],
           [0.76590786, 0.5103548 ],
           [0.51841799, 0.96110308],
           [0.2968005 , 0.37151262],
           [0.18772123, 0.01236941],
           [0.08074127, 0.85970689],
           [0.7384403 , 0.11111075],
           [0.44130922, 0.47833904],
           [0.15830987, 0.84998003],
           [0.87993703, 0.51473797],
           [0.27408646, 0.44660783],
           [0.41423502, 0.80047642],
           [0.29607993, 0.02039138],
           [0.62878791, 0.57261865],
           [0.57983781, 0.41138362],
           [0.5999292 , 0.9851368 ],
           [0.26581912, 0.80140153],
           [0.28468588, 0.0539621 ],
           [0.25358821, 0.19047777],
           [0.32756395, 0.45241885],
           [0.1441643 , 0.70294208],
           [0.16561286, 0.33204815],
           [0.96393053, 0.3599832 ],
           [0.96022672, 0.92147057],
           [0.18841466, 0.95363051],
           [0.02430656, 0.40768573],
           [0.20455555, 0.89857115],
           [0.69984361, 0.33025325],
           [0.77951459, 0.08273857],
           [0.02293309, 0.52671757],
           [0.57766286, 0.66084439],
           [0.00164217, 0.89298429],
           [0.51547261, 0.96515755],
           [0.63979518, 0.76993268],
           [0.9856244 , 0.75909912],
           [0.2590976 , 0.71004905],
           [0.80249689, 0.70155283],
           [0.87048309, 0.76733138],
           [0.92274961, 0.97434948],
           [0.00221421, 0.37371452],
           [0.46948837, 0.08305444],
           [0.98146874, 0.23964044],
           [0.3989448 , 0.22148276],
           [0.81373248, 0.3635998 ],
           [0.5464565 , 0.81031423],
           [0.77085409, 0.06009617],
           [0.48493107, 0.44974356],
           [0.02911156, 0.81317053],
           [0.08652569, 0.26423837],
           [0.11145381, 0.06340098],
           [0.25124511, 0.24210767],
           [0.96491529, 0.08507046],
           [0.63176605, 0.80777744],
           [0.8166602 , 0.17025008],
           [0.566082  , 0.19534463],
           [0.63535621, 0.81464201],
           [0.81190239, 0.81028553],
           [0.92668262, 0.58937388],
           [0.91262676, 0.91473434],
           [0.82481072, 0.05982164],
           [0.09420273, 0.96499664],
           [0.36104842, 0.57097522],
           [0.03550903, 0.30251811],
           [0.54635835, 0.82570583],
           [0.79614272, 0.65941728],
           [0.0511428 , 0.98650144],
           [0.18866774, 0.10748953],
           [0.36547777, 0.58091853],
           [0.24429087, 0.47282816],
           [0.79508747, 0.65227079],
           [0.35209494, 0.24185905],
           [0.63887768, 0.03113006],
           [0.49341505, 0.54423235],
           [0.58349974, 0.36471018],
           [0.93929935, 0.89233328],
           [0.94354008, 0.45934559],
           [0.11169243, 0.4186541 ]]) *100
    
    
    # figax = plt.subplots()
    # fig, ax = figax
    
    border_data_ = [34, 68, 37, 7, 41, 74, 42, 44, 18, 52, 54, 26, 62]
    
    # ax.set_aspect('equal')
    # mesh.plot(figax)
    pp = []
    projected_pp = []
    l = []
    for i in border_data_:
        pp.append(Point(Points[i][0], Points[i][1]))
        # ax.plot(Points[i][0] ,Points[i][1] ,'x', color = 'blue')
        P = projection_on_boundary(mesh, pp[-1])
        projected_pp.append(P)
    #     ax.plot(P.x, P.y, 'o', color = 'r')
    #     mesh.find_triangle(pp[-1]).plot(figax)
        # l.append(mesh.find_triangle(pp[i])[-1])
        
    # BM = Boundary_Mesh(mesh)
    
    # p1 = Point(1.4, 3)
    # print(p1)
    # print(projection_on_boundary(mesh, p1))
    
    data_deformation = np.array([[-1.53832353e-02,  9.99081718e-03],
           [-1.98750614e+00,  5.24696468e-01],
           [-3.67360187e-01, -1.62272457e-01],
           [-4.56454351e-01,  1.01593240e-01],
           [-4.01721875e-01,  1.34832090e-01],
           [-3.10920145e-02,  2.33269716e-02],
           [-9.50867448e-02,  3.12158440e-03],
           [-4.74415040e-01, -2.07533164e-02],
           [-1.13585798e-02,  1.11292884e-02],
           [-1.10875848e+00, -4.92713078e-01],
           [-1.09994169e-01,  7.04703635e-02],
           [-1.44304605e-02,  6.67030271e-03],
           [-4.48186287e-01,  1.11803515e-01],
           [-4.17153393e-02,  5.80839051e-03],
           [-3.64337759e-02,  3.13544105e-02],
           [-4.67963660e-01, -8.39650320e-02],
           [-1.83969841e-01,  8.87730282e-02],
           [-2.81826996e-01,  7.74679287e-02],
           [-4.07890194e-02,  4.42572908e-02],
           [-2.12169835e-02,  1.52004908e-02],
           [-3.49313558e-01, -6.24852489e-02],
           [-1.67714809e-01, -3.18901304e-02],
           [-5.15984703e-02,  1.72399155e-03],
           [-1.96043628e-02,  1.08876373e-02],
           [-6.13726553e-02,  1.15423421e-02],
           [-1.73495720e+00,  1.30203924e-01],
           [-4.42051428e-02,  1.38892933e-01],
           [-2.74428554e-02,  8.72233877e-03],
           [-2.26203663e-02,  3.03893565e-02],
           [-2.80809124e-02,  1.10995699e-02],
           [-8.00844244e-01,  9.43391638e-03],
           [-1.04667502e+00, -3.85416888e-01],
           [-2.13672675e-02,  2.31889493e-02],
           [-8.63557981e-02,  6.27845807e-02],
           [-1.47494914e-02,  1.47338618e-02],
           [-3.28058454e-02,  2.18953729e-02],
           [-5.81529739e-02,  6.20041209e-02],
           [-1.54054111e-02,  1.80842040e-01],
           [-3.73061808e-02,  2.80475571e-02],
           [-4.00782988e-02,  1.11845306e-01],
           [-1.98901688e-02,  1.20727131e-01],
           [-4.41124420e-02,  1.23983073e-01],
           [-4.45049174e-02,  4.34870305e-02],
           [-3.81449792e-01, -1.42254663e-01],
           [-5.66061544e+00,  4.08592874e-01],
           [-2.57121543e-01, -5.62805212e-02],
           [-1.28928762e+00,  2.19243213e-01],
           [-3.86231099e-02,  4.39568913e-02],
           [-7.80147683e-01, -3.83156176e-01],
           [-2.12218872e-01,  5.52330261e-02],
           [-1.07319626e-02,  1.64397821e-02],
           [-9.35205671e-02,  4.50174446e-02],
           [-2.55331570e-01,  1.39948357e-02],
           [-1.09481848e-01, -1.52735247e-02],
           [-1.52715990e+00,  1.26000524e+00],
           [-2.81091620e-02,  5.22783931e-02],
           [-2.38713950e+00, -1.18019438e-01],
           [-4.77768896e-01, -1.54069566e-01],
           [-2.49275165e-02,  5.07251085e-02],
           [-2.64424874e-02,  9.92386301e-02],
           [-1.19484521e-01,  1.11537669e-01],
           [-3.93738342e-02,  1.24082433e-01],
           [-8.24077043e-01,  5.06108375e-02],
           [-2.16933270e-02,  1.08133189e-02],
           [-8.28831919e-02,  3.16307157e-02],
           [-6.55763327e-02,  4.35646083e-02],
           [-2.78736760e-02,  3.49562952e-02],
           [-5.58387915e-02,  1.03053399e-01],
           [-3.17178883e-02,  7.60831719e-03],
           [-2.28613388e-01, -5.45467241e-03],
           [-7.69715961e-02,  2.89974168e-02],
           [-3.50464812e-02,  7.88956729e-03],
           [-1.24999348e-01,  1.23738316e-01],
           [-1.62006036e-01,  5.94973139e-03],
           [-5.50461532e-01, -5.62905263e-01],
           [-1.13996526e-01,  7.64117598e-02],
           [-4.47419810e-01,  1.53165748e-02],
           [-3.31808457e-02,  1.34484158e-01],
           [-8.25013814e-01,  8.84844258e-02],
           [-2.64316484e-02,  1.85501969e-02]])*100
    
    
    Points[border_data_], data_deformation[border_data_] # Point d'interpolation et donnees
    
    
    # Sample data points (replace with your own data)
    # Extract coordinates and values
    
    # x = Points[border_data_][:,0] 
    # y = Points[border_data_][:,1] 
    
    
    x = np.array([projected_pp[i].x for i in range(len(border_data_))])
    y = np.array([projected_pp[i].y for i in range(len(border_data_))])
    values_directionx = data_deformation[border_data_][:,0]
    values_directiony = data_deformation[border_data_][:,1]
    
    interp_function_x = LinearNDInterpolator(list(zip(x, y)), values_directionx)
    interp_function_y = LinearNDInterpolator(list(zip(x, y)), values_directiony)
    
    def f_x(x_eval, y_eval):
        # return rbf_x(x_eval, y_eval)
        return interp_function_x(x_eval, y_eval)
    
    def f_y(x_eval, y_eval):
        # return rbf_x(x_eval, y_eval)
        return interp_function_y(x_eval, y_eval)
    
    def Dirichlet_function(x_eval, y_eval):
        return np.array([f_x(x_eval, y_eval), f_y(x_eval, y_eval)])
    
    
    # plt.figure()
    grid_x, grid_y = np.meshgrid(np.linspace(0, 100, num=100),
                              np.linspace(0, 100, num=100))
    grid_values_x = f_x(grid_x, grid_y)

    # Plot the interpolated values
    # plt.imshow(grid_values_x, origin='lower', aspect='auto')
    # plt.colorbar()
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.title('Interpolated Function using LinearNDInterpolator')
    #plt.show()


    # plt.figure()
    # grid_values_y = f_y(grid_x, grid_y)
    # plt.imshow(grid_values_y, origin='lower', aspect='auto')
    # plt.colorbar()
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.title('Interpolated Function using LinearNDInterpolator')
    # plt.show()
    
    ### plot the boundary dirichlet function
    
    # plt.figure()
    # u_vals, v_vals = Dirichlet_function(grid_x, grid_y)
    # plt.quiver(grid_x, grid_y, u_vals, v_vals)
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.title('Vector Field Visualization')
    # plt.show()
    
    # Taking into account the Dirichlet condition 
    
    
    logger = log.Log('griffith_solver.log', level=log.INFO, console_output=True)
    logger.log_description(mesh_file=mesh,args=None)
    log_queue = logger._log_queue
    
    T = problem_data.lame_tensor_ice #Lame tensor
    
    XY_data = np.stack((x,y), axis = -1) 

    boundary_data = np.hstack( (XY_data, data_deformation[border_data_]) )
    
    
    boundary_displacement = problem_data.Boundary_Displacement_by_percussion(boundary_data = boundary_data)
    # boundary_displacement = problem_data.Constant_Displacement_On_Y(traction_coefficient= -100.)
    
        
    physical_data = problem_data.Physical_Data(T, 1., boundary_displacement, initial_fracture=None)
    
    
    ### Solution without fracture
    classical_solution = solver.Classical_Solution(mesh=mesh, physical_data=physical_data)
    print(classical_solution.energy)
    classical_solution.plot_displacement()
    classical_solution.plot_energy()
    
    ### Solution with fracture 
    boundary_point = [float(20), float(0)]
    time_discretization = None
    fracture_discretization = problem_data.Fracture_Discretization(angular_step = np.pi/4., boundary_point= boundary_point, lengh_step = 100 )
    discretization_data =  problem_data.Discretization_Data(mesh, time_discretization, fracture_discretization, tip_enrichement=False)
    
    solution = solver.smart_time_solver(discretization_data, physical_data, log_queue)
    # solution.plot_displacement()
    # solution.plot_energy()
    
    
    
    
    
    