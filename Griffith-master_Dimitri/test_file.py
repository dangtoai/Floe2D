#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 10:52:55 2023

@author: phandangtoai
"""
import time
# import argparse
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
from Func import * 


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
    
    
    figax = plt.subplots()
    fig, ax = figax
    
    border_data_ = [34, 68, 37, 7, 41, 74, 42, 44, 18, 52, 54, 26, 62]
    
    ax.set_aspect('equal')
    mesh.plot(figax)
    pp = []
    projected_pp = []
    l = []
    for i in border_data_:
        pp.append(Point(Points[i][0], Points[i][1]))
        ax.plot(Points[i][0] ,Points[i][1] ,'x', color = 'blue')
        P = projection_on_boundary(mesh, pp[-1])
        projected_pp.append(P)
        ax.plot(P.x, P.y, 'o', color = 'r')
        mesh.find_triangle(pp[-1]).plot(figax)
        # l.append(mesh.find_triangle(pp[i])[-1])
        
    # BM = Boundary_Mesh(mesh)
    
    # p1 = Point(1.4, 3)
    # print(p1)
    # print(projection_on_boundary(mesh, p1))
    
    data_deformation = np.array([[-3.41564832e-11,  2.57830424e-11],
           [-1.85632055e-05,  1.60722770e-05],
           [-1.68090776e-07,  5.46621101e-09],
           [-5.90340877e-07,  3.30167871e-07],
           [-1.73471728e-07,  2.98723199e-07],
           [-4.26390578e-09, -9.02010577e-10],
           [-9.10353676e-09, -9.09289866e-11],
           [-2.84013588e-07, -1.65241980e-08],
           [-2.36965031e-11,  3.13005177e-11],
           [-6.91180390e-06, -5.22930459e-06],
           [-1.28927169e-08,  1.10020130e-08],
           [-6.86573853e-11, -1.23645538e-11],
           [-5.31102186e-07,  2.64411533e-07],
           [-3.66286779e-10,  2.19537888e-10],
           [-2.27156699e-09,  1.03928866e-09],
           [-3.07803054e-07, -1.78753714e-08],
           [-3.24440611e-08,  7.07586515e-08],
           [-9.96135235e-08,  1.06302033e-07],
           [-3.16967317e-08,  7.67280417e-09],
           [-2.96822678e-10,  7.40313366e-11],
           [-7.33792120e-08, -1.67906928e-08],
           [-2.56531951e-08, -5.14857504e-09],
           [-1.39390521e-09, -5.69953196e-10],
           [-6.19866936e-11,  2.11769491e-11],
           [-1.19518179e-09, -2.75484635e-11],
           [-8.76883123e-06,  6.45705240e-06],
           [-1.68036509e-07,  1.26565855e-06],
           [-9.82863652e-10,  1.02475806e-10],
           [ 5.36953641e-11,  6.53346655e-10],
           [-7.62802765e-10,  4.30283587e-11],
           [-2.60321540e-06,  1.55478140e-06],
           [-3.50276836e-06, -5.16977010e-06],
           [ 1.24188888e-11,  3.70704689e-10],
           [-1.48122060e-08,  1.07737762e-08],
           [-1.57431613e-10,  1.62251879e-10],
           [-4.01763789e-09, -7.25233096e-10],
           [-2.97540446e-08,  2.03909694e-08],
           [ 1.02407351e-06,  5.99424169e-06],
           [-5.91996785e-10,  4.31133573e-10],
           [-7.16648130e-08,  1.30981994e-07],
           [ 9.35105804e-08,  3.66702893e-07],
           [-1.22356978e-07,  1.91494036e-07],
           [-1.64501331e-09,  3.08136194e-09],
           [-1.88512149e-07,  4.05692886e-08],
           [-1.62956887e-04, -4.23574112e-06],
           [-1.65706386e-07,  1.12943992e-08],
           [-4.65129630e-06, -3.46665979e-09],
           [-5.04706887e-09,  3.51773632e-09],
           [-1.06312720e-06, -1.68442299e-06],
           [-8.58842049e-08,  5.73161460e-08],
           [-1.04872465e-11,  6.68903821e-11],
           [-3.36958308e-09,  4.00421735e-09],
           [-3.26280235e-08,  1.49091369e-08],
           [-4.18891494e-09, -1.01954631e-09],
           [-9.01131036e-06, -3.64816986e-06],
           [-1.22843236e-09,  4.79960827e-09],
           [-3.69792557e-05, -6.81259045e-06],
           [-7.61860042e-07,  1.62195023e-07],
           [ 1.57917002e-09,  5.38664213e-09],
           [ 1.41696419e-08,  6.18326643e-08],
           [ 4.62336754e-07,  1.60572083e-06],
           [-1.20418035e-07,  1.60805077e-07],
           [-2.05413696e-06, -2.22421480e-06],
           [-5.06905268e-10,  2.08218776e-10],
           [-3.75964970e-09,  1.66969472e-09],
           [-1.66255316e-09,  3.42459711e-09],
           [-4.84872253e-10,  6.26695362e-11],
           [-6.28926858e-08,  6.58440377e-08],
           [-4.22673813e-09, -6.28026520e-11],
           [-1.19261106e-08, -2.00837065e-09],
           [-3.78260062e-09,  4.15184997e-10],
           [-3.03252284e-10, -6.68750610e-11],
           [-6.63710703e-08,  1.19405506e-07],
           [-2.53531429e-08,  1.05397728e-08],
           [-1.28482418e-06, -1.42335666e-06],
           [-1.27142917e-08,  2.40307161e-08],
           [-6.04920633e-07,  1.76406808e-07],
           [-3.72702007e-07,  1.12824000e-06],
           [-4.18852618e-07,  4.09280536e-06],
           [ 5.95624799e-11,  2.86838220e-10]])
    
    
    Points[border_data_], data_deformation[border_data_] # Point d'interpolation et donnees
    
    
    # Sample data points (replace with your own data)
    # Extract coordinates and values
    
    x = Points[border_data_][:,0] 
    # x = np.array([projected_pp[i].x for i in range(len(border_data_))])
    y = Points[border_data_][:,1] 
    # y = np.array([projected_pp[i].y for i in range(len(border_data_))])
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
    
    plt.figure()
    grid_x, grid_y = np.meshgrid(np.linspace(0, 100, num=100),
                             np.linspace(0, 100, num=100))
    grid_values_x = f_x(grid_x, grid_y)

    # Plot the interpolated values
    plt.imshow(grid_values_x, origin='lower', aspect='auto')
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Interpolated Function using Rbf')
    plt.show()

    # Evaluate the interpolated function at specific points
    result = f_x(0.5, 0.5)
    print(result)
    
    
    
    
    