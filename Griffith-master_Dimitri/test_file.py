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
import numpy as np
# from griffith.mesh import 
from Func import * 


if __name__ == '__main__':
    
    #define Lame tensor. 
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
    
    
    # if the oriented area is positive, then the point belongs to the interior of the triangle
    # print(" oriented area of triangle ",mesh.triangles[0].oriented_area())
    # print(" All oriented area of triangle ",mesh.triangles[0].all_oriented_area())
    
    
    # if not, one of the oriented area is negative 
    
    # triangle_test = Triangle(mesh.nodes[186],mesh.nodes[75], mesh.nodes[0], None)
    # print(triangle_test.oriented_area())
    # triangle_test = Triangle(mesh.nodes[75], mesh.nodes[268],mesh.nodes[0], None)
    # print(triangle_test.oriented_area())
    # triangle_test = Triangle(mesh.nodes[186], mesh.nodes[268], mesh.nodes[0], None)
    # print(triangle_test.oriented_area())
    
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
    l = []
    for i in border_data_:
        pp.append(Point(Points[i][0], Points[i][1]))
        ax.plot(Points[i][0] ,Points[i][1] ,'o')
        mesh.find_triangle(pp[-1])[-1].plot(figax)
        
        # l.append(mesh.find_triangle(pp[i])[-1])
        
        
    # for i in range(37):
        # l.append(mesh.find_triangle(pp[i])[-1])
        # print(mesh.find_triangle(pp[i])[-1])
    

    # ax.plot(point_test.x ,point_test.y ,'o')
    # point_test.plot(figax, 'o')
    # l[-1].plot(figax)
    # triangle.plot(figax)
    # mesh.adjacent_element(point, triangle).plot(figax)
    
    # for e in l[:]:
        # e.plot(figax)
    
    # BM = Boundary_Mesh(mesh)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    