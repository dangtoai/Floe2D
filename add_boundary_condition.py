#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 16:37:49 2023

@author: phandangtoai
"""

import csv
import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append("/Users/phandangtoai/Documents/Floe2D/Griffith-master_Dimitri")
from griffith import solver, log, problem_data
from griffith.mesh import Mesh, projection_on_boundary
from griffith.geometry import Point

mesh = Mesh("file.msh")
T = problem_data.lame_tensor_ice #Lame tensor

data = []
with open('boundary_data.csv', 'r', encoding='utf-8') as csv_file:
    lines = csv_file.readlines()
    lines_to_read = lines[1:-2]
    csv_reader = csv.reader(csv_file)
    csv_reader = csv.reader(lines_to_read)
    for row in csv_reader:
        print(row)
        data.append( [float(item) for item in row[0].split()] )

data = np.array(data)
print(data)
boundary_displacement = problem_data.Boundary_Displacement_by_percussion(boundary_data = data)
physical_data = problem_data.Physical_Data(T, 1., boundary_displacement, initial_fracture=None)
classical_solution = solver.Classical_Solution(mesh=mesh, physical_data=physical_data)
# mesh.plot()

# print(mesh.is_node_boundary(Node(1036.5, 94.0, None)))

# logger = log.Log('griffith_solver.log', level=log.INFO, console_output=True)
# logger.log_description(mesh_file=mesh,args=None)
# log_queue = logger._log_queue

projected_points = []

figax = plt.subplots()
fig, ax = figax

for i in range(len(data)):
    P = projection_on_boundary(mesh, Point(data[:,0][i], data[:,1][i]))
    projected_points.append(P)
    ax.plot(P.x, P.y, 'x', color = 'blue')
mesh.boundary_mesh.plot(figax)

#modify the boundary data

data[:,0] = np.array([P.x for P in projected_points ])
data[:,1] = np.array([P.y for P in projected_points ])
print(data)
boundary_displacement = problem_data.Boundary_Displacement_by_percussion(boundary_data = data)


physical_data = problem_data.Physical_Data(T, 1., boundary_displacement, initial_fracture=None)
classical_solution = solver.Classical_Solution(mesh=mesh, physical_data=physical_data)
# print(classical_solution.energy)
# classical_solution.plot_displacement()
# classical_solution.plot_energy()
