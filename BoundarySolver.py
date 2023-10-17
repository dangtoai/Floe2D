#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 17:18:33 2023

@author: phandangtoai
"""

import csv
import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append("/Users/phandangtoai/Documents/Floe2D/Griffith-master_Dimitri")
# import griffith.geometry as geo
from griffith.mesh import Mesh
# from griffith.geometry import Point, dist
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
from griffith import solver, log, problem_data
from numpy.linalg import norm

def line_coefficient(point1, point2):
    """
    compute the coefficient of line contains P1, P2
    return A,B,C of Ax+By=C
    """
    x1, y1 = point1
    x2, y2 = point2
    A = (y2-y1)/(x2-x1)
    B = y1 - A*x1
    return np.array([A,B])
    
def test_line(t, A, B):
    """
    return the line of equation y = At+B
    """
    return A*t + B

def intersection_line(A1, B1, A2, B2):
    """
    find the intersection of 2 line 
    y = A1*t +B1 and y = A2*t+B2"""
    M = np.array([[-A1, 1.], [-A2, 1.]])
    B = np.array([B1, B2])
    return np.linalg.solve(M, B) 

def boundary_edges_index(tri: Delaunay):
    """
    return the (i,j) index of all edges at the boundary 
    """
    All_edges = []
    for triangles in tri.simplices:
        for i in range(3):
            for j in range(i+1, 3):
                All_edges.append(
                    (min(triangles[i], triangles[j]), max(triangles[i], triangles[j])))
    S = set(All_edges)
    Border = ([e for e in S if All_edges.count(e) == 1])
    return Border

def boundary_triangles_index(tri: Delaunay):
    """
    return the (i,j,k) index of all triangles at the boundary 
    """
    boundary_edges = boundary_edges_index(tri)
    boundary_triangles = []
    for edge in boundary_edges:
        for triangle in tri.simplices:
            if set(edge).issubset(triangle):
                boundary_triangles.append(triangle)
    return boundary_triangles

def boundary_nodes_index(tri: Delaunay):
    """
    return the (i) index of all nodes at the boundary 
    """
    boundary_nodes = []
    boundary_edges = boundary_edges_index
    for i,j in boundary_edges:
        boundary_nodes.append(i)
        boundary_nodes.append(j)
    return list(set(boundary_nodes))

def cones_list(tri: Delaunay):
    """
    return the (i,j,k) index of all triangles at the boundary 
    the node k is in the interior of the mesh Delaunay
    such that ki and kj are 2 lines of the cone
    k is the head of the cone
    """
    boundary_edges = boundary_edges_index(tri)
    boundary_triangles = boundary_triangles_index(tri)
    triangles_out = []
    for triangle in boundary_triangles:
        for edge in boundary_edges:
            if set(edge).issubset(triangle):
                common_el = np.intersect1d(triangle, edge)
                not_common = np.setdiff1d(triangle, common_el)
                # print(common_el, not_common)
                triangle = np.append(not_common, common_el)
                triangles_out.append(triangle)
    return triangles_out

def inside_cone(head, base1, base2, point):
    """
    Verify if point is inside of the cone 
    
    """
    direction1 = base1 - head
    direction2 = base2 - head
    vector_point = point - head
    
    direction1 = direction1/norm(direction1)
    direction2 = direction2/norm(direction2)
    vector_point = vector_point/norm(vector_point)
    
    angle1 = np.arccos(np.dot(direction1, vector_point))
    angle2 = np.arccos(np.dot(direction2, vector_point))
    cone_angle = np.arccos(np.dot(direction1, direction2))
    
    return angle1<cone_angle and angle2<cone_angle

def P1_coefficient(Points, data):
    """
    if datas contains the value at points,
    return the affine function f that
    f(point[i]) = datas[i]
    f (x,y) = Ax+By+C
    """
    
    Matrix = np.array([[Points[0][0], Points[0][1], 1],
                      [Points[1][0], Points[1][1], 1],
                      [Points[2][0], Points[2][1], 1]])
    B_ = np.array(data)
    A, B, C = np.linalg.solve(Matrix, B_)

    return A,B,C
    

mesh = Mesh("file.msh")
data = []
with open('boundary_data.csv', 'r', encoding='utf-8') as csv_file:
    lines = csv_file.readlines()
    lines_to_read = lines[4:]
    csv_reader = csv.reader(csv_file)
    csv_reader = csv.reader(lines_to_read)
    for row in csv_reader:
        # print(row)
        data.append([float(item) for item in row[0].split()])

data = np.array(data)

xdata, ydata = data[:, 0], data[:, 1]
z1data = data[:,2]
z2data = data[:,3]
Points = np.array(list(zip(xdata, ydata)))

tri = Delaunay(Points)
# print(cones_list(tri))
# print(boundary_edge_index(tri))
# print(boundary_triangles_index(tri))
plt.triplot(Points[:,0], Points[:,1], cones_list(tri))
for i in range(27):
    plt.plot(Points[i][0], Points[i][1])
    plt.text(Points[i][0], Points[i][1], str(i), color = 'red')
# plt.plot(1009,123,'o')
# print(inside_cone(Points[14], Points[3], Points[7], np.array([1011,120])))
plt.plot(Points[:,0], Points[:,1], 'o')
plt.show()


interp_function_x = LinearNDInterpolator(list(zip(xdata, ydata)), z1data)
# interp_function_x = LinearNDInterpolator(list(zip(xdata, ydata)), z1data_new)

interp_function_y = LinearNDInterpolator(list(zip(xdata, ydata)), z2data)
def f_x(x_eval, y_eval):
    interpolated_value = interp_function_x(x_eval, y_eval)
    Cones = cones_list(tri) 
    if np.isnan(interpolated_value):
        for i,j,k in Cones:
            if inside_cone(Points[i], Points[j], Points[k], np.array([x_eval, y_eval])):
                # print(i,j,k)
                P_ = np.array([Points[i], Points[j], Points[k]])
                data_ = np.array([z1data[i], z1data[j], z1data[k]])
                A, B, C = P1_coefficient(P_, data_)
                interpolated_value = A*x_eval + B*y_eval + C
    return interpolated_value

def f_y(x_eval, y_eval):
    interpolated_value = interp_function_y(x_eval, y_eval)
    Cones = cones_list(tri) 
    if np.isnan(interpolated_value):
        for i,j,k in Cones:
            if inside_cone(Points[i], Points[j], Points[k], np.array([x_eval, y_eval])):
                # print(i,j,k)
                P_ = np.array([Points[i], Points[j], Points[k]])
                data_ = np.array([z1data[i], z1data[j], z1data[k]])
                A, B, C = P1_coefficient(P_, data_)
                interpolated_value = A*x_eval + B*y_eval + C
    return interpolated_value

grid_x, grid_y = np.meshgrid(np.linspace(min(xdata)-1, max(xdata)+1, num=10),
                            np.linspace(min(ydata)-1, max(ydata)+2, num=10))

grid_values_x = np.vectorize(f_x)(grid_x, grid_y)

figax = plt.subplots()
fig, ax = figax
# ax.plot(x_reordered, y_reordered, 'x')
mesh.boundary_mesh.plot(figax)
for i in range(27):
    ax.plot(Points[i][0], Points[i][1])
    ax.text(Points[i][0], Points[i][1], str(i), color = 'red')
contour_plot = ax.contourf(grid_x, grid_y, grid_values_x, cmap = 'viridis')
# ax.set_xlim(min(xdata)-1, max(xdata)+1)
# ax.set_ylim(min(ydata)-1, max(ydata)+2)
cbar = plt.colorbar(contour_plot)
cbar.set_label('Data Value')

logger = log.Log('griffith_solver.log', level=log.INFO, console_output=True)
logger.log_description(mesh_file=mesh,args=None)
log_queue = logger._log_queue

T = problem_data.lame_tensor_ice #Lame tensor
boundary_displacement = problem_data.Boundary_Displacement_by_percussion(boundary_data = data)
physical_data = problem_data.Physical_Data(T, 1., boundary_displacement, initial_fracture=None)
classical_solution = solver.Classical_Solution(mesh=mesh, physical_data=physical_data)
print(classical_solution.energy)
classical_solution.plot_displacement()
classical_solution.plot_energy()


