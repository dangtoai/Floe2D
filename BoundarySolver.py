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
import griffith.geometry as geo
from griffith.mesh import Mesh
from griffith.geometry import Point, dist
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
from griffith import solver, log, problem_data, fracture_iterators
from numpy.linalg import norm
from shapely.geometry import Polygon, LineString

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
    center = sum(tri.points)/tri.npoints
    for triangle in boundary_triangles:
        for edge in boundary_edges:
            if set(edge).issubset(triangle):
                common_el = np.intersect1d(triangle, edge)
                if np.cross(tri.points[common_el[0]] - center, tri.points[common_el[1]] - center)<0:
                    common_el = common_el[::-1]
                not_common = np.setdiff1d(triangle, common_el)
                triangle = np.append(not_common, common_el)
                triangles_out.append(triangle)
    return triangles_out

def inside_cone(head, base1, base2, point):
    # Calculate vectors from the  head to the rays and the point
    vector1 = base1 -  head
    vector2 = base2 -  head
    vector_point = point -  head

    # Calculate the angles between the vectors using the dot product
    angle1 = np.arctan2(vector1[1], vector1[0])
    angle2 = np.arctan2(vector2[1], vector2[0])
    point_angle = np.arctan2(vector_point[1], vector_point[0])

    # Normalize the angles to the range [0, 2*pi)
    angle1 = (angle1 + 2 * np.pi) % (2 * np.pi)
    angle2 = (angle2 + 2 * np.pi) % (2 * np.pi)
    point_angle = (point_angle + 2 * np.pi) % (2 * np.pi)

    # Check if the point angle is between the ray angles
    if angle1 <= angle2:
        return angle1 <= point_angle <= angle2
    else:
        return angle1 <= point_angle or point_angle <= angle2

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
    lines_to_read = lines[3:]
    csv_reader = csv.reader(csv_file)
    csv_reader = csv.reader(lines_to_read)
    for row in csv_reader:
        data.append([float(item) for item in row[0].split()])

data = np.array(data)

xdata, ydata = data[:, 0], data[:, 1]
z1data = data[:,2] 
z2data = data[:,3] 
Points = np.array(list(zip(xdata, ydata)))

tri = Delaunay(Points)


interp_function_x = LinearNDInterpolator(list(zip(xdata, ydata)), z1data)
interp_function_y = LinearNDInterpolator(list(zip(xdata, ydata)), z2data)

def f_x(x_eval, y_eval):
    interpolated_value = interp_function_x(x_eval, y_eval)
    Cones = cones_list(tri) 
    if np.isnan(interpolated_value):
        # print(x_eval, y_eval)
        interpolated_value = 0.
        l = 0                       # triangle count
        L = []                      # stock the list of triangles 
        for i,j,k in Cones:
            if inside_cone(Points[i], Points[j], Points[k], np.array([x_eval, y_eval])):
                l += 1
                L.append([i,j,k])
        # print(L)
        
        if l == 1:
            i0, j0, k0 = L[0]
            P_ = np.array([Points[i0], Points[j0], Points[k0]])
            data_ = np.array([z1data[i0], z1data[j0], z1data[k0]])
            A, B, C = P1_coefficient(P_, data_)
            interpolated_value += A*x_eval + B*y_eval + C
            # print(interpolated_value)
            return interpolated_value
        
        if l>=2: 
            i = 0
            while i < len(L):
                j = i+1
                while j < len(L):
                    intersection = np.intersect1d(L[i], L[j])
                    if len(intersection) == 0:
                        del L[j]
                    else: j += 1
                i += 1

            # create new data from the interpolation with the mesh's boundary
            # find the intersection between the cone and the mesh's boundary 
            if len(L) >= 2:
                Origine1 = Points[L[0][0]]
                # print(Points[L[0][0] )
                # print(L[0][0], Origine1)
                direction1 = (Points[L[0][-1]] - Points[L[0][0]]) * 1000
                data_coord1 = LineString([Origine1, Origine1 + direction1])
                data_coord1 = polygon.intersection(data_coord1)
                data_coord1 = np.array(list(data_coord1.coords)[-1])
                
                #compute the P1 approximation on the old network
                i1, j1, k1 = L[0]
                P1 = np.array([Points[i1], Points[j1], Points[k1]])
                data_ = np.array([z1data[i1], z1data[j1], z1data[k1]])
                A, B, C = P1_coefficient(P1, data_)
                new_data1   = A*data_coord1[0] + B*data_coord1[1] + C

                Origine2 = Points[L[1][0]]
                direction2 = (Points[L[1][-2]] - Points[L[1][0]]) * 1000
                data_coord2 = LineString([Origine2, Origine2 + direction2])
                data_coord2 = polygon.intersection(data_coord2)
                data_coord2 = np.array(list(data_coord2.coords)[-1])
                
                i_, j_, k_ = L[1]
                P2 = np.array([Points[i_], Points[j_], Points[k_]])
                data_ = np.array([z1data[i_], z1data[j_], z1data[k_]])
                A, B, C = P1_coefficient(P2, data_)
                new_data2   = A*data_coord2[0] + B*data_coord2[1] + C
                
                common_el = np.intersect1d(L[0], L[1])[0]
                P_ = np.array([Points[common_el], data_coord1, data_coord2 ])
                data_ = np.array([z1data[common_el], new_data1 , new_data2])
                
                # print(P_)
                # print(data_)
                A, B, C = P1_coefficient(P_, data_)
                # print(A,B,C)
                interpolated_value = A*x_eval + B*y_eval + C

    return interpolated_value

poly_coord = [p.array for p in mesh.boundary_mesh.points]
polygon = Polygon(poly_coord)

def f_y(x_eval, y_eval):
    interpolated_value = interp_function_y(x_eval, y_eval)
    Cones = cones_list(tri) 
    if np.isnan(interpolated_value):
        # print(x_eval, y_eval)
        interpolated_value = 0.
        l = 0                       # triangle count
        L = []                      # stock the list of triangles 
        for i,j,k in Cones:
            if inside_cone(Points[i], Points[j], Points[k], np.array([x_eval, y_eval])):
                l += 1
                L.append([i,j,k])
        # print(L)
        
        if l == 1:
            i0, j0, k0 = L[0]
            P_ = np.array([Points[i0], Points[j0], Points[k0]])
            data_ = np.array([z2data[i0], z2data[j0], z2data[k0]])
            A, B, C = P1_coefficient(P_, data_)
            interpolated_value += A*x_eval + B*y_eval + C
            # print(interpolated_value)
            return interpolated_value
        
        if l>=2: 
            i = 0
            while i < len(L):
                j = i+1
                while j < len(L):
                    intersection = np.intersect1d(L[i], L[j])
                    if len(intersection) == 0:
                        del L[j]
                    else: j += 1
                i += 1

            # create new data from the interpolation with the mesh's boundary
            # find the intersection between the cone and the mesh's boundary 
            if len(L) >= 2:
                Origine1 = Points[L[0][0]]
                # print(Points[L[0][0] )
                # print(L[0][0], Origine1)
                direction1 = (Points[L[0][-1]] - Points[L[0][0]]) * 1000
                data_coord1 = LineString([Origine1, Origine1 + direction1])
                data_coord1 = polygon.intersection(data_coord1)
                data_coord1 = np.array(list(data_coord1.coords)[-1])
                
                #compute the P1 approximation on the old network
                i1, j1, k1 = L[0]
                P1 = np.array([Points[i1], Points[j1], Points[k1]])
                data_ = np.array([z2data[i1], z2data[j1], z2data[k1]])
                A, B, C = P1_coefficient(P1, data_)
                new_data1   = A*data_coord1[0] + B*data_coord1[1] + C

                Origine2 = Points[L[1][0]]
                direction2 = (Points[L[1][-2]] - Points[L[1][0]]) * 1000
                data_coord2 = LineString([Origine2, Origine2 + direction2])
                data_coord2 = polygon.intersection(data_coord2)
                data_coord2 = np.array(list(data_coord2.coords)[-1])
                
                i_, j_, k_ = L[1]
                P2 = np.array([Points[i_], Points[j_], Points[k_]])
                data_ = np.array([z2data[i_], z2data[j_], z2data[k_]])
                A, B, C = P1_coefficient(P2, data_)
                new_data2   = A*data_coord2[0] + B*data_coord2[1] + C
                
                common_el = np.intersect1d(L[0], L[1])[0]
                P_ = np.array([Points[common_el], data_coord1, data_coord2 ])
                data_ = np.array([z2data[common_el], new_data1 , new_data2])
                
                # print(P_)
                # print(data_)
                A, B, C = P1_coefficient(P_, data_)
                # print(A,B,C)
                interpolated_value = A*x_eval + B*y_eval + C

    return interpolated_value

# Define the polygon vertices as a list of NumPy arrays

# mesh.boundary_mesh.plot()
# plt.triplot(Points[:,0], Points[:,1], cones_list(tri))
# for i in range(data.shape[0]):
#     plt.plot(Points[i][0], Points[i][1])
#     plt.text(Points[i][0], Points[i][1], str(i), color = 'red')


# def normfxy(x_eval, y_eval):
#     return norm(np.array([f_x(x_eval, y_eval), f_y(x_eval, y_eval)]))

def Dirichlet(x, y):
    return np.array([f_x(x, y), f_y(x, y)])

grid_x, grid_y = np.meshgrid(np.linspace(min(xdata)-1, max(xdata)+1, num = 100),
                            np.linspace(min(ydata)-1, max(ydata)+2, num = 100))

grid_values_x = np.vectorize(f_x)(grid_x, grid_y)
grid_values_y = np.vectorize(f_y)(grid_x, grid_y)
# grid_values_f = np.vectorize(normfxy)(grid_x, grid_y)

figax = plt.subplots()
fig, ax = figax
mesh.boundary_mesh.plot(figax)
for i in range(data.shape[0]):
    ax.plot(Points[i][0], Points[i][1], 'x')
    # ax.text(Points[i][0], Points[i][1], str(i), color = 'red')
contour_plot = ax.contourf(grid_x, grid_y, grid_values_x, cmap = 'viridis')
ax.set_xlim(min(xdata)-1, max(xdata)+1)
ax.set_ylim(min(ydata)-1, max(ydata)+2)
cbar = plt.colorbar(contour_plot)
cbar.set_label('Data Value')

figax = plt.subplots()
fig, ax = figax
mesh.boundary_mesh.plot(figax)
for i in range(data.shape[0]):
    ax.plot(Points[i][0], Points[i][1], 'x')
    # ax.text(Points[i][0], Points[i][1], str(i), color = 'red')
contour_plot = ax.contourf(grid_x, grid_y, grid_values_y, cmap = 'viridis')
ax.set_xlim(min(xdata)-1, max(xdata)+1)
ax.set_ylim(min(ydata)-1, max(ydata)+2)
cbar = plt.colorbar(contour_plot)
cbar.set_label('Data Value')

figax = plt.subplots()
fig, ax = figax
mesh.boundary_mesh.plot(figax)
ax.quiver(grid_x, grid_y, grid_values_x, grid_values_y , scale = 1.)

# figax = plt.subplots()
# fig, ax = figax
# mesh.boundary_mesh.plot(figax)
# for i in range(data.shape[0]):
    # ax.plot(Points[i][0], Points[i][1], 'x')
    # ax.text(Points[i][0], Points[i][1], str(i), color = 'red')
# contour_plot = ax.contourf(grid_x, grid_y, grid_values_f, cmap = 'viridis')
# ax.set_xlim(min(xdata)-1, max(xdata)+1)
# ax.set_ylim(min(ydata)-1, max(ydata)+2)
# cbar = plt.colorbar(contour_plot)
# cbar.set_label('Data Value')

# logger = log.Log('griffith_solver.log', level=log.INFO, console_output=True)
# logger.log_description(mesh_file=mesh,args=None)
# log_queue = logger._log_queue

T = problem_data.lame_tensor_ice #Lame tensor

# # boundary_displacement = problem_data.Constant_Displacement_On_Y(abscissa_mid= 1005.)
boundary_displacement = problem_data.Boundary_Displacement_by_percussion(boundary_data = data, Dirichlet_func = Dirichlet)
# # print(boundary_displacement.collision_point())
physical_data = problem_data.Physical_Data(T, problem_data.Constant_Toughness(10.), boundary_displacement, initial_fracture=None)

# print('start computation of classical energy')
# classical_solution = solver.Classical_Solution(mesh=mesh, physical_data=physical_data)
# # classical_solution.field.boundary_elements
# # solution = solver.smart_time_solver(discretization_data, physical_data, log_queue, args.number_processes)
# # Test = solver.Imposed_Fracture_Solution(mesh=mesh, physical_data=physical_data, fracture = None)

# # print(classical_solution.energy)
# classical_solution.plot_displacement()
# classical_solution.plot_energy()

# print('start computation of fractures')
# boundary_point = [1011.299987792969, 104.0]
# time_discretization = None
# fracture_discretization = problem_data.Fracture_Discretization(angular_step = np.pi/5., boundary_point= boundary_point, lengh_step = 100 )
# discretization_data =  problem_data.Discretization_Data(mesh, time_discretization, fracture_discretization, tip_enrichement=False)

# find fracture from every boundary admissible points. 
# f = fracture_iterators.Admissible_Boundary_Point(fracture_discretization, mesh.boundary_mesh)
# fig, ax = mesh.plot()
# for pointsegment in f:
#     point, segment = pointsegment
#     boundary_point = [point.x, point.y]
# # #     print(point)
#     fracture_discretization = problem_data.Fracture_Discretization(angular_step = np.pi/5., boundary_point= boundary_point, lengh_step = 100)
#     discretization_data =  problem_data.Discretization_Data(mesh, time_discretization, fracture_discretization, tip_enrichement=False)
#     solution = solver.smart_time_solver(discretization_data, physical_data, log_queue)
#     solution.plot()
    # point.plot((fig, ax), marker='+', color='red')
#     fracture_discretization = problem_data.Fracture_Discretization(angular_step = np.pi/20., boundary_point= [point.x, point.y], lengh_step = 100 )
#     discretization_data =  problem_data.Discretization_Data(mesh, time_discretization, fracture_discretization, tip_enrichement=False)
#     solution = solver.smart_time_solver(discretization_data, physical_data, log_queue)
    # solution.plot_displacement()
    # print(point, segment)

# solution = solver.smart_time_solver(discretization_data, physical_data, log_queue)
# solution.plot_displacement()
# solution.plot()


# solution = solver.solver_with_time_discretization(discretization_data, physical_data, log_queue, number_processes = None)
# solution.plot_displacement()
# solution.plot_energy()
