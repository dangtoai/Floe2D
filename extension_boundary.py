#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 11:41:27 2025

@author: phandangtoai
"""

import csv
import sys
import numpy as np
from numpy import pi, sin, cos
import matplotlib.pyplot as plt
sys.path.append("/Users/phandangtoai/Documents/Floe2D/Griffith-master_Dimitri")
# import griffith.geometry as geo
from griffith.mesh import Mesh
# from griffith.geometry import Point, dist
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
# from griffith import solver, log, problem_data, fracture_iterators
# from numpy.linalg import norm
from shapely.geometry import Polygon, LineString
from scipy.interpolate import Rbf

radius = 1.

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


mesh = Mesh("unitcircle.msh")
poly_coord = [p.array for p in mesh.boundary_mesh.points]
polygon = Polygon(poly_coord)


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

# plt.figure()
# plt.plot(xdata, ydata, 'o')
# for i, (x, y) in enumerate(zip(xdata, ydata)):
#     plt.text(x, y, str(i), fontsize=12, ha='right', va='bottom', color='red')

tri = Delaunay(np.column_stack((xdata, ydata)))

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")
ax.plot_trisurf(xdata, ydata, z1data, triangles=tri.simplices, cmap="coolwarm", edgecolor="black", alpha=0.8)
plt.title("$u_1$")

# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111, projection="3d")
# ax.plot_trisurf(xdata, ydata, z2data, triangles=tri.simplices, cmap="coolwarm", edgecolor="black", alpha=0.8)
# plt.title("$u_2$")


interp_function_x = LinearNDInterpolator(tri, z1data)
interp_function_y = LinearNDInterpolator(tri, z2data)


def f_x(x_eval, y_eval):
    interpolated_value = interp_function_x(x_eval, y_eval)
    Cones = cones_list(tri) 
    if np.isnan(interpolated_value):
        # print("problem at point ", x_eval, y_eval)
        interpolated_value = 0.
        l = 0                       # triangle count
        L = []                      # stock the list of triangles 
        for i,j,k in Cones:
            if inside_cone(Points[i], Points[j], Points[k], np.array([x_eval, y_eval])):
                l += 1
                L.append([i,j,k])
        
        
        if l == 1:
            i0, j0, k0 = L[0]
            P_ = np.array([Points[i0], Points[j0], Points[k0]])
            data_ = np.array([z1data[i0], z1data[j0], z1data[k0]])
            A, B, C = P1_coefficient(P_, data_)
            interpolated_value += A*x_eval + B*y_eval + C
            # print('value = ', interpolated_value)
            return interpolated_value
        
        # print(l, L)
        if l >= 2: 
            i = 0
            while i < len(L):
                j = i+1
                while j < len(L):
                    intersection = np.intersect1d(L[i], L[j])
                    if len(intersection) == 0:
                        del L[j]
                    else: j += 1
                i += 1
                
        #     # create new data from the interpolation with the mesh's boundary
        #     # find the intersection between the cone and the mesh's boundary 
            if len(L) >= 2:
                Origine1 = Points[L[0][0]]
                # print(Points[L[0][0]])
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

def Dirichlet(x, y):
    # rbfx = Rbf(xdata, ydata, z1data, function='multiquadric')
    # rbfy = Rbf(xdata, ydata, z2data, function='multiquadric')
    rbfx = Rbf(xdata, ydata, z1data, function='linear')
    rbfy = Rbf(xdata, ydata, z2data, function='linear')

    # return np.array([f_x(x, y), f_y(x, y)])
    return np.array([rbfx(x,y), rbfy(x,y)])


print(" test value at (1., 0.) = ", Dirichlet(1., 0.))


# taking 41 points on the right side of the circle and compute the Dirichlet boundary condition

approx_data = np.array([Dirichlet(radius * cos(theta), radius * sin(theta))
               for theta in np.linspace(-pi/2, pi/2, 41) ])
# print(approx_data)

new_xdata = [radius * cos(theta) for theta in np.linspace(-pi/2, pi/2, 41) ]
new_ydata = [radius * sin(theta) for theta in np.linspace(-pi/2, pi/2, 41) ]
new_z1data = approx_data[:,0] 
new_z2data = approx_data[:,1] 

tri = Delaunay(np.column_stack((new_xdata, new_ydata)))

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")
ax.plot_trisurf(new_xdata, new_ydata, new_z1data, triangles=tri.simplices, cmap="coolwarm", edgecolor="black", alpha=0.8)
plt.title("$u_1$")


plt.figure()
plt.quiver(new_xdata, new_ydata, new_z1data, new_z2data, angles='xy', scale_units='xy', scale=1)




# with open("data_circle_RBF.txt", "a") as file:  # 'a' mode appends instead of overwriting
#     for row in approx_data:
        
#         file.write(f"{row[0]} {row[1]}\n")  # Write each element in row format
        
        
# with open("data_circle_RBF.txt", "a") as file:  # 'a' mode appends instead of overwriting
#     for row in approx_data:
        
#         file.write(f"{row[0]} {row[1]}\n")  # Write each element in row format


# grid_x, grid_y = np.meshgrid(np.linspace(-1., 1., num = 50),
#                             np.linspace(-1., 1., num = 50))

# grid_values_x = np.vectorize(f_x)(grid_x, grid_y)
# grid_values_y = np.vectorize(f_y)(grid_x, grid_y)

# figax = plt.subplots()
# fig, ax = figax
# for i in range(data.shape[0]):
#     ax.plot(Points[i][0], Points[i][1], 'x')

# contour_plot = ax.contourf(grid_x, grid_y, grid_values_x, cmap = 'viridis')
# ax.set_xlim(-1.1, 1.1)
# ax.set_ylim(-1.1, 1.1)
# cbar = plt.colorbar(contour_plot)
# cbar.set_label('Data Value')

