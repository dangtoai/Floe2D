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
from griffith.geometry import Point, dist
from griffith.mesh import Edge, closest_point_on_segment
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon, Point
from Func import Unit_vect
# from numpy.linalg import norm

mesh = Mesh("file.msh")
T = problem_data.lame_tensor_ice #Lame tensor

data = []
with open('boundary_data.csv', 'r', encoding='utf-8') as csv_file:
    lines = csv_file.readlines()
    lines_to_read = lines[1:-2]
    csv_reader = csv.reader(csv_file)
    csv_reader = csv.reader(lines_to_read)
    for row in csv_reader:
        # print(row)
        data.append( [float(item) for item in row[0].split()] )

data = np.array(data)


projected_points = []

figax = plt.subplots()
fig, ax = figax
for i in range(len(data)):
    P = projection_on_boundary(mesh, Point(data[:,0][i], data[:,1][i]))
    projected_points.append(P)
    ax.plot(P.x, P.y, 'x', color = 'blue')
mesh.boundary_mesh.plot(figax)

#modify the boundary data
data[:, 0] = np.array([P.x for P in projected_points ])
data[:, 1] = np.array([P.y for P in projected_points ])

# Function to calculate the angle of a point with respect to the centroid
def angle_with_centroid(point, centroid):
    x, y = point
    cx, cy = centroid
    return np.arctan2(y - cy, x - cx)

# Function to rearrange points to form a proper polygon
def rearrange_polygon_points(xdata, ydata):
    # Calculate the centroid
    centroid = (np.mean(xdata), np.mean(ydata))
    # Create a list of points and calculate angles with respect to the centroid
    points = list(zip(xdata, ydata))
    angles = [angle_with_centroid(point, centroid) for point in points]
    # Sort the points based on the angles
    sorted_points = [point for _, point in sorted(zip(angles, points))]
    return zip(*sorted_points)

def rearrange_polygon_data(xdata, ydata, z1data, z2data):
    x_reordered, y_reordered = rearrange_polygon_points(xdata, ydata)
    # Rearrange z1data based on the new order of x and y
    z1_reordered = np.array([z1data[np.where((xdata == x) & (ydata == y))[0][0]] for x, y in zip(x_reordered, y_reordered)])
    # Rearrange z2data based on the new order of x and y
    z2_reordered = np.array([z2data[np.where((xdata == x) & (ydata == y))[0][0]] for x, y in zip(x_reordered, y_reordered)])

    return x_reordered, y_reordered, z1_reordered, z2_reordered

xdata, ydata = data[:, 0], data[:, 1]

x_reordered, y_reordered, z1_reordered, z2_reordered = rearrange_polygon_data(xdata, ydata, data[:,2], data[:,3])

xdata = np.array(x_reordered)
ydata = np.array(y_reordered)
z1data = np.array(z1_reordered)
z2data = np.array(z2_reordered)

# polygon = Polygon(np.column_stack((xdata, ydata)))

# interp_function_x = LinearNDInterpolator(list(zip(xdata, ydata)), z1data)
# interp_function_y = LinearNDInterpolator(list(zip(xdata, ydata)), z2data)

# def project_on_polygon(point, polygon):
#     if polygon.contains(Point(point)):
#         return point
#     nearest_point = polygon.exterior.interpolate(polygon.exterior.project(Point(point)))
#     return np.array([nearest_point.x, nearest_point.y])

# def f_x(x_eval, y_eval):
#     return np.interp(x_eval, xdata, z1data) + np.interp(y_eval, ydata, z1data)

# def f_y(x_eval,y_eval):
#     return np.interp(x_eval, xdata, z2data) + np.interp(y_eval, ydata, z2data)

# grid_x, grid_y = np.meshgrid(np.linspace(min(xdata)-1, max(xdata)+1, num=100),
#                           np.linspace(min(ydata)-1, max(ydata)+2, num=100))

# grid_values_x = np.vectorize(f_x)(grid_x, grid_y)
# grid_values_y = np.vectorize(f_y)(grid_x, grid_y)
# contour_plot = ax.contourf(grid_x, grid_y, grid_values_x, cmap = 'viridis')

# ax.set_xlim(min(xdata)-1, max(xdata)+1)
# ax.set_ylim(min(ydata)-1, max(ydata)+2)
# cbar = plt.colorbar(contour_plot)
# cbar.set_label('Data Value')
# mesh.boundary_mesh.plot(figax)


# figax = plt.subplots()
# fig, ax = figax
# # ax.plot(xdata, ydata)
# for i in range(len(data)):
#     P = projection_on_boundary(mesh, Point(data[:,0][i], data[:,1][i]))
#     projected_points.append(P)
#     ax.plot(P.x, P.y, 'x', color = 'blue')
# mesh.boundary_mesh.plot(figax)


# contour_plot = ax.contourf(grid_x, grid_y, grid_values_y, cmap = 'viridis')
# ax.set_xlim(min(xdata)-1, max(xdata)+1)
# ax.set_ylim(min(ydata)-1, max(ydata)+2)
# cbar = plt.colorbar(contour_plot)
# cbar.set_label('Data Value')
# mesh.boundary_mesh.plot(figax)


boundary_displacement = problem_data.Boundary_Displacement_by_percussion(boundary_data = data)
physical_data = problem_data.Physical_Data(T, 1., boundary_displacement, initial_fracture=None)
classical_solution = solver.Classical_Solution(mesh=mesh, physical_data=physical_data)
print(classical_solution.energy)
classical_solution.plot_displacement()
classical_solution.plot_energy()
