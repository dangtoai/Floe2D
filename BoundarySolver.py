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
from scipy.interpolate import LinearNDInterpolator
from Func import Unit_vect
from numpy.linalg import norm

def angle_with_centroid(point, centroid):
    x, y = point
    cx, cy = centroid
    return np.arctan2(y - cy, x - cx)

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

def linear_interpolation(x, y, x1, y1, data1, x2, y2, data2):
    # Compute the interpolated data
    interpolated_data = data1 + ((data2 - data1) * ((x - x1) + (y - y1)) / ((x2 - x1) + (y2 - y1)))
    return interpolated_data

mesh = Mesh("file.msh")

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

xdata, ydata = data[:, 0], data[:, 1]
z1data = data[:,2]
z2data = data[:,3]
x_reordered, y_reordered, z1_reordered, z2_reordered = rearrange_polygon_data(xdata, ydata, data[:,2], data[:,3])

# plt.figure()
# plt.plot(x_reordered, y_reordered, 'x')

n = len(xdata)
Points_inter = []
z1data_new = []
figax = plt.subplots()
fig, ax = figax
for i in range(n):
    i = i%n
    point1 = np.array([x_reordered[i-1], y_reordered[i-1]])
    point2 = np.array([x_reordered[i], y_reordered[i]])
    point3 = np.array([x_reordered[(i+1)%n], y_reordered[(i+1)%n]])
    point4 = np.array([x_reordered[(i+2)%n], y_reordered[(i+2)%n]])
    A1, B1 = line_coefficient(point1, point2)
    A2, B2 = line_coefficient(point3, point4)
    point_inter = intersection_line(A1, B1, A2, B2)
    
    coeff = norm(point_inter - point2)/norm(point1-point2)
    print(coeff)
    new_data = z1data[i-1] + coeff * (z1data[i]-z1data[i-1])
    z1data_new.append(new_data)
    Points_inter.append(point_inter)
    # print(point_inter)
    # t1 = np.linspace(x_reordered[0]-11, x_reordered[0]+11,20)
    # t2 = np.linspace(x_reordered[1]-11, x_reordered[1]+11,20)
    # plt.plot(t1, test_line(t1, A1, B1),'--' )
    # plt.plot(t2, test_line(t2, A2, B2),'--')
    ax.plot(point_inter[0], point_inter[1], 'o')
Points_inter = np.array(Points_inter)
z1data_new = np.array(z1data_new)

xdata = Points_inter[:,0]
ydata = Points_inter[:,1]
interp_function_x = LinearNDInterpolator(list(zip(xdata, ydata)), z1data_new)

# interp_function_x = LinearNDInterpolator(list(zip(xdata, ydata)), z1data)
def f_x(x_eval, y_eval):
    interpolated_value = interp_function_x(x_eval, y_eval)
    return interpolated_value

grid_x, grid_y = np.meshgrid(np.linspace(min(xdata)-1, max(xdata)+1, num=200),
                           np.linspace(min(ydata)-1, max(ydata)+2, num=200))

grid_values_x = np.vectorize(f_x)(grid_x, grid_y)

# figax = plt.subplots()
# fig, ax = figax
ax.plot(x_reordered, y_reordered)
mesh.boundary_mesh.plot(figax)
contour_plot = ax.contourf(grid_x, grid_y, grid_values_x, cmap = 'viridis')
# ax.set_xlim(min(xdata)-1, max(xdata)+1)
# ax.set_ylim(min(ydata)-1, max(ydata)+2)
cbar = plt.colorbar(contour_plot)
# cbar.set_label('Data Value')
# mesh.boundary_mesh.plot(figax)