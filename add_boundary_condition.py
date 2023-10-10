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
import griffith.geometry as geo
from griffith.geometry import Point, dist
from griffith.mesh import Edge, closest_point_on_segment
from griffith.mesh import Fracture, Mesh, Broken_Mesh, NotAdmissibleFracture, ShouldBeAdmissibleFracture, Broken_Mesh_Linear_Tip

from scipy.interpolate import LinearNDInterpolator
import shapely.geometry as sg
from shapely.geometry import Polygon, Point
from Func import Unit_vect
from griffith import fracture_iterators, log
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

# figax = plt.subplots()
# fig, ax = figax
# for i in range(len(data)):
#     P = projection_on_boundary(mesh, Point(data[:,0][i], data[:,1][i]))
#     projected_points.append(P)
#     ax.plot(P.x, P.y, 'x', color = 'blue')
# mesh.boundary_mesh.plot(figax)

#modify the boundary data
# data[:, 0] = np.array([P.x for P in projected_points ])
# data[:, 1] = np.array([P.y for P in projected_points ])

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

# def 

xdata, ydata = data[:, 0], data[:, 1]
z1data = data[:,2]
z2data = data[:,3]
x_reordered, y_reordered, z1_reordered, z2_reordered = rearrange_polygon_data(xdata, ydata, data[:,2], data[:,3])

xdata = np.array(x_reordered)
ydata = np.array(y_reordered)
z1data = np.array(z1_reordered)
z2data = np.array(z2_reordered)
# polygon = sg.Polygon(np.column_stack((xdata, ydata)))

# interp_function_x = LinearNDInterpolator(list(zip(xdata, ydata)), z1data)
# interp_function_y = LinearNDInterpolator(list(zip(xdata, ydata)), z2data)

# def project_on_polygon(point, polygon):
#     if polygon.contains(Point(point)):
#         return point
#     nearest_point = polygon.exterior.interpolate(polygon.exterior.project(Point(point)))
#     return np.array([nearest_point.x, nearest_point.y])

# def f_x(x_eval, y_eval):
#     interpolated_value = interp_function_x(x_eval, y_eval)
#     return interpolated_value
#     # return np.interp(x_eval, xdata, z1data) + np.interp(y_eval, ydata, z1data)

# def f_y(x_eval,y_eval):
#     interpolated_value = interp_function_y(x_eval, y_eval)
#     return interpolated_value

# grid_x, grid_y = np.meshgrid(np.linspace(min(xdata)-1, max(xdata)+1, num=100),
#                            np.linspace(min(ydata)-1, max(ydata)+2, num=100))

# grid_values_x = np.vectorize(f_x)(grid_x, grid_y)
# grid_values_y = np.vectorize(f_x)(grid_x, grid_y)

# contour_plot = ax.contourf(grid_x, grid_y, grid_values_x, cmap = 'viridis')

# figax = plt.subplot()
# fig, ax = figax
# contour_plot = ax.contourf(grid_x, grid_y, grid_values_y, cmap = 'viridis')
# ax.set_xlim(min(xdata)-1, max(xdata)+1)
# ax.set_ylim(min(ydata)-1, max(ydata)+2)
# cbar = plt.colorbar(contour_plot)
# cbar.set_label('Data Value')
# mesh.boundary_mesh.plot(figax)

# boundary_displacement = problem_data.Boundary_Displacement_by_percussion(boundary_data = data)
# physical_data = problem_data.Physical_Data(T, 1., boundary_displacement, initial_fracture=None)
# classical_solution = solver.Classical_Solution(mesh=mesh, physical_data=physical_data)
# print(classical_solution.energy)
# classical_solution.plot_displacement()
# classical_solution.plot_energy()

p_init = None #[1036.5, 94.0]
# time_discretization = None

# # look at fracture_iterators 
fracture_discretization = problem_data.Fracture_Discretization(angular_step = np.pi/10., boundary_point= p_init, lengh_step = 100 )


# discretization_data =  problem_data.Discretization_Data(mesh, time_discretization, fracture_discretization, tip_enrichement=False)
# f = fracture_iterators.Admissible_Boundary_Point(fracture_discretization, mesh.boundary_mesh)

# plot points admissibles on boundary 
# fig, ax = mesh.plot()
# for pointsegment in f:
#     point, segment = pointsegment
#     print(point)
#     point.plot((fig, ax), marker='+', color='red')

# plt.show()

# plot fracture admissible from a point on boundary 
f = fracture_iterators.Admissible_Fractures_From_Fixed_Boundary_Point(fracture_discretization, mesh, geo.Point(1015.0, 99.69999694824219))

def tuple_to_string(tple):
  result = "{}".format(tple[0])
  for i in tple[1:]:
    result += "-{}".format(i)
  return result

# # logger = log.Log(args.output_directory + '/' + 'admissible_fracture.log', level=log.INFO)
# # log_queue = logger._log_queue

def compute_fracture(ind, fracture):
  strind = tuple_to_string(ind)
  #import pudb
  #pudb.set_trace()
  try:
    frac_mesh = Broken_Mesh_Linear_Tip(fracture, mesh)
  except NotAdmissibleFracture as e:
    # log_queue.put(('DEBUG', '{} : {} #NOT ADMISSIBLE\n'.format(strind, fracture)))
    print(0)
  except ShouldBeAdmissibleFracture as e:
    print(0)
    # log_queue.put(('WARNING', '{} : {} #SHOULD BE ADMISSIBLE\n'.format(strind, fracture)))
  #except Exception as e:
    #log_queue.put(('ERROR', '{} : {} #ERROR\n'.format(strind, fracture)))
  except NotImplementedError:
      print(0)
    # log_queue.put(('WARNING', '{} : {} #NOT IMPLEMENTED\n'.format(strind, fracture)))
  else:
      print(0)
    # log_queue.put(('INFO', '{} : {}\n'.format(strind, fracture)))

do_fracture = compute_fracture

fig, ax = mesh.plot()
for i, fracture in enumerate(f):
    do_fracture([i], fracture)
    print(fracture)
    fracture.plot((fig,ax))
    p1, p2 = fracture.start_point, fracture.end_point
    p1.plot((fig, ax), marker='+', color='red')
    p2.plot((fig, ax), marker='+', color='red')

    
# logger = log.Log('griffith_solver.log', level=log.INFO, console_output=True)
# logger.log_description(mesh_file=mesh,args=None)
# log_queue = logger._log_queue

# solution = solver.smart_time_solver(discretization_data, physical_data, log_queue)
