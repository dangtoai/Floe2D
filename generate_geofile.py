#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 17:05:00 2023

@author: phandangtoai
"""
import sys
import scipy as sp
import csv
import numpy as np
from subprocess import run
sys.path.append(
    "/Users/phandangtoai/Documents/Floe2D/Griffith-master_Dimitri/")
from griffith.geometry import Point, dist
from griffith.mesh import Mesh

Input = input("Select a geometry to run the random process: ")
Input = int(Input)
print(f"the geometry number {Input}")
# 100 geometries of ice floes
mat = sp.io.loadmat("Biblio_Floes.mat")['G'][Input][0]
lc = 30.
GEO_SCRIPT = "/* the geometry file\n */ \n\n"
GEO_SCRIPT += f"lc = {lc}; // mesh precision \n"
Points = []
for i, point in enumerate(mat, start=1):
    GEO_SCRIPT += f"Point({i}) = {{{point[0]}, {point[1]}, 0 , lc}};\n"
    Points.append(Point(point[0], point[1]))

for i in range(1, len(mat)-1):
    GEO_SCRIPT += f"Line({i}) = {{{i}, {i+1}}};\n"

GEO_SCRIPT += f"Line({i+1}) = {{{i+1}, {1}}};\n"

# line_loop_command = f"Line Loop(1) = {{{NUMBERS_STRING}}};"

# define Dirichlet and Neumann part:
with open('boundary_data.csv', mode='r') as csv_file:
    csvreader = csv.reader(csv_file)
    # print(csvreader)
    for row in csvreader:
        print(row)

# numbers = [float(num) for num in row[0].split()]
# contact_point = np.array(numbers)
contact_point = np.array([1006.4074181139762, 119.04863615511762])

# define 2 Dirichlet's region: 1 near the collision area and 1 far from it
THRESHOLD = 4.  # TODO: Support of Dirichlet boundary depends on floe size

# 1st Dirichlet's part
Condition1 = np.array(
    [dist(Point(contact_point[0], contact_point[1]), p) for p in Points]) < THRESHOLD
index_dirichlet_points1 = np.where(Condition1)[0]
print(index_dirichlet_points1)

# 2nd Dirichlet's part
Condition2 = np.array([dist(Point(contact_point[0], contact_point[1]), p)
                      for p in Points]) > 7*THRESHOLD
index_dirichlet_points2 = np.where(Condition2)[0]
print(index_dirichlet_points2)

combined_condition = Condition1 | Condition2
index_dirichlet_lines = np.where(combined_condition)[0]
print(index_dirichlet_lines)

# index_dirichlet_lines = range( min(index_dirichlet_points1) - 1, max(index_dirichlet_points1) + 1 )
index_neumann_lines = [i for i in range(
    1, len(mat)) if i not in index_dirichlet_lines]

NUMBERS_STRING = ', '.join(str(i) for i in index_neumann_lines)
GEO_SCRIPT += f'Physical Line("N") = {{{NUMBERS_STRING}}};\n'

NUMBERS_STRING = ', '.join(str(i) for i in index_dirichlet_lines)
GEO_SCRIPT += f'Physical Line("D") = {{{NUMBERS_STRING}}};\n'

NUMBERS_STRING = ', '.join(str(i) for i in range(1, len(mat)+1))
GEO_SCRIPT += f"Line Loop({1}) = {{{NUMBERS_STRING}}};\n"

GEO_SCRIPT += f"Plane Surface({1}) = {{1}};\n"

GEO_SCRIPT += 'Physical Surface("S") = {1};\n'

with open('file.geo', 'w', encoding='utf-8') as geo_file:
    geo_file.write(GEO_SCRIPT)

run('gmsh -2 file.geo -o file.msh -format msh2', shell=True)
run('gmsh file.msh', shell = True)

mesh = Mesh("file.msh")
