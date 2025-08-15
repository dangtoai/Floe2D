#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 17:26:48 2025

@author: phandangtoai
"""

import numpy as np
from numpy import array
import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/phandangtoai/Documents/Floe2D/")
from mean_displacement_computation import max_displacement
from Func import Node, Floe, Angle_mat, Traction_mat
import os
# from math import sin


radius = 100. 
rho = 900.
Volume = 100**2
# Volume = np.pi * radius**2
Mass = Volume * rho

matrix = np.array([[32*Volume/(9*np.pi**2), 3*Volume/4],
                    [32*Volume/(9*np.pi**2), -3*Volume/4]])

### Young modulus E and Poisson ration nu

E, nu = 8.95*10**9, 0.295

### Lame parameters
lamb, mu = E*nu/ ((1+nu)*(1-2*nu)), E/(2*(1+nu))
    
### K,G springs stiffness
traction_stiff, torsion_stiff = np.linalg.solve(matrix, np.array([lamb, mu])) 


directories = [
    "data circle network v0 = (-1,0) m=1e2",
    # "data circle network v0 = (-1,0) m=150",
    "data circle network v0 = (-1,0) m=1e3",
    "data circle network v0 = (-1,0) m=1e4",
    "data circle network v0 = (-1,0) m=1e5",
    "data circle network v0 = (-1,0) m=1e6",
    "data circle network v0 = (-1,0) m=1e7",
    # "data circle network v0 = (-1,0) m=1e8",
]

# Filename to read in each directory
filename = "data_spring_1.txt"

all_data = []
Floes = []


i=0


List_K_eff = []
List_M_eff = []

for directory in directories:

    file_path = os.path.join(directory, filename)
    if not os.path.exists(file_path):
        print(f"Warning: File not found in {directory}")
        continue
    
    K_effs = []
    M_contact = []
    
    with open(file_path, "r") as file:
        current_data = []
        inside_floe = False

        for line in file:
            line = line.strip()

            if line.startswith('"T*'):
                # Start of new floe block
                current_data = []
                inside_floe = True
                continue

            if all(part == '-' for part in line.split()) or '---' in line:
                # End of floe block
                if inside_floe and current_data:
                    print(i)
                    data_array = np.array(current_data, dtype=float)
                    positions = data_array[:, :2]  # Keep only x, y
                    Nodes = [Node([x, y]) for x, y in positions]
                    Floes.append(Floe(Nodes))
                    
                    Angle_m = Angle_mat(Floes[-1])
                    Neighbors = Floes[-1].Neighbors()[Floes[-1].n-1]
                    M_contact.append([Mass/ Floes[-1].n])
                    Traction_Mat = Traction_mat(Floes[-1], Angle_m)
                    K_neighbors = [Traction_Mat[Floes[-1].n-1, i] for i in Neighbors]
                    K_eff = sum(K_neighbors)
                    K_effs.append(K_eff)
                    
                    i=i+1
                    
                inside_floe = False
                continue

            # Skip empty lines
            if not line:
                continue

            try:
                values = list(map(float, line.split()))
                current_data.append(values)
            except ValueError:
                print(f"Skipping invalid line: {line}")
                continue
    List_K_eff.append(np.mean(K_effs))
    List_M_eff.append(np.mean(M_contact))

M_P = np.array([1e2, 1e3, 1e4, 1e5, 1e6, 1e7])
M_eff_array = np.array(List_M_eff)
K_eff_array = np.array(List_K_eff)

Effective_sol = np.sqrt((M_eff_array+M_P)/K_eff_array)

Numerical_values = array(max_displacement)

# plt.figure()
# plt.loglog(M_P, Effective_sol, label = "estimation")
# plt.loglog(M_P, -Numerical_values[:,0], label = "numerical val")
# plt.legend()

plt.figure()
plt.semilogx(M_P, Effective_sol, label = "Estimation")
plt.semilogx(M_P, -Numerical_values[:,0], label = "Numerical val")
plt.legend()







