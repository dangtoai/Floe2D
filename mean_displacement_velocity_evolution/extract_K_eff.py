#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 00:43:49 2025

@author: phandangtoai
"""

import numpy as np
from numpy import array
import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/phandangtoai/Documents/Floe2D/")
from mean_displacement_computation import max_displacement
from Func import Node, Floe, Angle_mat, Traction_mat, Torsion_mat
import os
# from math import sin


radius = 100. 
rho = 900.
Volume = np.pi * radius**2
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
    "data circle network v0 = (-0.5,0) m=1e5",
    "data circle network v0 = (-1,0) m=1e5",
    "data circle network v0 = (-1.5,0) m=1e5",
    "data circle network v0 = (-2,0) m=1e5",
    "data circle network v0 = (-2.5,0) m=1e5",
]


# Filename to read in each directory
filename = "data_spring_1.txt"

all_data = []
Floes = []


i=0


List_K_eff = []
List_M_eff = []
List_K_eff_2 = []
List_M_eff_2 = []

for directory in directories:

    file_path = os.path.join(directory, filename)
    if not os.path.exists(file_path):
        print(f"Warning: File not found in {directory}")
        continue
    
    K_effs = []
    M_contact = []
    K_effs_ = []
    M_contact_ = []
    
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
                    f = Floe(Nodes)
                    Floes.append(f)
                    Angle_m = Angle_mat(f)
                    Traction_m = Traction_mat(f, Angle_m)
                    Torsion_m  = Torsion_mat(f)
                    M_eff, K_eff = f.effective(Traction_m, Torsion_m, 0)
                    M_eff2, K_eff2 = f.effective(Traction_m, Torsion_m, 1)
                    
                    # Append both first-order and second-order effective parameters
                    M_contact.append(M_eff)
                    K_effs.append(K_eff)
                    
                    # Store second-order effective constants
                    M_contact_.append(M_eff2)
                    K_effs_.append(K_eff2)
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
    List_K_eff_2.append(np.mean(K_effs_))
    List_M_eff_2.append(np.mean(M_contact_))

v0 = np.array([0.5, 1, 1.5, 2, 2.5])

# Compute bounds
lower_bounds = [np.sqrt(m / k) for m, k in zip(List_M_eff, List_K_eff)]*v0
upper_bounds = [np.sqrt(m / k) for m, k in zip(List_M_eff_2, List_K_eff_2)]*v0
Numerical_values = -array(max_displacement)
# # Plot both

plt.figure(figsize=(8, 5))
plt.plot(v0, lower_bounds, 'o-', label='Lower Bound (order 1)')
plt.plot(v0, upper_bounds, 'x--', label='Upper Bound (order 2)')
plt.plot(v0, Numerical_values[:,0], label = "Numerical val")
plt.xlabel('Simulation Index')
plt.ylabel('Estimated Max Displacement (m)')
plt.title('Displacement Bounds from Effective Parameters')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()









