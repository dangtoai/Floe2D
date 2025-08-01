#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  6 01:27:05 2025

@author: phandangtoai
"""

import subprocess
from mpi4py.futures import MPIPoolExecutor  # Import MPIPoolExecutor for dynamic task distribution

# Runs file.py 100 times

num_procs = 1
scripts = ["generate_network.py", 
           "deformation_computation.py", 
           "extension_geometry3.py"]



####
for _ in range(25):
    print("simulation number: ",_+1)
    for script in scripts:
        # print(f"Running {script} with {num_procs} processes...")
        # subprocess.run(["mpirun", "-np", str(num_procs), "python", script])
        # print(f"Finished {script}\n")
        subprocess.run(["python3", script])




print("All scripts completed successfully!")
# Open the file in read mode



