#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 16:54:02 2025

@author: phandangtoai
"""

import subprocess
# Runs file.py 100 times

num_procs = 6
scripts = ["Circle_geometry_simulation.py", 
           "deformation_computation_circle.py", 
           "extension_boundary.py"]

###
for _ in range(10):
    print("simulation number: ",_+1)
    for script in scripts:
        print(f"Running {script} with {num_procs} processes...")
        subprocess.run(["mpirun", "-np", str(num_procs), "python", script])
        print(f"Finished {script}\n")



print("All scripts completed successfully!")
# Open the file in read mode



