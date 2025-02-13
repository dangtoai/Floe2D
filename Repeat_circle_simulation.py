#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 16:54:02 2025

@author: phandangtoai
"""

import subprocess
# import numpy as np
# import csv
# import matplotlib.pyplot as plt 

# Runs file.py 100 times

for _ in range(100):
    subprocess.run(["python", "Circle_geometry_simulation.py"])
    subprocess.run(["python", "deformation_computation_circle.py"])
    subprocess.run(["python", "extension_boundary.py"])  # Runs file.py 100 times

# Open the file in read mode



