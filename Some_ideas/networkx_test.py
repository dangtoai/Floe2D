#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 16:30:05 2022

@author: phandangtoai
"""
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

G = nx.Graph()
l = [i for  i in range(9)]

G.add_nodes_from(l)
ubax1 = plt.subplot(121)

G.add_edges_from([(0, 1),
 (0, 3),
 (0, 5),
 (0, 8),
 (1, 2),
 (1, 4),
 (1, 5),
 (1, 8),
 (2, 4),
 (2, 5),
 (2, 6),
 (2, 7),
 (3, 5),
 (3, 7),
 (4, 6),
 (4, 8),
 (5, 7),
 (6, 7)])

nx.draw(G, with_labels=True, font_weight='bold')
for path in nx.all_simple_edge_paths(G, 8, 3, cutoff=(7)):
    if len(path) == 7: print(path)

#  idee pour la suite, chercher tous les chemin simple