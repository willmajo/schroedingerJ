#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Executable file which contains solver and visualizer"""

import schroed_solver
import schroed_visualize

# enter the directory of the problem
direc = "inputdata/asym_potential_well"

# set options for a better illustration of the solutions
bulge_factor = 0.2
lim_x = 0.5
lim_y = 0.2

# execute solver and visualizer
schroed_solver.solver(direc)
schroed_visualize.visualize(direc, bulge_factor, lim_x, lim_y)
