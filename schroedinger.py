#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 16:39:23 2018

@author: isi
"""

import schroed_solver
import schroed_visualize

direc = "inputdata/asym_potential_well"
bulge_factor = 0.2
lim_x = 0.5
lim_y = 0.2

schroed_solver.solver(direc)
schroed_visualize.visualize(direc, bulge_factor, lim_x, lim_y)
