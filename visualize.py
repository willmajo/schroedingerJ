#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Executable file which contains solver and visualizer"""

import argparse
from schroed_visualize import visualize


PARSER = argparse.ArgumentParser(description='Visualizes solution of the\
                                 solver')

MSG_DIREC = 'Set directory of the figures\
    (default value: ./inputdata/potential_well)'
PARSER.add_argument('-d', '--directory', default='./inputdata/potential_well',
                    help=MSG_DIREC)

MSG_BULGE_FACTOR = 'Set bulge factor for the wavefuncs\
                    (default value: 0.3)'
PARSER.add_argument('-bf', '--bulge_factor', default=0.3, type=float,
                    help=MSG_BULGE_FACTOR)

MSG_LIM_X = 'Set x limits of the figures\
                    (default value: 0.2)'
PARSER.add_argument('-xl', '--lim_x', default=0.2, type=float,
                    help=MSG_LIM_X)

MSG_LIM_Y = 'Set y limits of the figures\
                    (default value: 0.5)'
PARSER.add_argument('-yl', '--lim_y', default=0.5, type=float,
                    help=MSG_LIM_Y)

ARGS = PARSER.parse_args()
print("Directory: {}".format(ARGS.directory))
print("bulge factor: {:f}".format(ARGS.bulge_factor))
print("x limit: {:f}".format(ARGS.lim_x))
print("y limit: {:f}".format(ARGS.lim_y))

visualize(ARGS.directory, ARGS.bulge_factor, ARGS.lim_x, ARGS.lim_y)
