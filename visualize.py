#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Executable file which contains solver and visualizer"""

import argparse
from schroed_visualize import visualize


parser = argparse.ArgumentParser(description='Visualizes solution of the\
                                 solver')

msg_direc = 'Set directory of the figures\
    (default value: ./inputdata/potential_well)'
parser.add_argument('-d', '--directory', default='./inputdata/potential_well',
                    help=msg_direc)

msg_bulge_factor = 'Set bulge factor for the wavefuncs\
                    (default value: 0.3)'
parser.add_argument('-bf', '--bulge_factor', default=0.3, type=float,
                    help=msg_bulge_factor)

msg_lim_x = 'Set x limits of the figures\
                    (default value: 0.2)'
parser.add_argument('-xl', '--lim_x', default=0.2, type=float,
                    help=msg_lim_x)

msg_lim_y = 'Set y limits of the figures\
                    (default value: 0.5)'
parser.add_argument('-yl', '--lim_y', default=0.5, type=float,
                    help=msg_lim_y)

args = parser.parse_args()
print("Directory: {}".format(args.directory))
print("bulge factor: {:f}".format(args.bulge_factor))
print("x limit: {:f}".format(args.lim_x))
print("y limit: {:f}".format(args.lim_y))

visualize(args.directory, args.bulge_factor, args.lim_x, args.lim_y)
