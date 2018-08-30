#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Executable file which contains solver and visualizer"""

import argparse
from schroed_solver import solver

parser = argparse.ArgumentParser(description='Executes solver for the\
                                              Schroedinger equation')

msg_direc = 'Set directory of the figures\
    (default value: .)'
parser.add_argument('-d', '--directory', default='.',
                    help=msg_direc)

args = parser.parse_args()
print("Directory: {}".format(args.directory))

solver(args.directory)
