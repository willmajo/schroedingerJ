#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Executable file which contains solver and visualizer"""

import argparse
from schroed_solver import solver

PARSER = argparse.ArgumentParser(description='Executes solver for the\
                                              Schroedinger equation')

MSG = 'Set directory of the figures\
    (default value: ./inputdata/potential_well)'
PARSER.add_argument('-d', '--directory', default='./inputdata/potential_well',
                    help=MSG)

ARGS = PARSER.parse_args()
print("Directory: {}".format(ARGS.directory))

solver(ARGS.directory)
