#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 14:20:57 2018

@author: isi
"""

import os.path
import pytest
import numpy as np
from schroed_solver import solver


# creat parametrized tests for potenial and energies of every input example
@pytest.mark.parametrize("direc", [
    ("/Users/isi/Desktop/schroedinger_project/inputdata/asym_potential_well"),
    ("/Users/isi/Desktop/schroedinger_project/inputdata/potential_well"),
    ("/Users/isi/Desktop/schroedinger_project/inputdata/infinit_potential_well"),
    ("/Users/isi/Desktop/schroedinger_project/inputdata/harm_osz"),
    ("/Users/isi/Desktop/schroedinger_project/inputdata/double_potential_well"),
    ("/Users/isi/Desktop/schroedinger_project/inputdata/double_potential_wellc")])
def test_function(direc):
    """test interpolated discretised potential and energies"""
    # import values for potential and energy from the solver
    xypotential, energy = solver(direc)

    # open reference data for the potential
    filepot = open(os.path.join(direc, 'ref_potential.dat'), "r")
    linespot = filepot.readlines()
    filepot.close()

    # create array for reference potential to compare with calculated potential
    refpotential = np.zeros((int(len(xypotential)), 2))
    for ii in range(0, len(xypotential)):
        refpotential[ii, 0] = float(linespot[ii].split()[0])
        refpotential[ii, 1] = float(linespot[ii].split()[1])

    # checks whether reference and calculated values are the same
    assert np.all(np.abs(refpotential - xypotential) < 1e-10)

    # open reference data for the energy
    fileener = open(os.path.join(direc, 'ref_energies.dat'), "r")
    linesener = fileener.readlines()
    fileener.close()

    # create array for reference energies to compare with calculated energies
    refenergies = np.zeros((int(len(energy)), ))
    for ii in range(0, len(energy)):
        refenergies[ii] = float(linesener[ii].split()[0])

    # checks whether reference and calculated values are the same
    assert np.all(np.abs(refenergies - energy) < 1e-10)
