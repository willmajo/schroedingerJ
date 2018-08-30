#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Test interpolated potential and energies with reference data"""

import os.path
import pytest
import numpy as np
from schroed_solver import solver


# creat parametrized tests for potential of every input example
@pytest.mark.parametrize("direc", [
    ("./inputdata/asym_potential_well"),
    ("./inputdata/potential_well"),
    ("./inputdata/infinit_potential_well"),
    ("./inputdata/harm_osz"),
    ("./inputdata/double_potential_well"),
    ("./inputdata/double_potential_well_cubic")])
def test_potential(direc):
    """Read data for potential from reference files and outputfiles, and
    compare them with each other.

    Args:
        direc: directory of the outputfiles and reference files

    Returns:

    """
    # import values for potential from the solver
    xypotential = solver(direc)[0]

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


# creat parametrized tests for energy of every input example
@pytest.mark.parametrize("direc", [
    ("./inputdata/asym_potential_well"),
    ("./inputdata/potential_well"),
    ("./inputdata/infinit_potential_well"),
    ("./inputdata/harm_osz"),
    ("./inputdata/double_potential_well"),
    ("./inputdata/double_potential_well_cubic")])
def test_energy(direc):
    """Read data for energy from reference files and outputfiles, and
    compare them with each other.

    Args:
        direc: directory of the outputfiles and reference files

    Returns:

    """

    # import values for energy from the solver
    energy = solver(direc)[1]

    # open reference data for the energy
    fileener = open(os.path.join(direc, 'ref_energies.dat'), "r")
    linesener = fileener.readlines()
    fileener.close()

    # create array for reference energies to compare with calculated energies
    refenergies = np.zeros((int(len(energy)), ))
    for ii in range(0, len(energy)):
        refenergies[ii] = float(linesener[ii].split()[0])

    # checks whether reference and calculated values are the same
    assert np.all(np.abs(refenergies - energy) < 1e-2)
