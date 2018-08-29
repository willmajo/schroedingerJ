#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Visualize the data from the outputfiles for Schroedinger equation."""

import os.path
import numpy as np
import matplotlib.pyplot as plt


def visualize(direc, bulge_factor, lim_x, lim_y):
    """Read data from outputfiles and create plots for potential, eigenvalues,
    wavefunctions, expactation values and uncertainty for the position.

    Args:
        direc: directory of the outputfiles
        bulge_factor: bulge amplitude of the wavefunctions
        lim_x: set the limits for the x-axis with additional summand
        lim_y: set the limits for the y-axis with additional summand

    Returns:
    """
    # load data files
    xx = np.loadtxt(os.path.join(direc, "potential.dat"))[:, 0]
    potential = np.loadtxt(os.path.join(direc, "potential.dat"))[:, 1]
    energies = np.loadtxt(os.path.join(direc, "energies.dat"))
    wavefuncs = np.loadtxt(os.path.join(direc, "wavefuncs.dat"))[:, 1:]
    expvalues = np.loadtxt(os.path.join(direc, "expvalues.dat"))[:, 0]
    uncertainty = np.loadtxt(os.path.join(direc, "expvalues.dat"))[:, 1]

    # plotting potential, energies, wafefunctions and expactationvalues
    plt.subplot(1, 2, 1)
    plt.plot(xx, potential, color="black")

    bf = bulge_factor
    limx = lim_x
    limy = lim_y
    bb = len(energies)
    for ii in range(0, bb):
        energy = np.linspace(energies[ii], energies[ii], len(xx))
        plt.plot(xx, energy, color="gray")
        plt.plot(expvalues, energies, "xg")
        if ii % 2 == 0:
            plt.plot(xx, bf*wavefuncs[:, ii]+energies[ii], color="blue")
        else:
            plt.plot(xx, bf*wavefuncs[:, ii]+energies[ii], color="red")
    plt.xlim(xx.min() - limx, xx.max() + limx)
    plt.ylim(potential.min() - limy, energies.max() + limy)
    plt.xlabel("$x$ [Bohr]")
    plt.ylabel("Energy [Hartree]")
    plt.title("Potential, eigenstates, ($x$)")

    # plotting energies and uncertainty
    plt.subplot(1, 2, 2)
    bb = len(energies)
    for ii in range(0, bb):
        energy = np.linspace(energies[ii], energies[ii], len(xx))
        plt.plot(xx, energy, color="gray")
        plt.plot(uncertainty, energies, "+m")
        plt.yticks([])
    plt.xlim(0, uncertainty.max() + limx)
    plt.ylim(potential.min() - limy, energies.max() + limy)
    plt.xlabel("[Bohr]")
    plt.title("$\sigma_\mathrm{x}$")

    fname = os.path.join(direc, 'solution_schoedinger.pdf')
    plt.savefig(fname, format="pdf")
