#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Visualize the data from the outputfiles for Schroedinger equation."""

import os.path
import numpy as np
import matplotlib.pyplot as plt


def _visualize(direc, bulge_factor, lim_x, lim_y):
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
    energy = np.loadtxt(os.path.join(direc, "energies.dat"))
    wavefuncs = np.loadtxt(os.path.join(direc, "wavefuncs.dat"))[:, 1:]
    expvalues = np.loadtxt(os.path.join(direc, "expvalues.dat"))[:, 0]
    uncertainty = np.loadtxt(os.path.join(direc, "expvalues.dat"))[:, 1]

    # plotting potential, energies, wafefunctions and expactation values
    plt.subplot(1, 2, 1)
    plt.plot(xx, potential, color="black")

    for ii in range(0, int(len(energy))):
        energies = np.linspace(energy[ii], energy[ii], len(xx))
        plt.plot(xx, energies, color="gray")
        plt.plot(expvalues[ii], energy[ii], "xg")
        if ii % 2 == 0:
            plt.plot(xx, bulge_factor*wavefuncs[:, ii]+energy[ii],
                     color="blue")
        else:
            plt.plot(xx, bulge_factor*wavefuncs[:, ii]+energy[ii], color="red")
    plt.xlim(xx.min() - lim_x, xx.max() + lim_x)
    plt.ylim(potential.min() - lim_y, energy.max() + lim_y)
    plt.xlabel("$x$ [Bohr]")
    plt.ylabel("Energy [Hartree]")
    plt.title("Potential, eigenstates, ($x$)")

    # plotting energies and uncertainty
    plt.subplot(1, 2, 2)
    for ii in range(0, int(len(energy))):
        energies = np.linspace(energy[ii], energy[ii], len(xx))
        plt.plot(xx, energies, color="gray")
        plt.plot(uncertainty[ii], energy[ii], "+m")
        plt.yticks([])
    plt.xlim(0, uncertainty.max() + lim_x)
    plt.ylim(potential.min() - lim_y, energy.max() + lim_y)
    plt.xlabel("[Bohr]")
    plt.title("$\sigma_\mathrm{x}$")

    # save figure as solution_schroedinger.pdf
    fname = os.path.join(direc, 'solution_schoedinger.pdf')
    plt.savefig(fname, format="pdf")


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
    # private visualize function will be executed
    _visualize(direc, bulge_factor, lim_x, lim_y)
