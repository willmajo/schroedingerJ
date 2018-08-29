#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 16:52:05 2018

@author: jonaswillmann
"""

import os.path
import sys
import numpy as np
import scipy
from scipy.interpolate import interp1d


def solver(direc):
    '''import parameters and solve the one-dimensional time-independent
    Schrödinger equation'''
    # import input file and test if file exists
    file = 'schroedinger.inp'
    fn = os.path.join(direc, file)
    inp_all = open(fn, 'r')

    # import parameters from input file
    lines = inp_all.readlines()
    inp_mass = lines[0].split()[0]
    mass = float(inp_mass)
    inp_xmin = lines[1].split()[0]
    xmin = float(inp_xmin)
    inp_xmax = lines[1].split()[1]
    xmax = float(inp_xmax)
    inp_npoint = lines[1].split()[2]
    npoint = int(inp_npoint)
    inp_firsteigval = lines[2].split()[0]
    firsteigval = int(inp_firsteigval)
    inp_lasteigval = lines[2].split()[1]
    lasteigval = int(inp_lasteigval)
    interpoltype = lines[3].split()[0]
    inp_nrinterpolpoints = lines[4].split()[0]
    nrinterpolpoints = int(inp_nrinterpolpoints)
    len_pot = len(lines)-5
    xpot = np.zeros(len_pot)
    ypot = np.zeros(len_pot)
    for ii in range(5, len_pot+5):
        xpot[ii - 5] = float(lines[ii].split()[0])
        ypot[ii - 5] = float(lines[ii].split()[1])

    # read interpolation type and use interpolation for potential
    xx = np.linspace(xmin, xmax, npoint)
    if interpoltype == 'linear':
        pot = np.interp(xx, xpot, ypot)
    elif interpoltype == 'polynomial':
        degree = int(nrinterpolpoints - 1)
        coef = np.polyfit(xpot, ypot, degree)
        polf = np.poly1d(coef)
        pot = polf(xx)
    elif interpoltype == 'cspline':
        cubicf = interp1d(xpot, ypot, kind='cubic')
        pot = cubicf(xx)
    else:
        print('interpolation type not found')
        sys.exit(1)

    # write values for xPot and yPot into potential.dat file
    potential = np.array([xx, pot])
    xypotential = potential.T
    np.savetxt(os.path.join(direc, 'potential.dat'), xypotential)

    # formulate matrix-problem for the discretised Schroedinger equation
    matrix = np.zeros((npoint, npoint))
    delta = abs((xmax-xmin)/(npoint))
    aa = 1/(mass*(delta)**2)
    for ii in range(1, npoint):
        matrix[ii, ii-1] = -aa/2
    for ii in range(0, npoint):
        matrix[ii, ii] = aa+pot[ii]
    for ii in range(0, npoint-1):
        matrix[ii, ii+1] = -aa/2

    # compute eigenvalues and eigenvectors
    energy, wavefunc = scipy.linalg.eigh(matrix, eigvals=(int(firsteigval-1),
                                                          int(lasteigval-1)))

    # normalize wavefunctions
    deltavec = delta*np.ones((1, npoint))
    wavefunc_sq = wavefunc**2
    norm_sq = np.dot(deltavec, wavefunc_sq)
    norm = 1/(np.sqrt(norm_sq))
    norm_wavefunc = np.dot(wavefunc, np.diag(np.reshape(norm,
                                                        (len(energy), ))))

    # save eigenvalues and eigenvectors in energies.dat and wavefuncs.dat files
    np.savetxt(os.path.join(direc, 'energies.dat'), energy)
    wavefuncs = np.hstack((xx.reshape((npoint, 1)), norm_wavefunc))
    np.savetxt(os.path.join(direc, 'wavefuncs.dat'), wavefuncs)

    # compute expectation values and uncertainty for the position
    exp_value = np.zeros(lasteigval-firsteigval+1)
    exp_value_sq = np.zeros(lasteigval-firsteigval+1)
    uncert_x = np.zeros(lasteigval-firsteigval+1)
    for ii in range(firsteigval-1, lasteigval):
        exp_value[ii] = delta*np.sum(norm_wavefunc[:, ii]**2*xx)
        exp_value_sq[ii] = delta*np.sum(norm_wavefunc[:, ii]**2*xx**2)
        uncert_x[ii] = np.sqrt(exp_value_sq[ii]-exp_value[ii]**2)

    # save expectation values and uncertainty for the position in expvalues.dat
    expvalues = np.array([exp_value, uncert_x])
    datexpvalues = expvalues.T
    np.savetxt(os.path.join(direc, 'expvalues.dat'), datexpvalues)

    return xypotential, energy
