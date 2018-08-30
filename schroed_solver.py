#!/usr/bin/env python3
"""Solver for the one-dimensional time-independent Schroedinger equation for
any potential which can be read from an input file"""

import os.path
import sys
import numpy as np
import scipy
from scipy.interpolate import CubicSpline


def _inputdata(direc):
    """Import parameters from an input file.

    Args:
        direc: directory of the inputfiles

    Returns:
        xmin: minimum x-value
        xmax: maximum x-value
        npoint: number of x-values
        xpot: x-values for the potential before interpolation
        ypot: y-values for the potential before interpolation
        interpoltyp: type of interpolation
        nrinterpolpoints: number of interpolation points
        mass: particle mass
        firsteigval: first eigenvalue
        lasteigval: last eigenvalue
    """

    # open inputfile
    fn = os.path.join(direc, 'schroedinger.inp')
    inp_all = open(fn, 'r')

    # read parameters from inputfile
    lines = inp_all.readlines()
    mass = float(lines[0].split()[0])
    xmin = float(lines[1].split()[0])
    xmax = float(lines[1].split()[1])
    npoint = int(lines[1].split()[2])
    firsteigval = int(lines[2].split()[0])
    lasteigval = int(lines[2].split()[1])
    interpoltype = lines[3].split()[0]
    nrinterpolpoints = int(lines[4].split()[0])
    xpot = np.zeros(len(lines)-5)
    ypot = np.zeros(len(lines)-5)
    for ii in range(5, len(lines)):
        xpot[ii - 5] = float(lines[ii].split()[0])
        ypot[ii - 5] = float(lines[ii].split()[1])

    return (xmin, xmax, npoint, xpot, ypot, interpoltype, nrinterpolpoints,
            mass, firsteigval, lasteigval)


def _interp(xmin, xmax, npoint, xpot, ypot, interpoltype, nrinterpolpoints):
    """Interpolate given values for the potential.

    Args:
        xmin: minimum x-value
        xmax: maximum x-value
        npoint: number of x-values
        xpot: x-values for the potential before interpolation
        ypot: y-values for the potential before interpolation
        interpoltyp: type of interpolation
        nrinterpolpoints: number of interpolation points

    Returns:
        xx: vector for the position
        pot: interpolated potential
        xypotential: array with x- and y-values of the interpolated potential
    """

    # read interpolation type and interpolate the potential
    xx = np.linspace(xmin, xmax, npoint)
    if interpoltype == 'linear':
        pot = np.interp(xx, xpot, ypot)
    elif interpoltype == 'polynomial':
        degree = int(nrinterpolpoints - 1)
        coef = np.polyfit(xpot, ypot, degree)
        polf = np.poly1d(coef)
        pot = polf(xx)
    elif interpoltype == 'cspline':
        cubicf = CubicSpline(xpot, ypot, bc_type='natural')
        pot = cubicf(xx)
    else:
        print('interpolation type not found')
        sys.exit(1)

    xypotential = np.array([xx, pot]).T

    return xx, pot, xypotential


def _solve(xx, xmin, xmax, npoint, mass, firsteigval, lasteigval, pot):
    """A matrix for the discretised problem will be solved to receive the
    eigenvalues and wavefunctions. The wavefunctions will be normed.

    Args:
        xx: vector for the position
        xmin: minimum x-value
        xmax: maximum x-value
        npoint: number of x-values
        mass: particle mass
        firsteigval: first eigenvalue
        lasteigval: last eigenvalue
        pot: interpolated potential

    Returns:
        energy: array with the eigenvalues for the matrix of the problem
        wavefuncs: array with the eigenvectors for the matrix of the problem
        norm_wavefunc: normalized wavefunctions
    """

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
    energy, wavefunc = scipy.linalg.eigh(matrix)
    energy = energy[0:lasteigval-firsteigval+1]

    # normalize wavefunctions
    norm_wavefunc = 1/np.sqrt(delta)*wavefunc

    # put wavefunctions in an array with the x-values
    wavefuncs = np.hstack((xx.reshape((npoint, 1)), norm_wavefunc))
    wavefuncs = wavefuncs[:, firsteigval-1:lasteigval+1]

    return energy, wavefuncs, norm_wavefunc, delta


def _expval(xx, xmin, xmax, npoint, firsteigval, lasteigval, norm_wavefunc):
    """The expectation values and the uncertainty for the position will be
    calculated.

    Args:
        xx: vector for the position
        xmin: minimum x-value
        xmax: maximum x-value
        firsteigval: first eigenvalue
        lasteigval: last eigenvalue
        norm_wavefunc: normalized wavefunctions

    Returns:
        expvalues: array with expectation values and position-uncertainty
    """

    # compute expectation values and uncertainty for the position
    delta = abs((xmax-xmin)/(npoint))
    exp_value = np.zeros(lasteigval-firsteigval+1)
    exp_value_sq = np.zeros(lasteigval-firsteigval+1)
    uncert_x = np.zeros(lasteigval-firsteigval+1)
    for ii in range(firsteigval-1, lasteigval):
        exp_value[ii] = delta*np.sum(norm_wavefunc[:, ii]**2*xx)
        exp_value_sq[ii] = delta*np.sum(norm_wavefunc[:, ii]**2*xx**2)
        uncert_x[ii] = np.sqrt(exp_value_sq[ii]-exp_value[ii]**2)

    expvalues = np.array([exp_value, uncert_x]).T

    return expvalues


def solver(direc):
    """Solver function for the one-dimensional time-independent Schroedinger
    equation. All results will be saved in .dat files.

    Args:
        direc: directory of the inputfiles

    Returns:
        xypotential: array with x- and y-values of the interpolated potential
        energy: array with the eigenvalues for the matrix of the problem
    """

    # redefine the variables of the functions
    inputargs1 = _inputdata(direc)[0:7]
    inputargs2 = _inputdata(direc)[7:]
    xx = _interp(*inputargs1)[0]
    pot = _interp(*inputargs1)[1]
    norm_wavefunc = _solve(xx, *inputargs1[0:3], *inputargs2, pot)[2]

    # load results from functions
    xypotential = _interp(*inputargs1)[2]
    energy = _solve(xx, *inputargs1[0:3], *inputargs2, pot)[0]
    wavefuncs = _solve(xx, *inputargs1[0:3], *inputargs2, pot)[1]
    expvalues = _expval(xx, *inputargs1[0:3], *inputargs2[1:], norm_wavefunc)

    # save results in .dat files
    np.savetxt(os.path.join(direc, 'potential.dat'), xypotential)
    np.savetxt(os.path.join(direc, 'energies.dat'), energy)
    np.savetxt(os.path.join(direc, 'wavefuncs.dat'), wavefuncs)
    np.savetxt(os.path.join(direc, 'expvalues.dat'), expvalues)

    return xypotential, energy
