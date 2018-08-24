#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 14:10:55 2018

@author: jonaswillmann
"""

import numpy as np
import sys
import scipy
import scipy.linalg
fname = "infinit_potential_well.inp"


def input_schroed(fname):
    """import input file and test if file exists"""
    try:
        inp_all = open(fname, 'r')
    except FileNotFoundError:
        print("Input file {} not found".format(fname))
        sys.exit(1)
    else:
        inp_all = open(fname, 'r')
        return inp_all


"""import parameters from input file"""


lines = input_schroed(fname).readlines()
inp_mass = lines[0].split()[0]
mass = float(inp_mass)
inp_xMin = lines[1].split()[0]
xMin = float(inp_xMin)
inp_xMax = lines[1].split()[1]
xMax = float(inp_xMax)
inp_nPoint = lines[1].split()[2]
nPoint = float(inp_nPoint)
inp_firsteigval = lines[2].split()[0]
firsteigval = float(inp_firsteigval)
inp_lasteigval = lines[2].split()[1]
lasteigval = float(inp_lasteigval)
interpoltype = lines[3].split()[0]
inp_nrinterpolpoints = lines[4].split()[0]
nrinterpolpoints = float(inp_nrinterpolpoints)

len_pot = len(lines)-5
xPot = np.zeros(len_pot)
yPot = np.zeros(len_pot)
for ll in range(5, len_pot+5):
    xPot[ll - 5] = float(lines[ll].split()[0])
    yPot[ll - 5] = float(lines[ll].split()[1])
input_schroed(fname).close()


"""potential interpolation"""


if interpoltype == 'linear':
    x = np.linspace(xMin, xMax, nPoint)
    y = np.interp(x, xPot, yPot)
elif interpoltype == 'polynomial':
    degree = nrinterpolpoints - 1
    polycoef = np.polyfit(xPot, yPot, degree)
    func = np.poly1d(polycoef)
    x = np.linspace(xMin, xMax, nPoint)
    y = func(x)
elif interpoltype == 'cspline':
    cubfunc = np.interp1d(xPot, yPot, kind='cubic')
    x = np.linspace(xMin, xMax, nPoint)
    y = cubfunc(x)
else:
    print('interpolation type not found')


"""put values for xPot and yPot into potential.dat file"""
xy = np.array([x, y])
data = xy.T
np.savetxt('potential.dat', data)


Delta = abs((xMax-xMin)/nPoint)


def matrix_schroed(y):
    """formulate matrix-problem for the discretised Schroedinger equation"""
    N = int(nPoint)
    A = np.zeros((N, N))
    a = 1/(mass*Delta*Delta)

    for ii in range(1, N):
        A[ii, ii-1] = -a/2
    for ii in range(0, N):
        A[ii, ii] = a+y[ii]
    for ii in range(0, N-1):
        A[ii, ii+1] = -a/2

    return A


def solve_schroed(y):
    """compute eigenvalues and eigenvectors"""
    A = matrix_schroed(y)
    E, psi = scipy.linalg.eigh(A)

    return E, psi


E, psi = solve_schroed(y)


def norm_wavefuncs(psi):
    """normalize wavefunctions by dividung with 1/Delta"""
    N = int(nPoint)
    normpsi = np.zeros((N, N))
    for ii in range(0, N):
            normpsi[ii] = 1/np.sqrt(Delta)*psi[ii]

    return normpsi


"""save eigenvalues and eigenvector in file"""
normpsi = norm_wavefuncs(psi)
np.savetxt("energies.dat", E[int(firsteigval-1):int(lasteigval)])
wavefuncs = np.hstack((x.reshape((-1, 1)), normpsi))
np.savetxt("wavefuncs.dat", wavefuncs[:, int(firsteigval-1):int(lasteigval+1)])


def someshit(psi, x):
    exp_value = np.zeros(int(lasteigval-firsteigval)+1)
    exp_value_sq = np.zeros(int(lasteigval-firsteigval)+1)
    uncert_x = np.zeros(int(lasteigval-firsteigval)+1)
    for ii in range(int(lasteigval-firsteigval)+1):
        exp_value[ii] = Delta*np.sum(normpsi[ii]*x*normpsi[ii])
        exp_value_sq[ii] = Delta*np.sum(normpsi[ii]*x*x*normpsi[ii])
        uncert_x[ii] = np.sqrt(exp_value_sq[ii]-exp_value[ii]*exp_value[ii])

    return exp_value, exp_value_sq, uncert_x


exp_value, exp_value_sq, uncert_x = someshit(psi, x)
someshittyshit = np.array([exp_value, uncert_x])
fuckshit = someshittyshit.T
np.savetxt("expvalues.dat", fuckshit)
