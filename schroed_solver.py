#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 16:52:05 2018

@author: jonaswillmann
"""

import sys
import numpy as np
import scipy
import scipy.linalg


def solver(fn):
    '''import parameters and solve the one-dimensional time-independent
    Schrödinger equation'''
    # import input file and test if file exists
    try:
        inp_all = open(fn, 'r')
    except FileNotFoundError:
        print("Input file {} not found".format(fn))
        sys.exit(1)
    else:
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
    npoint = float(inp_npoint)
    inp_firsteval = lines[2].split()[0]
    firsteval = float(inp_firsteval)
    inp_lasteval = lines[2].split()[1]
    lasteval = float(inp_lasteval)
    interpoltype = lines[3].split()[0]
    inp_nrinterpolpoints = lines[4].split()[0]
    nrinterpolpoints = float(inp_nrinterpolpoints)
    len_pot = len(lines)-5
    xpot = np.zeros(len_pot)
    ypot = np.zeros(len_pot)
    for ll in range(5, len_pot+5):
        xpot[ll - 5] = float(lines[ll].split()[0])
        ypot[ll - 5] = float(lines[ll].split()[1])
    inp_all().close()

    # potential interpolation
    if interpoltype == 'linear':
        xx = np.linspace(xmin, xmax, npoint)
        yy = np.interp(xx, xpot, ypot)
    elif interpoltype == 'polynomial':
        degree = nrinterpolpoints - 1
        polycoef = np.polyfit(xpot, ypot, degree)
        func = np.poly1d(polycoef)
        xx = np.linspace(xmin, xmax, npoint)
        yy = func(xx)
    elif interpoltype == 'cspline':
        cubfunc = np.interp1d(xpot, ypot, kind='cubic')
        xx = np.linspace(xmin, xmax, npoint)
        yy = cubfunc(xx)
    else:
        print('interpolation type not found')

    # put values for xPot and yPot into potential.dat file
    xy = np.array([xx, yy])
    data = xy.T
    np.savetxt('potential.dat', data)
    hh = abs((xmax-xmin)/npoint)

    # formulate matrix-problem for the discretised Schroedinger equation
    nn = int(npoint)
    aa = np.zeros((nn, nn))
    bb = 1/(mass*hh*hh)
    for ii in range(1, nn):
        aa[ii, ii-1] = -bb/2
    for ii in range(0, nn):
        aa[ii, ii] = bb+yy[ii]
    for ii in range(0, nn-1):
        aa[ii, ii+1] = -bb/2

    # compute eigenvalues and eigenvectors
    ee, psi = scipy.linalg.eigh(aa)

    # normalize wavefunctions
    nn = int(npoint)
    normpsi = np.zeros((nn, nn))
    for ii in range(0, nn):
        normpsi[ii] = 1/np.sqrt(hh)*psi[ii]

    # save eigenvalues and eigenvector in file
    np.savetxt("energies.dat", ee[int(firsteval-1):int(lasteval)])
    wavefuncs = np.hstack((xx.reshape((-1, 1)), normpsi))
    np.savetxt("wavefuncs.dat", wavefuncs[:, int(firsteval-1):int(lasteval+1)])

    # compute expectation values and uncertainty for the position
    exp_value = np.zeros(int(lasteval-firsteval)+1)
    exp_value_sq = np.zeros(int(lasteval-firsteval)+1)
    uncert_x = np.zeros(int(lasteval-firsteval)+1)
    for ii in range(int(lasteval-firsteval)+1):
        exp_value[ii] = hh*np.sum(normpsi[ii]*xx*normpsi[ii])
        exp_value_sq[ii] = hh*np.sum(normpsi[ii]*xx*xx*normpsi[ii])
        uncert_x[ii] = np.sqrt(exp_value_sq[ii]-exp_value[ii]*exp_value[ii])

    # save expectation values and uncertainty for the position in file
    expvalues = np.array([exp_value, uncert_x])
    datexpvalues = expvalues.T
    np.savetxt("expvalues.dat", datexpvalues)
