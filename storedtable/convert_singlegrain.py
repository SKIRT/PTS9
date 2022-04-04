#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.convert_singlegrain Functions for converting single-grain dust mix properties
#                                               to stored table format
#
# The functions in this module can convert various kinds of single-grain dust mix properties to SKIRT stored table
# format. Single-grain dust mixes define the complete dust mix through the optical properties of just a single
# "representative" grain. These dust mixes cannot be used to calculate NLTE dust emission.
#
# The stored tables for basic optical properties produced by the functions in this module
# have a single wavelength axis and the following quantities:
#  - sigmaabs (m2/H): absorption cross section per hydrogen nucleon (count, not mass)
#  - sigmasca (m2/H): scattering cross section per hydrogen nucleon (count, not mass)
#  - g (1): scattering anisotropy parameter, defined as the average of the scattering angle \f$<cos\theta>\f$
#  - mu (kg/H): dust mass per hydrogen nucleon; this is wavelength indepedent so all items have the same value
#
# The stored tables for Mueller matrix coefficients produced by the functions in this module
# have two axes (wavelength and scattering angle) and four quantities, one for each Mueller matrix coefficient.
# Because the Mueller matrix coefficients are always used relative to the first coefficient, their values
# are provided in arbitrary units. We use IAU convebtions for the signs of the coefficients.
#

# -----------------------------------------------------------------

import glob
import numpy as np
from .io import writeStoredTable

# -----------------------------------------------------------------

## This function converts the column text data for the single-grain MeanInterstellarDustMix to stored table format.
# The function expects a single input file path and a single output file path.
# Refer to the notes in the input directory for more information on the data format.
def convertMeanInterstellarOpticalProps(inFilePaths, outFilePaths):
    # read the input file
    w,a,g,se = np.loadtxt(inFilePaths[0], usecols=(0,1,2,3), skiprows=80, unpack=True)

    # convert units from micron to m and from cm^2/H to m^2/H
    w *= 1e-6
    se *= 1e-4

    # calculate absorption and scattering cross section from total cross section and albedo
    sa = (1-a)*se
    ss = a*se

    # replace asymmetry parameter values in x-ray range by analytical approximation 1 - lambda in micron
    mask = w < 1e-8
    g[mask] = 1 - w[mask]*1e6

    # determine dust mass
    mu = np.zeros_like(w) + 1.870e-29     # in kg/H

    # write stored table -- reverse order because input file has decreasing wavelengths
    writeStoredTable(outFilePaths[0], ['lambda'], ['m'], ['log'], [w[::-1]],
                          ['sigmaabs','sigmasca','g','mu'], ['m2/H','m2/H','1','kg/H'], ['log','log','lin','lin'],
                          [sa[::-1],ss[::-1],g[::-1],mu[::-1]])

# -----------------------------------------------------------------

## This function converts the column text data for the single-grain MeanDraineLiDustMix to stored table format.
# The function expects a single input file path and a single output file path.
# Refer to the notes in the input directory for more information on the data format.
def convertMeanDraineLiOpticalProps(inFilePaths, outFilePaths):
    # read the input file
    w,sa,ss,g = np.loadtxt(inFilePaths[0], usecols=(0,1,2,5), unpack=True)

    # convert units from micron to m and from cm^2/H to m^2/H
    w *= 1e-6
    sa *= 1e-4
    ss *= 1e-4

    # determine dust mass
    mu = np.zeros_like(w) + (5.4e-4 + 5.4e-4 + 1.8e-4 + 2.33e-3 + 8.27e-3) * 1.67262178e-27

    # write stored table
    writeStoredTable(outFilePaths[0], ['lambda'], ['m'], ['log'], [w],
                          ['sigmaabs','sigmasca','g','mu'], ['m2/H','m2/H','1','kg/H'], ['log','log','lin','lin'],
                          [sa,ss,g,mu])

# -----------------------------------------------------------------

## This function converts the column text data for the single-grain MeanZubkoDustMix to stored table format.
# The function expects a single input file path and a single output file path.
# Refer to the notes in the input directory for more information on the data format.
def convertMeanZubkoOpticalProps(inFilePaths, outFilePaths):
    # read the input file
    w,se,a,g = np.loadtxt(inFilePaths[0], usecols=(0,3,4,5), unpack=True)

    # convert units from micron to m and from cm^2/H to m^2/H
    w *= 1e-6
    se *= 1e-4

    # calculate absorption and scattering cross section from total cross section and albedo
    sa = (1-a)*se
    ss = a*se

    # determine dust mass
    mu = np.zeros_like(w) + 1.44e-29    # in kg/H

    # write stored table
    writeStoredTable(outFilePaths[0], ['lambda'], ['m'], ['log'], [w],
                          ['sigmaabs','sigmasca','g','mu'], ['m2/H','m2/H','1','kg/H'], ['log','log','lin','lin'],
                          [sa,ss,g,mu])

# -----------------------------------------------------------------

## This function creates the power-laws for the single-grain MeanIvezicBenchmarkDustMix in stored table format.
# The function expects no input file path (because the output is analytically generated) and a single output file path.
# Refer to the notes in the input directory for more information.
def createMeanIvezicBenchmarkOpticalProps(inFilePaths, outFilePaths):

    # wavelength grid in micron and corresponding optical properties
    w = np.array((1e-3, 1., 1e4))
    sa = np.array((1., 1., 1e-4))   # 1/lambda    for lambda>1
    ss = np.array((1., 1., 1e-16))  # 1/lambda^4  for lambda>1
    g = np.array((0., 0., 0.))      # isotropic

    # convert wavelength units from micron to m, cross sections from arbitrary units to reasonable range
    w *= 1e-6
    sa *= 2e-26
    ss *= 2e-26

    # determine dust mass so that kappa values have reasonable order of magnitude
    mu = np.zeros_like(w) + 1.5e-29   # arbitrary value, in kg/H

    # write stored table
    writeStoredTable(outFilePaths[0], ['lambda'], ['m'], ['log'], [w],
                          ['sigmaabs','sigmasca','g','mu'], ['m2/H','m2/H','1','kg/H'], ['log','log','lin','lin'],
                          [sa,ss,g,mu])

# -----------------------------------------------------------------

## This function converts the column text data for the single-grain MeanPascucciBenchmarkDustMix to stored table format.
# The function expects a single input file path and a single output file path.
# Refer to the notes in the input directory for more information on the data format.
def convertMeanPascucciBenchmarkOpticalProps(inFilePaths, outFilePaths):
    # read the input file
    w,ss,se = np.loadtxt(inFilePaths[0], usecols=(0,1,2), comments=';', unpack=True)

    # determine absorption cross section from scattering cross section and total cross section
    sa = se-ss

    # convert wavelength units from micron to m, cross sections from arbitrary units to reasonable range
    w *= 1e-6
    sa *= 7e-13
    ss *= 7e-13

    # setup isotropic scattering
    g = np.zeros_like(w)

    # determine dust mass so that kappa values have reasonable order of magnitude
    mu = np.zeros_like(w) + 1.5e-29   # arbitrary value, in kg/H

    # write stored table
    writeStoredTable(outFilePaths[0], ['lambda'], ['m'], ['log'], [w],
                          ['sigmaabs','sigmasca','g','mu'], ['m2/H','m2/H','1','kg/H'], ['log','log','lin','lin'],
                          [sa,ss,g,mu])

# -----------------------------------------------------------------

## This function converts the column text file with basic optical properties for the single-grain
# MeanPinteBenchmarkDustMix to stored table format.
# The function expects a single input file path and a single output file path.
# Refer to the notes in the input directory for more information on the data format.
def convertMeanPinteBenchmarkOpticalProps(inFilePaths, outFilePaths):
    # read the input file
    w,ka,ks = np.loadtxt(inFilePaths[0], usecols=(0,2,3), unpack=True)

    # convert units from micron to m and from cm^2/g to m^2/kg
    w *= 1e-6
    ka *= 0.1
    ks *= 0.1

    # set arbitrary dust mass per hydrogen nucleon, and convert kappa's to cross sections per hydrogen nucleon
    mu = np.zeros_like(w) + 1.5e-29   # arbitrary value, in kg/H
    sa = ka * mu
    ss = ks * mu

    # set all scattering asymmetry values to the (only) published value, i.e. the value for 1 micron;
    # these values won't be used for the benchmark because the Mueller matrix is used instead
    g = np.zeros_like(w) + 0.6296066

    # write stored table
    writeStoredTable(outFilePaths[0], ['lambda'], ['m'], ['log'], [w],
                          ['sigmaabs','sigmasca','g','mu'], ['m2/H','m2/H','1','kg/H'], ['log','log','lin','lin'],
                          [sa,ss,g,mu])

# -----------------------------------------------------------------

## This function converts the column text file with Mueller matrix coefficients for the single-grain
# MeanPolarizedPinteBenchmarkDustMix to stored table format.
# The function expects a single input file path and a single output file path. The input file contains data
# for just a single wavelength (1 micron). To enable the regular interpolation scheme in SKIRT, we actually
# put three wavelengths (with the same data) at the extreme ends of the wavelength range in the stored table.
# Refer to the notes in the input directory for more information on the data format.
def convertMeanPinteBenchmarkMuellerMatrix(inFilePaths, outFilePaths):
    # set the wavelength grid to include wavelengths at the extreme ends of the wavelength range in opacity.dat
    w = np.array((0.1,1.,3000.))

    # read the input file
    theta,S11,S12,S33,S34 = np.loadtxt(inFilePaths[0], usecols=(0,1,2,3,4), unpack=True)

    # convert units from micron to m and from degrees to radians
    # (Sxx values are always used relative to S11 so we don't need to worry about units)
    w *= 1e-6
    theta *= np.pi/180.

    # convert input values from Legacy convention to IAU convention
    S34 = -S34

    # copy the Mueller matrix coefficients for each wavelength
    S11 = np.tile(S11,(len(w),1))
    S12 = np.tile(S12,(len(w),1))
    S33 = np.tile(S33,(len(w),1))
    S34 = np.tile(S34,(len(w),1))

    # write stored table
    writeStoredTable(outFilePaths[0], ['lambda','theta'], ['m','rad'], ['log','lin'], [w, theta],
                          ['S11','S12','S33','S34'], ['1','1','1','1'], ['lin','lin','lin','lin'], [S11,S12,S33,S34])

# -----------------------------------------------------------------

## This function converts the column text file with basic optical properties for the single-grain
# MeanTrustBenchmarkDustMix to stored table format.
# The function expects a single input file path and a single output file path.
# Refer to the notes in the input directory for more information on the data format.
def convertMeanTrustBenchmarkOpticalProps(inFilePaths, outFilePaths):
    # read the input file
    w,se,a,g = np.loadtxt(inFilePaths[0], usecols=(0,3,4,5), skiprows=4, unpack=True)

    # calculate absorption and scattering cross section from total cross section and albedo
    sa = (1-a)*se
    ss = a*se

    # convert units from micron to m and from cm^2/H to m^2/H
    w *= 1e-6
    sa *= 1e-4
    ss *= 1e-4

    # determine dust mass
    mu = np.zeros_like(w) + 1.434e-29   # in kg/H

    # write stored table
    writeStoredTable(outFilePaths[0], ['lambda'], ['m'], ['log'], [w],
                          ['sigmaabs','sigmasca','g','mu'], ['m2/H','m2/H','1','kg/H'], ['log','log','lin','lin'],
                          [sa,ss,g,mu])

# -----------------------------------------------------------------

## This function converts the set of column text files with Mueller matrix coefficients for the single-grain
# MeanPolarizedTrustBenchmarkDustMix to stored table format.
# The function expects a single input file path and a single output file path. The input file path is actually
# interpreted as a template where the asterisk is replaced by the appropriate code for the scattering angle theta.
# Refer to the notes in the input directory for more information on the data format.
def convertMeanTrustBenchmarkMuellerMatrix(inFilePaths, outFilePaths):

    # get a list of the input file names, in order of increasing theta
    paths = sorted(glob.glob(inFilePaths[0]))

    # determine the scattering angle grid from the file names
    theta = np.array([ float(path.split('_')[-1][0:3]) for path in paths ])

    # determine the wavelength grid from the first file
    w = np.loadtxt(paths[0], usecols=(0,))

    # allocate arrays for the Mueller matrix coefficients
    S11 = np.zeros((len(w), len(theta)))
    S12 = np.zeros((len(w), len(theta)))
    S33 = np.zeros((len(w), len(theta)))
    S34 = np.zeros((len(w), len(theta)))

    # read the Mueller matrix coefficients from each file
    for index in range(len(paths)):
        S11[:,index],S12[:,index],S33[:,index],S34[:,index] = np.loadtxt(paths[index], usecols=(1,2,3,4), unpack=True)

    # convert units from micron to m and from degrees to radians
    # (Sxx values are always used relative to S11 so we don't need to worry about units)
    w *= 1e-6
    theta *= np.pi/180.

    # write stored table
    writeStoredTable(outFilePaths[0], ['lambda','theta'], ['m','rad'], ['log','lin'], [w, theta],
                          ['S11','S12','S33','S34'], ['1','1','1','1'], ['lin','lin','lin','lin'], [S11,S12,S33,S34])

# -----------------------------------------------------------------
