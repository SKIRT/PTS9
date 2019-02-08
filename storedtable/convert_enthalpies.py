#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.convert_enthalpies Functions for converting enthalpies to stored table format
#
# The functions in this module can convert various kinds of enthalpy data to SKIRT stored table format.
# The stored tables produced by the functions in this module have:
# - a single axis: temperature T in K
# - a single quantity: specific enthalpy h specified as a value per unit of volume in J/m3
#
# Dividing the value listed in the table (in J/m3) by the bulk density of the material (in kg/m3)
# yields a specific enthalpy value per unit of mass (in J/kg).
#
# Providing the specific enthalpy per unit of volume implies that the information can be used for different
# bulk densities. However, some of the tables have been prepared assuming a specific bulk density value.
# In those cases, a proper dust model will be obtained only when using that particular bulk density value.
#
# In any case, the same bulk density value must be specified for the complete dust composition (i.e. also
# for calculating the dust mass from the size distribution) to obtain a consistent dust model.
#

# -----------------------------------------------------------------

import numpy as np
import pts.storedtable.io

# -----------------------------------------------------------------

## This function converts the column text files with enthalpies calculated according to Draine & Li (2001)
# to stored table format.
# The function expects a single input file path and a single output file path, plus a bulk density value in kg/m3
# posing as the second input file path. This value is used to convert the enthalpies from per-mass to per-volume.
# Refer to the notes in the input directory for more information on the data format.
def convertDraineEnthalpies(inFilePaths, outFilePaths):
    # get the arguments
    inFilePath = inFilePaths[0]
    outFilePath = outFilePaths[0]
    bulkdensity = float(inFilePaths[1].split("/")[-1])

    # read the input file
    T,h = np.loadtxt(inFilePath, usecols=(0,1), unpack=True)

    # convert from J/kg to J/m3
    h *= bulkdensity

    # write stored table
    pts.storedtable.io.writeStoredTable(outFilePath, ['T'], ['K'], ['log'], [T],  ['h'], ['J/m3'], ['log'], [h])

# -----------------------------------------------------------------

## This function converts the column text files with enthalpies for the TRUST dust model to stored table format.
# The function expects a single input file path and a single output file path, plus a bulk density value in kg/m3
# posing as the second input file path. This value is used to convert the enthalpies from per-mass to per-volume.
# Refer to the notes in the input directory for more information on the data format.
def convertTrustBenchmarkEnthalpies(inFilePaths, outFilePaths):
    # get the arguments
    inFilePath = inFilePaths[0]
    outFilePath = outFilePaths[0]
    bulkdensity = float(inFilePaths[1].split("/")[-1])

    # read the input file
    T,h = np.loadtxt(inFilePath, usecols=(0,1), skiprows=4, unpack=True)

    # convert units from erg/g to J/kg
    h *= 1e-4

    # convert from J/kg to J/m3
    h *= bulkdensity

    # write stored table
    pts.storedtable.io.writeStoredTable(outFilePath, ['T'], ['K'], ['log'], [T],  ['h'], ['J/m3'], ['log'], [h])

# -----------------------------------------------------------------

## This function converts "heat capacity" data files shipped with the DustEM code to enthalpies in stored table format.
# The function expects a sequence of input file paths and a sequence of corresponding output file paths.
# Refer to the notes in the input directory for more information on the input data format and the conversion process.
def convertDustemEnthalpies(inFilePaths, outFilePaths):
    for inFilePath, outFilePath in zip(inFilePaths, outFilePaths):
        # read the first two columns of the input file
        logTin,logCin = np.loadtxt(inFilePath, usecols=(0,1), skiprows=11, unpack=True)

        # interpolate the heat capacity values on a larger grid, to enable accurate integration
        logT = np.linspace(logTin[0], logTin[-1], 6000)
        logC = np.interp(logT, logTin, logCin)

        # perform the integration
        h = np.cumsum( np.log(10.) * 10.**(logC+logT) * (logT[1]-logT[0]) )

        # convert units from erg/cm3 to J/m3
        h *= 0.1

        # get the actual temperature grid
        T = 10.**logT

        # write stored table
        pts.storedtable.io.writeStoredTable(outFilePath, ['T'], ['K'], ['log'], [T],  ['h'], ['J/m3'], ['log'], [h])

# -----------------------------------------------------------------
