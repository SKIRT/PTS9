#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.convert_opticalprops Functions for converting optical properties to stored table format
#
# The functions in this module can convert various kinds of optical property data to SKIRT stored table format.
# The stored tables for basic optical properties produced by the functions in this module
# have a two axes:
#  - lambda (m): wavelength
#  - a (m): radius of the dust grain, which is assumed to be spherical
# and the following quantities, defined on those axes:
#  - Qabs (1): absorption efficiency
#  - Qsca (1): scattering efficiency
#  - g (1): scattering anisotropy parameter, defined as the average of the scattering angle \f$<cos\theta>\f$
#

# -----------------------------------------------------------------

import numpy as np
import glob
from .io import writeStoredTable
from .tokenizedfile import TokenizedFile

# -----------------------------------------------------------------

## This function converts optical properties listed in a single text file (in the format described below)
# to stored table format. Several options provide extra flexibility to support various existing formats.
#
# In addition to a single output file path for the stored table, the function expects five input argument
# strings, in order of occurrence: input-file-path, reverse, skip1, skip2, skip3. The last four arguments
# are booleans. The string "true" (case insensitive) represents True, anything else represents False.
#
# The input file should have a simple text format as described here. Any initial lines that start with a
# '#' character are considered to be part of a header and are thus ignored.
# The first number on the first non-header line specifies the number of grain size grid points \f$N_a\f$;
# the first number on the second non-header line specifies the number of wavelength grid points \f$N_\lambda\f$.
#
# Subsequently there must be \f$N_a\f$ blocks containing \f$N_\lambda+1\f$ lines each.
# The first number on the first line in each block provides the grain size value
# for this block (in micron); blocks must be in order of increasing grain size. Subsequent
# lines in the block provide information for the various wavelengths at this grain size. If
# the \em reverse flag is false, these lines are in order of increasing wavelength, otherwise
# they are in order of decreasing wavelength. All data lines in a block must have 4 to 7
# columns in the following order: \f$X_1\f$, \f$\lambda\f$, \f$X_2\f$, \f$Q^\text{abs}\f$,
# \f$Q^\text{sca}\f$, \f$X_3\f$, \f$g\f$, where the presence or absence of each of the dummy
# \f$X_i\f$ fields is indicated by the corresponding \em skip flag (true means the field
# is present, false means it is absent). The wavelength must be given in micron; the other
# three values are dimensionless.
#
# For all lines discussed above, any additional information at the end of the line is ignored.
#
def convertGenericOpticalProps(inFilePaths, outFilePaths):
    # convert options to Boolean
    reverse, skip1, skip2, skip3 = [ b.lower().endswith("/true") for b in inFilePaths[1:5] ]

    # open the input file and skip the header
    infile = TokenizedFile(open(inFilePaths[0]))
    infile.skipHeaderLines()

    # read the grid size
    Na = int(infile.next())
    infile.skipToEndOfLine()
    Nlambda = int(infile.next())
    infile.skipToEndOfLine()

    # allocate arrays
    w = np.zeros(Nlambda)
    a = np.zeros(Na)
    Qabs = np.zeros((Na,Nlambda))
    Qsca = np.zeros((Na,Nlambda))
    g = np.zeros((Na,Nlambda))

    # read the data blocks
    for i in range(Na):
        a[i] = float(infile.next())
        infile.skipToEndOfLine()

        for k in range(Nlambda):
            if skip1: infile.next()
            w[k] = float(infile.next())
            if skip2: infile.next()
            Qabs[i,k] = float(infile.next())
            Qsca[i,k] = float(infile.next())
            if skip3: infile.next()
            g[i,k] = float(infile.next())
            infile.skipToEndOfLine()

    # reverse the wavelengths if needed
    if reverse:
        w = np.flip(w, 0)
        Qabs = np.flip(Qabs, 1)
        Qsca = np.flip(Qsca, 1)
        g = np.flip(g, 1)

    # write stored table
    writeStoredTable(outFilePaths[0], ['a','lambda'], ['m','m'], ['log','log'], [a*1e-6,w*1e-6],
            ['Qabs','Qsca','g'], ['1','1','1'], ['log','log','lin'], [Qabs,Qsca,g])

# -----------------------------------------------------------------

## This function converts optical properties given in DustEM format to stored table format. It expects
# three input file paths and a single output file path. The input files should respectively contain:
#  - the wavelength grid \f$\lambda_k\f$, in increasing order;
#  - the grain sizes \f$a_i\f$ plus the absorption and scattering efficiencies
#    \f$Q_{k,i}^{\text{abs}}\f$ and \f$Q_{k,i}^{\text{sca}}\f$ on the two-dimensional grid;
#  - the grain sizes \f$a_i\f$ plus the scattering phase function asymmetry parameter
#    \f$g_{k,i}\f$ on the same grid; the grain sizes values must be identical in both files.
#
# The files should have a simple text format as described here. For all files, any initial
# lines that start with a # character are considered to be part of a header and are thus
# ignored.
#
# In the wavelength grid file, the first number on the first non-header line specifies the number of
# wavelength grid points \f$N_\lambda\f$. The remaining numbers in the file (separated by whitespace,
# possibly including line endings) specify the wavelength grid points \f$\lambda_k\f$ in (micron).
#
# The efficiencies file contains three blocks, each optionally starting with a header
# consisting of lines that start with a # character. In the first block, the first number on
# the first non-header line specifies the number of grain size grid points \f$N_a\f$. The
# second line has \f$N_a\f$ columns specifying the grain sizes \f$a_i\f$ (in micron), in
# increasing order. The second block contains \f$N_\lambda\f$ non-header lines. Each line
# specifies the \f$N_a\f$ absorption efficiencies \f$Q_{k,i}^{\text{abs}}\f$ corresponding to
# line (wavelength) \f$k\f$ and column (grain size) \f$i\f$. Similarly, the third block
# contains \f$N_\lambda\f$ non-header lines. Each line now specifies the \f$N_a\f$ scattering
# efficiencies \f$Q_{k,i}^{\text{sca}}\f$ corresponding to line (wavelength) \f$k\f$ and
# column (grain size) \f$i\f$.
#
# The scattering assymetry parameter file similarly contains two blocks. The first block
# describes the grain size grid and should be identical to the first block in the
# efficiencies file. The second block again contains \f$N_\lambda\f$ non-header lines. Each
# line now specifies the \f$N_a\f$ scattering assymetry parameters \f$g_{k,i}\f$
# corresponding to line (wavelength) \f$k\f$ and column (grain size) \f$i\f$. */
#
def convertDustemOpticalProps(inFilePaths, outFilePaths):

    # ------------ read wavelengths file ------------

    infile = TokenizedFile(open(inFilePaths[0]))

    # skip header lines, read the wavelength grid size, and read the wavelengths
    infile.skipHeaderLines()
    Nlambda = int(infile.next())
    w = np.array([ float(infile.next()) for k in range(Nlambda) ])

    # ------------ read efficiencies file ------------

    infile = TokenizedFile(open(inFilePaths[1]))

    # first block: skip header lines, read the grain size grid size, and read the grain sizes
    infile.skipHeaderLines()
    Na = int(infile.next())
    a = np.array([ float(infile.next()) for i in range(Na) ])

    # second block: skip header lines, and read the absorption efficiencies
    infile.skipHeaderLines()
    Qabs = np.array([ [ float(infile.next()) for i in range(Na) ] for k in range(Nlambda) ]).T

    # third block: skip header lines, and read the scattering efficiencies
    infile.skipHeaderLines()
    Qsca = np.array([ [ float(infile.next()) for i in range(Na) ] for k in range(Nlambda) ]).T

    # ------------ read scattering asymmetry file ------------

    infile = TokenizedFile(open(inFilePaths[2]))

    # first block: skip header lines, read the grain size grid size, and skip the grain sizes
    infile.skipHeaderLines()
    if Na != int(infile.next()):
        raise ValueError("Efficiencies and scattering asymmetry files have different number of grain sizes")
    for i in range(Na): infile.next()

    # second block: skip header lines, and read the scattering asymmetry parameters
    infile.skipHeaderLines()
    g = np.array([ [ float(infile.next()) for i in range(Na) ] for k in range(Nlambda) ]).T

    # ------------ write stored table ------------

    writeStoredTable(outFilePaths[0], ['a','lambda'], ['m','m'], ['log','log'], [a*1e-6,w*1e-6],
            ['Qabs','Qsca','g'], ['1','1','1'], ['log','log','lin'], [Qabs,Qsca,g])

# -----------------------------------------------------------------

## This function converts optical properties given in the format provided by Michiel Min to stored table format.
# In addition to a single output file path for the stored table, the function expects two input arguments: the
# input filename template, in which the asterisk will be replaced by the grain size to obtain the actual filenames,
# and the bulk density of the grain material.
#
# There is an input file for each grain size; the grain size (in micron) is encoded in the file name.
# Each file has the following columns:
#    wavelength (micron), kappa_ext (cm^2/g), kappa_abs (cm^2/g), kappa_scat (cm^2/g), asymmetry parameter.
# For solid spherical grains we have:  Q = kappa * rho_bulk * 4/3*a.
#
def convertMinOpticalProps(inFilePaths, outFilePaths):
    # discover the input files and parse the bulk density
    infiles = sorted(glob.glob(inFilePaths[0]))
    rhobulk = float(inFilePaths[1].split("/")[-1])

    # get the grain size grid from the filenames and convert to m
    a = np.array([ float(infile.split("/")[-1].split("_")[-1][0:-8]) for infile in infiles ]) * 1e-6

    # get the wavelength grid from the first file and convert to m
    w = np.loadtxt(infiles[0], usecols=(0,), unpack=True) * 1e-6

    # allocate arrays
    Qabs = np.zeros((len(w), len(a)))
    Qsca = np.zeros((len(w), len(a)))
    g = np.zeros((len(w), len(a)))

    # read all data into memory
    for i in range(len(a)):
        Qabs[:,i], Qsca[:,i], g[:,i] = np.loadtxt(infiles[i], usecols=(2,3,4), unpack=True)

    # convert kappa's in cm^2/g to m2/kg and then to Q's  (need a in m!)
    Qabs *= 0.1 * (4./3.*a*rhobulk)
    Qsca *= 0.1 * (4./3.*a*rhobulk)

    # write stored table
    writeStoredTable(outFilePaths[0], ['a','lambda'], ['m','m'], ['log','log'], [a,w],
            ['Qabs','Qsca','g'], ['1','1','1'], ['log','log','lin'], [Qabs.T,Qsca.T,g.T])

# -----------------------------------------------------------------

## This function converts optical properties given in STOKES format to stored table format. It expects
# a single input file paths and two output file paths: one for the basic optical properties and one for
# the Mueller matrix coeffients.
#
# The input file should have the text format as used by the STOKES code version 2.06. The file contains
# the absorption and scattering coefficients as a function of grain size and wavelength, and Mueller
# matrix coeffients as a function of grain size, wavelength, and scattering angle.
#
def convertStokesPolarizedOpticalProps(inFilePaths, outFilePaths):
    # open the input file
    infile = TokenizedFile(open(inFilePaths[0]))

    # skip header lines and read the grid sizes (sizes are given as n-1 in input file)
    for h in range(int(infile.next())): infile.skipLine()
    Na = int(infile.next()) + 1
    infile.skipToEndOfLine()
    Nlambda = int(infile.next()) + 1
    infile.skipToEndOfLine()
    Ntheta = int(infile.next()) + 1
    infile.skipToEndOfLine()
    for h in range(3): infile.skipLine()

    # allocate arrays
    w = np.zeros(Nlambda)
    a = np.zeros(Na)
    theta = np.zeros(Ntheta)
    Qabs = np.zeros((Na, Nlambda))
    Qsca = np.zeros((Na, Nlambda))
    g = np.zeros((Na, Nlambda))             # g remains zero because it is unused for polarized dust mixes
    S11 = np.zeros((Na, Nlambda, Ntheta))
    S12 = np.zeros((Na, Nlambda, Ntheta))
    S33 = np.zeros((Na, Nlambda, Ntheta))
    S34 = np.zeros((Na, Nlambda, Ntheta))

    # read the data
    for i in range(Na):
        a[i] = float(infile.next())
        infile.skipToEndOfLine()

        for k in range(Nlambda-1,-1,-1):
            if infile.next()!='#': raise ValueError("Expected line with wavelength column titles")
            infile.skipToEndOfLine()
            w[k] = float(infile.next())
            Qabs[i,k] = float(infile.next())
            Qsca[i,k] = float(infile.next())
            if infile.next()!='#': raise ValueError("Expected line with angle column titles")
            infile.skipToEndOfLine()

            for t in range(Ntheta):
                theta[t] = float(infile.next())
                S11[i,k,t] = float(infile.next())
                S12[i,k,t] = float(infile.next())
                S33[i,k,t] = float(infile.next())
                S34[i,k,t] = float(infile.next())

    # convert units
    w *= 1e-6
    a *= 1e-6
    theta *= np.pi / 180.

    # write stored table with absorption and scattering coefficients (plus dummy anisotropy parameter)
    writeStoredTable(outFilePaths[0], ['a','lambda'], ['m','m'], ['log','log'], [a,w],
            ['Qabs','Qsca','g'], ['1','1','1'], ['log','log','lin'], [Qabs,Qsca,g])

    # write stored table with Mueller matrix coeffients
    writeStoredTable(outFilePaths[1], ['a','lambda','theta'], ['m','m','rad'], ['log','log','lin'], [a,w,theta],
            ['S11','S12','S33','S34'], ['1','1','1','1'], ['lin','lin','lin','lin'], [S11,S12,S33,S34])

# -----------------------------------------------------------------
