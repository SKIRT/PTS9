#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.convert_sed Functions for converting SEDs and SED families to stored table format
#
# The functions in this module can convert various kinds of SED and SED family data to SKIRT stored table format.

# -----------------------------------------------------------------

import astropy.io.fits as fits
import glob
import gzip
import numpy as np
import scipy.io
from .io import writeStoredTable
from .tokenizedfile import TokenizedFile

# -----------------------------------------------------------------

## This function calculates an SED representing a simple model for the spectral energy distribution of a quasar;
# see Stalevski et al. (2012, MNRAS, 420, 2756–2772) and Schartmann et al. (2005, A&A, 437, 861–881).
# The SED is defined in the wavelength range between 0.001 \f$\mu\f$m and 1000 \f$\mu\f$m and is characterized by
#    \f[ S_\lambda \propto \begin{cases} \; \lambda^{1/5} &
#    \text{if $0.001~\mu{\text{m}}<\lambda<0.01~\mu{\text{m}}$} \\ \; \lambda^{-1} & \text{if
#    $0.01~\mu{\text{m}}<\lambda<0.1~\mu{\text{m}}$} \\ \; \lambda^{-3/2} & \text{if
#    $0.1~\mu{\text{m}}<\lambda<5~\mu{\text{m}}$} \\ \; \lambda^{-4} & \text{if
#    $5~\mu{\text{m}}<\lambda<1000~\mu{\text{m}}$.} \end{cases} \f]
#
# The function expects no input file names and a single output file name. There is no input file because the
# function calculates the SED analytically. Because SEDs are log-log interpolated in SKIRT, and the segments of
# the QuasarSED are power-laws, we can restrict the tabulation to the wavelength grid points separating the segments.
#
# The QuasarSED is arbitrarily normalized so that the value at 0.001 micron equals 1 W/m.
#
def createQuasarSED(inFilePaths, outFilePaths):
    # wavelength grid in micron
    w = np.array((0.001, 0.01, 0.1, 5., 1000.))
    L = np.zeros_like(w)

    # calculate SED at grid points, with arbitrary first value
    L[0] = 1.
    L[1] = L[0] * (w[1]/w[0])**0.2
    L[2] = L[1] * (w[2]/w[1])**-1.0
    L[3] = L[2] * (w[3]/w[2])**-1.5
    L[4] = L[3] * (w[4]/w[3])**-4.0

    # write stored table
    writeStoredTable(outFilePaths[0], ['lambda'], ['m'], ['log'], [w*1e-6],
                                           ['Llambda'], ['W/m'], ['log'], [L])

# -----------------------------------------------------------------

## This function converts a regular text column file with an optional header and two columns
# (wavelength in micron and specific luminosity in W/micron) to stored table file format.
# The function expects a sequence of input file paths and corresponding output file paths.
def convertTextSEDinWattsPerMicron(inFilePaths, outFilePaths):
    for inFilePath, outFilePath in zip(inFilePaths, outFilePaths):
        w, L = np.loadtxt(inFilePath, unpack=True)
        writeStoredTable(outFilePath, ['lambda'], ['m'], ['log'], [w*1e-6],
                                           ['Llambda'], ['W/m'], ['log'], [L*1e6])

# -----------------------------------------------------------------

## This function converts data representing the family of Bruzual & Charlot SEDs for single stellar populations,
# parameterized on metallicity and age (Bruzual & Charlot 2003, RAS 344, 1000-1026), to stored table format.
# The function expects a sequence of input file paths and corresponding output file paths. Each input file path
# is actually interpreted as a template where the string "MM" will be replaced by the appropriate metallicity code.
#
# The stored table output file contains specific luminosities in W/m normalized on an *initial* stellar population
# mass of 1 solar mass, in function of wavelength and parameterized on metallicity and age.
#
def convertBruzualCharlotSEDFamily(inFilePaths, outFilePaths):
    for inFilePath, outFilePath in zip(inFilePaths, outFilePaths):

        # read the age and wavelength grids from one of the files
        with gzip.open(inFilePath.replace("MM","22"), 'rt') as infile:
            tokens = TokenizedFile(infile)
            Nt = int(tokens.next())
            t = np.array([ float(tokens.next()) for p in range(Nt) ])
            while tokens.next()!="Padova": pass
            for i in range(3): tokens.skipLine()
            Nw = int(tokens.next())
            w = np.array([ float(tokens.next()) for k in range(Nw) ])

        # initialize the metallicity grid and the corresponding filename codes
        Z = np.array((0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05))
        NZ = len(Z)
        Zcode = [ "22", "32", "42", "52", "62", "72" ]

        # allocate hypercube for the luminosities
        L = np.zeros((Nw,NZ,Nt))    # indices k, m, p

        # read the luminosities from each of the files
        for m in range(NZ):
            with gzip.open(inFilePath.replace("MM",Zcode[m]), 'rt') as infile:
                tokens = TokenizedFile(infile)
                # skip the ages and wavelengths
                while tokens.next()!="Padova": pass
                for i in range(3): tokens.skipLine()
                for i in range(Nw+1): tokens.next()
                # iterate over ages
                for p in range(Nt):
                    # read luminosities for each age
                    if int(tokens.next()) != Nw:
                        raise ValueError("Number of luminosities does not match number of wavelengths")
                    for k in range(Nw):
                        L[k, m, p] = float(tokens.next())
                    # skip intervening dummy data
                    for i in range(int(tokens.next())): tokens.next()

        # write stored table, converting from input units
        # wavelengths (Angstrom), metallicities (dimensionless), ages (years), luminosities (Lsun/Angstrom)
        Angstrom = 1e-10
        Lsun = 3.839e26
        t[0] = 1  # avoid zero values in a logarithmically scaled axis
        writeStoredTable(outFilePath,
                              ['lambda','Z','t'], ['m','1','yr'], ['log','log','log'], [w*Angstrom,Z,t],
                              ['Llambda'], ['W/m'], ['log'], [L*(Lsun/Angstrom)])

# -----------------------------------------------------------------

## This function converts data representing the family of SEDs for simple stellar populations (SSPs) according to
# the model of Maraston (1998, MNRAS, 300, 872–892) to stored table format.
# The function expects a sequence of input file paths and corresponding output file paths. Each input file path
# is actually interpreted as a template where the string "*" will be replaced by the appropriate metallicity code.
#
# The stored table output file contains specific luminosities in W/m normalized on an *initial* stellar population
# mass of 1 solar mass, in function of wavelength and parameterized on metallicity and age.
#
def convertMarastonSEDFamily(inFilePaths, outFilePaths):
    for inFilePath, outFilePath in zip(inFilePaths, outFilePaths):

        # read and concatenate all column text files for this model
        # (all info is in the files; the file names themselves are thus redundant)
        tv, Zv, wv, Lv = np.concatenate([ np.loadtxt(path) for path in glob.glob(inFilePath) ]).T

        # determine the grid points for each axis as the sorted list of unique values
        # also get the indices in the unique array that can be used to reconstruct the original vector
        w, wi = np.unique(wv, return_inverse=True)
        Z, Zi = np.unique(Zv, return_inverse=True)
        t, ti = np.unique(tv, return_inverse=True)

        # allocate hypercube for the luminosities and copy the input luminosities to it
        # values that are not provided in the input remain at zero
        L = np.zeros((len(w),len(Z),len(t)))
        L[wi,Zi,ti] = Lv

        # convert from input units: Age(Gyr)  [Z/H]  lambda(AA)  L_{lambda}(erg/s/AA)
        w *= 1e-10
        Z = 0.02 * 10**Z
        t *= 1e9
        L *= 1e3

        # write stored table
        writeStoredTable(outFilePath,
                              ['lambda','Z','t'], ['m','1','yr'], ['log','log','log'], [w,Z,t],
                              ['Llambda'], ['W/m'], ['log'], [L])

# -----------------------------------------------------------------

## This function converts data representing the family of Starburst99 SEDs for single stellar populations
# (Leitherer et al. 1999, ApJS, 123, 3), assuming the Kroupa initial mass function (Kroupa 2001, MNRAS, 322, 231),
# and parameterized on metallicity and age, to stored table format. The function expects a single input file path
# pointing to the compressed FITS file containing the data, and a single output file path.
#
# The stored table output file contains specific luminosities in W/m normalized on an *initial* stellar population
# mass of 1 solar mass, in function of wavelength and parameterized on metallicity and age.
#
def convertStarburst99SEDFamily(inFilePaths, outFilePaths):
    # open the fits file
    hdul = fits.open(inFilePaths[0])

    # get the axes information
    table = hdul['AXES'].data
    w = table['lambda'][table['lambda']>0]
    Z = table['metallicity'][table['metallicity']>0]
    t = table['time'][table['time']>0]

    # get the luminosities (stored in file as log10)
    L = 10**hdul['SED'].data

    # write stored table
    writeStoredTable(outFilePaths[0],
                          ['lambda','Z','t'], ['m','1','yr'], ['log','log','log'], [w,Z,t],
                          ['Llambda'], ['W/m'], ['log'], [L])

# -----------------------------------------------------------------

## This function converts data representing a BPASS family of SED templates for stellar populations
# (Eldridge, Stanway et al, 2017, PASA 34, 58; Stanway & Eldridge, 2018, MNRAS, 479, 75) to stored table format.
# The function expects a sequence of input file paths and corresponding output file paths. Each input file path
# is actually interpreted as a template where the string "*" will be replaced by the appropriate metallicity code.
#
# The stored table output file contains specific luminosities in W/m normalized on an *initial* stellar population
# mass of 1 solar mass, as a function of wavelength and parameterized on metallicity and age.
#
def convertBpassSEDFamily(inFilePaths, outFilePaths):
    for inFilePath, outFilePath in zip(inFilePaths, outFilePaths):

        # read the wavelength grid from one of the files (and convert from Angstrom to meter)
        w = np.loadtxt(inFilePath.replace("*", "001"), usecols=0) * 1e-10

        # initialize the metallicity grid and the corresponding filename codes
        Z = np.array((1e-5, 1e-4, 0.001, 0.002, 0.003, 0.004, 0.006, 0.008, 0.01, 0.014, 0.02, 0.03, 0.04))
        Zcode = ["em5", "em4", "001", "002", "003", "004", "006", "008", "010", "014", "020", "030", "040"]

        # initialize the age grid (in years)
        t = np.array([10**(6+c/10) for c in range(51)])

        # allocate hypercube for the luminosities
        L = np.zeros((len(w), len(Z), len(t)))

        # read the luminosities from each of the files
        for m in range(len(Z)):
            L[:, m, :] = np.loadtxt(inFilePath.replace("*", Zcode[m]), usecols=range(1, 52))

        # convert from input units: Lsun/Angstrom normalized to 1e6 Msun
        #          to output units: W/m normalized to 1 Msun
        L *= 3.848e26 * 1e10 * 1e-6

        # write stored table
        writeStoredTable(outFilePath,
                              ['lambda','Z','t'], ['m','1','yr'], ['log','log','log'], [w,Z,t],
                              ['Llambda'], ['W/m'], ['log'], [L])

# -----------------------------------------------------------------

## This function converts data representing the family of MAPPINGS III starburst template SEDs, parameterized on
# metallicity, compactness, ISM pressure and PDR covering factor, as described in Groves et al. (2008, ApJS,176,438),
# to stored table format. The function expects a single input file path pointing to the IDL save file containing
# the data, and a single output file path.
#
# The stored table output file contains specific luminosities in W/m scaled to a star formation rate
# of 1 Msun/yr (over 10^7 yrs).
#
# The input data must be converted for SKIRT as follows:
#   - lambda: micron --> * 1e-6  (to m)
#   - Z: relative to solar, with Zsun = 0.0122 as in Asplund et al. (2005) --> Z = Zrel * Zsun
#   - logC: dimensionless --> no conversion
#   - P: given as log10(P/k) ('ponk') with k Boltzmann's constant, and in K/cm^3 --> P = k * 1e6 * 10**ponk (to Pa)
#   - fPDR: dimensionless fraction -> no conversion
#   - Llambda: erg/s/Hz  -->  * 1e-7 * c / lambda**2   (to W/m)
#
def convertMappingsSEDFamily(inFilePaths, outFilePaths):
    # load the input file
    data = scipy.io.readsav(inFilePaths[0])

    # get the axis grid points (explicitly use float() because some are defined as strings)
    w = np.array(list(map(float,data['wave_micron']))) * 1e-6
    Z = np.array(list(map(float,data['metal']))) * 0.0122
    logC = np.array(list(map(float,data['cparam'])))
    P = 1.3806488e-23 * 1e6 * 10**np.array(list(map(float,data['ponk'])))
    fPDR = np.array((0.,1.))

    # get the luminosities hypercube
    # - convert to W/m
    # - rearrange axes from input order (lambda, fPDR, P, logC, Z) to desired order (lambda, Z, logC, P, fPDR)
    L = np.moveaxis(data['model_spectra'].T * 1e-7 * 2.99792458e8 / w**2, -1, 0)

    # reverse the wavelength axis so that wavelengths are in increasing order
    w = np.flip(w,0)
    L = np.flip(L,0)

    # remove the three overly prominent carbon emission lines from the dust continuum:
    #    [CII] line at 158 micron and [CI] lines at 370 and 609 micron
    for wline in (158e-6, 370e-6, 609e-6):
        i = np.argmin(np.abs(w-wline))
        w = np.delete(w, i, 0)
        L = np.delete(L, i, 0)

    # write stored table
    writeStoredTable(outFilePaths[0],
          ['lambda','Z','logC','P','fPDR'], ['m','1','1','Pa','1'], ['log','log','lin','log','lin'], [w,Z,logC,P,fPDR],
          ['Llambda'], ['W/m'], ['lin'], [L])   # use linear interpolation to make fPDR interpolation work

# -----------------------------------------------------------------

## This function converts data representing the Castelli-Kurucz family of stellar atmosphere models (i.e. SEDs)
# to stored table format. The function expects a single input file path and the corresponding output file path.
# The input file path is actually interpreted as a template where the asterisks will be replaced by the appropriate
# metallicity and temperature codes.
#
# The stored table output file contains specific fluxes in W/m/m2. To obtain the specific luminosity for a particular
# star, multiply this value by the stellar surface, i.e. 4 pi R^2.
#
def convertCastelliKuruczSEDFamily(inFilePaths, outFilePaths):
    # return metallicity value for a given file path
    def toM(path):
        name = path.split('/')[-1]
        return (-0.1 if name[2]=='m' else 0.1)*(float(name[3:5]))
    # return temperature value for a given file path
    def toT(path):
        return float(path.split('/')[-1].split('_')[1].split('.')[0])
    # return logg value for a given column name
    def toG(name):
        return 0.1*float(name[1:])

    # get a list of all the file paths and file names
    paths = sorted(glob.glob(inFilePaths[0]))

    # get the metallicity grid and the temperature grid from the filenames
    M = np.unique([ toM(path) for path in paths ])
    T = np.unique([ toT(path) for path in paths ])

    # get the wavelength grid and the gravity grid from the first file
    data = fits.open(paths[0])[1].data  # the tables are in the first extension
    w = data['WAVELENGTH']
    g = np.unique([ toG(data.columns[i].name) for i in range(1, len(data.columns)) ])

    # allocate the flux array
    F = np.zeros((len(w), len(M), len(T), len(g)))     # indices k, m, t, n

    # read the data from each file
    for path in paths:
        # get the metallicity index and the temperature index from the filename and the grid
        m = M.tolist().index(toM(path))
        t = T.tolist().index(toT(path))

        # loop over each column in the file (skipping the first wavelength column)
        data = fits.open(path)[1].data
        for i in range(1, len(data.columns)):
            # get the gravity index from the column name and the grid
            name = data.columns[i].name
            n = g.tolist().index(toG(name))
            # get the fluxes
            F[:,m,t,n] = data[name]

    # convert units
    w *= 1e-10          # from Angstrom to m
    Z = 0.02 * 10**M    # from [M/H] to fraction
    g = 10**(g-2)       # from log_g (in cm/s2) to g (in m/s2)
    F *= 1e7            # from erg/s/cm2/A to W/m2/m

    # write stored table
    writeStoredTable(outFilePaths[0],
                          ['lambda','Z','Teff','g'], ['m','1','K','m/s2'], ['log','log','log','log'], [w,Z,T,g],
                          ['Flambda'], ['W/m2/m'], ['log'], [F])

# -----------------------------------------------------------------

## This function converts data representing a family of SEDs for single stellar populations generated by the FSPS code,
# parameterized on metallicity and age, to stored table format. For info on FSPS, see https://github.com/cconroy20/fsps
# and Conroy, Gunn, & White (2009, ApJ) and Conroy & Gunn (2010, ApJ).
# The function expects a sequence of input file paths and corresponding output file paths. Each input file path
# is actually interpreted as a template where the string "*" will be replaced by some metallicity identifier.
# The metallicity identifiers in the filenames are used only for locating the files; the actual metallicity value is
# obtained from the file contents, and the files are automatically sorted in order of increasing metallicity.
#
# The stored table output file contains specific luminosities in W/m normalized on an *initial* stellar population
# mass of 1 solar mass, in function of wavelength and parameterized on metallicity and age.
#
def convertFSPSSEDFamily(inFilePaths, outFilePaths):
    for inFilePath, outFilePath in zip(inFilePaths, outFilePaths):

        # get the list of input file paths
        inpaths = np.array(glob.glob(inFilePath))

        # read the wavelength grid and the age grid from one of the files
        with open(inpaths[0]) as infile:
            # skip header comments and read the number of ages and number of wavelengths
            tokens = TokenizedFile(infile)
            tokens.skipHeaderLines()
            Nt = int(tokens.next())
            Nw = int(tokens.next())

            # read the wavelengths (all on a single line)
            w = np.array([ float(tokens.next()) for i in range(Nw) ])

            # read the ages (on alternating lines)
            def readlogt(tokens):
                result = tokens.next()
                tokens.skipLine()
                tokens.skipLine()
                return float(result)
            logt = np.array([ readlogt(tokens) for i in range(Nt) ])

        # read the metallicities from the first line of each file)
        def readlogZ(inpath):
            with open(inpath) as infile:
                tokens = TokenizedFile(infile)
                tokens.next()   # skip hash character
                tokens.next()   # skip label
                return float(tokens.next())
        logZ = np.array([ readlogZ(inpath) for inpath in inpaths ])
        NZ = len(logZ)

        # sort the metallicity values and the corresponding filepaths in increasing order of metallicity
        ind = np.argsort(logZ)
        logZ = logZ[ind]
        inpaths = inpaths[ind]

        # allocate hypercube for the luminosities
        L = np.zeros((Nw,NZ,Nt))    # indices k, m, p

        # read the luminosities from each of the files
        for m in range(NZ):
            with open(inpaths[m]) as infile:
                tokens = TokenizedFile(infile)

                # skip the header and wavelengths
                tokens.skipHeaderLines()
                tokens.skipLine()
                tokens.skipLine()

                # iterate over ages, and read luminosities for each age
                for p in range(Nt):
                    tokens.skipLine()  # skip intermediate dummy line
                    for k in range(Nw):
                        L[k, m, p] = float(tokens.next())

        # write stored table, converting from input units
        # wavelengths (Angstrom), metallicities (log Z/Zsol), ages (log years), luminosities (Lsun/Hz)
        Angstrom = 1e-10
        c = 2.99792458e8
        Lsun = 3.839e26
        Zsun = 0.0190   # see MANUAL.pdf
        w *= Angstrom
        Z = Zsun * 10**logZ
        t = 10**logt
        L = Lsun * c * (L.T/w**2).T
        writeStoredTable(outFilePath,
                              ['lambda','Z','t'], ['m','1','yr'], ['log','log','log'], [w,Z,t],
                              ['Llambda'], ['W/m'], ['log'], [L])

# -----------------------------------------------------------------
