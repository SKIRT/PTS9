#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.convert_band Functions for writing filter transmission curves to stored table format
#
# The functions in this module can load PTS built-in or user-provided broadband filters
# and output their transmission curves to SKIRT stored table format.

# -----------------------------------------------------------------

import glob
import pathlib
import xml.etree.ElementTree
import numpy as np
import astropy.units as u
import pts.band as bnd
from .io import writeStoredTable

# -----------------------------------------------------------------

## This function outputs the transmission curves for all PTS built-in broadbands to stored table files.
# The function ignores the input file paths and it expects a single output file path. The latter is interpreted
# as a template where the asterisk will be replaced by the appropriate PTS band name.
def writeBroadBands(inFilePaths, outFilePaths):

    # loop over all builtin bands
    for name in bnd.builtinBandNames():
        band = bnd.BroadBand(name)

        # construct the output file path based on the band ID
        outFilePath = outFilePaths[0].replace("*",name)

        # get the transmission curve
        w, T = band.transmissionCurve()

        # write stored table
        writeStoredTable(outFilePath, ['lambda'], ['m'], ['lin'], [w.to_value(u.m)],
                                           ['T'], ['1'], ['lin'], [T.to_value(u.m**(-1))])

# -----------------------------------------------------------------

## This function converts a set of filter response curves given in SVO XML format to stored table files.
# The function expects a single input file path template with wildcards expanding to the set of input file names,
# and a single output file path template in which the asterisk will be replaced by the appropriate band name.
# This band name is derived from the input file name by removing the filename extension, converting to uppercase,
# and replacing any dots by underscores.
def convertSvoResponseCurve(inFilePaths, outFilePaths):

    # loop over the list of input files
    for inFilePath in glob.glob(inFilePaths[0]):

        # load the XML tree
        with open(inFilePath, 'r') as bandfile: bandtree = xml.etree.ElementTree.parse(bandfile)

        # get the response curve, converting wavelengths from Angstrom to meter
        values = np.array([ float(el.text) for el in bandtree.findall(".//RESOURCE/TABLE/DATA/TABLEDATA[1]/TR/TD") ])
        if len(values) < 4: raise ValueError("Transmission table not found in '{}'".format(inFilePath))
        wavelengths, transmissions = np.reshape(values, (-1, 2)).T
        wavelengths *= 1e-10

        # convert from response curve to transmission curve, with arbitrary scaling
        transmissions *= wavelengths

        # normalize the transmission curve to unity
        transmissions /= np.trapz(x=wavelengths, y=transmissions)

        # construct the output file path based on the input file name
        name = pathlib.Path(inFilePath).stem.upper().replace(".", "_")
        outFilePath = outFilePaths[0].replace("*",name)

        # write stored table
        writeStoredTable(outFilePath, ['lambda'], ['m'], ['lin'], [wavelengths],
                                      ['T'], ['1'], ['lin'], [transmissions])

# -----------------------------------------------------------------
