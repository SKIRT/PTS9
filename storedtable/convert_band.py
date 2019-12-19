#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.convert_band Functions for writing filter transmission curves to stored table format
#
# The functions in this module can load PTS built-in broadband filters
# and output their transmission curves to SKIRT stored table format.

# -----------------------------------------------------------------

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
