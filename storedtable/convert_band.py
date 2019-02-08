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

import pts.band.broadband as bb
import pts.storedtable.io

# -----------------------------------------------------------------

## This function outputs the transmission curves for all PTS built-in broadbands to stored table files.
# The function ignores the input file paths and it expects a single output file path. The latter is interpreted
# as a template where the asterisk will be replaced by the appropriate PTS band name.
def writeBroadBands(inFilePaths, outFilePaths):

    # loop over all builtin bands
    for name in bb.builtinBandNames():
        band = bb.BroadBand(name)

        # construct the output file path based on the band ID
        outFilePath = outFilePaths[0].replace("*",name)

        # get the transmission curve and convert from micron to m
        w, T = band.transmissionCurve()
        w *= 1e-6
        T *= 1e6

        # write stored table
        pts.storedtable.io.writeStoredTable(outFilePath, ['lambda'], ['m'], ['lin'], [w], ['T'], ['1'], ['lin'], [T])

# -----------------------------------------------------------------
