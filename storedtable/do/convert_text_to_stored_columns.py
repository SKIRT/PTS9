#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.do.convert_text_to_stored_columns Convert column text file to SKIRT stored columns format
#
# This script converts a column text file to a SKIRT stored columns file representing the same data.
#
# Provide the name or (relative or absolute) file path of the text file to be converted as the single
# positional argument.
# The output file is placed next to the input file with the ".scol" filename extension.
#

# -----------------------------------------------------------------

def do( filepath : (str,"name or path of the column text file to be converted"),
        names : (str,"white-space-separated list of column names"),
        units : (str,"white-space-separated list of unit strings"),
        ) -> "convert column text file to SKIRT stored columns format":

    import numpy as np
    import pts.storedtable as stab
    import pts.utils as ut

    # get the file paths
    inpath = ut.absPath(filepath)
    outpath = inpath.with_suffix('.scol')

    # load the text file contents
    values = np.loadtxt(inpath, unpack=True)

    # construct list of names and unit string, allowing commas as separators
    names = names.replace(",", " ").split()
    units = units.replace(",", " ").split()

    # save the binary file
    stab.writeStoredColumns(outpath, names, units, values)

# -----------------------------------------------------------------
