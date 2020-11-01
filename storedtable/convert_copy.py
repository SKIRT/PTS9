#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.convert_copy Function to "convert" stored table files by simply copying them
#
# The function in this module has the same signature as the functions for converting some type of data
# into SKIRT stored table format, however, it expects its input files to already have the correct format.

# -----------------------------------------------------------------

import shutil

# -----------------------------------------------------------------

## This function expects a sequence of input file paths and corresponding output file paths.
# It simply copies each of the input files to corresponding output file path.
# In other words, the function has the same signature as the functions for converting some type of data
# into SKIRT stored table format, however, it expects its input files to already have the correct format.
#
def copyWithoutConversion(inFilePaths, outFilePaths):
    for inFilePath, outFilePath in zip(inFilePaths, outFilePaths):
        shutil.copyfile(inFilePath, outFilePath)

# -----------------------------------------------------------------
