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
import pathlib

# -----------------------------------------------------------------

## This function expects a sequence of input file paths and corresponding output file paths.
# If there are multiple input paths and a single output path, the function uses the output path
# as a directory path and preserve the input filenames. Otherwise it simply copies each of the
# input files to corresponding output file path.
# In other words, the function has the same signature as the functions for converting some type of data
# into SKIRT stored table format, however, it expects its input files to already have the correct format.
#
def copyWithoutConversion(inFilePaths, outFilePaths):
    # if there are multiple input paths and a single output path,
    # use the output path as a directory path and preserve the input filenames
    if len(inFilePaths) > 1 and len(outFilePaths) == 1:
        outDirPath = pathlib.Path(outFilePaths[0])
        for inFilePath in inFilePaths:
            outFilePath = outDirPath / pathlib.Path(inFilePath).name
            shutil.copyfile(inFilePath, outFilePath)

    # otherwise use the full input and output paths, which allows renaming the file
    else:
        for inFilePath, outFilePath in zip(inFilePaths, outFilePaths):
            shutil.copyfile(inFilePath, outFilePath)

# -----------------------------------------------------------------
