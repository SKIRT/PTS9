#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.do.construct_skirt_resources Construct all or part of the SKIRT resources from original data
#
# This script constructs all or part of the SKIRT resources (in SKIRT stored table format) from the original data
# residing in the \c SKIRT/Resources9/OriginalData directory hierarchy, according to the instructions embedded
# in that same hierarchy.
# See the pts.storedtable.conversionspec module for more information.
#
# Provide the name of a (possibly nested) subdirectory as the single argument to restrict execution to that
# portion of the directory tree, or the string "." to construct all resources in the directory hierarchy.
#

# -----------------------------------------------------------------

def do( subDirectory : (str,"name of subdirectory to process or . for all"),
        ) -> "construct SKIRT resources from original data":

    import pts.storedtable.conversionspec
    import pts.utils.path as pp

    # get the paths to the top-level SKIRT resources input and output directories
    parentPath = pp.projectParent()
    inputPath = parentPath / "Resources9" / "OriginalData"
    outputPath = parentPath / "Resources9" / "StoredTables"

    # create the conversion suite
    specs = pts.storedtable.conversionspec.createConversionSpecs(str(inputPath), str(outputPath),
                                                                 subDirectory if subDirectory!="." else None)

    # perform the conversion(s)
    specs.perform()

# -----------------------------------------------------------------
