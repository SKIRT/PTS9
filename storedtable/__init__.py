#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
#  Package initialization file
# -----------------------------------------------------------------

## \package pts.storedtable Support for handling SKIRT stored table and stored columns formats
#
# This package contains classes and functions to help convert data used by SKIRT as a built-in
# resource or as an input file from their original format (whatever that may be) to the proprietary
# "SKIRT stored table" or "SKIRT stored columns" file formats.
#

from .conversionspec import createConversionSpecs, ConversionSpecs, ConversionSpec
from .io import listStoredTableInfo, readStoredTable, writeStoredTable, \
                listStoredColumnsInfo, readStoredColumns, writeStoredColumns
from .tokenizedfile import TokenizedFile
