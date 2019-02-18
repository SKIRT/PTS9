#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.simulation.text Handling SKIRT text input/output files
#
# This module offers functions for reading information from SKIRT text output files (including column text files,
# log file and convergence reports) and for writing SKIRT text input files (column text files).
#

# -----------------------------------------------------------------

import numpy as np
import pts.utils as ut
from .units import unit as smunit

# -----------------------------------------------------------------

## This function extracts a numeric value from a text file such as a SKIRT log file or a file produced by
# SKIRT's SpatialGridConvergenceProbe class. The value is returned as an astropy quantity with appropriate units,
# assuming that a valid unit string has been found in the text file.
#
# The file path is interpreted as described for the pts.utils.absPath() function.
# The text line containing the value is located by the \em trigger (a text string that must occur on a line
# before the one containing the value) and the \em header (a text string that must occur on the line containing
# the field). In fact, the trigger can consist of multiple subtriggers, separated by forward slashes, that
# must be triggered in sequence before the function starts looking for the header.
# If no conforming line is found, the function raises an error.
#
# If the selected line has a section between parentheses at the end, this section (including parentheses) is removed.
# Then the line is split in segments on white space. If the last segment is a valid representation of a floating point
# number, this number is returned as a dimensionless quantity. Otherwise, the last two segments represent the value and
# the unit string, respectively. If there are problems converting these representations, the function raises an error.
#
# For example, the following call will return the total gridded dust mass from a spatial grid convergence file:
#
#     getQuantityFromFile("xxx_convergence.dat", "Dust/Total mass", "Gridded")
#
# and the following call will return the total metallic mass for the first particle medium from a log file:
#
#     getQuantityFromFile("xxx_log.txt", "ParticleMedium", "Total metallic mass")
#
def getQuantityFromFile(path, trigger, header):
    path = ut.absPath(path)

    triggers = trigger.split("/")
    triggered = 0
    with open(path) as infile:
        for line in infile:
            # select the line
            if triggered < len(triggers) and triggers[triggered] in line: triggered += 1
            elif triggered==len(triggers) and header in line:
                # remove section between parentheses
                if line.strip().endswith(')'):
                    line = line[:line.rfind('(')]
                segments = line.split()
                if len(segments) >= 2:
                    try:
                        return float(segments[-1]) << smunit("")
                    except ValueError: pass
                    return float(segments[-2]) << smunit(segments[-1])

    raise ValueError("Quantity '{}' not found in text file".format(header))

# -----------------------------------------------------------------

## This function reads the header of a SKIRT column text file, and returns a list of the column descriptions,
# in order of occurrence. If the specified file does not have a column header in SKIRT format, the function
# returns an empty list. The file path is interpreted as described for the pts.utils.absPath() function.
def getColumnDescriptions(path):
    path = ut.absPath(path)

    # parse the header
    descriptions = []
    with open(path) as infile:
        for line in infile:
            # skip empty lines
            line = line.strip()
            if len(line) > 0:
                # handle end-of-header
                if not line.startswith("#"): break;
                # remove hash character and skip non-column header lines
                line = line[1:].strip()
                if line.lower().startswith("column") and ":" in line:
                    # extract the description
                    colon = line.find(":")
                    if line.endswith(')'):
                        left = line.rfind("(")
                        description = line[colon+1:left].strip()
                    else:
                        description = line[colon+1:].strip()
                    # add the description to the list
                    descriptions.append(description)
    return descriptions

## This function loads data from a SKIRT column text file, using the metadata in the file's header to assign
# appropriate units to each column, and optionally, to locate columns by name.
#
# \note The specified file \em must have a column header in SKIRT format, because the information in the header
# is used to determine the appropriate units for each column. If the file does not have a conforming column header,
# use the numpy.loadtxt() function directly.
#
# The file path is interpreted as described for the pts.utils.absPath() function.
#
# The columns to be loaded can be specified in one of the following ways:
#  - A sequence (list or tuple) of integers specifying zero-based column indices.
#  - A single string including a comma-separated list of one-based column indices and/or description fragments;
#    each fragment must be contained in exactly one column description in the file header.
#  - None or omitted: causes all columns to be loaded.
#
# The function always returns a list of astropy quantities, even if only one column has been loaded. The units for
# each item are derived from the metadata in the file header. If the file contains a single row, each item in the
# returned list is a scalar astropy quantity. Otherwise each item in the list is an astropy quantity array, and
# all arrays have the same length.
#
def loadColumns(path, columns=None):
    path = ut.absPath(path)

    # parse the header
    header = []
    with open(path) as infile:
        for line in infile:
            # skip empty lines
            line = line.strip()
            if len(line) > 0:
                # handle end-of-header
                if not line.startswith("#"): break;
                # remove hash character and skip non-column header lines
                line = line[1:].strip()
                if line.lower().startswith("column") and ":" in line:
                    # extract the description and the unit
                    colon = line.find(":")
                    if line.endswith(')'):
                        left = line.rfind("(")
                        description = line[colon+1:left].strip()
                        unit = line[left+1:-1].strip()
                    else:
                        description = line[colon+1:].strip()
                        unit = ""
                    # add the column metadata to the internal list
                    header.append( (description, unit) )

    # construct a list of zero-based indices for columns to be loaded
    if columns is None:
        # None --> all columns
        usecols = list(range(len(header)))
    elif isinstance(columns, str):
        # single string --> comma-separated list of one-based column indices and/or description fragments
        usecols = []
        for col in columns.split(","):
            # string representing an integer -> one-based index
            if _representsInteger(col):
                col = int(col)
                if col<1 or col>len(header):
                    raise ValueError("one-based column index is out of range: {}".format(col))
                usecols.append(col-1)
            # other string --> match column description
            else:
                usecols.append(_indexForDescriptionInHeader(col, header))
    else:
        # sequence --> zero-based column indices
        usecols = []
        for col in columns:
            col = int(col)
            if col<0 or col>=len(header):
                raise ValueError("zero-based column index is out of range: {}".format(col))
            usecols.append(col)

    # load the data and assign units
    data = np.loadtxt(path, usecols=usecols, ndmin=2, unpack=True)
    return [coldata << smunit(header[colindex][1]) for colindex, coldata in zip(usecols, data)]


## This helper function returns true is the specified value can be converted to an integer, and False otherwise.
def _representsInteger(value):
    try:
        int(value)
        return True
    except ValueError:
        return False

## This helper function searches the header metadata for a column with a description containing the specified string.
# If there is exactly one such column, the function returns its zero-based index. Otherwise, it raises an error.
def _indexForDescriptionInHeader(spec, header):
    count = 0
    result = -1
    for index,(description,unit) in enumerate(header):
        if spec.strip().lower() in description.lower():
            count += 1
            result = index

    if count == 0:
        raise ValueError("Column specification does not match header descriptions: {}".format(spec))
    if count > 1:
        raise ValueError("Column specification matches multiple header descriptions: {}".format(spec))
    return result

# -----------------------------------------------------------------

## This function writes data to a SKIRT column text file, generating the metadata in the file's header from the
# information provided to the function, and automatically converting the data to the specified units.
#
# The function expects the following arguments:
#  - \em path: the file path to the output file, interpreted as described for the pts.utils.absPath() function.
#  - \em quantities: a sequence of astropy quantity arrays of the same length, one for each column
#    (if the output has only a single row, each item in sequence is a scalar scalar astropy quantity).
#  - \em units: a string containing a comma-separated list of valid SKIRT unit strings.
#  - \em descriptions: a string containing a comma-separated list of quantity descriptions.
#
# The arguments must specify the same number of quantities, units, and descriptions. The units and descriptions
# will appear in the header of the output file. The quantities are converted to the requested units before being
# written to the file.
#
def saveColumns(path, quantities, units, descriptions):

    # split descriptions and units into segments
    descriptions = [ s.strip() for s in descriptions.split(",") ]
    units = [ s.strip() for s in units.split(",") ]

    # verify length of input sequences
    if len(quantities) != len(units) or len(quantities) != len(descriptions):
        raise ValueError("Number of units or descriptions does not match number of quantities")

    # convert the quantities to the requested units
    quantities = [quantity.to(smunit(unit)) for quantity, unit in zip(quantities, units)]

    # open the file
    path = ut.absPath(path)
    with open(path, 'wt') as outfile:

        # write the header
        for col, (description, unit) in enumerate(zip(descriptions, units)):
            outfile.write("# column {}: {} ({})\n".format(col+1, description, unit))

        # write the data
        np.savetxt(outfile, np.stack(quantities).T, fmt="%1.9e")

# -----------------------------------------------------------------
