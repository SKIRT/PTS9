#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.io Input/output functions for files in the SKIRT stored table and stored column format
#
# The functions in this module allow reading and writing data and metadata from and to a disk file
# in SKIRT stored table format or in SKIRT stored column format.
#
# SKIRT stored table format
# -------------------------
#
# The SKIRT stored table file format is intended for tabulating one or more physical quantities
# as a function of on one or more independent variables. It has been designed for representing
# built-in SKIRT resources but it can be used for multi-dimensional tables in other contexts.
#
# A SKIRT stored table file contains binary data with a straightforward layout. The file is essentially
# a sequence of 8-byte data items. A data item can have one of three types:
#  - string: 1 to 8 printable and non-whitespace 7-bit ASCII characters, padded with spaces if needed.
#  - unsigned integer: 64-bit integer in little-endian byte order
#  - floating point: 64-bit double (IEEE 754) in little-endian byte order
#
# The overall layout is as follows:
#  - SKIRT name/version tag
#  - Endianness tag
#  - numAxes
#  - axisName (x numAxes)
#  - axisUnit (x numAxes)
#  - axisScale (x numAxes)
#  - [ numPoints  axisPoint (x numPoints) ] (x numAxes)
#  - numQuantities
#  - quantityName (x numQuantities)
#  - quantityUnit (x numQuantities)
#  - quantityScale (x numQuantities)
#  - value (x numQuantities x numPoints1 x ... x numPointsN)
#  - end-of-file tag
#
# The values are ordered so that the quantity values for a particular point are next to each other,
# the first axis index varies most rapidly, and the last axis index varies least rapidly.
#
# SKIRT stored columns format
# ---------------------------
#
# The SKIRT stored columns file format stores a set of data columns and is intended as a faster
# replacement for large regular text column data files, without the benefit of being human readable.
# The format does not support non-leaf rows, so it cannot be used for adaptive mesh data.
#
# A SKIRT stored columns file contains binary data with a straightforward layout. It can contain
# the same type of data items as described for stored tables above. The layout is similar but not
# identical to stored tables:
#
#  - SKIRT name/version tag
#  - Endianness tag
#  - 0  (to differentiate from stored table format, which stores the nr of axes here)
#  - numRows
#  - numColumns
#  - columnName (x numColumns)
#  - columnUnit (x numColumns)
#  - value (x numColumns x numRows)
#  - end-of-file tag
#
# The values are ordered so that the column values for a particular row are next to each other.
#

# -----------------------------------------------------------------

import functools
import logging
import numpy as np
import struct
import pts.simulation as sm
import pts.utils as ut

# -----------------------------------------------------------------

## This function lists relevant metadata about the specified SKIRT stored table file. The printed information
# includes the names, units and ranges for each of the axes, and the names and units for each of the quantities.
# The table file path is interpreted as described for the pts.utils.absPath() function.
def listStoredTableInfo(tableFilePath):
    inpath = ut.absPath(tableFilePath)

    # open the file
    with open(inpath, 'rb') as infile:
        # verify header tags
        if stringFromFile(infile) != "SKIRT X" or intFromFile(infile) != 0x010203040A0BFEFF:
            raise ValueError("File does not have SKIRT stored table format: {}".format(inpath))

        # get the axes metadata and grids
        numAxes = intFromFile(infile)
        axisNames = [ stringFromFile(infile) for i in range(numAxes) ]
        axisUnits = [ stringFromFile(infile) for i in range(numAxes) ]
        axisScales = [ stringFromFile(infile) for i in range(numAxes) ]
        axisGrids = [ arrayFromFile(infile, (intFromFile(infile),)) for i in range(numAxes) ]

        # get the quantities metadata
        numQuantities = intFromFile(infile)
        quantityNames = [ stringFromFile(infile) for i in range(numQuantities) ]
        quantityUnits = [ stringFromFile(infile) for i in range(numQuantities) ]
        quantityScales = [ stringFromFile(infile) for i in range(numQuantities) ]

    # print file path
    logging.info("Stored table file: {}".format(inpath))

    # print axes information
    for i in range(numAxes):
        logging.info("  axis {}: {} ({}) with {} points from {:.3e} to {:.3e} on {} scale"  \
           .format(i, axisNames[i], axisUnits[i], len(axisGrids[i]), axisGrids[i][0], axisGrids[i][-1], axisScales[i]))

    # print quantities information
    for i in range(numQuantities):
        logging.info("  quantity {}: {} ({}) on {} scale"  \
                     .format(i, quantityNames[i], quantityUnits[i], quantityScales[i]))

# -----------------------------------------------------------------

## This function reads the specified SKIRT stored table file and returns a dictionary with the following structure:
#  - axisNames: list of axis names, in order of occurrence
#  - axisUnits: list of corresponding units
#  - axisScales: list of corresponding scales
#  - quantityNames: list of quantity names, in order of occurrence
#  - quantityUnits: list of corresponding units
#  - quantityScales: list of corresponding scales
#  - for each axis name: array with grid points as an astropy quantity with the appropriate unit
#  - for each quantity name: array with values as an astropy quantity with the appropriate unit
# The table file path is interpreted as described for the pts.utils.absPath() function.
#
def readStoredTable(tableFilePath):
    inpath = ut.absPath(tableFilePath)

    # open the file
    with open(inpath, 'rb') as infile:
        # verify header tags
        if stringFromFile(infile) != "SKIRT X" or intFromFile(infile) != 0x010203040A0BFEFF:
            raise ValueError("File does not have SKIRT stored table format: {}".format(inpath))

        # get the axes metadata and grids
        numAxes = intFromFile(infile)
        axisNames = [ stringFromFile(infile) for i in range(numAxes) ]
        axisUnits = [ stringFromFile(infile) for i in range(numAxes) ]
        axisScales = [ stringFromFile(infile) for i in range(numAxes) ]
        axisGrids = [ arrayFromFile(infile, (intFromFile(infile),)) for i in range(numAxes) ]

        # get the quantities metadata
        numQuantities = intFromFile(infile)
        quantityNames = [ stringFromFile(infile) for i in range(numQuantities) ]
        quantityUnits = [ stringFromFile(infile) for i in range(numQuantities) ]
        quantityScales = [ stringFromFile(infile) for i in range(numQuantities) ]

        # get the quantity values
        shapeValues = tuple([ numQuantities ] + [ len(axisGrid) for axisGrid in axisGrids ])
        values = arrayFromFile(infile, shapeValues)

        # verify the trailing tag
        if stringFromFile(infile) != "STABEND":
            raise ValueError("File does not have the proper trailing tag: {}".format(inpath))

    # construct the dictionary that will be returned, adding basic metadata
    d = dict(axisNames=axisNames, axisUnits=axisUnits, axisScales=axisScales,
             quantityNames=quantityNames, quantityUnits=quantityUnits, quantityScales=quantityScales)

    # add axis grids
    for i in range(numAxes):
        d[axisNames[i]] = axisGrids[i] << sm.unit(axisUnits[i])

    # add quantities information
    for i in range(numQuantities):
        d[quantityNames[i]] = values[i] << sm.unit(quantityUnits[i])

    return d

# -----------------------------------------------------------------

## This function writes the specified data and metadata to a disk file in SKIRT stored table format.
#
# \note This low-level function does \em not support astropy quantities; it expects plain numpy arrays with values
#       presented in the proper units (i.e. as declared by the unit strings).
#
# The function expects the following arguments:
#  - outFilePath: the output file's absolute or relative path, including filename and '.stab' extension
#  - axisNames: a sequence of axis name strings
#  - axisUnits: a sequence of axis unit strings (should match internal SKIRT units, i.e. usually SI units)
#  - axisScales: a sequence of axis scale strings ('lin' or 'log') controlling interpolation behavior
#  - axisGrids: a sequence of numpy arrays specifying the grid points for each axis in the correct units
#  - quantityNames: a sequence of quantity name strings
#  - quantityUnits: a sequence of quantity unit strings (should match internal SKIRT units, i.e. usually SI units)
#  - quantityScales: a sequence of quantity scale strings ('lin' or 'log') controlling interpolation behavior
#  - quantityValues: a sequence of numpy arrays specifying the values for each quantity in the correct units;
#                     each array has indices in the same order and with the same range as the specified axes
#
# The sequences must be nonempty and have the same number of items within each group (axis and quantity).
# A name or unit string must contain 1 to 8 printable and non-whitespace 7-bit ASCII characters.
# The output table file path is interpreted as described for the pts.utils.absPath() function.
#
def writeStoredTable(outFilePath, axisNames, axisUnits, axisScales, axisGrids,
                                  quantityNames, quantityUnits, quantityScales, quantityValues):
    outpath = ut.absPath(outFilePath)

    # verify some of the requirements/restrictions on the specified data
    if outpath.suffix != ".stab":
        raise ValueError("Stored table filename extension is not '.stab': {}".format(outpath))
    numAxes = len(axisNames)
    if numAxes<1 or numAxes>9 or numAxes!=len(axisUnits) or numAxes!=len(axisGrids):
        raise ValueError("Mismatch in number of axes")
    shapeValues = tuple([ len(axisGrid) for axisGrid in axisGrids ])
    numQuantities = len(quantityNames)
    if numQuantities<1 or numQuantities>9 or numQuantities!=len(quantityUnits) or numQuantities!=len(quantityValues):
        raise ValueError("Mismatch in number of quantities")
    for valueArray in quantityValues:
        if valueArray.shape != shapeValues:
            raise ValueError("Mismatch in number of values")

    # open the output file
    with open(outpath, 'wb') as out:
        # write the SKIRT and endianness tags
        out.write(b"SKIRT X\n")
        intToFile(out, 0x010203040A0BFEFF)

        # write the axes metadata
        intToFile(out, numAxes)
        for axisName in axisNames: stringToFile(out, axisName)
        for axisUnit in axisUnits: stringToFile(out, axisUnit)
        for axisScale in axisScales: stringToFile(out, axisScale)

        # write the grid for each axis
        for axisGrid in axisGrids:
            intToFile(out, len(axisGrid))
            arrayToFile(out, axisGrid)

        # write the quantities metadata
        intToFile(out, numQuantities)
        for quantityName in quantityNames: stringToFile(out, quantityName)
        for quantityUnit in quantityUnits: stringToFile(out, quantityUnit)
        for quantityScale in quantityScales: stringToFile(out, quantityScale)

        # write the values
        arrayToFile(out, np.stack(quantityValues))

        # write the EOF tag
        out.write(b"STABEND\n")

    # report file creation to user
    logging.info("Created stored table file: {}".format(outpath))

# -----------------------------------------------------------------

## This function lists relevant metadata about the specified SKIRT stored columns file. The printed information
# includes the names and units and ranges for each of the columns and the number of rows.
# The stored columns file path is interpreted as described for the pts.utils.absPath() function.
def listStoredColumnsInfo(columnsFilePath):
    inpath = ut.absPath(columnsFilePath)

    # open the file
    with open(inpath, 'rb') as infile:
        # verify header tags
        if stringFromFile(infile) != "SKIRT X" or intFromFile(infile) != 0x010203040A0BFEFF \
                        or intFromFile(infile) != 0:
            raise ValueError("File does not have SKIRT stored columns format: {}".format(inpath))

        # get the number of columns and rows
        numRows = intFromFile(infile)
        numColumns = intFromFile(infile)

        # get the column metadata
        columnNames = [ stringFromFile(infile) for i in range(numColumns) ]
        columnUnits = [ stringFromFile(infile) for i in range(numColumns) ]

    # print file path
    logging.info("Stored columns file: {}".format(inpath))

    # print columns information
    for i in range(numColumns):
        logging.info("  column {}: {} ({})".format(i+1, columnNames[i], columnUnits[i]))
    logging.info("  nr of rows: {}".format(numRows))

# -----------------------------------------------------------------

## This function reads the specified SKIRT stored columns file and returns a dictionary with the following structure:
#  - columnNames: list of column names, in order of occurrence
#  - columnUnits: list of corresponding units
#  - for each column name: array with values (one per row) as an astropy quantity with the appropriate unit
# The stored columns file path is interpreted as described for the pts.utils.absPath() function.
#
def readStoredColumns(columnsFilePath):
    inpath = ut.absPath(columnsFilePath)

    # open the file
    with open(inpath, 'rb') as infile:
        # verify header tags
        if stringFromFile(infile) != "SKIRT X" or intFromFile(infile) != 0x010203040A0BFEFF \
                        or intFromFile(infile) != 0:
            raise ValueError("File does not have SKIRT stored table format: {}".format(inpath))

        # get the number of columns and rows
        numRows = intFromFile(infile)
        numColumns = intFromFile(infile)

        # get the column metadata
        columnNames = [ stringFromFile(infile) for i in range(numColumns) ]
        columnUnits = [ stringFromFile(infile) for i in range(numColumns) ]

        # get the data values
        values = arrayFromFile(infile, (numColumns, numRows))

        # verify the trailing tag
        if stringFromFile(infile) != "SCOLEND":
            raise ValueError("File does not have the proper trailing tag: {}".format(inpath))

    # construct the dictionary that will be returned, adding basic metadata
    d = dict(columnNames=columnNames, columnUnits=columnUnits)

    # add data values
    for i in range(numColumns):
        d[columnNames[i]] = values[i] << sm.unit(columnUnits[i])

    return d

# -----------------------------------------------------------------

## This function writes the specified data and metadata to a disk file in SKIRT stored columns format.
#
# \note This low-level function does \em not support astropy quantities; it expects plain numpy arrays with values
#       presented in the proper units (i.e. as declared by the unit strings).
#
# The function expects the following arguments:
#  - outFilePath: the output file's absolute or relative path, including filename and '.scol' extension
#  - columnNames: a sequence of column name strings
#  - columnUnits: a sequence of column unit strings (should match supported SKIRT import units)
#  - values: a single 2D numpy array with shape (numColumns, numRows) or
#            a sequence of numColumns 1D numpy arrays with length numRows
#            specifying the values for each column in the correct units
#
# A name or unit string must contain 1 to 8 printable and non-whitespace 7-bit ASCII characters.
# There must be at least one column; the number of rows may be zero.
# The output columns file path is interpreted as described for the pts.utils.absPath() function.
#
def writeStoredColumns(outFilePath, columnNames, columnUnits, values):
    outpath = ut.absPath(outFilePath)

    # assemble all values in single array
    values = np.array(values)
    numColumns, numRows = values.shape

    # verify some of the requirements/restrictions on the specified data
    if outpath.suffix != ".scol":
        raise ValueError("Stored columns filename extension is not '.scol': {}".format(outpath))
    numColumns = len(columnNames)
    if numColumns<1 or numColumns!=len(columnNames) or numColumns!=len(columnUnits):
        raise ValueError("Mismatch in number of columns")

    # open the output file
    with open(outpath, 'wb') as out:
        # write the SKIRT and endianness tags
        out.write(b"SKIRT X\n")
        intToFile(out, 0x010203040A0BFEFF)
        intToFile(out, 0)

        # write the rows and columns metadata
        intToFile(out, numRows)
        intToFile(out, numColumns)
        for columnName in columnNames: stringToFile(out, columnName)
        for columnUnit in columnUnits: stringToFile(out, columnUnit)

        # write the values
        arrayToFile(out, values)

        # write the EOF tag
        out.write(b"SCOLEND\n")

    # report file creation to user
    logging.info("Created stored columns file: {}".format(outpath))

# -----------------------------------------------------------------

## This helper function reads an 8-byte little-endian unsigned integer from the specified (open) file and returns
# the corresponding integer value.
def intFromFile(fd):
    return struct.unpack('<Q', fd.read(8))[0]

## This helper function reads an 8-byte name or unit string from the specified (open) file and returns the
# resulting variable-sized string after stripping trailing white space.
def stringFromFile(fd):
    return str(fd.read(8).strip(),'utf8')

## This helper function reads a sequence of doubles from the specified (open) file to a numpy array with the given
# shape. The doubles should appear in the file in the required order (first index varies most rapidly).
def arrayFromFile(fd, shape):
    count = functools.reduce(lambda x, y: x*y, shape)
    a = np.fromfile(fd, dtype='<f8', count=count)
    return np.reshape(a, tuple(reversed(shape))).T

# -----------------------------------------------------------------

## This helper function writes the specified unsigned integer value to the given (open) file
# as a little-endian 8-byte sequence.
def intToFile(fd, n):
    fd.write(struct.pack('<Q', n))

## This helper function writes the specified name or unit string to the given (open) file as an 8-byte sequence,
# after checking that the string adheres to the proper requirements.
def stringToFile(fd, s):
    s = str(s).strip()
    if len(s)<1 or len(s)>8:
        raise ValueError("Too few or too many characters in name or unit string")
    for c in s:
        if not ( 32 < ord(c) < 127 ):
            raise ValueError("Non-printable or 8-bit ASCII character in name or unit string: " + s)
    fd.write(bytes(s.ljust(8),'utf8'))

## This helper function writes the specified numpy array to the given (open) file as a sequence of doubles
# in the required order (first index varies most rapidly).
def arrayToFile(fd, a):
    a.T.astype(dtype='<f8', order='C', casting='safe', copy=False).tofile(fd)

# -----------------------------------------------------------------
