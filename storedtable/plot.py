#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.plot Plot data from a SKIRT stored table file
#
# This module offer functions to create plots from the data in a SKIRT stored table file.
#

# -----------------------------------------------------------------

import logging
import numpy as np
import matplotlib.pyplot as plt
import pts.utils.path
from pts.storedtable.io import readStoredTable

# -----------------------------------------------------------------

## This function creates a plot of a particular 1D curve from the data in a specified SKIRT stored table file.
# In addition to the file path of a SKIRT stored table file, the function accepts the following optional arguments
# to configure the plot:
#  - horAxis: zero-based index of the table axis on the horizontal plot axis (default = 0)
#  - verAxis: zero-based index of the table quantity on the vertical plot axis (default = 0)
#  - axis0, axis1, axis2, axis3, axis4: ordinate value for the table axis with the indicated zero-based index
#    (default = arithmetic or geometric mean of the axis range; ignored for axes not in the table)
#
# Thus, by default, the script plots the first table quantity as a function of the first table axis,
# with half-way values for the other axes, if any.
#
# The table file path and the plot file path are interpreted according to the rules described for
# the pts.utils.path.absolute() function.
# If no plot path is given, the figure is not saved and it is left open so that is displayed in notebooks.
#
def plotStoredTableCurve(tableFilePath, horAxis=0, verAxis=0, *,
                         axis0=None, axis1=None, axis2=None, axis3=None, axis4=None,
                         plotFilePath=None, figsize=(8,6)):

    # load the complete stored table
    table = readStoredTable(tableFilePath)

    # get info on horizontal axis
    horName = table['axisNames'][horAxis]
    horUnit = table['axisUnits'][horAxis]
    horScale = table['axisScales'][horAxis]
    horGrid = table[horName]

    # get info on vertical axis
    verName = table['quantityNames'][verAxis]
    verUnit = table['quantityUnits'][verAxis]
    verScale = table['quantityScales'][verAxis]

    # get the appropriate slice from the values hypercube
    index = []
    for axisName,axisScale,axisValue in zip(table['axisNames'], table['axisScales'], (axis0,axis1,axis2,axis3,axis4)):
        if axisName == horName:
            index.append(Ellipsis)
        else:
            axisGrid = table[axisName]
            if axisValue is None:
                if axisScale == 'log': axisValue = np.sqrt(axisGrid[0] * axisGrid[-1])
                else: axisValue = (axisGrid[0] + axisGrid[-1]) / 2
            index.append((np.abs(axisGrid - axisValue)).argmin())
    verValues = table[verName][tuple(index)]

    # create the plot
    plt.figure(figsize=figsize)
    if horScale == 'log': plt.xscale('log')
    if verScale == 'log': plt.yscale('log')
    plt.plot(horGrid, verValues)
    plt.xlabel("{} ({})".format(horName, horUnit))
    plt.ylabel("{} ({})".format(verName, verUnit))

    # if a filepath is provided, save the figure; otherwise leave it open
    if plotFilePath is not None:
        plotpath = pts.utils.path.absolute(plotFilePath)
        plt.savefig(plotpath, bbox_inches='tight', pad_inches=0.25)
        plt.close()
        logging.info("Created {}".format(plotpath))

# ----------------------------------------------------------------------
