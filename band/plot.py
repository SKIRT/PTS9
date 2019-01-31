#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.band.plot Plot built-in broadbands in a given wavelength range or of a given family
#
# This module offers a function to create a plot of the transmission curves for all built-in broadbands
# with a pivot wavelength in a given wavelength range, and/or with a band name that contains specified segments.
#

# -----------------------------------------------------------------

import logging
import matplotlib.pyplot as plt
import pts.utils.path
import pts.band.broadband as bb

# -----------------------------------------------------------------

## This function creates a plot of the transmission curves for all built-in broadbands that satisfy the
# all of the specified selection criteria:
#  - \em minWavelength (float): if specified, the pivot wavelength must exceed this value
#  - \em maxWavelength (float): if specified, the pivot wavelength must be lower than this value
#  - \em nameSegments (string or iterable of strings): if specified, the band name must contain
#    at least one of these segments
#
# The plot file path is interpreted according to the rules described for the pts.utils.path.absolute() function.
# If no plot path is given, the figure is not saved and it is left open so that is displayed in notebooks.
#
def plotBuiltinBands(plotFilePath=None, figsize=(20,6), minWavelength=1e-6, maxWavelength=1e6, nameSegments=None):

    # load all bands that satisfy the specified criteria
    bands = [ bb.BroadBand(name) for name in bb.builtinBandNames() ]
    bands = [ band for band in bands if minWavelength <= band.pivotWavelength() <= maxWavelength ]
    if nameSegments is not None:
        if isinstance(nameSegments, str): nameSegments = [ nameSegments ]
        bands = [ band for band in bands if any([s.lower() in band.name().lower().split("_") for s in nameSegments]) ]

    # sort the remaining bands on pivot wavelength
    bands = sorted(bands, key=bb.BroadBand.pivotWavelength)
    logging.info("Plotting {} built-in bands...".format(len(bands)))

    # setup the figure
    figure = plt.figure(figsize=figsize)
    colors = ('r','g','b','c','m','y')

    # loop over bands
    labelpos = 0.25
    colorindex = 0
    for band in bands:
        wavelengths, transmissions = band.transmissionCurve()
        transmissions /= transmissions.max()  # normalize to a maximum of 1
        plt.plot(wavelengths, transmissions, color=colors[colorindex])
        labelpos += 0.05
        if labelpos > 0.69: labelpos = 0.25
        plt.text(band.pivotWavelength(), labelpos, band.name(), horizontalalignment='center', fontsize='x-small',
                 color=colors[colorindex], backgroundcolor='w')
        colorindex = (colorindex + 1) % len(colors)

    # set axis details
    plt.xscale('log')
    plt.grid(True, axis='y')

    # add axis labels and a legend
    plt.xlabel(r"$\lambda\,(\mu \mathrm{m})$", fontsize='large')
    plt.ylabel("Transmission", fontsize='large')

    # if a filepath is provided, save the figure; otherwise leave it open
    if plotFilePath is not None:
        plotpath = pts.utils.path.absolute(plotFilePath)
        plt.savefig(plotpath, bbox_inches='tight', pad_inches=0.25)
        plt.close()
        logging.info("Created {}".format(plotpath))

# ----------------------------------------------------------------------
