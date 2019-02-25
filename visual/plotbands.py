#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.plotbands Plot built-in broadbands in a given wavelength range or of a given family
#
# This module offers a function to create a plot of the transmission curves for all built-in broadbands
# with a pivot wavelength in a given wavelength range, and/or with a band name that contains specified segments.
#

# -----------------------------------------------------------------

import astropy.units as u
import logging
import matplotlib.pyplot as plt
import pts.band as bnd
import pts.simulation as sm
import pts.utils as ut

# -----------------------------------------------------------------

## This function creates a plot of the transmission curves for all built-in broadbands that satisfy the
# all of the specified selection criteria:
#  - \em minWavelength (astropy quantity): if specified, the pivot wavelength must exceed this value
#  - \em maxWavelength (astropy quantity): if specified, the pivot wavelength must be lower than this value
#  - \em nameSegments (string or iterable of strings): if specified, the band name must contain
#    at least one of these segments
#
# By default, the figure is saved as FigBuiltinBands.pdf in the current directory. This can be overridden with the
# out* arguments as described for the pts.utils.savePath() function. In interactive mode (see the
# pts.utils.interactive() function), the figure is not saved and it is left open so that is displayed in notebooks.
#
def plotBuiltinBands(minWavelength=1e-6*u.micron, maxWavelength=1e6*u.micron, nameSegments=None, *,
                     outDirPath=None, outFileName=None, outFilePath=None, figSize=(20, 6), interactive=None):

    # load all bands that satisfy the specified criteria
    bands = [ bnd.BroadBand(name) for name in bnd.builtinBandNames() ]
    bands = [ band for band in bands if minWavelength <= band.pivotWavelength() <= maxWavelength ]
    if nameSegments is not None:
        if isinstance(nameSegments, str): nameSegments = [ nameSegments ]
        bands = [ band for band in bands if any([s.lower() in band.name().lower().split("_") for s in nameSegments]) ]

    # sort the remaining bands on pivot wavelength
    bands = sorted(bands, key=bnd.BroadBand.pivotWavelength)
    logging.info("Plotting {} built-in bands...".format(len(bands)))

    # setup the figure
    plt.figure(figsize=figSize)
    colors = ('r','g','b','c','m','y')

    # loop over bands
    labelpos = 0.25
    colorindex = 0
    for band in bands:
        wavelengths, transmissions = band.transmissionCurve()
        wavelengths <<= u.micron              # convert to micron
        transmissions /= transmissions.max()  # normalize to a maximum of 1
        plt.plot(wavelengths.value, transmissions.value, color=colors[colorindex])
        labelpos += 0.05
        if labelpos > 0.69: labelpos = 0.25
        plt.text(band.pivotWavelength().to_value(wavelengths.unit), labelpos, band.name(),
                 horizontalalignment='center', fontsize='x-small', color=colors[colorindex], backgroundcolor='w')
        colorindex = (colorindex + 1) % len(colors)

    # set axis details
    plt.xscale('log')
    plt.grid(True, axis='y')

    # add axis labels and a legend
    plt.xlabel(r"$\lambda$" + sm.latexForUnit(wavelengths), fontsize='large')
    plt.ylabel("Transmission", fontsize='large')

    # if not in interactive mode, save the figure; otherwise leave it open
    if not ut.interactive(interactive):
        saveFilePath = ut.savePath("FigBuiltinBands.pdf", (".pdf",".png"),
                                   outDirPath=outDirPath, outFileName=outFileName, outFilePath=outFilePath)
        plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
        plt.close()
        logging.info("Created {}".format(saveFilePath))

# ----------------------------------------------------------------------
