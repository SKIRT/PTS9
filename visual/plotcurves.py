#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.plotcurves Plot SEDs or other SKIRT output as a function of wavelength
#
# This module offers functions to plot spectral distributions or other SKIRT output as a function of wavelength.
# Each function accepts a list of Instrument or Probe instances of a certain type, and plots the corresponding output.
#

# -----------------------------------------------------------------

import logging
import matplotlib.pyplot as plt
import pts.simulation as sm
import pts.utils as ut

# -----------------------------------------------------------------

## This function creates a plot of the SEDs produced by one or more SKIRT simulations. Specifically, the function
# accepts a sequence of Simulation and/or Instrument instances (or a single instance of either of these types),
# and it plots the total flux density for each of the instruments that actually produced an "*_sed.dat" file.
#
# If there is only a single SED, and the corresponding output file includes the individual flux components,
# the function plots these as well.
#
# If there are multiple SEDs, the plot uses the flux flavor and units of the first SED and the data for all
# other SEDs are converted, if needed.
#
# If there is no SED in the list, the function does nothing.
#
# The plot file path is interpreted as described for the pts.utils.absPath() function.
# If no plot path is given, the figure is not saved and it is left open so that is displayed in notebooks.
#
def plotSeds(instruments, minWavelength=None, maxWavelength=None, decadesFlux=None,
             *, plotFilePath=None, figsize=(8, 6)):
    # get the (instrument, output file path) tuples
    instr_paths = sm.instrumentOutFilePaths(instruments, "sed.dat")
    if len(instr_paths) < 1:
        return

    # setup the figure
    plt.figure(figsize=figsize)

    # if there is a single output file, and it has components, plot the components
    if len(instr_paths) == 1 and any(["transparent" in col for col in sm.getColumnDescriptions(instr_paths[0][1])]):
        instrument, filepath = instr_paths[0]
        # load the columns (we assume that all components have the same units)
        wave, tot, tra, dirpri, scapri, dirsec, scasec = sm.loadColumns(filepath, (0, 1, 2, 3, 4, 5, 6))
        waveUnit = wave.unit
        fluxUnit = tot.unit
        fluxMax = max(tot.max(), tra.max())
        # plot the various components
        label = "{} {} ".format(instrument.prefix(), instrument.name())
        plt.plot(wave.value, tot.value, color='k', ls='solid', label=label + "total")
        plt.plot(wave.value, tra.value, color='b', ls='dotted', label=label + "transparent")
        plt.plot(wave.value, (dirpri + scapri).value, color='b', ls='solid', label=label + "primary")
        plt.plot(wave.value, (dirsec + scasec).value, color='r', ls='solid', label=label + "secondary")

    # otherwise loop over all SEDs
    else:
        colors = ('r', 'g', 'b', 'c', 'm', 'y')
        colorindex = 0
        first = True
        for instrument, filepath in instr_paths:
            # load the total flux; first time remember units; thereafter convert units
            wave, flux = sm.loadColumns(filepath, (0, 1))
            if first:
                waveUnit = wave.unit
                fluxUnit = flux.unit
                fluxMax = flux.max()
                first = False
            else:
                wave <<= waveUnit
                flux = sm.convertToFlavor(wave, flux, fluxUnit)
                fluxMax = max(fluxMax, flux.max())
            # plot
            plt.plot(wave.value, flux.value, color=colors[colorindex],
                     label="{} {} total".format(instrument.prefix(), instrument.name()))
            # advance color index
            colorindex = (colorindex + 1) % len(colors)

    # set axis details and add a legend
    plt.xscale('log')
    plt.yscale('log')
    if minWavelength is not None:
        plt.xlim((minWavelength << waveUnit).value, None)
    if maxWavelength is not None:
        plt.xlim(None, (maxWavelength << waveUnit).value)
    if decadesFlux is not None:
        plt.ylim(fluxMax.value * 10 ** (0.2 - decadesFlux), fluxMax.value * 10 ** 0.2)
    plt.xlabel(r"$\lambda$" + sm.latexForUnit(waveUnit), fontsize='large')
    plt.ylabel(sm.latexForSpectralFlux(fluxUnit) + sm.latexForUnit(fluxUnit), fontsize='large')
    plt.legend(loc='best')

    # if a filepath is provided, save the figure; otherwise leave it open
    if plotFilePath is not None:
        plotpath = ut.absPath(plotFilePath)
        plt.savefig(plotpath, bbox_inches='tight', pad_inches=0.25)
        plt.close()
        logging.info("Created SED plot {}".format(plotpath))

# ----------------------------------------------------------------------
