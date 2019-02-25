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
# By default, the figure is saved in the output directory of the first instrument in the list,
# using a name starting with the corresponding simulation prefix, possibly including the instrument name,
# and ending with ".pdf". This can be overridden with the out* arguments as described for the
# pts.utils.savePath() function. In interactive mode (see the pts.utils.interactive() function),
# the figure is not saved and it is left open so that is displayed in notebooks.
#
def plotSeds(simulation, minWavelength=None, maxWavelength=None, decades=None, *,
             outDirPath=None, outFileName=None, outFilePath=None, figSize=(8, 6), interactive=None):

    # get the (instrument, output file path) tuples
    instr_paths = sm.instrumentOutFilePaths(simulation, "sed.dat")
    if len(instr_paths) < 1:
        return

    # setup the figure
    plt.figure(figsize=figSize)

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
    if decades is not None:
        plt.ylim(fluxMax.value * 10 ** (-decades), fluxMax.value * 10 ** 0.2)
    plt.xlabel(r"$\lambda$" + sm.latexForUnit(waveUnit), fontsize='large')
    plt.ylabel(sm.latexForSpectralFlux(fluxUnit) + sm.latexForUnit(fluxUnit), fontsize='large')
    plt.legend(loc='best')

    # if not in interactive mode, save the figure; otherwise leave it open
    if not ut.interactive(interactive):
        # use the first instrument output path; if there are multiple instruments, remove the instrument name
        defSaveFilePath = instr_paths[0][1]
        if len(instr_paths) > 1:
            defSaveFilePath = defSaveFilePath.with_name(instr_paths[0][0].prefix() + "_sed.pdf")
        saveFilePath = ut.savePath(defSaveFilePath, (".pdf",".png"),
                                   outDirPath=outDirPath, outFileName=outFileName, outFilePath=outFilePath)
        plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
        plt.close()
        logging.info("Created {}".format(saveFilePath))

# ----------------------------------------------------------------------

## This function creates a plot of the behavior of the sources in a SKIRT simulations. Specifically, it plots
# the luminosity of and the number of photon packets launched by each source as a function of wavelength.
# Both curves are normalized so that their maximum value is equal to one. This allows a user to evaluate
# whether the photon packet distribution (as a function of wavelength) is appropriate for the source spectra.
#
# The function accepts a single Simulation instance and it assumes that the simulation includes exactly one
# LuminosityProbe and exactly one LaunchedPacketsProbe. If this is not the case, the function does nothing.
#
# By default, the figure is saved in the simulation output directory, using a name starting with the simulation prefix
# and ending with "sources.pdf". This can be overridden with the out* arguments as described for the
# pts.utils.savePath() function. In interactive mode (see the pts.utils.interactive() function),
# the figure is not saved and it is left open so that is displayed in notebooks.
#
def plotSources(simulation, minWavelength=None, maxWavelength=None, decades=None, *,
                outDirPath=None, outFileName=None, outFilePath=None, figSize=(8, 6), interactive=None):

    # find the required probes
    probes = simulation.probes()
    lumiProbes = [ probe for probe in probes if probe.type() == "LuminosityProbe" ]
    packProbes = [ probe for probe in probes if probe.type() == "LaunchedPacketsProbe" ]
    if len(lumiProbes) != 1 or len(packProbes) != 1:
        return
    lumiProbe = lumiProbes[0]
    packProbe = packProbes[0]

    # load the luminosities
    lumiFilePath = lumiProbe.outFilePaths("luminosities.dat")[0]
    descriptions = sm.getColumnDescriptions(lumiFilePath)
    columns = list(range(len(descriptions)))
    del columns[1]  # remove the "specific luminosity column"
    lumiWave, lumiTot, *lumiFracs = sm.loadColumns(lumiFilePath, columns)

    # load the number of launched packets
    packFilePath = packProbe.outFilePaths("launchedpackets.dat")[0]
    descriptions = sm.getColumnDescriptions(packFilePath)
    columns = list(range(len(descriptions)))
    packWave, packTot, *packSplits = sm.loadColumns(packFilePath, columns)
    packWave <<= lumiWave.unit

    # setup the figure
    plt.figure(figsize=figSize)
    label = "{} ".format(simulation.prefix())

    # plot the total
    lumiMax = lumiTot.max()
    packMax = packTot.max()
    plt.plot(lumiWave.value, lumiTot/lumiMax, color='k', ls='solid', label=label + "total")
    plt.plot(packWave.value, packTot/packMax, color='k', ls='dashed')

    # loop over all sources
    colors = ('r', 'g', 'b', 'c', 'm', 'y')
    colorindex = 0
    sourceindex = 1
    for lumiFrac, packSplit in zip(lumiFracs, packSplits):
        plt.plot(lumiWave.value, lumiFrac*lumiTot/lumiMax, color=colors[colorindex], ls='solid',
                 label=label + str(sourceindex))
        plt.plot(packWave.value, packSplit/packMax, color=colors[colorindex], ls='dashed')
        # advance color and source index
        colorindex = (colorindex + 1) % len(colors)
        sourceindex += 1

    # set axis details and add a legend
    plt.xscale('log')
    plt.yscale('log')
    if minWavelength is not None:
        plt.xlim((minWavelength << lumiWave.unit).value, None)
    if maxWavelength is not None:
        plt.xlim(None, (maxWavelength << lumiWave.unit).value)
    if decades is not None:
        plt.ylim(10 ** (-decades), 10 ** 0.2)
    plt.xlabel(r"$\lambda$" + sm.latexForUnit(lumiWave.unit), fontsize='large')
    plt.ylabel(r"Normalized $L$ and $N_\mathrm{pp}$", fontsize='large')
    plt.legend(loc='best')

    # if not in interactive mode, save the figure; otherwise leave it open
    if not ut.interactive(interactive):
        saveFilePath = ut.savePath(simulation.outFilePath("sources.pdf"), (".pdf",".png"),
                                   outDirPath=outDirPath, outFileName=outFileName, outFilePath=outFilePath)
        plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
        plt.close()
        logging.info("Created {}".format(saveFilePath))

# ----------------------------------------------------------------------
