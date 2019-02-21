#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.plotdensitycuts Plot planar cuts through the medium density in a SKIRT simulation
#
# This module offers functions to plot planar cuts through the medium density in a SKIRT simulation,
# allowing a visual comparison between the theoretical and spatially gridded medium density distribution.
#

# -----------------------------------------------------------------

import logging
import matplotlib.colors
import matplotlib.pyplot as plt
import pts.simulation as sm
import pts.utils as ut

# -----------------------------------------------------------------


## This function creates a plot of media density cuts generated by the DefaultMediaDensityCutsProbe in a SKIRT
# simulation, allowing a visual comparison between the theoretical and spatially gridded medium density distribution.
#
# The function accepts a single Simulation instance and it assumes that the simulation includes exactly one
# DefaultMediaDensityCutsProbe. If this is not the case, the function does nothing.
#
# The plot file path is interpreted as described for the pts.utils.absPath() function.
# If no plot path is given, the figure is not saved and it is left open so that is displayed in notebooks.
#
def plotDefaultMediaDensityCuts(simulation, medium="dust", cut="xy", decades=5,
                                *, plotFilePath=None, figSize=(18, 6)):

    # find the relevant probe
    cutProbes = [ probe for probe in simulation.probes() if probe.type() == "DefaultMediaDensityCutsProbe" ]
    if len(cutProbes) != 1:
        return
    cutProbe = cutProbes[0]

    # load the theoretical and gridded cuts for the requested medium and cut plane
    tPaths = cutProbe.outFilePaths("{}_t_{}.fits".format(medium,cut))
    gPaths = cutProbe.outFilePaths("{}_g_{}.fits".format(medium,cut))
    if len(tPaths) != 1 or len(gPaths) != 1:
        return
    tFrame = sm.loadFits(tPaths[0])
    gFrame = sm.loadFits(gPaths[0])

    # determine the range of the x and y axes
    xgrid, ygrid = sm.getFitsAxes(tPaths[0])

    # determine the range of density values to display and clip the data arrays
    vmax = max(tFrame.max(), gFrame.max())
    vmin = vmax / 10**decades
    tFrame[ tFrame<vmin ] = vmin
    gFrame[ gFrame<vmin ] = vmin

    # setup the figure
    fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=figSize)

    # plot the cuts and a color bar
    normalizer = matplotlib.colors.LogNorm(vmin.value, vmax.value)
    extent = (xgrid[0].value, xgrid[-1].value, ygrid[0].value, ygrid[-1].value)
    im = ax1.imshow(tFrame.value, norm=normalizer, cmap='gnuplot', extent=extent,
                    aspect='equal', interpolation='bicubic', origin='lower')
    fig.colorbar(im, ax=(ax1,ax2)).ax.set_ylabel("density" + sm.latexForUnit(tFrame.unit), fontsize='large')
    ax2.imshow(gFrame.value, norm=normalizer, cmap='gnuplot', extent=extent,
               aspect='equal', interpolation='bicubic', origin='lower')

    # set axis details
    ax1.set_xlim(xgrid[0].value, xgrid[-1].value)
    ax1.set_ylim(ygrid[0].value, ygrid[-1].value)
    ax2.set_xlim(xgrid[0].value, xgrid[-1].value)
    ax2.set_ylim(ygrid[0].value, ygrid[-1].value)
    ax1.set_xlabel(cut[0] + sm.latexForUnit(xgrid.unit), fontsize='large')
    ax1.set_ylabel(cut[-1] + sm.latexForUnit(ygrid.unit), fontsize='large')
    ax2.set_xlabel(cut[0] + sm.latexForUnit(xgrid.unit), fontsize='large')
    ax2.set_ylabel(cut[-1] + sm.latexForUnit(ygrid.unit), fontsize='large')

    # if a filepath is provided, save the figure; otherwise leave it open
    if plotFilePath is not None:
        plotpath = ut.absPath(plotFilePath)
        plt.savefig(plotpath, bbox_inches='tight', pad_inches=0.25)
        plt.close()
        logging.info("Created {} density {} cut plot {}".format(medium, cut, plotpath))

# ----------------------------------------------------------------------
