#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.plotvectorcuts Plot planar cuts or projections for a vector quantity from a SKIRT simulation
#
# The function in this module creates plots of the planar cuts or projections for a vector quantity
# produced by one of the relevant probes in a SKIRT simulation.
#

# -----------------------------------------------------------------

import logging
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.lines
import matplotlib.patches
import numpy as np
import pts.simulation as sm
import pts.utils as ut

# -----------------------------------------------------------------

## This function creates plots of the planar cuts or projections for a vector quantity produced by one of the relevant
# probes in a SKIRT simulation. Specifically, the function accepts a single Simulation instance and it assumes that
# the simulation includes one or more instances of one or more of the probe types listed in \em probeTypes,
# with an associated probe form that produces a planar cut (DefaultCutsForm, PlanarCutsForm) or planar
# projection (ParallelProjectionForm, AllSkyProjectionForm). If this is not the case, the function does nothing.
#
# The figure displays an arrow for each bin of Nx by Ny pixels. The bin size can be specified as an argument.
# The orientation and length of the arrow indicate respectively the direction and magnitude of the vector
# projected on the image plane and averaged over the bin. The color of the arrow scales with the vector
# component orthogonal to the image plane, also averaged over the bin. Vectors pointing away from the
# observer and are red-ish and vectors pointing towards the observer and are blue-ish.
#
# By default, the figures are saved in the simulation output directory with a filename that includes the simulation
# prefix, the probe name, and the medium component or type indicator, and has the ".pdf" filename extension.
# This can be overridden with the out* arguments as described for the pts.utils.savePath() function.
# In interactive mode (see the pts.utils.interactive() function), the figures are not saved and are left open
# for display in notebooks.
#
# The function takes the following arguments:
#   - simulation: the Simulation instance for which to plot the vector quantity.
#   - probeTypes: a sequence of strings with probe types (SKIRT class names) for which to plot the vector quantity.
#   - binSize: the number of pixels in each bin, in horizontal and vertical directions; default is 32 x 32 pixels.
#   - outDirPath: string or Pathlib.Path object that overrides the output directory.
#   - outFileName: string or Pathlib.Path object that overrides the output file name.
#   - outFilePath: string or Pathlib.Path object that overrides the output file path.
#   - figSize: the horizontal and vertical size of the output figure in inch; default is 6x6 inch.
#   - interactive: whether to leave figures open (True) or save them to file (False).
#
def plotVectorCuts(simulation, probeTypes, binSize=(32,32), *,
                   outDirPath=None, outFileName=None, outFilePath=None, figSize=(6,6), interactive=None):

    # find the relevant probes
    probes = simulation.probes(probeTypes,
                    ("DefaultCutsForm", "PlanarCutsForm", "ParallelProjectionForm", "AllSkyProjectionForm"))

    # iterate over the output file paths for each probe
    for probe, path in sm.probeOutFilePaths(probes, ".fits"):

        # load vector data cube with shape (nx, ny, 3)
        vs = sm.loadFits(path)

        # load the axes grids
        xgrid, ygrid, dummygrid = sm.getFitsAxes(path)
        xmin = xgrid[0].value
        xmax = xgrid[-1].value
        ymin = ygrid[0].value
        ymax = ygrid[-1].value
        extent = (xmin, xmax, ymin, ymax)

        # determine binning configuration
        binX = binSize[0]
        orLenX = vs.shape[0]
        dropX =  orLenX % binX
        startX = dropX//2
        binY = binSize[1]
        orLenY = vs.shape[1]
        dropY = orLenY % binY
        startY = dropY//2

        # construct arrays with central bin positions in pixel coordinates
        posX = np.arange(startX - 0.5 + binX / 2.0, orLenX - dropX + startX - 0.5, binX)
        posY = np.arange(startY - 0.5 + binY / 2.0, orLenY - dropY + startY - 0.5, binY)

        # perform the actual binning, while splitting in vector components
        vx = np.zeros((len(posX),len(posY)))
        vy = np.zeros((len(posX),len(posY)))
        vz = np.zeros((len(posX),len(posY)))
        for x in range(len(posX)):
            for y in range(len(posY)):
                vx[x,y] = np.mean(vs[startX+binX*x : startX+binX*(x+1), startY+binY*y : startY+binY*(y+1) , 0].value)
                vy[x,y] = np.mean(vs[startX+binX*x : startX+binX*(x+1), startY+binY*y : startY+binY*(y+1) , 1].value)
                vz[x,y] = np.mean(vs[startX+binX*x : startX+binX*(x+1), startY+binY*y : startY+binY*(y+1) , 2].value)

        # start the figure
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=figSize)

        # configure the axes
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect('equal')

        # configure the axes labels
        if "AllSky" in probe.formType():
            xlabel = r"$\phi$"
            ylabel = r"$\theta$"
        elif "_xy" in path.stem:
            xlabel = "x"
            ylabel = "y"
        elif "_xz" in path.stem:
            xlabel = "x"
            ylabel = "z"
        elif "_yz" in path.stem:
            xlabel = "y"
            ylabel = "z"
        else:
            xlabel = "horizontal"
            ylabel = "vertical"
        ax.set_xlabel(xlabel + sm.latexForUnit(xgrid), fontsize='large')
        ax.set_ylabel(ylabel + sm.latexForUnit(ygrid), fontsize='large')

        # determine a characteristic 'large' field strength in the image plane
        vmax = np.percentile(np.sqrt(vx**2 + vy**2), 99.0)
        if vmax==0: vmax=1      # guard against all zeros

        # determine the scaling so that the longest arrows do not to overlap with neighboring arrows
        lengthScale = 2 * vmax * max(float(len(posX))/figSize[0], float(len(posY))/figSize[1])

        # determine the quiver key label derived from the filename's last or one but last segment
        segments = path.stem.split("_")
        symbol = segments[-1]
        if symbol == "xy" or symbol == "xz" or symbol == "yz":
            symbol = segments[-2]
        if len(symbol)>1:
            symbol = r"$\{}$".format(symbol)
        key = "{}={:.3g}{}".format(symbol, vmax, sm.latexForUnit(vs))

        # determine the color scheme for the component orthogonal to image plane
        vzmax = np.abs(vz).max()
        if vzmax==0: vzmax=1      # guard against all zeros
        normalizer = matplotlib.colors.Normalize(-vzmax, vzmax)

        # plot the vector field (scale positions to data coordinates)
        X,Y = np.meshgrid(xmin + posX * (xmax - xmin) / orLenX, ymin + posY * (ymax - ymin) / orLenY, indexing='ij')
        quiverPlot = ax.quiver(X,Y, vx, vy, vz, cmap='jet', norm=normalizer, pivot='middle', units='inches',
                               angles='xy', scale=lengthScale, scale_units='inches',
                               width=0.015, headwidth=2.5, headlength=2, headaxislength=2, minlength=0.8)
        ax.quiverkey(quiverPlot, 0.75, -0.08*abs((xmax-xmin)/(ymax-ymin)), vmax, key, coordinates='axes', labelpos='E')

        # if not in interactive mode, save the figure; otherwise leave it open
        if not ut.interactive(interactive):
            saveFilePath = ut.savePath(path.with_suffix(".pdf"), (".pdf", ".png"),
                                       outDirPath=outDirPath, outFileName=outFileName, outFilePath=outFilePath)
            plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
            plt.close()
            logging.info("Created {}".format(saveFilePath))

# -----------------------------------------------------------------
