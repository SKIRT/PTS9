#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.plotmagneticfield Plot planar cuts through the magnetic field in a SKIRT input model
#
# The function in this module creates PDF vector maps for the planar magnetic field cuts through the SKIRT input model
# produced by the DefaultMagneticFieldCutsProbe and PlanarMagneticFieldCutsProbe probes.
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

## This function creates PDF vector maps for the magnetic field cuts through the SKIRT input model produced by the
# relevant probes. Specifically, the function accepts a single Simulation instance and it assumes that the simulation
# includes one or more DefaultMagneticFieldCutsProbe or PlanarMagneticFieldCutsProbe instances. If this is not the
# case, the function does nothing.
#
# The figure displays an arrow for each bin of Nx by Ny pixels. The bin size can be specified as an argument.
# The orientation and length of the arrow indicate resepctively the direction and strength of the magnetic
# field projected on the cut plane and averaged over the bin. The color of the arrow scales with the magnetic
# field component orthogonal to the cut plane, also averaged over the bin. Vectors pointing away from the
# observer and are red-ish and vectors pointing towards the observer and are blue-ish.
#
# By default, the figures are saved in the simulation output directory with a filename that includes the simulation
# prefix and the probe name, and has the ".pdf" filename extension. The output directory can be overridden as
# described for the pts.utils.savePath() function. In interactive mode (see the pts.utils.interactive() function),
# the figures are not saved and are left open for display in notebooks.
#
# The function takes the following arguments:
#   - simulation: the Simulation instance for which to plot the cuts.
#   - binSize: the number of pixels in each bin, in horizontal and vertical directions.
#   - outDirPath: string or Pathlib.Path object that specifies (overrides) the output directory.
#   - figSize: the horizontal and vertical size of the output figure in inch; default is 6x6 inch.
#   - interactive: whether to leave figures open (True) or save them to file (False).
#
def plotMagneticFieldCuts(simulation, *, binSize=(32,32), outDirPath=None, figSize=(6, 6), interactive=None):

    # find the relevant probes
    probes = [ probe for probe in simulation.probes() \
               if probe.type() in ("DefaultMagneticFieldCutsProbe", "PlanarMagneticFieldCutsProbe") ]

    # iterate over them
    for probe in probes:
        for cut in ("xy", "xz", "yz"):

            # load magnetic field for the this probe and cut
            paths = probe.outFilePaths("{}.fits".format(cut))
            if len(paths) == 1:

                # load data cube with shape (nx, ny, 3)
                Bs = sm.loadFits(paths[0])

                # load the axes grids
                xgrid, ygrid, dummygrid = sm.getFitsAxes(paths[0])
                xmin = xgrid[0].value
                xmax = xgrid[-1].value
                ymin = ygrid[0].value
                ymax = ygrid[-1].value
                extent = (xmin, xmax, ymin, ymax)

                # determine binning configuration
                binX = binSize[0]
                orLenX = Bs.shape[0]
                dropX =  orLenX % binX
                startX = dropX//2
                binY = binSize[1]
                orLenY = Bs.shape[1]
                dropY = orLenY % binY
                startY = dropY//2

                # construct arrays with central bin positions in pixel coordinates
                posX = np.arange(startX - 0.5 + binX / 2.0, orLenX - dropX + startX - 0.5, binX)
                posY = np.arange(startY - 0.5 + binY / 2.0, orLenY - dropY + startY - 0.5, binY)

                # perform the actual binning, while splitting in vector components
                Bx = np.zeros((len(posX),len(posY)))
                By = np.zeros((len(posX),len(posY)))
                Bz = np.zeros((len(posX),len(posY)))
                for x in range(len(posX)):
                    for y in range(len(posY)):
                        Bx[x,y] = np.mean(Bs[startX+binX*x : startX+binX*(x+1), startY+binY*y : startY+binY*(y+1) , 0].value)
                        By[x,y] = np.mean(Bs[startX+binX*x : startX+binX*(x+1), startY+binY*y : startY+binY*(y+1) , 1].value)
                        Bz[x,y] = np.mean(Bs[startX+binX*x : startX+binX*(x+1), startY+binY*y : startY+binY*(y+1) , 2].value)

                # start the figure
                fig, ax = plt.subplots(ncols=1, nrows=1, figsize=figSize)

                # configure the axes
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)
                ax.set_xlabel(cut[0] + sm.latexForUnit(xgrid), fontsize='large')
                ax.set_ylabel(cut[-1] + sm.latexForUnit(ygrid), fontsize='large')
                ax.set_aspect('equal')

                # determine a characteristic 'large' field strength in the cut plane
                Bmax = np.percentile(np.sqrt(Bx**2 + By**2), 99.0)
                if Bmax==0: Bmax=1      # guard against all zeros

                # determine the scaling so that the longest arrows do not to overlap with neighboring arrows
                lengthScale = 2 * Bmax * max(float(len(posX))/figSize[0], float(len(posY))/figSize[1])
                key = "{:.3g}{}".format(Bmax, sm.latexForUnit(Bs))

                # determine the color scheme for the component orthogonal to cut plane
                Bzmax = np.abs(Bz).max()
                if Bzmax==0: Bzmax=1      # guard against all zeros
                normalizer = matplotlib.colors.Normalize(-Bzmax, Bzmax)

                # plot the vector field (scale positions to data coordinates)
                X,Y = np.meshgrid(xmin + posX * (xmax - xmin) / orLenX, ymin + posY * (ymax - ymin) / orLenY, indexing='ij')
                quiverPlot = ax.quiver(X,Y, Bx, By, Bz, cmap='jet', norm=normalizer, pivot='middle', units='inches',
                                       angles='xy', scale=lengthScale, scale_units='inches',
                                       width=0.015, headwidth=2.5, headlength=2, headaxislength=2, minlength=0.8)
                ax.quiverkey(quiverPlot, 0.8, -0.08, Bmax, key, coordinates='axes', labelpos='E')

                # if not in interactive mode, save the figure; otherwise leave it open
                if not ut.interactive(interactive):
                    saveFilePath = ut.savePath(simulation.outFilePath("{}_B_{}.pdf".format(probe.name(),cut)),
                                               (".pdf", ".png"), outDirPath=outDirPath)
                    plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
                    plt.close()
                    logging.info("Created {}".format(saveFilePath))

# -----------------------------------------------------------------
