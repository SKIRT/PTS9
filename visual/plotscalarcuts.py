#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.plotscalarcuts Plot planar cuts or projections for a scalar quantity from a SKIRT simulation
#
# The function in this module creates plots of the planar cuts or projections for a scalar quantity
# produced by one of the relevant probes in a SKIRT simulation.
#

# -----------------------------------------------------------------

import logging
import matplotlib.colors
import matplotlib.pyplot as plt
import pts.simulation as sm
import pts.utils as ut

# -----------------------------------------------------------------

## This function creates plots of the planar cuts or projections for a scalar quantity produced by one of the relevant
# probes in a SKIRT simulation. Specifically, the function accepts a single Simulation instance and it assumes that
# the simulation includes one or more instances of one or more of the probe types listed in \em probeTypes,
# with an associated probe form that produces a planar cut (DefaultCutsForm, PlanarCutsForm) or planar
# projection (ParallelProjectionForm, AllSkyProjectionForm). If this is not the case, the function does nothing.
#
# For \em decades=0, the color scale is linear. For \em decades>0, the color scale is logarithmic with the given
# dynamic range in dex.
#
# If the output for a given medium component or medium type (dust, gas, or electrons) includes two or three files
# with names that differ only by orientation labels ("_xy", "_xz", and/or "_yz"), the maps for these files are
# included in a single plot and share the same color scale.
#
# By default, the figure is saved in the simulation output directory with a filename that includes the simulation
# prefix, the probe name, and the medium component or type indicator, and has the ".pdf" filename extension.
# This can be overridden with the out* arguments as described for the pts.utils.savePath() function.
# In interactive mode (see the pts.utils.interactive() function), the figure is not saved and it is left open
# so that is displayed in notebooks.
#
def plotScalarCuts(simulation, probeTypes, decades=0, *,
                   outDirPath=None, outFileName=None, outFilePath=None, figSize=None, interactive=None):

    # find the relevant probes
    probes = simulation.probes(probeTypes,
                    ("DefaultCutsForm", "PlanarCutsForm", "ParallelProjectionForm", "AllSkyProjectionForm"))

    # find the list of sorted output file paths and map them to the corresponding probe
    pathdict = { path: probe for probe, path in sm.probeOutFilePaths(probes, ".fits") }
    paths = sorted(pathdict.keys())

    # iterate over the paths, grouping 1, 2 or 3 "related" paths
    index = 0
    while index < len(paths):
        pathgroup = [ paths[index] ]
        index += 1
        if index < len(paths) and _related(pathgroup[0], paths[index]):
            pathgroup.append(paths[index])
            index += 1
            if index < len(paths) and _related(pathgroup[1], paths[index]):
                pathgroup.append(paths[index])
                index += 1

        # load the data frames and the range of the x and y axes
        # (there can be one to three frames depending on symmetries)
        numframes = len(pathgroup)
        frames = [ sm.loadFits(path) for path in pathgroup ]
        grids = [ sm.getFitsAxes(path) for path in pathgroup ]

        # determine the range of values to display and clip the data arrays
        vmax = max([ frame.max() for frame in frames ])
        vmin = vmax - vmax  # to retain units
        if decades>0:
            vmin = vmax / 10**decades
            for frame in frames:
                frame[ frame<vmin ] = vmin

        # setup the figure depending on the number of frames
        fig, axes = plt.subplots(ncols=numframes, nrows=1, figsize=figSize if figSize is not None else (8*numframes,6))
        if numframes==1: axes = [axes]

        # plot the frames and set axis details for each
        for ax, path, frame, (xgrid, ygrid) in zip(axes, pathgroup, frames, grids):
            # set normalizer (logarithmic normalizer crashes if all values are zero)
            if decades>0 and vmax>0:
                normalizer = matplotlib.colors.LogNorm(vmin.value, vmax.value)
            else:
                normalizer = matplotlib.colors.Normalize(vmin.value, vmax.value)
            # plot the frame
            extent = (xgrid[0].value, xgrid[-1].value, ygrid[0].value, ygrid[-1].value)
            im = ax.imshow(frame.value.T, norm=normalizer, cmap='gnuplot', extent=extent,
                           aspect='equal', interpolation='bicubic', origin='lower')
            # set axis limits
            ax.set_xlim(xgrid[0].value, xgrid[-1].value)
            ax.set_ylim(ygrid[0].value, ygrid[-1].value)
            # set axis labels
            if "AllSky" in pathdict[path].formType():
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
            ax.set_xlabel(xlabel + sm.latexForUnit(xgrid.unit), fontsize='large')
            ax.set_ylabel(ylabel + sm.latexForUnit(ygrid.unit), fontsize='large')

        # add a color bar with a label derived from the filename's last or one but last segment
        segments = path.stem.split("_")
        symbol = segments[-1]
        if symbol == "xy" or symbol == "xz" or symbol == "yz":
            symbol = segments[-2]
        if len(symbol)>1:
            symbol = r"$\{}$".format(symbol)
        fig.colorbar(im, ax=axes).ax.set_ylabel(symbol + sm.latexForUnit(frame.unit), fontsize='large')

        # if not in interactive mode, save the figure; otherwise leave it open
        if not ut.interactive(interactive):
            basepath = pathgroup[0].with_suffix(".pdf")
            if numframes>1:
                basepath = basepath.with_stem(_reduce(basepath))
            saveFilePath = ut.savePath(basepath, (".pdf", ".png"),
                                       outDirPath=outDirPath, outFileName=outFileName, outFilePath=outFilePath)
            plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
            plt.close()
            logging.info("Created {}".format(saveFilePath))

# ----------------------------------------------------------------------

# This private function returns true if the stems of the two paths are identical
# except for variations of "_xy", "_xz" and "_yz".
def _related(path1, path2):
    cut1 = "_xy"
    cut2 = "_xz"
    cut3 = "_yz"
    stem1 = path1.stem
    stem2 = path2.stem
    if (cut1 in stem1 or cut2 in stem1 or cut3 in stem1) and (cut1 in stem2 or cut2 in stem2 or cut3 in stem2):
        stem1 = stem1.replace(cut1, "?", 1).replace(cut2, "?", 1).replace(cut3, "?", 1)
        stem2 = stem2.replace(cut1, "?", 1).replace(cut2, "?", 1).replace(cut3, "?", 1)
        if stem1 == stem2:
            return True
    return False

# This private function returns the stem of the path after removing "_xy", "_xz" and "_yz"
def _reduce(path):
    return path.stem.replace("_xy", "", 1).replace("_xz", "", 1).replace("_yz", "", 1)

# ----------------------------------------------------------------------
