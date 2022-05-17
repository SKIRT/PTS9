#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.plotconvergencecuts Plot planar cuts through the input and gridded medium density in a SKIRT simulation
#
# The function in this module creates plots of the planar cuts through the input and gridded medium density
# produced by the convergence cuts probe in a SKIRT simulation. This allows a visual comparison between the medium
# density distribution of the input model and the corresponding spatially discretized version used in the simulation.
#

# -----------------------------------------------------------------

import logging
import matplotlib.colors
import matplotlib.pyplot as plt
import pts.simulation as sm
import pts.utils as ut

# -----------------------------------------------------------------

## This function creates plots of the planar cuts through the input and gridded medium density produced by the
# convergence cuts probe in a SKIRT simulation. This allows a visual comparison between the medium density distribution
# of the input model and the corresponding spatially discretized version used in the simulation.

# The function accepts a single Simulation instance and it assumes that the simulation includes one or more
# ConvergenceCutsProbe instances. If this is not the case, the function does nothing. If the simulation also includes
# one or more ConvergenceInfoProbe instances, information about the total input and gridded mass is extracted from
# the output of those probes and displayed on the plot.
#
# A single plot is produced for each medium (dust, electrons, gas) including panels for both input and gridded density
# (rows) and for all orientations (columns) depending on the symmetry of the model. All panels share the same color
# scale. The dynamic range of the logarithmic color scale is specified by \em decades.
#
# By default, the figures are saved in the simulation output directory with a filename that includes the simulation
# prefix, the probe name, and the medium type indicator, and has the ".pdf" filename extension.
# This can be overridden with the out* arguments as described for the pts.utils.savePath() function.
# In interactive mode (see the pts.utils.interactive() function), the figures are not saved and are left open
# for display in notebooks.
#
# The function takes the following arguments:
#   - simulation: the Simulation instance for which to plot the vector quantity.
#   - decades: the dynamic range of a logarithmic color scale in dex, or zero for linear color scale; default is 5.
#   - outDirPath: string or Pathlib.Path object that overrides the output directory.
#   - outFileName: string or Pathlib.Path object that overrides the output file name.
#   - outFilePath: string or Pathlib.Path object that overrides the output file path.
#   - figSize: the horizontal and vertical size of the output figure in inch;
#              default is 7*num-orientations x 10 inch.
#   - interactive: whether to leave figures open (True) or save them to file (False).
#
def plotConvergenceCuts(simulation, decades=5, *,
                        outDirPath=None, outFileName=None, outFilePath=None, figSize=None, interactive=None):

    # find the convergence cuts probes and iterate over them
    for probe in simulation.probes("ConvergenceCutsProbe"):

        # iterate over the medium types for which output has been produced
        for medium in ("dust", "elec", "gas"):
            if len(probe.outFilePaths("{}_t_xy.fits".format(medium))) > 0:

                # determine the cuts for which output has been produced
                cuts = []
                for cut in ("xy", "xz", "yz"):
                    if len(probe.outFilePaths("{}_t_{}.fits".format(medium,cut))) > 0:
                        cuts.append(cut)

                # load the input and gridded cuts for the requested medium
                numcuts = len(cuts)
                paths = [ probe.outFilePaths("{}_t_{}.fits".format(medium,cut))[0] for cut in cuts ] + \
                        [ probe.outFilePaths("{}_g_{}.fits".format(medium,cut))[0] for cut in cuts ]
                frames = [ sm.loadFits(path) for path in paths ]
                grids = [ sm.getFitsAxes(path) for path in paths ]

                # determine the range of values to display and clip the data arrays
                vmax = max([ frame.max() for frame in frames ])
                vmin = vmax - vmax  # to retain units
                if decades>0:
                    vmin = vmax / 10**decades
                    for frame in frames:
                        frame[ frame<vmin ] = vmin

                # setup the figure depending on the number of cuts
                fig, axes = plt.subplots(ncols=numcuts, nrows=2, figsize=figSize if figSize is not None else (7*numcuts,10))
                if numcuts>1: axes = list(axes[0]) + list(axes[1])
                else: axes = list(axes)

                # set normalizer (logarithmic normalizer crashes if all values are zero)
                if decades>0 and vmax>0:
                    normalizer = matplotlib.colors.LogNorm(vmin.value, vmax.value)
                else:
                    normalizer = matplotlib.colors.Normalize(vmin.value, vmax.value)

                # plot the frames and set axis details for each
                for ax, cut, frame, (xgrid, ygrid) in zip(axes, cuts+cuts, frames, grids):
                    # plot the frame
                    extent = (xgrid[0].value, xgrid[-1].value, ygrid[0].value, ygrid[-1].value)
                    im = ax.imshow(frame.value.T, norm=normalizer, cmap='gnuplot', extent=extent,
                                   aspect='equal', interpolation='bicubic', origin='lower')
                    # set axis limits
                    ax.set_xlim(xgrid[0].value, xgrid[-1].value)
                    ax.set_ylim(ygrid[0].value, ygrid[-1].value)
                    # set axis labels
                    ax.set_xlabel(cut[0] + sm.latexForUnit(xgrid.unit), fontsize='large')
                    ax.set_ylabel(cut[-1] + sm.latexForUnit(ygrid.unit), fontsize='large')

                # add a color bar
                symbol = r"$\rho$" if medium == "dust" else "n"
                cbarax = fig.colorbar(im, ax=axes, shrink=0.5, anchor=(0,0)).ax
                cbarax.set_ylabel(symbol + sm.latexForUnit(frame.unit), fontsize='large')

                # find the path to an appropriate convergence info data file, if present
                infopath = _convergenceInfoPath(simulation, probe)
                if infopath is not None:
                    # get the input and gridded mass for the appropriate medium
                    m_type = "Dust"
                    if medium=="elec": m_type = "Electrons"
                    if medium=="gas": m_type = "Gas"
                    m_in = sm.getQuantityFromFile(infopath, m_type+"/Total mass", "Input")
                    m_gr = sm.getQuantityFromFile(infopath, m_type+"/Total mass", "Gridded")

                    # add info on input and gridded mass
                    keys = dict(transform=cbarax.transAxes, horizontalalignment='center', fontsize='large')
                    cbarax.text(2, 1.92, m_type, **keys)
                    cbarax.text(2, 1.80, r"$M_\mathrm{input}$", **keys)
                    cbarax.text(2, 1.72, r"${}$".format(m_in.value) + sm.latexForUnit(m_in), **keys)
                    cbarax.text(2, 1.60, r"$M_\mathrm{gridded}$", **keys)
                    cbarax.text(2, 1.52, r"${}$".format(m_gr.value) + sm.latexForUnit(m_gr), **keys)
                    if m_in > 0:
                        cbarax.text(2, 1.44, r"$({:4.2} \%)$".format(100 * (m_gr/m_in-1)), **keys)

                # if not in interactive mode, save the figure; otherwise leave it open
                if not ut.interactive(interactive):
                    basepath = simulation.outFilePath("{}_{}.pdf".format(probe.name(), medium))
                    saveFilePath = ut.savePath(basepath, (".pdf", ".png"),
                                               outDirPath=outDirPath, outFileName=outFileName, outFilePath=outFilePath)
                    plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
                    plt.close()
                    logging.info("Created {}".format(saveFilePath))

# ----------------------------------------------------------------------

## This private function returns the path to the output file of the first convergence info probe that has the same
# \em probeAfter property value as the requesting convergence cuts probe, or None if there is no such probe.
def _convergenceInfoPath(simulation, cutsprobe):
    # get the probeAter value of the cuts probe
    requested = cutsprobe.getStringAttribute("probeAfter")

    # loop over all convergence info probes
    for infoprobe, infopath in sm.probeOutFilePaths(simulation.probes("ConvergenceInfoProbe"), ".dat"):
        # get the probeAter value of the info probe
        candidate = infoprobe.getStringAttribute("probeAfter")

        # return if this is a hit
        if candidate == requested:
            return infopath

    # indicate failure
    return None

# ----------------------------------------------------------------------
