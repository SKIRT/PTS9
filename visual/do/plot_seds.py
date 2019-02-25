#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.do.plot_seds Plot the SEDs produced by one or more SKIRT simulations
#
# This script creates plots of the SEDs produced by one or more SKIRT simulations. It takes the following arguments:
#  - \em simDirPath (positional string argument): the path to the SKIRT simulation output directory,
#                                                 or "." for the current directory
#  - \em prefix (string): the prefix of the simulation to handle; by default handles all simulations in the directory
#  - \em wmin (float): smallest wavelength on the horizontal axis, in micron; default is 0.1 micron
#  - \em wmax (float): largest wavelength on the horizontal axis, in micron; default is 1000 micron
#  - \em dex (float): if specified, the number of decades to be plotted on the vertical axis; default is 5
#
# In all cases, the plot file is placed next to the "*_sed.dat" file(s) being represented. The filename starts
# with the simulation prefix and ends with "sed.pdf".
#

# -----------------------------------------------------------------

def do( simDirPath : (str,"SKIRT simulation output directory"),
        prefix : (str,"SKIRT simulation prefix") = "",
        wmin : (float,"smallest wavelength on the horizontal axis, in micron") = 0.1,
        wmax : (float,"largest wavelength on the horizontal axis, in micron") = 1000.,
        dex : (float,"number of decades to be plotted on the vertical axis") = 5,
        ) -> "plot the SEDs produced by one or more SKIRT simulations":

    import astropy.units as u
    import pts.simulation as sm
    import pts.visual as vis

    for sim in sm.createSimulations(simDirPath, prefix if len(prefix)>0 else None):
        vis.plotSeds(sim, minWavelength=wmin * u.micron, maxWavelength=wmax * u.micron, decades=dex)

# ----------------------------------------------------------------------
