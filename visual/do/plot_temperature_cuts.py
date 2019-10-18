#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.do.plot_temperature_cuts Plot planar cuts through the dust temperature in one or more
# SKIRT simulations
#
# This script creates plots of the dust temperature cuts generated by the DefaultDustTemperatureCutsProbe
# or the PlanarDustTemperatureCutsProbe in a SKIRT simulation.
#
# If the simulation does not include any DefaultDustTemperatureCutsProbe or PlanarDustTemperatureCutsProbe instances,
# the script does nothing.
#
# The script takes the following arguments:
#  - \em simDirPath (positional string argument): the path to the SKIRT simulation output directory,
#                                                 or "." for the current directory
#  - \em prefix (string): the prefix of the simulation to handle; by default handles all simulations in the directory
#
# In all cases, the plot file is placed next to the simulation output file(s) being handled. The filename includes
# the simulation prefix, the probe name, and the medium indicator, and has the ".pdf" filename extension.
#

# -----------------------------------------------------------------

def do( simDirPath : (str,"SKIRT simulation output directory"),
        prefix : (str,"SKIRT simulation prefix") = "",
        ) -> "plot the planar cuts through the temperature in one or more SKIRT simulations":

    import pts.simulation as sm
    import pts.visual as vis

    for sim in sm.createSimulations(simDirPath, prefix if len(prefix)>0 else None):
        vis.plotDefaultDustTemperatureCuts(sim)

# ----------------------------------------------------------------------
