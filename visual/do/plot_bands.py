#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.do.plot_bands Plot built-in broadband transmission curves
#
# This script creates a plot of the transmission curves for all built-in broadbands that satisfy the
# all of the specified selection criteria:
#  - \em wmin (float): if specified, the pivot wavelength must exceed this value
#  - \em wmax (float): if specified, the pivot wavelength must be lower than this value
#  - \em names (string with comma-separated segments): if specified, the band name must contain
#    at least one of these segments
#
# The plot file is named "FigBuiltinBands.pdf" and is placed in the current working directory.
#

# -----------------------------------------------------------------

def do( wmin : (float,"smallest pivot wavelength to be plotted, in micron") = 0.001,
        wmax : (float,"largest pivot wavelength to be plotted, in micron") = 20000.,
        names : (str,"band name segments for bands to be plotted, comma separated") = "",
        ) -> "plot built-in broadbands in a given wavelength range":

    import astropy.units as u
    import pts.visual as vis
    vis.plotBuiltinBands(minWavelength=wmin * u.micron, maxWavelength=wmax * u.micron,
                         nameSegments=names.split(',') if len(names)>0 else None)

# ----------------------------------------------------------------------
