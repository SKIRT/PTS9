#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.do.plot_bands Plot built-in broadband transmission curves
#
# This script creates a plot of the transmission curves for all built-in broadbands that satisfy
# all of the specified selection criteria:
#  - \em names (string with comma-separated segments): if specified, the band name must contain
#    at least one of these segments
#  - \em wmin (float): if specified, the pivot wavelength must exceed this value
#  - \em wmax (float): if specified, the pivot wavelength must be lower than this value
#
# The plot file is named "FigBuiltinBands.pdf" and is placed in the current working directory.
#

# -----------------------------------------------------------------

def do( names : (str,"band name segments for bands to be plotted, comma separated") = "",
        wmin : (float,"smallest pivot wavelength to be plotted, in micron") = 1e-99,
        wmax : (float,"largest pivot wavelength to be plotted, in micron") = 1e99,
        ) -> "plot built-in broadbands":

    import astropy.units as u
    import pts.visual as vis
    vis.plotBuiltinBands(wmin << u.micron, wmax << u.micron, names)

# ----------------------------------------------------------------------
