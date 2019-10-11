#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.do.plot_spectral_resolution Plot the spectral resolution of a wavelength axis
#
# This script creates a plot of the spectral resolution \f$R=\lambda/\Delta\lambda\f$ of a wavelength grid
# loaded from a SKIRT stored table (".stab") or a SKIRT output file (".dat") that includes a wavelength axis.
#
# Provide the name or (relative or absolute) file path of a SKIRT stored table file or a SKIRT output file
# as the single required argument, and use optional arguments to further configure the plot:
#  - \em wmin (float): smallest wavelength on the horizontal axis, in micron; default is 0.1 micron
#  - \em wmax (float): largest wavelength on the horizontal axis, in micron; default is 1000 micron
#
# The plot file has the same name as the input file but with the ".pdf" filename extension and is placed
# in the current working directory.
#

# -----------------------------------------------------------------

def do( filepath : (str,"name or path of SKIRT stored table or text column file"),
        wmin : (float,"smallest wavelength on the horizontal axis, in micron") = 0.1,
        wmax : (float,"largest wavelength on the horizontal axis, in micron") = 1000.,
        ) -> "plot the spectral resolution of a wavelength axis":

    import astropy.units as u
    import pts.visual as vis

    vis.plotSpectralResolution(inFilePath=filepath, minWavelength=wmin * u.micron, maxWavelength=wmax * u.micron)

# -----------------------------------------------------------------
