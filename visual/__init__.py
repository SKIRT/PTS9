#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
#  Package initialization file
# -----------------------------------------------------------------

## \package pts.visual Facilities for visualizing SKIRT-related data through plots and images
#
# This package includes facilities for visualizing SKIRT-related data through plots (e.g. spectra)
# and images (e.g. data frames).
#

from .makergbimages import makeRGBImages, makeConvolvedRGBImages
from .makewavelengthmovie import makeWavelengthMovie
from .moviefile import MovieFile
from .plotbands import plotBuiltinBands
from .plotcurves import plotSeds, plotSources, plotSpectralResolution
from .plotgrids import plotGrids
from .plotmagneticfield import plotMagneticField
from .plotpolarization import plotPolarization
from .plotscalarcuts import plotScalarCuts
from .plotstoredtable import plotStoredTableCurve, plotStoredTableInteractive
from .plotvelocity import plotVelocity
from .rgbimage import RGBImage
