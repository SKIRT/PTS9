#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
#  Package initialization file
# -----------------------------------------------------------------

## \package pts.band Facilities for representing broadband filters, including transmission curve data
#
# This package includes facilities for representing broadband filters, including transmission curve data
# for a set of standard bands, options to produce basic plots for these curves, and the capability to
# convolve SEDs or frame data cubes with the transmission curves.
#

from .broadband import BroadBand, builtinBands, builtinBand
