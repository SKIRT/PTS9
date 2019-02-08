#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
#  Package initialization file
# -----------------------------------------------------------------

## \package pts Python toolkit for working with SKIRT (PTS)
#
# The Python toolkit for working with SKIRT (PTS) includes a set of Python packages offering functionality
# related to working with SKIRT, the advanced 3D Monte Carlo radiative transfer code also developed
# by the Astronomical Observatory at the Ghent University.
#
#
# Developer Guidelines
# --------------------
#
# This information will move elsewhere over time...
#
# PTS code uses one of three import styles:
#
# import style           | comments
# -----------------------|---------
# import module          | default style; each reference must include full module name
# import module as local | style with local name for frequently-used modules as listed below
# from module import *   | exceptional style for pulling all symbols into the local namespace
#
# The following frequently-used modules are imported with a local name:
#
#     import astropy.io.fits as fits
#     import astropy.units as u
#     import lxml.etree as etree
#     import matplotlib.pyplot as plt
#     import numpy as np
#
#     import pts.band.broadband as bb
#     import pts.simulation.simulation as sim
#     import pts.simulation.skifile as ski
#     import pts.simulation.skirt as sk
#     import pts.simulation.units as su
#     import pts.utils.error as pe
#     import pts.utils.path as pp
#
#

