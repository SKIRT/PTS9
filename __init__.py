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
# This information may be moved elsewhere over time...
#
# __Organization of PTS packages__
#
# Each PTS package (directory just inside the top-level pts directory) exposes all public functions and classes
# (i.e. those intended for use outside of the package) at the package level. The functionality is implemented
# in various modules (python source files) residing inside the package. The initialization file for each package
# explicitly places the public names into the package namespace using explicit imports.
#
# __Import styles used in PTS code__
#
# Default style for importing external packages (including standard-library packages):
#
#     import some.package           # each reference must include full package name
#
# External packages imported with a local name:
#
#     import astropy.io.fits as fits
#     import astropy.units as u
#     import lxml.etree as etree
#     import matplotlib.pyplot as plt
#     import numpy as np
#
# Importing other PTS packages (or same package from within do subdirectory):
#
#     import pts.admin as adm
#     import pts.band as bnd
#     import pts.do as do
#     import pts.simulation as sm
#     import pts.storedtable as stab
#     import pts.utils as ut
#     import pts.visual as vis
#
# Importing symbols from within the same package, including initialization file:
#
#     from .module import name      # default style is to use explicit import
#     from .module import *         # exceptional style, for example in conversionspec.py
#


