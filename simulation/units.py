#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.simulation.units Handling SKIRT units
#
# This module offers functions for converting SKIRT units to astropy units, and for handling conversions that
# may not be trivial to perform using pure astropy.
#

# -----------------------------------------------------------------

import astropy.units as u
import warnings

# -----------------------------------------------------------------

## This function returns an astropy unit for any of the following "unit-like" inputs:
#  - a string representing a valid SKIRT unit: returns the corresponding astropy unit.
#  - a string representing a valid astropy unit in generic format: returns the corresponding astropy unit,
#    except that 'A' is usually interpreted as Angstrom instead of Ampere.
#  - an astropy quantity: returns the unit of the quantity.
#  - an astropy unit: simply returns the unit itself.
#
# While the function is mainly intended for converting SKIRT units to astropy units, supporting the other
# input types is useful when the function is invoked from some of the other functions in this module.
#
def unit(unitlike):
    if isinstance(unitlike, str):
        # handle SKIRT-specific pressure equivalence: 1 K/m3 <==> k_B Pa
        if unitlike=="K/m3":
            unitlike = "1.3806488e-23 Pa"
        # replace SKIRT-specific symbol for Angstrom
        elif 'A' in unitlike:
            unitlike = "/".join([ s if s!='A' else "Angstrom" for s in unitlike.split('/') ])
        # parse string ignoring warnings about multiple divisions, as in "W/m2/sr"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=u.UnitsWarning)
            return u.Unit(unitlike)
    if isinstance(unitlike, u.Quantity):
        return unitlike.unit
    if isinstance(unitlike, u.UnitBase):
        return unitlike
    raise ValueError("Unsupported unit-like type: {}".format(type(unitlike)))

## This function returns a latex-formatted string representation for the "unit-like" input,
# enclosed in square brackets and preceded by a short space. The input argument
# is interpreted as described for the unit() function in this module.
def latex(unitlike):
    return r"$\,$[" + unit(unitlike).to_string("latex_inline") + r"]"

# -----------------------------------------------------------------
