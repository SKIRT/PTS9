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
#  - a string representing a valid astropy unit in generic format: returns the corresponding astropy unit.
#  - an astropy quantity: returns the unit of the quantity.
#  - an astropy unit: simply returns the unit itself.
#
# While the function is mainly intended for converting SKIRT units to astropy units, supporting the other
# input types is useful when the function is invoked directly or indirectly from some other context.
#
def unit(unitlike):
    if isinstance(unitlike, str):
        # handle SKIRT-specific pressure equivalence: 1 K/m3 <==> k_B Pa
        if unitlike=="K/m3":
            unitlike = "1.3806488e-23 Pa"
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
def latexForUnit(unitlike):
    return r"$\,$[" + unit(unitlike).to_string("latex_inline") + r"]"

## This function returns a latex-formatted string representation for a flux density (\f$\lambda F_\lambda\f$,
# \f$F_\lambda\f$, or \f$F_\nu\f$) or a surface brightness (\f$\lambda f_\lambda\f$, \f$f_\lambda\f$, or \f$f_\nu\f$)
# that has the units represented by the "unit-like" input (the latex string does not include the unit itself).
# The input argument is interpreted as described for the unit() function in this module. If no match is found,
# the function returns an empty string.
def latexForSpectralFlux(unitlike):
    un = unit(unitlike)
    if un.is_equivalent(unit("W/m2")): return r"$\lambda\,F_\lambda$"
    if un.is_equivalent(unit("W/m2/m")): return r"$F_\lambda$"
    if un.is_equivalent(unit("W/m2/Hz")): return r"$F_\nu$"
    if un.is_equivalent(unit("W/m2/sr")): return r"$\lambda\,f_\lambda$"
    if un.is_equivalent(unit("W/m2/sr/m")): return r"$f_\lambda$"
    if un.is_equivalent(unit("W/m2/sr/Hz")): return r"$f_\nu$"
    return ""

## This function returns a latex-formatted string representation for the mean intensity or spectral radiance
# (\f$\lambda J_\lambda\f$, \f$J_\lambda\f$, or \f$J_\nu\f$)
# that has the units represented by the "unit-like" input (the latex string does not include the unit itself).
# The input argument is interpreted as described for the unit() function in this module. If no match is found,
# the function returns an empty string.
def latexForSpectralRadiance(unitlike):
    un = unit(unitlike)
    if un.is_equivalent(unit("W/m2/sr")): return r"$\lambda\,J_\lambda$"
    if un.is_equivalent(unit("W/m2/sr/m")): return r"$J_\lambda$"
    if un.is_equivalent(unit("W/m2/sr/Hz")): return r"$J_\nu$"
    return ""

## This function returns a latex-formatted string representation for the spectral luminosity
# (\f$\lambda L_\lambda\f$, \f$L_\lambda\f$, or \f$L_\nu\f$)
# that has the units represented by the "unit-like" input (the latex string does not include the unit itself).
# The input argument is interpreted as described for the unit() function in this module. If no match is found,
# the function returns an empty string.
def latexForSpectralLuminosity(unitlike):
    un = unit(unitlike)
    if un.is_equivalent(unit("W")): return r"$\lambda\,L_\lambda$"
    if un.is_equivalent(unit("W/m")): return r"$L_\lambda$"
    if un.is_equivalent(unit("W/Hz")): return r"$L_\nu$"
    return ""

# -----------------------------------------------------------------
