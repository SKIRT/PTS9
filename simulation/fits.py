#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.simulation.fits Handling SKIRT FITS output files
#
# This module offers functions for reading frames and data cubes from SKIRT FITS output files with support for units.
#

# -----------------------------------------------------------------

import astropy.io.fits as fits
import numpy as np
import pts.utils as ut
from .units import unit as smunit

# -----------------------------------------------------------------

## This function loads the complete contents from a FITS file produced by SKIRT, and returns the resulting data frame
# or data cube as a two- or three-dimensional astropy quantity array. The indices of the array are in the order
# [x, y] for a data frame and [x, y, z] for a data cube, where x and y are usually spatial axes and z is often the
# wavelength axis. The appropriate units for the array values are obtained from the FITS file header written by SKIRT.
#
# The file path is interpreted as described for the pts.utils.absPath() function.
#
def loadFits(path):
    path = ut.absPath(path)
    data, header = fits.getdata(path, header=True)
    return data.T.astype(float) << smunit(header['BUNIT'])

# -----------------------------------------------------------------

## This function retrieves the axis grid points from a FITS file produced by SKIRT, and returns the result as a
# tuple of astropy quantity arrays: (x, y) for data frames or (x, y, z) for data cubes.
# The function assumes that x and y represent spatial axes with a regular, linear grid. For a data cube, the
# function expects the z-axis grid points (often wavelengths) to be listed in a binary table extension.
#
# The file path is interpreted as described for the pts.utils.absPath() function.
#
def getFitsAxes(path):
    # open the file
    path = ut.absPath(path)
    with fits.open(path) as hdul:
        # build x and y grids
        h = hdul[0].header
        x = _grid(h['NAXIS1'], h['CRPIX1'], h['CRVAL1'], h['CDELT1'], h['CUNIT1'])
        y = _grid(h['NAXIS2'], h['CRPIX2'], h['CRVAL2'], h['CDELT2'], h['CUNIT2'])
        # if there are only two axes, we're done
        if int(h['NAXIS']) == 2:
            return x,y

        # if there are three axis, read the z grid from the table extension
        hdu = hdul["Z-axis coordinate values"]
        z = hdu.data["GRID_POINTS"].astype(float) << smunit(hdu.header["TUNIT1"])
        return x,y,z

## This helper function constructs a regular linear grid based on the information for an axis in the FITS header
def _grid(npix, crpix, crval, cdelt, cunit):
    npix = int(npix)
    crpix = float(crpix)
    crval = float(crval)
    cdelt = float(cdelt)
    start = crval - cdelt * (crpix-1)
    stop = crval + cdelt * (npix-crpix)
    return np.linspace(start, stop, npix) << smunit(cunit)

# -----------------------------------------------------------------
