#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.simulation.fits Handling SKIRT FITS output files
#
# This module offers functions for reading frames and data cubes from SKIRT FITS output files with support for units,
# and for writing frames and data cubes to FITS files with the same structure as SKIRT output files.
#

# -----------------------------------------------------------------

import datetime
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

## This function writes a frame or data cube to a FITS file that resembles a SKIRT output file.
# It expects the following arguments:
#  - outFilePath: the output file's absolute or relative path, including filename and '.fits' extension,
#    interpreted as described for the pts.utils.absPath() function.
#  - data: a two- or three-dimensional astropy quantity array; the indices of the array are in the order [x, y] for
#    a data frame and [x, y, z] for a data cube, where usually x and y are spatial axes and z the wavelength axis.
#  - xaxis (optional): a 1D astropy quantity array listing the grid points for the x-axis.
#  - yaxis (optional): a 1D astropy quantity array listing the grid points for the y-axis.
#  - zaxis (optional; used only if \em data is 3D): a 1D astropy quantity array listing the grid points for the z-axis.
#
# The grid points for the x- and y-axis are assumed to have a regular linear distribution. The grid points for the
# z-axis can have any distribution; all points are included in the FITS file as a table.
#
# The function includes information in the FITS header fields about the data and axes in the same way as SKIRT does
# for its output files. If the information for a given axis is missing, the corresponding information is simply not
# written to the FITS file. The information for the z-axis is written only if the \em data array has three dimensions.
#
# The current implementation does \em not include information on the line of sight and distance from the model.
#
def writeFits(outFilePath, data, xaxis=None, yaxis=None, zaxis=None):
    outpath = ut.absPath(outFilePath)

    # verify some of the requirements/restrictions on the specified data
    if outpath.suffix != ".fits":
        raise ValueError("FITS filename extension is not '.fits': {}".format(outpath))
    numAxes = len(data.shape)
    if numAxes<2 or numAxes>3:
        raise ValueError("The data array has {} dimensions rather than 2 or 3".format(numAxes))
    if xaxis is not None and len(xaxis) != data.shape[0]:
        raise ValueError("The x-grid has {} grid points while the data array has {}".format(len(xaxis), data.shape[0]))
    if yaxis is not None and len(yaxis) != data.shape[1]:
        raise ValueError("The y-grid has {} grid points while the data array has {}".format(len(yaxis), data.shape[1]))
    if numAxes == 3 and zaxis is not None and len(zaxis) != data.shape[2]:
        raise ValueError("The z-grid has {} grid points while the data array has {}".format(len(zaxis), data.shape[2]))

    # store the data in a primary HDU, and put it in the HDU list
    hdu = fits.PrimaryHDU(data.value.T.astype(np.float32))
    hdul = fits.HDUList([hdu])

    # add the basic header fields
    hdr = hdu.header
    hdr["BSCALE"] = (1, "Array value scale")
    hdr["BZERO"]  = (0, "Array value offset")
    hdr['DATE'] = (datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), "Date and time of creation")
    hdr['ORIGIN'] = ("Python Toolkit SKIRT", "Astronomical Observatory, Ghent University")
    hdr['BUNIT'] = (data.unit.to_string(), "Physical unit of the array values")

    # add the header fields for the x and y axes, if present
    if xaxis is not None:
        hdr['CRPIX1'] = ((len(xaxis)+1)/2, "X-axis coordinate system reference pixel")
        hdr['CRVAL1'] = ((xaxis[0]+xaxis[-1]).value/2, "Coordinate value at X-axis reference pixel")
        hdr['CDELT1'] = ((xaxis[1]-xaxis[0]).value, "Coordinate increment along X-axis")
        hdr['CUNIT1'] = (xaxis.unit.to_string(), "Physical unit of the X-axis")
        hdr['CTYPE1'] = (" ", "Linear X coordinates")
    if yaxis is not None:
        hdr['CRPIX2'] = ((len(yaxis)+1)/2, "Y-axis coordinate system reference pixel")
        hdr['CRVAL2'] = ((yaxis[0]+yaxis[-1]).value/2, "Coordinate value at Y-axis reference pixel")
        hdr['CDELT2'] = ((yaxis[1]-yaxis[0]).value, "Coordinate increment along Y-axis")
        hdr['CUNIT2'] = (yaxis.unit.to_string(), "Physical unit of the Y-axis")
        hdr['CTYPE2'] = (" ", "Linear Y coordinates")

    # store the full z-axis as a table, if present
    if numAxes == 3 and zaxis is not None:
        hdr['CUNIT3'] = (zaxis.unit.to_string(), "Physical unit of the Z-axis")

        # create and add extension HDU containing the z grid points as a table
        zcol = fits.Column(name="GRID_POINTS", array=zaxis.value, unit=zaxis.unit.to_string(), format='D')
        zhdu = fits.BinTableHDU.from_columns([zcol], name="Z-axis coordinate values")
        hdul.append(zhdu)

    # write the file
    hdul.writeto(outpath, overwrite=True)

# -----------------------------------------------------------------
