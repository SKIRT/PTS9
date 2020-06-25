#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.simulation.do.healpix_to_planck Create a Planck compatible all-sky FITS
# file based on the raw output of a HEALPixSkyInstrument.
#
# This script converts the given HEALPixSkyInstrument output FITS file to a Planck compatible
# all-sky FITS file. The image data in the original FITS file is moved to a FITS table, so that
# the resulting FITS file can be analysed using the same tools as were used for Planck.
#
# The function takes the following arguments:
#  - \em inputFileName (string): the name of the file containing the raw HEALPix data cube
#  - \em outputFileName (string): the name of the generated Planck compatible FITS file
#

# -----------------------------------------------------------------


def do(
    inputFileName: (
        str,
        "HEALPixSkyInstrument output image containing the raw HEALPix pixel data",
    ),
    outputFileName: (
        str,
        "output file containing a Planck compatible version of the data cube",
    ),
) -> "create a Planck compatible FITS file from a raw HEALPixSkyInstrument output file":

    import numpy as np
    import astropy.io.fits as fits
    import pts.simulation.healpix as healpix
    import pathlib
    from pts.utils.error import UserError
    import logging

    inputPath = pathlib.Path(inputFileName)
    if not inputPath.exists():
        raise UserError("Input file does not exist!")
    if not inputPath.suffix == ".fits":
        raise UserError("Input file is not a FITS file!")

    outputPath = pathlib.Path(outputFileName)
    if outputPath.exists():
        raise UserError("Output file already exists!")
    if not outputPath.suffix == ".fits":
        raise UserError("Output file is not a FITS file!")

    inputFile = fits.open(inputPath)
    logging.info("Opened input file {0}".format(inputPath))

    hData = healpix.HEALPixGrid(inputFile[0].data)
    hData.printInfo()
    j, i = hData.getPixelIndices()
    validCube = hData.getPixel(j, i)
    logging.info("Extracted actual data cube of shape {0}".format(validCube.shape))

    wav = inputFile[1].data[:]

    if len(validCube.shape) == 1:
        validCube.reshape((1, validCube.shape[0]))

    cols = []
    for iCol in range(len(wav)):
        col = fits.Column(
            name="{0:.0f} mu".format(wav[iCol][0]),
            format="D",
            array=validCube[iCol],
            unit=inputFile[0].header["BUNIT"],
        )
        cols.append(col)
    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.header["PIXTYPE"] = "HEALPIX"
    hdu.header["ORDERING"] = "RING"
    hdu.header["NSIDE"] = hData._NSide

    outputHDUL = fits.HDUList([fits.PrimaryHDU(), hdu])

    outputHDUL.writeto(outputPath)
    logging.info("Wrote output file {0}".format(outputPath))
