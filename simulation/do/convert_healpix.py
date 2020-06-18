#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.simulation.do.convert_healpix Create a projected image based on a HEALPixSkyInstrument output cube
#
# This script creates a projected image based on the raw output of a HEALPixSkyInstrument.
#
# The function takes the following arguments:
#  - \em inputFileName (string): the name of the file containing the raw HEALPix data cube
#  - \em outputFileName (string): the name of the generated projected image
#  - \em projection (string): the projection to use. Accepted values are Mollweide (for the Mollweide projection)
#    and HammerAitoff (for the Hammer-Aitoff projection); default is Mollweide.
#  - \em nPixelX (int): Number of pixels in the vertical direction for the generated image. The number of pixels
#    in the horizontal direction is twice this value. Default is 250.
#  - \em thetaCenter (float): Zenith angle of the central position relative to the original crosshair of the
#    HEALPixSkyInstrument (in degrees). Corresponds to a vertical rotation that goes downward in the right half
#    of the image. Default is 0.
#  - \em phiCenter (float): Azimuth angle of the central position relative to the original crosshair of the
#    HEALPixSkyInstrument (in degrees). Corresponds to a horizontal rotation to the right. Default is 0.
#

# -----------------------------------------------------------------


def do(
    inputFileName: (
        str,
        "HEALPixSkyInstrument output image containing the raw HEALPix pixel data",
    ),
    outputFileName: (
        str,
        "output image containing an all sky projection of the input HEALPix grid",
    ),
    projection: (str, "projection to use (Mollweide/HammerAitoff)") = "Mollweide",
    nPixelY: (
        int,
        "number of pixels to use in the vertical direction for the output image",
    ) = 250,
    thetaCenter: (
        float,
        "zenith angle of the central position relative to the original crosshair (in degrees),"
        " corresponding to a vertical rotation that goes downward in the right half of the image",
    ) = 0.0,
    phiCenter: (
        float,
        "azimuth angle of the central position relative to the original crosshair (in degrees),"
        " corresponding to a horizontal rotation to the right",
    ) = 0.0,
) -> "create a projected image based on a HEALPixSkyInstrument output cube":

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

    if nPixelY <= 0:
        raise UserError("Invalid number of pixels: {0}!".format(nPixelY))

    inputFile = fits.open(inputPath)
    logging.info("Opened input file {0}".format(inputPath))

    hData = healpix.HEALPixGrid(inputFile[0].data)
    hData.printInfo()
    outputData = hData.getProjectionMap(nPixelY, projection, thetaCenter, phiCenter)

    outputHDUL = fits.HDUList(
        [fits.PrimaryHDU(outputData, header=inputFile[0].header), inputFile[1]]
    )

    outputHDUL.writeto(outputPath)
    logging.info("Wrote output file {0}".format(outputPath))
