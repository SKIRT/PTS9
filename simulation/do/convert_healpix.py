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
) -> "create a projected image based on a HEALPixSkyInstrument output cube":

    import numpy as np
    import astropy.io.fits as fits
    import pts.simulation.healpix as healpix
    import pathlib
    from pts.utils.error import UserError

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

    inputData = inputFile[0].data
    numWav = len(inputData)
    outputData = np.zeros((numWav, nPixelY, 2 * nPixelY))
    for i in range(numWav):
        outputData[i] = healpix.getProjectionMap(inputData[i], nPixelY, projection)

    outputHDUL = fits.HDUList(
        [fits.PrimaryHDU(outputData, header=inputFile[0].header), inputFile[1]]
    )

    outputHDUL.writeto(outputPath)
