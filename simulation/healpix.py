#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.simulation.healpix Handling of HEALPix grids
#
# This module offers support for projecting HEALPix grids output by the HEALPixSkyInstrument.
#

# -----------------------------------------------------------------

import numpy as np
from pts.utils.error import UserError

# -----------------------------------------------------------------


## This function returns the index arrays corresponding to the given input zenith and azimuth angles for a
# HEALPix grid in SKIRT's modified ring order with the given Nside parameter.
#
# For each pair of input angles, the returned pair of output indices points to the pixel that contains
# those angles. Using those indices to access the raw HEALPix data will return the corresponding pixel
# value.
#
# Note that the code here is (and always should be) identical to the implementation used internally by
# SKIRT to guarantee consistency. We do however make use of NumPy arrays and NumPy array indexing for
# increased efficiency.
#
def getHEALPixIndices(theta, phi, Nside):
    z = np.cos(theta)
    za = np.abs(z)
    tt = np.mod(2.0 * phi / np.pi, 4.0)

    # empty ring (j) and pixel-in-ring (i) index arrays
    j = np.zeros(za.shape, dtype=int)
    i = np.zeros(za.shape, dtype=int)

    # we need different expressions for HEALPix pixels in the equatorial plane
    # (za <= 2./3.) and in the polar caps (za > 2./3.)
    # to exploit NumPy efficiency, we use index arrays

    # equatorial plane
    ieq = za <= 2.0 / 3.0
    t1 = Nside * (0.5 + tt[ieq])
    t2 = 0.75 * Nside * z[ieq]
    jp = np.floor(t1 - t2)
    jm = np.floor(t1 + t2)
    j[ieq] = Nside + 1 + jp - jm
    kshift = 1 - np.bitwise_and(j[ieq], 1)
    t1 = jp + jm + kshift + 1 + 7 * Nside
    i[ieq] = np.bitwise_and(np.right_shift(t1.astype(int), 1), 4 * Nside - 1)
    j[ieq] += Nside - 2

    # polar caps
    ipol = ~ieq
    ntt = np.floor(tt[ipol])
    tp = tt[ipol] - ntt
    tmp = np.zeros(tp.shape)
    # deal with small values of sin(theta)
    tmp[za[ipol] >= 0.99] = (
        Nside
        * np.sin(theta[ipol & (za >= 0.99)])
        / np.sqrt((1.0 + za[ipol & (za >= 0.99)]) / 3.0)
    )
    tmp[za[ipol] < 0.99] = Nside * np.sqrt(3.0 * (1 - za[ipol & (za < 0.99)]))
    jp = np.floor(tp * tmp)
    jm = np.floor((1.0 - tp) * tmp)
    j[ipol] = jp + jm + 1
    i[ipol] = np.floor(tt[ipol] * j[ipol])
    # distinguish between north and south pole
    j[ipol & (z > 0)] -= 1
    j[ipol & (z < 0)] = 4 * Nside - j[ipol & (z < 0)] - 1

    return j, i


## This function returns the projection of the given HEALPix data map, using the given
# number of vertical pixels for the resulting projection image, and the given projection
# transformation. The optional parameters thetaCenter and phiCenter allow to select a
# different central position than the original crosshair of the HEALPixSkyInstrument.
#
# The function first sets up the linear coordinates for the projection image and then converts
# them to the appropriate zenith and azimuth angles using the appropriate projection. These angles
# are then rotated according to thetaCenter and phiCenter so that they point to the correct angles
# on the rotated HEALPix sphere. Finally, the angles are passed on to getHEALPixIndices() to get
# the corresponding HEALPix pixels.
#
# Because of the nature of the Mollweide projection, some pixels in the four corners of the output
# image are unused. The values in these pixels are set to zero.
#
# The angles thetaCenter and phiCenter lead to vertical and horizontal rotations respectively. The
# horizontal rotation is simply to the right. The vertical rotation goes upward in the left part of
# the image and downward in the right part. Pixels that rotate through the pole move from the left
# part to the right part. Due to the projection, thetaCenter can go through a full 360 degrees
# rotation before the projection image is the same; thetaCenter angles larger than 180 degrees
# correspond to projection images for smaller thetaCenter values that are seen upside down.
#
def getProjectionMap(
    HEALPixCube, nPixelY, projection="Mollweide", thetaCenter=0.0, phiCenter=0.0
):
    # derive the HEALPix parameters from the image size
    Nside = HEALPixCube.shape[1] // 4

    # convert the azimuth angle from degrees to radians
    phiCenter *= np.pi / 180.0
    thetaCenter *= np.pi / 180.0

    # make sure we oversample the HEALPix pixels initially, we will resample after the projection
    if nPixelY < 2 * Nside:
        nPixelYHigh = 2 * Nside
        nPixelYHigh += nPixelY - nPixelYHigh % nPixelY
        resFac = nPixelYHigh // nPixelY
    else:
        nPixelYHigh = nPixelY
        resFac = 1

    # initialize the image pixels
    nPixelXHigh = 2 * nPixelYHigh
    yMin = -1.0 + 1.0 / nPixelYHigh
    dY = 2.0 / nPixelYHigh
    xMin = -2.0 + 2.0 / nPixelXHigh
    dX = 4.0 / nPixelXHigh
    y, x = np.mgrid[yMin:1.0:dY, xMin:2.0:dX]

    # create an empty pixel map
    image = np.zeros(x.shape)

    if projection == "Mollweide":
        # compute the Mollweide angles for each pixel
        temp = np.arcsin(y)
        theta = np.arcsin((2.0 * temp + np.sin(2.0 * temp)) / np.pi) + 0.5 * np.pi
        phi = np.pi + 0.5 * np.pi * x / np.cos(temp)

        # filter out invalid angles in the corners
        iValid = (0.5 * x) ** 2 + (y) ** 2 < 1.0
        theta = theta[iValid]
        phi = phi[iValid]
    elif projection == "HammerAitoff":
        x *= np.sqrt(2.0)
        y *= np.sqrt(2.0)
        # compute the Hammer-Aitoff angles for each pixel
        temp = np.sqrt(1.0 - (0.25 * x) ** 2 - (0.5 * y) ** 2)
        theta = np.arcsin(temp * y) + 0.5 * np.pi
        phi = 2.0 * np.arctan(0.5 * temp * x / (2.0 * temp ** 2 - 1.0)) + np.pi

        # filter out invalid angles in the corners
        iValid = (0.5 * x) ** 2 + (y) ** 2 < 2.0
        theta = theta[iValid]
        phi = phi[iValid]
    else:
        raise UserError(
            "Unknown projection ("
            + projection
            + ")! Possible values are Mollweide, HammerAitoff."
        )

    # invert the azimuth angle to account for the SKIRT orientation convention
    phi = 2.0 * np.pi - phi

    # we have successfully applied the projection to obtain the angular coordinates of each
    # pixel of the projection image on the unit sphere, centred on (theta, phi) = (0, 0)
    # now we need to transform these angles to the unit sphere centred on (thetaCenter, phiCenter)
    # we do this by applying two rotations:
    #  - a rotation over phiCenter along the z axis. This rotation only affects phi
    phi += phiCenter
    #  - a rotation over thetaCenter along the (rotated) x axis. This rotation needs to be done
    #    in Cartesian space
    sinTheta = np.sin(theta)
    x = sinTheta * np.cos(phi)
    y = sinTheta * np.sin(phi)
    z = np.cos(theta)
    sinThetaCenter = np.sin(thetaCenter)
    cosThetaCenter = np.cos(thetaCenter)
    temp = y * cosThetaCenter + z * sinThetaCenter
    z = -y * sinThetaCenter + z * cosThetaCenter
    y = temp
    phi = np.arctan2(y, x)
    theta = np.arccos(z)

    # make sure the final angles are within range for the HEALPix grid
    phi = np.mod(phi, 2.0 * np.pi)
    theta = np.mod(theta, np.pi)

    j, i = getHEALPixIndices(theta, phi, Nside)

    # copy the HEALPix pixel values into the corresponding image pixels (only for valid pixels)
    image[iValid] = HEALPixCube[j, i]

    # resample the image onto the desired resolution
    if resFac > 1:
        temp = image.reshape(
            (image.shape[0] // resFac, resFac, image.shape[1] // resFac, resFac)
        )
        image = np.sum(temp, axis=(1, 3))

    return image


# -----------------------------------------------------------------
