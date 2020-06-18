#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.simulation.healpix Handling of HEALPix grids
#
# This module offers functionality to work with HEALPix grids output by the HEALPixSkyInstrument.
#

# -----------------------------------------------------------------

import numpy as np
from pts.utils.error import UserError
import logging

## This class represents a HEALPix grid. It contains a reference to the raw
# HEALPixSkyInstrument data cube and some additional geometrical information
# about the grid that is derived from it.
#
# It offer functionality to access, manipulate and project the HEALPix pixels.
#
class HEALPixGrid:

    ## Constructor.
    #
    # Responsible for initialising the geometrical information about the HEALPix grid.
    #
    def __init__(self, HEALPixCube):
        self._HEALPixCube = HEALPixCube
        self._NSide = HEALPixCube.shape[-1] // 4
        self._order = int(np.log2(self._NSide))

    ## Print some information about the HEALPix grid to the console.
    def printInfo(self):
        logging.info(
            "HEALPix grid has side length {0} and order {1}.".format(
                self._NSide, self._order
            )
        )
        if len(self._HEALPixCube) == 3:
            logging.info(
                "It contains {0} variables per pixel.".format(
                    self._HEALPixCube.shape[0]
                )
            )

    ## This function returns the angle arrays corresponding to the given input index arrays.
    #
    # The index array j contains the ring indices, the index array i contains the pixel-in-ring
    # indices. The return arrays contain the centers of the corresponding pixels.
    #
    # This function is the inverse of getHEALPixIndices(). No checks are done on the input values
    # to ensure that they are actually valid indices.
    #
    def getHEALPixAngles(self, j, i):

        # make sure the input values are NumPy arrays
        if not j is np.ndarray:
            j = np.array(j)
        if not i is np.ndarray:
            i = np.array(i)

        # reserve empty result arrays
        theta = np.zeros(j.shape)
        phi = np.zeros(i.shape)

        # cache the total number of pixels that is used below
        Ntot = 12 * self._NSide ** 2

        # now compute the angles
        # we treat the poles and equatorial region separately

        # north pole
        iNPol = j < self._NSide - 1
        z = 1.0 - 4.0 * (j[iNPol] + 1.0) ** 2 / Ntot
        theta[iNPol] = np.arccos(z)
        phi[iNPol] = 0.5 * np.pi * (i[iNPol] + 0.5) / (j[iNPol] + 1.0)

        # equator
        iEq = (j >= self._NSide - 1) & (j < 3 * self._NSide)
        z = 2.0 * (2 * self._NSide - j[iEq] - 1.0) / (3 * self._NSide)
        theta[iEq] = np.arccos(z)
        phi[iEq] = (
            (i[iEq] + 0.5 * np.bitwise_and(j[iEq], 1)) * 1.5 * np.pi / (3 * self._NSide)
        )

        # south pole
        iSPol = j >= 3 * self._NSide
        z = 4.0 * (4 * self._NSide - j[iSPol] - 1.0) ** 2 / Ntot - 1.0
        theta[iSPol] = np.arccos(z)
        phi[iSPol] = 0.5 * np.pi * (i[iSPol] + 0.5) / (4 * self._NSide - j[iSPol] - 1.0)

        return theta, phi

    ## This function returns the index arrays corresponding to the given input zenith and azimuth angles for the
    # HEALPix grid in SKIRT's modified ring order.
    #
    # For each pair of input angles, the returned pair of output indices points to the pixel that contains
    # those angles. Using those indices to access the raw HEALPix data will return the corresponding pixel
    # value.
    #
    # Note that the code here is (and always should be) identical to the implementation used internally by
    # SKIRT to guarantee consistency. We do however make use of NumPy arrays and NumPy array indexing for
    # increased efficiency.
    #
    def getHEALPixIndices(self, theta, phi):

        # make sure the input angles are NumPy arrays
        if not theta is np.ndarray:
            theta = np.array(theta)
        if not phi is np.ndarray:
            phi = np.array(phi)

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
        t1 = self._NSide * (0.5 + tt[ieq])
        t2 = 0.75 * self._NSide * z[ieq]
        jp = np.floor(t1 - t2)
        jm = np.floor(t1 + t2)
        j[ieq] = self._NSide + 1 + jp - jm
        kshift = 1 - np.bitwise_and(j[ieq], 1)
        t1 = jp + jm + kshift + 1 + 7 * self._NSide
        i[ieq] = np.bitwise_and(np.right_shift(t1.astype(int), 1), 4 * self._NSide - 1)
        j[ieq] += self._NSide - 2

        # polar caps
        ipol = ~ieq
        ntt = np.floor(tt[ipol])
        tp = tt[ipol] - ntt
        tmp = np.zeros(tp.shape)
        # deal with small values of sin(theta)
        tmp[za[ipol] >= 0.99] = (
            self._NSide
            * np.sin(theta[ipol & (za >= 0.99)])
            / np.sqrt((1.0 + za[ipol & (za >= 0.99)]) / 3.0)
        )
        tmp[za[ipol] < 0.99] = self._NSide * np.sqrt(3.0 * (1 - za[ipol & (za < 0.99)]))
        jp = np.floor(tp * tmp)
        jm = np.floor((1.0 - tp) * tmp)
        j[ipol] = jp + jm + 1
        i[ipol] = np.floor(tt[ipol] * j[ipol])
        # distinguish between north and south pole
        j[ipol & (z > 0)] -= 1
        j[ipol & (z < 0)] = 4 * self._NSide - j[ipol & (z < 0)] - 1

        return j, i

    ## This function returns the projection of the HEALPix data map, using the given
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
        self, nPixelY, projection="Mollweide", thetaCenter=0.0, phiCenter=0.0
    ):

        # convert the azimuth angle from degrees to radians
        phiCenter *= np.pi / 180.0
        thetaCenter *= np.pi / 180.0

        # make sure we oversample the HEALPix pixels initially, we will resample after the projection
        if nPixelY < 2 * self._NSide:
            nPixelYHigh = 2 * self._NSide
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
        if len(self._HEALPixCube.shape) == 2:
            image = np.zeros(x.shape)
        elif len(self._HEALPixCube.shape) == 3:
            image = np.zeros((self._HEALPixCube.shape[0], x.shape[0], x.shape[1]))
        else:
            raise ValueError(
                "Expected HEALPixCube of dimension 2 or 3, but instead got cube with shape {0}".format(
                    self._HEALPixCube.shape
                )
            )

        logging.info(
            "Creating a projected data cube with shape {0}".format(image.shape)
        )

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

        j, i = self.getHEALPixIndices(theta, phi)

        # copy the HEALPix pixel values into the corresponding image pixels (only for valid pixels)
        if len(self._HEALPixCube.shape) == 2:
            image[iValid] = self._HEALPixCube[j, i]
        elif len(self._HEALPixCube.shape) == 3:
            image[:, iValid] = self._HEALPixCube[:, j, i]

        # resample the image onto the desired resolution
        if resFac > 1:
            if len(self._HEALPixCube.shape) == 2:
                temp = image.reshape(
                    (image.shape[0] // resFac, resFac, image.shape[1] // resFac, resFac)
                )
                image = np.sum(temp, axis=(1, 3))
            elif len(self._HEALPixCube.shape) == 3:
                temp = image.reshape(
                    (
                        image.shape[0],
                        image.shape[1] // resFac,
                        resFac,
                        image.shape[2] // resFac,
                        resFac,
                    )
                )
                image = np.sum(temp, axis=(2, 4))
            logging.info("Resampled projected cube to shape {0}".format(image.shape))

        return image

    ## This function returns a new HEALPixGrid that contains a degraded copy of this one.
    # The degraded copy combines 2*2^degradeFactor pixels into a single pixel using the
    # given operator (default is sum) to combine pixel values. This effectively changes
    # the order of the HEALPixGrid by -degradeFactor.
    #
    def degrade(self, degradeFactor, operator=np.sum):

        # compute the properties of the degraded grid
        newOrder = self._order - degradeFactor
        newNSide = 1 << newOrder

        # allocate the degraded data cube
        if len(self._HEALPixCube.shape) == 2:
            newHEALPixCube = np.zeros((4 * newNSide - 1, 4 * newNSide))
        elif len(self._HEALPixCube.shape) == 3:
            newHEALPixCube = np.zeros(
                (self._HEALPixCube.shape[0], 4 * newNSide - 1, 4 * newNSide)
            )

        # create the degraded grid
        newHEALPixGrid = HEALPixGrid(newHEALPixCube)

        # generate index arrays for the entire old HEALPix grid
        # we first overdo it a bit
        j, i = np.mgrid[0 : 4 * self._NSide - 1 : 1, 0 : 4 * self._NSide : 1]
        # now mask out the top and bottom triangles for the missing pixels in the polar regions
        select = np.ones(j.shape, dtype=bool)
        select[i + 1 > 4 * (j + 1)] = False
        select[i + 1 > 4 * (4 * self._NSide - j - 1)] = False
        j = j[select]
        i = i[select]

        # compute the angles for the centres of the old pixels
        theta, phi = self.getHEALPixAngles(j, i)
        # now get the pixels that contain these angles
        newj, newi = newHEALPixGrid.getHEALPixIndices(theta, phi)
        # what we want to do at this point, is something like
        #  newHEALPixCube[newj, newi] += self._HEALPixCube[j, i]
        # (assuming our operator is a summation)
        # this does not work, and that is very sensible, since this is like
        # trying to add to the same variable with multiple threads in a shared
        # memory context (and for all we know, that could be exactly what happens
        # inside the NumPy library)
        # we can however bypass this issue by using our knowledge of the grid:
        # we know that every new pixel corresponds to exactly 2*2^degradeFactor
        # old pixels, so if we could do the summation 2*2^degradeFactor times,
        # and only sum the corresponding fraction of pixels in one go, then
        # we are done
        # this only works if the new pixels are ordered sensibly, which they
        # will not be in our RING pixel ordering (the first 4 indices for
        # example *always* correspond to different pixels, no matter what
        # the resolution of the HEALPix grid is)
        # we can however find a good ordering by computing the pixel index
        # for our new pixels, and argsorting on this. We don't even care
        # about properly accounting for the corners, as we are only
        # interested in the relative index ordering, not in retrieving the
        # correct RING pixel index
        # So:
        #  - compute the RING pixel indices within the degraded grid
        pixel = newj * 4 * newNSide + newi
        #  - argsort these. If we index the j and i arrays with the pixSort
        #    array, then 2*2^degradeFactor consecutive elements in the
        #    resulting arrays will correspond to the same pixel in the degraded
        #    grid
        pixSort = np.argsort(pixel)
        #  - now do the summations
        #    we could do this using a loop, but that would make it harder to
        #    apply our operator (not all operators can be applied sequentially)
        #    and would possibly be inefficient
        #    instead, we do some clever reshaping and directly put the new pixel
        #    values into their corresponding new pixels
        numOldPerNew = 2 * (1 << degradeFactor)
        # we do need the new pixel indices (in the right order) for this to work
        jsmall, ismall = np.mgrid[0 : 4 * newNSide - 1 : 1, 0 : 4 * newNSide : 1]
        select = np.ones(jsmall.shape, dtype=bool)
        select[ismall + 1 > 4 * (jsmall + 1)] = False
        select[ismall + 1 > 4 * (4 * newNSide - jsmall - 1)] = False
        jsmall = jsmall[select]
        ismall = ismall[select]
        # now get the old pixels sorted per new pixel index and reshaped
        # appropriately so that we can simply apply our operator to the last axis
        if len(self._HEALPixCube.shape) == 2:
            sortedPixels = self._HEALPixCube[j[pixSort], i[pixSort]].reshape(
                (12 * newNSide ** 2, numOldPerNew)
            )
            newHEALPixCube[jsmall[select], ismall[select]] = operator(
                sortedPixels, axis=1
            )
        elif len(self._HEALPixCube.shape) == 3:
            sortedPixels = self._HEALPixCube[:, j[pixSort], i[pixSort]].reshape(
                (self._HEALPixCube.shape[0], 12 * newNSide ** 2, numOldPerNew)
            )
            newHEALPixCube[:, jsmall, ismall] = operator(sortedPixels, axis=2)

        return newHEALPixGrid


# -----------------------------------------------------------------
