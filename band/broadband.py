#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.band.broadband Class for representing broadband filters
#
# This module contains the BroadBand class. An instance of this class represents
# a broadband filter with a given transmission curve.

# -----------------------------------------------------------------

import pathlib
import numpy as np
import astropy.units as u
import pts.utils as ut
import pts.simulation as sm
import pts.storedtable as stab

# -----------------------------------------------------------------

## An instance of the BroadBand class represents a broadband filter with a given transmission curve.
# An instance can be constructed from one of the "built-in" band definitions (see below),
# as a uniform filter with constant transmission over a specified wavelength range,
# or by providing a custom wavelength grid and corresponding transmission curve.
#
# Refer to the documentation of the constructor for this class for information on how to construct
# one of these types of BroadBand objects.
#
# Built-in bands
# --------------
#
# This class relies on SKIRT resource files to provide built-in broadband definitions. These resource files are
# written in the SKIRT stored table format and have names that end with "_BroadBand.stab". The first part of
# the file name identifies the particular band being represented.
#
# The class recursively searches for broadband resource files in the contents of two directories, if they exist:
# the ~/SKIRT/resources directory and the ~/PTS/resources directory. In other words, if SKIRT has been installed
# and all relevant resource packs have been installed, all broadbands available to SKIRT are available to this class
# as well. In case SKIRT is not installed on the same system (and in the same project directory), one can manually
# download the aprpopriate resource packs from the SKIRT web site, unzip the archives, and place the files in the
# ~/PTS/resources directory.
#
# For a list of available broadbands, see the documentation of the SKIRT %BroadBand class, or use the
# "pts list_bands" command.
#
# Band properties and operations
# ------------------------------
#
# This class offers operations such as determining the band's pivot wavelength and convolving an SED
# or a frame data cube with the band's transmission curve.
#
# Refer to the appendix in Camps et al. 2016 (MNRAS 462, 1057-1075) for a brief introduction of
# the relevant concepts and a derivation of the corresponding formulas. For the purposes of this
# class, a band is defined through its arbitrarily scaled transmission curve \f$T(\lambda)\f$.
# The relative transmission at a given wavelength can then be written as \f[
# T_\text{rel}(\lambda) = \frac{ T(\lambda) } { \max[ T(\lambda) ] }. \f]
#
# Given a spectral energy distribution \f$L_\lambda(\lambda)\f$, the mean specific luminosity
# over the band \f$\left<L_\lambda\right>\f$ can then be obtained through \f[
# \left<L_\lambda\right> = \frac{ \int L_\lambda(\lambda)T(\lambda) \,\mathrm{d}\lambda } { \int
# T(\lambda) \,\mathrm{d}\lambda }. \f]
#
# The pivot wavelength of the band is defined as the wavelength at which the mean specific
# luminosity can be properly converted between wavelength and frequency representations. It is
# given by \f[ \lambda_\mathrm{pivot} = \sqrt{ \frac{ \int T(\lambda) \,\mathrm{d}\lambda } {
# \int T(\lambda) \,\mathrm{d}\lambda/\lambda^2 } }. \f]
#
# The effective width of the band is defined as the horizontal size of a rectangle with height
# equal to the maximum transmission and with the same area as the one covered by the band's
# transmission curve. It is given by \f[ W_\mathrm{eff} = \frac{ \int T(\lambda)
# \,\mathrm{d}\lambda } { \max[ T(\lambda) ] }. \f]
#
# As set forth by Camps et al. 2016, for energy measuring devices (bolometers) the total system
# transmission \f$T(\lambda)\f$ is usually given by the instrument designers, and it can be used
# directly for this class. For photon counters (including all instruments in the UV, optical and
# near infrared), the total system response \f$R(\lambda)\f$ is usually given instead. The
# transmission curve needed for this class can be derived from the response curve through \f[
# T(\lambda) = \lambda\,R(\lambda). \f]
#
# Also, to further simplify the above formulas, the constructor normalizes the transmission curve
# to unity, i.e. \f[ \int T(\lambda) \,\mathrm{d}\lambda = 1. \f]
#
class BroadBand:
    ## The file paths for all detected built-in band definitions (i.e. ending with "_BroadBand.stab")
    _bandpaths = set()

    ## Flag becomes True as soon as bandpaths have been added for the SKIRT and PTS resource directories
    _added = False

    ## The constructor creates a BroadBand instance in one of the following three ways, depending on the type of
    # the \em bandspec argument:
    #
    #   - if \em bandspec is a string, it looks for a built-in broadband filter based on the text segments
    #     in the string, and then loads the transmission curve from the corresponding resource file.
    #     To specify a band, it suffices to include two or more segments of its name (case insensitive) that
    #     uniquely identify the band, seperated by an underscore or a space. For example, to select the
    #     HERSCHEL_PACS_100 band, one could enter "Herschel 100", "PACS 100", or "HERSCHEL_PACS_100".
    #
    #   - if \em bandspec is a tuple with two numbers, a bolometer-type band is constructed with a uniform
    #     transmission curve in the indicated (min,max) wavelength range (expressed in micron).
    #
    #   - if \em bandspec is a two-dimensional numpy.ndarray, a bolometer-type band is constructed using
    #     \em bandspec[:,0] as wavelength grid and \em bandspec[:,1] as transmission curve. The wavelengths are
    #     assumed to be expressed in micron, while the transmission curve has arbitrary scaling.
    #     The array layout used here is identical to what one would obtain when applying numpy.loadtxt() to the
    #     input file for a SKIRT FileBand.
    #
    # If \em bandspec has a different type, or if its string contents does not unambiguously match a built-in band
    # name, the constructor raises an exception.
    #
    # \note During the first part of construction, wavelengths are initialized in micron but without explicit units,
    # and transmissions have arbitrary scale. Thus, these are plain numpy arrays, not astropy quantities. At the
    # end of construction, the appropriate astropy units are attached.
    #
    def __init__(self, bandspec):
        # construction of built-in band
        if isinstance(bandspec, str):
            self._ensureBuiltinBands()

            # get the corresponding built-in band path
            bandpath = self._bandPathFromSpec(bandspec)
            self._bandname = bandpath.stem[:-10]

            # load wavelengths and normalized transmissions from the stored table, including proper astropy units
            table = stab.readStoredTable(bandpath)
            self._wavelengths = table['lambda']
            self._transmissions = table['T']

        # construction of uniform band
        elif isinstance(bandspec, tuple) and len(bandspec)==2:
            self._bandname = "Uniform"
            self._wavelengths = np.array(list(map(float, bandspec)))
            self._transmissions = np.ones(2)
            self._normalize()

        # construction of a custom band
        elif isinstance(bandspec, np.ndarray) and len(bandspec.shape)==2:
            self._bandname = "Custom"
            self._wavelengths = bandspec[:,0].astype(float)
            self._transmissions = bandspec[:,1].astype(float)
            self._normalize()

        else: raise ValueError("Unsuppported type of band specification '{}'".format(bandspec))

    ## This function matches the given band specification with all detected built-in band paths.
    # If there is a single match, the corresponding path is returned. If there is no match,
    # or if there are multiple matches, an exception is raised.
    # A band specification must include two or more segments of the band name (case insensitive),
    # separated by whitespace and/or underscores.
    def _bandPathFromSpec(self, bandspec):
        specsegments = bandspec.upper().replace("_"," ").split()
        result = None
        for bandpath in self._bandpaths:
            namesegments = bandpath.stem[:-10].upper().replace("_"," ").split()
            if all([ (specsegment in namesegments) for specsegment in specsegments ]):
                if result is not None:
                    raise ValueError("Band specification '{}' matches multiple band names, "
                                     "including {} and {}".format(bandspec, result.stem[:-10], bandpath.stem[:-10]))
                result = bandpath
        if result is None:
            raise ValueError("Band specification '{}' matches no band names".format(bandspec))
        return result

    ## This function normalizes the transmission curve after it has been loaded during construction,
    # and attaches the appropriate astropy units (micron and 1/micron, respectively).
    def _normalize(self):
        # normalize the transmission curve to unity
        self._transmissions /= np.trapz(x=self._wavelengths, y=self._transmissions)

        # attach the appropriate astropy units
        self._wavelengths <<= u.micron
        self._transmissions <<= u.micron**(-1)

    # -----------------------------------------------------------------

    ## This function adds any band definitions detected in the SKIRT or PTS resources
    # directories to the list of built-in bands. It does so only the first time it is called.
    @classmethod
    def _ensureBuiltinBands(cls):
        if not cls._added:
            cls._addBuiltinBands(ut.skirtResourcesPath())
            cls._addBuiltinBands(ut.ptsResourcesPath())
            cls._added = True

    ## This function recursively searches the contents of the specified directory and adds
    # files with a name ending in "_BroadBand.stab" to the list of built-in bands.
    @classmethod
    def _addBuiltinBands(cls, directory):
        if directory is not None:
            for path in pathlib.Path(directory).glob("**/*_BroadBand.stab"):
                if path.is_file():
                    cls._bandpaths.add(path)

    ## This function returns an iterable over all built-in band names, in arbitrary order.
    @classmethod
    def builtinBandNames(cls):
        cls._ensureBuiltinBands()
        return [bandpath.stem[:-10] for bandpath in cls._bandpaths]

    # -----------------------------------------------------------------

    ## This function returns the built-in band's name, or "Uniform" or "Custom"
    # for uniform or custom bands respectively.
    def name(self):
        return self._bandname

    ## This function returns a tuple of astropy quantities specifying the wavelength range for this band.
    # The transmission may be zero in some (usually small) intervals within the range, but it is guaranteed to be zero
    # outside the range.
    def wavelengthRange(self):
        return self._wavelengths[0], self._wavelengths[-1]

    ## This function returns the minimum wavelength for the band, as an astropy quantity.
    def minWavelength(self):
        return self._wavelengths[0]

    ## This function returns the maximum wavelength for the band, as an astropy quantity.
    def maxWavelength(self):
        return self._wavelengths[-1]

    ## This function returns the pivot wavelength as defined in the class header, as an astropy quantity.
    #
    # \note The implementation here does purposefully not use a call to numpy.trapz() to exactly reproduce the
    # numerical quadrature rule implemented in the Band class in SKIRT. This to avoid significant differences in
    # pivot wavelength for a uniform band consisting of only two wavelength points. Both the implementation here
    # and in SKIRT are bad approximations for the real pivot wavelength in this case, but at least they will be
    # the same, which is important when identifying bands in SKIRT output.
    def pivotWavelength(self):
        dlam = self._wavelengths[1:] - self._wavelengths[:-1]
        lam = 0.5 * (self._wavelengths[1:] + self._wavelengths[:-1])
        T = 0.5 * (self._transmissions[1:] + self._transmissions[:-1])
        return 1. / np.sqrt((T*dlam/lam**2).sum())

    ## This function returns the effective band width (a wavelength interval) as defined in the class header
    # as an astropy quantity.
    def effectiveWidth(self):
        return 1. / self._transmissions.max()

    ## This function returns a tuple containing a copy of the wavelength and transmission arrays representing
    # the transmission curve for this band, normalized to unity. Both arrays are astropy quantities.
    def transmissionCurve(self):
        return self._wavelengths.copy(), self._transmissions.copy()

    # -----------------------------------------------------------------

    ## This function calculates and returns the band-averaged value \f$\left<F_\lambda\right>\f$ for a given
    # spectral energy distribution \f$F_\lambda(\lambda)\f$. For a band with transmission curve \f$T(\lambda)\f$,
    # \f[ \left<F_\lambda\right> = \frac{ \int F_\lambda(\lambda)T(\lambda) \,\mathrm{d}\lambda }
    # { \int T(\lambda) \,\mathrm{d}\lambda }. \f]
    #
    # The function in fact accepts a distribution (over a range of wavelengths) of various spectral quantities,
    # including flux density, surface brightness, spectral radiance, or spectral luminosity of any "flavor"
    # (neutral, per wavelength, or per frequency) and in arbitrary units. For the purposes of this function,
    # these quantities are generically referred to as "flux". The incoming fluxes are converted to an equivalent
    # "per-wavelength" flavor, the convolution is calculated according to the equation above, and the result
    # is converted back to the flavor and units of the incoming fluxes, or to the optionally specified flavor
    # and/or units.
    #
    # The function accepts the following arguments:
    # - \em wavelengths: an astropy quantity array specifying the wavelengths \f$\lambda_\ell\f$, in increasing order,
    #   on which the fluxes are sampled. The convolution is performed on a wavelength grid that combines the grid
    #   points given here with the grid points on which the band transmission curve is defined.
    # - \em fluxes: an astropy quantity array specifying the flux distribution(s) \f$F_\lambda(\lambda_\ell)\f$.
    #   This can be an array with the same length as \em wavelengths, or a multi-dimensional
    #   array where the last dimension has the same length as \em wavelengths.
    #   The returned result will have the shape of \em fluxes minus the last (or only) dimension.
    # - \em numWavelengths: an integer specifying the approximate number of wavelength grid points used in the
    #   convolution calculation, or None (or omitted) to indicate no limit. The incoming flux distribution is resampled
    #   to all wavelength points in the combined convolution grid. Convolving large data cubes may thus consume
    #   a prohibitive amount of memory and/or computation time unless the number of wavelengths is limited.
    # - \em flavor: the flavor and/or units to which the result must be converted. If missing or None, the result
    #   is converted back to the flavor and units of the incoming fluxes. If specified, this can be one of the
    #   strings "neutral", "wavelength", or "frequency", or an explicit astropy unit instance that is compatible
    #   with the incoming flux type (density, brightness, luminosity), after flavor conversion.
    #
    def convolve(self, wavelengths, fluxes, *, numWavelengths=None, flavor=None):

        # convert fluxes to per-wavelength flavor
        if flavor is None: flavor = fluxes.unit
        fluxes = sm.convertToFlavor(wavelengths, fluxes, "wavelength")

        # get the involved wavelength grids, in micron, without units
        wa = wavelengths.to_value(u.micron)
        wb = self._wavelengths.to_value(u.micron)

        # create a combined wavelength grid (without units), restricted to the overlapping interval
        w1 = wa[ (wa>=wb[0]) & (wa<=wb[-1]) ]
        w2 = wb[ (wb>=wa[0]) & (wb<=wa[-1]) ]
        w = np.unique(np.hstack((w1, w2)))
        if len(w) < 2: return 0

        # downsample the grid if so requested
        if numWavelengths is not None:
            numWavelengths = max(2, numWavelengths)
            if len(w) > numWavelengths:
                downfactor = len(w) // numWavelengths
                w = w[::downfactor]

        # log-log interpolate transmission and fluxes (without units) on the combined wavelength grid
        T = np.exp(np.interp(np.log(w), np.log(wb), _log(self._transmissions.to_value(1/u.micron)), left=0., right=0.))
        # use SciPy: NumPy interpolation doesn't support 1D interpolation of multi-dimensional arrays
        from scipy.interpolate import interp1d
        F = np.exp(interp1d(np.log(wa), _log(fluxes.value), copy=False, bounds_error=False, fill_value=0.)(np.log(w)))

        # perform integration, re-assign stripped per-wavelength units and convert back to original units
        convolved = np.trapz(x=w, y=F*T) << fluxes.unit
        return sm.convertToFlavor(self.pivotWavelength(), convolved, flavor)

# -----------------------------------------------------------------

## This helper function returns the natural logarithm for positive values, and a large negative number
# (but not infinity) for zero or negative values.
def _log(X):
    zeromask = X <= 0
    logX = np.empty(X.shape)
    logX[zeromask] = -750.  # the smallest (in magnitude) negative value x for which np.exp(x) returns zero
    logX[~zeromask] = np.log(X[~zeromask])
    return logX

# -----------------------------------------------------------------

## This function returns an iterable over all built-in band names, in arbitrary order.
def builtinBandNames():
    return BroadBand.builtinBandNames()

# -----------------------------------------------------------------
