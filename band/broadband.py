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
# Broadcast instances are usually obtained by loading built-in band definitions using one of the class functions
# builtinBands() or builtinBand() as opposed to directly invoking the constructor. However, for specific use cases,
# a Broadcast instance can also be constructed directly, for example to obtain a uniform filter with constant
# transmission over a specified wavelength range. For more information, refer to the documentation of the constructor.
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
# download the appropriate resource packs from the SKIRT web site, unzip the archives, and place the files in the
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

    ## As mentioned in the class header, built-in Broadcast instances can be obtained using one of the class functions
    # builtinBands() or builtinBand() as opposed to directly invoking the constructor. For more information, refer to
    # the doumentation of these functions.
    #
    # The constructor supports the following use cases, depending on the type of the \em bandspec argument:
    #
    #   - If \em bandspec is a string or a pathlib.Path, it should refer to a file in SKIRT stored table format
    #     that contains a properly normalized band definition (i.e., the same format as built-in band definitions).
    #     The BroadBand instance is constructed by loading the contents of this file.
    #
    #   - if \em bandspec is a tuple with two numbers, a bolometer-type band is constructed with a uniform
    #     transmission curve in the indicated (min,max) wavelength range expressed in micron.
    #
    #   - if \em bandspec is a two-dimensional numpy array, a custom bolometer-type band is constructed using
    #     \em bandspec[:,0] as wavelength grid and \em bandspec[:,1] as transmission curve. The wavelengths
    #     must be expressed in micron, while the transmission curve has arbitrary scaling.
    #     The array layout used here is identical to what one would obtain when applying numpy.loadtxt()
    #     to the input file for a SKIRT FileBand.
    #
    # If \em bandspec has a different type, or if the string/path does not refer to a file in the proper format,
    # the constructor raises an exception.
    #
    # \note As indicated above, the constructor expects the wavelengths specifying uniform and custom bands
    # as plain numbers expressed in micron. In contrast, the wavelengths and transmissions stored and returned
    # by this class are astropy quantities that include poper units. For uniform and custom bands, these units
    # will be micron (wavelengths) and 1/micron (transmissions). For built-in bands, these units will be
    # m (wavelengths) and 1/m (transmissions). Clients of this class should thus properly honor the astropy
    # units of any values returned by this class.
    #
    def __init__(self, bandspec):
        # construction of band from stored table
        if isinstance(bandspec, (str, pathlib.Path)):
            # get the corresponding band path
            bandpath = ut.absPath(bandspec)
            self._bandname = bandpath.stem
            if self._bandname.endswith("_BroadBand"):
                self._bandname = self._bandname[:-10]

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

        else: raise ValueError("Unsupported type of band specification '{}'".format(bandspec))

    ## This private function normalizes the transmission curve for uniform and custom bands,
    # and attaches the appropriate astropy units (micron and 1/micron, respectively).
    def _normalize(self):
        # normalize the transmission curve to unity
        self._transmissions /= np.trapz(x=self._wavelengths, y=self._transmissions)

        # attach the appropriate astropy units
        self._wavelengths <<= u.micron
        self._transmissions <<= u.micron**(-1)

    # -----------------------------------------------------------------

    ## This private function adds any band definitions detected in the SKIRT or PTS resources
    # directories to the list of built-in bands. It does so only the first time it is called.
    @classmethod
    def _ensureBuiltinBands(cls):
        if not cls._added:
            cls._addBuiltinBands(ut.skirtResourcesPath())
            cls._addBuiltinBands(ut.ptsResourcesPath())
            cls._added = True

    ## This private function recursively searches the contents of the specified directory and adds
    # files with a name ending in "_BroadBand.stab" to the list of built-in bands.
    @classmethod
    def _addBuiltinBands(cls, directory):
        if directory is not None:
            for path in pathlib.Path(directory).glob("**/*_BroadBand.stab"):
                if path.is_file():
                    cls._bandpaths.add(path)

    ## This function returns a list of BroadBand instances that match the specified criteria. The list is in
    # arbitrary order an can be empty.
    #
    # - \em nameSegments is a string specifying broadband name segments seperated by an underscore, a comma or a space.
    #   A built-in broad-band matches as soon as its name contains one or more of the specified segments.
    #   The comparison is case-insensitive. An empty string (the default) matches all bands.
    #   To specify a given band, it suffices to include two or more segments of its name that
    #   uniquely identify the band, seperated by an underscore or a space. For example, to select the
    #   HERSCHEL_PACS_100 band, one could enter "Herschel 100", "PACS 100", or "HERSCHEL_PACS_100".
    #
    # - \em minWavelength and \em maxWavelength are astropy quantities specifying a wavelength range.
    #   A built-in broad-band matches as soon as its pivot wavelength is inside the specified range.
    #   The default values specify an unlimited wavelength range.
    #
    @classmethod
    def builtinBands(cls, nameSegments="", minWavelength=0<<u.m, maxWavelength=1e99<<u.m):
        cls._ensureBuiltinBands()

        # construct a list with all bands that match the specified name segments
        result = []
        specsegments = nameSegments.upper().replace("_"," ").replace(","," ").split()
        for bandpath in cls._bandpaths:
            bandsegments = bandpath.stem[:-10].upper().replace("_"," ").split()
            if all([ (specsegment in bandsegments) for specsegment in specsegments ]):
                result.append(BroadBand(bandpath))

        # remove bands outside of the specified wavelength range
        result = [ band for band in result if minWavelength <= band.pivotWavelength() <= maxWavelength ]
        return result

    ## This function returns the single BroadBand instance that matches the specified criteria. It raises an error
    # if multiple bands or no band match the criteria. The arguments are the same as those for builtinBands().
    @classmethod
    def builtinBand(cls, nameSegments="", minWavelength=0<<u.m, maxWavelength=1e99<<u.m):
        result = cls.builtinBands(nameSegments, minWavelength, maxWavelength)
        if len(result) > 1:
            raise ValueError("Name segments '{}' match multiple band names, "
                             "including {} and {}".format(nameSegments, result[0].name(), result[1].name()))
        if len(result) == 0:
            raise ValueError("Band specification matches no bands")
        return result[0]

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

## This private helper function returns the natural logarithm for positive values, and a large negative number
# (but not infinity) for zero or negative values.
def _log(X):
    zeromask = X <= 0
    logX = np.empty(X.shape)
    logX[zeromask] = -750.  # the smallest (in magnitude) negative value x for which np.exp(x) returns zero
    logX[~zeromask] = np.log(X[~zeromask])
    return logX

# -----------------------------------------------------------------

## This function returns a list of BroadBand instances that match the specified criteria. It merely calls the
# BroadBand.builtinBands() function; see the documentation of that function for more information.
def builtinBands(nameSegments="", minWavelength=0<<u.m, maxWavelength=1e99<<u.m):
    return BroadBand.builtinBands(nameSegments, minWavelength, maxWavelength)

## This function returns the single BroadBand instance that matches the specified criteria. It merely calls the
# BroadBand.builtinBand() function; see the documentation of that function for more information.
def builtinBand(nameSegments="", minWavelength=0<<u.m, maxWavelength=1e99<<u.m):
    return BroadBand.builtinBand(nameSegments, minWavelength, maxWavelength)

# -----------------------------------------------------------------
