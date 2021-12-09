#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.band.broadband Class for representing broadband filters
#
# This module contains the BroadBand class. An instance of this class represents either a built-in
# "standard" broadband filter (including transmission curve data) or a uniform filter over a
# specified wavelength range (i.e. with a constant transmission curve in that range).
#

# -----------------------------------------------------------------

import astropy.constants as const
import astropy.units as u
import numpy as np
import xml.etree.ElementTree
import pts.utils as ut
import pts.simulation as sm

# -----------------------------------------------------------------

## An instance of the BroadBand class represents either a built-in standard broadband filter
# (including transmission curve data), a uniform filter over a specified wavelength range
# (with a constant transmission curve in that range), or a custom filter based on a given
# wavelength grid and transmission curve. Refer to the constructor of this class for
# information on how to construct one of these two types of BroadBand objects.
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
# T(\lambda) = \lambda\,R(\lambda). \f] When applicable, this operation is performed by the
# constructor just after loading the response curve data from file.
#
# Also, to further simplify the above formulas, the constructor normalizes the transmission curve
# to unity, i.e. \f[ \int T(\lambda) \,\mathrm{d}\lambda = 1. \f]
#
# Built-in bands
# --------------
#
# The table below lists the built-in bands available at the time of writing. The first column lists
# the complete band name; the second column indicates the corresponding pivot wavelength.
#
#   | Band name | Pivot wavelength (micron)
#   |-----------|--------------------------
#   | 2MASS_2MASS_J      | 1.2393
#   | 2MASS_2MASS_H      | 1.6495
#   | 2MASS_2MASS_KS     | 2.1639
#   | ALMA_ALMA_10       | 349.89
#   | ALMA_ALMA_9        | 456.2
#   | ALMA_ALMA_8        | 689.59
#   | ALMA_ALMA_7        | 937.98
#   | ALMA_ALMA_6        | 1244.4
#   | ALMA_ALMA_5        | 1616
#   | ALMA_ALMA_4        | 2100.2
#   | ALMA_ALMA_3        | 3043.4
#   | EUCLID_VIS_VIS     | 0.71032
#   | EUCLID_NISP_Y      | 1.0808
#   | EUCLID_NISP_J      | 1.3644
#   | EUCLID_NISP_H      | 1.7696
#   | GALEX_GALEX_FUV    | 0.15351
#   | GALEX_GALEX_NUV    | 0.23008
#   | GENERIC_JOHNSON_U  | 0.35247
#   | GENERIC_JOHNSON_B  | 0.44169
#   | GENERIC_JOHNSON_V  | 0.55251
#   | GENERIC_JOHNSON_R  | 0.6899
#   | GENERIC_JOHNSON_I  | 0.8739
#   | GENERIC_JOHNSON_J  | 1.2431
#   | GENERIC_JOHNSON_M  | 5.0122
#   | HERSCHEL_PACS_70   | 70.77
#   | HERSCHEL_PACS_100  | 100.8
#   | HERSCHEL_PACS_160  | 161.89
#   | HERSCHEL_SPIRE_250 | 252.55
#   | HERSCHEL_SPIRE_350 | 354.27
#   | HERSCHEL_SPIRE_500 | 515.35
#   | IRAS_IRAS_12       | 11.41
#   | IRAS_IRAS_25       | 23.609
#   | IRAS_IRAS_60       | 60.409
#   | IRAS_IRAS_100      | 101.14
#   | JCMT_SCUBA2_450    | 449.3
#   | JCMT_SCUBA2_850    | 853.81
#   | PLANCK_HFI_857     | 352.42
#   | PLANCK_HFI_545     | 545.55
#   | PLANCK_HFI_353     | 839.3
#   | PLANCK_HFI_217     | 1367.6
#   | PLANCK_HFI_143     | 2130.7
#   | PLANCK_HFI_100     | 3001.1
#   | PLANCK_LFI_70      | 4303
#   | PLANCK_LFI_44      | 6845.9
#   | PLANCK_LFI_30      | 10674
#   | RUBIN_LSST_U       | 0.368
#   | RUBIN_LSST_G       | 0.47823
#   | RUBIN_LSST_R       | 0.62178
#   | RUBIN_LSST_I       | 0.75323
#   | RUBIN_LSST_Z       | 0.86851
#   | RUBIN_LSST_Y       | 0.97301
#   | SLOAN_SDSS_U       | 0.35565
#   | SLOAN_SDSS_G       | 0.47025
#   | SLOAN_SDSS_R       | 0.61756
#   | SLOAN_SDSS_I       | 0.749
#   | SLOAN_SDSS_Z       | 0.89467
#   | SPITZER_IRAC_I1    | 3.5508
#   | SPITZER_IRAC_I2    | 4.496
#   | SPITZER_IRAC_I3    | 5.7245
#   | SPITZER_IRAC_I4    | 7.8842
#   | SPITZER_MIPS_24    | 23.759
#   | SPITZER_MIPS_70    | 71.985
#   | SPITZER_MIPS_160   | 156.43
#   | SWIFT_UVOT_UVW2    | 0.20552
#   | SWIFT_UVOT_UVM2    | 0.22462
#   | SWIFT_UVOT_UVW1    | 0.25804
#   | SWIFT_UVOT_U       | 0.34628
#   | SWIFT_UVOT_B       | 0.43496
#   | SWIFT_UVOT_V       | 0.54254
#   | TNG_OIG_U          | 0.37333
#   | TNG_OIG_B          | 0.43978
#   | TNG_OIG_V          | 0.53729
#   | TNG_OIG_R          | 0.63902
#   | TNG_NICS_J         | 1.2758
#   | TNG_NICS_H         | 1.6265
#   | TNG_NICS_K         | 2.2016
#   | UKIRT_UKIDSS_Z     | 0.88264
#   | UKIRT_UKIDSS_Y     | 1.0314
#   | UKIRT_UKIDSS_J     | 1.2501
#   | UKIRT_UKIDSS_H     | 1.6354
#   | UKIRT_UKIDSS_K     | 2.2058
#   | WISE_WISE_W1       | 3.3897
#   | WISE_WISE_W2       | 4.6406
#   | WISE_WISE_W3       | 12.568
#   | WISE_WISE_W4       | 22.314
#
# Built-in band implementation
# ----------------------------
#
# The hard-coded class dictionary \c _bandinfo maps the band name for all built-in bands to the information
# needed to load the corresponding transmission curve. The value tuple includes the following items:
#    - the directory name (relative to this package's \c data directory), which also indicates the file format.
#    - the name of the file from which to load the data, located in the specified directory.
#    - additional information, depending on the format, or None if the format does not need extra information.
#
# __SVO__
#
# The band descriptions in the SVO directory were downloaded from the SVO filter profile service web site
# http://svo2.cab.inta-csic.es/theory/fps in XML format.
#
# The last item of the value tuple in the \c _bandinfo dictionary is a Boolean flag indicating
# the type of instrument (because this information cannot be derived from the SVO file contents):
# the flag is True for photon counters and False for bolometers.
#
# __JCMT__
#
# The transmission curves for the SCUBA-2 450 and 850 micron instruments on the JCMT telescope were converted to
# a simple text column format from 'model450' and 'model850' data
# downloaded at http://www.eaobservatory.org/jcmt/instrumentation/continuum/scuba-2/filters/
#
# The JCMT instruments are bolometers.
#
# __PLANCK__
#
# The PLANCK transmission curves are provided as a text column file for each of the following bands:
#
# | Band name | Frequency (GHz)| Wavelength (micron) |
# |---------|-----|-------|
# | LFI 030 |  28 | 10600 |
# | LFI 044 |  44 |  6810 |
# | LFI 070 |  70 |  4260 |
# | HFI 100 | 100 |  3000 |
# | HFI 143 | 142 |  2100 |
# | HFI 217 | 222 |  1380 |
# | HFI 353 | 352 |   850 |
# | HFI 545 | 545 |   550 |
# | HFI 857 | 857 |   350 |
#
# The PLANCK instruments are bolometers.
#
# __ALMA__
#
# The ALMA transmission curves are provided as a single text column file for all bands. ALMA sensitivity depends
# on the amount of water in the atmosphere, quantified by the precipitable water vapour (pwv) in mm. However,
# because we use only transmission curves normalized in each band, the effect of the pwv is minimal for our purposes.
# Somewhat arbitrarily, we opted to use the data for pwv = 0.2 mm.
#
# The ALMA instruments are bolometers.
#
# The last item of the value tuple in the \c _bandinfo dictionary specifies the wavelength range for each ALMA band
# supported here, so that we can select the appropriate portion of the transmission curve (provdided in a single file).
#
# | ALMA BAND | Frequency range (GHz) | Wavelength range (mm) | Wavelength range (micron) |
# |----|-----------|-------------|-------------|
# |  1 |  31 - 45  | unsupported | unsupported |
# |  2 |  67 - 90  | unsupported | unsupported |
# |  3 |  84 - 116 | 3.57 - 2.59 | 2590 - 3570 |
# |  4 | 125 - 163 | 2.40 - 1.84 | 1840 - 2400 |
# |  5 | 162 - 211 | 1.84 - 1.42 | 1420 - 1840 |
# |  6 | 211 - 275 | 1.42 - 1.09 | 1090 - 1420 |
# |  7 | 275 - 373 | 1.09 - 0.80 |  800 - 1090 |
# |  8 | 385 - 500 | 0.78 - 0.60 |  600 - 780  |
# |  9 | 602 - 720 | 0.50 - 0.42 |  420 - 500  |
# | 10 | 787 - 950 | 0.38 - 0.32 |  320 - 380  |
#
#
# __EUCLID and RUBIN__
#
# The Euclid and Rubin response curves are provided as simple text column files with wavelength in Angstrom.
# Some files include a lot of leading and/or trailing wavelength points with zero or very small transmission;
# our load function removes these points.
#
# All Euclid and Rubin instruments are photon counters.
#
class BroadBand:
    ## The band information dictionary described in the class header
    _bandinfo = {
        # SVO
        "2MASS_2MASS_H": ("SVO", "2MASS.2MASS.H.xml", True),
        "2MASS_2MASS_J": ("SVO", "2MASS.2MASS.J.xml", True),
        "2MASS_2MASS_KS": ("SVO", "2MASS.2MASS.Ks.xml", True),
        "GALEX_GALEX_FUV": ("SVO", "GALEX.GALEX.FUV.xml", True),
        "GALEX_GALEX_NUV": ("SVO", "GALEX.GALEX.NUV.xml", True),
        "GENERIC_JOHNSON_B": ("SVO", "Generic.Johnson.B.xml", True),
        "GENERIC_JOHNSON_I": ("SVO", "Generic.Johnson.I.xml", True),
        "GENERIC_JOHNSON_J": ("SVO", "Generic.Johnson.J.xml", True),
        "GENERIC_JOHNSON_M": ("SVO", "Generic.Johnson.M.xml", True),
        "GENERIC_JOHNSON_R": ("SVO", "Generic.Johnson.R.xml", True),
        "GENERIC_JOHNSON_U": ("SVO", "Generic.Johnson.U.xml", True),
        "GENERIC_JOHNSON_V": ("SVO", "Generic.Johnson.V.xml", True),
        "HERSCHEL_PACS_100": ("SVO", "Herschel.Pacs.green.xml", False),
        "HERSCHEL_PACS_160": ("SVO", "Herschel.Pacs.red.xml", False),
        "HERSCHEL_PACS_70": ("SVO", "Herschel.Pacs.blue.xml", False),
        "HERSCHEL_SPIRE_250": ("SVO", "Herschel.SPIRE.PSW_ext.xml", False),
        "HERSCHEL_SPIRE_350": ("SVO", "Herschel.SPIRE.PMW_ext.xml", False),
        "HERSCHEL_SPIRE_500": ("SVO", "Herschel.SPIRE.PLW_ext.xml", False),
        "IRAS_IRAS_100": ("SVO", "IRAS.IRAS.100mu.xml", True),
        "IRAS_IRAS_12": ("SVO", "IRAS.IRAS.12mu.xml", True),
        "IRAS_IRAS_25": ("SVO", "IRAS.IRAS.25mu.xml", True),
        "IRAS_IRAS_60": ("SVO", "IRAS.IRAS.60mu.xml", True),
        "SLOAN_SDSS_G": ("SVO", "SLOAN.SDSS.g.xml", True),
        "SLOAN_SDSS_I": ("SVO", "SLOAN.SDSS.i.xml", True),
        "SLOAN_SDSS_R": ("SVO", "SLOAN.SDSS.r.xml", True),
        "SLOAN_SDSS_U": ("SVO", "SLOAN.SDSS.u.xml", True),
        "SLOAN_SDSS_Z": ("SVO", "SLOAN.SDSS.z.xml", True),
        "SPITZER_IRAC_I1": ("SVO", "Spitzer.IRAC.I1.xml", True),
        "SPITZER_IRAC_I2": ("SVO", "Spitzer.IRAC.I2.xml", True),
        "SPITZER_IRAC_I3": ("SVO", "Spitzer.IRAC.I3.xml", True),
        "SPITZER_IRAC_I4": ("SVO", "Spitzer.IRAC.I4.xml", True),
        "SPITZER_MIPS_160": ("SVO", "Spitzer.MIPS.160mu.xml", True),
        "SPITZER_MIPS_24": ("SVO", "Spitzer.MIPS.24mu.xml", True),
        "SPITZER_MIPS_70": ("SVO", "Spitzer.MIPS.70mu.xml", True),
        "SWIFT_UVOT_B": ("SVO", "Swift.UVOT.Bband.xml", True),
        "SWIFT_UVOT_UVM2": ("SVO", "Swift.UVOT.UVM2.xml", True),
        "SWIFT_UVOT_UVW1": ("SVO", "Swift.UVOT.UVW1.xml", True),
        "SWIFT_UVOT_UVW2": ("SVO", "Swift.UVOT.UVW2.xml", True),
        "SWIFT_UVOT_U": ("SVO", "Swift.UVOT.Uband.xml", True),
        "SWIFT_UVOT_V": ("SVO", "Swift.UVOT.Vband.xml", True),
        "TNG_NICS_H": ("SVO", "TNG.NICS.H.xml", True),
        "TNG_NICS_J": ("SVO", "TNG.NICS.J.xml", True),
        "TNG_NICS_K": ("SVO", "TNG.NICS.K.xml", True),
        "TNG_OIG_B": ("SVO", "TNG.OIG.B.xml", True),
        "TNG_OIG_R": ("SVO", "TNG.OIG.R.xml", True),
        "TNG_OIG_U": ("SVO", "TNG.OIG.U.xml", True),
        "TNG_OIG_V": ("SVO", "TNG.OIG.V.xml", True),
        "UKIRT_UKIDSS_H": ("SVO", "UKIRT.UKIDSS.H.xml", True),
        "UKIRT_UKIDSS_J": ("SVO", "UKIRT.UKIDSS.J.xml", True),
        "UKIRT_UKIDSS_K": ("SVO", "UKIRT.UKIDSS.K.xml", True),
        "UKIRT_UKIDSS_Y": ("SVO", "UKIRT.UKIDSS.Y.xml", True),
        "UKIRT_UKIDSS_Z": ("SVO", "UKIRT.UKIDSS.Z.xml", True),
        "WISE_WISE_W1": ("SVO", "WISE.WISE.W1.xml", True),
        "WISE_WISE_W2": ("SVO", "WISE.WISE.W2.xml", True),
        "WISE_WISE_W3": ("SVO", "WISE.WISE.W3.xml", True),
        "WISE_WISE_W4": ("SVO", "WISE.WISE.W4.xml", True),

        # JCMT
        "JCMT_SCUBA2_450": ("JCMT", "scuba2_450_transmission.dat", None),
        "JCMT_SCUBA2_850": ("JCMT", "scuba2_850_transmission.dat", None),

        # PLANCK
        "PLANCK_LFI_30": ("PLANCK", "LFI_BANDPASS_F030.txt", None),
        "PLANCK_LFI_44": ("PLANCK", "LFI_BANDPASS_F044.txt", None),
        "PLANCK_LFI_70": ("PLANCK", "LFI_BANDPASS_F070.txt", None),
        "PLANCK_HFI_100": ("PLANCK", "HFI_BANDPASS_F100.txt", None),
        "PLANCK_HFI_143": ("PLANCK", "HFI_BANDPASS_F143.txt", None),
        "PLANCK_HFI_217": ("PLANCK", "HFI_BANDPASS_F217.txt", None),
        "PLANCK_HFI_353": ("PLANCK", "HFI_BANDPASS_F353.txt", None),
        "PLANCK_HFI_545": ("PLANCK", "HFI_BANDPASS_F545.txt", None),
        "PLANCK_HFI_857": ("PLANCK", "HFI_BANDPASS_F857.txt", None),

        # ALMA
        "ALMA_ALMA_3":  ("ALMA", "alma-0-2000-02.dat", (2590, 3570)),
        "ALMA_ALMA_4":  ("ALMA", "alma-0-2000-02.dat", (1840, 2400)),
        "ALMA_ALMA_5":  ("ALMA", "alma-0-2000-02.dat", (1420, 1840)),
        "ALMA_ALMA_6":  ("ALMA", "alma-0-2000-02.dat", (1090, 1420)),
        "ALMA_ALMA_7":  ("ALMA", "alma-0-2000-02.dat", (800, 1090)),
        "ALMA_ALMA_8":  ("ALMA", "alma-0-2000-02.dat", (600, 780)),
        "ALMA_ALMA_9":  ("ALMA", "alma-0-2000-02.dat", (420, 500)),
        "ALMA_ALMA_10": ("ALMA", "alma-0-2000-02.dat", (320, 380)),

        # EUCLID
        "EUCLID_VIS_VIS":  ("EUCLID", "Euclid_VIS.txt", True),
        "EUCLID_NISP_Y":  ("EUCLID", "Euclid_NISP_Y.txt", True),
        "EUCLID_NISP_J":  ("EUCLID", "Euclid_NISP_J.txt", True),
        "EUCLID_NISP_H":  ("EUCLID", "Euclid_NISP_H.txt", True),

        # RUBIN
        "RUBIN_LSST_U":  ("EUCLID", "lsst_u.txt", True),
        "RUBIN_LSST_G":  ("EUCLID", "lsst_g.txt", True),
        "RUBIN_LSST_R":  ("EUCLID", "lsst_r.txt", True),
        "RUBIN_LSST_I":  ("EUCLID", "lsst_i.txt", True),
        "RUBIN_LSST_Z":  ("EUCLID", "lsst_z.txt", True),
        "RUBIN_LSST_Y":  ("EUCLID", "lsst_y.txt", True),
    }

    ## The constructor creates a BroadBand instance in one of the following three ways, depending on the type of
    # the \em bandspec argument:
    #
    #   - if \em bandspec is a string, it looks for a built-in standard broadband filter based on the text segments
    #     in the string, and then loads the corresponding transmission curve from built-in resources
    #     (residing in this package's \c data directory).
    #     To specify a band, it suffices to include just two segments of its name (case insensitive) that
    #     uniquely identify the band, seperated by an underscore or a space. For example, to select the
    #     HERSCHEL_PACS_100 band, one could enter "Herschel 100", "PACS 100", or "HERSCHEL_PACS_100".
    #
    #   - if \em bandspec is a tuple with two numbers, a bolometer-type band is constructed with a uniform
    #     transmission curve in the indicated (min,max) wavelength range (expressed in micron).
    #
    #   - if \em bandspec is a two-dimensional numpy.ndarray, a bolometer-type band is constructed using
    #     \em bandspec[:,0] as wavelength grid and \em bandspec[:,1] as transmission curve. The wavelengths are
    #     assumed to be expressed in micron, while the transmission curves are assumed dimensionless.
    #     The array layout used here is identical to what one would obtain when applying numpy.loadtxt() to the
    #     input file for a SKIRT FileBand.
    #
    # If \em bandspec has a different type, or if its string contents does not unambiguously match a built-in band
    # name, the constructor raises an exception.
    #
    # \note During the first part of construction, _wavelengths are initialized in micron but without explicit units,
    # and _transmissions have arbitrary scale. Thus, these are plain numpy arrays, not astropy quantities. At the
    # end of construction, the _normalize() function attaches the appropriate units.
    #
    def __init__(self, bandspec):
        # construction of built-in band
        if isinstance(bandspec, str):
            # get the built-in band name and format
            self._bandname = self._bandNameFromSpec(bandspec)
            form = self._bandinfo[self._bandname][0]

            # invoke the loader for this format; this sets _wavelengths, _transmissions, and _photoncounter
            if form=="SVO": self._loadSVO()
            elif form=="JCMT": self._loadJCMT()
            elif form=="PLANCK": self._loadPLANCK()
            elif form=="ALMA": self._loadALMA()
            elif form=="EUCLID": self._loadEUCLID()
            else: raise ValueError("Unsuppported built-in band format '{}'".format(form))

        # construction of uniform band
        elif isinstance(bandspec, tuple) and len(bandspec)==2:
            self._bandname = "Uniform"
            self._wavelengths = np.array(list(map(float, bandspec)))
            self._transmissions = np.ones(2)
            self._photoncounter = False

        # construction of a custom band
        elif isinstance(bandspec, np.ndarray) and len(bandspec.shape)==2:
            self._bandname = "Custom"
            self._wavelengths = bandspec[:,0].astype(float)
            self._transmissions = bandspec[:,1].astype(float)
            self._photoncounter = False

        else: raise ValueError("Unsuppported type of band specification '{}'".format(bandspec))

        # normalize the transmission curve and attach the appropriate astropy units
        self._normalize()

    ## This function matches the given band specification with all built-in band names, i.e. the keys of the
    # _bandinfo dictionary. If there is a single match, the corresponding key is returned. If there is no match,
    # or if there are multiple matches, an exception is raised.
    # To match a band name, a band specification must include just two segments of the band name (case insensitive).
    # On both sides of the match, segments are separated by whitespace and/or underscores.
    def _bandNameFromSpec(self, bandspec):
        specsegments = bandspec.upper().replace("_"," ").split()
        result = None
        for bandname in self._bandinfo.keys():
            namesegments = bandname.upper().replace("_"," ").split()
            if all([ (specsegment in namesegments) for specsegment in specsegments ]):
                if result is not None:
                    raise ValueError("Band specification '{}' matches multiple band names, "
                                     "including {} and {}".format(bandspec, result, bandname))
                result = bandname
        if result is None:
            raise ValueError("Band specification '{}' matches no band names".format(bandspec))
        return result

    ## This function loads the transmission curve for the SVO format. It expects _bandname to be set,
    # and it sets _wavelengths, _transmissions, and _photoncounter.
    def _loadSVO(self):
        # get band info
        subdir, filename, self._photoncounter = self._bandinfo[self._bandname]
        filepath = ut.dataPath(BroadBand) / subdir / filename

        # load the XML tree
        with open(filepath, 'r') as bandfile: bandtree = xml.etree.ElementTree.parse(bandfile)

        # get the transmission table (converting wavelengths from Angstrom to micron)
        values = np.array([ float(el.text) for el in bandtree.findall(".//RESOURCE/TABLE/DATA/TABLEDATA[1]/TR/TD") ])
        if len(values) < 4: raise ValueError("Transmission table not found for band '{}'".format(self._bandname))
        self._wavelengths, self._transmissions = np.reshape(values, (-1, 2)).T
        self._wavelengths *= 1e-4

    ## This function loads the transmission curve for the JCMT format. It expects _bandname to be set,
    # and it sets _wavelengths, _transmissions, and _photoncounter.
    def _loadJCMT(self):
        subdir, filename, dummy = self._bandinfo[self._bandname]
        self._photoncounter = False
        filepath = ut.dataPath(BroadBand) / subdir / filename
        self._wavelengths, self._transmissions = np.loadtxt(filepath, unpack=True)

    ## This function loads the transmission curve for the PLANCK format. It expects _bandname to be set,
    # and it sets _wavelengths, _transmissions, and _photoncounter.
    def _loadPLANCK(self):
        # get band info
        subdir, filename, dummy = self._bandinfo[self._bandname]
        self._photoncounter = False
        filepath = ut.dataPath(BroadBand) / subdir / filename

        # load text columns and process depending on instrument type
        if "LFI" in self._bandname:
            frequencies, transmissions = np.loadtxt(filepath, skiprows=1, usecols=(0, 1), unpack=True)
            frequencies *= 1e9  # from GHz to Hz
            wavelengths = const.c.si.value / frequencies
            wavelengths *= 1e6  # from m to micron
        else:
            wavenumbers, transmissions, uncertainties = np.loadtxt(filepath, skiprows=3, usecols=(0, 1, 2), unpack=True)

            # only keep the rows where the transmission value is higher than 10 x the uncertainty
            mask = transmissions > 10. * uncertainties
            wavenumbers = wavenumbers[mask]
            transmissions = transmissions[mask]

            # only keep the rows where the transmission is above 1/5000 of the peak transmission
            peak = np.max(transmissions)
            mask = transmissions > peak / 5000.
            wavenumbers = wavenumbers[mask]
            transmissions = transmissions[mask]

            # convert to wavelengths
            wavenumbers *= 1e2  # from 1/cm to 1/m
            wavelengths = 1. / wavenumbers
            wavelengths *= 1e6  # from m to micron

        # reverse order
        self._wavelengths = np.flipud(wavelengths)
        self._transmissions = np.flipud(transmissions)

    ## This function loads the transmission curve for the ALMA format. It expects _bandname to be set,
    # and it sets _wavelengths, _transmissions, and _photoncounter.
    def _loadALMA(self):
        subdir, filename, (wmin, wmax) = self._bandinfo[self._bandname]
        self._photoncounter = False
        filepath = ut.dataPath(BroadBand) / subdir / filename

        # load text columns
        frequencies, transmissions = np.loadtxt(filepath, usecols=(0, 1), unpack=True)
        frequencies *= 1e9  # from GHz to Hz
        wavelengths = const.c.si.value / frequencies
        wavelengths *= 1e6  # from m to micron

        # only keep the rows inside the wavelength range for this band
        mask = (wmin < wavelengths) * (wavelengths < wmax)
        wavelengths = wavelengths[mask]
        transmissions = transmissions[mask]

        # reverse order
        self._wavelengths = np.flipud(wavelengths)
        self._transmissions = np.flipud(transmissions)

    ## This function loads the transmission curve for the EUCLID format. It expects _bandname to be set,
    # and it sets _wavelengths, _transmissions, and _photoncounter.
    def _loadEUCLID(self):
        # get band info
        subdir, filename, self._photoncounter = self._bandinfo[self._bandname]

        # load text columns and convert from Angstrom to micron
        filepath = ut.dataPath(BroadBand) / subdir / filename
        wavelengths, transmissions = np.loadtxt(filepath, unpack=True)
        wavelengths *= 1e-4

        # remove leading and trailing values that are close to zero
        zeroT = transmissions.max() / 500.
        first = max(0, np.argmax(transmissions > zeroT) - 1)
        last = len(transmissions) - np.argmax(transmissions[::-1] > zeroT) + 1
        self._wavelengths = wavelengths[first:last]
        self._transmissions = transmissions[first:last]

    ## This function normalizes the transmission curve after it has been loaded during construction,
    # and attaches the appropriate astropy units (micron and 1/micron, respectively)
    # It expects _bandname, _wavelengths, _transmissions, and _photoncounter to be set.
    def _normalize(self):
        # for photon counters, convert from response curve to transmission curve (with arbitrary scaling)
        if self._photoncounter:
            self._transmissions *= self._wavelengths

        # normalize the transmission curve to unity
        self._transmissions /= np.trapz(x=self._wavelengths, y=self._transmissions)

        # attach the appropriate astropy units
        self._wavelengths <<= u.micron
        self._transmissions <<= u.micron**(-1)

    # -----------------------------------------------------------------

    ## This function returns an iterable over all built-in band names, in arbitrary order.
    @classmethod
    def builtinBandNames(cls):
        return cls._bandinfo.keys()

    ## This function returns the band's name, i.e. one of the keys of the _bandinfo dictionary, or "Uniform"
    # for uniform bands.
    def name(self):
        return self._bandname

    ## This function returns a tuple of astropy quantities specifying the wavelength range for this band in micron.
    # The transmission may be zero in some (usually small) intervals within the range, but it is guaranteed to be zero
    # outside the range.
    def wavelengthRange(self):
        return self._wavelengths[0], self._wavelengths[-1]

    ## This function returns the minimum wavelength for the band, as an astropy quantity in micron.
    def minWavelength(self):
        return self._wavelengths[0]

    ## This function returns the maximum wavelength for the band, as an astropy quantity in micron.
    def maxWavelength(self):
        return self._wavelengths[-1]

    ## This function returns the pivot wavelength as defined in the class header, as an astropy quantity in micron.
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
    # as an astropy quantity in micron.
    def effectiveWidth(self):
        return 1. / self._transmissions.max()

    ## This function returns a tuple containing a copy of the wavelength and transmission arrays representing
    # the transmission curve for this band. For photon counters the response curve has already been converted to a
    # transmission curve as described in the class header. Both arrays are astropy quantities. The wavelengths
    # are in micron and the transmissions in 1/micron (because of the normalization to unity).
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
        wb = self._wavelengths.value

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
        T = np.exp(np.interp(np.log(w), np.log(wb), _log(self._transmissions.value), left=0., right=0.))
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
