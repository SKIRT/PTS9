#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.band.do.list_bands List all broadbands built into PTS
#
# This script lists the names and pivot wavelengths of all built-in broadbands that satisfy
# all of the specified selection criteria:
#  - \em names (string with comma-separated segments): if specified, the band name must contain
#    at least one of these segments
#  - \em wmin (float): if specified, the pivot wavelength must exceed this value
#  - \em wmax (float): if specified, the pivot wavelength must be lower than this value
#
# The bands are sorted on pivot wavelength within each band family.
#

# -----------------------------------------------------------------

def do( names : (str,"band name segments for bands to be listed, comma separated") = "",
        wmin : (float,"smallest pivot wavelength to be listed, in micron") = 1e-99,
        wmax : (float,"largest pivot wavelength to be listed, in micron") = 1e99,
        ) -> "list built-in broadbands":

    import logging
    import astropy.units as u
    import pts.band as bnd

    # load selected bands
    bands = bnd.builtinBands(names, wmin << u.micron, wmax << u.micron)
    if len(bands) == 0:
        logging.info("There are no matching built-in bands")
    else:
        # sort the bands on pivot wavelength within each band family
        bands = sorted(bands, key=bnd.BroadBand.pivotWavelength)
        bands = sorted(bands, key=lambda band: "_".join(band.name().split('_')[:2]))

        # get the length of the longest name, with a minimum to accomodate the title
        maxlen = max([ len(band.name()) for band in bands ] + [10])

        # list band info
        logging.info("| {} | Pivot wavelength".format("Band name".ljust(maxlen)))
        logging.info("|-{}-|-----------------".format("-"*maxlen))
        for band in bands:
            logging.info("| {:s} | {:#.5g}".format(band.name().ljust(maxlen), band.pivotWavelength().to(u.micron)))
        logging.info("|-{}-|-----------------".format("-"*maxlen))

# -----------------------------------------------------------------
