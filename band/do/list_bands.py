#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.band.do.list_bands List all broadbands built into PTS
#
# This script lists the names of all broadbands built into PTS with their corresponding pivot wavelength.
# The bands are sorted on pivot wavelength within each band family.
#

# -----------------------------------------------------------------

def do() -> "list all built-in broadbands":

    import logging
    import pts.band as bnd

    # load all bands
    bands = [ bnd.BroadBand(name) for name in bnd.builtinBandNames() ]
    logging.info("There are {} built-in bands:".format(len(bands)))

    # sort them on pivot wavelength within each band family
    bands = sorted(bands, key=bnd.BroadBand.pivotWavelength)
    bands = sorted(bands, key=lambda b: b.name().split('_')[0])

    # list band info
    logging.info("| Band name          | Pivot wavelength")
    logging.info("|--------------------|-----------------")
    for band in bands:
        logging.info("| {:18s} | {:1.5g}".format(band.name(), band.pivotWavelength()))
    logging.info("|--------------------|-----------------")

# -----------------------------------------------------------------
