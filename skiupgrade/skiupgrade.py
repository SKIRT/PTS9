#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.skiupgrade.skiupgrade Contains the upgradeSkiFile function for upgrading SKIRT parameter files
#
# The upgradeSkiFile function in this module allows upgrading SKIRT parameter files (\em ski files)
# between versions of SKIRT 9.

# -----------------------------------------------------------------

import logging
import pts.simulation as sm
import pts.utils as ut

# -----------------------------------------------------------------

## This function
# Specifically, the script iterates over all files with the .ski filename extension and performs as follows for each:
# - if the file is not a SKIRT parameter file, or it is intended for a SKIRT version before version 9, a warning
#   is issued and the file is otherwise ignored.
# - if the file is a SKIRT 9 parameter file and it is up to date with the latest version of SKIRT 9, a informational
#   message is issued and the file is otherwise ignored.
# - if the file is a SKIRT 9 parameter file that needs upgrading, it is first copied to a backup version with a
#   filename including a time stamp, and it is then replaced by a version upgraded to the latest version of SKIRT 9.

def upgradeSkiFile(inpath, *, backup=True, replace=True):
    inpath = ut.absPath(inpath)
    logging.info("Processing ski file {}".format(inpath))

    #outpath = inpath if replace else inpath.with_name(inpath.stem+"_upgraded.ski")

# -----------------------------------------------------------------
