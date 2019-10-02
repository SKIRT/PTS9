#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.skiupgrade.do.upgrade_ski_files Upgrade ski files in a given directory to the latest version of SKIRT
#
# This script upgrades all SKIRT 9 ski files in a given directory to be appropriate for the latest version of SKIRT 9.
# Specifically, the script iterates over all files in the directory with the .ski filename extension and performs
# as follows for each:
# - if the file is not a SKIRT parameter file, or it is intended for a SKIRT version before version 9, a warning
#   is issued and the file is otherwise ignored.
# - if the file is a SKIRT 9 parameter file and it is up to date with the latest version of SKIRT 9, an informational
#   message is issued and the file is otherwise ignored.
# - if the file is a SKIRT 9 parameter file that needs upgrading, it is first copied to a backup version with a
#   filename including a time stamp, and it is then replaced by a version upgraded to the latest version of SKIRT 9.
#
# The script takes a single positional string argument specifying the path to the directory containing the ski files to
# be upgraded, or "." (a single period) for the current directory.
#
# See the pts.skiupgrade.skiupgrade module for more information on the ski file upgrade process.
#

# -----------------------------------------------------------------

def do( skiDirPath : (str,"directory containing the ski files to be upgraded"),
        ) -> "upgrade ski files in a given directory to the latest version of SKIRT 9":

    import pts.skiupgrade
    import pts.utils as ut

    for skipath in ut.absPath(skiDirPath).glob("*.ski"):
        pts.skiupgrade.upgradeSkiFile(skipath)

# -----------------------------------------------------------------
