#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.admin.do.remove_pyc Remove all compiled Python files from the local PTS repository
#
# This script removes all compiled Python files (\c *.pyc and \c __pycache__ directories) from the
# local PTS repository. Although these files are ignored by git, removing them can be useful to
# clean the local repository after one or more Python source files have been relocated or renamed.
#

# -----------------------------------------------------------------

def do() -> "remove all compiled Python files from the local PTS repository":
    import logging
    import pts.utils as ut

    # get the path to the top-level pts directory
    ptsdir = ut.ptsPath()

    # remove all .pyc files
    for pyc in ptsdir.rglob("*.pyc"):
        logging.info("Removing {!s}".format(pyc))
        pyc.unlink()

    # remove all __pycache__ directories
    for pyc in ptsdir.rglob("__pycache__"):
        logging.info("Removing {!s}".format(pyc))
        pyc.rmdir()

# -----------------------------------------------------------------
