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

## This function upgrades the specified SKIRT 9 parameter file (\em ski file) so that it becomes appropriate for the
# latest SKIRT 9 version. The upgrade process supports all ski files created by the SKIRT project version 9 tools
# (including the SKIRT command line Q&A and the graphical MakeUp wizard) since SKIRT 9 was publicly released.
# Ski files created for older SKIRT versions (such as SKIRT 7 and 8) cannot be upgraded to SKIRT 9 automatically.
#
# The function accepts three arguments:
# - inpath: the absolute or relative path to the ski file to be handled; the filename extension should be ".ski"
# - backup: if the input ski file needs upgrading and \em backup is True (the default), a copy of the input ski file
#           is created with a filename including a time stamp and ending with "_backupski.xml"
# - replace: if the input ski file needs upgrading and \em replace is True (the default), the input ski file is
#            overwritten by the upgraded version; otherwise it is saved as a new file using a similar filename ending
#            with "_upgradedski.xml"
#
def upgradeSkiFile(inpath, *, backup=True, replace=True):
    # load the ski file
    inpath = ut.absPath(inpath)
    try:
        ski = sm.SkiFile(inpath)
    except SyntaxError:
        logging.error("File does not contain well-formed XML: {}".format(inpath))
        return

    # verify the ski file format version
    try:
        version = ski.getStringAttribute("/skirt-simulation-hierarchy", "format")
    except ValueError:
        logging.error("XML file does not have ski file format: {}".format(inpath))
        return
    if version!="9":
        logging.error("Ski file is older than version 9: {}".format(inpath))
        return

    # perform the upgrade in the XML tree in memory, keeping track of whether the contents was actually changed
    changed = "up" in inpath.name  ## stub

    # save the upgraded version if needed
    if changed:
        if backup:
            inpath.rename(inpath.with_name(inpath.stem + "_" + ut.timestamp() + "_backupski.xml"))
        if replace:
            ski.saveTo(inpath)
            logging.warning("Ski file UPGRADED:  {}".format(inpath))
        else:
            outpath = inpath.with_name(inpath.stem + "_upgradedski.xml")
            ski.saveTo(outpath)
            logging.warning("Ski file UPGRADED:  {} --> {}".format(inpath, outpath.name))
    else:
        logging.info("Ski file unchanged: {}".format(inpath))

# -----------------------------------------------------------------
