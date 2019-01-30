#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
# Main module for the do package
# -----------------------------------------------------------------

## \package pts.do.__main__ Execute one of the PTS command scripts from the command line
#
# Assuming that the PYTHONPATH environment variable has been set to include the ~PTS/pts directory,
# this __main__ script gets executed automatically when this module is specified on the python command line:
#
#     python -m pts.do
#
# This script configures logging for use from the command line and then simply invokes the
# doWithCommandLineArguments() function of the command module situated in the admin package.
# See there for more information.
#

# -----------------------------------------------------------------

# initialize PTS for use from the Terminal command line
import pts.admin.initialize
pts.admin.initialize.initializePTS()

# invoke the do() function
import pts.admin.command
pts.admin.command.doWithCommandLineArguments()

# -----------------------------------------------------------------
