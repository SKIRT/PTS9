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
# This script simply invokes the perform() function of the command module situated in the admin package,
# so that the bulk of the code is contained in file with a regular name (i.e. without underscores).
#

# -----------------------------------------------------------------

# import the command module and invoke the perform() function
from pts.admin.command import perform
perform()

# -----------------------------------------------------------------
