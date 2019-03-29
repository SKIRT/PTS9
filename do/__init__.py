#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
#  Package initialization file
# -----------------------------------------------------------------

## \package pts.do Exposing PTS functionality to the command line
#
# The Python scripts residing in the (optional) \c do sub-directory of each PTS package are intended to be invoked
# from the command line, exposing (part of) the package's functionality for interactive use in a terminal session.
# The \c __main__ module in this package enables executing any of these scripts from the command line,
# without having to specify the exact location.
#

from .command import doWithCommandLineArguments, listCommands, CommandScript
from .initialize import initializePTS
