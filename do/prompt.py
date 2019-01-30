#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.prompt Access PTS command scripts and functions from the interactive Python prompt
#
# While being imported, this module configures PTS logging for use from the interactive Python prompt.
# It also defines a function to perform command scripts from the interactive Python prompt.
#

# -----------------------------------------------------------------

# initialize PTS for use from the Terminal command line
import pts.admin.initialize
pts.admin.initialize.initializePTS(prompt=True)

# -----------------------------------------------------------------

## This function performs a command script as if it were invoked from the command line. The
# The \em commandline argument should not include the program part (i.e. "pts"); the first token is
# the specification of the command script. The \em commandline argument can be one of the following:
#   - A single string providing the complete command; in this case the tokens are extracted by splitting
#     on whitespace so a token cannot include whitespace.
#   - An iterable of strings specifying the individual tokens of the command; in this case a token can
#     include whitespace.
#
def do(commandline):

    # build the list of arguments
    if isinstance(commandline,str): arguments = commandline.split()
    else: arguments = list(commandline)

    # perform the command
    import pts.admin.command
    pts.admin.command.doWithCommandLineArguments(arguments)

# -----------------------------------------------------------------
