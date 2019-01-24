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

import logging
from pts.admin import command

# -----------------------------------------------------------------

## This class overrides the logging formatter to adjust the format for the message level.
class CommandLineLoggingFormatter(logging.Formatter):
    def format(self, record):
        line = super().format(record)
        line = line.replace("<DEBUG>",    ".", 1)
        line = line.replace("<INFO>",     " ", 1)
        line = line.replace("<WARNING>",  "!", 1)
        line = line.replace("<ERROR>",    "* ERROR:", 1)
        line = line.replace("<CRITICAL>", "* CRITICAL ERROR:", 1)
        return line

# -----------------------------------------------------------------

# configure logging facilities for use from the command line
formatter = CommandLineLoggingFormatter(fmt='%(asctime)s.%(msecs)03d <%(levelname)s> %(message)s',
                                        datefmt='%m/%d/%Y %H:%M:%S')
handler = logging.StreamHandler()
handler.setFormatter(formatter)
logging.root.addHandler(handler)
logging.root.setLevel(logging.DEBUG)

# invoke the do() function
command.doWithCommandLineArguments()

# -----------------------------------------------------------------
