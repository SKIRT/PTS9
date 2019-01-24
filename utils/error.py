#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.utils.error Custom exceptions for use in PTS
#
# This module currently contains just a single custom exception class for use in PTS
# when a fatal problem is most likely caused by a user error as opposed to a programming error.
#

# -----------------------------------------------------------------

## An instance of this class should be raised by PTS scripts when a fatal problem
# is most likely caused by a user error as opposed to a programming error. This allows top-level code
# (such as the code supporting the execution of PTS scripts from the command line) to catch these
# errors and report them in a user-friendly manner (for example, show usage help instead of a stack trace).
class UserError(Exception):
    ## The constructor accepts the error message (a string) that will be reported to the user
    def __init__(self, message):
        self.message = message

# -----------------------------------------------------------------
