#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.utils.config Configuration utilities
#
# This module offers utilities related to the current configuration of PTS.
#

# -----------------------------------------------------------------

## This global variable is set by the setInteractive() function in this module, which is typically invoked
# during initialization of PTS. A value of true means that PTS operates in interactive mode, where a user
# expects, for example, visualization to occur in a window or a notebook cell rather than being saved to a file.
_interactive = False

# -----------------------------------------------------------------

## This function sets the PTS interactive mode to the specified \em mode value (i.e.\ True or False).
def setInteractive(mode):
    global _interactive
    _interactive = bool(mode)  # guarantee that the value is either True or False

## If the \em mode argument is omitted or has a value of None, this function returns the current interactive mode
# (as set by the most recent call to the setInteractive() function). If the \em mode argument is specified and has
# value other than None, its value is returned (after conversion to bool).
def interactive(mode=None):
    global _interactive
    if mode is None: return _interactive
    else: return bool(mode)    # guarantee that the value is either True or False

# -----------------------------------------------------------------
