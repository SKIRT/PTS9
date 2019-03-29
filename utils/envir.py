#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.utils.envir System environment utilities
#
# This module offers system environment related utilities specific to PTS.
#

# -----------------------------------------------------------------

import datetime
import getpass
import socket

# -----------------------------------------------------------------

## This function returns a string representing the current time and date in the format "YYYY-MM-DD--hh-mm-ss".
# This format ensures proper collation (strings sort in date/time order), is easy to read for a human user,
# and can be used as part of a filename on any platform (since there are no nasty characters).
def timestamp():
    return datetime.datetime.now().strftime("%Y-%m-%d--%H-%M-%S")

# -----------------------------------------------------------------

## This function returns the name of the currently logged-in user.
def username():
    return getpass.getuser()

# -----------------------------------------------------------------

## This function returns the name of the local host.
def hostname():
    return socket.getfqdn()

# -----------------------------------------------------------------
