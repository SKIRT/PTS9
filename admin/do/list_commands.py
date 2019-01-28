#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.admin.do.list_commands List all PTS command scripts, per package
#
# This script lists all available PTS command scripts, per package. Packages are sorted alphabetically, and
# commands are sorted alphabetically within each package.
#

# -----------------------------------------------------------------

def do() -> "list all PTS command scripts":

    from pts.admin import command
    command.listCommands()

# -----------------------------------------------------------------
