#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.do.list_stored_table_info List basic metadata for a SKIRT stored table file
#
# This script lists relevant metadata about the specified SKIRT stored table file. The printed information
# includes the names, units and ranges for each of the axes, and the names and units for each of the quantities.
#
# Provide the name or (relative or absolute) file path of a SKIRT stored table file as the single argument.
#

# -----------------------------------------------------------------

def do( filepath : (str,"name or path of SKIRT stored table file"),
        ) -> "list basic metadata for a SKIRT stored table file":

    import pts.storedtable as stab

    stab.listStoredTableInfo(filepath)

# -----------------------------------------------------------------
