#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.do.list_stored_table_or_columns_info List basic metadata for a SKIRT stored table or stored columns file
#
# This script lists relevant metadata about the specified SKIRT stored table or stored columns file.
# For a stored table file, the printed information includes the names, units and ranges for each of the axes,
# and the names and units for each of the quantities. For a stored columns file, the printed information
# includes the names and units and ranges for each of the columns and the number of rows.
#
# Provide the name or (relative or absolute) file path of a SKIRT stored table/columns file as the single argument.
#

# -----------------------------------------------------------------

def do( filepath : (str,"name or path of SKIRT stored table or stored columns file"),
        ) -> "list basic metadata for a SKIRT stored table or stored columns file":

    import pts.storedtable as stab
    from pts.utils.error import UserError

    if filepath.endswith(".stab"):
        stab.listStoredTableInfo(filepath)
    elif filepath.endswith(".scol"):
        stab.listStoredColumnsInfo(filepath)
    else:
        raise UserError("Filename extension should be .stab (for stored table) or .scol (for stored colunns)")

# -----------------------------------------------------------------
