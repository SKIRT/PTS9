#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.do.plot_stored_table Plot a selected curve from data in a SKIRT stored table file
#
# This script plots a particular 1D curve from the data in a specified SKIRT stored table file.
#
# Provide the name or (relative or absolute) file path of a SKIRT stored table file as the single required argument,
# and use optional arguments to further configure the plot:
#  - hor: zero-based index of the table axis on the horizontal plot axis (default = 0)
#  - ver: zero-based index of the table quantity on the vertical plot axis (default = 0)
#  - ax0, ax1, ax2, ax3, ax4: ordinate value for the table axis with the indicated zero-based index
#    (default = arithmetic or geometric mean of the axis range; ignored for axes not in the table)
#
# Thus, by default, the script plots the first table quantity as a function of the first table axis,
# with half-way values for the other axes, if any.
#
# The plot file is named "FigStoredTable.pdf" and is placed in the current working directory.
#

# -----------------------------------------------------------------

def do( filepath : (str,"name or path of SKIRT stored table file"),
        hor : (int,"zero-based index of the table axis on the horizontal plot axis") = 0,
        ver : (int,"zero-based index of the table quantity on the vertical plot axis") = 0,
        ax0 : (float,"ordinate value for the table axis with index 0") = None,
        ax1 : (float,"ordinate value for the table axis with index 1") = None,
        ax2 : (float,"ordinate value for the table axis with index 2") = None,
        ax3 : (float,"ordinate value for the table axis with index 3") = None,
        ax4 : (float,"ordinate value for the table axis with index 4") = None,
        ) -> "plot a curve from data in a SKIRT stored table file":

    import pts.visual as vis

    vis.plotStoredTableCurve(tableFilePath=filepath, horAxis=hor, verAxis=ver,
                             axis0=ax0, axis1=ax1, axis2=ax2, axis3=ax3, axis4=ax4)

# -----------------------------------------------------------------
