#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.admin.do.try_do Test basic PTS command line functionality.
#
# You can use this script to test basic PTS command line functionality. It simply prints the values
# of its argument (a single positional argument and three optional arguments).
#
# The script can also serve as a starting template for developing other scripts.
#

# -----------------------------------------------------------------

def do( aFixedString : (str,"first and only positional argument"),
        aString : (str,"optional string argument") = "PTS is great",
        aFloat : (float,"optional float argument") = 3.14,
        anInteger : (int,"optional integer argument") = 7,
        ) -> "try the PTS command mechanism":

    print("Command line arguments are:")
    print("  Fixed string:    {}".format(aFixedString))
    print("  Optional string: {}".format(aString))
    print("  Float number:    {}".format(aFloat))
    print("  Integer number:  {}".format(anInteger))

# -----------------------------------------------------------------
