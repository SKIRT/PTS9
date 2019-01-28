#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.admin.do.try_do Test basic PTS command line functionality
#
# You can use this script to test basic PTS command line functionality. It simply outputs the values
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

    import logging
    logging.info("Command line arguments are:")
    logging.info("  Fixed string:    {}".format(aFixedString))
    logging.info("  Optional string: {}".format(aString))
    logging.info("  Float number:    {}".format(aFloat))
    logging.info("  Integer number:  {}".format(anInteger))

# -----------------------------------------------------------------
