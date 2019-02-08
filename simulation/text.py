#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.simulation.text Handling SKIRT text input/output files
#
# This module offers functions for reading information from SKIRT text output files (including column text files,
# log file and convergence reports) and for writing SKIRT text input files (column text files).
#

# -----------------------------------------------------------------

import pts.simulation.units as su
import pts.utils.path as pp

# -----------------------------------------------------------------

## This function extracts a numeric value from a text file such as a SKIRT log file or a file produced by
# SKIRT's SpatialGridConvergenceProbe class. The value is returned as an astropy quantity with appropriate units,
# assuming that a valid unit string has been found in the text file.
#
# The text line containing the value is located by a trigger (a text string that must occur on a line
# before the one containing the value) and a header (a text string that must occur on the line containing
# the field). In fact, the trigger can consist of multiple subtriggers, separated by forward slashes, that
# must be triggered in sequence before the function starts looking for the header.
# If no conforming line is found, the function raises an error.
#
# If the selected line has a section between parentheses at the end, this section (including parentheses) is removed.
# Then the line is split in segments on white space. If the last segment is a valid representation of a floating point
# number, this number is returned as a dimensionless quantity. Otherwise, the last two segments represent the value and
# the unit string, respectively. If there are problems converting these representations, the function raises an error.
#
# For example, the following call will return the total gridded dust mass from a spatial grid convergence file:
#
#     getQuantityFromFile("xxx_convergence.dat", "Dust/Total mass", "Gridded")
#
# and the following call will return the total metallic mass for the first particle medium from a log file:
#
#     getQuantityFromFile("xxx_log.txt", "ParticleMedium", "Total metallic mass")
#
def getQuantityFromFile(path, trigger, header):
    path = pp.absolute(path)
    triggers = trigger.split("/")
    triggered = 0
    with open(path) as infile:
        for line in infile:
            # select the line
            if triggered < len(triggers) and triggers[triggered] in line: triggered += 1
            elif triggered==len(triggers) and header in line:
                # remove section between parentheses
                if line.strip().endswith(')'):
                    line = line[:line.rfind('(')]
                segments = line.split()
                if len(segments) >= 2:
                    try:
                        return float(segments[-1]) << su.unit("")
                    except ValueError: pass
                    return float(segments[-2]) << su.unit(segments[-1])

    raise ValueError("Quantity '{}' not found in text file".format(header))

# -----------------------------------------------------------------
