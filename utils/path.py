#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.utils.path Special paths related to PTS
#
# This module allows retrieving some special paths related to PTS, such as for example
# the path to the pts repository.
#

# -----------------------------------------------------------------

import inspect
import pathlib

# -----------------------------------------------------------------

## This function returns the absolute path to the top-level directory of the pts repository
# as a pathlib.Path object.
def pts():
    return pathlib.Path(inspect.getfile(inspect.currentframe())).parent.parent

## This function returns the absolute path to the directory containing the SKIRT and PTS
# project directories as a pathlib.Path object. In the canonical SKIRT/PTS developer
# directory structure, this is the grandparent of the pts repository.
def projectParent():
    return pts().parent.parent

# -----------------------------------------------------------------
