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

## This function returns the absolute path to the \c data subdirectory for the PTS package
# in which the given object is defined.
def data(object):
    packagedir = pathlib.Path(inspect.getfile(object)).parent
    if packagedir.name == "do": packagedir = packagedir.parent
    return packagedir / "data"

## This function returns the absolute path to the top-level directory of the pts repository
# as a pathlib.Path object.
def pts():
    return pathlib.Path(inspect.getfile(pts)).parent.parent

## This function returns the absolute path to the directory containing the SKIRT and PTS
# project directories as a pathlib.Path object. In the canonical SKIRT/PTS developer
# directory structure, this is the grandparent of the pts repository.
def projectParent():
    return pts().parent.parent

# -----------------------------------------------------------------
