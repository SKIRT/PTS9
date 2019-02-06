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

## This function returns the absolute path to the SKIRT executable built as part of the project structure including
# this PTS source file, or None if no executable is found. The function assumes that a SKIRT release version is built
# under a \c ~/SKIRT9 or \c SKIRT directory in the directory returned by the projectParent() function. The \c SKIRT9
# directory is tried first.
def skirt():
    project = projectParent()
    path = project / "SKIRT9" / "release" / "SKIRT" / "main" / "skirt"
    if path.is_file(): return path
    path = project / "SKIRT" / "release" / "SKIRT" / "main" / "skirt"
    if path.is_file(): return path
    return None

## This function returns the absolute canonical path corresponding to the given path. The input path may
# be specified as a string or as a pathlib.Path object. If the input path is relative, it is interpreted
# relative to the current working directory. If it starts with a tilde, the tilde is expanded to the home
# directory of the current user. In all cases, the path is streamlined by removing and "." and ".." segments
# and following symbolic links.
def absolute(inpath):
    path = pathlib.Path(inpath).expanduser()
    if not path.is_absolute(): path = pathlib.Path.cwd() / path
    return path.resolve()

# -----------------------------------------------------------------
