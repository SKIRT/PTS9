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
def dataPath(object):
    packagedir = pathlib.Path(inspect.getfile(object)).parent
    if packagedir.name == "do": packagedir = packagedir.parent
    return packagedir / "data"

## This function returns the absolute path to the top-level directory of the pts repository
# as a pathlib.Path object.
def ptsPath():
    return pathlib.Path(inspect.getfile(ptsPath)).parent.parent

## This function returns the absolute path to the PTS resources directory as part of the project structure including
# this PTS source file, or None if no such directory is found.
def ptsResourcesPath():
    path = ptsPath().parent / "resources"
    if path.is_dir(): return path
    return None

## This function returns the absolute path to the directory containing the SKIRT and PTS
# project directories as a pathlib.Path object. In the canonical SKIRT/PTS developer
# directory structure, this is the grandparent of the pts repository.
def projectParentPath():
    return ptsPath().parent.parent

## This function returns the absolute path to the SKIRT executable built as part of the project structure including
# this PTS source file, or None if no executable is found. The function assumes that a SKIRT release version is built
# under a \c ~/SKIRT9 or \c SKIRT directory in the directory returned by the projectParent() function. The \c SKIRT9
# directory is tried first.
def skirtPath():
    project = projectParentPath()
    path = project / "SKIRT9" / "release" / "SKIRT" / "main" / "skirt"
    if path.is_file(): return path
    path = project / "SKIRT" / "release" / "SKIRT" / "main" / "skirt"
    if path.is_file(): return path
    return None

## This function returns the absolute path to the SKIRT resources directory as part of the project structure including
# this PTS source file, or None if no such directory is found. The function assumes that SKIRT is installed
# under a \c ~/SKIRT9 or \c SKIRT directory in the directory returned by the projectParent() function. The \c SKIRT9
# directory is tried first.
def skirtResourcesPath():
    project = projectParentPath()
    path = project / "SKIRT9" / "resources"
    if path.is_dir(): return path
    path = project / "SKIRT" / "resources"
    if path.is_dir(): return path
    return None

## This function returns the absolute canonical path corresponding to the given path. The input path may
# be specified as a string or as a pathlib.Path object. If the input path is relative, it is interpreted
# relative to the current working directory. If it starts with a tilde, the tilde is expanded to the home
# directory of the current user. In all cases, the path is streamlined by removing and "." and ".." segments
# and following symbolic links.
def absPath(inpath):
    path = pathlib.Path(inpath).expanduser()
    if not path.is_absolute(): path = pathlib.Path.cwd() / path
    return path.resolve()

## This function returns the absolute canonical path corresponding to the file path that will be used to save
# a result to file, such as a PDF plot or a PNG image. The returned path is derived from the input arguments:
#  - \em defFilePath: required argument; specifies the default save file path and is used in case the other arguments
#    do not fully specify the path. This path is typically constructed by the immediate caller, e.g. a plot function.
#  - \em outDirPath: optional argument; overrides the directory path in \em defFilePath.
#  - \em outFileName: optional argument; overrides the filename portion of the path in \em defFilePath. It is valid
#    for both \em outDirPath and \em outFileName to be specified, in which case \em defFilePath is ignored. However,
#    in this case it may be more natural to use the \em outFilePath argument instead.
#  - \em outFilePath: optional argument; overrides the complete path, which means that \em defFilePath, \em outDirPath
#    and \em outFileName are ignored. The three optional arguments are typically (but not necessarily) supplied by
#    the caller of the immediate caller.
#  - \em suffix: required argument; string or sequence of strings listing the allowed suffixes for the returned path
#    (including the leading "." in each suffix).
#    If the path resulting from the above rules does not already has one of the specified suffixes, the suffix
#    of the path will be replaced by the first (or only) specified suffix.
#
# The three input paths may be specified as a string or as a pathlib.Path object. If the input path is relative,
# it is interpreted relative to the current working directory. If it starts with a tilde, the tilde is expanded
# to the home directory of the current user.
def savePath(defFilePath, suffix, *, outDirPath=None, outFileName=None, outFilePath=None):

    # construct the path based on defFilePath, outDirPath, outFileName, and outFilePath
    result = None
    if outFilePath is not None: result = absPath(outFilePath)
    elif outDirPath is not None and outFileName is not None: result = absPath(outDirPath) / outFileName
    elif outFileName is not None: result = absPath(defFilePath).with_name(outFileName)
    elif outDirPath is not None: result = absPath(outDirPath) / absPath(defFilePath).name
    else: result = absPath(defFilePath)

    # adjust the suffix as requested
    if isinstance(suffix, str): suffix = [ suffix ]
    if not (result.suffix in suffix):
        result = result.with_suffix(suffix[0])

    return result

# -----------------------------------------------------------------
