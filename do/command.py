#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.command Execute one of the PTS command scripts from the command line

# -----------------------------------------------------------------

import sys
import inspect
from pathlib import Path

# -----------------------------------------------------------------

## This top-level function locates the PTS command script corresponding to the first command line argument,
# and then invokes it with the remaining command line arguments.
def perform():
    script = findScript()
    exec("from {} import do; do()".format(script))

# -----------------------------------------------------------------

## This helper function locates the PTS command script corresponding to the first command line argument
# and returns its import path as a string.
def findScript():
    # split the command specification into a package name (if any) and a script name
    commandspec = sys.argv[1].replace("_","/").split('/')
    if len(commandspec) == 1:
        packagename = ""
        scriptname = commandspec[0]
    elif len(commandspec) == 2:
        packagename = commandspec[0].lower()
        scriptname = commandspec[1].lower()
    else:
        raise ValueError("invalid command spec")

    # get the path to the top-level pts directory
    ptsdir = Path(inspect.getfile(inspect.currentframe())).parent.parent

    # make a list of all pts package directories containing a "do" directory
    packagenames = [x.name for x in ptsdir.iterdir() if x.is_dir() and (x/"do").is_dir() ]

    # shorten the list to package names that match the specified package name, if any
    packagenames = [pname for pname in packagenames if pname.lower().startswith(packagename) ]

    # make a list of all script names in the eligible packages
    scriptnames = [ (pname,x.stem) for pname in packagenames for x in (ptsdir/pname/"do").iterdir() if x.suffix==".py" ]

    # shorten the list to script names that match the specified script name
    scriptnames = [ (pname,sname) for (pname,sname) in scriptnames if sname.lower().startswith(scriptname) ]

    # verify that the list has exactly one script
    if len(scriptnames) < 1:
        raise ValueError("command not found")
    elif len(scriptnames) > 1:
        raise ValueError("command is ambiguous")

    # format the absolute import name for the script
    return "pts.{}.do.{}".format(*scriptnames[0])

# -----------------------------------------------------------------
