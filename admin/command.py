#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.admin.command Execute one of the PTS command scripts from the command line
#
# This module enables command scripts to be executed from the command line to expose PTS functionality
# to a terminal session. The top-level perform() function is invoked from the \_\_main\_\_ module in the
# \c do package when that package is specified on the python command line. Typical usage is as follows:
#
#     alias pts='export PYTHONPATH=~/SKIRT/PTS9 ; python -m pts.do'
#     ...
#     pts try me
#
# Command scripts are located in the \c do subdirectory of each PTS package (an immediate subdirectory
# of the top-level \c pts directory). To access a script, the first argument on the command line
# (after the -m option) must specify its package and script names separated by a forward slash.
# The remaining command line arguments are passed to the script.
#
# Several shortcuts are allowed for the package and script names, as long as these shortcuts do not
# result in any ambiguities. The package name can be omitted, and if the script name contains underscores,
# each of the segments between the underscores can be used instead of the full script name. Also,
# any name can be shortened to an initial subset of the name.
#
# For example, assume that a command script called \c try_do.py resides in the \c do subdirectory of
# the \c admin package, and  that this script accepts a single string argument. Assuming that none of
# the other command scripts match the same shortcut, this script can be accessed by any of the
# following commands,
#
#     pts admin/try_do me
#     pts ad/tr me
#     pts tr "you and me"
#     pts do "you and me"
#
# A PTS command script must define a top-level function called \c do() with a number of positional
# and/or optional parameters as required by the script. These parameters must be declared in a specific
# way as described below. The perform() function in this module inspects the function signature to
# obtain information about the number of parameters and the semantics of each parameter. This information
# is then used to create a command line argument parser (including help text), and to translate
# command line arguments to function argument values.
#
# For example, the try_do script used in the previous example could have a function declaration as follows:
#
#     def do( aFixedString : (str,"first and only positional argument"),
#             aString : (str,"optional string argument") = "PTS is great",
#             aFloat : (float,"optional float argument") = 3.14,
#             anInteger : (int,"optional integer argument") = 7,
#             ) -> "try the PTS command mechanism":
#
# The following rules apply:
#   - positional arguments are indicated by the lack of a default value; optional arguments have a default value.
#   - optional arguments are translated to command line options with the same name prefixed with "--".
#   - each argument is annotated by a tuple specifying the argument type and a help text.
#   - supported argument types are str, int and float (specified wthout quotes); no other types are allowed.
#   - the function's return value annotation briefly describes the purpose of the script.
#   - the "-h" and "--help" options are automatically added.
#
# With the above definition of the do() function, requesting help from the command line results in:
#
#     $ pts try --help
#     usage: pts [-h] [--aString <str>] [--aFloat <float>] [--anInteger <int>]
#     admin/try_do aFixedString
#
#     admin/try_do: try the PTS command mechanism
#
#     positional arguments:
#     admin/try_do       packagename/scriptname
#     aFixedString       first and only positional argument
#
#     optional arguments:
#     -h, --help         show this help message and exit
#     --aString <str>    optional string argument (default: PTS is great)
#     --aFloat <float>   optional float argument (default: 3.14)
#     --anInteger <int>  optional integer argument (default: 7)
#
# And here's an example of providing optional arguments on the command line:
#
#     $ pts try me --aFloat 8.3 --anInteger 17
#     Starting admin/try_do...
#     Command line arguments are:
#       Fixed string:    me
#       Optional string: PTS is great
#       Float number:    8.3
#       Integer number:  17
#     Finished admin/try_do...
#

# -----------------------------------------------------------------

import argparse
import importlib
import inspect
import pathlib
import sys

# -----------------------------------------------------------------

## This function locates the PTS command script corresponding to the first command line argument,
# and then invokes it with the remaining command line arguments.
def perform():
    # split the command specification into a package name (if any) and a script name
    commandspec = sys.argv[1].split('/')
    if len(commandspec) == 1:
        packagename = ""
        scriptname = commandspec[0]
    elif len(commandspec) == 2:
        packagename = commandspec[0].lower()
        scriptname = commandspec[1].lower()
    else:
        raise ValueError("invalid command spec")

    # get the path to the top-level pts directory
    ptsdir = pathlib.Path(inspect.getfile(inspect.currentframe())).parent.parent

    # make a list of all pts package directories containing a "do" directory
    packagenames = [x.name for x in ptsdir.iterdir() if x.is_dir() and (x/"do").is_dir() ]

    # shorten the list to package names that match the specified package name, if any
    packagenames = [pname for pname in packagenames if pname.lower().startswith(packagename) ]

    # make a list of all script names in the eligible packages
    scriptnames = [ (pname,x.stem) for pname in packagenames for x in (ptsdir/pname/"do").iterdir() if x.suffix==".py" ]

    # shorten the list to script names that match the specified script name
    if '_' in scriptname:
        scriptnames = [ (pname,sname) for (pname,sname) in scriptnames if sname.lower().startswith(scriptname) ]
    else:
        scriptnames = [ (pname,sname) for (pname,sname) in scriptnames \
                                if any([ segm.startswith(scriptname) for segm in sname.lower().split('_') ]) ]

    # verify that the list has exactly one script
    if len(scriptnames) < 1:
        raise ValueError("command not found")
    elif len(scriptnames) > 1:
        raise ValueError("command is ambiguous")

    # construct the command script object and perform its function
    script = CommandScript(*scriptnames[0])
    script.performWithCommandLineArguments()

# -----------------------------------------------------------------

## This public function lists all available PTS commands, per package.
def listCommands():
    # get the path to the top-level pts directory
    ptsdir = pathlib.Path(inspect.getfile(inspect.currentframe())).parent.parent

    # loop over all pts package directories containing a "do" directory
    for pname in sorted([ x.name for x in ptsdir.iterdir() if x.is_dir() and (x/"do").is_dir() ]):
        print("Package {}:".format(pname))

        # loop over all command scripts in this package
        for sname in sorted([ x.stem for x in (ptsdir/pname/"do").iterdir() if x.suffix==".py" ]):
            script = CommandScript(pname, sname)
            print("  {}: {}".format(script.name(), script.description()))

# -----------------------------------------------------------------

## The CommandScript class encapsulates a particular PTS command script. It offers functions to
# obtain information such as the script's name and description, and to execute its do() function
# using the current command line arguments.#
class CommandScript:

    ## The constructor imports the script specified through its package and script names, and obtains
    # the signature of its do() function.
    def __init__(self, packagename, scriptname):
        # remember script identification
        self._name = "{}/{}".format(packagename, scriptname)

        # get and remember the script's do() function and its signature
        importpath = "pts.{}.do.{}".format(packagename, scriptname)
        self._dofunction = vars(importlib.import_module(importpath))["do"]
        self._signature = inspect.signature(self._dofunction)

    ## This function returns the script's command name, i.e. its package and script names joined by a forward slash.
    def name(self):
        return self._name

    ## This function returns the script's description, i.e. the return value annotation of its do() function.
    def description(self):
        return self._signature.return_annotation

    ## This function inspects the arguments offered by the script's do() function, retrieves the
    # appropriate values from the corresponding command line arguments using the standard argument parser,
    # and finally invokes the function with these arguments.
    def performWithCommandLineArguments(self):
        # initialize an argument parser based on the function parameters
        parser = argparse.ArgumentParser(prog="pts",
                                         description="{}: {}".format(self._name, self._signature.return_annotation))
        parser.add_argument(self._name, type=str, help="packagename/scriptname")
        for p in self._signature.parameters.values():
            if p.default == inspect.Parameter.empty:  # no default ==> positional argument
                parser.add_argument(p.name, type=p.annotation[0], help=p.annotation[1])
            else:
                parser.add_argument("--"+p.name, default=p.default, type=p.annotation[0],
                                    help="{} (default: {})".format(p.annotation[1],p.default),
                                    metavar="<{}>".format(p.annotation[0].__name__))

        # make the parser process the command line arguments
        args = vars(parser.parse_args())
        del args[self._name]

        # invoke the do function in the command script with the appropriate arguments
        print("Starting {}..." .format(self._name))
        self._dofunction(**args)
        print("Finished {}." .format(self._name))

# -----------------------------------------------------------------

