#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.admin.do.list_dependencies List external package dependencies for PTS
#
# This script lists all non-PTS packages on which the PTS code directly depends, including standard
# packages and third-party packages, and indicates whether these packages are currently installed or not.
# To accomplish this feat, the script reads all PTS source files looking for import statements of non-PTS packages.
# It then obtains a list of installed packages through a standard-library function.
#
# Note that the listed external packages in turn may have their additional dependencies, which are not checked.
#

# -----------------------------------------------------------------

def do() -> "list external package dependencies for PTS":
    import logging
    import pkgutil
    import re
    import pts.utils as ut

    # ----- find dependencies -----

    # get the path to the top-level pts directory
    ptsdir = ut.ptsPath()

    # initialize the set of dependencies (i.e. package names)
    packages = set()

    # compile the regular expression matching the module name in an import statement
    regex = re.compile(r"(from\s+(?P<module1>(\w|\.)+)\s+import)|(import\s+(?P<module2>(\w|\.)+)(\s|\Z))")

    # loop over all lines in all Python source files in PTS
    for pypath in ptsdir.rglob("*.py"):
        with open(pypath, 'r') as pyfile:
            for line in pyfile:

                # if the line contains the word "import", we need to look further
                if "import" in line:

                    # remove comments and split in statements
                    for statement in line.split('#',2)[0].split(';'):

                        # parse the module token from the statement
                        match = regex.search(statement)
                        if match:
                            module = match.group('module1')
                            if not module: module = match.group('module2')

                            # use the portion before the first dot as the package name
                            #   - empty package means relative import and thus internal to PTS
                            #   - "pts" package means internal to PTS
                            package = module.split('.',2)[0]
                            if len(package) > 0 and package != "pts":
                                packages.add(package)

    # ----- find availability -----

    # build a set of available package names
    available = set([ module.name for module in pkgutil.iter_modules() ])
    available.update(("sys","time"))   # add these as special packages that are always installed

    # output a sorted list of packages with their availability
    logging.info("PTS depends on {} packages:".format(len(packages)))
    for package in sorted(packages):
        if package in available:
            logging.info("  {} -- installed".format(package))
        else:
            logging.warning("  {} -- NOT INSTALLED".format(package))

# -----------------------------------------------------------------
