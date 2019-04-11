#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.admin.do.create_resource_archives Create publishable archives for SKIRT 9 resources
#
# This script creates archives (ZIP files) containing SKIRT 9 resources that can be copied to the public server.
#
# SKIRT 9 resources are bundled into one or more archives. The core archive is required, the others are optional.
# In other words, SKIRT 9 users must download at least the core resource archive for the code to be functional.
# Each archive is versioned, and the download procedure at the user's side attempts to synchronize the downloaded
# version with the version required by the currently checked out version of the SKIRT 9 C++ code. It is important
# to maintain older versions of the archives on the public server (at least for a while) so that users can continue
# to employ them (for example, because a newer version of SKIRT has a compatibility issue or a show-stopping bug).
#
# This script expects the following directory structure:
#
#     <project-parent-directory>
#         PTS9
#         SKIRT9
#         Resources9
#             StoredTables
#                 TableDirA
#                 TableDirB
#                 TableDirC
#                 ...
#             Publish
#                 Archives
#                 Definitions
#                     Core
#                         include.txt
#                         version.txt
#                         history.txt
#                     ExtraA
#                         include.txt
#                         version.txt
#                         history.txt
#                     ExtraB
#                     ...
#
# An archive can be created for each subdirectory of the Definitions directory. Each of these subdirectories should
# contain three basic text files:
#
#   - include.txt: a list of StoredTables subdirectories to be included in the archive (one directory name per line)
#   - version.txt: the current version of this archive (a single integer number on the first line)
#   - history.txt: a human readable summary of the version history for this archive
#
# The archive will contain the (recursive) contents of the StoredTables subdirectories listed in the include.txt file,
# in addition to the version.txt and history.txt files, placed at the top level. The archive will be named
# "SKIRT9_Resources_name_vxx.zip", where name is replaced by the name of the corresponding Definitions subdirectory,
# and xx is replaced by the version number. The archive will be placed in the Resources9/Publish/Archives directory,
# replacing any existing file with the same name (but leaving files with earlier version numbers untouched).
#
# By default, the script creates an archive for each Definitions subdirectory. Specifying the \c --name=dirname option
# restricts creation to the specified subdirectory.
#

# -----------------------------------------------------------------

def do( name : (str,"create only the resource archive with this name") = "",
        ) -> "create publishable archives for SKIRT 9 resources":

    import logging
    import pathlib
    import zipfile
    import pts.utils as ut

    # get the paths to the key involved directories
    publishdir = ut.projectParentPath() / "Resources9" / "Publish"
    archivedir = publishdir / "Archives"
    definisdir = publishdir / "Definitions"
    stablesdir = publishdir.parent / "StoredTables"

    # create a list of archives to be created
    if len(name)>0:
        archivenames = [ name ]
    else:
        archivenames = [ path.name for path in definisdir.glob("*") if path.is_dir() ]

    # create each archive
    for name in archivenames:
        # get the archive version
        with open(definisdir/name/"version.txt") as versionfile:
            version = int(versionfile.readline().strip())

        # get list of stored table subdirectories to include
        with open(definisdir/name/"include.txt") as includefile:
            subdirnames = [ line.strip() for line in includefile if len(line.strip())>0 ]

        # create the archive
        rootname = "SKIRT9_Resources_{}".format(name)
        archivepath = archivedir / (rootname + "_v{}.zip".format(version))
        logging.info("Creating archive {}...".format(archivepath))
        with zipfile.ZipFile(archivepath, mode='w', compression=zipfile.ZIP_DEFLATED,
                             compresslevel=6) as ziparchive:

            # include the top-level files
            for sourcename in ("version.txt", "history.txt"):
                logging.info("  Including {}".format(sourcename))
                ziparchive.write(definisdir/name/sourcename, arcname=pathlib.Path(rootname)/sourcename)

            # include all contents of the required subdirectories
            for subdirname in subdirnames:
                for source in (stablesdir/subdirname).rglob("*"):
                    if not source.name.startswith("."):
                        sourcename = source.relative_to(stablesdir)
                        logging.info("  Including {}".format(sourcename))
                        ziparchive.write(source, arcname=rootname/sourcename)

# -----------------------------------------------------------------
