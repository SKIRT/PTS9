#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.test.functional Contains the SkirtTestSuite class for performing a suite of SKIRT test cases
#
# An instance of the SkirtTestSuite class in this module represents a suite of functional SKIRT test cases, stored as
# a nested structure of files and directories according to a specific layout, and provides facilities to
# perform the tests, verify the results, and prepare a summary test report.

# -----------------------------------------------------------------

import filecmp
import logging
import multiprocessing
import time
import pts.do
import pts.simulation as sm
import pts.utils as ut

# -----------------------------------------------------------------

## An instance of the SkirtTestSuite class represents a suite of SKIRT test cases, stored as
# a nested structure of files and directories according to a specific layout, and provides facilities to
# perform the tests, verify the results, and prepare a summary test report.
#
# A test suite consists of a set of independent test cases (i.e. test cases can be executed in arbitrary order)
#
# Each test case in a test suite is defined by a collection of files and directories as follows:
#  - a directory with arbitrary name containing all test case files and directories, called the "case directory";
#    a test suite is named after this directory
#  - immediately inside the case directory there is:
#    - exactly one \em ski file with an arbitrary name (with the \c .ski filename extension) specifying the simulation
#      to be performed for the test case
#    - a directory named \c in containing the input files for the simulation, if any
#    - a directory named \c ref containing the reference files for the test, i.e. a copy of the output files
#      generated by a correct simulation run
#    - a directory named \c out to receive the actual output files when the test is performed; the contents
#      of this directory are automatically removed and recreated when running the test case
#    - everything else is ignored, as long as there are no additional files with a \c .ski filename extension
#
# A test suite is defined by a collection of files and directories as follows:
#  - a directory directly or indirectly containing all test cases, called the "suite directory";
#    a test suite is named after this directory
#  - each ski file directly or indirectly contained in the suite directory defines a test case that
#    must adhere to the description above (no other ski files in the same directory, special directories
#    next to the \em ski file, etc.)
#
# For example, a test suite may be structured with nested sub-suites as follows (where each \c CaseN directory
# contains a ski file plus \c ref, \c in, and \c out directories):
#
# \verbatim
# SKIRT Tests
#   SPH simulations
#       Case1
#       Case2
#   Geometries
#     Radial
#         Case1
#         Case2
#     Cylindrical
#         Case1
#         Case2
#         Case3
#     Full 3D
#         Case1
#         Case2
#   Instruments
# \endverbatim
#
# It is also allowed to nest test cases inside another test case, but this is not recommended.
#
class SkirtTestSuite:

    ## The first argument accepted by the constructor indicates the selection of functional tests to be performed:
    #  - "." (a single period): perform all test cases in the standard suite.
    #  - "testcase" (the name of a test case directory): perform all test cases with that name.
    #  - "subsuite" (the name of an intermediate directory in the hierarchy): perform the test cases in all sub-suites
    #    with that name.
    #  - "parentsubsuite/testcase" or "parentsubsuite/subsuite": perform the indicated test case(s) or sub-suite(s)
    #    that reside immediately inside the indicated parent sub-suite; this can disambiguate items with the same name.
    #
    # In addition, the constructor accepts an optional argument specifying the path of the directory containing the
    # complete functional test suite. If specified, the path is interpreted as described for the pts.utils.absPath()
    # function. If omitted, the default path is pts.utils.projectParentPath()/"Functional9".
    def __init__(self, subSuite=".", *, suitePath=None):

        # set the top-level suite path and remember the skirt path
        self._suitePath = ut.absPath(suitePath) if suitePath is not None else ut.projectParentPath()/"Functional9"

        # find all matching sub-suite paths
        if subSuite is None or subSuite == "" or "." in subSuite:
            subSuitePaths = [ self._suitePath ]
        else:
            subSuitePaths = self._suitePath.rglob(subSuite)

        # find all valid test cases in any of the sub-suite paths and remember the paths to the corresponding ski files
        skiPathSet = set()
        for subSuitePath in subSuitePaths:
            for skiPath in subSuitePath.rglob("[!.]*.ski"):
                if len(list(skiPath.parent.glob("[!.]*.ski"))) == 1:
                    skiPathSet.add(skiPath)
        self._skiPaths = sorted(skiPathSet)

        # abort if there are no valid test cases
        if len(self._skiPaths) == 0:
            raise ut.UserError("No valid test cases found for sub-suite specification: '{}'".format(subSuite))

    ## This function returns the number of test cases in this test suite
    def size(self):
        return len(self._skiPaths)

    ## This function prepares the contents of all test case directories in the sub-suite for performing the tests.
    # Specifically, it creates \c in, \c out and \c ref directories next to the ski file, if they don't exist,
    # and it removes all files from the \c out directory (without touching any of its subdirectories, which should
    # not be present anyway).
    def clean(self):
        # loop over all ski files in the sub-suite
        for skipath in self._skiPaths:

            # create test case subdirectories as needed
            casedir = skipath.parent
            (casedir/"in").mkdir(exist_ok=True)
            (casedir/"out").mkdir(exist_ok=True)
            (casedir/"ref").mkdir(exist_ok=True)

            # remove files from the "out" subdirectory
            for path in (casedir/"out").glob("*"):
                if path.is_file(): path.unlink()

    ## This function "endorses" the current output of all test case directories in the sub-suite by
    # replacing the contents of the "ref" directory by the contents of the "out" directory.
    # \note This function destroys the current reference output for test cases; USE WITH CARE!
    def endorse(self):
        # loop over all ski files in the sub-suite
        for skipath in self._skiPaths:
            casedir = skipath.parent

            # remove the contents of the ref directory
            for path in (casedir/"ref").glob("*"):
                if path.is_file(): path.unlink()

            # move the contents of the out directory to the ref directory
            for path in (casedir/"out").glob("*"):
                if path.is_file():
                    target = casedir/"ref"/path.name
                    path.replace(target)

    ## This function performs all tests in the test suite, verifies the results, and prepares a summary test report.
    # The function accepts two optional arguments:
    #  - \em skirtPath: the path to the skirt executable. If specified, the path is interpreted as described for
    #    the pts.utils.absPath() function. If omitted, the default path is used as described for the constructor
    #    of the pts.simulation.skirt.Skirt class.
    #  - \em visual: a string describing the visualizations to be created for the output of each test case, if any.
    #    If specified, the string is a semicolon-separated list of PTS visualization commands as they would be entered
    #    at the command line, omitting the simulation output path argument (which will be inserted automatically).
    def perform(self, *, skirtPath=None, visual=None):
        # create a SKIRT execution context for each core on the host computer
        numCores = min(24, len(self._skiPaths), max(multiprocessing.cpu_count(), 1))
        skirts = [ sm.Skirt(skirtPath) for core in range(numCores) ]

        # inform the user of the fact that the tests are being initiated
        skirtversion = skirts[0].version()
        skirtpath = skirts[0].path()
        logging.info("Using {}".format(skirtversion))
        logging.info("With path {}".format(skirtpath))
        logging.info("Performing {} functional test case(s) in {} parallel processes..." \
                     .format(len(self._skiPaths), numCores))

        # prepare the test case directories
        self.clean()

        # initialize statistics
        statistics = { }

        # initialize the list of remaining test cases and the list of simulations being processed
        cases = self._skiPaths.copy()
        sims = [ ]
        numTotal = len(cases)
        numDone = 0

        # main "event loop" runs as long as there are test cases to be processed
        while True:
            # start new skirt executions where possible, and remove the started test cases from the list
            for skirt in skirts:
                if not skirt.isRunning() and len(cases)>0:
                    skipath = cases.pop(0)
                    sims.append(skirt.execute(skiFilePath=skipath, skiRelative=True, inDirPath="in", outDirPath="out",
                                              numProcesses=1, numThreadsPerProcess=1, console='silent', wait=False))

            # if a simulation is ready, report on it and remove it from the list
            # we only report on one simulation at a time to ensure that new skirt executions are launched frequently
            oldNumDone = numDone
            for simindex, sim in enumerate(sims):
                if not sim.isRunning():
                    numDone += 1
                    status = reportTestCase(sim)
                    if not status in statistics: statistics[status] = 0
                    statistics[status] += 1
                    testname = sim.skiFilePath().relative_to(self._suitePath).parent
                    level = logging.INFO if status == "Succeeded" else logging.WARNING
                    logging.log(level, "{:3d} -- {}: {}".format(numDone, testname, status))
                    del sims[simindex]

                    # create visualizations if requested
                    if visual is not None and status != "Crashed":
                        for command in visual.split(";"):
                            arguments = command.split()
                            if len(arguments)>0:
                                arguments.insert(1, str(sim.outDirPath()))
                                pts.do.doWithCommandLineArguments(arguments)
                    break

            if numDone == numTotal:
                break

            # if no simulation was ready, sleep for a while to avoid consuming lots of CPU
            if numDone == oldNumDone:
                time.sleep(1)

        # log test summary
        logging.info("Summary for {} test case(s):".format(numTotal))
        for status in sorted(statistics.keys()):
            level = logging.INFO if status == "Succeeded" else logging.WARNING
            logging.log(level, "  {}: {}".format(status, statistics[status]))

        # create a consolidated test report at the top suite level
        reportdir = self._suitePath / "_testreports"
        reportname = "{}_testreport_{}.txt".format(ut.timestamp(), numTotal)
        reportdir.mkdir(exist_ok=True)
        with open(reportdir / reportname, "wt") as confile:
            # write header and summary
            confile.write("Using {}\n".format(skirtversion))
            confile.write("With path {}\n".format(skirtpath))
            confile.write("Summary for {} test case(s):\n".format(numTotal))
            for status in sorted(statistics.keys()):
                confile.write("  {}: {}\n".format(status, statistics[status]))
            confile.write("---------------\n")

            # copy reports for all performed test cases
            for skipath in self._skiPaths:
                testname = skipath.relative_to(self._suitePath).parent
                confile.write("{}: ".format(testname))
                with open(skipath.parent/"out"/"_testreport.txt", "rt") as reportfile:
                    for line in reportfile:
                        confile.write(line)

            # write footer
            confile.write("---------------\n")

# -----------------------------------------------------------------

## This function verifies the results for the test case represented by the specified simulation, creates
# a test report in the output directory, and returns an overall status string: "Crashed", "Failed", or "Succeeded".
def reportTestCase(sim):
    # if the simulation has not properly finished by now, it must have crashed
    if sim.status() != "Finished":
        status = "Crashed"
    else:
        # verify the results
        extra, missing, differing = findDifferences(sim)
        if len(extra) + len(missing) + len(differing) > 0:
            status = "Failed"
        else:
            status = "Succeeded"

    # create the report file
    with open(sim.outDirPath() / "_testreport.txt", "wt") as reportfile:
        reportfile.write(status + "\n")
        if status == "Failed":
            if len(extra)>0:
                reportfile.write("  Extra files:\n")
                for name in extra:
                    reportfile.write("    {}\n".format(name))
            if len(missing)>0:
                reportfile.write("  Missing files:\n")
                for name in missing:
                    reportfile.write("    {}\n".format(name))
            if len(differing)>0:
                reportfile.write("  Differing files:\n")
                for name in differing:
                    reportfile.write("    {}\n".format(name))
    return status

## This function looks for relevant differences between the contents of the output and reference directories of the
# test case represented by the specified simulation. The function returns three lists of file names, which should all
# be empty for the test case to be successful:
#   - a list of extra files: files present in the output directory but not in the reference directory.
#   - a list of missing files: files present in the reference directory but not in the output directory.
#   - a list of differing files: files present in both directories but which have a different contents.
#
# The function only handles files of the following types:
#   - "*.dat": the contents of these files must be identical.
#   - "*.fits": the contents of these files must be identical except for the information in one header line
#               (allowing a different time stamp)
#
# Any other files are ignored. Specifically, this means that the following file types are \em not compared:
#   - the SKIRT log file ("*.log")
#   - the parameter file generated by SKIRT ("*_parameters.xml")
#   - visualization files generated by PTS ("*.pdf", "*.png")
def findDifferences(sim):

    # get the output and reference paths
    outpath = sim.outDirPath()
    refpath = outpath.parent / "ref"

    # ensure that neither has any subdirectories
    if len([ path for path in outpath.glob("*") if path.is_dir() ]) > 0:
            raise ValueError("Test case 'out' directory has subdirectories")
    if len([ path for path in refpath.glob("*") if path.is_dir() ]) > 0:
            raise ValueError("Test case 'ref' directory has subdirectories")

    # create list of extra, missing and common files of types ".dat" and ".fits"
    dircomp = filecmp.dircmp(outpath, refpath)
    extra = [ name for name in dircomp.left_only if name.endswith((".dat",".fits")) ]
    missing = [ name for name in dircomp.right_only if name.endswith((".dat",".fits")) ]
    common = [ name for name in dircomp.common if name.endswith((".dat",".fits")) ]

    # compare common files, focusing on those that aren't trivially equal.
    matches, mismatches, errors = filecmp.cmpfiles(outpath, refpath, common, shallow=False)
    differing = [ name for name in mismatches + errors
                    if not name.endswith(".fits") or not _areEqualFitsFiles(outpath/name, refpath/name) ]

    # return the lists of extra, missing, and differing files
    return extra, missing, differing

## This helper function returns True if the specified FITS files are equal except for the information
#  in a single header record, and False otherwise.
def _areEqualFitsFiles(filePath1, filePath2):
    return _areEqualLists(_readBlocks(filePath1), _readBlocks(filePath2), 1)

## This helper function returns True if the specified lists are equal except for the contents of
#  a number of items not larger than \em allowedDiffs, and False otherwise.
def _areEqualLists(list1, list2, allowedDiffs):

    # the lists must have the same length, which must be at least 2 (to avoid everything being read into 1 line)
    length = len(list1)
    if length < 2 or length != len(list2): return False

    # compare the lists item by item
    diffs = 0
    for index in range(length):
        if list1[index] != list2[index]:
            # verify against allowed number of differences
            diffs += 1
            if diffs > allowedDiffs: return False

    # no relevant differences
    return True

## This helper function reads blocks of 80 bytes from the specified binary file into a list of bytes objects,
# and returns the list.
def _readBlocks(filepath):
    result = []
    with open(filepath, 'rb') as f:
        while True:
            block = f.read(80)
            if len(block) == 0: break
            result.append(block)
    return result

# -----------------------------------------------------------------
