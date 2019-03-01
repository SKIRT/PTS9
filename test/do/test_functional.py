#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.test.do.test_functional Perform all or a selection of the standard functional SKIRT tests
#
# This script performs all or a selection of the standard functional SKIRT test cases.
# The test case definitions are expected to reside in the \c SKIRT/Functional9 directory hierarchy.
# See the pts.test.functional module for more information.
#
# The script takes a single positional string argument, which can be one of the following:
#  - "." (a single period): perform all test cases in the standard suite.
#  - "testcase" (the name of a test case directory): perform all test cases with that name.
#  - "subsuite" (the name of an intermediate directory in the hierarchy): perform the test cases in all sub-suites
#    with that name.
#  - "parentsubsuite/testcase" or "parentsubsuite/subsuite": perform the indicated test case(s) or sub-suite(s)
#    that reside immediately inside the indicated parent sub-suite; this can disambiguate items with the same name.
#

# -----------------------------------------------------------------

def do( subSuite : (str,"name of sub-suite or test case to be performed or '.' to perform all"),
        ) -> "perform all or a selection of the standard functional SKIRT tests":

    import pts.test
    pts.test.SkirtTestSuite(subSuite=subSuite).perform()

# -----------------------------------------------------------------
