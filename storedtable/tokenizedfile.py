#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.tokenizedfile Read a text file token by token
#
# The class in this module allows reading a text file token by token (where tokens are separated by whitespace).
#

# -----------------------------------------------------------------

## This class allows reading a text file token by token (where tokens are separated by whitespace),
# while still allowing to skip complete lines where needed.
class TokenizedFile:

    ## The constructor accepts an open file object. The file must be opened in text mode.
    def __init__(self, file):
        self.file = file
        self.line = []

    ## This function returns the next token in the file, skipping white space (including line endings) as needed.
    def next(self):
        while len(self.line)==0:
            line = self.file.readline()
            if not line:
                raise StopIteration
            self.line = line.split()
        return self.line.pop(0)

    ## This function skips any remaining tokens on the current line (without advancing to the next line).
    def skipToEndOfLine(self):
            self.line = []

    ## This function skips tokens as follows. If the current line has one or more remaining tokens,
    # these tokens are skipped. Otherwise, all tokens on the next line are skipped.
    def skipLine(self):
        if len(self.line)>0:
            self.line = []
        else:
            self.file.readline()

    ## This function skips any lines that start with a '#' character. If the current line has remaining tokens,
    # it is skipped only if the first remaining token starts with a '#'.
    def skipHeaderLines(self):
        if len(self.line)>0 and not self.line[0].startswith('#'): return
        while True:
            line = self.file.readline()
            if not line.startswith('#'): break
        self.line = line.split()

# -----------------------------------------------------------------
