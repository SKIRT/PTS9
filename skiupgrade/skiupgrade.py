#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.skiupgrade.skiupgrade Contains the upgradeSkiFile function for upgrading SKIRT parameter files
#
# The upgradeSkiFile function in this module allows upgrading SKIRT parameter files (\em ski files)
# between versions of SKIRT 9.

# -----------------------------------------------------------------

import logging
import pts.simulation as sm
import pts.utils as ut

# -----------------------------------------------------------------

## This function upgrades the specified SKIRT 9 parameter file (\em ski file) so that it becomes appropriate for the
# latest SKIRT 9 version. The upgrade process supports all ski files created by the SKIRT project version 9 tools
# (including the SKIRT command line Q&A and the graphical MakeUp wizard) since SKIRT 9 was publicly released.
# Ski files created for older SKIRT versions (such as SKIRT 7 and 8) cannot be upgraded to SKIRT 9 automatically.
# The function logs an appropriate error message when an unsupported file is specified.
#
# The function accepts three arguments:
# - inpath: the absolute or relative path to the ski file to be handled; the filename extension should be ".ski"
# - backup: if the input ski file needs upgrading and \em backup is True (the default), a copy of the input ski file
#           is created with a filename including a time stamp and ending with "_backupski.xml"
# - replace: if the input ski file needs upgrading and \em replace is True (the default), the input ski file is
#            overwritten by the upgraded version; otherwise it is saved as a new file using a similar filename ending
#            with "_upgradedski.xml"
#
def upgradeSkiFile(inpath, *, backup=True, replace=True):
    # load the ski file
    inpath = ut.absPath(inpath)
    try:
        ski = sm.SkiFile(inpath)
    except SyntaxError:
        logging.error("File does not contain well-formed XML: {}".format(inpath))
        return

    # verify the ski file format version
    try:
        version = ski.getStringAttribute("/skirt-simulation-hierarchy", "format")
    except ValueError:
        logging.error("XML file does not have ski file format: {}".format(inpath))
        return
    if version!="9":
        logging.error("Ski file is older than version 9: {}".format(inpath))
        return

    # perform the upgrade in the XML tree in memory, keeping track of whether the contents has actually changed
    changed = False
    for condition,templates in _getUpgradeDefinitions():
        changed |= ski.transformIf(condition, templates)

    # save the upgraded version if needed
    if changed:
        if backup:
            inpath.rename(inpath.with_name(inpath.stem + "_" + ut.timestamp() + "_backupski.xml"))
        if replace:
            ski.saveTo(inpath)
            logging.warning("Ski file UPGRADED:  {}".format(inpath))
        else:
            outpath = inpath.with_name(inpath.stem + "_upgradedski.xml")
            ski.saveTo(outpath)
            logging.warning("Ski file UPGRADED:  {} --> {}".format(inpath, outpath.name))
    else:
        logging.info("Ski file unchanged: {}".format(inpath))

# -----------------------------------------------------------------

## This private function returns a sequence of 2-tuples, each defining the XPath condition and XSLT template
# for a single modification to the ski file format.
#
# Using XSLT is extremely flexible, but unfortunately the XSLT language is fairly obscure an thus has a steep
# learning curve. To alleviate this problem to some extent, we provide a set of functions that generate
# upgrade definitions for specific types of changes, such as, for example, changing the name of a property.
# That way, as soon as the XSLT sheet for a particular type of change has been developed, it can be more easily
# reused for other, similar changes. As a result, this function consists of a sequence of calls to definition
# generators. The git hash and date listed for each (set of) calls identifies the change in the SKIRT 9 code
# requiring the corresponding ski file modification(s).
#
# New generators will have to be added as the need arises. Examples of specific XSLT templates can be found
# in the previous version of PTS in the source file <tt>~/PTS/pts/core/prep/upgradeskifile.py</tt>.
#
def _getUpgradeDefinitions():
    return [

        # SKIRT 9 master git xxxx (14 jan 2020): replace single bulk velocity with velocity field
        # !!! We just remove the zero-velocity components without upgrading non-zero velocities
        # !!! A true upgrade is hard to implement and would be applicable to very few ski files
        _removeScalarPropertyWithValue("GeometricMedium", "velocityX", "0 "),
        _removeScalarPropertyWithValue("GeometricMedium", "velocityY", "0 "),
        _removeScalarPropertyWithValue("GeometricMedium", "velocityZ", "0 "),
    ]

# -----------------------------------------------------------------

## Generates the definition for unconditionally adding a scalar property to a given type.
def _addScalarProperty(typeName, propName, propStringValue):
    return ('''//{0}[not(@{1})]'''.format(typeName, propName),
            '''
            <xsl:template match="//{0}[not(@{1})]">
                <xsl:element name="{0}">
                    <xsl:apply-templates select="@*"/>
                    <xsl:attribute name="{1}">
                        <xsl:value-of select="{2}"/>
                    </xsl:attribute>
                    <xsl:apply-templates select="node()"/>
                </xsl:element>
            </xsl:template>
            '''.format(typeName, propName, propStringValue))

## Generates the definition for removing a scalar property with a value starting with a given string from a given type.
def _removeScalarPropertyWithValue(typeName, oldPropName, oldValue):
    return ('''//{0}[starts-with(@{1},'{2}')]'''.format(typeName, oldPropName, oldValue),
            '''
            <xsl:template match="//{0}[starts-with(@{1},'{2}')]/@{1}">
            </xsl:template>
            '''.format(typeName, oldPropName, oldValue))

## Generates the definition for changing the name of a scalar property for a given type.
def _changeScalarPropertyName(typeName, oldPropName, newPropName):
    return ('''//{0}/@{1}'''.format(typeName, oldPropName),
            '''
            <xsl:template match="//{0}/@{1}">
                <xsl:attribute name="{2}">
                    <xsl:value-of select="."/>
                </xsl:attribute>
            </xsl:template>
            '''.format(typeName, oldPropName, newPropName))

## Generates the definition for changing the name of a compound property for a given type.
def _changeCompoundPropertyName(typeName, oldPropName, newPropName):
    return ('''//{0}/{1}'''.format(typeName, oldPropName),
            '''
            <xsl:template match="//{0}/{1}">
                <xsl:element name="{2}">
                    <xsl:apply-templates select="@*|node()"/>
                </xsl:element>
            </xsl:template>
            '''.format(typeName, oldPropName, newPropName))

# -----------------------------------------------------------------
