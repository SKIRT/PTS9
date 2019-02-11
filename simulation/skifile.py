#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.simulation.skifile Reading from and adjusting a SKIRT parameter file.
#
# An instance of the SkiFile class in this module allows reading from and adjusting an existing
# SKIRT parameter file (\em ski file).
#

# -----------------------------------------------------------------

import datetime
import lxml.etree as etree
import pts.simulation.units as su
import pts.utils.path as pp

# -----------------------------------------------------------------

## An instance of the SkiFile class represents a particular existing SKIRT parameter file (\em ski file).
# There are functions to read information from the ski file, and to adjust specific items in its contents.
# The class offers two types of functions.
#
#  - "Specific" functions get or set a particular piece of
#    information, such as for example the number of photon packets launched by the simulation. These
#    functions encapsulate explicit knowledge about the ski file structure, and should be updated when that
#    structure changes (preferably in such a way that older ski files continue to be supported).
#  - "Generic" functions allow getting or setting information addressed through a subset of the XPath syntax.
#    These functions require the caller to know (the relevant portion of) the structure of the ski file,
#    and put the maintenance burden in the event of changes on the side of the caller.
#
# Updates made to a SkiFile instance do \em not affect the underlying file; use the saveto() function to save
# the updated contents of a SkiFile instance to another file (or to replace the original file if so desired).
# A SkiFile class instance is always constructed from an existing ski file; creating a new ski file from scratch
# is not supported. To create a new ski file, start SKIRT in interactive mode (without any arguments).
#
class SkiFile:
    # ---------- Constructing and saving -----------------------------

    ## The constructor loads the contents of the specified ski file into a new SkiFile instance.
    # The path may be specified as a string or a pathlib.Path object.
    # It may be absolute, relative to a user's home folder, or relative to the current working directory.
    # The filename of the ski file \em must end with ".ski" or with "_parameters.xml".
    #
    def __init__(self, skiFilePath):
        # get the absolute path and verify the file name
        self._path = pp.absolute(skiFilePath)
        if not self._path.name.lower().endswith((".ski", "_parameters.xml")):
            raise ValueError("Invalid filename extension for ski file")

        # load the XML tree from the ski file (remove blank text to avoid confusing the pretty printer when saving)
        self._tree = etree.parse(str(self._path), parser=etree.XMLParser(remove_blank_text=True))

    ## This function returns the absolute file path of the ski file represented by this SkiFile instance.
    def skiFilePath(self):
        return self._path

    ## This function saves the (possibly updated) contents of the SkiFile instance into the specified file.
    # The file path may be specified as a string or a pathlib.Path object.
    # It may be absolute, relative to a user's home folder, or relative to the current working directory.
    # The name of the file \em must end with ".ski" or with "_parameters.xml".
    # Saving to and thus replacing the ski file from which this
    # SkiFile instance was originally constructed is allowed, but often not the intention.
    def saveTo(self, saveFilePath):
        # get the absolute path and verify the file name
        path = pp.absolute(saveFilePath)
        if not path.name.lower().endswith(".ski"):
            raise ValueError("Invalid filename extension for ski file")

        # update the producer and time attributes on the root element
        root = self._tree.getroot()
        root.set("producer", "Python toolkit for SKIRT (SkiFile class)")
        root.set("time", datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S"))

        # serialize the XML tree
        self._tree.write(str(path), encoding="UTF-8", xml_declaration=True, pretty_print=True)

    # ---------- Generic functions -----------------------------

    ## This function returns an attribute value or element tag as a string. The first argument is an XPath expression
    # relative to the document root that selects exactly one element. The second argument is an attribute name string
    # or None. In the latter case, the element tag (corresponding to the SKIRT property or class name) is returned.
    # If the element is not found, there are multiple elements, or the selected element does not have the specified
    # attribute, an error is raised.
    def getStringAttribute(self, xpath, attribute):
        # get the element(s) and verify there is exactly one
        elems = self._tree.xpath(xpath)
        if len(elems) == 0: raise ValueError("Ski file has no element for xpath '{}'".format(xpath))
        if len(elems) > 1: raise ValueError("Ski file has multiple elements for xpath '{}'".format(xpath))
        elem = elems[0]
        if not etree.iselement(elem): raise ValueError("Ski file xpath '{}' returns a non-element".format(xpath))

        # get the attribute value or the element tag
        value = elem.get(attribute) if attribute is not None else elem.tag
        if not isinstance(value, str): raise ValueError("Ski file element for xpath '{}' has no attribute '{}'" \
                                                        .format(xpath, attribute))
        return value

    ## This function returns a list including the string values of a particular attribute for a number of elements.
    # The first argument is an XPath expression relative to the document root that selects one or more elements.
    # The second argument is the attribute name (the attribute must be present for all selected elements).
    # If the attribute argument is None, the element tags (corresponding to SKIRT class names) are returned instead.
    # If no elements are found, the function returns an empty list. If one or more of the selected element does
    # not have the specified attribute, an error is raised.
    def getStringAttributes(self, xpath, attribute):
        values = []
        # loop over the element(s)
        for elem in self._tree.xpath(xpath):
            if not etree.iselement(elem): raise ValueError("Ski file xpath '{}' returns a non-element".format(xpath))
            # get the attribute value or the element tag
            value = elem.get(attribute) if attribute is not None else elem.tag
            if not isinstance(value, str): raise ValueError("Ski file element for xpath '{}' has no attribute '{}'" \
                                                            .format(xpath, attribute))
            values.append(value)
        return values

    ## This function sets an attribute value from a string. The first argument is an XPath expression relative to
    # the document root that selects exactly one element. The second argument is the attribute name. The third
    # argument is the string value. If the element did not have an attribute with the specified name, it is added.
    # If the element is not found or there are multiple elements, an error is raised.
    def setStringAttribute(self, xpath, attribute, value):
        # get the element(s) and verify there is exactly one
        elems = self._tree.xpath(xpath)
        if len(elems) == 0: raise ValueError("Ski file has no element for xpath '{}'".format(xpath))
        if len(elems) > 1: raise ValueError("Ski file has multiple elements for xpath '{}'".format(xpath))
        elem = elems[0]
        if not etree.iselement(elem): raise ValueError("Ski file xpath '{}' returns a non-element".format(xpath))

        # set the attribute value
        if not isinstance(value, str):
            raise ValueError("Ski file attribute value has type '' instead of string".format(type(value)))
        elem.set(attribute, value)

    ## This function returns a Boolean attribute value as a true Python Boolean value. The string value is
    # considered to represent True if it contains "true", "t", "yes", "y" or "1" (case insensitive),
    # and False otherwise.
    # The function arguments are the same as those for the getStringAttribute() function.
    def getBoolAttribute(self, xpath, attribute):
        value = self.getStringAttribute(xpath, attribute).strip().lower()
        return value == "true" or value == "t" or value == "yes" or value == "y" or value == "1"

    ## This function sets an attribute value from a Boolean value.
    # The function arguments are the same as those for the setStringAttribute() function, except that the
    # value argument is interpreted as a Python Boolean.
    def setBoolAttribute(self, xpath, attribute, value):
        self.setStringAttribute(xpath, attribute, 'true' if value else 'false')

    ## This function returns a numeric attribute value as an integer number.
    # The function arguments are the same as those for the getStringAttribute() function.
    def getIntAttribute(self, xpath, attribute):
        return int(self.getStringAttribute(xpath, attribute))

    ## This function sets an attribute value from an integer number.
    # The function arguments are the same as those for the setStringAttribute() function, except that the
    # value argument must be convertible to int.
    def setIntAttribute(self, xpath, attribute, value):
        self.setStringAttribute(xpath, attribute, str(int(value)))

    ## This function returns a numeric attribute value (representing a dimensionless quantity)
    # as a floating point number.
    # The function arguments are the same as those for the getStringAttribute() function.
    def getFloatAttribute(self, xpath, attribute):
        return float(self.getStringAttribute(xpath, attribute))

    ## This function sets an attribute value representing a dimensionless quantity from a floating point number.
    # The function arguments are the same as those for the setStringAttribute() function, except that the
    # value argument must be convertible to float.
    def setFloatAttribute(self, xpath, attribute, value):
        self.setStringAttribute(xpath, attribute, _prettyStringForFloat(value))

    ## This function returns an attribute value representing a physical quantity with given units as an astropy
    # scalar quantity. If the attribute value does not include a unit string, the function raises an error.
    # The function arguments are the same as those for the getStringAttribute() function.
    def getQuantityAttribute(self, xpath, attribute):
        value = self.getStringAttribute(xpath, attribute)
        segments = value.split()
        if len(segments) != 2:
            raise ValueError("Ski file quantity attribute has no units or invalid format: '{}'".format(value))
        return float(segments[0]) * su.unit(segments[1])

    ## This function sets an attribute value representing a physical quantity with given units.
    # The function arguments are the same as those for the setStringAttribute() function, except that the
    # value argument must be an astropy scalar quantity with appropriate units, and that a skirt unit string
    # can be specified (i.e. the units to be used in the ski file). If the \em skirtUnit argument is omitted,
    # the skirt unit string is taken from the current attribute value. Thus, in this case, the attribute must
    # already be present and include a unit string. In any case, the specified quantity is converted to the
    # units used in the ski file before being written to the ski file.
    def setQuantityAttribute(self, xpath, attribute, value, skirtUnit=None):
        if skirtUnit is None:
            oldvalue = self.getStringAttribute(xpath, attribute)
            segments = oldvalue.split()
            if len(segments) != 2:
                raise ValueError("Ski file quantity attribute has no units or invalid format: '{}'".format(oldvalue))
            skirtUnit = segments[1]
        newvalue = _prettyStringForFloat(value.to(su.unit(skirtUnit)).value) + " " + skirtUnit
        self.setStringAttribute(xpath, attribute, newvalue)

     # ---------- Specific functions -----------------------------

    ## This function returns the number of photon packets launched per simulation segment for primary sources
    def numPrimaryPackets(self):
        return self.getFloatAttribute("//MonteCarloSimulation", "numPackets")

    ## This function sets the number of photon packets launched per simulation segment for primary sources
    def setNumPrimaryPackets(self, value):
        self.setFloatAttribute("//MonteCarloSimulation", "numPackets", value)

    ## This function returns the names of the instruments in the instrument system as a list of strings,
    # in their order of appearance in the ski file.
    def instrumentNames(self):
        return self.getStringAttributes("//InstrumentSystem/instruments/*", "instrumentName")

    ## This function returns the names of the probes in the probe system as a list of strings,
    # in their order of appearance in the ski file.
    def probeNames(self):
        return self.getStringAttributes("//ProbeSystem/probes/*", "probeName")

# -----------------------------------------------------------------

## This function returns a "pretty" string representation for a floating point number suitable for inclusion
# in a ski file. The formatting proceeds in the same way as when SKIRT outputs numbers in a ski file.
def _prettyStringForFloat(value):
    # start with a regular semi-smart conversion
    s = "{:1.10g}".format(float(value))

    # remove leading zeroes and the + sign in the exponent
    s = s.replace("e-0","e-").replace("e+0","e").replace("e+","e")

    # replace 4 or more trailing zeroes by exponent
    zeroes = len(s) - len(s.rstrip('0'))
    if zeroes > 3:
        s = s[:-zeroes] + "e" + str(zeroes)

    return s

# -----------------------------------------------------------------
