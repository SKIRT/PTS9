#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.moviefile Creating a movie file from a sequence of still images
#
# An instance of the MovieFile class in this module allows creating a movie file from a sequence of still images.

# -----------------------------------------------------------------

import pts.utils as ut
import subprocess

# -----------------------------------------------------------------
#  MovieFile class
# -----------------------------------------------------------------

## An instance of the MovieFile class allows creating a movie file from a sequence of still images. Use the constructor
# to specify the filename and the frame size, then insert the frames one at a time with the addFrame() function,
# and finally flush all movie information to the file with the close() function. The resulting movie file uses MPEG-4
# encoding and must have the .mp4 filename extension.
#
# This class requires \c ffmpeg (https://ffmpeg.org) to be installed and in the PATH. The Anaconda Python distribution
# offers a package called "ffmpeg" that does just this. In other cases, the installation procedure might differ.
#
class MovieFile:
    ## The constructor launches \c ffmpeg to create a movie output file (replacing any existing file with the same
    # name). When the constructor finishes, \c ffmpeg is waiting for movie frames to be sent to it.
    #
    # The constructor accepts the following arguments:
    # - filepath: the absolute or relative path of the movie output file, which \em must end with the .mp4 extension.
    # - shape: the number of movie frame pixels in the x and y directions; the default value is (800,600).
    # - rate: the number of frames per second in the output movie; the default value is 24 fps.
    #
    def __init__(self, outFilePath, shape=(800,600), rate=24):
        outFilePath = ut.absPath(outFilePath)
        assert outFilePath.suffix.lower() == ".mp4"

        # remember the frame shape
        self._shape = shape

        # ensure that we have access rights to create the output file (since we ignore any messages from ffmpeg)
        open(outFilePath,'w').close()

        # construct the first part of the command line for raw video input
        cmdline = [ 'ffmpeg',               # path to executable
                    '-v', 'quiet',          # be less verbose
                    '-y',                   # overwrite output file if it exists
                    '-f', 'rawvideo',       # input format: raw, uncompressed data stream
                    '-pix_fmt', 'rgba',     # input pixel format (3 channels plus dummy alpha channel, 8 bits each)
                    '-s', '{:1d}x{:1d}'.format(*shape), # frame size (pixels)
                    '-r', '{:1d}'.format(rate),         # frame rate (frames per second)
                    '-i', '-',              # the input comes from a pipe
                    '-an',                  # there is no audio
                    '-vcodec', 'mpeg4',     # output encoding
                    outFilePath ]

        # launch ffmpeg; pipe the input from this process and pipe any messages to the null device
        self._p = subprocess.Popen(cmdline, stdin=subprocess.PIPE)
                                   #stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    # This function writes a byte sequence containing the pixel data for a single movie frame to \c ffmpeg.
    # The frame must be specified as a pts.visual.rgbimage.RGBImage instance.
    # The shape of the image must match the shape specified when constructing the movie.
    def addFrame(self, image):
        assert image.shape() == self._shape
        self._p.stdin.write(image.bytesArray())

   # This function flushes all movie information to the file.
   # The function \em must be called for the movie file to be playable.
    def close(self):
        if self._p:
            self._p.stdin.close()
            self._p.wait()
            self._p = None

# -----------------------------------------------------------------
