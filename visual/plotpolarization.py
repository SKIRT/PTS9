#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.plotpolarization Plot polarization maps from SKIRT simulation output
#
# The function in this module creates PDF polarization maps for the output of polarization-enabled SKIRT simulations
# that include instruments generating surface brightness frames or data cubes (*_total.fits and *_stokes*.fits files).
#

# -----------------------------------------------------------------

import logging
import matplotlib.pyplot as plt
import matplotlib.lines
import matplotlib.colors
import numpy as np
import pts.simulation as sm
import pts.utils as ut

# -----------------------------------------------------------------

## This function creates PDF polarization maps for the output of polarization-enabled SKIRT simulations
# that include instruments generating surface brightness frames or data cubes (*.fits files).
# Specifically, the function accepts a sequence of Simulation and/or Instrument instances (or a single instance of
# either of these types), and it creates polarization maps for each of the instruments that actually produced
# both *_total.fits and *_stokes*.fits files. Other instruments are silently ignored.
#
# By default, the figures are saved in the output directory of the corresponding instrument, using a name starting
# with the simulation prefix and instrument name, and ending with ".pdf". The output directory can be overridden.
# In interactive mode (see the pts.utils.interactive() function), the figures are not saved and are left open
# for display in notebooks.
#
# The function takes the following arguments:
#   - simulation: a sequence of Simulation and/or Instrument instances, or a single instance of either of these types.
#   - plotLinMap: whether to plot a linear polarization map (binned degree and angle); True or False.
#   - plotDegMap: whether to plot a linear polarization degree map; True or False.
#   - plotDegAvg: whether to plot the y-axis averaged linear polarization degree for all x pixels; True or False.
#   - plotCirMap: whether to plot a circular polarization map; True, False or None for automatic.
#   - wavelength: an astropy quantity or a sequence of astropy quantities specifying the wavelength(s) for which to
#     create the plot, or the string 'all' to create a plot for all wavelengths in the data cube, or None to use the
#     first frame in the data cube.
#   - binSize: the number of pixels in each bin, in horizontal and vertical directions.
#   - polScale: (degree, length) to set the scale of the polarization degree; use Nones for automatic scale.
#   - decades: number of surface brightness decades (dex) shown in the background image; default is 5.
#   - outDirPath: string or Pathlib.Path object that specifies (overrides) the output directory.
#   - figSize: the horizontal and vertical size of the output figure in inch; default is 6x6 inch.
#   - interactive: whether to leave figures open (True) or save them to file (False).
#
def plotPolarization(simulation, *, plotLinMap=True, plotDegMap=False, plotDegAvg=False, plotCirMap=False,
                     wavelength=None, binSize=(7,7), polScale=(None,None), decades=5,
                     outDirPath=None, figSize=(8, 6), interactive=None):

    # loop over all applicable instruments
    for instrument, filepath in sm.instrumentOutFilePaths(simulation, "stokesQ.fits"):
        # form the simulation/instrument name
        insname = "{}_{}".format(instrument.prefix(), instrument.name())

        # get the file paths for the frames/data cubes
        filepathI = instrument.outFilePaths("total.fits")[0]
        filepathQ = instrument.outFilePaths("stokesQ.fits")[0]
        filepathU = instrument.outFilePaths("stokesU.fits")[0]
        filepathV = instrument.outFilePaths("stokesV.fits")[0]

        # load datacubes with shape (nx, ny, nlambda)
        Is = sm.loadFits(filepathI)
        Qs = sm.loadFits(filepathQ)
        Us = sm.loadFits(filepathU)
        Vs = sm.loadFits(filepathV)

        # load the axes grids (assuming all files have the same axes)
        xgrid, ygrid, wavegrid = sm.getFitsAxes(filepathI)
        xmin = xgrid[0].value
        xmax = xgrid[-1].value
        ymin = ygrid[0].value
        ymax = ygrid[-1].value
        extent = (xmin, xmax, ymin, ymax)

        # determine binning configuration
        binX = binSize[0]
        orLenX = Is.shape[0]
        dropX =  orLenX % binX
        startX = dropX//2
        binY = binSize[1]
        orLenY = Is.shape[1]
        dropY = orLenY % binY
        startY = dropY//2

        # construct arrays with central bin positions scaled to actual axes values
        posX = np.arange(startX - 0.5 + binX / 2.0, orLenX - dropX + startX - 0.5, binX)
        posY = np.arange(startY - 0.5 + binY / 2.0, orLenY - dropY + startY - 0.5, binY)
        posX = xmin + posX * (xmax-xmin)/orLenX
        posY = ymin + posY * (ymax-ymin)/orLenY

        # determine the appropriate wavelength index or indices
        if wavelength is None:
            indices = [0]
        elif wavelength == 'all':
            indices = range(Is.shape[2])
        else:
            if not isinstance(wavelength, (list,tuple)): wavelength = [wavelength]
            indices = instrument.wavelengthIndices(wavelength)

        # loop over all requested wavelength indices
        for index in indices:
            wave = wavegrid[index]
            wavename = "{:09.4f}".format(wave.to_value(sm.unit("micron")))
            wavelatex = r"$\lambda={:.4g}$".format(wave.value) + sm.latexForUnit(wave)

            # extract the corresponding frame, and transpose to (y,x) style for compatibility with legacy code
            I = Is[:,:,index].T.value
            Q = Qs[:,:,index].T.value
            U = Us[:,:,index].T.value
            V = Vs[:,:,index].T.value

            # if plotCirMap is None, we'll detect circular polarization for this frame during binning
            plotCirMapThis = plotCirMap

            # perform the actual binning
            binnedI = np.zeros((len(posY),len(posX)))
            binnedQ = np.zeros((len(posY),len(posX)))
            binnedU = np.zeros((len(posY),len(posX)))
            binnedV = np.zeros((len(posY),len(posX)))
            for x in range(len(posX)):
                for y in range(len(posY)):
                    binnedI[y,x] = np.sum(I[startY+binY*y : startY+binY*(y+1) , startX+binX*x : startX+binX*(x+1)])
                    binnedQ[y,x] = np.sum(Q[startY+binY*y : startY+binY*(y+1) , startX+binX*x : startX+binX*(x+1)])
                    binnedU[y,x] = np.sum(U[startY+binY*y : startY+binY*(y+1) , startX+binX*x : startX+binX*(x+1)])
                    binnedV[y,x] = np.sum(V[startY+binY*y : startY+binY*(y+1) , startX+binX*x : startX+binX*(x+1)])
                    # we only need circular polarization relative to total Intensity
                    binnedV[y,x] /= binnedI[y,x]
                    # if needed, detect circular polarization
                    if plotCirMapThis is None:
                        stdV = np.nanstd(V[startY+binY*y : startY+binY*(y+1) , startX+binX*x : startX+binX*(x+1)]/
                                         I[startY+binY*y : startY+binY*(y+1) , startX+binX*x : startX+binX*(x+1)])
                        if stdV == 0: stdV = 1 # if we only have one polarized pixel in the bins, stdV is undefined
                        if np.abs(binnedV[y,x])/stdV > 2: plotCirMapThis = True

            # -----------------------------------------------------------------

            # plot a linear polarization map
            if plotLinMap:
                fig, ax = plt.subplots(ncols=1, nrows=1, figsize=figSize)

                # configure the axes
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)
                ax.set_xlabel("x" + sm.latexForUnit(xgrid), fontsize='large')
                ax.set_ylabel("y" + sm.latexForUnit(ygrid), fontsize='large')

                # determine intensity range for the background image, ignoring pixels with outrageously high flux
                Ib = I.copy()
                highmask = Ib > 1e6*np.nanmedian(np.unique(Ib))
                vmax = np.nanmax(Ib[~highmask])
                Ib[highmask] = vmax
                vmin = vmax / 10**decades

                # plot the background image and the corresponding color bar
                normalizer = matplotlib.colors.LogNorm(vmin, vmax)
                cmap = plt.get_cmap('PuRd')
                cmap.set_under('w')
                backPlot = ax.imshow(Ib, norm=normalizer, cmap=cmap, extent=extent,
                                     aspect='equal', interpolation='bicubic', origin='lower')
                cbarlabel = sm.latexForSpectralFlux(Is) + sm.latexForUnit(Is) + " @ " +  wavelatex
                plt.colorbar(backPlot, ax=ax).set_label(cbarlabel, fontsize='large')

                # compute the linear polarization degree
                degreeLD = np.sqrt(binnedQ**2 + binnedU**2)
                degreeLD[degreeLD>0] /= binnedI[degreeLD>0]

                # determine a characteristic 'high' degree of polarization in the frame
                # (this has to be done before degreeLD contains 'np.NaN')
                charDegree = np.percentile(degreeLD, 99.0)
                if not 0<charDegree<1:
                    charDegree = np.nanmax((np.nanmax(degreeLD), 0.0001))

                # remove pixels with minuscule polarization
                degreeLD[degreeLD<charDegree/50] = np.NaN

                # determine the scaling so that the longest arrows do not to overlap with neighboring arrows
                degreeScale, lengthScale = polScale
                if degreeScale is None:
                    degreeScale = _roundUp(charDegree, slack=0.0001)
                if lengthScale is None:
                    lengthScale = 1/(degreeScale*2.2)/max(float(len(posX))/figSize[0], float(len(posY))/figSize[1])
                key = "{:.3g}%".format(100 * degreeScale)

                # compute the polarization angle
                angle = 0.5 * np.arctan2(binnedU, binnedQ) # angle from North through East while looking at the sky

                # create the polarization vector arrays
                xPolarization = -degreeLD*np.sin(angle) #For angle = 0: North & x=0, For angle = 90deg: West & x=-1
                yPolarization =  degreeLD*np.cos(angle) #For angle = 0: North & y=1, For angle = 90deg: West & y=0

                # plot the vector field
                X,Y = np.meshgrid(posX, posY)
                quiverPlot = ax.quiver(X,Y, xPolarization, yPolarization, pivot='middle', units='inches',
                                        angles = 'xy', scale = 1/lengthScale, scale_units = 'inches', headwidth=0,
                                        headlength=1, headaxislength=1, minlength = 0.8, width = 0.02)
                ax.quiverkey(quiverPlot, 0.85, 0.02, degreeScale, key, coordinates='axes', labelpos='E')

                # if not in interactive mode, save the figure; otherwise leave it open
                if not ut.interactive(interactive):
                    saveFilePath = ut.savePath(filepath, ".pdf", outDirPath=outDirPath,
                                               outFileName = "{}_{}_pollinmap.pdf".format(insname, wavename))
                    plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
                    plt.close()
                    logging.info("Created {}".format(saveFilePath))

            # -----------------------------------------------------------------

            # plot a linear polarization degree map
            if plotDegMap:
                fig, ax = plt.subplots(ncols=1, nrows=1, figsize=figSize)

                # configure the axes
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)
                ax.set_xlabel("x" + sm.latexForUnit(xgrid), fontsize='large')
                ax.set_ylabel("y" + sm.latexForUnit(ygrid), fontsize='large')

                # calculate polarization degree for each pixel, in percent
                degreeHD = np.sqrt(Q**2 + U**2)
                degreeHD[degreeHD>0] /= I[degreeHD>0]
                degreeHD *= 100

                # set degree to zero for pixels with very low intensity
                Icutoff = np.nanmax(I[I < 1e6*np.nanmedian(np.unique(I))]) / 10**decades
                degreeHD[I < Icutoff] = 0

                # plot the image and the corresponding color bar
                normalizer = matplotlib.colors.Normalize(vmin=0, vmax=np.percentile(degreeHD, 99))
                backPlot = ax.imshow(degreeHD, norm=normalizer, cmap='plasma', extent=extent,
                                     aspect='equal', interpolation='bicubic', origin='lower')
                plt.colorbar(backPlot, ax=ax).set_label("Linear polarization degree (%)" + " @ " +  wavelatex,
                                                        fontsize='large')

                # if not in interactive mode, save the figure; otherwise leave it open
                if not ut.interactive(interactive):
                    saveFilePath = ut.savePath(filepath, ".pdf", outDirPath=outDirPath,
                                               outFileName = "{}_{}_poldegmap.pdf".format(insname, wavename))
                    plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
                    plt.close()
                    logging.info("Created {}".format(saveFilePath))

            # -----------------------------------------------------------------

            # plot the y-axis averaged linear polarization degree
            if plotDegAvg:
                # construct the plot
                fig, ax = plt.subplots(ncols=1, nrows=1, figsize=figSize)
                degreeHD = np.sqrt(np.average(Q, axis = 0)**2 + np.average(U, axis=0)**2)
                degreeHD /= np.average(I,axis = 0)
                ax.plot(xgrid.value, degreeHD*100)
                ax.set_xlim(xmin, xmax)
                ax.set_title("{}   {}".format(insname, wavelatex), fontsize='large')
                ax.set_xlabel("x" + sm.latexForUnit(xgrid), fontsize='large')
                ax.set_ylabel('Average linear polarization degree (%)', fontsize='large')

                # if not in interactive mode, save the figure; otherwise leave it open
                if not ut.interactive(interactive):
                    saveFilePath = ut.savePath(filepathI, ".pdf", outDirPath=outDirPath,
                                               outFileName = "{}_{}_poldegavg.pdf".format(insname, wavename))
                    plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
                    plt.close()
                    logging.info("Created {}".format(saveFilePath))

            # -----------------------------------------------------------------

            # plot a circular polarization map
            if plotCirMapThis:
                pass

# -----------------------------------------------------------------

## This helper function rounds up a given number (away from zero) to one significant digit, in a logarithmic fashion.
# The \em slack argument allows a number just larger than the nearest limit to be clamped to that limit.
def _roundUp(x, slack = 0.):
    x *= 1.-slack
    signX = np.sign(x)
    x = np.abs(x)
    dec = 10**(np.floor(np.log10(x)))
    mant = x/dec
    if mant <= 1: return signX*1*dec
    if mant <= 2: return signX*2*dec
    if mant <= 4: return signX*4*dec
    if mant <= 5: return signX*5*dec
    return signX*10*dec
