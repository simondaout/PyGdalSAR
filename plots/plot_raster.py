#!/usr/bin/env python2
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Mathieu Volat 
############################################

"""\
plot\_raster.py
-------------
Display and Cut image file 

Usage: plot\_raster.py --infile=<path> [--cpt=<values>] [--wrap=<values>] [<ibeg>] [<iend>] [<jbeg>] [<jend>] [--format=<value>]
Options:
-h --help           Show this screen.
--infile PATH       File to be cut
--ibeg VALUE        Ligne numbers bounded the cutting zone [default: 0]
--iend VALUE        Ligne numbers bounded the cutting zone [default: nlign]
--jbeg VALUE        Column numbers bounded the cutting zone [default: 0]
--jend VALUE        Column numbers bounded the cutting zone [default: ncol]
--cpt               Indicate colorscale for phase
--wrap  VALUE       Wrapped phase between value for unwrapped files [default: no]
--format VALUE        Format input files: ROI_PAC, GAMMA, GTIFF [default: ROI_PAC]
"""

from __future__ import print_function
import os, sys

# docopt (command line parser)
from nsbas import docopt

import numpy as np
from numpy.lib.stride_tricks import as_strided

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import setp
from osgeo import gdal
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gamma as gm

# Initialize a matplotlib figure
fig, ax = plt.subplots(1,figsize=(8,10))

# read arguments
arguments = docopt.docopt(__doc__)
infile = arguments["--infile"]
# Open phiset (image)
ds = gdal.Open(infile, gdal.GA_ReadOnly)

if arguments["--cpt"] is  None:
    # cmap=cm.jet 
    cmap=cm.rainbow
else:
    cmap=arguments["--cpt"]

if arguments["--format"] ==  None:
    sformat = 'ROI_PAC'
else:
    sformat = arguments["--format"]

ds_extension = os.path.splitext(sys.argv[1])[1]

if sformat == "ROI_PAC": 
    if (ds_extension == ".unw" or ds_extension ==".hgt"):
        phase_band = ds.GetRasterBand(2)
        amp_band = ds.GetRasterBand(1)
    else:
        phase_band = ds.GetRasterBand(1)

    # Attributes
    print("> Driver:   ", ds.GetDriver().ShortName)
    print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
    print("> Datatype: ", gdal.GetDataTypeName(phase_band.DataType))
    nlines, ncols = ds.RasterYSize, ds.RasterXSize

elif sformat == "GTIFF":
    phase_band = ds.GetRasterBand(1)
    # Attributes
    print("> Driver:   ", ds.GetDriver().ShortName)
    print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
    print("> Datatype: ", gdal.GetDataTypeName(phase_band.DataType))
    nlines, ncols = ds.RasterYSize, ds.RasterXSize

elif sformat == 'GAMMA':
    par_file =  '.par'
    nlines,ncols = gm.readpar()
    phi = gm.readgamma(infile)

if arguments["<ibeg>"] ==  None:
    ibeg = 0
else:
    ibeg = int(arguments["<ibeg>"])
if arguments["<iend>"] ==  None:
    iend = nlines
else:
    iend = int(arguments["<iend>"])
if arguments["<jbeg>"] ==  None:
    jbeg = 0
else:
    jbeg = int(arguments["<jbeg>"])
if arguments["<jend>"] ==  None:
    jend = ncols
else:
    jend = int(arguments["<jend>"])


if sformat == "ROI_PAC":
    if (ds_extension == ".unw" or ds_extension ==".hgt"):
        phi = phase_band.ReadAsArray(0, 0,
                               ds.RasterXSize, ds.RasterYSize,
                               ds.RasterXSize, ds.RasterYSize)
        amp = amp_band.ReadAsArray(0, 0,
                               ds.RasterXSize, ds.RasterYSize,
                               ds.RasterXSize, ds.RasterYSize)
        
        # cut image
        cutphi = as_strided(phi[ibeg:iend,jbeg:jend])
        cutamp = as_strided(amp[ibeg:iend,jbeg:jend])

        # Unwrapped inteferogram
        hax = ax.imshow(cutamp, cm.Greys,vmax=np.percentile(cutamp, 95))

        if arguments["--wrap"] is not None:
            cutphi =  np.mod(cutphi+float(arguments["--wrap"]),2*float(arguments["--wrap"]))-float(arguments["--wrap"])
            vmax=float(arguments["--wrap"])
            vmin = -vmax 
        else:
            vmax=np.nanmean(cutphi)+2*np.nanstd(cutphi)
            vmin=np.nanmean(cutphi)-2*np.nanstd(cutphi)

        cax = ax.imshow(cutphi, cmap, interpolation=None,vmax=vmax,vmin=vmin,alpha=0.6)

    else:
        phi = phase_band.ReadAsArray(0, 0,
                               ds.RasterXSize, ds.RasterYSize,
                               ds.RasterXSize, ds.RasterYSize)
        
        cutphi = as_strided(phi[ibeg:iend,jbeg:jend]) 
      
        if ds_extension == ".slc":
            # SLC, display amplitude
            cutphi = np.absolute(cutphi)
            cax = ax.imshow(cutphi, cm.Greys_r, vmax=np.percentile(cutphi, 95))
        elif ds_extension in [".int", ".flat"]:

            # Wrapped interferogram, display computed phase
            cutamp = np.absolute(cutphi)
            cutphi = np.angle(cutphi)
            hax = ax.imshow(cutamp, cm.Greys_r)
            cax = ax.imshow(cutphi, cmap,alpha=0.6)
            # hax = ax.imshow(cutamp, cm.Greys,vmax=1,vmin=0.)
        elif ds_extension ==  ".trans":
            # Geocoding lookup table
            cax = ax.imshow(cutphi)

elif sformat == "GTIFF":
    phi = phase_band.ReadAsArray(0, 0,
           ds.RasterXSize, ds.RasterYSize,
           ds.RasterXSize, ds.RasterYSize)

    cutphi = as_strided(phi[ibeg:iend,jbeg:jend]) 
    cax = ax.imshow(cutphi, cmap,alpha=0.6, interpolation=None)
elif sformat == "GAMMA":
    cutphi = as_strided(phi[ibeg:iend,jbeg:jend]) 
    cax = ax.imshow(cutphi, cmap,alpha=0.6, interpolation=None)

divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

try:
    del ds
except:
    pass

# Display the data
fig.canvas.set_window_title(sys.argv[1])
##ax.set_rasterized(True)
fig.savefig('{}.eps'.format(infile), format='EPS',dpi=150)
plt.show()
