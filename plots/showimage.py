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
myshowimage.py
-------------
Display and Cut image file 

Usage: myshowimage.py --infile=<path> [<ibeg>] [<iend>] [<jbeg>] [<jend>]
Options:
-h --help           Show this screen.
--infile PATH       File to be cut
--ibeg VALUE        Ligne numbers bounded the cutting zone [default: 0]
--iend VALUE        Ligne numbers bounded the cutting zone [default: nlign]
--jbeg VALUE        Column numbers bounded the cutting zone [default: 0]
--jend VALUE        Column numbers bounded the cutting zone [default: ncol]
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

# Initialize a matplotlib figure
fig, ax = plt.subplots(1,figsize=(10,11))

# read arguments
arguments = docopt.docopt(__doc__)
infile = arguments["--infile"]
# Open phiset (image)
ds = gdal.Open(infile, gdal.GA_ReadOnly)

# Ok, assume this is a roipac image, so we can rely on the file extension
# to know how to display the phi (in the real world, we would probably write
# different programs)
ds_extension = os.path.splitext(sys.argv[1])[1]

# Get the band that have the phi we want
if ds_extension == ".unw":
    phase_band = ds.GetRasterBand(2)
    amp_band = ds.GetRasterBand(1)
else:
	phase_band = ds.GetRasterBand(1)

# Attributes
print("> Driver:   ", ds.GetDriver().ShortName)
print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
print("> Datatype: ", gdal.GetDataTypeName(phase_band.DataType))

if arguments["<ibeg>"] ==  None:
    ibeg = 0
else:
    ibeg = int(arguments["<ibeg>"])
if arguments["<iend>"] ==  None:
    iend = ds.RasterYSize
else:
    iend = int(arguments["<iend>"])
if arguments["<jbeg>"] ==  None:
    jbeg = 0
else:
    jbeg = int(arguments["<jbeg>"])
if arguments["<jend>"] ==  None:
    jend = ds.RasterXSize
else:
    jend = int(arguments["<jend>"])

# Read phi in numpy array, resampled by a factor of 4
#   Resampling is nearest neighbour with that function: hardly acceptable,
#   but easy for a simple demo!
if ds_extension == ".unw":
	phi = phase_band.ReadAsArray(0, 0,
                           ds.RasterXSize, ds.RasterYSize,
                           ds.RasterXSize, ds.RasterYSize)
	amp = amp_band.ReadAsArray(0, 0,
                           ds.RasterXSize, ds.RasterYSize,
                           ds.RasterXSize, ds.RasterYSize)
	
	# cut image
	cutphi = as_strided(phi[ibeg:iend,jbeg:jend])
	cutamp = as_strided(amp[ibeg:iend,jbeg:jend])

else:
	phase = phase_band.ReadAsArray(0, 0,
                           ds.RasterXSize, ds.RasterYSize,
                           ds.RasterXSize, ds.RasterYSize)
	
	cutphase = as_strided(phase[ibeg:iend,jbeg:jend])	

#phi = np.nan_to_num(phi)
#amp = np.nan_to_num(amp)
#kk = np.flatnonzero(phi==0)
#phi[kk] = np.float('NaN')


if ds_extension == ".slc":
    # SLC, display amplitude
    cutphi = np.absolute(phi)
    cax = ax.imshow(cutphi, cm.Greys_r, vmax=np.percentile(cutphi, 95))
elif ds_extension in [".int", ".flat"]:
    # Wrapped interferogram, display computed phase
    cutamp = np.absolute(cutphase)
    cutphi = np.angle(cutphase)
    # hax = ax.imshow(cutamp, cm.Greys_r,vmax=np.percentile(cutamp, 90))
    hax = ax.imshow(cutamp, cm.Greys,vmax=1,vmin=0.)
    # cax = ax.imshow(cutphi, cm.gist_rainbow, interpolation='bicubic',alpha=0.8)

elif ds_extension == ".hgt":
    # DEM in radar geometry
    cax = ax.imshow(cutphi, cm.Greys_r, vmax=np.percentile(cutphi, 95))
elif ds_extension == ".unw":
    # Unwrapped inteferogram
    # hax = ax.imshow(cutamp, cm.Greys,vmax=1,vmin=0.)
    #hax = ax.imshow(cutamp, cm.Greys,vmax=np.percentile(cutamp, 95))
    vmax=np.nanmean(cutphi)+2*np.nanstd(cutphi)
    vmin=np.nanmean(cutphi)-2*np.nanstd(cutphi)
    vmax=8
    vmin=-8
    cutwrap = np.mod(cutphi,2*np.pi)-np.pi
    cax = ax.imshow(cutphi, cm.gist_rainbow, interpolation='bicubic',vmax=vmax,vmin=vmin,alpha=0.5)
    # cax = ax.imshow(cutwrap, cm.gist_rainbow, interpolation='bicubic',alpha=0.7,vmax=np.pi,vmin=-np.pi)
elif ds_extension ==  ".trans":
    # Geocoding lookup table
    cax = ax.imshow(cutphi)

setp( ax.get_xticklabels(), visible=False)
cbar = fig.colorbar(cax, orientation='vertical',aspect=9,fraction=0.02,pad=0.06)
# cbar = fig.colorbar(hax, orientation='vertical',aspect=10,fraction=0.02,pad=0.01)

# Close dataset (if this was a new image, it would be written at this moment)
# This will free memory for the display stuff
del ds

# Display the data
fig.canvas.set_window_title(sys.argv[1])
##ax.set_rasterized(True)
fig.savefig('{}.eps'.format(infile), format='EPS',dpi=150)
plt.show()
