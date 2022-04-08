#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Simon DAOUT 
############################################

"""\
cut_int.py
-------------
Cut int file 

Usage: cut_int.py --infile=<path> --outfile=<path> --plot=<yes/no> [<ibeg>] [<iend>] [<jbeg>] [<jend>]

Options:
-h --help           Show this screen.
--infile PATH       File to be cut
--outfile PATH      Cutting file
--plot VALUE        Display interferograms
--ibeg VALUE        Ligne numbers bounded the cutting zone [default: 0]
--iend VALUE        Ligne numbers bounded the cutting zone [default: nlign]
--jbeg VALUE        Column numbers bounded the cutting zone [default: 0]
--jend VALUE        Column numbers bounded the cutting zone [default: ncol]
"""

# lib
import numpy as np
from numpy.lib.stride_tricks import as_strided

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import *

import os, sys

from osgeo import gdal
import time

# docopt (command line parser)
import docopt

# read arguments
arguments = docopt.docopt(__doc__)
infile = arguments["--infile"]
outfile = arguments["--outfile"]
if arguments["--plot"] ==  None:
    plot = 'no'
else:
    plot = arguments["--plot"]

gdal.UseExceptions()
# Open phiset (image)
ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])

# Ok, assume this is a roipac image, so we can rely on the file extension
# to know how to display the phi (in the real world, we would probably write
# different programs)
ds_extension = os.path.splitext(infile)[1]

# Get the band that have the phi we want
ds_band1 = ds.GetRasterBand(1)

# Attributes
print("> Driver:   ", ds.GetDriver().ShortName)
print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
print("> Datatype: ", gdal.GetDataTypeName(ds_band1.DataType))

# Read phi in numpy array, resampled by a factor of 4
#   Resampling is nearest neighbour with that function: hardly acceptable,
#   but easy for a simple demo!
pha = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
amp = np.absolute(pha)
phi = np.angle(pha)
nlign,ncol = ds.RasterYSize, ds.RasterXSize

if arguments["<ibeg>"] ==  None:
    ibeg = 0
else:
    ibeg = int(arguments["<ibeg>"])
if arguments["<iend>"] ==  None:
    iend = nlign
else:
    iend = int(arguments["<iend>"])
if arguments["<jbeg>"] ==  None:
    jbeg = 0
else:
    jbeg = int(arguments["<jbeg>"])
if arguments["<jend>"] ==  None:
    jend = ncol
else:
    jend = int(arguments["<iend>"])

# cut
cutamp = np.copy(amp)
cutamp[0:ibeg,:] = 0.
cutamp[iend:nlign,:] = 0.
cutamp[:,0:jbeg] = 0.
cutamp[:,jend:ncol] = 0.
    
# create new GDAL imaga
drv = gdal.GetDriverByName("roi_pac")
dst_ds = drv.CreateCopy(outfile, ds)
dst_band1 = dst_ds.GetRasterBand(1)

# write in the output file 
cutout = cutamp * exp(1j*phi)
dst_band1.WriteArray(cutout,0,0)

# close image
del dst_ds
del ds

if plot=="yes":
    # Plot
    fig = plt.figure(0,figsize=(8,9))
    ax = fig.add_subplot(1,2,1)
    hax = ax.imshow(amp,cmap=cm.Greys_r,vmin=0, vmax=1)
    
    # mask for plot 
    phi[phi==0] = float('NaN')
    masked_array = np.ma.array(phi, mask=np.isnan(phi))
    cax = ax.imshow(masked_array,cmap=cm.gist_rainbow,alpha=.8,interpolation='nearest')
    ax.set_title('In Phase')

    ax = fig.add_subplot(1,2,2)
    hax = ax.imshow(np.absolute(cutout),cmap=cm.Greys_r,vmin=0, vmax=1)
    phi = np.angle(cutout)
    phi[cutamp==0] = float('NaN')
    masked_array = np.ma.array(phi, mask=np.isnan(phi))
    cax = ax.imshow(masked_array,cmap=cm.gist_rainbow,alpha=.8,interpolation='nearest')
    ax.set_title('Cut Phase')

    #cbar = fig.colorbar(cax, orientation='vertical',aspect=5)

    # Display the phi
    fig.canvas.set_window_title(sys.argv[1])
    plt.show()
