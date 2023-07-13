#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Simon DAOUT (Oxford)
############################################

"""\
cut_unw.py
-------------
Cut unw file 

Usage: cut_unw.py --infile=<path> --outfile=<path> --plot=<yes/no> [<ibeg>] [<iend>] [<jbeg>] [<jend>]

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
plot = arguments["--plot"]

gdal.UseExceptions()
# Open dataset (image)
ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])

# Ok, assume this is a roipac image, so we can rely on the file extension
# to know how to display the data (in the real world, we would probably write
# different programs)
ds_extension = os.path.splitext(infile)[1]

# Get the band that have the data we want
ds_band1 = ds.GetRasterBand(1)
ds_band2 = ds.GetRasterBand(2)

# Attributes
print("> Driver:   ", ds.GetDriver().ShortName)
print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
print("> Datatype: ", gdal.GetDataTypeName(ds_band2.DataType))

# Read data in numpy array, resampled by a factor of 4
#   Resampling is nearest neighbour with that function: hardly acceptable,
#   but easy for a simple demo!
amp = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
data = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
nlign,ncol = ds.RasterYSize, ds.RasterXSize

# Close dataset (if this was a new image, it would be written at this moment)
# This will free memory for the display stuff
try:
    del ds
    print('... waiting ...')
    time.sleep(0.5)
except:
    sys.exit(1)

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
cutamp[0:ibeg] = 0
cutamp[iend:nlign] = 0
cutamp[0:jbeg] = 0
cutamp[jend:ncol] = 0
    
cutdata = np.copy(data)
cutdata[0:ibeg] = 0
cutdata[iend:nlign] = 0
cutdata[0:jbeg] = 0
cutdata[jend:ncol] = 0

# create new GDAL imaga
# recupérer driver ROI_PAC
drv = gdal.GetDriverByName("roi_pac")
# Creer l'image, avec le chemin, largeur, hauteur
dst_ds = drv.Create(outfile, ncol, nlign, 2, gdal.GDT_Float32)
# Sortir les 2 bandes (float)
dst_band1 = dst_ds.GetRasterBand(1)
dst_band2 = dst_ds.GetRasterBand(2)

# write in the output file 
# dst_band2.WriteArray(array2, 0, y) pour ajouter une ligne
# dst_band1.WriteArray(data, 0, 0) pour ajouter un tableau 2d
dst_band1.WriteArray(cutamp,0,0)
dst_band2.WriteArray(cutdata,0,0)

# close image
del dst_ds

# save output
#fid2 = open(outfile,'wb')
#cutdata.flatten().astype('f4').tofile(fid2)

if plot=="yes":
    #ax.imshow(data, cm.gist_rainbow, vmax=np.percentile(data, 98))
    vmax = np.percentile(data, 98)
    # Plot
    fig = plt.figure(0,figsize=(8,9))
    ax = fig.add_subplot(1,2,1)
    cax = ax.imshow(data,cmap=cm.gist_rainbow,vmax=vmax)
    ax.set_title(infile)
    setp( ax.get_xticklabels(), visible=False)

    ax = fig.add_subplot(1,2,2)
    cax = ax.imshow(cutdata,cmap=cm.gist_rainbow,vmax=vmax)
    ax.set_title(outfile)
    setp( ax.get_xticklabels(), visible=False)

    cbar = fig.colorbar(cax, orientation='vertical',aspect=5)


    # Display the data
    fig.canvas.set_window_title(sys.argv[1])
    plt.show()
