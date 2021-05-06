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
add_rmg.py
-------------
Add, substract, and plot RMG (BIL) unwrapped files.

Usage: add_rmg.py --infile=<path> --outfile=<path> [--add=<path>] [--remove=<path>] [--plot=<yes/no>]

Options:
-h --help           Show this screen.
--infile PATH       Input interferogram
--outfile PATH      Output interferogram
--add PATH          Model to be added [Optional]
--remove PATH       Model to be removed [Optional]
--plot VALUE        Display [default:no]
"""

import gdal
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import *
import docopt

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning) 

# read arguments
arguments = docopt.docopt(__doc__)
infile = arguments["--infile"]
outfile = arguments["--outfile"]
if arguments["--add"] == None:
    add = 'no'
else:
    add = arguments["--add"]
if arguments["--remove"] == None:
    remove = 'no'
else:
    remove = arguments["--remove"]
if arguments["--add"] == None and arguments["--remove"] == None:
    print("Nothing to do, add or remove a file!")
    sys.exit()

if arguments["--plot"] == None:
    plot = 'no'
else:
    plot = arguments["--plot"]

# print(infile)
# print(add)
# print remove

gdal.UseExceptions()
# Open int
print(infile)
ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])
ds_band1 = ds.GetRasterBand(1)
ds_band2 = ds.GetRasterBand(2)

# # Attributes
# print("> Driver:   ", ds.GetDriver().ShortName)
# print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
# print("> Datatype: ", gdal.GetDataTypeName(ds_band2.DataType))
# print

amp = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
phi = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
nlign,ncol = ds.RasterYSize, ds.RasterXSize
temp = np.copy(phi)
kk = np.nonzero(temp==0.0)
temp[kk]=float('NaN')

del ds

if remove is not "no":
    ds = gdal.OpenEx(remove, allowed_drivers=["ROI_PAC"])
    # print("> Driver:   ", ds.GetDriver().ShortName)
    # print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
    # print("> Datatype: ", gdal.GetDataTypeName(ds_band2.DataType))
    # print
    ds_band2 = ds.GetRasterBand(2)
    remphi = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:nlign,:ncol]
    temp -= remphi
    del ds

#Open new model
if add is not "no":
    ds = gdal.OpenEx(add, allowed_drivers=["ROI_PAC"])
    # print("> Driver:   ", ds.GetDriver().ShortName)
    # print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
    # print("> Datatype: ", gdal.GetDataTypeName(ds_band2.DataType))
    # print
    ds_band2 = ds.GetRasterBand(2)
    addphi = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:nlign,:ncol]
    temp += addphi
    del ds

kk = np.nonzero(temp > -10000)
out = np.zeros((nlign,ncol))
out[kk]=temp[kk]

# Create new file
drv = gdal.GetDriverByName("roi_pac")
# Creer l'image, avec le chemin, largeur, hauteur
dst_ds = drv.Create(outfile, ncol, nlign, 2, gdal.GDT_Float32)
# Sortir les 2 bandes (float)
dst_band1 = dst_ds.GetRasterBand(1)
dst_band2 = dst_ds.GetRasterBand(2)
# Write new data
dst_band1.WriteArray(amp,0,0)
dst_band2.WriteArray(out,0,0)
# close image
del dst_ds

# plot
if plot=="yes":
  if remove is not "no" and add is not "no":
    vmax = np.percentile(phi, 98)
    fig = plt.figure(0,figsize=(8,9))
    # data
    ax = fig.add_subplot(1,4,1)
    cax = ax.imshow(phi,cmap=cm.gist_rainbow,vmax=vmax)
    ax.set_title('infile')
    setp( ax.get_xticklabels(), visible=False)
    # remove
    ax = fig.add_subplot(1,4,2)
    cax = ax.imshow(remphi,cmap=cm.gist_rainbow,vmax=vmax)
    ax.set_title('remove')
    setp( ax.get_xticklabels(), visible=False)
    # model
    ax = fig.add_subplot(1,4,3)
    cax = ax.imshow(addphi,cmap=cm.gist_rainbow,vmax=vmax)
    ax.set_title('add')
    setp( ax.get_xticklabels(), visible=False)
    # output
    #vmax = np.percentile(out, 98)
    ax = fig.add_subplot(1,4,4)
    cax = ax.imshow(out,cmap=cm.gist_rainbow,vmax=vmax)
    ax.set_title('outfile')
    setp( ax.get_xticklabels(), visible=False)
  
  elif add is not "no":
    vmax = np.percentile(phi, 98)
    fig = plt.figure(0,figsize=(8,9))
    # data
    ax = fig.add_subplot(1,3,1)
    cax = ax.imshow(phi,cmap=cm.gist_rainbow,vmax=vmax)
    ax.set_title('infile')
    setp( ax.get_xticklabels(), visible=False)
    # model
    ax = fig.add_subplot(1,3,2)
    cax = ax.imshow(addphi,cmap=cm.gist_rainbow,vmax=vmax)
    ax.set_title('add')
    setp( ax.get_xticklabels(), visible=False)
    # output
    #vmax = np.percentile(out, 98)
    ax = fig.add_subplot(1,3,3)
    cax = ax.imshow(out,cmap=cm.gist_rainbow,vmax=vmax)
    ax.set_title('outfile')
    setp( ax.get_xticklabels(), visible=False)

  else:
    vmax = np.percentile(phi, 98)
    fig = plt.figure(0,figsize=(8,9))
    # data
    ax = fig.add_subplot(1,3,1)
    cax = ax.imshow(phi,cmap=cm.gist_rainbow,vmax=vmax)
    ax.set_title('infile')
    setp( ax.get_xticklabels(), visible=False)
    # remove
    ax = fig.add_subplot(1,3,2)
    cax = ax.imshow(remphi,cmap=cm.gist_rainbow,vmax=vmax)
    ax.set_title('remove')
    setp( ax.get_xticklabels(), visible=False)
    # output
    #vmax = np.percentile(out, 98)
    ax = fig.add_subplot(1,3,3)
    cax = ax.imshow(out,cmap=cm.gist_rainbow,vmax=vmax)
    ax.set_title('outfile')
    setp( ax.get_xticklabels(), visible=False)

plt.show()
