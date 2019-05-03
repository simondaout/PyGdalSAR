#!/usr/bin/env python2
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Mathieu Volat 
#                 Simon Daout
############################################

"""\
plot\_raster.py
-------------
Display and Cut image file (.unw/.int/.r4/.tiff)

Usage: plot\_raster.py --infile=<path> [--cpt=<values>] [<ibeg>] [<iend>] [<jbeg>] [<jend>] \
[--format=<value>] [--lectfile=<value>] [--rad2mm=<value>] [--title=<value>] [--wrap=<values>] 
       plot\_raster.py --infile=<path> [--cpt=<values>] [<ibeg>] [<iend>] [<jbeg>] [<jend>] \
[--format=<value>] [--parfile=<path>] [--lectfile=<value>] [--rad2mm=<value>] [--title=<value>] [--wrap=<values>] [--vmin=<value>] [--vmax=<value>]


Options:
-h --help             Show this screen.
--infile=<file>       Raster to be displayed 
--ibeg=<value>        Ligne numbers bounded the cutting zone [default: 0]
--iend=<value>        Ligne numbers bounded the cutting zone [default: nlign]
--jbeg=<value>        Column numbers bounded the cutting zone [default: 0]
--jend=<value>        Column numbers bounded the cutting zone [default: ncol]
--format=<value>      Format input files: ROI_PAC, GAMMA, GTIFF [default: ROI_PAC]
--cpt==<value>        Indicate colorscale for phase
--wrap=<value>        Wrapped phase between value for unwrapped files 
--lectfile=<file>     Path of the lect.in file for r4 format
--parfile=<file>      Path of the .par file of GAMMA
--rad2mm=<value>      Convert data [default: 1]
--tile=<value>        Title plot 
--vmax                Max colorscale [default: 98th percentile]
--vmin                Min colorscale [default: 2th percentile]
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

if arguments["--lectfile"] ==  None:
    lecfile = "lect.in"
else:
    lecfile = arguments["--lectfile"]

if arguments["--rad2mm"] ==  None:
        rad2mm = 1
else:
        rad2mm = np.float(arguments["--rad2mm"])

ds_extension = os.path.splitext(infile)[1]

if sformat == "ROI_PAC": 
    if (ds_extension == ".unw" or ds_extension ==".hgt"):
        phase_band = ds.GetRasterBand(2)
        amp_band = ds.GetRasterBand(1)
        # Attributes
        print("> Driver:   ", ds.GetDriver().ShortName)
        print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
        print("> Datatype: ", gdal.GetDataTypeName(phase_band.DataType))
        nlines, ncols = ds.RasterYSize, ds.RasterXSize
    elif (ds_extension == ".r4" or ds_extension == ""):
        fid = open(infile, 'r')
        ncols, nlines = map(int, open(lecfile).readline().split(None, 2)[0:2])
        phi = np.fromfile(fid,dtype=np.float32)[:nlines*ncols].reshape((nlines,ncols))
        print("> Driver:   REAL4  band file")
        print("> Size:     ", ncols,'x',nlines,'x')
        print("> Datatype: FLOAT32")
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
    import gamma as gm
    if arguments["--parfile"] ==  None:
        par_file =  arguments["--parfile"]
        nlines,ncols = gm.readpar(par=par_file)
        phi = gm.readgamma(infile,par=par_file)
    else:
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
        cutphi = as_strided(phi[ibeg:iend,jbeg:jend])*rad2mm
        cutamp = as_strided(amp[ibeg:iend,jbeg:jend])

        # Unwrapped inteferogram
        hax = ax.imshow(cutamp, cm.Greys,vmax=np.percentile(cutamp, 95))

    elif (ds_extension == ".r4" or ds_extension == ""):

        cutphi = as_strided(phi[ibeg:iend,jbeg:jend])*rad2mm 

    else:
        phi = phase_band.ReadAsArray(0, 0,
                               ds.RasterXSize, ds.RasterYSize,
                               ds.RasterXSize, ds.RasterYSize)
        
        cutphi = as_strided(phi[ibeg:iend,jbeg:jend])*rad2mm 
        
        if ds_extension == ".slc":
            # SLC, display amplitude
            cutphi = np.absolute(cutphi)
        elif ds_extension in [".int", ".flat"]:
            # Wrapped interferogram, display computed phase
            cutamp = np.absolute(cutphi)
            cutphi = np.angle(cutphi)
            hax = ax.imshow(cutamp, cm.Greys_r)

elif sformat == "GTIFF":
    phi = phase_band.ReadAsArray(0, 0,
           ds.RasterXSize, ds.RasterYSize,
           ds.RasterXSize, ds.RasterYSize)
    cutphi = as_strided(phi[ibeg:iend,jbeg:jend])*rad2mm

elif sformat == "GAMMA":
    cutphi = as_strided(phi[ibeg:iend,jbeg:jend])*rad2mm


if arguments["--wrap"] is not None:
    cutphi =  np.mod(cutphi+float(arguments["--wrap"]),2*float(arguments["--wrap"]))-float(arguments["--wrap"])
    vmax=float(arguments["--wrap"])
    vmin = -vmax
elif (arguments["--vmax"] is not None) or (arguments["--vmin"] is not None):
    if arguments["--vmax"] is not  None:
        vmax = np.float(arguments["--vmax"])
    if arguments["--vmin"] is not  None:
        vmin = np.float(arguments["--vmin"])
else:
    vmax = np.nanpercentile(cutphi,95)
    vmin = np.nanpercentile(cutphi,5)

cax = ax.imshow(cutphi, cmap, interpolation=None,vmax=vmax,vmin=vmin,alpha=0.6)

divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
if arguments["--title"] ==  None:
    fig.canvas.set_window_title(infile)
else:
    fig.canvas.set_window_title(arguments["--title"])

try:
    del ds
except:
    pass

# Display the data
##ax.set_rasterized(True)
# fig.savefig('{}.eps'.format(infile), format='EPS',dpi=150)
plt.show()
