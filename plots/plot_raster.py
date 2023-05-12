#!/usr/bin/env python3
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
[--format=<value>] [--lectfile=<value>] [--rad2mm=<value>] [--title=<value>] [--wrap=<values>] [--band=<values>] 
       plot\_raster.py --infile=<path> [--cpt=<values>] [<ibeg>] [<iend>] [<jbeg>] [<jend>] \
[--format=<value>] [--parfile=<path>] [--lectfile=<value>] [--rad2mm=<value>] [--title=<value>] [--wrap=<values>] [--vmin=<value>] [--vmax=<value>] [--cols=<values>] [--lines=<values>] [--band=<values>]
       plot\_raster.py --infile=<path> [--cpt=<values>] [--crop=<values>] \
[--format=<value>] [--lectfile=<value>] [--rad2mm=<value>] [--title=<value>] [--wrap=<values>] [--vmin=<value>] [--vmax=<value>] [--band=<values>]  


Options:
-h --help             Show this screen.
--infile=<file>       Raster to be displayed 
--ibeg=<value>        Ligne numbers bounded the cutting zone [default: 0]
--iend=<value>        Ligne numbers bounded the cutting zone [default: nlines]
--jbeg=<value>        Column numbers bounded the cutting zone [default: 0]
--jend=<value>        Column numbers bounded the cutting zone [default: ncols]
--crop=<values>          Crop option with smoothing of boundaries (same as ibeg,iend..) [default: 0,ncols,0,nlines]
--format=<value>      Format input files: ROI_PAC, GAMMA, GTIFF [default: ROI_PAC]
--cpt==<value>        Indicate colorscale for phase
--wrap=<value>        Wrapped phase between value for unwrapped files 
--lectfile=<file>     Path of the lect.in file for r4 format
--parfile=<file>      Path of the .par file of GAMMA
--rad2mm=<value>      Convert data [default: 1]
--tile=<value>        Title plot 
--band=<values>      Select band number [default: 1] 
--vmax                Max colorscale [default: 98th percentile]
--vmin                Min colorscale [default: 2th percentile]
--cols VALUE         Add crosses on pixel column numbers (eg. 200,400,450)
--lines VALUE        Add crosses on pixel lines numbers  (eg. 1200,1200,3000)
"""

import os, sys

# docopt (command line parser)
try:
    from nsbas import docopt
except:
    import docopt

import numpy as np
from numpy.lib.stride_tricks import as_strided

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import setp
from osgeo import gdal
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

# read arguments
arguments = docopt.docopt(__doc__)
infile = arguments["--infile"]

if arguments["--cpt"] is  None:
    # cmap=cm.jet 
    try:
        from matplotlib.colors import LinearSegmentedColormap
        cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
        cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
        cmap = cmap.reversed()
    except:
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
    rad2mm = float(arguments["--rad2mm"])

if arguments["--band"] ==  None:
    band = 1
else:
    band = int(arguments["--band"])

if (arguments["--cols"] is not None and arguments["--lines"] is not None):
    ipix = list(map(int,arguments["--cols"].replace(',',' ').split()))
    jpix = list(map(int,arguments["--lines"].replace(',',' ').split()))
    if len(jpix) != len(ipix):
      raise Exception("ncols and nlines lists are not the same size")

ds_extension = os.path.splitext(infile)[1]
if (ds_extension == ".tif" or ds_extension ==".tiff"):
     sformat = 'GTIFF'

if sformat == "ROI_PAC": 
    ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])
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
        ncols, nlines = list(map(int, open(lecfile).readline().split(None, 2)[0:2]))
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
    ds = gdal.Open(infile, gdal.GA_ReadOnly)
    phase_band = ds.GetRasterBand(band)
    # Attributes
    print("> Driver:   ", ds.GetDriver().ShortName)
    print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
    print("> Datatype: ", gdal.GetDataTypeName(phase_band.DataType))
    nlines, ncols = ds.RasterYSize, ds.RasterXSize

elif sformat == 'GAMMA':
    from parsers import gamma as gm
    if ds_extension == ".diff":
        if arguments["--parfile"] !=  None:
            par_file =  arguments["--parfile"]
            nlines,ncols = gm.readpar(par=par_file)
            phi = gm.readgamma_int(infile,par=par_file)
        else:
            nlines,ncols = gm.readpar()
            phi = gm.readgamma_int(infile)
    else:
        if arguments["--parfile"] !=  None:
            par_file =  arguments["--parfile"]
            nlines,ncols = gm.readpar(par=par_file)
            phi = gm.readgamma(infile,par=par_file)
        else:
            nlines,ncols = gm.readpar()
            phi = gm.readgamma(infile)


# crop
if arguments["--crop"] ==  None:
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
    crop = [0,jend,0,iend]
else:
    crop = list(map(float,arguments["--crop"].replace(',',' ').split()))
jbeg,jend,ibeg,iend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])

# Initialize a matplotlib figure
fig = plt.figure(1,figsize=(12,8))

if sformat == "ROI_PAC":
    if (ds_extension == ".unw" or ds_extension ==".hgt"):
        ax = fig.add_subplot(1,2,1)
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
        hax = ax.imshow(cutamp, cm.Greys_r,vmin=0, vmax=np.percentile(cutamp, 90))
        divider = make_axes_locatable(ax)
        c = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(hax, cax=c)

        ax = fig.add_subplot(1,2,2)
    elif (ds_extension == ".r4" or ds_extension == ""):
        ax = fig.add_subplot(1,1,1)
        cutphi = as_strided(phi[ibeg:iend,jbeg:jend])*rad2mm 

    else:
        ax = fig.add_subplot(1,2,1)
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
            hax = ax.imshow(cutamp, cm.Greys_r, vmin=0, vmax=np.percentile(cutamp, 90))
            divider = make_axes_locatable(ax)
            c = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(hax, cax=c)
            
            ax = fig.add_subplot(1,2,2)
elif sformat == "GTIFF":
    ax = fig.add_subplot(1,1,1)
    phi = phase_band.ReadAsArray(0, 0,
           ds.RasterXSize, ds.RasterYSize,
           ds.RasterXSize, ds.RasterYSize)
    cutphi = as_strided(phi[ibeg:iend,jbeg:jend])*rad2mm

elif sformat == "GAMMA":
    ax = fig.add_subplot(1,1,1)
    if ds_extension == ".diff":
        cutphi = as_strided(phi[ibeg:iend,jbeg:jend])*rad2mm 
        cutamp = np.absolute(cutphi)
        cutphi = np.angle(cutphi)
        hax = ax.imshow(cutamp, cm.Greys_r)
    else:
        cutphi = as_strided(phi[ibeg:iend,jbeg:jend])*rad2mm

if arguments["--wrap"] is not None:
    cutphi =  np.mod(cutphi+float(arguments["--wrap"]),2*float(arguments["--wrap"]))-float(arguments["--wrap"])
    vmax=float(arguments["--wrap"])
    vmin = -vmax
elif (arguments["--vmax"] is not None) or (arguments["--vmin"] is not None):
    if arguments["--vmax"] is not  None:
        vmax = float(arguments["--vmax"])
        if arguments["--vmin"] is not  None:
          vmin = float(arguments["--vmin"])
        else:
          vmin = -vmax
else:
    #print(np.nanmax(cutphi),np.nanmin(cutphi))
    vmax = np.nanpercentile(cutphi,98)
    vmin = np.nanpercentile(cutphi,2)

# replace 0 by nan
try:
    cutphi[cutphi==0] = float('NaN')
except:
    pass
masked_array = np.ma.array(cutphi, mask=np.isnan(cutphi))

cax = ax.imshow(masked_array, cmap, interpolation='nearest',vmax=vmax,vmin=vmin)
#cax = ax.imshow(masked_array, cmap, interpolation='none',vmax=vmax,vmin=vmin)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

# plot optional crosses
if (arguments["--cols"] is not None and arguments["--lines"] is not None):
    for i in range(len(ipix)):
      ax.scatter(ipix[i]-jbeg,jpix[i]-ibeg,marker='x',color='black',s=150.)

#if arguments["--title"] ==  None:
#    fig.canvas.set_window_title(infile)
#else:
#    fig.canvas.set_window_title(arguments["--title"])

try:
    del ds
except:
    pass

# Display the data
fig.tight_layout()
# ax.set_rasterized(True)
try:
    fig.savefig('{}.pdf'.format(infile), format='PDF',dpi=180)
except:
    pass
plt.show()
