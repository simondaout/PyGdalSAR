#!/usr/bin/env python2
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
dem2slope.py
-------------
Convert dem to slope file in the LOS geeometry knowing los heading and lattitude mean

Usage: 
    dem2slope.py --infile=<path> --los=<value> --heading=<value>  [--out_slope=<path>] [--out_slopelos=<path>] [--lat=<value>] [--red=<value>]
    dem2slope.py --infile=<path> --incidence=<path> [--out_slope=<path>] [--out_slopelos=<path>] [--lat=<value>] [--az_res=<value>] [--rg_res=<value>]

dem2slope.py -h | --help

Options:
-h --help           Show this screen
--infile PATH       DEM file 
--out_slope FILE    Output slope file ['Slope.r4']
--out_slopelos FILE Output slope in the LOS file ['LOSslope.r4']
--los VALUE         Mean Los angle
--heading VALUE.    Mean Heading angle
--lat VALUE.        Average latitude in deg.. if raster is given in WGS84 [default: None]
--az_res VALUE.     Resolutiom azimuth pixel size (size pixel) [default: 30m]
--rg_res VALUE.     Resolutiom range pixel size (size pixel) [default: 30m]
--incidence=<file>  Path to incidence file .unw 
"""

import gdal, os
import scipy.ndimage
import docopt
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
np.warnings.filterwarnings('ignore')

# read arguments
arguments = docopt.docopt(__doc__)
infile = arguments["--infile"]
if arguments["--out_slope"] == None:
    arguments["--out_slope"] = 'Slope.r4'
if arguments["--out_slopelos"] == None:
    arguments["--out_slopelos"] = 'LOSslope.r4'

if arguments["--lat"] == None:
    lat = 0
else:
    lat = np.deg2rad(float(arguments["--lat"]))
if arguments["--az_res"] == None:
    az_res = 30
else:
    az_res = float(arguments["--az_res"])
if arguments["--rg_res"] == None:
    rg_res = 30
else:
    rg_res = float(arguments["--rg_res"])

# read input
ds_extension = os.path.splitext(infile)[1]
if (ds_extension == ".r4" or ds_extension == ""):
    fid = open(infile, 'r')
    ncols, nlines = map(int, open('lect.in').readline().split(None, 2)[0:2])
    topo = np.fromfile(fid,dtype=np.float32)[:nlines*ncols].reshape((nlines,ncols))
    fid.close()
else:
    ds = gdal.Open(infile, gdal.GA_ReadOnly)
    ds_band = ds.GetRasterBand(1)
    topo = ds_band.ReadAsArray()
    ncols, nlines = ds.RasterYSize, ds.RasterXSize

look, heading = np.zeros((nlines,ncols)),np.zeros((nlines,ncols))
if arguments["--incidence"] == None:
    look.fill(np.deg2rad(float(arguments["--los"])))
    heading.fill(np.deg2rad(float(arguments["--heading"])))
else:
    ds = gdal.Open(arguments["--incidence"], gdal.GA_ReadOnly)
    ds_band1 = ds.GetRasterBand(1)
    ds_band2 = ds.GetRasterBand(2)
    look[:,:] = np.deg2rad(ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize))[:nlines,:ncols]
    # in roipac incidence file, band 2 is the angle between north and horizontal LOS
    heading[:,:] = np.deg2rad(90 - ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize))[:nlines,:ncols]

toposmooth = scipy.ndimage.filters.gaussian_filter(topo,2.)
Py, Px = np.gradient(toposmooth,az_res,rg_res*np.cos(lat))
slope = np.sqrt(Px**2+Py**2)
slopelos = (np.cos(heading)*Px+np.sin(heading)*Py)/np.sin(look)

if (ds_extension == ".r4" or ds_extension == ""):
    fid1 = open(arguments["--out_slopelos"],'wb')
    slopelos.flatten().astype('float32').tofile(fid1)
    fid1.close()
    fid2 = open(arguments["--out_slope"],'wb')
    slope.flatten().astype('float32').tofile(fid2)
    fid2.close()
else:
    # create output file
    drv = gdal.GetDriverByName('GTiff')
    ds2 = drv.CreateCopy(arguments["--out_slopelos"],ds)
    ds2_band = ds2.GetRasterBand(1)
    ds2_band.WriteArray(slopelos)
    del ds2. drv
    drv = gdal.GetDriverByName('GTiff')
    ds2 = drv.CreateCopy(arguments["--out_slope"],ds)
    ds2_band = ds2.GetRasterBand(1)
    ds2_band.WriteArray(slope)
    del ds2, drv

# ibeg,iend=600,900
# jbeg,jend=900,1200

# initiate figure depl
fig = plt.figure(1,figsize=(9,8))
# plot topo
ax = fig.add_subplot(2,2,1)
cmap = cm.gist_earth
cmap.set_bad('white')
# cax = ax.imshow(toposmooth[ibeg:iend,jbeg:jend],cmap=cmap,vmax=np.nanpercentile(toposmooth,98),vmin=np.nanpercentile(toposmooth,2))
cax = ax.imshow(toposmooth,cmap=cmap,vmax=np.nanpercentile(toposmooth,98),vmin=np.nanpercentile(toposmooth,2))
ax.set_title('DEM',fontsize=6)
fig.colorbar(cax, orientation='vertical',aspect=10)

ax = fig.add_subplot(2,2,2)
cmap = cm.jet
cmap.set_bad('white')
# cax = ax.imshow(Px[ibeg:iend,jbeg:jend],cmap=cmap,vmax=np.nanpercentile(Px,98),vmin=np.nanpercentile(Px,2))
cax = ax.imshow(Px,cmap=cmap,vmax=np.nanpercentile(Px,98),vmin=np.nanpercentile(Px,2))
ax.set_title('Gradient in West',fontsize=6)
fig.colorbar(cax, orientation='vertical',aspect=10)

ax = fig.add_subplot(2,2,3)
cmap = cm.jet
cmap.set_bad('white')
# cax = ax.imshow(Py[ibeg:iend,jbeg:jend],cmap=cmap,vmax=np.nanpercentile(Py,98),vmin=np.nanpercentile(Py,2))
cax = ax.imshow(Py,cmap=cmap,vmax=np.nanpercentile(Py,98),vmin=np.nanpercentile(Py,2))
ax.set_title('Gradient in North',fontsize=6)
fig.colorbar(cax, orientation='vertical',aspect=10)

ax = fig.add_subplot(2,2,4)
cmap = cm.jet
cmap.set_bad('white')
# cax = ax.imshow(slopelos[ibeg:iend,jbeg:jend],cmap=cmap,vmax=np.nanpercentile(slopelos,98),vmin=np.nanpercentile(slopelos,2))
cax = ax.imshow(slopelos,cmap=cmap,vmax=np.nanpercentile(slopelos,98),vmin=np.nanpercentile(slopelos,2))
ax.set_title('Slope Toward the Satellite',fontsize=6)
fig.colorbar(cax, orientation='vertical',aspect=10)

plt.show()




