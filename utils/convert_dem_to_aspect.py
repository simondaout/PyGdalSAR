#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy.ndimage
import scipy.optimize as opt
import scipy.linalg as lst
from osgeo import gdal
gdal.UseExceptions()
from math import sin, cos
from numpy.lib.stride_tricks import as_strided
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
from sys import argv,exit,stdin,stdout
import getopt
import os, math
from os import path

def compute_slope_aspect(path):
   
    #### LOAD DEM
    ds = gdal.Open(path,gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    topo = band.ReadAsArray()
    ncols, nlines = ds.RasterYSize, ds.RasterXSize
    filtrer = scipy.ndimage.gaussian_filter(topo,2.)
    gt = ds.GetGeoTransform()
    projref = ds.GetProjectionRef()
    drv = gdal.GetDriverByName('GTiff')

    # Get middle latitude
    data1 = ds.GetGeoTransform()   
    lats = data1[3] + (np.arange(nlines) * data1[5])
    lat_moy = np.mean(lats)
    
    # Get resolution depending on whether the 
    res = data1[1]*40075e3/360
    print("Spatial resolution in deg: {}, in meter: {}".format(data1[1],res))
    if res<1 or res>500:
        print("Spatial resolution seems unrealistic. Exit!")
        exit()  
    
    # Calcul gradient
    Py, Px = np.gradient(filtrer, res, res*np.cos(np.deg2rad(lat_moy)))
    Px = Px.astype(float); Py = Py.astype(float)
    # Slope
    slope = np.arctan(np.sqrt(Py**2+Px**2))
    # smooth slope to avoid crazy values 
    # Line of max slope
    aspect = np.arctan2(Py,-Px)
    
    # Create aspect and slope files
    dst = drv.Create('slope.tif', nlines, ncols, 1, gdal.GDT_Float32)
    bandt = dst.GetRasterBand(1)
    bandt.WriteArray(np.rad2deg(slope))
    dst.SetGeoTransform(gt)
    dst.SetProjection(projref)
    bandt.FlushCache()

    dst = drv.Create('aspect.tif', nlines, ncols, 1, gdal.GDT_Float32)
    bandt = dst.GetRasterBand(1)
    bandt.WriteArray(np.rad2deg(aspect))
    dst.SetGeoTransform(gt)
    dst.SetProjection(projref)
    bandt.FlushCache()

    # Plot DEM, Slope, Py and aspect
    fig = plt.figure(1,figsize=(11,7))
    cmap = cm.terrain
    # Plot topo
    ax = fig.add_subplot(2,2,1)
    cax = ax.imshow(filtrer,cmap=cmap,vmax=np.nanpercentile(filtrer,98),vmin=np.nanpercentile(filtrer,2))
    ax.set_title('DEM',fontsize=6)
    fig.colorbar(cax, orientation='vertical')

    ax = fig.add_subplot(2,2,2)
    cax = ax.imshow(np.rad2deg(slope),cmap=cmap,vmax=np.nanpercentile(np.rad2deg(slope),98),vmin=np.nanpercentile(np.rad2deg(slope),2))
    ax.set_title('Slope',fontsize=6)
    fig.colorbar(cax, orientation='vertical')

    ax = fig.add_subplot(2,2,3)
    cax = ax.imshow(Py,cmap=cmap,vmax=np.nanpercentile(Py,98),vmin=np.nanpercentile(Py,2))
    ax.set_title('Py',fontsize=6)
    fig.colorbar(cax, orientation='vertical')

    ax = fig.add_subplot(2,2,4)
    cax = ax.imshow(np.rad2deg(aspect),cmap=cmap,vmax=np.nanpercentile(aspect,98),vmin=np.nanpercentile(aspect,2))
    ax.set_title('aspect',fontsize=6)
    fig.colorbar(cax, orientation='vertical')
    
    fig.savefig('dem_slope_aspect.png',format='PNG',dpi=300)

    if PLOT:
        plt.show()      

    return -aspect, slope

def usage():
  print('convert_dem_to_aspect.py demfile [-v] [-h]')
  print('-v Verbose mode. Show more information about the processing')
  print('-h Show this screen')

try:
    opts,args = getopt.getopt(argv[1:], "h", ["help"])
except:
    print("for help use --help")
    exit()

for o in argv:
    if o in ("-h","--help"):
       usage()
       exit()
    if o in ("-v","--verbose"):
      level = 'debug'

PLOT=True
demfile = argv[1]
compute_slope_aspect(demfile)
