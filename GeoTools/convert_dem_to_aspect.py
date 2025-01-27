#!/usr/bin/env python3
# -*- coding: utf-8 -*-

print()
print()
print('Author: Simon DAOUT')
print()
print()

import numpy as np
import scipy.ndimage
import scipy.optimize as opt
import scipy.linalg as lst
from osgeo import gdal, osr
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
from os import path, sys

def compute_slope_aspect(path):
   
    #### LOAD DEM
    ds = gdal.Open(path,gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    topo = band.ReadAsArray()
    ncols, nlines = ds.RasterYSize, ds.RasterXSize
    gt = ds.GetGeoTransform()
    projref = ds.GetProjectionRef()
    drv = gdal.GetDriverByName('GTiff')

    # filter over a two pixels windows to remove outliers
    filtrer = scipy.ndimage.gaussian_filter(topo,sigma = (2., 2.))
    
    # Get middle latitude
    data = ds.GetGeoTransform()   
    lats = data[3] + (np.arange(nlines) * data[5])
    lat_moy = np.mean(lats)
    print("Metadata:", data)
    print("Average lattitude:", lat_moy)    

    # Get resolution depending on the coordinate reference system 
    res_x = data[1]
    res_y = data[5]
    srs = osr.SpatialReference(wkt=ds.GetProjection()) 
    #if abs(res_x) < 0.01: # si inf. à 0.01m, certainement en degré, mais serait plus propre de vérifier le SCR
    if srs.IsGeographic():
        res_x = (res_x*40075e3/360) # conversion degre to meter 
        res_y = (res_y*40075e3/360)*np.cos(np.deg2rad(lat_moy)) 
        print("Geographic coordinates")
        print("dx: {}, dy: {}".format(res_x, res_y))
        dsm_dy, dsm_dx = np.gradient(filtrer, res_x, res_y) 
    else:
        print("Cartographic coordinates")
        print("dx: {}, dy: {}".format(res_x, res_y))
        dsm_dy, dsm_dx = np.gradient(filtrer, res_x, res_y) # d'abord resolution en y (négative) puis x
        # dsm_dy est la pente dans la direction des lignes
        # dsm_dx est la pente dans la direction des colonnes
 
    # Compute slope and downslope direction
    slope = np.sqrt(dsm_dx**2 + dsm_dy**2)
    aspect = np.arctan2(-dsm_dy, dsm_dx)    # clockwise slope direction from 0 to 360
  
    #crop_slice = (slice(450, 850), slice(500, 900))
    #aspect_crop = np.rad2deg(aspect[crop_slice]) 
    #fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    #cax= ax.imshow(aspect_crop, cmap='Greys_r', origin='upper', alpha=0.5, vmax=180, vmin=0)
    #ax.set_xticks([])
    #ax.set_yticks([])
    #divider = make_axes_locatable(ax)
    #c = divider.append_axes("right", size="5%", pad=0.05)
    #plt.colorbar(cax, cax=c) 
    #plt.show()
    #sys.exit()

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

    # Plot DEM, Slope, gradients and aspect
    fig = plt.figure(1,figsize=(11,7))
    cmap = cm.terrain
    # Plot topo
    ax = fig.add_subplot(2,2,1)
    cax = ax.imshow(filtrer, cmap='Greys_r', vmax=np.nanpercentile(filtrer,98), vmin=np.nanpercentile(filtrer,2))
    ax.set_title('DEM',fontsize=6)
    fig.colorbar(cax, orientation='vertical')

    ax = fig.add_subplot(2,2,2)
    cax = ax.imshow(np.rad2deg(slope), cmap='Greys_r', vmax=np.nanpercentile(np.rad2deg(slope),90), vmin=np.nanpercentile(np.rad2deg(slope),10))
    ax.set_title('Slope',fontsize=6)
    fig.colorbar(cax, orientation='vertical')

    ax = fig.add_subplot(2,2,3)
    cax = ax.imshow(dsm_dy,cmap='Greys_r', vmax=np.nanpercentile(dsm_dy,98), vmin=np.nanpercentile(dsm_dy,2))
    ax.set_title('Gradient in y',fontsize=6)
    fig.colorbar(cax, orientation='vertical')

    ax = fig.add_subplot(2,2,4)
    cax = ax.imshow(np.rad2deg(aspect), cmap='Greys_r', vmax=180, vmin=-180)
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
