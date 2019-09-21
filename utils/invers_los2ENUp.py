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

import sys,os
import numpy as np
import scipy.optimize as opt
import scipy.linalg as lst
import gdal
gdal.UseExceptions()

# plot
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap

wdir = '/home/cometraid14/daouts/work/tibet/qinghai/data/insar/'

# ASC T172
infile1 = wdir + "T172/T172_LOSVel_mmyr_nan.grd_s90_cropV.tiff"
lookf1 = wdir +"T172/T172_look_s90_cropV.tiff"
headf1 = wdir +"T172/T172_head_s90_cropV.tiff"

print('Read ascending track: {}'.format(infile1))

ds = gdal.Open(infile1,gdal.GA_ReadOnly)
band = ds.GetRasterBand(1)
los1 = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
band.FlushCache()
nlines,ncols = ds.RasterYSize, ds.RasterXSize
print('Number of pixel: {}'.format(len(los1.flatten())))
print('Nlines: {}, Ncols: {}'.format(nlines,ncols))
del ds, band


ds = gdal.Open(lookf1,gdal.GA_ReadOnly)
gt = ds.GetGeoTransform()
proj = ds.GetProjectionRef()
driver = gdal.GetDriverByName('GTiff')
band = ds.GetRasterBand(1)
look1 = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
band.FlushCache()
del ds, band

ds = gdal.Open(headf1,gdal.GA_ReadOnly)
band = ds.GetRasterBand(1)
head1 = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
del ds, band

# convert head, look to angle phi, theta in rad
# theta: vertical angle between LOS and vertical
theta1 = np.deg2rad(90.-look1)
# phi: horizontal angle between LOS and EAST
phi1 = np.deg2rad(-90-head1)

# compute proj
proj1=[np.cos(theta1)*np.cos(phi1),
                np.cos(theta1)*np.sin(phi1),
                np.sin(theta1)
                ]

print('Average LOS projection to east, north, up: {0:.5f} {1:.5f} {2:.5f}'.\
    format(np.nanmean(proj1[0]),np.nanmean(proj1[1]),np.nanmean(proj1[2])))
print()

# Descending T004
infile2 = wdir +"T004/T004_LOSVel_mmyr_nan_s90_cropV.tiff"
lookf2 = wdir + "T004/T004_look_s90_cropV.tiff"
headf2 = wdir +"T004/T004_head_s90_cropV.tiff"

print('Read descending track: {}'.format(infile2))

ds = gdal.Open(infile2,gdal.GA_ReadOnly)
band = ds.GetRasterBand(1)
los2 = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
band.FlushCache()
del ds, band

ds = gdal.Open(lookf2,gdal.GA_ReadOnly)
band = ds.GetRasterBand(1)
look2 = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
band.FlushCache()
del ds, band

ds = gdal.Open(headf2,gdal.GA_ReadOnly)
band = ds.GetRasterBand(1)
head2 = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
band.FlushCache()
del ds, band

# convert head, look to angle phi, theta in rad
theta2 = np.deg2rad(90.-look2)
phi2 = np.deg2rad(-90-head2)

# compute proj
proj2=[np.cos(theta2)*np.cos(phi2),
                np.cos(theta2)*np.sin(phi2),
                np.sin(theta2)
                ]

print('Average LOS projection to east, north, up: {0:.5f} {1:.5f} {2:.5f}'.\
    format(np.nanmean(proj2[0]),np.nanmean(proj2[1]),np.nanmean(proj2[2])))
print()

# prepare output matrix
east, north, up = np.ones((nlines,ncols))*np.float('NaN'), np.ones((nlines,ncols))*np.float('NaN'), \
np.ones((nlines,ncols))*np.float('NaN')
seast, snorth, sup = np.ones((nlines,ncols))*np.float('NaN'), np.ones((nlines,ncols))*np.float('NaN'), \
np.ones((nlines,ncols))*np.float('NaN')

# loop over each line and cols
for i in range(nlines):
    for j in range(ncols):

        if not np.isnan(los1[i,j]) and not np.isnan(los2[i,j]): 
            # build data matrix 
            data = np.zeros((2))
            data[0] = los1[i,j]
            data[1] = los2[i,j]

            # Build G matrix
            G = np.zeros((2,3))
            G[0,0] = proj1[0][i,j]
            G[0,1] = proj1[1][i,j]
            G[0,2] = proj1[2][i,j]
            G[1,0] = proj2[0][i,j]
            G[1,1] = proj2[1][i,j]
            G[1,2] = proj2[2][i,j]

            # build uncertainty matrix
            rms = np.ones((2))

            # Inversion
            x0 = lst.lstsq(G,data)[0]
            _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
            _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
            pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1.e-9)[0]
            east[i,j] = pars[0]; north[i,j] = pars[1]; up[i,j] = pars[2];
            # print('east, north, up: {0:.5f} {1:.5f} {2:.5f}'.\
            # format(ve,vn,vup)) 

            # try:
            # print(np.dot(G.T,G))
            # varx = np.linalg.inv(np.dot(G.T,G))
            # res = np.sum(pow((d-np.dot(G,pars)),2))
            # scale = 1./(G.shape[0]-G.shape[1])
            # sigmam = np.sqrt(scale*res*np.diag(varx)) 
            # except:
            #     sigmam = np.ones((G.shape[1]))*float('NaN')

            # print('sigma east, north, up: {0:.5f} {1:.5f} {2:.5f}'.\
            # format(sigmam[0],sigmam[1],sigmam[2])) 
            # sys.exit()

cm_locs = '/home/comethome/jdd/ScientificColourMaps5/by_platform/python/'
cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
fig=plt.figure(figsize=(14,12))
vmax=3
vmin=-9

# plot los1
ax = fig.add_subplot(2,3,1)
cax = ax.imshow(los1,cmap=cmap,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('LOS Asc.')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

# plot los2
ax = fig.add_subplot(2,3,2)
cax = ax.imshow(los2,cmap=cmap,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('LOS Desc.')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

vmax=10
vmin=-10

# plot east
ax = fig.add_subplot(2,3,4)
cax = ax.imshow(east,cmap=cmap,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('Ue')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

# plot north
ax = fig.add_subplot(2,3,5)
cax = ax.imshow(north,cmap=cmap,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('Un')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

# plot up
ax = fig.add_subplot(2,3,6)
cax = ax.imshow(up,cmap=cmap,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('Uz')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
fig.tight_layout()
fig.savefig('decomposition.png', format='PNG',dpi=150)

plt.show()

# Save output files
ds = driver.Create('Ue.tif', ncols, nlines, 1, gdal.GDT_Float32)
band = ds.GetRasterBand(1)
band.WriteArray(east)
ds.SetGeoTransform(gt)
ds.SetProjection(proj)
band.FlushCache()

ds = driver.Create('Un.tif', ncols, nlines, 1, gdal.GDT_Float32)
band = ds.GetRasterBand(1)
band.WriteArray(north)
ds.SetGeoTransform(gt)
ds.SetProjection(proj)
band.FlushCache()

ds = driver.Create('Uz.tif', ncols, nlines, 1, gdal.GDT_Float32)
band = ds.GetRasterBand(1)
band.WriteArray(up)
ds.SetGeoTransform(gt)
ds.SetProjection(proj)
band.FlushCache()