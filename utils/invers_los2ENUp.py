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

# output = 'lake'
# # ASC T099
# infile1 = wdir + "T099/T099_LOSVel_mmyr_nan_s90_cropL.tiff"
# lookf1 = wdir +"T099/T099_look_s90_cropL.tiff"
# headf1 = wdir +"T099/T099_head_s90_cropL.tiff"
# sigf1 = wdir +"T099/T099_LOS_sigmaVel_mmyr_nan_s90_cropL.tiff"
# # Descending T004
# infile2 = wdir +"T004/T004_LOSVel_mmyr_nan_s90_cropL.tiff"
# lookf2 = wdir + "T004/T004_look_s90_cropL.tiff"
# headf2 = wdir +"T004/T004_head_s90_cropL.tiff"
# sigf2 = wdir + "T004/T004_LOS_sigmaVel_mmyr_nan_s90_cropL.tiff"


output = 'glacier'
# ASC T172
infile1 = wdir + "T172/T172_LOSVel_mmyr_nan.grd_s90_cropV.tiff"
lookf1 = wdir +"T172/T172_look_s90_cropV.tiff"
headf1 = wdir +"T172/T172_head_s90_cropV.tiff"
sigf1 = wdir +"T172/T172_LOS_sigVel_mmyr_nan.grd_s90_cropV.tiff"
# Descending T004
infile2 = wdir +"T004/T004_LOSVel_mmyr_nan_s90_cropV.tiff"
lookf2 = wdir + "T004/T004_look_s90_cropV.tiff"
headf2 = wdir +"T004/T004_head_s90_cropV.tiff"
sigf2 = wdir + "T004/T004_LOS_sigmaVel_mmyr_nan_s90_cropV.tiff"
# Descending T099
infile3 = wdir +"T099/T099_LOSVel_mmyr_nan_s90_cropV.tiff"
lookf3 = wdir + "T099/T099_look_s90_cropV.tiff"
headf3 = wdir +"T099/T099_head_s90_cropV.tiff"
sigf3 = wdir + "T099/T099_LOS_sigmaVel_mmyr_nan_s90_cropV.tiff"


################################

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

ds = gdal.Open(sigf1,gdal.GA_ReadOnly)
band = ds.GetRasterBand(1)
sig1 = -band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
sig1[np.isnan(los1)] = np.float('NaN') 
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

################################

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

ds = gdal.Open(sigf2,gdal.GA_ReadOnly)
band = ds.GetRasterBand(1)
sig2 = -band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
sig2[np.isnan(los2)] = np.float('NaN') 
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


################################

print('Read ascending track: {}'.format(infile3))

ds = gdal.Open(infile3,gdal.GA_ReadOnly)
band = ds.GetRasterBand(1)
los3 = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
band.FlushCache()
nlines,ncols = ds.RasterYSize, ds.RasterXSize
print('Number of pixel: {}'.format(len(los3.flatten())))
print('Nlines: {}, Ncols: {}'.format(nlines,ncols))
del ds, band


ds = gdal.Open(lookf3,gdal.GA_ReadOnly)
gt = ds.GetGeoTransform()
proj = ds.GetProjectionRef()
driver = gdal.GetDriverByName('GTiff')
band = ds.GetRasterBand(1)
look3 = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
band.FlushCache()
del ds, band

ds = gdal.Open(sigf3,gdal.GA_ReadOnly)
band = ds.GetRasterBand(1)
sig3 = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
sig3[np.isnan(los3)] = np.float('NaN') 
band.FlushCache()
del ds, band

ds = gdal.Open(headf3,gdal.GA_ReadOnly)
band = ds.GetRasterBand(1)
head3 = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
del ds, band

# convert head, look to angle phi, theta in rad
# theta: vertical angle between LOS and vertical
theta3 = np.deg2rad(90.-look3)
# phi: horizontal angle between LOS and EAST
phi3 = np.deg2rad(-90-head3)

# compute proj
proj3=[np.cos(theta3)*np.cos(phi3),
                np.cos(theta3)*np.sin(phi3),
                np.sin(theta3)
                ]

print('Average LOS projection to east, north, up: {0:.5f} {1:.5f} {2:.5f}'.\
    format(np.nanmean(proj3[0]),np.nanmean(proj3[1]),np.nanmean(proj3[2])))
print()

################################
# plot DATA

cm_locs = '/home/comethome/jdd/ScientificColourMaps5/by_platform/python/'
cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
cmap_r = cmap.reversed()
fig=plt.figure(0, figsize=(14,12))
vmax=3
vmin=-9

# plot los1
ax = fig.add_subplot(4,3,1)
cax = ax.imshow(los1,cmap=cmap,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('LOS Asc.')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

# plot los2
ax = fig.add_subplot(4,3,2)
cax = ax.imshow(los2,cmap=cmap,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('LOS Desc.')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

# plot los2
ax = fig.add_subplot(4,3,3)
cax = ax.imshow(los3,cmap=cmap,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('LOS Desc.')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

vmax = 1
vmin = 0

# plot SIGMA los1
ax = fig.add_subplot(4,3,4)
cax = ax.imshow(sig1,cmap=cmap_r,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('Sigma LOS Asc.')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

# plot SIGMA los2
ax = fig.add_subplot(4,3,5)
cax = ax.imshow(sig2,cmap=cmap_r,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('sigma LOS Desc.')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

ax = fig.add_subplot(4,3,6)
cax = ax.imshow(sig3,cmap=cmap_r,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('sigma LOS Desc.')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

vmax=10
vmin=-10

# plt.show()

################################

# prepare output matrix
east, north, up = np.ones((nlines,ncols))*np.float('NaN'), np.ones((nlines,ncols))*np.float('NaN'), \
np.ones((nlines,ncols))*np.float('NaN')
seast, snorth, sup = np.ones((nlines,ncols))*np.float('NaN'), np.ones((nlines,ncols))*np.float('NaN'), \
np.ones((nlines,ncols))*np.float('NaN')


## least-square 2 views 
# # loop over each line and cols
# for i in range(nlines):
#     for j in range(ncols):

#         if not np.isnan(los1[i,j]) and not np.isnan(los2[i,j]): 
#             # build data matrix 
#             data = np.zeros((2))
#             data[0] = los1[i,j]
#             data[1] = los2[i,j]

#             # Build G matrix
#             G = np.zeros((2,2))
#             G[0,0] = proj1[0][i,j]
#             G[0,1] = proj1[2][i,j]
#             G[1,0] = proj2[0][i,j]
#             G[1,1] = proj2[2][i,j]

#             # build uncertainty matrix
#             rms = np.zeros((2))
#             rms[0] = sig1[i,j]
#             rms[1] = sig1[i,j]

#             # Inversion
#             x0 = lst.lstsq(G,data)[0]
#             _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
#             _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
            
#             pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1.e-9)[0]
#             east[i,j] = pars[0]; up[i,j] = pars[1]
#             # print('east, north, up: {0:.5f} {1:.5f} {2:.5f}'.\
#             # format(ve,vn,vup)) 

#             Cd = np.diag(rms**2, k = 0)
#             sigmam = np.linalg.inv(np.dot(np.dot(G.T,np.linalg.inv(Cd)),G))
#             seast[i,j], sup[i,j] = sigmam[0,0], sigmam[1,1]
#             # A = (np.dot(G,G.T) + Cd)
#             # A1 = np.linalg.inv(A)
#             # sigmam = 1 - np.dot(G.T, A1)

#             # print('sigma east: {0:.5f}, sigma up: {1:.5f}'.\
#             # format(sigmam[0,0],sigmam[1,1])) 
#             # sys.exit()

## least-square 3 views 
# loop over each line and cols
for i in range(nlines):
    for j in range(ncols):

        if (not np.isnan(los1[i,j]) or not np.isnan(los3[i,j])) and not np.isnan(los2[i,j]):  
            # build data matrix 
            data = np.zeros((3))
            data[0] = los1[i,j]
            data[1] = los2[i,j]
            data[2] = los3[i,j]

            # Build G matrix
            G = np.zeros((3,2))
            G[0,0] = proj1[0][i,j]
            G[0,1] = proj1[2][i,j]
            G[1,0] = proj2[0][i,j]
            G[1,1] = proj2[2][i,j]
            G[2,0] = proj3[0][i,j]
            G[2,1] = proj3[2][i,j]

            # build uncertainty matrix
            rms = np.zeros((3))
            rms[0] = sig1[i,j]
            rms[1] = sig2[i,j]
            rms[2] = sig3[i,j]

            if np.isnan(los3[i,j]) or np.isnan(sig3[i,j]):
                data = np.delete(data,2,0); G = np.delete(G,2,0); rms = np.delete(rms,2,0)
            if np.isnan(los1[i,j]) or np.isnan(sig1[i,j]):
                data =np.delete(data,0,0); G = np.delete(G,0,0); rms = np.delete(rms,0,0)

            # Inversion
            x0 = lst.lstsq(G,data)[0]
            _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
            _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
            
            pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1.e-9)[0]
            east[i,j] = pars[0]; up[i,j] = pars[1]
            # print('east, north, up: {0:.5f} {1:.5f} {2:.5f}'.\
            # format(ve,vn,vup)) 

            Cd = np.diag(rms**2, k = 0)
            sigmam = np.linalg.inv(np.dot(np.dot(G.T,np.linalg.inv(Cd)),G))
            seast[i,j], sup[i,j] = sigmam[0,0], sigmam[1,1]
            # A = (np.dot(G,G.T) + Cd)
            # A1 = np.linalg.inv(A)
            # sigmam = 1 - np.dot(G.T, A1)

            # print('sigma east: {0:.5f}, sigma up: {1:.5f}'.\
            # format(sigmam[0,0],sigmam[1,1])) 
            # sys.exit()


################################

cm_locs = '/home/comethome/jdd/ScientificColourMaps5/by_platform/python/'
cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
cmap_r = cmap.reversed()
fig=plt.figure(1, figsize=(14,12))
vmax=3
vmin=-9

# plot los1
ax = fig.add_subplot(4,3,1)
cax = ax.imshow(los1,cmap=cmap,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('LOS Asc.')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

# plot los2
ax = fig.add_subplot(4,3,2)
cax = ax.imshow(los2,cmap=cmap,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('LOS Desc.')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

# plot los2
ax = fig.add_subplot(4,3,3)
cax = ax.imshow(los3,cmap=cmap,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('LOS Desc.')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

vmax = 1
vmin = 0

# plot SIGMA los1
ax = fig.add_subplot(4,3,4)
cax = ax.imshow(sig1,cmap=cmap_r,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('Sigma LOS Asc.')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

# plot SIGMA los2
ax = fig.add_subplot(4,3,5)
cax = ax.imshow(sig2,cmap=cmap_r,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('sigma LOS Desc.')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

ax = fig.add_subplot(4,3,6)
cax = ax.imshow(sig3,cmap=cmap_r,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('sigma LOS Desc.')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

vmax=10
vmin=-10

# plot east
ax = fig.add_subplot(4,3,7)
cax = ax.imshow(east,cmap=cmap,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('Ue')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

# plot north
ax = fig.add_subplot(4,3,8)
cax = ax.imshow(north,cmap=cmap,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('Un')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

# plot up
ax = fig.add_subplot(4,3,9)
cax = ax.imshow(up,cmap=cmap,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('Uz')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
fig.tight_layout()
fig.savefig('decomposition.png', format='PNG',dpi=150)

vmax = 1
vmin = 0

# plot sigma east
ax = fig.add_subplot(4,3,10)
cax = ax.imshow(seast,cmap=cmap_r,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('Ue')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

# plot sigma north
ax = fig.add_subplot(4,3,11)
cax = ax.imshow(snorth,cmap=cmap_r,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('Un')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

# plot sigma up
ax = fig.add_subplot(4,3,12)
cax = ax.imshow(sup,cmap=cmap_r,vmax=vmax,vmin=vmin,interpolation=None)
ax.set_title('Uz')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)


fig.tight_layout()
fig.savefig('decomposition_{}.pdf'.format(output), format='PDF',dpi=150)

# Save output files
ds = driver.Create('Ue_{}.tif'.format(output), ncols, nlines, 1, gdal.GDT_Float32)
band = ds.GetRasterBand(1)
band.WriteArray(east)
ds.SetGeoTransform(gt)
ds.SetProjection(proj)
band.FlushCache()

ds = driver.Create('Un_{}.tif'.format(output), ncols, nlines, 1, gdal.GDT_Float32)
band = ds.GetRasterBand(1)
band.WriteArray(north)
ds.SetGeoTransform(gt)
ds.SetProjection(proj)
band.FlushCache()

ds = driver.Create('Uz_{}.tif'.format(output), ncols, nlines, 1, gdal.GDT_Float32)
band = ds.GetRasterBand(1)
band.WriteArray(up)
ds.SetGeoTransform(gt)
ds.SetProjection(proj)
band.FlushCache()

plt.show()
