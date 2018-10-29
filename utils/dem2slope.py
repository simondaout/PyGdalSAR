#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-


"""\
dem2slope.py
-------------

Usage: dem2slope.py --infile=<path> --outfile=<path> --los=<value> --heading=<value> --lat=<value>
dem2slope.py -h | --help

Options:
-h --help           Show this screen
--infile PATH       DEM file 
--outfile PATH      Output file 
--los VALUE         Mean Los angle
--heading VALUE.    Mean Heading angle
--lat VALUE.    	Average latitude
"""

import gdal
import scipy.ndimage
from nsbas import docopt
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm

# read arguments
arguments = docopt.docopt(__doc__)
infile = arguments["--infile"]
outfile = arguments["--outfile"]
look = float(arguments["--los"])*np.pi/180
heading = float(arguments["--heading"])*np.pi/180
lat = float(arguments["--lat"])

# read input
ds = gdal.Open(infile, gdal.GA_ReadOnly)
ds_band = ds.GetRasterBand(1)
topo = ds_band.ReadAsArray()

# create output file
drv = gdal.GetDriverByName('GTiff')
ds2 = drv.CreateCopy(outfile,ds)
ds2_band = ds2.GetRasterBand(1)

topo = ds_band.ReadAsArray()
toposmooth = scipy.ndimage.filters.gaussian_filter(topo,.1)
Py, Px = np.gradient(toposmooth,90,90*np.cos(lat*np.pi/180))
slope = np.sqrt(Px**2+Py**2)
slopelos = (np.cos(heading)*Px+np.sin(heading)*Py)/np.sin(look)
ds2_band.WriteArray(slopelos)
del ds2

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




