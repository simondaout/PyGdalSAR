#!/usr/bin/env python2

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided

import os
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

l = 1
ncol, nlines, N = 2125, 2810,75

cubei = np.fromfile('depl_cumule',dtype=np.float32)
cube = as_strided(cubei[:nlines*ncol*N])
kk = np.flatnonzero(np.logical_or(cube==0, cube>9999))
cube[kk] = float('NaN')
data = cube.reshape((nlines,ncol,N))[700:1300,1180:1780,l]*4.4563
del cubei, cube

cubei = np.fromfile('cube_era5_flat',dtype=np.float32)
cube = as_strided(cubei[:nlines*ncol*N])
kk = np.flatnonzero(np.logical_or(cube==0, cube>9999))
cube[kk] = float('NaN')
model = cube.reshape((nlines,ncol,N))[700:1300,1180:1780,l]*4.4563
del cubei, cube

cubei = np.fromfile('depl_cumule-era5',dtype=np.float32)
cube = as_strided(cubei[:nlines*ncol*N])
kk = np.flatnonzero(np.logical_or(cube==0, cube>9999))
cube[kk] = float('NaN')
flat = cube.reshape((nlines,ncol,N))[700:1300,1180:1780,l]*4.4563
del cubei, cube

from matplotlib.colors import LinearSegmentedColormap
cm_locs = '/home/comethome/jdd/ScientificColourMaps5/by_platform/python/'
cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
cmap = cmap.reversed()

vmax = np.nanpercentile(flat, 98)
vmin = np.nanpercentile(flat, 2)

fig = plt.figure(1,figsize=(10,6))

ax = fig.add_subplot(1,3,1)
cax = ax.imshow(data,cmap=cmap,vmax=vmax,vmin=vmin)
plt.setp( ax.get_xticklabels(), visible=False)
plt.setp( ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(1,3,2)
cax = ax.imshow(model,cmap=cmap,vmax=vmax,vmin=vmin)
plt.setp( ax.get_xticklabels(), visible=False)
plt.setp( ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(1,3,3)
cax = ax.imshow(flat,cmap=cmap,vmax=vmax,vmin=vmin)
plt.setp( ax.get_xticklabels(), visible=False)
plt.setp( ax.get_yticklabels(), visible=False)

fig.savefig('era5-correction.eps', format='EPS',dpi=150)
plt.show()