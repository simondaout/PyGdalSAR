#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from osgeo import gdal
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pylab import *
import docopt
import os

from matplotlib.colors import LinearSegmentedColormap
cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
cmap = cmap.reversed()

file1 = 'stack_amp_da_2rlks.unw'
file2 = '20161014/grad_radar_2rlks.hgt'
outfile='cutfile_2rlks.unw'

ds = gdal.OpenEx(file1, allowed_drivers=["ROI_PAC"])
ds_band = ds.GetRasterBand(2) 
std1 = ds_band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
std = 1./ (std1+1e-5)
maxstd = 20

ds = gdal.OpenEx(file2, allowed_drivers=["ROI_PAC"])
ds_band1 = ds.GetRasterBand(1)
ds_band2 = ds.GetRasterBand(2)
dem = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize) 
graddem = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize) 
maxgrad = np.nanmax(graddem)

out = ((graddem / maxgrad)*0.5 + std[:ds.RasterYSize,:ds.RasterXSize]/maxstd) * maxgrad

# création du nouveau fichier
drv = gdal.GetDriverByName('ROI_PAC')
dst_ds = drv.Create(outfile, ds.RasterXSize, ds.RasterYSize, 2, gdal.GDT_Float32)
dst_band1 = dst_ds.GetRasterBand(1)
dst_band2 = dst_ds.GetRasterBand(2)
dst_band1.WriteArray(dem)
dst_band2.WriteArray(out)
del dst_ds, ds

fig = plt.figure(figsize=(16,4))
fig.subplots_adjust(hspace=0.5)
ax1 = fig.add_subplot(1,3,1)
cax = ax1.imshow(std/maxstd, cmap = cmap, vmin=0.1, vmax=0.3 , extent=None,interpolation='nearest')
divider = make_axes_locatable(ax1)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
setp( ax1.get_xticklabels(), visible=False)
ax1.set_title('normalized std')

ax2 = fig.add_subplot(1,3,2)
cax = ax2.imshow(graddem/maxgrad, cmap = cmap, extent=None,interpolation='nearest')
setp( ax2.get_xticklabels(), visible=False)
ax2.set_title('normalized grad DEM')
divider = make_axes_locatable(ax2)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

ax3 = fig.add_subplot(1,3,3)
cax = ax3.imshow(out, cmap = cmap, vmax=np.nanpercentile(out,90) ,extent=None,interpolation='nearest')
setp( ax3.get_xticklabels(), visible=False)
divider = make_axes_locatable(ax3)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
fig.tight_layout()

plt.show()
