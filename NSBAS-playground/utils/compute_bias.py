#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Compute bias
========================

Usage:
    compute_bias.py --long_baseline_int=<path> --list_short_baselines=<path> --radar_file=<file> 

Options:
  --long_baseline_int=<path>       Absolute path to long baseline interferogram 
  --list_short_baselines=<path>    Absolute path to list of short baselines ifgs
  --radar_file FILE   Path to the radar.hgt file
  -h --help           Show this screen
"""

import docopt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
# gdal
from osgeo import gdal
import os
gdal.UseExceptions()

arguments = docopt.docopt(__doc__)
radar=arguments["--radar_file"]
driver = gdal.GetDriverByName("roi_pac")
ds = gdal.OpenEx(radar, allowed_drivers=["ROI_PAC"])
nlines,ncols= ds.RasterYSize, ds.RasterXSize

int_list=arguments["--list_short_baselines"]
date_1,date_2=np.loadtxt(int_list,comments="#",unpack=True,dtype='i,i')
kmax=len(date_1)

# open short baselines
los = np.zeros((nlines,ncols,kmax))
for i in range((kmax)):
    date1,date2 = date_1[i], date_2[i]
    infile = str(date1) + '-' + str(date2) + '_pre_inv.unw'
    ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])
    ds_band2 = ds.GetRasterBand(2)
    los_map = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:nlines,:ncols]
    los_map[los_map==0] = float('NaN') 
    los[:ds.RasterYSize,:ds.RasterXSize,i] = np.copy(los_map) 
    #plt.imshow(los[:,:,i])
    #plt.show()
sum_short = np.nansum(los,axis=2)
del los

# open long baseline
ds = gdal.OpenEx(arguments["--long_baseline_int"], allowed_drivers=["ROI_PAC"])
ds_band2 = ds.GetRasterBand(2)
long_los = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:nlines,:ncols]   
long_los[long_los==0] = float('NaN')

# compute diff
bias = sum_short - long_los[:ds.RasterYSize,:ds.RasterXSize] 
cst = np.nanmean(bias)
bias = bias - cst
sum_short = sum_short -cst

#Plot
try:
    from matplotlib.colors import LinearSegmentedColormap
    cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
    cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
except:
    cmap=cm.rainbow
cmap_r = cmap.reversed()

fig = plt.figure(figsize=(12,5))
ax = fig.add_subplot(1,3,1)
cax =ax.imshow(sum_short,vmax=np.nanpercentile(long_los,98),vmin=np.nanpercentile(long_los,2),cmap=cmap,interpolation=None)
ax.set_title('Sum short baselines')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
plt.setp( ax.get_xticklabels(), visible=False)
plt.setp( ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(1,3,2)
cax = ax.imshow(long_los,vmax=np.nanpercentile(long_los,98),vmin=np.nanpercentile(long_los,2),cmap=cmap,interpolation=None)
ax.set_title('Long baselines ifg')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
plt.setp( ax.get_xticklabels(), visible=False)
plt.setp( ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(1,3,3)
cax =ax.imshow(bias,vmax=np.nanpercentile(bias,85),vmin=np.nanpercentile(bias,15),cmap=cmap,interpolation=None)
ax.set_title('Biases')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
plt.setp( ax.get_xticklabels(), visible=False)
plt.setp( ax.get_yticklabels(), visible=False)

fig.tight_layout()
fig.savefig('biases.pdf', format='PDF',dpi=150)
plt.show()


 
