#!/usr/bin/env python
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
csv2geo.py
-------------
Convert .csv file to geotif.

Usage: csv2geo.py --csv_file=<file> --ref_file=<file>
geo2r4.py -h | --help

Options:
-h --help           Show this screen
--csv_file=<file>   .csv file to convert to raster
--ref_file=<file>   ref raster file
"""

import gdal
import numpy as np
import docopt
import os, sys
from matplotlib import pyplot as plt
import matplotlib.cm as cm

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# read input parameters and arguments
arguments = docopt.docopt(__doc__)

print('Read infile:', arguments["--csv_file"])
#lat,lon,height,defo,coh  = np.loadtxt(arguments["--csv_file"],comments="#", delimiter=',', unpack=True, dtype='f,f,f,f,f')
latp,lonp,defo,coh  = np.loadtxt(arguments["--csv_file"],comments="#", unpack=True, dtype='f,f,f,f')

# open ref
print('Read infile:', arguments["--ref_file"])
ds = gdal.Open(arguments["--ref_file"])
ds_geo = ds.GetGeoTransform()
# (95.201531714, 0.001666666667, 0.0, 39.02513845, 0.0, -0.001666666667)
proj = ds.GetProjection()
driver = gdal.GetDriverByName('GTiff')

pix_az, pix_rg = np.indices((ds.RasterYSize,ds.RasterXSize))
lats,lons = ds_geo[3]+ds_geo[5]*pix_az, ds_geo[0]+ds_geo[1]*pix_rg

# initiate new band
los = np.zeros((ds.RasterYSize,ds.RasterXSize))

print('Convert csv to raster')
for (lat,i) in zip(lats.flatten(),pix_az.flatten()):
   for (lon,j) in zip(lons.flatten(),pix_rg.flatten()): 
     m = np.float('NaN')
     lat_min,lat_max = lat + ds_geo[5], lat
     lon_min,lon_max = lon, lon + ds_geo[1]
     #print(i,j) 
     #print(lat_min,lat_max)
     #print(lon_min,lon_max)
     
     index = np.nonzero(
     np.logical_and(latp>lat_min,
     np.logical_and(latp<lat_max,
     np.logical_and(lonp<lon_max,
     lonp>lon_min))))
     
     # mean
     #m = np.nanmean(defo[index])
     # max coherence
     idx_cohm = coh[index].index(max(coh[index]))
     m = defo[index][idx_cohm]
     los[i,j] = m 

outfile = os.path.splitext(csv_file)[0] + '.tiff'
print('Convert output file:', outfile)
dst_ds = driver.Create(outfile, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Float32)
dst_band2 = dst_ds.GetRasterBand(1)
dst_band2.WriteArray(los,0,0)
dst_ds.SetGeoTransform(ds_geo)
dst_ds.SetProjection(proj)
dst_band2.FlushCache()
del dst_ds

print('Plot results')
fig = plt.figure(figsize=(16,4))
fig.subplots_adjust(hspace=0.5)

try:
   from matplotlib.colors import LinearSegmentedColormap
   cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
   cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
   cmap = cmap.reversed()
except:
   cmap=cm.rainbow

vmax = np.nanpercentile(los,98)
vmin = np.nanpercentile(los,2)

ax1 = fig.add_subplot(1,2,1)
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
m = cm.ScalarMappable(norm=norm,cmap=cmap)
facel=m.to_rgba(defo)
cax = ax1.scatter(lats,lons,s = 1., marker='o', color=facel) 
divider = make_axes_locatable(ax1)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
plt.setp( ax1.get_xticklabels(), visible=False)
ax1.set_title(csv_file)

ax2 = fig.add_subplot(1,2,2)
cax = ax2.imshow(los, cmap = cmap, vmax=vmax, vmin=vmin, extent=None,interpolation='nearest')
plt.setp( ax2.get_xticklabels(), visible=False)
divider = make_axes_locatable(ax2)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
ax2.set_title(outfile)

plt.show()
