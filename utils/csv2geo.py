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

from osgeo import gdal
import numpy as np
import docopt
import os, sys
from os import environ
import matplotlib

if environ["TERM"].startswith("screen"):
    matplotlib.use('Agg')

from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# read input parameters and arguments
arguments = docopt.docopt(__doc__)

print('Read .csv file:', arguments["--csv_file"])
#lat,lon,height,defo,coh  = np.loadtxt(arguments["--csv_file"],comments="#", delimiter=',', unpack=True, dtype='f,f,f,f,f')
latp,lonp,defo,coh  = np.loadtxt(arguments["--csv_file"],comments="#", unpack=True, dtype='f,f,f,f')

print('Plot data')
fig = plt.figure(figsize=(16,4))
fig.subplots_adjust(hspace=0.5)

try:
   from matplotlib.colors import LinearSegmentedColormap
   cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
   cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
   cmap = cmap.reversed()
except:
   cmap=cm.rainbow

vmax = np.nanpercentile(defo,98)
vmin = np.nanpercentile(defo,2)

ax1 = fig.add_subplot(1,2,1)
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
m = cm.ScalarMappable(norm=norm,cmap=cmap)
facel=m.to_rgba(defo[::4])
cax = ax1.scatter(latp[::4],lonp[::4],s = 1., marker='o', color=facel) 
ax1.set_title(arguments["--csv_file"])

#plt.show()

# open ref
print('Read reference raster file:', arguments["--ref_file"])
ds = gdal.Open(arguments["--ref_file"])
ds_geo = ds.GetGeoTransform()
proj = ds.GetProjection()
driver = gdal.GetDriverByName('GTiff')
print("> Driver:     ", ds.GetDriver().ShortName)
print("> Size:       ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
print("> Origin:     ", ds_geo[0], ds_geo[3] )
print("> Pixel Size: ", ds_geo[1], ds_geo[5])

# extract lines,cols,lat,lon raster
pix_lin, pix_col = np.indices((ds.RasterYSize,ds.RasterXSize))
lats,lons = ds_geo[3]+ds_geo[5]*pix_lin, ds_geo[0]+ds_geo[1]*pix_col

# initiate new band
los = np.zeros((ds.RasterYSize,ds.RasterXSize))

print('Convert csv to raster...this may take a while...')
for (lat,lon,i,j) in zip(lats.flatten(),lons.flatten(),pix_lin.flatten(),pix_col.flatten()):
#   for (lon,j) in zip(lons.flatten(),pix_col.flatten()): 
    
     #print(lat,lon)
     #print(i,j)
     if ((i % 5) == 0) and (j==0):
         print('Processing line: {}'.format(i))

     # initialisation
     m = float('NaN')
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
     
     #m = np.nanmean(defo[index])
     m =  np.nansum(defo[index]*coh[index]) / np.nansum(coh[index])
     #try:
       #m = defo[index][np.argmax(coh[index])]
     #except:
      # pass
     los[i,j] = m 

outfile = os.path.splitext(arguments["--csv_file"])[0] + '.tiff'
print('Convert output file:', outfile)
dst_ds = driver.Create(outfile, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Float32)
dst_band2 = dst_ds.GetRasterBand(1)
dst_band2.WriteArray(los,0,0)
dst_ds.SetGeoTransform(ds_geo)
dst_ds.SetProjection(proj)
dst_band2.FlushCache()
del dst_ds

print('Plot results')

ax2 = fig.add_subplot(1,2,2)
cax = ax2.imshow(los, cmap = cmap, vmax=vmax, vmin=vmin, extent=None,interpolation='nearest')
divider = make_axes_locatable(ax2)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
ax2.set_title(outfile)

plt.show()
