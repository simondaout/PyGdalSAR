#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Simon DAOUT (CRPG-ENSG)
############################################


"""\
plot_geo.py
-------------
Plot georeferenced file in RMG or Tif format with optional DEM or ARCGISIMAGE

Usage: plot_geots.py --infile=<value> [--vmin=<value>] [--vmax=<value>] \
[--geocrop=<values>] [--wrap=<values>] [--cpt=<values>] [--dem=<values>] \
[--coeff=<values>] [--plot=<yes/no>] [--rad2mm=<value>] [--outfile=<name>] \
[--shapefile=<file>]

Options:
-h --help           Show this screen.
--infile VALUE
--geocrop VALUE     Crop in geo coordiantes [default: latmin,latmax,lonmin,lonmax]
--vmax              Max colorscale [default: 90th percentile]
--vmin              Min colorscale [default: 10th percentile]
--cpt               Indicate colorscale
--wrap  VALUE       Wrapped phase between value [default: no]
--dem               Path to DEM file
--coeff             Optional value to scale data 
--plot              Display results [default: yes] 
--rad2mm            Scaling value between input data (rad) and desired output [default: -4.4563]
--outfile           Give a name to output file
--shapefile         Plot shape file
"""

print()
print()
print('Author: Simon Daout')
print()
print()

# gdal
from osgeo import gdal, osr
gdal.UseExceptions()
import os
# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from pylab import *
import contextily as ctx
from mpl_toolkits.axes_grid1 import make_axes_locatable
import shapefile as shp  # Requires the pyshp package
from shapely.geometry import box
import geopandas as gpd


import docopt, os
arguments = docopt.docopt(__doc__)

infile = arguments["--infile"]

if arguments["--cpt"] is  None:
    try:
        from matplotlib.colors import LinearSegmentedColormap
        cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
        cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
        cmap = cmap.reversed()
    except:
        cmap=cm.rainbow
else:
    cmap=arguments["--cpt"]

if arguments["--coeff"] is  None:
  vel2disp=1
else:
  vel2disp=float(arguments["--coeff"])

if arguments["--plot"] ==  None:
    plot = 'yes'
else:
    plot = arguments["--plot"]

if arguments["--rad2mm"] ==  None:
    rad2mm = -4.4563 # toward = postive
else:
    rad2mm = float(arguments["--rad2mm"])  

ds_extension = os.path.splitext(infile)[1]
# print(ds_extension)  

if ds_extension == ".unw":
  ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])
  ds_band1 = ds.GetRasterBand(1)
  ds_band2 = ds.GetRasterBand(2)
  los = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)*rad2mm*vel2disp
  amp = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
  los[amp==0] = float('NaN')

if (ds_extension == ".tif") or (ds_extension == ".tiff"):
  ds = gdal.Open(infile, gdal.GA_ReadOnly)
  ds_band1 = ds.GetRasterBand(1)
  los = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)*rad2mm*vel2disp

los[los==255] = float('NaN')

ds_geo=ds.GetGeoTransform()
print('Read infile:', infile)
pix_az, pix_rg = np.indices((ds.RasterYSize,ds.RasterXSize))
lat,lon = ds_geo[3]+ds_geo[5]*pix_az, ds_geo[0]+ds_geo[1]*pix_rg
minx,maxx,maxy,miny = ds_geo[0], ds_geo[0]+ ds_geo[1]*ds.RasterXSize, ds_geo[3], ds_geo[3]+ds_geo[5]*ds.RasterYSize  
print('Original coordinates:', miny,maxy,minx,maxx)
# # print lat
# sys.exit()

if arguments["--geocrop"] is not  None:
    geocrop = list(map(float,arguments["--geocrop"].replace(',',' ').split()))
    latbeg,latend,lonbeg,lonend = float(geocrop[0]),float(geocrop[1]),float(geocrop[2]),float(geocrop[3])
else:
    latbeg,latend,lonbeg,lonend = miny,maxy,minx,maxx
print('Cooridnates plot:', latbeg,latend,lonbeg,lonend)

los[los==0.]=float('NaN')
kk = np.nonzero(np.logical_or(np.logical_or(~np.isnan(los), np.abs(los)<999.),los==0.0))
mprim = los[kk]

if arguments["--vmax"] ==  None:
  vmax = np.nanpercentile(mprim, 90)
else:
  vmax = float(arguments["--vmax"])

if arguments["--vmin"] ==  None:
  vmin = np.nanpercentile(mprim, 10)
else:
  vmin = float(arguments["--vmin"])

if arguments["--wrap"] is not None: 
  los = np.mod(los+float(arguments["--wrap"]),2*float(arguments["--wrap"]))-float(arguments["--wrap"])
  # los = los - nanmean(los)
  vmax=float(arguments["--wrap"])
  vmin=-vmax

# plot diplacements maps
fig = plt.figure(1,figsize=(8,5))
fig.subplots_adjust(wspace=0.001)
ax = fig.add_subplot(1,1,1)
basename = os.path.splitext(infile)[0]
cdem = cm.Greys
cdem.set_bad('white')
# cmap.set_bad('white')
masked_array = np.ma.array(los, mask=np.isnan(los))

# définir la limite de la carte
ax.axis((lonbeg,lonend,latbeg,latend))

# Affichez la carte de fond
ctx.add_basemap(ax,crs=4326 ,source=ctx.providers.Esri.WorldShadedRelief)

if arguments["--dem"] is not None:
   ds2 = gdal.Open(arguments["--dem"], gdal.GA_ReadOnly)
   ds2_geo = ds2.GetGeoTransform()
   ds2_band = ds2.GetRasterBand(1)
   dem = ds2_band.ReadAsArray(0, 0, ds2.RasterXSize, ds2.RasterYSize)
   dminx,dmaxx,dmaxy,dminy = ds2_geo[0], ds2_geo[0]+ ds2_geo[1]*ds2.RasterXSize, ds2_geo[3], ds2_geo[3]+ds2_geo[5]*ds2.RasterYSize  
   hax = ax.imshow(dem, extent=(dminx,dmaxx,dminy,dmaxy), cmap=cdem,\
    vmax=255,vmin=1,zorder=3)

cax = ax.imshow(masked_array,extent=(minx,maxx,miny,maxy),cmap=cmap,\
     vmax=vmax,vmin=vmin, alpha=0.7, zorder=4,interpolation='none')

if arguments["--shapefile"] is not None and os.path.exists(arguments["--shapefile"]):
    shape = gpd.read_file(arguments["--shapefile"])
    shape = shape.to_crs(4326)
    shape.plot(ax=ax,facecolor='none', color='none',edgecolor='black',zorder=10)
    
ax.set_xlim([lonbeg,lonend])
ax.set_ylim([latbeg,latend])
ax.set_xticks(np.linspace(lonbeg,lonend,3))
ax.set_yticks(np.linspace(latbeg,latend,3))

if arguments["--outfile"] ==  None:
    outfile = basename+'.pdf'
else:
    outfile = arguments["--outfile"]+'.pdf'
del ds, ds_band1
plt.suptitle(outfile)
try:
  divider = make_axes_locatable(ax)
  c = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(cax, cax=c)
except:
  pass
fig.savefig(outfile, format='PDF',dpi=180)

if plot == 'yes':
    plt.show()
