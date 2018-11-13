#!/usr/bin/env python2
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
plot_geo.py
-------------
Plot georeferenced map 

Usage: plot_geots.py --infile=<value> [--vmin=<value>] [--vmax=<value>] \
[--geocrop=<values>] [--wrap=<yes/no>] [--cpt=<values>] [--dem=<values>] \
[--coeff=<values>] [--plot=<yes/no>] [--rad2mm=<value>]

Options:
-h --help           Show this screen.
--infile VALUE
--geocrop VALUE     Crop in geo coordiantes [default: latmin,latmax,lonmin,lonmax]
--vmax              Max colorscale [default: 90th percentile]
--vmin              Min colorscale [default: 10th percentile]
--cpt				Indicate colorscale
--wrap 	VALUE		Wrapped phase between value [default: no]
--dem 				path to DEM file
--coeff				optional value to scale data 
--plot              Display results [default: yes] 
--rad2mm            Scaling value between input data (rad) and desired output [default: -4.4563]
"""

# gdal
import gdal
gdal.UseExceptions()
import os
# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided

import matplotlib as mpl
# mpl.rcParams['backend'] = 'TkAgg' 
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from pylab import *
from mpl_toolkits.basemap import Basemap

import docopt
arguments = docopt.docopt(__doc__)

infile = arguments["--infile"]

if arguments["--cpt"] is  None:
	# cmap=cm.jet 
	cmap=cm.rainbow
else:
	cmap=arguments["--cpt"]

if arguments["--coeff"] is  None:
	vel2disp=1
else:
	vel2disp=np.float(arguments["--coeff"])

if arguments["--plot"] ==  None:
    plot = 'yes'
else:
    plot = arguments["--plot"]

if arguments["--rad2mm"] ==  None:
    rad2mm = -4.4563 # toward = postive
else:
    rad2mm = float(arguments["--rad2mm"])  

ds_extension = os.path.splitext(infile)[1]

if ds_extension == ".unw":
	ds = gdal.Open(infile, gdal.GA_ReadOnly)
	ds_band1 = ds.GetRasterBand(1)
	ds_band2 = ds.GetRasterBand(2)
	los = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)*rad2mm*vel2disp
	amp = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
	los[amp==0] = np.float('NaN')
if ds_extension == ".tif":
	ds = gdal.Open(infile, gdal.GA_ReadOnly)
	ds_band1 = ds.GetRasterBand(1)
	los = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)*rad2mm*vel2disp

ds_geo=ds.GetGeoTransform()
print 'Read infile:', infile
pix_az, pix_rg = np.indices((ds.RasterYSize,ds.RasterXSize))
lat,lon = ds_geo[3]+ds_geo[5]*pix_az, ds_geo[0]+ds_geo[1]*pix_rg
minx,maxx,maxy,miny = ds_geo[0], ds_geo[0]+ ds_geo[1]*ds.RasterXSize, ds_geo[3], ds_geo[3]+ds_geo[5]*ds.RasterYSize  
# print minx,maxx,miny,maxy
# # print lat
# sys.exit()

if arguments["--geocrop"] is not  None:
    geocrop = map(float,arguments["--geocrop"].replace(',',' ').split())
    latbeg,latend,lonbeg,lonend = float(geocrop[0]),float(geocrop[1]),float(geocrop[2]),float(geocrop[3])
else:
    latbeg,latend,lonbeg,lonend = minx,maxx,maxy,miny

los[los==0.]=np.float('NaN')
kk = np.nonzero(np.logical_or(np.logical_or(~np.isnan(los), np.abs(los)<999.),los==0.0))
mprim = los[kk]

if arguments["--vmax"] ==  None:
	vmax = np.nanpercentile(mprim, 90)
else:
	vmax = np.float(arguments["--vmax"])

if arguments["--vmin"] ==  None:
	vmin = np.nanpercentile(mprim, 10)
else:
	vmin = np.float(arguments["--vmin"])

if arguments["--wrap"] is not None:	
	los = np.mod(los+float(arguments["--wrap"]),2*float(arguments["--wrap"]))-float(arguments["--wrap"])
	# los = los - nanmean(los)
	vmax=float(arguments["--wrap"])
	vmin=-vmax

# plot diplacements maps
fig = plt.figure(1,figsize=(9,7))
fig.subplots_adjust(wspace=0.001)
ax = fig.add_subplot(1,1,1)
basename = os.path.splitext(infile)[0]

cdem = cm.Greys
cdem.set_bad('white')

# cmap.set_bad('white')
masked_array = np.ma.array(los, mask=np.isnan(los))


m = Basemap(
    projection='cyl',\
    llcrnrlon=lonbeg, \
    llcrnrlat=latbeg, \
    urcrnrlon=lonend, \
    urcrnrlat=latend, \
    resolution='i',
    ax=ax,
    )

# m.drawparallels(np.arange(latbeg,latend,0.25),linewidth=0.25,zorder=1)
# m.drawmeridians(np.arange(lonbeg,lonend,0.25),linewidth=0.25,zorder=1)
# m.drawparallels(np.arange(latbeg,latend,.25),labels=[1,0,0,0],zorder=1)
# m.drawmeridians(np.arange(lonbeg,lonend,.25),labels=[0,0,0,1],zorder=1)
m.fillcontinents(color='white',lake_color='lightblue',zorder=1)


if arguments["--dem"] is not None:
   ds2 = gdal.Open(arguments["--dem"], gdal.GA_ReadOnly)
   ds2_geo = ds2.GetGeoTransform()
   # print ds2_geo
   ds2_band = ds2.GetRasterBand(1)
   dem = ds2_band.ReadAsArray(0, 0, ds2.RasterXSize, ds2.RasterYSize)
   dminx,dmaxx,dmaxy,dminy = ds2_geo[0], ds2_geo[0]+ ds2_geo[1]*ds2.RasterXSize, ds2_geo[3], ds2_geo[3]+ds2_geo[5]*ds2.RasterYSize  
   print dminx,dmaxx,dmaxy,dminy
   print np.nanpercentile(dem,98), np.nanpercentile(dem,2)
   # hax = ax.imshow(dem, extent=(dminx,dmaxx,dminy,dmaxy), cmap=cdem,\
   #  vmax=np.nanpercentile(dem,98),vmin=np.nanpercentile(dem,2),zorder=1)
   hax = ax.imshow(dem, extent=(dminx,dmaxx,dminy,dmaxy), cmap=cdem,\
    vmax=255,vmin=1,zorder=2)
else:
	# pass
   m.arcgisimage(service='World_Shaded_Relief', xpixels = 1000,zorder=2)
   # hax = ax.imshow(amp, extent=(minx,maxx,miny,maxy), cmap=cdem,\
   #  vmax=4500,vmin=2000,alpha=1.,zorder=1)


cax = ax.imshow(masked_array,extent=(minx,maxx,miny,maxy),cmap=cmap,\
     vmax=vmax,vmin=vmin,alpha=0.45, zorder=3) 
ax.set_title(basename,fontsize=6)


outfile = basename+'.pdf'

del ds, ds_band1, ds_band2
fig.tight_layout()
plt.suptitle(outfile)
try:
	fig.colorbar(cax, orientation='vertical',aspect=10)
except:
	pass

# fig.savefig(basename+'.tiff', format='tiff',dpi=180)
fig.savefig(outfile, format='PDF',dpi=180)

if plot == 'yes':
    plt.show()
