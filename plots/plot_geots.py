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
plot_geots.py
-------------
Plot georeferenced time series maps additional cropping, conversion options 

Usage: plot_geots.py --vmin=<value> --vmax=<value> [--lectfile=<path>] \
[--images=<path>]  [--crop=<values>] [--geocrop=<values>] [--dem=<values>] \
[--cpt=<values>] [--plot=<yes/no>] [--rad2mm=<value>]

Options:
-h --help           Show this screen.
--lectfile PATH     Path of the lect.in file [default: lect_geo.in]
--images PATH       Path to image_retuenues file [default: images_retenues]
--crop VALUE        Crop in radar coordinates [default: 0,nlign,0,ncol]
--crop VALUE        Crop in geo coordiantes [default: 0,nlign,0,ncol]
--vmax              Max colorscale 
--vmin              Min colorscale 
--cpt               Indicate colorscale
--dem               path to DEM file
--rad2mm                Scaling value between input data (rad) and desired output [default: -4.4563]
--plot              Display results [default: yes] 
"""

# gdal
import gdal,osr
gdal.UseExceptions()

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from pylab import *
from mpl_toolkits.basemap import Basemap  

import docopt
arguments = docopt.docopt(__doc__)

if arguments["--lectfile"] ==  None:
   lecfile = "lect_geo.in"
else:
   lecfile = arguments["--lectfile"]

if arguments["--images"] ==  None:
   fimages = "images_retenues"
else:
   fimages = arguments["--images"]

if arguments["--cpt"] is  None:
    cmap=cm.rainbow
else:
    cmap=arguments["--cpt"]

if arguments["--plot"] ==  None:
    plot = 'yes'
else:
    plot = arguments["--plot"]

if arguments["--rad2mm"] ==  None:
    rad2mm = -4.4563 # toard = postive
else:
    rad2mm = float(arguments["--rad2mm"])

vmin = np.float(arguments["--vmin"])
vmax = np.float(arguments["--vmax"])

# read lect.in 
ncol, nlign = map(int, open(lecfile).readline().split(None, 2)[0:2])

if arguments["--crop"] ==  None:
    ibeg,iend,jbeg,jend = 0,nlign,0,ncol
else:
    crop = map(float,arguments["--crop"].replace(',',' ').split())
    ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])

if arguments["--geocrop"] is not  None:
    geocrop = map(float,arguments["--geocrop"].replace(',',' ').split())
    latbeg,latend,lonbeg,lonend = float(geocrop[0]),float(geocrop[1]),float(geocrop[2]),float(geocrop[3])

# load images_retenues file
nb,idates,dates,base=np.loadtxt(fimages, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
# nb images
N=len(dates)
print 'Number images: ', N

# plot diplacements maps
# fig = plt.figure(1,figsize=(14,10))
fig = plt.figure(1,figsize=(18,16))
fig.subplots_adjust(wspace=0.001)


# LOOK for Nan on the last date
infile = 'geo_'+str(idates[-1])+'_'+str(N-1)+'.unw'
ds = gdal.Open(infile, gdal.GA_ReadOnly)
ds_band2 = ds.GetRasterBand(2)
los = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)*rad2mm
index = np.nonzero(los==0)

for l in xrange((N)):  
#for l in xrange((1)): 
#    l=26

    infile = 'geo_'+str(idates[l])+'_'+str(l)+'.unw'
    rscfile = 'geo_'+str(idates[l])+'_'+str(l)+'.unw.rsc'

    ds = gdal.Open(infile, gdal.GA_ReadOnly)
    ds_band1 = ds.GetRasterBand(1)
    ds_band2 = ds.GetRasterBand(2)
    los = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)*rad2mm
    amp = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
    los[index] = np.float('NaN')

    ds_geo=ds.GetGeoTransform()
    print 'Read infile:', infile
    pix_az, pix_rg = np.indices((ds.RasterYSize,ds.RasterXSize))
    lat,lon = ds_geo[3]+ds_geo[5]*pix_az, ds_geo[0]+ds_geo[1]*pix_rg
    minx,maxx,maxy,miny = ds_geo[0], ds_geo[0]+ ds_geo[1]*ds.RasterXSize, ds_geo[3], ds_geo[3]+ds_geo[5]*ds.RasterYSize  

    # ax = fig.add_subplot(1,1,1)
    ax = fig.add_subplot(4,int(N/4)+1,l+1)

    masked_array = np.ma.array(los, mask=np.isnan(los))
    # cmap.set_bad('white')
    # cmap[...,-1] = 0
    cdem = cm.Greys
    cdem.set_bad('white')
    
    if arguments["--geocrop"] is not None:
        
        # check if proejected
        prj = ds.GetProjection()
        srs=osr.SpatialReference(wkt=prj)
        if srs.GetAttrValue('projcs') is not None:
            projection = 'tmerc'
        else:
            projection = 'cyl'

        m = Basemap(
            projection=projection,\
            llcrnrlon=lonbeg, \
            llcrnrlat=latbeg, \
            urcrnrlon=lonend, \
            urcrnrlat=latend, \
            resolution='i',
            ax=ax,
            )


    if arguments["--dem"] is not None:
        ds2 = gdal.Open(arguments["--dem"], gdal.GA_ReadOnly)
        ds2_geo = ds2.GetGeoTransform()
        # print ds2_geo
        ds2_band = ds2.GetRasterBand(1)
        dem = ds2_band.ReadAsArray(0, 0, ds2.RasterXSize, ds2.RasterYSize)
        dminx,dmaxx,dmaxy,dminy = ds2_geo[0], ds2_geo[0]+ ds2_geo[1]*ds2.RasterXSize, ds2_geo[3], ds2_geo[3]+ds2_geo[5]*ds2.RasterYSize  
        print dminx,dmaxx,dmaxy,dminy
        hax = ax.imshow(dem, extent=(dminx,dmaxx,dminy,dmaxy), cmap=cdem,\
        vmax=np.nanpercentile(dem,98),vmin=np.nanpercentile(dem,2),zorder=1)
    else:
        #hax = ax.imshow(amp, extent=(minx,maxx,miny,maxy), cmap=cdem,\
        #vmax=4500,vmin=2000,alpha=1.,zorder=3)
        m.arcgisimage(service='World_Shaded_Relief', xpixels = 1000,zorder=1)
        

    cax = ax.imshow(masked_array,extent=(minx,maxx,miny,maxy),cmap=cmap,\
        vmax=vmax,vmin=vmin,alpha=1, zorder=2)    
    m.drawparallels(np.arange(latbeg,latend,.5),linewidth=0.05,dashes=[1, 0],zorder=3)
    m.drawmeridians(np.arange(lonbeg,lonend,.5),linewidth=0.05,dashes=[1, 0],zorder=3)
    #m.drawparallels(np.arange(latbeg,latend,.5),labels=[1,0,0,0],zorder=3)
    #m.drawmeridians(np.arange(lonbeg,lonend,.5),labels=[0,0,0,1],zorder=3)
    ax.set_title(idates[l],fontsize=6)

    del ds, ds_band1, ds_band2

fig.tight_layout()
plt.suptitle('Time series maps')
fig.colorbar(cax, orientation='vertical',aspect=10)
fig.savefig('geots.pdf', format='PDF',dpi=300)

if plot == 'yes':
    plt.show()
