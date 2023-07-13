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
plot_geots.py
-------------
Plot georeferenced time series maps additional cropping, conversion options 

Usage: plot_geots.py [--vmin=<value>] [--vmax=<value>] [--wrap=<values>] [--lectfile=<path>] \
[--images=<path>]  [--crop=<values>] [--geocrop=<values>]  \
[--cpt=<values>] [--plot=<yes/no>] [--rad2mm=<value>] [--output=<path>]

Options:
-h --help           Show this screen.
--lectfile PATH     Path of the lect.in file [default: lect_geo.in]
--images PATH       Path to image_retuenues file [default: images_retenues]
--crop VALUE        Crop in radar coordinates [default: 0,nlign,0,ncol]
--crop VALUE        Crop in geo coordiantes [default: 0,nlign,0,ncol]
--vmax              Max colorscale 
--vmin              Min colorscale 
--wrap  VALUE       Wrapped phase between value [default: no]
--cpt               Indicate colorscale
--rad2mm            Scaling value between input data (rad) and desired output [default: -4.4563]
--plot              Display results [default: yes]
--output            Optinional saving as mp4
"""


from osgeo import gdal
import osr, os
gdal.UseExceptions()

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from pylab import *
from mpl_toolkits.basemap import Basemap  
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

if arguments["--vmax"] is not  None:
    vmax = float(arguments["--vmax"])

if arguments["--vmin"] is not  None:
    vmin = float(arguments["--vmin"])

if arguments["--wrap"] is not None:
    vmax=float(arguments["--wrap"])
    vmin=-vmax

# read lect.in 
ncol, nlign = map(int, open(lecfile).readline().split(None, 2)[0:2])

if arguments["--crop"] ==  None:
    ibeg,iend,jbeg,jend = 0,nlign,0,ncol
else:
    crop = map(float,arguments["--crop"].replace(',',' ').split())
    ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])



# load images_retenues file
nb,idates,dates,base=np.loadtxt(fimages, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
# nb images
N=len(dates)
print('Number images: ', N)

# plot diplacements maps
fig, ax = plt.subplots(1)

# LOOK for Nan on the last date
infile = 'geo_'+str(idates[-1])+'_'+str(N-1)+'.unw'
ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])
ds_band2 = ds.GetRasterBand(2)
los = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)*rad2mm
index = np.nonzero(los==0)

# plot basemap
ds_geo=ds.GetGeoTransform()
pix_az, pix_rg = np.indices((ds.RasterYSize,ds.RasterXSize))
lat,lon = ds_geo[3]+ds_geo[5]*pix_az, ds_geo[0]+ds_geo[1]*pix_rg
minx,maxx,maxy,miny = ds_geo[0], ds_geo[0]+ ds_geo[1]*ds.RasterXSize, ds_geo[3], ds_geo[3]+ds_geo[5]*ds.RasterYSize  

if arguments["--geocrop"] is not  None:
    geocrop = list(map(float,arguments["--geocrop"].replace(',',' ').split()))
    latbeg,latend,lonbeg,lonend = float(geocrop[0]),float(geocrop[1]),float(geocrop[2]),float(geocrop[3])
else:
    latbeg,latend,lonbeg,lonend = miny,maxy,minx,maxx
    # print  latbeg,latend,lonbeg,lonend 

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
m.arcgisimage(service='World_Shaded_Relief', xpixels = 1000,zorder=1)


def f(i):
    global idates

    # plt.title(idates[i])
    infile = 'geo_'+str(idates[i])+'_'+str(i)+'.unw'
    rscfile = 'geo_'+str(idates[i])+'_'+str(i)+'.unw.rsc'
    print('Read image {}: {}'.format(i,infile))

    ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])
    ds_band2 = ds.GetRasterBand(2)
    los = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)*rad2mm
    los[index] = float('NaN')

    if arguments["--wrap"] is not None:
        los = np.mod(los+float(arguments["--wrap"]),2*float(arguments["--wrap"]))-float(arguments["--wrap"])
        vmax=float(arguments["--wrap"])
        vmin=-vmax
   
    # masked_array = np.ma.array(los, mask=np.isnan(los))     
    return los

# Initialize
i = 0
im = ax.imshow(f(i),extent=(minx,maxx,miny,maxy),cmap=cmap,\
        vmax=vmax,vmin=vmin,alpha=0.8, zorder=2)    
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=c)
# fig.colorbar(im, orientation='vertical',aspect=10)
m.drawparallels(np.arange(latbeg,latend,.5),linewidth=0.05,dashes=[1, 0],zorder=3)
m.drawmeridians(np.arange(lonbeg,lonend,.5),linewidth=0.05,dashes=[1, 0],zorder=3)
m.drawparallels(np.arange(latbeg,latend,.5),labels=[1,0,0,0],zorder=3)
m.drawmeridians(np.arange(lonbeg,lonend,.5),labels=[0,0,0,1],zorder=3)
# plt.show()
# sys.exit()

# Animation update function
def updatefig(frame):
    global N
    frame = frame % N 
    im.set_array(f(frame))
    return im,

# Animate
ani = animation.FuncAnimation(fig, updatefig, frames=range(1,2*N),  interval=600, blit=True)

if arguments["--output"] is not None:
    base = os.path.splitext(arguments["--output"])[0]
    ani.save(base+'.mp4')

if plot == 'yes':
    plt.show()
