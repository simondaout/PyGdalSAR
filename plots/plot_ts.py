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
plot_ts.py
-------------
Plot a time series file (cube in binary format) 

Usage: plot_ts.py --infile=<path> [--vmin=<value>] [--vmax=<value>] [--lectfile=<path>] [--imref=<value>] \
[--list_images=<path>] [--crop=<values>] [--cpt=<values>] [--dateslim=<values>] 

Options:
-h --help           Show this screen.
--infile PATH       path to time series (depl_cumule)
--lectfile PATH     Path of the lect.in file [default: lect.in]
--imref VALUE       Reference image number [default: 1]
--list_images PATH  Path to image_retuenues file [default: images_retenues]
--crop VALUE        Crop option [default: 0,nlines,0,ncol]
--vmax              Max colorscale [default: 98th percentile]
--vmin              Min colorscale [default: 2th percentile]
--dateslim              Datemin,Datemax time series
"""

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided

import os
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable


# scipy
import scipy
import scipy.optimize as opt
import scipy.linalg as lst
from osgeo import gdal
from datetime import datetime as datetimes
import datetime
import docopt
arguments = docopt.docopt(__doc__)


def date2dec(dates):
    dates  = np.atleast_1d(dates)
    times = []
    for date in dates:
        x = datetimes.strptime('{}'.format(date),'%Y%m%d')
        dec = float(x.strftime('%j'))/365.1
        year = float(x.strftime('%Y'))
        times.append(year + dec)
    return times

infile = arguments["--infile"]
if arguments["--lectfile"] ==  None:
   lecfile = "lect.in"
else:
   lecfile = arguments["--lectfile"]

if arguments["--imref"] !=  None:
    if arguments["--imref"] < 1:
        print('--imref must be between 1 and Nimages')
    else:
        imref = int(arguments["--imref"]) - 1

if arguments["--cpt"] is  None:
    # cmap=cm.jet 
    try:
        from matplotlib.colors import LinearSegmentedColormap
        cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
        cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
        cmap = cmap.reversed()
    except:
        cmap=cm.rainbow
else:
    cmap=arguments["--cpt"]

# lect cube
try:
    ds = gdal.Open(infile)
except:
    passs
if not ds:
  print('.hdr file time series cube {0}, not found, open {1}'.format(infile,lecfile))
  # read lect.in 
  ncol, nlines = map(int, open(lecfile).readline().split(None, 2)[0:2])
else:
  ncol, nlines = ds.RasterXSize, ds.RasterYSize

if arguments["--crop"] ==  None:
    crop = [0,nlines,0,ncol]
else:
    crop = list(map(float,arguments["--crop"].replace(',',' ').split()))
ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])

if arguments["--list_images"] ==  None:
    listim = "images_retenues"
else:
    listim = arguments["--list_images"]

nb,idates,dates,base=np.loadtxt(listim, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
N = len(dates)

if arguments["--dateslim"] is not  None:
    dmin,dmax = arguments["--dateslim"].replace(',',' ').split()
    datemin = date2dec(dmin)
    datemax = date2dec(dmax)
else:
    datemin, datemax = int(np.min(dates)), int(np.max(dates))+1
    dmax = str(datemax) + '0101'
    dmin = str(datemin) + '0101'

cubei = np.fromfile(infile,dtype=np.float32)
cube = as_strided(cubei[:nlines*ncol*N])
kk = np.flatnonzero(np.logical_or(cube==9990, cube==9999))
cube[kk] = float('NaN')
cube[cube==0] = float('NaN')
print('Number of line in the cube: ', cube.shape)
maps = cube.reshape((nlines,ncol,N))
print('Reshape cube: ', maps.shape)
if arguments["--imref"] !=  None:
    cst = np.copy(maps[:,:,imref])
    for l in range((N)):
        maps[:,:,l] = maps[:,:,l] - cst
        if l != imref:
            index = np.nonzero(maps[:,:,l]==0.0)
            maps[:,:,l][index] = float('NaN')

# clean dates
indexd = np.flatnonzero(np.logical_and(dates<datemax,dates>datemin))
nb,idates,dates,base = nb[indexd],idates[indexd],dates[indexd],base[indexd]
# new number of dates
N = len(dates)
maps = as_strided(maps[:,:,indexd])
print('Reshape cube: ', maps.shape)

if arguments["--vmax"] ==  None:
    vmax = np.nanpercentile(maps, 98)*4.4563
else:
    vmax = float(arguments["--vmax"])

if arguments["--vmin"] ==  None:
    vmin = np.nanpercentile(maps, 2)*4.4563 
else:
    vmin = float(arguments["--vmin"])

# plot diplacements maps
fig = plt.figure(1,figsize=(14,10))
fig.subplots_adjust(wspace=0.001)

# vmax = np.abs([np.nanmedian(maps[:,:,-1]) + 1.*np.nanstd(maps[:,:,-1]),\
#     np.nanmedian(maps[:,:,-1]) - 1.*np.nanstd(maps[:,:,-1])]).max()
# vmin = -vmax

# for l in range((N)):
# i = 1
# listdate=[1,4,8,74]
# for l in listdate:
for l in range((N)): 
    d = as_strided(maps[ibeg:iend,jbeg:jend,l])*4.4563
    ax = fig.add_subplot(4,int(N/4)+1,l+1)
    # ax = fig.add_subplot(1,4,i)
    # i = i+1
    cmap.set_bad('white')
    cax = ax.imshow(d,cmap=cmap,vmax=vmax,vmin=vmin)
    ax.set_title(idates[l],fontsize=6)
    plt.setp( ax.get_xticklabels(), visible=False)
    plt.setp( ax.get_yticklabels(), visible=False)

plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)
plt.suptitle('Time series maps')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
fig.savefig('maps_clean.pdf', format='PDF',dpi=150)

plt.show()

