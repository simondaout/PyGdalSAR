#!/usr/bin/env python3
# -*- coding: utf-8 -*-
################################################################################
#
# NSBAS - New Small Baseline Chain
#
################################################################################
# Author        : Simon DAOUT (Oxford)
################################################################################


"""\
correct_flatr_unw.py
-------------
Add ramp correction to unrapped file

usage: correct_flatr_unw.py --infile=<path> --param=<path> --outfile=<path> [--rlook_factor=<path>] [--plot=<yes/no>]

--infile=<path>           Unwrapped IFG to be corrected from ramp in range
--param=<path>            Parameter text file .flatr containing containing polynomial fit in range (.flatr) or azimuth (.flata)
--outfile=<outfile>       Prefix name $prefix$date1-$date2$suffix_$rlookrlks.unw
--rlook_factor=<value>    Look factor between wrapped correction and unwrapped file [default: 1]
--plot=<yes/no>           If yes, plot figures for each ints [default: no]
"""


from os import path, environ
import os
import matplotlib
if environ["TERM"].startswith("screen"):
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from datetime import datetime
import datetime
import time
import numpy as np
from numpy.lib.stride_tricks import as_strided
from osgeo import gdal
gdal.UseExceptions()
from nsbas import docopt
import shutil

# read arguments
arguments = docopt.docopt(__doc__)

infile = arguments["--infile"]
rscfile = infile + '.rsc'
outfile = arguments["--outfile"]
outrsc = outfile + '.rsc'

if arguments["--rlook_factor"] ==  None:
    factor = 1
else:
    factor = int(arguments["--rlook_factor"])

if arguments["--plot"] ==  None:
    plot = 'no'
else:
    plot = str(arguments["--plot"])

ds = gdal.Open(infile, gdal.GA_ReadOnly)
# Get the band that have the data we want
ds_band1 = ds.GetRasterBand(1)
ds_band2 = ds.GetRasterBand(2)

los_map = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
coh_map = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
print()
print('Nlign:{}, Ncol:{}:'.format(ds.RasterYSize, ds.RasterXSize))

param = arguments["--param"]
extension = os.path.splitext(param)[1]
if extension == '.flatr':
    print('Add back range correction...')
    print('Open range file:', param)
    rg_a, rg_b, rg_c, rg_d, rg_f, rg_g = np.loadtxt(param,comments="#",usecols=(0,1,2,3,4,5),unpack=True,dtype='f,f,f,f,f,f')
    rg = np.tile(np.arange(1,ds.RasterXSize+1)*factor,ds.RasterYSize).reshape(ds.RasterYSize, ds.RasterXSize)
    corr_map = rg_a*rg + rg_b*rg**2 + rg_c*rg**3 + rg_d**rg**4 + rg_f*rg**5 + rg_g*rg**6
    print('Add range ramp %f r, %f r**2  + %f r**3 + %f r**4 + %f r**5 + %f r**6'%(rg_a, rg_b, rg_c, rg_d, rg_f, rg_g))
if extension == '.flata':
    print('Add back azimutal correction...')
    print('Open azimutal file:', param)
    az_a, az_b, az_c, az_d, az_f, az_g = np.loadtxt(azfile,comments="#",usecols=(0,1,2,3,4,5),unpack=True,dtype='f,f,f,f,f,f')
    corr_map = az_a*az + az_b*az**2 + az_c*az**3 + az_d**az**4 + az_f*az**5 + az_g*az**6
    print('Add azimuthal ramp %f az, %f az**2  + %f az**3 + %f az**4 + %f az**5 + %f az**6'%(az_a, az_b, az_c, az_d, az_f, az_g))

flatlos = np.copy(los_map) - corr_map
flatlos[los_map==0] = 0
flatlos[np.isnan(los_map)] = np.float('NaN')

nfigure=0
fig = plt.figure(nfigure,figsize=(11,8))
vmax,vmin = np.nanpercentile(flatlos,95), np.nanpercentile(flatlos,5)

ax = fig.add_subplot(1,3,1)
cax = ax.imshow(los_map,cmap=cm.gist_rainbow,vmax=vmax,vmin=vmin)
ax.set_title('{}'.format(infile))
setp( ax.get_xticklabels(), visible=None)
fig.colorbar(cax, orientation='vertical',aspect=10)

ax = fig.add_subplot(1,3,2)
cax = ax.imshow(corr_map,cmap=cm.gist_rainbow,vmax=vmax,vmin=vmin)
ax.set_title('Unwrapped Ramp')
setp( ax.get_xticklabels(), visible=None)
fig.colorbar(cax, orientation='vertical',aspect=10)

ax = fig.add_subplot(1,3,3)
cax = ax.imshow(flatlos,cmap=cm.gist_rainbow,vmax=vmax,vmin=vmin)
ax.set_title('{0}'.format(outfile))
setp( ax.get_xticklabels(), visible=None)
fig.colorbar(cax, orientation='vertical',aspect=10)
if extension == '.flata':
    fig.savefig('add_az_corrrection.png', format='PNG',dpi=150)
if extension == '.flatr':
    fig.savefig('add_rg_corrrection.png', format='PNG',dpi=150)

if plot=='yes':
    plt.show()

drv = gdal.GetDriverByName("roi_pac")
dst_ds = drv.Create(outfile, ds.RasterXSize, ds.RasterYSize, 2, gdal.GDT_Float32)
dst_band1 = dst_ds.GetRasterBand(1)
dst_band2 = dst_ds.GetRasterBand(2)
dst_band1.WriteArray(coh_map,0,0)
dst_band2.WriteArray(flatlos,0,0)
shutil.copy(rscfile,outrsc)
