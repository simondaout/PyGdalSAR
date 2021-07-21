#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT
################################################################################

"""
count_unw_pixel.py
========================

This script count the number of unw ifgs for a given pixel for various temporal baselines thresholds 

Usage:
        count_unw_pixel.py [--int_path=<path>] [--int_list=<path>] [--dates_list=<path>]  [--lectfile=<path>]  [--prefix=<value>] [--suffix=<value>]  [--Bc=<values>] [--count=<yes/no>] 

Options:
  --int_path=<dir>    path to interferograms directory [default: LN_DATA]
  --dates_list=<file>   Path to text file containing date,date_dec,bt,bp [default: list_dates]
  --int_list=<file>     Text file containing list of interferograms dates in two colums, $data1 $date2 [default: list_pair]
  --lectfile=<path>       Path to the lect.in file. Simple text file containing width and length and number of images of the time series cube (output of invers_pixel). By default the program will try to find an .hdr file. [default: lect.in].
  --prefix=<value>    Prefix name $prefix$date1-$date2$suffix.unw [default: '']
  --suffix=<vaue>     Suffix name $prefix$date1-$date2$suffix.unw [default: _pre_inv]
  --Bc=<value>        Critical temporal baseline (eg. 0.1).
  --count<yes/no>     if no, open count_pixel_cube and do not count number of ifg 
  -h --help           Show this screen
"""

import glob, math, os, sys
from osgeo import gdal
import numpy as np
import docopt
gdal.UseExceptions()
import shutil
from datetime import datetime
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy.lib.stride_tricks import as_strided

arguments = docopt.docopt(__doc__)
if arguments["--int_path"] == None:
    int_path = 'LN_DATA'
else:
    int_path=arguments["--int_path"]
if arguments["--int_list"] == None:
    int_list = 'list_pair'
else:
    int_list=arguments["--int_list"]
if arguments["--dates_list"] == None:
    baseline = 'list_dates'
else:
    baseline=arguments["--dates_list"]

if arguments["--prefix"] == None:
  prefix = "" 
else:
  prefix=arguments["--prefix"]
if arguments["--suffix"] == None:
  suffix = '_pre_inv'
else:
  suffix=arguments["--suffix"]
if (arguments["--Bc"] == None):
   bc = 0.1
else:
   bc = float(arguments["--Bc"])
if arguments["--lectfile"] ==  None:
   arguments["--lectfile"] = "lect.in"

run = True
if arguments["--count"] ==  'no':
   run = False
    

ncols, nlines = list(map(int, open(arguments["--lectfile"]).readline().split(None, 2)[0:2]))

def date2dec(dates):
    dates  = np.atleast_1d(dates)
    times = []
    for date in dates:
        x = datetime.strptime('{}'.format(date),'%Y%m%d')
        dec = float(x.strftime('%j'))/365.1
        year = float(x.strftime('%Y'))
        # print date,dec,year
        times.append(year + dec)
    return times

intdir = os.path.join(arguments["--int_path"]) + '/'

# read int
date_1,date_2=np.loadtxt(int_list,comments="#",unpack=True,dtype='i,i')
kmax=len(date_1)
print("number of interferogram: ",kmax)
# open baseline.rsc
im,imd,bt,bp=np.loadtxt(baseline,comments="#",usecols=(0,1,2,3),unpack=True,dtype='i,f,f,f')
print("image list=",baseline)
nmax=len(imd)
print("number of image: ",nmax)

if run is True:
  # loop over all images
  count = np.zeros((nlines,ncols,nmax))
  for n in range(nmax):
    print()
    print('Open image {}'.format(im[n]))
    date = imd[n]
    # loop over all ifgs
    for k in range(kmax):
        date1 = date2dec(date_1[k])[0]
        date2 = date2dec(date_2[k])[0]
        # check if Bt > Bc
        Bt = date2 - date1
        if (im[n] == date_1[k] or date == date_2[k]) and (Bt >= bc):
                # Open ifg
                folder =  'int_'+ str(date_1[k]) + '_' + str(date_2[k]) + '/'
                infile= intdir + folder + str(prefix) + str(date_1[k]) + '-' + str(date_2[k]) + str(suffix) + '.unw'
                print('Open ifg {}'.format(str(prefix) + str(date_1[k]) + '-' + str(date_2[k]) + str(suffix) + '.unw'))
                ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])
                ds_band2 = ds.GetRasterBand(2)
                los_map = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
                # loop over all lines
                for j in range(0,nlines): # maps[i,j]
                    index = np.flatnonzero(np.logical_and(~np.isnan(los_map[j,:]),los_map[j,:] != 0))
                    count[j,index,n] +=  1
                del ds, ds_band2, los_map
else:
  print('Open count_pixel_cube')  
  count = np.fromfile('count_pixel_cube',dtype=np.float32).reshape(nlines,ncols,nmax)

# loop over count
countf = np.zeros((nlines,ncols))
for n in range(nmax):
    index = np.nonzero(count[:,:,n])
    countf[index] += 1

# plot results
fig = plt.figure(0)
ax = fig.add_subplot(1,1,1)
cax = ax.imshow(countf,cmap=cm.jet,vmax=np.percentile(countf,98))
ax.set_title('Number of images with bt>{}'.format(bc))
plt.setp( ax.get_xticklabels(), visible=False)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
fig.savefig('{}.eps'.format('count_pixels_all_{}'.format(bc)), format='EPS',dpi=150)

# save r4
fid2 = open('count_pixel_all.r4','wb')
countf.flatten().astype('float32').tofile(fid2)

fig = plt.figure(1,figsize=(14,10))
fig.subplots_adjust(wspace=0.01)
for n in range((nmax)):
  ax = fig.add_subplot(10,int(nmax/10)+1,n+1)
  cax = ax.imshow(count[:,:,n],cmap=cm.jet,vmax=np.percentile(count,98))
  ax.set_title(im[n],fontsize=6)
  plt.setp(ax.get_xticklabels(), visible=False)
  plt.setp(ax.get_yticklabels(), visible=False)
fig.savefig('{}.eps'.format('count_pixels_{}'.format(bc)), format='EPS',dpi=150)
plt.show()

# save cube
fid = open('count_pixel_cube', 'wb')
count.flatten().astype('float32').tofile(fid)


