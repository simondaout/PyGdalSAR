#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT
################################################################################

"""
count_unw_pixel.py
========================

This script count the number of unw pixels for qll ifgs with a given temporal baselines thresholds 

Usage:
        count_unw_pixel.py [--int_path=<path>] [--int_list=<path>] [--dates_list=<path>]  [--lectfile=<path>]  [--prefix=<value>] [--suffix=<value>] [--count=<yes/no>] 

Options:
  --int_path=<dir>    path to interferograms directory [default: LN_DATA]
  --dates_list=<file>   Path to text file containing date,date_dec,bt,bp [default: list_dates]
  --int_list=<file>     Text file containing list of interferograms dates in two colums, $data1 $date2 [default: list_pair]
  --lectfile=<path>       Path to the lect.in file. Simple text file containing width and length and number of images of the time series cube (output of invers_pixel). By default the program will try to find an .hdr file. [default: lect.in].
  --prefix=<value>    Prefix name $prefix$date1-$date2$suffix.unw [default: '']
  --suffix=<vaue>     Suffix name $prefix$date1-$date2$suffix.unw [default: _pre_inv]
  --count<yes/no>     if no, open count_pixel_cube and do not count number of ifg 
  -h --help           Show this screen
"""

from os import environ
import glob, math, os, sys
from osgeo import gdal
import numpy as np
import docopt
gdal.UseExceptions()
import shutil
from datetime import datetime
import matplotlib
if environ["TERM"].startswith("screen"):
        matplotlib.use('Agg')
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

summ = np.zeros((nlines,ncols))
count = np.zeros((nlines,ncols))
for k in range(kmax):
     date1 = date2dec(date_1[k])[0]
     date2 = date2dec(date_2[k])[0]
     # check if Bt > Bc
     Bt = date2 - date1
     # Open ifg
     print('Open ifg {}'.format(str(prefix) + str(date_1[k]) + '-' + str(date_2[k]) + str(suffix) + '.unw'))
     try:
        infile = intdir +  str(date_1[k]) + '-' + str(date_2[k]) + '_pre_inv.unw'
        ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])
     except:
        folder =  'int_'+ str(date_1[k]) + '_' + str(date_2[k]) + '/'
        infile= intdir + str(prefix) + str(date_1[k]) + '-' + str(date_2[k]) + str(suffix) + '.unw'
        ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])
     ds_band2 = ds.GetRasterBand(2)
     los_map = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
     # loop over all lines
     for j in range(0,nlines): # maps[i,j]
         index = np.flatnonzero(np.logical_and(~np.isnan(los_map[j,:]),los_map[j,:] != 0))
         summ[j,index] +=  Bt
         count[j,index] +=  1
     del ds, ds_band2, los_map
average = summ/count
  
# save r4
fid2 = open('average_bt.r4','wb')
average.flatten().astype('float32').tofile(fid2) 

# plot results
fig = plt.figure(0)
ax = fig.add_subplot(1,1,1)
cax = ax.imshow(average,cmap=cm.jet,vmax=np.percentile(average,98))
ax.set_title('Average Bt for all pixels')
plt.setp( ax.get_xticklabels(), visible=False)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
fig.savefig('{}.eps'.format('average_bt'), format='EPS',dpi=150)

