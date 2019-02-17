#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
################################################################################
#
# NSBAS - New Small Baseline Chain
#
################################################################################
# Author        : Simon DAOUT (Kiel)
################################################################################


"""\
correct_ramp_unw.py
-------------
Add ramps corrections to unrapped files

usage: correct_ramp_unw.py --int_list=<path> --int_path=<path> \
--prefix=<value> --suffix=<value> --rlook=<value>  \
[--suffix_range_file=<path>] [--preffix_range_file=<path>] \
[--suffix_azimuth_file=<path>] [--preffix_azimuth_file=<path>] \
[--suffix_output=<value>] [--rlook_factor=<path>] [--plot=<yes/no>]

--int_list PATH         Text file containing list of interferograms dates in two colums, $data1 $date2
--int_path PATh         Absolute path to interferograms directory
--prefix VALUE          Prefix name $prefix$date1-$date2$suffix_$rlookrlks.unw
--suffix value          Suffix name $prefix$date1-$date2$suffix_$rlookrlks.unw
--rlook value           look int. $prefix$date1-$date2$suffix_$rlookrlks.unw
--suffix_range_file     Suffix of .flatr range file $prefix_range_file $date1-$date2$suffix_range_file .flatr
--preffix_range_file    Preffix of .flatr range file $prefix_range_file $date1-$date2$suffix_range_file .flatr
--suffix_azimuth_file   Suffix of .flatz azimuth file $prefix_azimuth_file$date1-$date2$suffix_azimuth_file.flatz
--preffix_azimuth_file  Preffix of .flatz azimuth file $prefix_azimuth_file$date1-$date2$suffix_azimuth_file.flatz
--suffix_output value   Suffix output file name $prefix$date1-$date2$suffix_output [default: '']
--rlook_factor  value   Look factor between wrapped correction and unwrapped file 
--plot yes/no           If yes, plot figures for each ints [default: no]
"""


# system
from os import path, environ
import os

# plot
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

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided

# gdal
import gdal
gdal.UseExceptions()

from nsbas import docopt
import shutil

# read arguments
arguments = docopt.docopt(__doc__)

int_list=arguments["--int_list"]
int_path=arguments["--int_path"]
prefix=arguments["--prefix"]
suffix=arguments["--suffix"]
rlook=arguments["--rlook"]

if arguments["--suffix_output"] == None:
    suffout = ''
else:
    suffout = arguments["--suffix_output"]

if arguments["--suffix_azimuth_file"] == None:
    sufaz = ''
else:
    sufaz = arguments["--suffix_azimuth_file"]

if arguments["--preffix_azimuth_file"] == None:
    preaz = ''
else:
    preaz = arguments["--preffix_azimuth_file"]


if arguments["--suffix_range_file"] == None:
    sufrg = ''
else:
    sufrg = arguments["--suffix_range_file"]

if arguments["--preffix_range_file"] == None:
    prerg = ''
else:
    prerg = arguments["--preffix_range_file"]

if arguments["--suffix_azimuth_file"] == None and arguments["--preffix_azimuth_file"] == None and \
arguments["--suffix_range_file"] == None and arguments["--preffix_range_file"] == None:
    print 'No input files for ramp corrections'
    sys.exit()

if arguments["--rlook_factor"] ==  None:
    factor = 1
else:
    factor = int(arguments["--rlook_factor"])
wlook = str(int(rlook)/factor) 


if arguments["--plot"] ==  None:
    plot = 'no'
else:
    plot = str(arguments["--plot"])

print
# read int
date_1,date_2=np.loadtxt(int_list,comments="#",unpack=True,dtype='i,i')
kmax=len(date_1)
print "number of interferogram: ",kmax

for kk in xrange((kmax)):
    date1, date2 = date_1[kk], date_2[kk]
    idate = str(date1) + '-' + str(date2) 
    folder = int_path + 'int_'+ str(date1) + '_' + str(date2) + '/'
    rscfile=folder + prefix + str(date1) + '-' + str(date2) + suffix + '_' + rlook + 'rlks.unw.rsc'
    infile=folder + prefix + str(date1) + '-' + str(date2) + suffix + '_' + rlook + 'rlks.unw'

    if arguments["--suffix_azimuth_file"] == None and arguments["--preffix_azimuth_file"] == None : 
        az_a, az_b, az_c, az_d, az_f, az_g = 0, 0, 0, 0, 0, 0
    else:
        azfile = folder + preaz + str(idate) + sufaz + '_' + wlook + 'rlks.flatr'
        print 'Open azimuth file:', azfile
        az_a, az_b, az_c, az_d, az_f, az_g = np.loadtxt(azfile,comments="#",usecols=(0,1,2,3,4,5),unpack=True,dtype='f,f,f,f,f,f')

    if arguments["--suffix_range_file"] == None and arguments["--preffix_range_file"] == None:
        rg_a, rg_b, rg_c, rg_d, rg_f, rg_g = 0, 0, 0, 0, 0, 0
    else:
        rgfile = folder + prerg + str(idate) + sufrg + '_' + wlook + 'rlks.flatr'
        print 'Open range file:', rgfile
        rg_a, rg_b, rg_c, rg_d, rg_f, rg_g = np.loadtxt(rgfile,comments="#",usecols=(0,1,2,3,4,5),unpack=True,dtype='f,f,f,f,f,f')

    ds = gdal.Open(infile, gdal.GA_ReadOnly)
    # Get the band that have the data we want
    ds_band1 = ds.GetRasterBand(1)
    ds_band2 = ds.GetRasterBand(2)

    los_map = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
    coh_map = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
    print
    print 'Nlign:{}, Ncol:{}, int:{}:'.format(ds.RasterYSize, ds.RasterXSize, idate)

    corr_map = np.ones((ds.RasterYSize, ds.RasterXSize))
    az = np.tile(np.arange(1,ds.RasterYSize+1)*factor,ds.RasterXSize).reshape(ds.RasterXSize,ds.RasterYSize).T
    rg = np.tile(np.arange(1,ds.RasterXSize+1)*factor,ds.RasterYSize).reshape(ds.RasterYSize, ds.RasterXSize)

    rg_corr = rg_a*rg + rg_b*rg**2 + rg_c*rg**3 + rg_d**rg**4 + rg_f*rg**5 + rg_g*rg**6
    print 'Add range ramp %f r, %f r**2  + %f r**3 + %f r**4 + %f r**5 + %f r**6'%(rg_a, rg_b, rg_c, rg_d, rg_f, rg_g)

    az_corr = az_a*az + az_b*az**2 + az_c*az**3 + az_d**az**4 + az_f*az**5 + az_g*az**6
    print 'Add azimuthal ramp %f az, %f az**2  + %f az**3 + %f az**4 + %f az**5 + %f az**6'%(az_a, az_b, az_c, az_d, az_f, az_g)

    # print rg
    # print 
    # print rg_corr
    # print

    corr_map = rg_corr + az_corr

    flatlos = los_map - corr_map
    flatlos[los_map==0] = 0
    flatlos[np.isnan(los_map)] = np.float('NaN')

    wrapfile =  folder + prerg + str(date1) + '-' + str(date2) + sufrg + '_' + wlook + 'rlks.int'
    ds2 = gdal.Open(wrapfile, gdal.GA_ReadOnly)
    ds2_band = ds2.GetRasterBand(1)
    wrapphi = ds2_band.ReadAsArray(0, 0, ds2.RasterXSize, ds2.RasterYSize, \
        ds2.RasterXSize, ds2.RasterYSize)

    rg = np.tile(np.arange(1,ds2.RasterXSize+1),ds2.RasterYSize).reshape(ds2.RasterYSize, ds2.RasterXSize)
    rg_corr = rg_a*rg + rg_b*rg**2 + rg_c*rg**3 + rg_d**rg**4 + rg_f*rg**5 + rg_g*rg**6
    wrapcorr = (np.cos(rg_corr) + np.sin(rg_corr)*1j)
    
    flatwrap = wrapphi*wrapcorr
    flatwrap[wrapphi==0] = 0
    flatwrap[np.isnan(wrapphi)] = np.float('NaN')

    corwrapfile =  folder + prerg + str(date1) + '-' + str(date2) + sufrg + '_flatrange' + '_' + wlook + 'rlks.int'
    ds3 = gdal.Open(corwrapfile, gdal.GA_ReadOnly)
    ds3_band = ds3.GetRasterBand(1)
    corwrapphi = np.angle(ds3_band.ReadAsArray(0, 0, ds3.RasterXSize, ds3.RasterYSize, \
        ds3.RasterXSize, ds3.RasterYSize))

    nfigure=0
    fig = plt.figure(nfigure,figsize=(11,8))
    vmax,vmin = np.nanpercentile(flatlos,95), np.nanpercentile(flatlos,5)

    ax = fig.add_subplot(2,3,1)
    cax = ax.imshow(np.angle(wrapphi),cmap=cm.gist_rainbow,vmax=np.pi,vmin=-np.pi)
    ax.set_title('Wrapped phase {}'.format(idate))
    setp( ax.get_xticklabels(), visible=None)
    fig.colorbar(cax, orientation='vertical',aspect=10)

    ax = fig.add_subplot(2,3,2)
    cax = ax.imshow(np.angle(wrapcorr),cmap=cm.gist_rainbow,vmax=np.pi,vmin=-np.pi)
    ax.set_title('Wrapped Ramp')
    setp( ax.get_xticklabels(), visible=None)
    fig.colorbar(cax, orientation='vertical',aspect=10)

    ax = fig.add_subplot(2,3,3)
    cax = ax.imshow(np.angle(flatwrap),cmap=cm.gist_rainbow,vmax=np.pi,vmin=-np.pi)
    ax.set_title('Wrapped phase+Ramp')
    setp( ax.get_xticklabels(), visible=None)
    fig.colorbar(cax, orientation='vertical',aspect=10)

    ax = fig.add_subplot(2,3,4)
    cax = ax.imshow(corwrapphi-np.angle(flatwrap),cmap=cm.gist_rainbow,vmax=np.pi,vmin=-np.pi)
    ax.set_title('Check null')
    setp( ax.get_xticklabels(), visible=None)
    fig.colorbar(cax, orientation='vertical',aspect=10)

    ax = fig.add_subplot(2,3,5)
    cax = ax.imshow(np.angle(np.cos(corr_map) + np.sin(corr_map)*1j),cmap=cm.gist_rainbow,vmax=np.pi,vmin=-np.pi)
    ax.set_title('Unwrapped Ramp')
    setp( ax.get_xticklabels(), visible=None)
    fig.colorbar(cax, orientation='vertical',aspect=10)

    ax = fig.add_subplot(2,3,6)
    cax = ax.imshow(flatlos,cmap=cm.gist_rainbow,vmax=vmax,vmin=vmin)
    ax.set_title('Unrapped phase + Ramp')
    setp( ax.get_xticklabels(), visible=None)
    fig.colorbar(cax, orientation='vertical',aspect=10)

    fig.savefig(folder+'add_ramp_corrrection.png', format='PNG',dpi=150)

    if plot=='yes':
        plt.show()

    outfile = folder + prefix + str(date1) + '-' + str(date2) + suffout + '_' + rlook + 'rlks.unw'  
    outrsc = folder + prefix + str(date1) + '-' + str(date2) + suffout + '_' + rlook + 'rlks.unw.rsc' 
    
    drv = gdal.GetDriverByName("roi_pac")
    dst_ds = drv.Create(outfile, ds.RasterXSize, ds.RasterYSize, 2, gdal.GDT_Float32)
    dst_band1 = dst_ds.GetRasterBand(1)
    dst_band2 = dst_ds.GetRasterBand(2)
    dst_band1.WriteArray(coh_map,0,0)
    dst_band2.WriteArray(flatlos,0,0)
    shutil.copy(rscfile,outrsc)

    del dst_ds, ds
    del los_map, coh_map, flatlos, corr_map, rg_corr, az_corr

    plt.close('all')
