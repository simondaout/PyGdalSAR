#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Simon Daout 
############################################

"""\
plot\_RMS_ts.py
-------------
Create jpg for all RMS files output from time series analysis

Usage: plot\_RMS_ts.py [--lectfile=<path>] [--outputdir=<path>]  [--vmin=<value>] [--vmax=<value>] [--dates_list=<path>] [--int_list=<path>] 

Options:
-h --help             Show this screen.
--lectfile=<file>     Path of the lect.in file for r4 format [default:lect.in]
--int_list=<file>     Text file containing list of interferograms dates in two colums, $data1 $date2 [default: interf_pair.rsc]
--dates_list=<file>   Path to text file containing date,bp,bt,doppler_frq,date_dec [default: baseline.rsc]
  --outputdir=<dir>   Output directory [default: ./plots]
--vmax                Max colorscale [default: 98th percentile]
--vmin                Min colorscale [default: 2th percentile]
"""

import os, sys

# docopt (command line parser)
try:
    from nsbas import docopt
except:
    import docopt

import numpy as np
from numpy.lib.stride_tricks import as_strided

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import setp
from osgeo import gdal
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

# read arguments
arguments = docopt.docopt(__doc__)
if arguments["--lectfile"] ==  None:
    lecfile = "lect.in"
else:
    lecfile = arguments["--lectfile"]
if arguments["--int_list"] == None:
    int_list = 'list_pair'
else:
    int_list=arguments["--int_list"]
if arguments["--dates_list"] == None:
    baseline = 'list_dates'
else:
    baseline=arguments["--dates_list"]

# A more predictable makedirs
def makedirs(name):
    if os.path.exists(name):
        return
    os.makedirs(name)

if arguments["--outputdir"] == None:
    plotdir = './plots/'
else:
    plotdir = os.path.join(arguments["--outputdir"]) + '/'
# Create output directories
makedirs(plotdir)

try:
    from matplotlib.colors import LinearSegmentedColormap
    cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
    cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
    cmap = cmap.reversed()
except:
    cmap=cm.jet

# read int
date_1,date_2=np.loadtxt(int_list,comments="#",unpack=True,dtype='i,i')
kmax=len(date_1)
print("number of interferogram: ",kmax)

# read dates
im,imd,bt,bp =np.loadtxt(baseline,comments="#",usecols=(0,1,2,3),unpack=True,dtype='i,f,f,f')
print("image list=",baseline)
nmax=len(imd)
print("number of image: ",nmax)

def plot_raster(infile,n):
    fig = plt.figure(n,figsize=(12,8))
    print('Plot', infile)
    fid = open(infile, 'r')
    ncols, nlines = list(map(int, open(lecfile).readline().split(None, 2)[0:2]))
    phi = np.fromfile(fid,dtype=np.float32)[:nlines*ncols].reshape((nlines,ncols))
    ax = fig.add_subplot(1,1,1)
    masked_array = np.ma.array(phi, mask=np.isnan(phi))

    if (arguments["--vmax"] is not None) and (arguments["--vmin"] is not None):
        if arguments["--vmax"] is not  None:
          vmax = float(arguments["--vmax"])
          if arguments["--vmin"] is not  None:
            vmin = float(arguments["--vmin"])
          else:
            vmin = -vmax
    else:
        vmax = np.nanpercentile(phi,98)
        vmin = np.nanpercentile(phi,2)

    cax = ax.imshow(masked_array, cmap, interpolation='nearest',vmax=vmax,vmin=vmin)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)
    fig.tight_layout()
    fig.savefig(plotdir+'/{}.png'.format(infile), format='PNG',dpi=90)
    plt.close('all')
    del fig 

n = 0
for k in range(kmax):
    infile = 'RMSpixel_{}_{}'.format(date_1[k],date_2[k])
    plot_raster(infile,n)
    n = n + 1

for k in range(nmax):
    infile = 'RMSpixel_{}'.format(im[k])
    plot_raster(infile,n)
    n = n + 1

# plot RMSpixel
plot_raster('RMSpixel',n)

