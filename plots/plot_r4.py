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
plot_r4.py
-------------

Usage: 	plot_r4.py --infile=<path>
		plot_r4.py --infile=<path> [--lectfile=<path>] [--vmin=<value>] [--vmax=<value>] \
[--wrap=<yes/no>] [--rad2mm=<value>] [--crop=<values>] [--title=<value>]

Options:
-h --help             Show this screen.
--infile=<file>       File to be cut
--lectfile=<file>     Path of the lect.in file [default: lect.in]
--vmax=<value> 	      Max colorscale (default: 90th percentile)
--vmin=<value>        Min colorscale (default: 10th percentile)
--wrap=<value> 	      Wrapped phase [default: no]
--rad2mm=<value>      Convert data [default: 1]
--tile=<value>        Title plot 
--crop=<value>        Define a region of interest [default: 0,nlign,0,ncol]
"""

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

# docopt (command line parser)
import docopt

# read arguments
arguments = docopt.docopt(__doc__)

infile = arguments["--infile"]
if arguments["--lectfile"] ==  None:
    lecfile = "lect.in"
else:
    lecfile = arguments["--lectfile"]

# read lect.in 
ncol, nlign = map(int, open(lecfile).readline().split(None, 2)[0:2])

if arguments["--crop"] ==  None:
    crop = [0,nlign,0,ncol]
else:
    crop = map(float,arguments["--crop"].replace(',',' ').split())
ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])

# Program
fid = open(infile, 'r')
m = np.fromfile(fid,dtype=np.float32)[:nlign*ncol].reshape((nlign,ncol))[ibeg:iend,jbeg:jend]
kk = np.nonzero(np.logical_or(np.logical_or(~np.isnan(m), np.abs(m)<999.),m!=0.0))
mprim = m[kk]

if arguments["--vmax"] ==  None:
	vmax = np.nanpercentile(mprim, 99)
else:
	vmax = np.float(arguments["--vmax"])

if arguments["--vmin"] ==  None:
	vmin = np.nanpercentile(mprim, 1)
else:
	vmin = np.float(arguments["--vmin"])
# print vmax,vmin

if arguments["--wrap"] !=  None:	
	m = np.mod(m,2*np.pi)-np.pi 
	vmax=np.pi
	vmin=-np.pi

if arguments["--rad2mm"] ==  None:
        rad2mm = 1
else:
        rad2mm = np.float(arguments["--rad2mm"])

# Plot
m, vmin, vmax = m*rad2mm, vmin*rad2mm, vmax*rad2mm

fig = plt.figure(0,figsize=(8,10))
ax = fig.add_subplot(1,1,1)
cax = ax.imshow(m,cmap=cm.jet,vmax=vmax,vmin=vmin)
if arguments["--title"] ==  None:
	ax.set_title(infile)
else:
	ax.set_title(arguments["--title"])
setp( ax.get_xticklabels(), visible=False)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

fig.savefig('{}.eps'.format(infile), format='EPS',dpi=150)
plt.show()
