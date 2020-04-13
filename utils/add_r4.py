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
add_r4.py
-------------
Add, substract, compute the amplitude (infile1**2,infile2**2) ar the phase (arctan(infile1/infile2)) from two real4 files.

Usage: add_r4.py --infile1=<path> --infile2=<path> [--outfile=<path>] [--lectfile=<path>] --sign=<+/-/amp/dphi> [--rad2mm=<value>] [--vmin=<value>] [--vmax=<value>] 
add_r4.py -h | --help

Options:
-h --help           Show this screen.
--infile1 PATH      
--infile2 PATH
--outfile PATH
--lectfile PATH     Path of the lect.in file [default: lect.in]
--sign VALUE        amp=sqrt(infile1**2,infile2**2), dphi=arctan(infile1/infile2)
--rad2mm=<value>      Convert data [default: 1]
--vmax                Max colorscale [default: 98th percentile]
--vmin                Min colorscale [default: 2th percentile]
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pylab import *
# docopt (command line parser)
import docopt

# read arguments
arguments = docopt.docopt(__doc__)
infile1 = arguments["--infile1"]
infile2 = arguments["--infile2"]
sign = arguments["--sign"]

if arguments["--lectfile"] ==  None:
    lecfile = "lect.in"
else:
    lecfile = arguments["--lectfile"]

# read lect.in
ncol, nlign = map(int, open(lecfile).readline().split(None, 2)[0:2])

if arguments["--outfile"] is not  None:
    outfile = arguments["--outfile"]

if arguments["--rad2mm"] ==  None:
        rad2mm = 1
else:
        rad2mm = np.float(arguments["--rad2mm"])

# read r4
in1 = np.fromfile(infile1,dtype=np.float32) 
in2 = np.fromfile(infile2,dtype=np.float32)
in1[in1==0] = np.float('NaN')
in2[in2==0] = np.float('NaN')

# center to 0
#in1 = in1 - np.nanmean(in1)
#in2 = in2 - np.nanmean(in2)

# reshape
in11 = in1.reshape((nlign,ncol)) 
in22 = in2.reshape((nlign,ncol))
#out= in11+in22
if sign=='+':
    out= in11 + in22
if sign=='-':
    out= in11 - in22
if sign=='amp':
    out = np.sqrt(in11**2+in22**2)
if sign=='dphi':
    out = np.arctan2(in11,in22)

# save outfile in the convention ROIPAC
if arguments["--outfile"] is not  None:
  out.flatten().astype('float32').tofile(outfile)

# convert
in11, in22, out = in11*rad2mm, in22*rad2mm, out*rad2mm

# plot maps
if arguments["--vmax"] is not  None:
    vmax = np.float(arguments["--vmax"])
else:
    vmax = np.max( [np.nanpercentile(in11,95),np.abs(np.nanpercentile(in11,5))] )
if arguments["--vmin"] is not  None:
    vmin = np.float(arguments["--vmin"])
else:
    vmin = -vmax

try:
        from matplotlib.colors import LinearSegmentedColormap
        cm_locs = '/home/comethome/jdd/ScientificColourMaps5/by_platform/python/'
        cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
        cmap = cmap.reversed()
except:
        cmap=cm.jet

fig = plt.figure(10, figsize=(10,8))
fig.subplots_adjust(hspace=0.5)
ax1 = fig.add_subplot(1,3,1)
cax = ax1.imshow(in11, cmap = cmap, vmax=vmax, vmin=vmin, extent=None)
divider = make_axes_locatable(ax1)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
setp( ax1.get_xticklabels(), visible=False)
ax1.set_title(infile1)

ax2 = fig.add_subplot(1,3,2)
cax = ax2.imshow(in22, cmap = cmap, vmax=vmax, vmin=vmin, extent=None)
setp( ax2.get_xticklabels(), visible=False)
ax2.set_title(infile2)
divider = make_axes_locatable(ax2)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

#vmax = np.max( [np.nanpercentile(out,95),np.abs(np.nanpercentile(out,5))] )

ax3 = fig.add_subplot(1,3,3)
cax = ax3.imshow(out, cmap = cmap, vmax=vmax, vmin=vmin, extent=None)
setp( ax3.get_xticklabels(), visible=False)
ax3.set_title(outfile)
divider = make_axes_locatable(ax3)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
fig.tight_layout()

if arguments["--outfile"] is not  None:
  fig.savefig('{}.pdf'.format('outfile'), format='PDF',dpi=150)

plt.show()
