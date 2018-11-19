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

Usage: add_r4.py --infile1=<path> --infile2=<path> --outfile=<path> --nlign=<value> --ncol=<value> --sign=<+/-/amp/dphi> 
add_r4.py -h | --help

Options:
-h --help           Show this screen.
--infile1 PATH      
--infile2 PATH
--outfile PATH
--nlign VALUE       
--ncol VALUE        
--sign VALUE        amp=sqrt(infile1**2,infile2**2), dphi=arctan(infile1/infile2)
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from pylab import *
# docopt (command line parser)
import docopt

# read arguments
arguments = docopt.docopt(__doc__)
ncol = int(arguments["--ncol"])
nlign = int(arguments["--nlign"])
infile1 = arguments["--infile1"]
infile2 = arguments["--infile2"]
outfile = arguments["--outfile"]
sign = arguments["--sign"]

# read r4
in1 = np.fromfile(infile1,dtype=np.float32) 
in2 = np.fromfile(infile2,dtype=np.float32)

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

# save outfile in the convention ROIPACK
out2 = out
out2.flatten().astype('float32').tofile(outfile)

# remove NaN to compute colorbar
index = np.nonzero(np.isnan(out))
#out[index] = 0

# plot maps
vmax = np.nanmean(in11) +2* np.nanstd(in11)
vmin = np.nanmean(in11) - 2*np.nanstd(in11)
fig = plt.figure(10, figsize=(10,8))
fig.subplots_adjust(hspace=0.5)
ax1 = fig.add_subplot(1,3,1)
ax1.imshow(in11, cmap = cm.jet, vmax=vmax, vmin=vmin, extent=None)
setp( ax1.get_xticklabels(), visible=False)
ax1.set_title(infile1)

ax2 = fig.add_subplot(1,3,2)
ax2.imshow(in22, cmap = cm.jet, vmax=vmax, vmin=vmin, extent=None)
setp( ax2.get_xticklabels(), visible=False)
ax2.set_title(infile2)
fig = plt.figure(10, figsize=(10,8))
fig.subplots_adjust(hspace=0.5)
ax1 = fig.add_subplot(1,3,1)
cax =ax1.imshow(in11, cmap = cm.jet, vmax=vmax, vmin=vmin, extent=None)
setp( ax1.get_xticklabels(), visible=False)
ax1.set_title(infile1)
cbar = fig.colorbar(cax, shrink =0.2 , aspect = 10)

ax2 = fig.add_subplot(1,3,2)
cax = ax2.imshow(in22, cmap = cm.jet, vmax=vmax, vmin=vmin, extent=None)
setp( ax2.get_xticklabels(), visible=False)
ax2.set_title(infile2)
cbar = fig.colorbar(cax, shrink = 0.2, aspect = 10)

vmax = np.nanmean(out) + np.nanstd(out)
vmin = np.nanmean(out) - np.nanstd(out)
ax3 = fig.add_subplot(1,3,3)
cax =ax3.imshow(out, cmap = cm.jet, vmax=vmax, vmin=0, extent=None)
setp( ax3.get_xticklabels(), visible=False)
ax3.set_title(outfile)
cbar = fig.colorbar(cax, shrink = 0.2, aspect=10)
fig.tight_layout()

plt.show()
