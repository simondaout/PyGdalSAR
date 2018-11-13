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
extend_r4.py
-------------
Cut r4 file 

Usage: extend_r4.py --infile=<path> --outfile=<path> [<ibeg>] [<iend>] [<jbeg>] [<jend>] [--lectfile=<path>]

Options:
-h --help           Show this screen.
--infile PATH       File to be extend
--outfile PATH      Cutting file
--ibeg VALUE        Ligne numbers bounded the cutting zone [default: 0]
--iend VALUE        Ligne numbers bounded the cutting zone [default: nlign]
--jbeg VALUE        Column numbers bounded the cutting zone [default: 0]
--jend VALUE        Column numbers bounded the cutting zone [default: ncol]
--lectfile PATH     Path of the lect.in file [default: lect.in]
"""

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from pylab import *

# docopt (command line parser)
import docopt

# read arguments
arguments = docopt.docopt(__doc__)
infile = arguments["--infile"]
outfile = arguments["--outfile"]

if arguments["--lectfile"] ==  None:
    lecfile = "lect.in"
else:
    lecfile = arguments["--lectfile"]

# read lect.in 
ncol, nlign = map(int, open(lecfile).readline().split(None, 2)[0:2])

if arguments["<ibeg>"] ==  None:
    ibeg = 0
else:
    ibeg = int(arguments["<ibeg>"])
if arguments["<iend>"] ==  None:
    iend = nlign
else:
    iend = int(arguments["<iend>"])
if arguments["<jbeg>"] ==  None:
    jbeg = 0
else:
    jbeg = int(arguments["<jbeg>"])
if arguments["<jend>"] ==  None:
    jend = ncol
else:
    jend = int(arguments["<iend>"])

# Program
fid = open(infile, 'r')
m = np.fromfile(fid,dtype=np.float32)
map = m.reshape((iend-ibeg,jbeg-jend))
extendmap = np.ones((nlign,ncol))*float('NaN')
extendmap[ibeg:iend,jbeg:jend] = map

# save output
fid2 = open(outfile,'wb')
extendmap.flatten().astype('float32').tofile(fid2)

vmax = np.nanmean(m) + 2*np.nanstd(m)
vmin = np.nanmean(m) - 2*np.nanstd(m)

# Plot
fig = plt.figure(0,figsize=(8,9))
ax = fig.add_subplot(1,2,1)
cax = ax.imshow(map,cmap=cm.jet,vmax=vmax,vmin=vmin, interpolation='bicubic')
ax.set_title(infile)
setp( ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(1,2,2)
cax = ax.imshow(extendmap,cmap=cm.jet,vmax=vmax,vmin=vmin, interpolation='bicubic')
ax.set_title(outfile)
setp( ax.get_xticklabels(), visible=False)

cbar = fig.colorbar(cax, orientation='vertical',aspect=5)

plt.show()
