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
cut_r4.py
-------------
Cut r4 file 

Usage: cut_r4.py --infile=<path> --outfile=<path> [--lectfile=<path>] [<ibeg>] [<iend>] [<jbeg>] [<jend>]

Options:
-h --help           Show this screen.
--infile PATH       File to be cut
--outfile PATH      Cutting file
--lectfile PATH     Path of the lect.in file [default: lect.in]
--ibeg VALUE        Ligne numbers bounded the cutting zone [default: 0]
--iend VALUE        Ligne numbers bounded the cutting zone [default: nlign]
--jbeg VALUE        Column numbers bounded the cutting zone [default: 0]
--jend VALUE        Column numbers bounded the cutting zone [default: ncol]
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
    jend = int(arguments["<jend>"])

# Program
fid = open(infile, 'r')
m = np.fromfile(fid,dtype=np.float32)   
print len(m)/ncol
map = m[:nlign*ncol].reshape((nlign,ncol)) 
print ibeg,iend,jbeg,jend
cutmap = as_strided(map[ibeg:iend,jbeg:jend])

# save output
fid2 = open(outfile,'wb')
cutmap.flatten().astype('float32').tofile(fid2)

# 
kk = flatnonzero(m==0)
m[kk] = float('NaN')

#print np.mean(m), np.percentile(abs(m),98)
#vmax=1
#vmin=0
vmax = np.nanmean(cutmap) + 2*np.nanstd(cutmap)
vmin = np.nanmean(cutmap) - 2*np.nanstd(cutmap)

# Plot
fig = plt.figure(0,figsize=(8,9))
ax = fig.add_subplot(1,2,1)
cax = ax.imshow(map,cmap=cm.jet,vmax=vmax,vmin=vmin, interpolation='bicubic')
ax.set_title(infile)
setp( ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(1,2,2)
cax = ax.imshow(cutmap,cmap=cm.jet,vmax=vmax,vmin=vmin, interpolation='bicubic')
ax.set_title(outfile)
setp( ax.get_xticklabels(), visible=False)

cbar = fig.colorbar(cax, orientation='vertical',aspect=5)
fig.savefig('{}.eps'.format(outfile), format='EPS',dpi=150)

plt.show()
