#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT 
################################################################################

"""\
clean_r4.py
-------------
Clean a r4 file given an other r4 file (mask) and a threshold on this mask

Usage: clean_r4.py --infile=<path> --mask=<path> --threshold=<value> --outfile=<path> [--lectfile=<path>] [--scale=<value>] 

Options:
-h --help           Show this screen.
--infile PATH       r4 file to clean
--outfile PATH      output file
--mask PATH         r4 file used as mask
--threshold VALUE   threshold value on mask file (Keep pixel with mask > threshold)
--scale VALUE       scale the mask [default:1]
--lectfile PATH     Path of the lect.in file [default: lect.in]
"""

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided


import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from pylab import *

import docopt
arguments = docopt.docopt(__doc__)
infile = arguments["--infile"]
outfile = arguments["--outfile"]
maskf = arguments["--mask"]
seuil = float(arguments["--threshold"])
if arguments["--lectfile"] ==  None:
   lecfile = "lect.in"
else:
   lecfile = arguments["--lectfile"]
if arguments["--scale"] ==  None:
   scale = 1.
else:
   scale = float(arguments["--scale"])

# read lect.in 
ncol, nlign = map(int, open(lecfile).readline().split(None, 2)[0:2])
fid = open(infile, 'r')
m = np.fromfile(fid,dtype=np.float32).reshape((nlign,ncol))
# m[m==0.0] = np.float('NaN')

perc=99.9
maxlos,minlos=np.nanpercentile(m,perc),np.nanpercentile(m,(100-perc))
print maxlos,minlos
kk = np.nonzero(
    np.logical_or(m<minlos,m>maxlos))
m[kk] = np.float('NaN')

# mask
fid2 = open(maskf, 'r')
mask = np.fromfile(fid2,dtype=np.float32).reshape((nlign,ncol))
mask =  mask*scale
kk = np.nonzero(np.logical_or(mask==0, mask>seuil))
mf = np.copy(m)
mf[kk] = np.float('NaN')
#plt.imshow(mf)
#plt.show()
#sys.exit()

#save output
fid3 = open(outfile,'wb')
mf.flatten().astype('float32').tofile(fid3)

# Plot
vmax = np.nanmean(mf) + 2*np.nanstd(mf)
vmin = np.nanmean(mf) - 2*np.nanstd(mf)

fig = plt.figure(0,figsize=(12,8))
ax = fig.add_subplot(1,3,1)
# hax = ax.imshow(mask, cm.Greys, vmin=0, vmax=seuil)
cax = ax.imshow(m, cm.jet,vmin=vmin, vmax=vmax)
ax.set_title(infile)
setp( ax.get_xticklabels(), visible=False)
cbar = fig.colorbar(cax, orientation='vertical',aspect=9)

ax = fig.add_subplot(1,3,2)
cax = ax.imshow(mask, cm.jet, vmin=0, vmax=seuil)
ax.set_title(maskf)
setp( ax.get_xticklabels(), visible=False)
cbar = fig.colorbar(cax, orientation='vertical',aspect=9)

ax = fig.add_subplot(1,3,3)
cax = ax.imshow(mf, cm.jet, vmin=vmin, vmax=vmax)
setp( ax.get_xticklabels(), visible=False)
setp( ax.get_xticklabels(), visible=False)
cbar = fig.colorbar(cax, orientation='vertical',aspect=9)
ax.set_title(outfile)

fig.tight_layout()
plt.show()
