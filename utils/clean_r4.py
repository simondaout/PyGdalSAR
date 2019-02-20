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
clean_r4.py
-------------
Clean a r4 file given an other r4 file (mask) and a threshold on this mask

Usage: clean_r4.py --infile=<path> --outfile=<path>  [--mask=<path>] [--threshold=<value>] \
[--perc=<value>] [--clean=<values>] [--buff=<value>] [--lectfile=<path>] [--scale=<value>] 

Options:
-h --help           Show this screen.
--infile PATH       r4 file to clean
--outfile PATH      output file
--mask PATH         r4 file used as mask
--threshold VALUE   threshold value on mask file (Keep pixel with mask > threshold)
--scale VALUE       scale the mask [default:1]
--lectfile PATH     Path of the lect.in file [default: lect.in]
--clean VALUE       Crop option with smoothing of boundaries [default: 0,ncol,0,nlign]
--buff VALUE        Number of pixels for crop smothing  (default: 50)
--perc VALUE        Percentile of hidden LOS pixel for the estimation and clean outliers [default:99.9]
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
if arguments["--mask"] ==  None:
    maskf = 'no'
else:
    maskf = arguments["--mask"]
if arguments["--threshold"] ==  None:
    seuil = -np.inf
else:
    seuil = float(arguments["--threshold"])
if arguments["--lectfile"] ==  None:
   lecfile = "lect.in"
else:
   lecfile = arguments["--lectfile"]
if arguments["--scale"] ==  None:
   scale = 1.
else:
   scale = float(arguments["--scale"])

if arguments["--perc"] ==  None:
    perc = 99.9
else:
    perc = float(arguments["--perc"])

# read lect.in 
ncol, nlign = map(int, open(lecfile).readline().split(None, 2)[0:2])
fid = open(infile, 'r')
m = np.fromfile(fid,dtype=np.float32).reshape((nlign,ncol))
# m[m==0.0] = np.float('NaN')

# clean
maxlos,minlos=np.nanpercentile(m,perc),np.nanpercentile(m,(100-perc))
print 'Clean outliers outside:', maxlos,minlos
kk = np.nonzero(
    np.logical_or(m<minlos,m>maxlos))
m[kk] = np.float('NaN')
mf = np.copy(m)

# mask
if maskf is not 'no':
    fid2 = open(maskf, 'r')
    mask = np.fromfile(fid2,dtype=np.float32).reshape((nlign,ncol))
    mask =  mask*scale
    kk = np.nonzero(np.logical_or(mask==0, mask>seuil)) 
    mf[kk] = np.float('NaN')
else:
    mask = np.zeros((nlign,ncol))

# crop
if arguments["--clean"] ==  None:
    crop = [0,ncol,0,nlign]
else:
    crop = map(float,arguments["--clean"].replace(',',' ').split())
ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])

if iend-ibeg<ncol and jend-jbeg<nlign:
    if arguments["--buff"] ==  None:
        buf = 50
    else:
        buf = int(arguments["--buff"])
    # ibeg, iend = 250, 900
    # jbeg, jend = 2000, 2900 
    ibeg1,iend1 = 0, ibeg
    ibeg2,iend2 = iend,ncol
    jbeg1,jend1 = 0, jbeg
    jbeg2,jend2 = jend, nlign

    # top 
    if jend1 > 0:
        for j in xrange(buf):
          for i in xrange(ibeg,iend):
                mf[jend1+j,i] = mf[jend1+j,i]*(np.float(j+1)/buf)
    
    #bottom
    if jbeg2 < nlign:
        for j in xrange(buf):
          for i in xrange(ibeg,iend):
            # print jbeg2-(buf-j)
            mf[jbeg2-(buf-j),i] = mf[jbeg2-(buf-j),i] - mf[jbeg2-(buf-j),i]*(np.float(j+1)/buf)

    # left
    if iend1 > 0:
        for j in xrange(buf):
          for i in xrange(jbeg,jend):   
            mf[i,iend1+(buf-(j+1))] = mf[i,iend1+(buf-(j+1))] - mf[i,iend1+(buf-(j+1))]*(np.float(j+1)/buf)
        # sys.exit()

    # right
    if ibeg2 < ncol:
        for j in xrange(buf):
          for i in xrange(jbeg,jend):
            mf[i,ibeg2-(buf-j)] = mf[i,ibeg2-(buf-j)] - mf[i,ibeg2-(buf-j)]*(np.float(j+1)/buf)

    mf[jbeg1:jend1,:] = 0.0
    mf[jbeg2:jend2,:] = 0.0
    mf[:,ibeg1:iend1] = 0.0
    mf[:,ibeg2:iend2] = 0.0

# mf[jbeg+buf:jend-buf,ibeg+buf:iend-buf] = np.float('NaN')

#plt.imshow(mf)
#plt.show()
#sys.exit()

#save output
fid3 = open(outfile,'wb')
mf.flatten().astype('float32').tofile(fid3)

# Plot
vmax = np.nanpercentile(m,98) 
# vmax= 2.
vmin = -vmax

fig = plt.figure(0,figsize=(12,8))
ax = fig.add_subplot(1,4,1)
# hax = ax.imshow(mask, cm.Greys, vmin=0, vmax=seuil)
cax = ax.imshow(m, cm.RdBu,vmin=vmin, vmax=vmax)
ax.set_title(infile)
setp( ax.get_xticklabels(), visible=False)
cbar = fig.colorbar(cax, orientation='vertical',aspect=9)

ax = fig.add_subplot(1,4,2)
cax = ax.imshow(mask, cm.RdBu, vmin=0, vmax=seuil)
ax.set_title('MASK')
setp( ax.get_xticklabels(), visible=False)
cbar = fig.colorbar(cax, orientation='vertical',aspect=9)

ax = fig.add_subplot(1,4,3)
cax = ax.imshow(mf[jbeg-buf:jend+buf,ibeg-buf:iend+buf], cm.RdBu, vmin=vmin, vmax=vmax)
setp( ax.get_xticklabels(), visible=False)
setp( ax.get_xticklabels(), visible=False)
cbar = fig.colorbar(cax, orientation='vertical',aspect=9)
ax.set_title(outfile)

ax = fig.add_subplot(1,4,4)
cax = ax.imshow(mf[jbeg-buf:jend+buf,ibeg-buf:iend+buf]-m[jbeg-buf:jend+buf,ibeg-buf:iend+buf], cm.RdBu, vmin=vmin, vmax=vmax)
setp( ax.get_xticklabels(), visible=False)
setp( ax.get_xticklabels(), visible=False)
cbar = fig.colorbar(cax, orientation='vertical',aspect=9)
ax.set_title('Check diff')

# ax = fig.add_subplot(1,4,3)
# cax = ax.imshow(mf, cm.RdBu, vmin=vmin, vmax=vmax)
# setp( ax.get_xticklabels(), visible=False)
# setp( ax.get_xticklabels(), visible=False)
# cbar = fig.colorbar(cax, orientation='vertical',aspect=9)
# ax.set_title(outfile)

# ax = fig.add_subplot(1,4,4)
# cax = ax.imshow(mf-m, cm.RdBu, vmin=vmin, vmax=vmax)
# setp( ax.get_xticklabels(), visible=False)
# setp( ax.get_xticklabels(), visible=False)
# cbar = fig.colorbar(cax, orientation='vertical',aspect=9)
# ax.set_title('Check diff')


fig.tight_layout()
plt.show()
