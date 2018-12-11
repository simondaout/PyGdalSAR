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

Usage: clean_r4.py --infile=<path> --mask=<path> --threshold=<value> --outfile=<path> [--perc=<value>] [--clean=<values>] [--lectfile=<path>] [--scale=<value>] 

Options:
-h --help           Show this screen.
--infile PATH       r4 file to clean
--outfile PATH      output file
--mask PATH         r4 file used as mask
--threshold VALUE   threshold value on mask file (Keep pixel with mask > threshold)
--scale VALUE       scale the mask [default:1]
--lectfile PATH     Path of the lect.in file [default: lect.in]
--clean VALUE       Clean option [default: 0,ncol,0,nlign]
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

if arguments["--clean"] ==  None:
    mask = [0,0,0,0]
else:
    print arguments["--clean"]
    mask = map(float,arguments["--clean"].replace(',',' ').split())
ibeg,iend,jbeg,jend = int(mask[0]),int(mask[1]),int(mask[2]),int(mask[3])

if arguments["--perc"] ==  None:
    perc = 99.9
else:
    perc = float(arguments["--perc"])

# read lect.in 
ncol, nlign = map(int, open(lecfile).readline().split(None, 2)[0:2])
fid = open(infile, 'r')
m = np.fromfile(fid,dtype=np.float32).reshape((nlign,ncol))
# m[m==0.0] = np.float('NaN')

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

# ibeg, iend = 250, 900
# jbeg, jend = 2000, 2900 
mcrop = np.nanpercentile(mf[jbeg:jend,ibeg:iend],85)

ibeg1,iend1 = 0, ibeg
ibeg2,iend2 = iend,ncol
jbeg1,jend1 = 0, jbeg
jbeg2,jend2 = jend, nlign

buf = 50
# top 
# print mcrop
for j in xrange(buf):
  # print jend1+j, (np.float(j+1)/buf)
  for i in xrange(ibeg,iend):
	# mf[jend1+j,:] = mcrop*(np.float(j+1)/buf)
	# print jend1+j, mcrop*(np.float(j+1)/buf)
    mf[jend1+j,i] = mf[jend1+j,i]*(np.float(j+1)/buf)

#bottom
# print mcrop
for j in xrange(buf):
  for i in xrange(ibeg,iend):
	# print jbeg2-(buf-j), mcrop - mcrop*(np.float(j+1)/buf)
	# mf[jbeg2-(buf-j),:] = mcrop - mcrop*(np.float(j+1)/buf)
	mf[jbeg2-(buf-j),i] = mf[jbeg2-(buf-j),i] - mf[jbeg2-(buf-j),i]*(np.float(j+1)/buf)

# left
# print mcrop
for j in xrange(buf):
  for i in xrange(ibeg,iend):	
	# print iend1+(buf-j), mcrop - mcrop*(np.float(j+1)/buf)
	# mf[jbeg+buf:jend-buf,iend1+(buf-j)] = mcrop - mcrop*(np.float(j+1)/buf)
	mf[i,iend1+(buf-j)] = mf[i,iend1+(buf-j)] - mf[i,iend1+(buf-j)]*(np.float(j+1)/buf)

# right
# print mcrop
for j in xrange(buf):
  for i in xrange(ibeg,iend):
	# print ibeg2-(buf-j),  mcrop - mcrop*(np.float(j+1)/buf)
	# mf[jbeg+buf:jend-buf,ibeg2-(buf-j)] = mcrop - mcrop*(np.float(j+1)/buf)
	mf[i,ibeg2-(buf-j)] = mf[i,ibeg2-(buf-j)] - mf[i,ibeg2-(buf-j)]*(np.float(j+1)/buf)


mf[jbeg1:jend1,:] = 0.0
mf[jbeg2:jend2,:] = 0.0
mf[:,ibeg1:iend1] = 0.0
mf[:,ibeg2:iend2] = 0.0
mf[jbeg1:jend1,ibeg1:iend1] = 0.0
mf[jbeg2:jend2,ibeg2:iend2] = 0.0



#plt.imshow(mf)
#plt.show()
#sys.exit()

#save output
fid3 = open(outfile,'wb')
mf.flatten().astype('float32').tofile(fid3)

# Plot
vmax = np.nanpercentile(m,99) 
vmin = -vmax

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
