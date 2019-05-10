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
[--perc=<value>] [--clean=<values>] [--buff=<value>] [--lectfile=<path>] [--scale=<value>] \
[--ramp=<yes/no>] [--ref=<jstart,jend>]

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
--ramp<yes/no>      If yes estimate a quadratic ramp in range
--ref=<jstart,jend> Set to zero displacements from jstart to jend
"""

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided
import scipy.optimize as opt
import scipy.linalg as lst

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

if arguments["--ramp"] == None:
    ramp = 'no'
else:
    ramp = arguments["--ramp"]

if arguments["--ref"] == None:
    refstart, refend = None,None
else:
    ref = map(int,arguments["--ref"].replace(',',' ').split())
    refstart,refend = ref[0], ref[1]

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

if ramp=='yes':
    index = np.nonzero(~np.isnan(mf))
    temp = np.array(index).T
    mi = m[index].flatten()
    az = temp[:,0]; rg = temp[:,1]

    G=np.zeros((len(mi),4))
    G[:,0] = rg**2
    G[:,1] = rg
    G[:,2] = az
    G[:,3] = 1

    x0 = lst.lstsq(G,mi)[0]
    _func = lambda x: np.sum(((np.dot(G,x)-mi))**2)
    _fprime = lambda x: 2*np.dot(G.T, (np.dot(G,x)-mi))
    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0)[0]    

    pars = np.dot(np.dot(np.linalg.inv(np.dot(G.T,G)),G.T),mi)
    a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
    print 'Remove ramp %f x**2 %f x  + %f y + %f'%(a,b,c,d)

    G=np.zeros((len(mask.flatten()),4))
    for i in xrange(nlign):
        G[i*ncol:(i+1)*ncol,0] = np.arange((ncol))**2
        G[i*ncol:(i+1)*ncol,1] = np.arange((ncol))
        G[i*ncol:(i+1)*ncol,2] = i
    G[:,3] = 1
    temp = (mf.flatten() - np.dot(G,pars))
    mf=temp.reshape(nlign,ncol)

if (refstart is not None) and (refend is not None):
    cst = np.nanmean(mf[refstart:refend,:])
    mf = mf - cst

# crop
if arguments["--clean"] ==  None:
    crop = [0,ncol,0,nlign]
else:
    crop = map(float,arguments["--clean"].replace(',',' ').split())
ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])
if iend>ncol:
    iend=ncol
if jend>nlign:
    jend=nlign

if iend-ibeg<ncol or jend-jbeg<nlign:
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

    mf[jbeg1:jend1,:] = np.float('NaN')
    mf[jbeg2:jend2,:] = np.float('NaN')
    mf[:,ibeg1:iend1] = np.float('NaN')
    mf[:,ibeg2:iend2] = np.float('NaN')
#plt.imshow(mf)
#plt.show()
#sys.exit()

#save output
fid3 = open(outfile,'wb')
mf.flatten().astype('float32').tofile(fid3)

# Plot
vmax = np.nanpercentile(mf,98) 
# vmax= 2.
vmin = -vmax

fig = plt.figure(0,figsize=(12,8))
ax = fig.add_subplot(1,4,1)
# hax = ax.imshow(mask, cm.Greys, vmin=0, vmax=seuil)
cax = ax.imshow(mf, cm.RdBu,vmin=vmin, vmax=vmax)
ax.set_title(infile)
setp( ax.get_xticklabels(), visible=False)
cbar = fig.colorbar(cax, orientation='vertical',aspect=9)

ax = fig.add_subplot(1,4,2)
cax = ax.imshow(mask, cm.RdBu, vmin=0, vmax=seuil)
ax.set_title('MASK')
setp( ax.get_xticklabels(), visible=False)
cbar = fig.colorbar(cax, orientation='vertical',aspect=9)

if iend-ibeg<ncol and jend-jbeg<nlign:
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
else:
    ax = fig.add_subplot(1,4,3)
    cax = ax.imshow(mf, cm.RdBu, vmin=vmin, vmax=vmax)
    setp( ax.get_xticklabels(), visible=False)
    setp( ax.get_xticklabels(), visible=False)
    cbar = fig.colorbar(cax, orientation='vertical',aspect=9)
    ax.set_title(outfile)


fig.tight_layout()
plt.show()
