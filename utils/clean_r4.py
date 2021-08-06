#!/usr/bin/env python3
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
[--perc=<value>] [--crop=<values>] [--buff=<value>] [--lectfile=<path>] [--scale=<value>] [--scale_mask=<value>]\
[--ramp=<lin/quad/cub/no>] [--ref=<jstart,jend,istart,iend>] [--removeNAN=<yes/no>] [--cst=<value>] [--reverse=<yes/no>]

Options:
-h --help           Show this screen.
--infile PATH       r4 file to clean
--outfile PATH      output file
--mask PATH         r4 file used as mask
--threshold VALUE   threshold value on mask file (Keep pixel with mask > threshold)
--reverse yes/no    if yes mask values above mask threshold (default: no)
--scale VALUE       scale data [default:1]
--scale_mask VALUE  scale the mask [default:1]
--lectfile PATH     Path of the lect.in file [default: lect.in]
--crop VALUE       Crop option with smoothing of boundaries [default: 0,ncol,0,nlign]
--buff VALUE        Number of pixels for crop smothing  (default: 50)
--perc VALUE        Percentile of hidden LOS pixel for the estimation and clean outliers [default:99.9]
--ramp<lin/quad/cub/no>      Correct the map from ramp in range and azimuth
--ref=<jstart,jend,istart,iend> Set to zero displacements from jstart to jend
--removeNAN         replace NaN by 0
--cst               Add constante to map
"""

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided
import scipy.optimize as opt
import scipy.linalg as lst

import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
if arguments["--scale_mask"] ==  None:
   scale_mask = 1.
else:
   scale_mask = float(arguments["--scale_mask"])

if arguments["--removeNAN"] ==  None:
   setzero = 'no'
else:
   setzero = arguments["--removeNAN"]

if arguments["--perc"] ==  None:
    perc = 99.9
else:
    perc = float(arguments["--perc"])

if arguments["--cst"] ==  None:
    shift = 0
else:
    shift = float(arguments["--cst"])

if (arguments["--ramp"] == None) or (arguments["--ramp"] == 'no'):
    ramp = 'no'
elif arguments["--ramp"] == 'lin':
    ramp = 'lin'
elif arguments["--ramp"] == 'quad':
    ramp = 'quad'
elif arguments["--ramp"] == 'cub':
    ramp = 'cub'
elif arguments["--ramp"] == 'yes':
    print('ramp==yes is depricated. Use lin/quad/cub/no. Set ramp to lin')
    ramp = 'lin'
else:
    print('ramp argument not recognized. Exit!')
    sys.exit()
if arguments["--reverse"] == 'yes':
   reverse = True
else:
   reverse = False

# read lect.in 
ncol, nlign = list(map(int, open(lecfile).readline().split(None, 2)[0:2]))

if arguments["--ref"] == None:
    lin_start, lin_jend, col_start, col_jend = None,None,None,None
else:
    ref = list(map(int,arguments["--ref"].replace(',',' ').split()))
    try:
        lin_start,lin_end, col_start, col_end = ref[0], ref[1], ref[2], ref[3]
    except:
        lin_start,lin_end = ref[0], ref[1]
        col_start, col_end = 0, ncol

fid = open(infile, 'r')
m = np.fromfile(fid,dtype=float32).reshape((nlign,ncol))
m[m==0.0] = float('NaN')

mf = np.copy(m)

# clean for ramp
maxlos,minlos=np.nanpercentile(m,90),np.nanpercentile(m,10)
print('Clean outliers for ramp estimation outside:', maxlos,minlos)
kk = np.nonzero(
    np.logical_or(m<minlos,m>maxlos))
m_ramp = np.copy(m)
m_ramp[kk] = float('NaN')

# mask
if maskf != 'no':
    fid2 = open(maskf, 'r')
    mask = np.fromfile(fid2,dtype=float32)[:nlign*ncol].reshape((nlign,ncol))
    mask =  mask*scale_mask
    if reverse:
      kk = np.nonzero(np.logical_or(mask==0, mask<seuil)) 
    else:
      kk = np.nonzero(np.logical_or(mask==0, mask>seuil)) 
    mf[kk] = float('NaN')
else:
    mask = np.zeros((nlign,ncol))

if ramp=='lin':
    index = np.nonzero(~np.isnan(m_ramp))
    temp = np.array(index).T
    mi = m[index].flatten()
    az = temp[:,0]; rg = temp[:,1]

    G=np.zeros((len(mi),3))
    G[:,0] = rg
    G[:,1] = az
    G[:,2] = 1

    x0 = lst.lstsq(G,mi)[0]
    _func = lambda x: np.sum(((np.dot(G,x)-mi))**2)
    _fprime = lambda x: 2*np.dot(G.T, (np.dot(G,x)-mi))
    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0)[0]    

    #pars = np.dot(np.dot(np.linalg.inv(np.dot(G.T,G)),G.T),mi)
    a = pars[0]; b = pars[1]; c = pars[2]
    print('Remove ramp %f x  + %f y + %f'%(a,b,c))

    G=np.zeros((len(mask.flatten()),3))
    for i in range(nlign):
        G[i*ncol:(i+1)*ncol,0] = np.arange((ncol))
        G[i*ncol:(i+1)*ncol,1] = i
    G[:,2] = 1
    temp = (mf.flatten() - np.dot(G,pars))
    mf=temp.reshape(nlign,ncol)
    mf_ramp=temp.reshape(nlign,ncol)

if ramp=='quad':
    index = np.nonzero(~np.isnan(m_ramp))
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
    #pars = np.dot(np.dot(np.linalg.inv(np.dot(G.T,G)),G.T),mi)
    a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3] 
    print('Remove ramp %f x**2 %f x  + %f y + %f'%(a,b,c,d))

    G=np.zeros((len(mask.flatten()),4))
    for i in range(nlign):
        G[i*ncol:(i+1)*ncol,0] = np.arange((ncol))**2
        G[i*ncol:(i+1)*ncol,1] = np.arange((ncol))
        G[i*ncol:(i+1)*ncol,2] = i
    G[:,3] = 1
    temp = (mf.flatten() - np.dot(G,pars))
    mf=temp.reshape(nlign,ncol)
    mf_ramp=temp.reshape(nlign,ncol)

elif ramp=='cub':
    index = np.nonzero(~np.isnan(m_ramp))
    temp = np.array(index).T
    mi = m[index].flatten()
    az = temp[:,0]; rg = temp[:,1]

    G=np.zeros((len(mi),6))
    G[:,0] = rg**2
    G[:,1] = rg
    G[:,2] = az**3
    G[:,3] = az**2
    G[:,4] = az
    G[:,5] = 1

    x0 = lst.lstsq(G,mi)[0]
    _func = lambda x: np.sum(((np.dot(G,x)-mi))**2)
    _fprime = lambda x: 2*np.dot(G.T, (np.dot(G,x)-mi))
    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0)[0]    
    #pars = np.dot(np.dot(np.linalg.inv(np.dot(G.T,G)),G.T),mi)
    a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
    print('Remove ramp %f x**2 + %f x  + %f y**3 + %f y**2 + %f y + %f'%(a,b,c,d,e,f))

    G=np.zeros((len(mask.flatten()),6))
    for i in range(nlign):
        G[i*ncol:(i+1)*ncol,0] = np.arange((ncol))**2
        G[i*ncol:(i+1)*ncol,1] = np.arange((ncol))
        G[i*ncol:(i+1)*ncol,2] = i**3
        G[i*ncol:(i+1)*ncol,3] = i**2
        G[i*ncol:(i+1)*ncol,4] = i
    G[:,5] = 1
    temp = (mf.flatten() - np.dot(G,pars))
    mf=temp.reshape(nlign,ncol)
    mf_ramp=temp.reshape(nlign,ncol)

else:
    mf_ramp= np.copy(mf)

if (arguments["--ref"] != None) :
    cst = np.nanmean(mf[lin_start:lin_end,col_start:col_end])
    mf = mf - cst
    mf_ramp = mf_ramp - cst

# crop
if arguments["--crop"] ==  None:
    crop = [0,ncol,0,nlign]
else:
    crop = list(map(float,arguments["--crop"].replace(',',' ').split()))
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
        for j in range(buf):
          for i in range(ibeg,iend):
                mf[jend1+j,i] = mf[jend1+j,i]*(float(j+1)/buf)
    
    #bottom
    if jbeg2 < nlign:
        for j in range(buf):
          for i in range(ibeg,iend):
            # print(jbeg2-(buf-j))
            mf[jbeg2-(buf-j),i] = mf[jbeg2-(buf-j),i] - mf[jbeg2-(buf-j),i]*(float(j+1)/buf)

    # left
    if iend1 > 0:
        for j in range(buf):
          for i in range(jbeg,jend):   
            mf[i,iend1+(buf-(j+1))] = mf[i,iend1+(buf-(j+1))] - mf[i,iend1+(buf-(j+1))]*(float(j+1)/buf)
        # sys.exit()

    # right
    if ibeg2 < ncol:
        for j in range(buf):
          for i in range(jbeg,jend):
            mf[i,ibeg2-(buf-j)] = mf[i,ibeg2-(buf-j)] - mf[i,ibeg2-(buf-j)]*(float(j+1)/buf)

    mf[jbeg1:jend1,:] = float('NaN')
    mf[jbeg2:jend2,:] = float('NaN')
    mf[:,ibeg1:iend1] = float('NaN')
    mf[:,ibeg2:iend2] = float('NaN')
#plt.imshow(mf)
#plt.show()
#sys.exit()

# apply  scale and manual cst
mf = scale*(mf - shift)

# clean based on perc
maxlos,minlos=np.nanpercentile(mf,perc),np.nanpercentile(mf,(100-perc))
print('Clean outliers outside:', maxlos,minlos)
kk = np.nonzero(
    np.logical_or(mf<minlos,mf>maxlos))
mf[kk] = float('NaN')

# Plot
vmax = np.max([np.nanpercentile(mf,99),np.nanpercentile(m,99),np.nanpercentile(mf,1),np.nanpercentile(m,1)])
vmin = -vmax

if setzero == 'yes':
    print('HELLO')
    # replace "NaN" to 0
    mf[np.isnan(mf)] = 0.

#save output
fid3 = open(outfile,'wb')
mf.flatten().astype('float32').tofile(fid3)
fid3.close()

fig = plt.figure(0,figsize=(16,6))
ax = fig.add_subplot(1,4,1)
cax = ax.imshow(m, cm.RdBu,vmin=vmin, vmax=vmax)
ax.set_title(infile)
setp( ax.get_xticklabels(), visible=False)
#cbar = fig.colorbar(cax, orientation='vertical',aspect=9)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

ax = fig.add_subplot(1,4,2)
cax = ax.imshow(mf_ramp, cm.RdBu, vmin=vmin, vmax=vmax)
ax.set_title('DATA - RAMP')
setp( ax.get_xticklabels(), visible=False)
#cbar = fig.colorbar(cax, orientation='vertical',aspect=9)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

ax = fig.add_subplot(1,4,3)
cax = ax.imshow(mask, cm.RdBu, vmin=0, vmax=np.nanpercentile(mask,99))
ax.set_title('MASK')
setp( ax.get_xticklabels(), visible=False)
#cbar = fig.colorbar(cax, orientation='vertical',aspect=9)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

#if iend-ibeg<ncol and jend-jbeg<nlign:
#    ax = fig.add_subplot(1,4,4)
#    cax = ax.imshow(mf[jbeg-buf:jend+buf,ibeg-buf:iend+buf], cm.RdBu, vmin=vmin, vmax=vmax)
#    setp( ax.get_xticklabels(), visible=False)
#    #cbar = fig.colorbar(cax, orientation='vertical',aspect=9)
#    ax.set_title(outfile)
#    divider = make_axes_locatable(ax)
#    c = divider.append_axes("right", size="5%", pad=0.05)
#    plt.colorbar(cax, cax=c)
#else:

ax = fig.add_subplot(1,4,4)
cax = ax.imshow(mf, cm.RdBu, vmin=vmin, vmax=vmax)
setp( ax.get_xticklabels(), visible=False)
setp( ax.get_xticklabels(), visible=False)
#cbar = fig.colorbar(cax, orientation='vertical',aspect=9)
ax.set_title(outfile)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

fig.tight_layout()
plt.show()
