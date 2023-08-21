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
clean_ts.py
-------------
Clean a time series file (cube in binary format) given an other real4 file (mask) and a threshold on this mask

Usage: clean_ts.py [--infile=<path> ] [--outfile=<path>]  [--mask=<path>] [--threshold=<value>]  \
[--perc=<value>] [--vmin=<value>] [--vmax=<value>] [--rampmask=<yes/no>] \
[--flatten_mask=<path>] [--lectfile=<path>] [--scale=<value>] [--imref=<value>] \
[--images=<path>] [--clean=<values>]  [--crop=<values>] [--clean_demerr=<path>] 

Options:
-h --help           Show this screen.
--infile PATH       path to time series [default: depl_cumule]
--outfile PATH      output file [default: depl_cumule_clean]
--mask PATH         r4 file used as mask [default: None]
--flatten_mask PATH output r4 flatten mask [default: None]
--rampmask VALUE    flatten mask [default: yes]
--threshold VALUE   threshold value on mask file (Keep pixel with mask < threshold)
--scale VALUE       scale the mask [default:1]
--lectfile PATH     Path of the lect.in file [default: lect.in]
--imref VALUE       Reference image number [default: 1]
--images PATH       Path to image_retuenues file [default: images_retenues]
--clean VALUE       Clean option [default: 0,0,0,0]
--crop VALUE        Crop option [default: 0,nlines,0,ncol]
--vmax              Max colorscale [default: 98th percentile]
--vmin              Min colorscale [default: 2th percentile]
--perc VALUE        Percentile of hidden LOS pixel for the estimation and clean outliers [default:99.9]
--clean_demerr      Path to dem error file
"""

print()
print()
print('Author: Simon Daout') 
print()
print()

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided

import os
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from osgeo import gdal
import shutil

# scipy
import scipy
import scipy.optimize as opt
import scipy.linalg as lst

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning) 

import docopt
arguments = docopt.docopt(__doc__)
if arguments["--infile"] ==  None:
    infile = 'depl_cumule'
else:
    infile = arguments["--infile"]
if arguments["--outfile"] ==  None:
    outfile = 'depl_cumule_clean'
else:
    outfile = arguments["--outfile"]    
maskf = arguments["--mask"]
if arguments["--threshold"] == None:
    seuil = 1
else:
    seuil = float(arguments["--threshold"])
if arguments["--rampmask"] ==  None:
   rampmask = "no"
else:
   rampmask = arguments["--rampmask"]
if arguments["--lectfile"] ==  None:
   lecfile = "lect.in"
else:
   lecfile = arguments["--lectfile"]
if arguments["--scale"] ==  None:
   scale = 1.
else:
   scale = float(arguments["--scale"])
if arguments["--imref"] ==  None:
    imref = 0
elif arguments["--imref"] < 1:
    print('--imref must be between 1 and Nimages')
else:
    imref = int(arguments["--imref"]) - 1

# lect cube
ds = gdal.Open(infile)
if not ds:
  print('.hdr file time series cube {0}, not found, open {1}'.format(infile,lecfile))
  ncol, nlines = list(map(int, open(lecfile).readline().split(None, 2)[0:2]))
else:
  hdr = infile + '.hdr'
  print('Read ', hdr)
  ncol, nlines = ds.RasterXSize, ds.RasterYSize
  driver =  ds.GetDriver()
  
if arguments["--clean"] ==  None:
    mask = [0,0,0,0]
else:
    print(arguments["--clean"])
    mask = map(float,arguments["--clean"].replace(',',' ').split())
mibeg,miend,mjbeg,mjend = int(mask[0]),int(mask[1]),int(mask[2]),int(mask[3])

if arguments["--crop"] ==  None:
    crop = [0,nlines,0,ncol]
else:
    crop = map(float,arguments["--crop"].replace(',',' ').split())
ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])

if arguments["--perc"] ==  None:
    perc = 99.9
else:
    perc = float(arguments["--perc"])

if arguments["--clean_demerr"] ==  None:
    demf = 'no'
    dem = np.zeros((iend-ibeg,jend-jbeg))
else:
    demf = arguments["--clean_demerr"]
    extension = os.path.splitext(demf)[1]
    if extension == ".tif":
        ds = gdal.Open(demf, gdal.GA_ReadOnly)
        dem = ds.GetRasterBand(1).ReadAsArray()
    else:
        dem = np.fromfile(demf,dtype=np.float32)[:nlines*ncol].reshape((nlines,ncol))[ibeg:iend,jbeg:jend]

# load images_retenues file
fimages='images_retenues'
nb,idates,dates,base=np.loadtxt(fimages, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
# nb images
N=len(dates)
print('Number images: ', N)

# open mask file
mask = np.zeros((iend-ibeg,jend-jbeg))
if maskf is not None:
    fid2 = open(maskf, 'r')
    mask = np.fromfile(fid2,dtype=np.float32)[:nlines*ncol].reshape((nlines,ncol))[ibeg:iend,jbeg:jend]
    mask =  mask*scale
    kk = np.nonzero(np.logical_or(mask==0.0, mask>999.))
    mask[kk] = float('NaN')
    # mask[kk] = 999

if rampmask=='yes':

    index = np.nonzero( ~np.isnan(mask))
    temp = np.array(index).T
    maski = mask[index].flatten()

    az = temp[:,0]; rg = temp[:,1]
    G=np.zeros((len(maski),5))
    G[:,0] = rg**2
    G[:,1] = rg
    G[:,2] = az**2
    G[:,3] = az
    G[:,4] = 1

    x0 = lst.lstsq(G,maski)[0]
    _func = lambda x: np.sum(((np.dot(G,x)-maski))**2)
    _fprime = lambda x: 2*np.dot(G.T, (np.dot(G,x)-maski))
    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]    

    pars = np.dot(np.dot(np.linalg.inv(np.dot(G.T,G)),G.T),maski)
    a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
    print('Remove ramp mask %f x**2 %f x  + %f y**2 + %f y + %f for : %s'%(a,b,c,d,e,maskf))

    G=np.zeros((len(mask.flatten()),5))
    for i in range(iend-ibeg):
        G[i*ncol:(i+1)*ncol,0] = np.arange((ncol))**2
        G[i*ncol:(i+1)*ncol,1] = np.arange((ncol))
        G[i*ncol:(i+1)*ncol,2] = i**2
        G[i*ncol:(i+1)*ncol,3] = i
    G[:,4] = 1
    temp = (mask.flatten() - np.dot(G,pars))
    maskflat=temp.reshape(iend-ibeg,jend-jbeg)
    # maskflat = (temp - np.nanmin(temp)).reshape(nlines,ncol)

else:
    # pass
    maskflat = np.copy(mask)

if maskf is not None:
    vmax = np.nanpercentile(maskflat,92)
    vmin= np.nanpercentile(maskflat,8)

    kk = np.nonzero(maskflat>seuil)
    spacial_mask = np.copy(maskflat)
    spacial_mask[kk] = float('NaN')
    nfigure=0
    fig = plt.figure(0,figsize=(11,4))
    ax = fig.add_subplot(1,3,1)
    cax = ax.imshow(mask,cmap=cm.jet,vmax=vmax,vmin=0,interpolation='none')
    ax.set_title('RMSpixel')
    plt.setp( ax.get_xticklabels(), visible=False)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    ax = fig.add_subplot(1,3,2)
    cax = ax.imshow(maskflat,cmap=cm.jet,vmax=vmax,vmin=0,interpolation='none')
    ax.set_title('flat RMSpixel')
    plt.setp( ax.get_xticklabels(), visible=False)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    ax = fig.add_subplot(1,3,3)
    cax = ax.imshow(spacial_mask,cmap=cm.jet,vmax=vmax,vmin=0,interpolation='none')
    ax.set_title('Mask')
    plt.setp( ax.get_xticklabels(), visible=False)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)
    del spacial_mask, mask
    #plt.show()
    #sys.exit()

if arguments["--flatten_mask"] != None:
    # save clean ts
    fid = open(arguments["--flatten_mask"], 'wb')
    maskflat.flatten().astype('float32').tofile(fid)
    fid.close()

cubei = np.fromfile(infile,dtype=np.float32)
cube = as_strided(cubei[:nlines*ncol*N])
kk = np.flatnonzero(np.logical_or(cube==9990, cube==9999))
cube[kk] = float('NaN')
print('Number of line in the cube: ', cube.shape)

maps = cube.reshape((nlines,ncol,N))[ibeg:iend,jbeg:jend,:N]
print('Reshape cube: ', maps.shape)

# clean 
cst = np.copy(maps[:,:,imref])
for l in range((N)):
    d = as_strided(maps[:,:,l])
    
    _d=np.copy(d)
    _d[d==0] = float('NaN')
    maxlos,minlos=np.nanpercentile(_d,perc),np.nanpercentile(_d,(100-perc))
    # print('Min, Max LOS:', minlos, maxlos)
    # set at NaN zero values for all dates
    kk = np.nonzero(
    np.logical_or(d<minlos,
    np.logical_or(d>maxlos,
        np.isnan(maskflat)
        )))
    
    d[kk] = float('NaN')
    # carefull stupid unit 
    maps[:,:,l] = maps[:,:,l] - cst - dem*(base[l]-base[imref])/100.
    if l != imref:
        index = np.nonzero(d==0.0)
        d[index] = float('NaN')
    maps[mibeg:miend,mjbeg:mjend,l] = float('NaN')

# _cube=np.copy(cube)
# _cube[cube==0] = float('NaN')
# maxlos,minlos=np.nanpercentile(_cube,perc),np.nanpercentile(_cube,(100-perc))
# print('Min, Max LOS:', minlos, maxlos)

# save clean ts
#if not ds:
fid = open(outfile, 'wb')
maps.flatten().astype('float32').tofile(fid)
fid.close()
wf = open("lect_clean.in", "w")
wf.write("%i %i %i\n" % (jend-jbeg, iend-ibeg, N))
wf.close()
#else:
#    dst_ds = driver.Create(outfile, jend-jbeg, iend-ibeg, N, gdal.GDT_Float32)
#    for i in range(N):
#        dst_band = dst_ds.GetRasterBand(i+1)
#        dst_band.WriteArray(maps[:,:,i],0,0)
#        dst_band.FlushCache()

# plot diplacements maps
fig = plt.figure(1,figsize=(14,10))
fig.subplots_adjust(wspace=0.001)

if arguments["--vmax"] ==  None:
    vmax = np.nanpercentile(maps, 98)
else:
    vmax = float(arguments["--vmax"])

if arguments["--vmin"] ==  None:
    vmin = np.nanpercentile(maps, 2)
else:
    vmin = float(arguments["--vmin"])

for l in range((N)):
    d = as_strided(maps[:,:,l])
    ax = fig.add_subplot(4,int(N/4)+1,l+1)
    cmap = cm.jet
    cmap.set_bad('white')
    cax = ax.imshow(d,cmap=cm.jet,vmax=vmax,vmin=vmin,interpolation='none')
    ax.set_title(idates[l],fontsize=6)
    plt.setp( ax.get_xticklabels(), visible=False)
    plt.setp( ax.get_yticklabels(), visible=False)

plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)
# fig.tight_layout()
plt.suptitle('Time series maps')
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
plt.show()

