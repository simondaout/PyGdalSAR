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
clean_raster.py
-------------
Clean a raster file given an other r4 file (mask) and a threshold on this mask

Usage: clean_raster.py --infile=<path> --outfile=<path>  [--mask=<path>] [--threshold=<value>] \
[--perc=<value>] [--crop=<values>] [--buff=<value>] [--lectfile=<path>] [--scale=<value>] [--scale_mask=<value>] [--ramp=<lin/quad/cub/4/no>] [--ref=<jstart,jend,istart,iend>] [--removeNAN=<yes/no>] [--cst=<value>] [--absolute=<yes/no>]  [--reverse=<yes/no>] [--filter=<HP/LP>] [--fwindsize=<value>]

Options:
-h --help           Show this screen.
--infile PATH       r4 file to clean
--outfile PATH      output file
--mask PATH         r4 file used as mask
--threshold VALUE   threshold value on mask file (Keep pixel with mask > threshold)
--reverse yes/no    if yes mask values above mask threshold (default: no)
--absolute yes/no    if yes, take the aboslute values of the mask (default: no)
--scale VALUE       scale data [default:1]
--scale_mask VALUE  scale the mask [default:1]
--lectfile PATH     Path of the lect.in file [default: lect.in]
--crop VALUE       Crop option with smoothing of boundaries [default: 0,ncols,0,nlines]
--buff VALUE        Number of pixels for crop smothing  (default: 50)
--perc VALUE        Percentile of hidden LOS pixel for the estimation and clean outliers [default:100]
--ramp<lin/quad/cub/no>      Correct the map from ramp in range and azimuth
--ref=<jstart,jend,istart,iend> Set to zero displacements from jstart to jend
--removeNAN         replace NaN by 0
--cst               Add constante to map
--filter=<HP/LP> Apply a high pass (HP) or a low pass (LP) gaussian filter to the image
--fwindsize=<value> Filter window size (default: 16)
"""

print()
print()
print('Author: Simon DAOUT')
print()
print()

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided
import scipy.optimize as opt
import scipy.linalg as lst
from scipy import ndimage

import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
from pylab import *
import os

import docopt
arguments = docopt.docopt(__doc__)
infile = arguments["--infile"]
outfile = arguments["--outfile"]
if arguments["--mask"] ==  None:
    maskf = 'no'
else:
    maskf = arguments["--mask"]
if arguments["--threshold"] ==  None:
    seuil = 1e10
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
    perc = 100
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
elif arguments["--ramp"] == '4':
    ramp = '4'
elif arguments["--ramp"] == 'yes':
    print('ramp==yes is depricated. Use lin/quad/cub/4/no. Set ramp to lin')
    ramp = 'lin'
else:
    print('ramp argument not recognized. Exit!')
    sys.exit()
if arguments["--reverse"] == 'yes':
   reverse = True
else:
   reverse = False
if arguments["--absolute"] == 'yes':
   absolute = True
else:
   absolute = False

if arguments["--fwindsize"] == None:
    arguments["--fwindsize"] = 16

ds_extension = os.path.splitext(infile)[1]
if (ds_extension == ".tif" or ds_extension ==".tiff" or ds_extension ==".grd"):
    from osgeo import gdal
    sformat = 'GRID'
    ds = gdal.Open(infile, gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    # Attributes
    print("> Driver:   ", ds.GetDriver().ShortName)
    print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
    print("> Datatype: ", gdal.GetDataTypeName(band.DataType))
    nlines, ncols = ds.RasterYSize, ds.RasterXSize
    gt = ds.GetGeoTransform()
    proj = ds.GetProjectionRef()
    driver = gdal.GetDriverByName('GTiff')
    m = band.ReadAsArray(0, 0,
                               ds.RasterXSize, ds.RasterYSize,
                               ds.RasterXSize, ds.RasterYSize)
else:
    sformat = 'R4'
    # read lect.in 
    ncols, nlines = list(map(int, open(lecfile).readline().split(None, 2)[0:2]))
    fid = open(infile, 'r')
    m = np.fromfile(fid,dtype=float32).reshape((nlines,ncols))
#    m[m==0.0] = float('NaN')

if arguments["--ref"] == None:
    lin_start, lin_jend, col_start, col_jend = None,None,None,None
else:
    ref = list(map(int,arguments["--ref"].replace(',',' ').split()))
    try:
        lin_start,lin_end, col_start, col_end = ref[0], ref[1], ref[2], ref[3]
    except:
        lin_start,lin_end = ref[0], ref[1]
        col_start, col_end = 0, ncols

# mask
if maskf != 'no':
  ds_extension = os.path.splitext(maskf)[1]
  if (ds_extension == ".tif" or ds_extension ==".tiff" or ds_extension ==".grd"):
    from osgeo import gdal
    ds = gdal.Open(maskf, gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    # Attributes
    print("> Driver:   ", ds.GetDriver().ShortName)
    print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
    print("> Datatype: ", gdal.GetDataTypeName(band.DataType))
    mask = band.ReadAsArray(0, 0,
                               ds.RasterXSize, ds.RasterYSize,
                               ds.RasterXSize, ds.RasterYSize)
  else:
    fid2 = open(maskf, 'r')
    mask = np.fromfile(fid2,dtype=float32)[:nlines*ncols].reshape((nlines,ncols))
    mask =  mask*scale_mask
    mask[np.isnan(mask)] = 0
else:
    mask = np.zeros((nlines,ncols))*float('nan')

mf = np.copy(m)
# apply mask
if absolute:
      mask = np.abs(mask)
if reverse:
      kk = np.nonzero(np.logical_or(mask==0, mask<seuil)) 
else:
      kk = np.nonzero(np.logical_or(mask==0, mask>seuil)) 
mf[kk] = float('NaN')

# apply  scale and manual cst
mf = scale*(mf - shift)

# apply high pass filter
if arguments["--filter"] == 'HP':
    # make array with nans replaced by zeros and filter it
    m_filter_vals = np.copy(mf)
    m_filter_vals[np.isnan(mf)] = 0.
    m_lp_vals = ndimage.gaussian_filter(m_filter_vals, int(arguments["--fwindsize"]))
    # make same size array full of ones, but set to zero where there is a nan in mf
    m_filter_ones = 0*np.copy(mf)+1
    m_filter_ones[np.isnan(mf)] = 0.
    m_lp_ones = ndimage.gaussian_filter(m_filter_ones, int(arguments["--fwindsize"]))
    # find the ratio to make coefficients sum to one near nan values
    m_lp = m_lp_vals/m_lp_ones
    # subtract to obtain hp filter
    mf = mf - m_lp

elif arguments["--filter"] == 'LP':
    m_filter = np.copy(mf)
    index = np.isnan(mf)
    m_filter[index] = 0.
    mf = ndimage.gaussian_filter(m_filter, int(arguments["--fwindsize"]))
    mf[index] = float('nan')

if ramp != 'no':
    # clean for ramp
    maxlos,minlos=np.nanpercentile(m,99.5),np.nanpercentile(m,.5)
    print('Clean outliers for ramp estimation outside:', maxlos,minlos)
    kk = np.nonzero(
        np.logical_or(m<minlos,m>maxlos))
    m_ramp = np.copy(mf)
    m_ramp[kk] = float('NaN')

if ramp == 'lin':
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
    for i in range(nlines):
        G[i*ncols:(i+1)*ncols,0] = np.arange((ncols))
        G[i*ncols:(i+1)*ncols,1] = i
    G[:,2] = 1

    # apply ramp correction
    temp = np.dot(G,pars)
    ramp=temp.reshape(nlines,ncols)
    mf = mf - ramp
    del temp

elif ramp=='quad':
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
    for i in range(nlines):
        G[i*ncols:(i+1)*ncols,0] = np.arange((ncols))**2
        G[i*ncols:(i+1)*ncols,1] = np.arange((ncols))
        G[i*ncols:(i+1)*ncols,2] = i
    G[:,3] = 1
    # apply ramp correction
    temp = np.dot(G,pars)
    ramp=temp.reshape(nlines,ncols)
    mf = mf - ramp
    del temp

elif ramp=='cub':
    index = np.nonzero(~np.isnan(m_ramp))
    temp = np.array(index).T
    mi = m[index].flatten()
    az = temp[:,0]; rg = temp[:,1]

    G=np.zeros((len(mi),7))
    G[:,0] = rg**3
    G[:,1] = rg**2
    G[:,2] = rg
    G[:,3] = az**3
    G[:,4] = az**2
    G[:,5] = az
    G[:,6] = 1

    x0 = lst.lstsq(G,mi)[0]
    _func = lambda x: np.sum(((np.dot(G,x)-mi))**2)
    _fprime = lambda x: 2*np.dot(G.T, (np.dot(G,x)-mi))
    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0)[0]    
    #pars = np.dot(np.dot(np.linalg.inv(np.dot(G.T,G)),G.T),mi)
    a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6] 
    print('Remove ramp %f x**3 + %f x**2 + %f x  + %f y**3 + %f y**2 + %f y + %f'%(a,b,c,d,e,f,g))

    G=np.zeros((len(mask.flatten()),7))
    for i in range(nlines):
        G[i*ncols:(i+1)*ncols,0] = np.arange((ncols))**3
        G[i*ncols:(i+1)*ncols,1] = np.arange((ncols))**2
        G[i*ncols:(i+1)*ncols,2] = np.arange((ncols))
        G[i*ncols:(i+1)*ncols,3] = i**3
        G[i*ncols:(i+1)*ncols,4] = i**2
        G[i*ncols:(i+1)*ncols,5] = i
    G[:,6] = 1
    # apply ramp correction
    temp = np.dot(G,pars)
    ramp=temp.reshape(nlines,ncols)
    mf = mf - ramp
    del temp

elif ramp=='4':
    index = np.nonzero(~np.isnan(m_ramp))
    temp = np.array(index).T
    mi = m[index].flatten()
    az = temp[:,0]; rg = temp[:,1]

    G=np.zeros((len(mi),9))
    G[:,0] = rg**3
    G[:,1] = rg**2
    G[:,2] = rg
    G[:,3] = az**3
    G[:,4] = az**2
    G[:,5] = (az*rg)**2
    G[:,6] = (az*rg)
    G[:,7] = az
    G[:,8] = 1

    x0 = lst.lstsq(G,mi)[0]
    _func = lambda x: np.sum(((np.dot(G,x)-mi))**2)
    _fprime = lambda x: 2*np.dot(G.T, (np.dot(G,x)-mi))
    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0)[0]    
    #pars = np.dot(np.dot(np.linalg.inv(np.dot(G.T,G)),G.T),mi)
    a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]; h = pars[7]; i = pars[8] 
    print('Remove ramp %f x**3 + %f x**2 + %f x  + %f y**3 + %f y**2 + %f xy**2 + %f xy + %f y + %f'%(a,b,c,d,e,f,g,h,i))

    G=np.zeros((len(mask.flatten()),9))
    for i in range(nlines):
        G[i*ncols:(i+1)*ncols,0] = np.arange((ncols))**3
        G[i*ncols:(i+1)*ncols,1] = np.arange((ncols))**2
        G[i*ncols:(i+1)*ncols,2] = np.arange((ncols))
        G[i*ncols:(i+1)*ncols,3] = i**3
        G[i*ncols:(i+1)*ncols,4] = i**2
        G[i*ncols:(i+1)*ncols,5] = (i*np.arange((ncols)))**2
        G[i*ncols:(i+1)*ncols,6] = i*np.arange((ncols))
        G[i*ncols:(i+1)*ncols,7] = i
    G[:,8] = 1
    # apply ramp correction
    temp = np.dot(G,pars)
    ramp=temp.reshape(nlines,ncols)
    mf = mf - ramp
    del temp

else:
    ramp= np.zeros(np.shape(m))

if (arguments["--ref"] != None) :
    cst = np.nanmean(mf[lin_start:lin_end,col_start:col_end])
    mf = mf - cst
    ramp = ramp - cst

# crop
if arguments["--crop"] ==  None:
    crop = [0,ncols,0,nlines]
else:
    crop = list(map(float,arguments["--crop"].replace(',',' ').split()))
ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])
if iend>ncols:
    iend=ncols
if jend>nlines:
    jend=nlines

if iend-ibeg<ncols or jend-jbeg<nlines:
    if arguments["--buff"] ==  None:
        buf = 50
    else:
        buf = int(arguments["--buff"])
    # ibeg, iend = 250, 900
    # jbeg, jend = 2000, 2900 
    ibeg1,iend1 = 0, ibeg
    ibeg2,iend2 = iend,ncols
    jbeg1,jend1 = 0, jbeg
    jbeg2,jend2 = jend, nlines

    # top 
    if jend1 > 0:
        for j in range(buf):
          for i in range(ibeg,iend):
                mf[jend1+j,i] = mf[jend1+j,i]*(float(j+1)/buf)
    
    #bottom
    if jbeg2 < nlines:
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
    if ibeg2 < ncols:
        for j in range(buf):
          for i in range(jbeg,jend):
            mf[i,ibeg2-(buf-j)] = mf[i,ibeg2-(buf-j)] - mf[i,ibeg2-(buf-j)]*(float(j+1)/buf)

    mf[jbeg1:jend1,:] = float('NaN')
    mf[jbeg2:jend2,:] = float('NaN')
    mf[:,ibeg1:iend1] = float('NaN')
    mf[:,ibeg2:iend2] = float('NaN')

# clean based on perc
maxlos,minlos=np.nanpercentile(mf,perc),np.nanpercentile(mf,(100-perc))
print('Clean outliers outside {}-{} with perc:{}:'.format(maxlos,minlos,perc))
kk = np.nonzero(
    np.logical_or(mf<minlos,mf>maxlos))
#mf[kk] = float('NaN')
mf[kk] = 0

# Plot
vmax = np.nanpercentile(mf,98)
vmin = np.nanpercentile(mf,2)

if setzero == 'yes':
    # replace "NaN" to 0
    mf[np.isnan(mf)] = 0.
    #mf[np.isnan(m)] = np.nan

#save output
if sformat == 'R4':
    fid3 = open(outfile,'wb')
    mf.flatten().astype('float32').tofile(fid3)
    fid3.close()
else:
   dst_ds = driver.Create(outfile, ncols, nlines, 1, gdal.GDT_Float32)
   dst_band2 = dst_ds.GetRasterBand(1)
   dst_band2.WriteArray(mf,0,0)
   dst_ds.SetGeoTransform(gt)
   dst_ds.SetProjection(proj)
   dst_band2.FlushCache()
   del dst_ds

try:
  from matplotlib.colors import LinearSegmentedColormap
  cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
  cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
  cmap = cmap.reversed()
except:
  cmap=cm.rainbow

fig = plt.figure(0,figsize=(12,6))
ax = fig.add_subplot(2,2,1)
cax = ax.imshow(m, cmap,vmin=vmin, vmax=vmax, interpolation='nearest')
ax.set_title(infile)
setp( ax.get_xticklabels(), visible=False)
#cbar = fig.colorbar(cax, orientation='vertical',aspect=9)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

ax = fig.add_subplot(2,2,2)
cax = ax.imshow(mask, cmap, vmin=0, vmax=np.nanpercentile(mask,99),interpolation='nearest')
ax.set_title('MASK')
setp( ax.get_xticklabels(), visible=False)
#cbar = fig.colorbar(cax, orientation='vertical',aspect=9)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

ax = fig.add_subplot(2,2,3)
cax = ax.imshow(ramp, cmap, vmin=vmin, vmax=vmax,interpolation='nearest')
ax.set_title('DATA - RAMP')
setp( ax.get_xticklabels(), visible=False)
#cbar = fig.colorbar(cax, orientation='vertical',aspect=9)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

ax = fig.add_subplot(2,2,4)
cax = ax.imshow(mf, cmap, vmin=vmin, vmax=vmax,interpolation='nearest')
setp( ax.get_xticklabels(), visible=False)
setp( ax.get_xticklabels(), visible=False)
#cbar = fig.colorbar(cax, orientation='vertical',aspect=9)
ax.set_title(outfile)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

fig.tight_layout()

try:
    fig.savefig('{}_clean.pdf'.format(infile), format='PDF',dpi=180)
except:
    pass

plt.show()
