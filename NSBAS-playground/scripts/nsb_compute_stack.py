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
compute_stack.py
-------------
Compute a stack from a list of interferograms

Usage: compute_stack.py --homedir=<path> --int_path=<path> --int_list=<file> [--suffix=<value>] [--prefix=<value>] [--output=<file>] [--rlook=<value>] [--ref=<jstart,jend,istart,iend>] [--ramp=<yes/no>]  --radar_file=<file>
compute_stack.py -h | --help

Options:
-h --help           Show this screen
--homedir PATH      Path to the Home directory
--int_path PATH     Path to the interferograms from homedir
--int_list FILE     List of interferograms
--suffix VALUE      Suffix interferograms [default: ""]
--prefix VALUE      Prefix interferograms [default: ""]
--rlook VALUE       Look interferograms [default: 2]
--radar_file FILE   Path to the radar.hgt file
--ref=<jstart,jend,istart,iend> Set to zero displacements from jstart to jend
--ramp<yes/no>      Correct the map from ramp in range and azimuth before stack computation
--output
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
# gdal
from osgeo import gdal
import os
gdal.UseExceptions()
import docopt
import scipy.optimize as opt
import scipy.linalg as lst

print()
print()
print('Author: Simon DAOUT')
print('Please cite:')
print('Daout, S., Sudhaus, H., Kausch, T., Steinberg, A., & Dini, B. (2019). Interseismic and postseismic shallow creep of the North Qaidam Thrust faults detected with a multitemporal InSAR analysis. Journal of Geophysical Research: Solid Earth, 124(7), 7259-7279.')
print()
print()

# read arguments
arguments = docopt.docopt(__doc__)
arguments = docopt.docopt(__doc__)
home=os.path.abspath(arguments["--homedir"])  + '/'
int_path=os.path.join(home, arguments["--int_path"]) + '/'
int_list = os.path.join(home,arguments["--int_list"])
radar=os.path.join(home, arguments["--radar_file"])

if arguments["--prefix"] == None:
    prefix = ''
else:
    prefix=arguments["--prefix"]
if arguments["--suffix"] == None:
    suffix = ''
else:
    suffix=arguments["--suffix"]
if arguments["--rlook"] == None:
    rlook = '_2'
else:
    rlook = '_' + arguments["--rlook"] 
if arguments["--output"] == None:
    output = stack.r4
else:
    output=arguments["--output"]
if arguments["--ref"] == None:
    lin_start, lin_jend, col_start, col_jend = None,None,None,None
else:
    ref = list(map(int,arguments["--ref"].replace(',',' ').split()))
    try:
        lin_start,lin_end, col_start, col_end = ref[0], ref[1], ref[2], ref[3]
    except:
        lin_start,lin_end = ref[0], ref[1]
        col_start, col_end = 0, ncol

def remove_ramp(los,nlign,ncol):
        index = np.nonzero(np.logical_and(np.logical_and(~np.isnan(los),los<np.percentile(los,80)),los>np.percentile(los,20)))
        temp = np.array(index).T
        mi = los[index].flatten()
        az = temp[:,0]; rg = temp[:,1]

        G=np.zeros((len(mi),5))
        G[:,0] = rg**2
        G[:,1] = az**2
        G[:,2] = rg
        G[:,3] = az
        G[:,4] = 1

        x0 = lst.lstsq(G,mi)[0]
        _func = lambda x: np.sum(((np.dot(G,x)-mi))**2)
        _fprime = lambda x: 2*np.dot(G.T, (np.dot(G,x)-mi))
        pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0)[0]
        a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; f = pars[4]
        print('Remove ramp %f x**2 +  %f y**2  +  %f x  + %f y + %f'%(a,b,c,d,f))

        G=np.zeros((len(los.flatten()),5))
        for i in range(nlign):
           G[i*ncol:(i+1)*ncol,0] = np.arange((ncol))**2
           G[i*ncol:(i+1)*ncol,1] = i**2
           G[i*ncol:(i+1)*ncol,2] = np.arange((ncol))
           G[i*ncol:(i+1)*ncol,3] = i
        G[:,4] = 1
        mf = (los.flatten() - np.dot(G,pars)).reshape(nlign,ncol)
        return mf, index

# list of dates
date1,date2=np.loadtxt(int_list,comments="#",unpack=True,dtype='i,i')
kmax = len(date1)

# radar file
driver = gdal.GetDriverByName("roi_pac")
ds = gdal.OpenEx(radar, allowed_drivers=["ROI_PAC"])
nlines,ncols= ds.RasterYSize, ds.RasterXSize

alllos = np.zeros((nlines,ncols,kmax))
allcor = np.zeros((nlines,ncols,kmax))

for i in range((kmax)):

    interf1,interf2 = date1[i],date2[i]
    los_map = np.zeros((nlines,ncols))
    cor_map = np.zeros((nlines,ncols))

    folder= int_path  + 'int_' + str(interf1) + '_' + str(interf2) + '/'    
    infile=folder + str(prefix) + str(interf1) + '-' + str(interf2) + str(suffix) + str(rlook) + 'rlks.unw'

    ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])
    ds_band2 = ds.GetRasterBand(2)
    ds_band1 = ds.GetRasterBand(1)
    print('Nlines:{}, Ncol:{}, int:{}-{}'.format(ds.RasterYSize, ds.RasterXSize,interf1,interf2))

    los_map[:ds.RasterYSize,:ds.RasterXSize] = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:nlines,:ncols]
    cor_map[:ds.RasterYSize,:ds.RasterXSize] = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:nlines,:ncols]
    
    #fig = plt.figure(0,figsize=(16,6))
    #ax = fig.add_subplot(1,3,1)
    #cax = ax.imshow(los_map,cm.jet)
    #los_temp = np.copy(los_map)  
    
    if arguments["--ramp"] == 'yes':
      los_map,index2 = remove_ramp(los_map,ds.RasterYSize,ds.RasterXSize)
    
    if (arguments["--ref"] != None) :
          los_ref = los_map[lin_start:lin_end,col_start:col_end].flatten()
          cor_ref = cor_map[lin_start:lin_end,col_start:col_end].flatten()
          index = np.nonzero(los_ref != 0)
          amp_ref = cor_ref[index]
          cst = np.nansum(los_ref[index]*amp_ref) / np.nansum(amp_ref)
          print(cst)
          los_map = los_map - cst

    #ax = fig.add_subplot(1,3,2)
    #cax = ax.imshow(los_map,cm.jet)
    #los_temp[index2] = np.nan
    #ax = fig.add_subplot(1,3,3)
    #cax = ax.imshow(los_temp,cm.jet)
    #plt.show()

    alllos[:,:,i] = los_map
    allcor[:,:,i] = cor_map

alllos[alllos==0] = float('NaN')
sumlos = np.nanmean(alllos,axis=2)

#sumlos=np.zeros((nlines,ncols))
#for i in range(0,nlines,1):
#    for j in range(0,ncols,1):
#        k=np.flatnonzero(alllos[i,j,:] != 0)
#        if len(k) > 0:
#            #sumlos[i,j]=np.sum(alllos[i,j,k[:]])/len(k)
#            sumlos[i,j]=np.sum(alllos[i,j,k[:]]*allcor[i,j,k[:]])/np.sum(allcor[i,j,k[:]])
               
sumlos.astype('float32').tofile(home + output)

# cmap=cm.jet
try:
        from matplotlib.colors import LinearSegmentedColormap
        cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
        cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
        cmap = cmap.reversed()
except:
        cmap=cm.rainbow

fig = plt.figure(figsize=(12,5))
ax = fig.add_subplot(1,1,1)
ax.set_title('colorMap')
plt.imshow(sumlos,vmax=np.nanpercentile(sumlos,98),vmin=np.nanpercentile(sumlos,2),cmap=cmap)
ax.set_aspect('equal')

cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)
plt.colorbar(orientation='vertical')
plt.show()
