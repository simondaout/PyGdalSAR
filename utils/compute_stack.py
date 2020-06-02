#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 11:42:19 2017

@author: thorek, daouts
"""

import numpy as np
import matplotlib.pyplot as plt

# gdal
import gdal
gdal.UseExceptions()

home = '/data3/work/tibet/qinghai/processing/T455/'
suffix= '_flatrange_flatcutbox_'
preffix = 'filt_coh_'
nlign,ncol=4082,1420

# home = '/data3/work/tibet/qinghai/processing/T319/'
# suffix= '_flatrange_flatcutbox_invert_'
# preffix = 'filt_coh_cor_'
# rlook = 4
# nlign,ncol=4078,1420

# home = '/data3/T047/'
# suffix= '_flatrange_flatcutbox_'
# preffix = 'filt_coh_'
# nlign,ncol=4077,1420

int_list= home + 'interf_pair_withEq1.rsc'
output = 'phs_sum_eq1'
# int_list= home + 'interf_pair_withEq2.rsc'
# output = 'phs_sum_eq2'
# int_list= home + 'interf_pair_withEq12.rsc'
# output = 'phs_sum_eq12'
rlook = 4



date1,date2=np.loadtxt(int_list,comments="#",unpack=True,dtype='i,i')
kmax = len(date1)

alllos = np.zeros((nlign,ncol,kmax))
allcor = np.zeros((nlign,ncol,kmax))

for i in xrange((kmax)):

    interf1,interf2 = date1[i],date2[i]
    los_map = np.zeros((nlign,ncol))
    cor_map = np.zeros((nlign,ncol))

    folder=home + 'int/' + 'int_' + str(interf1) + '_' + str(interf2) + '/'    
    infile=folder + str(preffix) + str(interf1) + '-' + str(interf2) + str(suffix) + str(rlook) + 'rlks.unw'

    ds = gdal.Open(infile, gdal.GA_ReadOnly)
    ds_band2 = ds.GetRasterBand(2)
    ds_band1 = ds.GetRasterBand(1)
    print 'Nlign:{}, Ncol:{}, int:{}-{}'.format(ds.RasterYSize, ds.RasterXSize,interf1,interf2)

    los_map[:ds.RasterYSize,:ds.RasterXSize] = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:nlign,:ncol]
    # cor_map[:ds.RasterYSize,:ds.RasterXSize] = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:nlign,:ncol]

    alllos[:,:,i] = los_map
    # allcor[:,:,i] = cor_map


sumlos=np.zeros((nlign,ncol))
count=np.zeros((nlign,ncol))

# maxcoh = np.nanmax(allcor[:,:,:])
# allcor[:,:,:] = allcor[:,:,:]/maxcoh
# allcor[np.isnan(allcor)] = 0.0

for i in xrange(0,nlign,1):
    for j in xrange(0,ncol,1):
        k=np.flatnonzero(alllos[i,j,:])
        if len(k) > 1:
            sumlos[i,j]=np.sum(alllos[i,j,k[:]])/len(k)
            # sumlos[i,j]=np.sum(alllos[i,j,k[:]]*allcor[i,j,k[:]])/np.sum(allcor[i,j,k[:]])
               
    
folder=home + '/int/stacks/'
sumlos.astype('float32').tofile(folder + output)

fig = plt.figure(figsize=(12,5))

ax = fig.add_subplot(1,1,1)
ax.set_title('colorMap')
plt.imshow(sumlos,cmap='plasma',vmin=-10, vmax=10)
ax.set_aspect('equal')

cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)
plt.colorbar(orientation='vertical')
plt.show()