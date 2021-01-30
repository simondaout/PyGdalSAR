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


import gdal,sys
import pymc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import *
import matplotlib.ticker as ticker

ds_amp = gdal.Open("NETCDF:"+sys.argv[1]+':Band1', gdal.GA_ReadOnly)
ds_phi = gdal.Open("NETCDF:"+sys.argv[1]+':Band2', gdal.GA_ReadOnly)

amp = ds_amp.ReadAsArray(0, 0, ds_amp.RasterXSize, ds_amp.RasterYSize)
phi = ds_phi.ReadAsArray(0, 0, ds_phi.RasterXSize, ds_phi.RasterYSize)

def hdi(trace, cred_mass=0.95):
    hdi_min, hdi_max = pymc.utils.calc_min_interval(np.sort(trace), 1.0-cred_mass)
    return hdi_min, hdi_max

# wrapped phase
phi = np.fmod(phi,2*np.pi).flatten()
kk = np.nonzero(phi<0)
phi[kk] += 2*np.pi

#phimin,phimax = hdi(phi)

fig = plt.figure(0, figsize=(9,6))
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlim([0,2*np.pi])
histo = ax1.hist(phi,normed=True,range=(0,2*np.pi),bins=100,histtype='step',color='black')
plt.axvline(x=np.median(phi), c="red")
plt.xlabel('Median phase shift: {:0.3f} (rad) '.format(np.median(phi)))
start,end = ax1.get_xlim()
ax1.set_xticks(np.arange(start,end , (end-start)/12))
ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.legend(loc='best')
fig.savefig('histo_phitemp.eps', format='EPS')
plt.show()
