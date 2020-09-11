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
compute_stack.py
-------------
Compute a stack from a list of interferograms

Usage: compute_stack.py --homedir=<path> --int_path=<path> --int_list=<file> [--suffix=<value>] [--prefix=<value>] [--output=<file>] [--rlook=<value>] --radar_file=<file>
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
--output
"""

import numpy as np
import matplotlib.pyplot as plt
# gdal
import gdal,os
gdal.UseExceptions()
import docopt

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

# list of dates
date1,date2=np.loadtxt(int_list,comments="#",unpack=True,dtype='i,i')
kmax = len(date1)

# radar file
driver = gdal.GetDriverByName("roi_pac")
ds = gdal.Open(radar, gdal.GA_ReadOnly)
nlines,ncols= ds.RasterYSize, ds.RasterXSize

alllos = np.zeros((nlines,ncols,kmax))
allcor = np.zeros((nlines,ncols,kmax))

for i in xrange((kmax)):

    interf1,interf2 = date1[i],date2[i]
    los_map = np.zeros((nlines,ncols))
    cor_map = np.zeros((nlines,ncols))

    folder= int_path  + 'int_' + str(interf1) + '_' + str(interf2) + '/'    
    infile=folder + str(prefix) + str(interf1) + '-' + str(interf2) + str(suffix) + str(rlook) + 'rlks.unw'

    ds = gdal.Open(infile, gdal.GA_ReadOnly)
    ds_band2 = ds.GetRasterBand(2)
    ds_band1 = ds.GetRasterBand(1)
    print 'Nlign:{}, Ncol:{}, int:{}-{}'.format(ds.RasterYSize, ds.RasterXSize,interf1,interf2)

    los_map[:ds.RasterYSize,:ds.RasterXSize] = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:nlines,:ncols]
    # cor_map[:ds.RasterYSize,:ds.RasterXSize] = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:nlines,:ncols]

    alllos[:,:,i] = los_map
    # allcor[:,:,i] = cor_map


sumlos=np.zeros((nlines,ncols))
count=np.zeros((nlines,ncols))

# maxcoh = np.nanmax(allcor[:,:,:])
# allcor[:,:,:] = allcor[:,:,:]/maxcoh
# allcor[np.isnan(allcor)] = 0.0

for i in xrange(0,nlines,1):
    for j in xrange(0,ncols,1):
        k=np.flatnonzero(alllos[i,j,:])
        if len(k) > 1:
            sumlos[i,j]=np.sum(alllos[i,j,k[:]])/len(k)
            # sumlos[i,j]=np.sum(alllos[i,j,k[:]]*allcor[i,j,k[:]])/np.sum(allcor[i,j,k[:]])
               
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
plt.imshow(sumlos,cmap=cmap)
ax.set_aspect('equal')

cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)
plt.colorbar(orientation='vertical')
plt.show()
