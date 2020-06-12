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
plot_hist_dphi.py
-------------

Usage: plot_hist_dphi.py --file1=<path> --file2=<path> [--lectfile=<path>] [--rad2mm=<value>] [--vmin=<value>] [--vmax=<value>] 
plot_hist_dphi_r4.py -h | --help

Options:
-h --help           Show this screen.
--file1 PATH
--file2 PATH
--lectfile PATH       Path of the lect.in file [default: lect.in]
--rad2mm=<value>      Convert data [default: 1]
--vmax                Max colorscale [default: 98th percentile]
--vmin                Min colorscale [default: 2th percentile]
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import scipy.stats as stats
from scipy.optimize import curve_fit
import docopt,os,gdal

# read arguments
arguments = docopt.docopt(__doc__)
if arguments["--rad2mm"] ==  None:
        rad2mm = 1
else:
        rad2mm = np.float(arguments["--rad2mm"])
file1 = arguments["--file1"]
file2 = arguments["--file2"]
ds_extension = os.path.splitext(file1)[1]

if (ds_extension == ".r4" or ds_extension == ""):
    if arguments["--lectfile"] ==  None:
        lecfile = "lect.in"
    else:
        lecfile = arguments["--lectfile"]
    cols, lines = map(int, open(lecfile).readline().split(None, 2)[0:2])
    phi1 = np.fromfile(file1,dtype=np.float32)[:lines*cols].reshape((lines,cols))*rad2mm 
    phi2 = np.fromfile(file2,dtype=np.float32)[:lines*cols].reshape((lines,cols))*rad2mm 
#elif (ds_extension == ".tif" or ds_extension == ".tiff"):
else:
    ds1 = gdal.Open(file1, gdal.GA_ReadOnly)
    ds2 = gdal.Open(file2, gdal.GA_ReadOnly)
    phi1 = ds1.GetRasterBand(1).ReadAsArray(0, 0,
           ds1.RasterXSize, ds1.RasterYSize,
           ds1.RasterXSize, ds1.RasterYSize)
    phi2 = ds2.GetRasterBand(1).ReadAsArray(0, 0,
           ds2.RasterXSize, ds2.RasterYSize,
           ds2.RasterXSize, ds2.RasterYSize)

data = phi1 - phi2
data[data==0] = np.float('NaN')
data[phi1==0] = np.float('NaN')
data[phi2==0] = np.float('NaN')

# plot maps
if arguments["--vmax"] is not None:
    vmax = np.float(arguments["--vmax"])
else:
    vmax = np.max( [np.nanpercentile(data,98),np.abs(np.nanpercentile(data,2))] )
if arguments["--vmin"] is not None:
    vmin = np.float(arguments["--vmin"])
else:
    vmin = -vmax
try:
    from matplotlib.colors import LinearSegmentedColormap
    cm_locs = '/home/cometsoft/PyGdalSAR/contrib/python/colormaps/'
    cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
    cmap = cmap.reversed()
except:
    cmap=cm.jet

fig = plt.figure(figsize=(16,4))
fig.subplots_adjust(hspace=0.5)
ax1 = fig.add_subplot(1,3,1)
cax = ax1.imshow(phi1, cmap = cmap, vmax=vmax, vmin=vmin, extent=None,interpolation='nearest')
divider = make_axes_locatable(ax1)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
plt.setp( ax1.get_xticklabels(), visible=False)
ax1.set_title(file1)

ax2 = fig.add_subplot(1,3,2)
cax = ax2.imshow(phi2, cmap = cmap, vmax=vmax, vmin=vmin, extent=None,interpolation='nearest')
plt.setp( ax2.get_xticklabels(), visible=False)
divider = make_axes_locatable(ax2)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
ax2.set_title(file2)

#vmax = np.max( [np.nanpercentile(out,95),np.abs(np.nanpercentile(out,5))] )
ax3 = fig.add_subplot(1,3,3)
cax = ax3.imshow(data, cmap = cmap, vmax=vmax, vmin=vmin, extent=None,interpolation='nearest')
plt.setp( ax3.get_xticklabels(), visible=False)
divider = make_axes_locatable(ax3)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
fig.tight_layout()
fig.savefig('{}.pdf'.format('dphi'), format='PDF',dpi=150)

## plot histogram difference

# # Calculate mean and standard deviation
index = np.nonzero(
            np.logical_and(~np.isnan(data),
            np.logical_and(data>np.nanpercentile(data, 2.), data<np.nanpercentile(data, 98.)
            )))

diff = data[index]
diff_med = np.nanmedian(diff)
diff_std = np.nanstd(diff)

opts = {'c':'red', 'linestyle':'--'}
fig, ax1  = plt.subplots(1, 1, figsize=(6,4))
sns.distplot(diff, norm_hist=True, hist=True, color="dodgerblue", ax=ax1)
#ax1.set_xlim([vmin,vmax])
ax1.axvline(x=diff_med,alpha=0.4, **opts)
ax1.set_xlabel('Median: {:.3f} STD: {:.3f}'.format(diff_med, diff_std))
ax1.set_ylabel('Norm. PDF')
ax1.set_ylim(bottom=0)

fig.savefig('{}.pdf'.format('hist_dphi'), format='PDF',dpi=150)

plt.show()
