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
plot_hist_dphi_r4.py
-------------

Usage: plot_hist_dphi_r4.py --infile=<path> --demfile=<path> [--lectfile=<path>] [--rad2mm=<value>] [--vmin=<value>] [--vmax=<value>] 
plot_hist_dphi_r4.py -h | --help

Options:
-h --help           Show this screen.
--infile PATH      
--demfile PATH
--lectfile PATH     Path of the lect.in file [default: lect.in]
--rad2mm=<value>      Convert data [default: 1]
--vmax                Max colorscale [default: 98th percentile]
--vmin                Min colorscale [default: 2th percentile]
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
from scipy.optimize import curve_fit
import docopt

# read arguments
arguments = docopt.docopt(__doc__)
infile = arguments["--infile"]
demfile = arguments["--demfile"]
if arguments["--lectfile"] ==  None:
    lecfile = "lect.in"
else:
    lecfile = arguments["--lectfile"]
# read lect.in
cols, lines = map(int, open(lecfile).readline().split(None, 2)[0:2])
if arguments["--rad2mm"] ==  None:
        rad2mm = 1
else:
        rad2mm = np.float(arguments["--rad2mm"])

# read r4
data = np.fromfile(infile,dtype=np.float32)[:lines*cols].reshape((lines,cols))*rad2mm 
dem = np.fromfile(demfile,dtype=np.float32)[:lines*cols].reshape((lines,cols))
data[data==0] = np.float('NaN')
dem[np.isnan(data)] = np.float('NaN')

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
    topomap=LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"batlow.txt"))
except:
    cmap=cm.jet
    topomap=cm.terrain

# plot hist and phase topo diff
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))

cax = ax1.imshow(data, cmap = cmap, vmax=vmax, vmin=vmin, extent=None)
divider = make_axes_locatable(ax1)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
plt.setp( ax1.get_xticklabels(), visible=False)
ax1.set_title(infile)

cax = ax2.imshow(dem, cmap = topomap, extent=None)
plt.setp( ax2.get_xticklabels(), visible=False)
ax2.set_title(demfile)
divider = make_axes_locatable(ax2)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

# # Calculate mean and standard deviation
index = np.nonzero(
            np.logical_and(~np.isnan(data),
            np.logical_and(data>np.nanpercentile(data, 2.), data<np.nanpercentile(data, 98.)
            )))

diff = data[index]
diff_mean = np.nanmedian(diff)
diff_std = np.nanstd(diff)
dem_clean = dem[index]

opts = {'c':'red', 'linestyle':'--'}
# plot hist and phase topo diff
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
sns.distplot(diff, norm_hist=True, hist=True, color="dodgerblue", ax=ax1)
#ax1.set_xlim([vmin,vmax])
ax1.axvline(x=diff_mean,alpha=0.4, **opts)
ax1.set_xlabel('Mean: {:.3f} STD: {:.3f}'.format(diff_mean, diff_std))
ax1.set_ylabel('Norm. PDF')
ax1.set_ylim(bottom=0)

f = ax1.lines[0]
xf = f.get_xydata()[:,0]
yf = f.get_xydata()[:,1]

ax1.fill_between(xf, yf, color="dodgerblue", alpha=0.5, where=(xf>(diff_mean-1*diff_std)) & (xf<(diff_mean+1*diff_std)))

def linear_f(x, a, b):
    return a*x + b

interval = 20
ax2.scatter(dem_clean[::interval], diff[::interval], color='dodgerblue', alpha=0.05,
        marker = 's', s = 5, edgecolor = 'none',rasterized=True,
        )

popt, pcov = curve_fit(linear_f, dem_clean, diff)
ax2.plot(dem_clean, linear_f(dem_clean,*popt), color='black', lw=2, 
        label="{:.6f} x + {:.2f}".format(*popt))
ax2.set_ylabel('$\phi$')
ax2.set_xlabel('Elevation')
ax2.set_xlim([np.percentile(dem_clean,2),np.percentile(dem_clean,98)])
plt.legend(loc='best')

fig.savefig('{}.pdf'.format('hist_dphi'), format='PDF',dpi=150)

plt.show()
