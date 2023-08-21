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
        rad2mm = float(arguments["--rad2mm"])
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
    cols, lines = ds2.RasterXSize, ds2.RasterYSize

data = phi1 - phi2
data[data==0] = float('NaN')
data[phi1==0] = float('NaN')
data[phi2==0] = float('NaN')

# plot maps
if arguments["--vmax"] is not None:
    vmax = float(arguments["--vmax"])
else:
    vmax = np.max( [np.nanpercentile(data,98),np.abs(np.nanpercentile(data,2))] )
if arguments["--vmin"] is not None:
    vmin = float(arguments["--vmin"])
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
ax1 = fig.add_subplot(2,2,1)
cax = ax1.imshow(phi1, cmap = cmap, vmax=vmax, vmin=vmin, extent=None,interpolation='nearest')
divider = make_axes_locatable(ax1)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
plt.setp( ax1.get_xticklabels(), visible=False)
ax1.set_title(file1)

ax2 = fig.add_subplot(2,2,2)
cax = ax2.imshow(phi2, cmap = cmap, vmax=vmax, vmin=vmin, extent=None,interpolation='nearest')
plt.setp( ax2.get_xticklabels(), visible=False)
divider = make_axes_locatable(ax2)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
ax2.set_title(file2)

#vmax = np.max( [np.nanpercentile(out,95),np.abs(np.nanpercentile(out,5))] )
ax3 = fig.add_subplot(2,2,3)
cax = ax3.imshow(data, cmap = cmap, vmax=vmax, vmin=vmin, extent=None,interpolation='nearest')
plt.setp( ax3.get_xticklabels(), visible=False)
divider = make_axes_locatable(ax3)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
ax3.set_title('Diff')

## remove ramp
ramp = 'yes'
if ramp == 'yes':
    import scipy.optimize as opt
    import scipy.linalg as lst

    index = np.nonzero(~np.isnan(data))
    temp = np.array(index).T
    mi = data[index].flatten()
    az = temp[:,0]; rg = temp[:,1]

    G=np.zeros((len(mi),3))
    G[:,0] = rg
    G[:,1] = az
    G[:,2] = 1

    x0 = lst.lstsq(G,mi)[0]
    _func = lambda x: np.sum(((np.dot(G,x)-mi))**2)
    _fprime = lambda x: 2*np.dot(G.T, (np.dot(G,x)-mi))
    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0)[0]
    a = pars[0]; b = pars[1]; c = pars[2]
    print('Remove ramp %f x  + %f y + %f'.format(a,b,c))
    
    G=np.zeros((len(data.flatten()),3))
    for i in range(lines):
        G[i*cols:(i+1)*cols,0] = np.arange((cols))
        G[i*cols:(i+1)*cols,1] = i
    G[:,2] = 1
    temp = (data.flatten() - np.dot(G,pars))
    data=temp.reshape(lines,cols)

#vmax = np.max( [np.nanpercentile(out,95),np.abs(np.nanpercentile(out,5))] )
ax4 = fig.add_subplot(2,2,4)
cax = ax4.imshow(data, cmap = cmap, vmax=vmax, vmin=vmin, extent=None,interpolation='nearest')
plt.setp( ax4.get_xticklabels(), visible=False)
divider = make_axes_locatable(ax4)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
ax4.set_title('Diff-Ramp')
fig.tight_layout()

fig.savefig('{}.pdf'.format('dphi'), format='PDF',dpi=150)

# save output 
if (ds_extension == ".r4" or ds_extension == ""):
    data.flatten().astype('float32').tofile('diff.r4')
else:
    outfile = 'diff.tiff'
    ds_geo = ds1.GetGeoTransform()
    proj = ds1.GetProjection()
    driver = gdal.GetDriverByName('GTiff')
    print('Convert output file:', outfile)
    dst_ds = driver.Create(outfile, ds1.RasterXSize, ds1.RasterYSize, 1, gdal.GDT_Float32)
    dst_band2 = dst_ds.GetRasterBand(1)
    dst_band2.WriteArray(data,0,0)
    dst_ds.SetGeoTransform(ds_geo)
    dst_ds.SetProjection(proj)
    dst_band2.FlushCache()
    del dst_ds

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
