#!/usr/bin/env python3

"""
compute_dphi_dz.py
-------------------
Reads in raster files - Velocities and topo
Plots topo on x axis vs  vel on y axis, pixel by pixel
Inputs must be the same size and the same number of pixels

.. Author:
    Roberta Wilkinson, University of Oxford; April 2023.
    Repurposes code from Simon Daout's plot_raster.py.
.. Last Modified:
    10th May 2023

Usage:
    vert-vs-topo.py --data=<path> --topo=<path> --outfile=<path> [--ymin=<value>] [--ymax=<value>] [--tp=<value>] [<ibeg>] [<iend>] [<jbeg>] [<jend>]
                    [--latmin=<val>] [--latmax=<val>] [--lonmin=<val>] [--lonmax=<val>] [--plot=<value>]

Options:
    -h --help             Show this screen.
    --data=PATH           Path to data velocity or LOS file.
    --topo=PATH           Path to topography file.
    --outfile=PATH        Output file name or path.
    --ymin=<value>        Minimum value for y-axis in scatter plot.
    --ymax=<value>        Maximum value for y-axis in scatter plot.
    --tp=<value>          Set transparency for scatter plot points [default: 0.1].
    --ibeg=<value>        Line number (smallest y-value) bounded the cutting zone [default: 0].
    --iend=<value>        Line number (biggest y-value) bounded the cutting zone [default: nlines].
    --jbeg=<value>        Column number (smallest x-value) bounded the cutting zone [default: 0].
    --jend=<value>        Column number (biggest x-value) bounded the cutting zone [default: ncols].
    --plot=<value>        if True, display figures  [default: True].
"""

import os, sys, logging
import docopt
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
from numpy.lib.stride_tricks import as_strided
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Rectangle

args = docopt.docopt(__doc__)

data = args["--data"]
topo = args["--topo"]
outfile = args["--outfile"]
tp = float(args["--tp"]) if args["--tp"] else 0.1
plot = args["--plot"]  if args["--plot"] else True

yflag = False
if args["--ymin"] and args["--ymax"]:
    ymin = float(args["--ymin"])
    ymax = float(args["--ymax"])
    yflag = True

sformat = 'GTIFF'
band = 1

ds = gdal.Open(data, gdal.GA_ReadOnly)
phase_band = ds.GetRasterBand(band)
nlines, ncols = ds.RasterYSize, ds.RasterXSize

dstopo = gdal.Open(topo, gdal.GA_ReadOnly)
phase_bandtopo = dstopo.GetRasterBand(band)

ibeg = int(args["<ibeg>"]) if args["<ibeg>"] else 0
iend = int(args["<iend>"]) if args["<iend>"] else nlines
jbeg = int(args["<jbeg>"]) if args["<jbeg>"] else 0
jend = int(args["<jend>"]) if args["<jend>"] else ncols

# Override with lat/lon if provided
if args["--latmin"] and args["--latmax"] and args["--lonmin"] and args["--lonmax"]:
    latmin, latmax = float(args["--latmin"]), float(args["--latmax"])
    lonmin, lonmax = float(args["--lonmin"]), float(args["--lonmax"])
    # Convert lat/lon to pixel coords
    def latlon_to_pixel(lat, lon, gt):
        px = int((lon - gt[0]) / gt[1])
        py = int((lat - gt[3]) / gt[5])
        return py, px
    ibeg, jbeg = latlon_to_pixel(latmax, lonmin, geotrans)
    iend, jend = latlon_to_pixel(latmin, lonmax, geotrans)

phi = phase_band.ReadAsArray()
phi[phi==0] = np.nan
cutphi = phi[ibeg:iend, jbeg:jend]
cutphi[cutphi==0] = np.nan

topo = phase_bandtopo.ReadAsArray()
topo = topo.astype(float)
topo[topo==0] = np.nan
cuttopo = topo[ibeg:iend, jbeg:jend]
cuttopo[cuttopo==0] = np.nan

xtopo = np.ravel(cuttopo)
yvel = np.ravel(cutphi)

# linear regression
index = ~np.isnan(yvel)
xtopop = xtopo[index]
yvelp = yvel[index]
G = np.zeros((len(xtopop), 2))
G[:, 0] = xtopop
G[:, 1] = 1
x = np.linalg.lstsq(G, yvelp, rcond=None)[0]
print('Linear regression %f x + %f' % (x[0], x[1]))
ysynth = np.dot(G, x)

fig2 = plt.figure(2, figsize=(12, 8))
axcompare = fig2.add_subplot(1, 1, 1)
axcompare.scatter(xtopo, yvel, s=0.5, alpha=tp, color='black')
axcompare.plot(xtopop, ysynth, 'r-', label='y =  {:.4f}x + {:.3f}'.format(x[0], x[1]))

if yflag:
    axcompare.set_ylim(ymin, ymax)

axcompare.set_xlabel('z (m)', fontsize=14)
axcompare.set_ylabel('Velocity (mm/yr)', fontsize=14)
axcompare.set_title(outfile, fontsize=14)

# Corrected legend
axcompare.legend(loc='best')  # Automatically find the best location
fig2.canvas.manager.set_window_title('Cropped data')
plt.tight_layout()
fig2.savefig(f'{outfile}_scatterplot.pdf', format='PDF', dpi=180)

#####################################################################################
# PLOT TIFF MAPS
#####################################################################################
try:
    from matplotlib.colors import LinearSegmentedColormap
    cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
    cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
    cmap = cmap.reversed()
    cmaptopo = cm.gray
except:
    cmap=cm.rainbow

# Initialize a matplotlib figure
fig = plt.figure(1,figsize=(15,8))

####### Plot the data velocity - Big context map
axcontext = fig.add_subplot(1,3,1)
# Put the values of the data velocity array into phi
phi = phase_band.ReadAsArray(0, 0,
       ds.RasterXSize, ds.RasterYSize,
       ds.RasterXSize, ds.RasterYSize)

vmaxcontext = np.nanpercentile(phi,98)
vmincontext = np.nanpercentile(phi,2)

# replace 0 by nan
try:
    phi[phi==0] = float('NaN')
except:
    pass
masked_arraycontext = np.ma.array(phi, mask=np.isnan(phi))

caxcontext = axcontext.imshow(masked_arraycontext, cmap, interpolation='nearest',vmax=vmaxcontext,vmin=vmincontext)
divider = make_axes_locatable(axcontext)
ccontext = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(caxcontext, cax=ccontext)

axcontext.set_title('Data', fontsize=14)

# Add the subplot outline
deli=abs(iend-ibeg)
delj=abs(jend-jbeg)
axcontext.add_patch(Rectangle((jbeg, ibeg), delj, deli, fc='none',color="black",linewidth=2))

####### Plot the data velocity - cropped map
ax = fig.add_subplot(1,3,2)
# Trim the maps
cutphi = as_strided(phi[ibeg:iend,jbeg:jend])

vmax = np.nanpercentile(cutphi,98)
vmin = np.nanpercentile(cutphi,2)

# replace 0 by nan
try:
    cutphi[cutphi==0] = float('NaN')
except:
    pass
masked_array = np.ma.array(cutphi, mask=np.isnan(cutphi))

cax = ax.imshow(masked_array, cmap, interpolation='nearest',vmax=vmax,vmin=vmin)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

ax.set_title('Data (cropped)', fontsize=14)

###### Plot the topography
axtopo = fig.add_subplot(1,3,3)

# Put the values of the data velocity array into phi
phitopo = phase_bandtopo.ReadAsArray(0, 0,
       dstopo.RasterXSize, dstopo.RasterYSize,
       dstopo.RasterXSize, dstopo.RasterYSize)

cutphitopo = as_strided(phitopo[ibeg:iend,jbeg:jend])

vmaxtopo = np.nanpercentile(cutphitopo,98)
vmintopo = np.nanpercentile(cutphitopo,2)

# replace 0 by nan
try:
    cutphitopo[cutphitopo==0] = float('NaN')
except:
    pass
masked_arraytopo = np.ma.array(cutphitopo, mask=np.isnan(cutphitopo))

caxtopo = axtopo.imshow(masked_arraytopo, cmaptopo, interpolation='nearest',vmax=vmaxtopo,vmin=vmintopo)
#cax = ax.imshow(masked_array, cmap, interpolation='none',vmax=vmax,vmin=vmin)
divider = make_axes_locatable(axtopo)
ctopo = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(caxtopo, cax=ctopo)

axtopo.set_title('Topography (cropped)', fontsize=14)

############################

fig.canvas.manager.set_window_title('Data and topo')

try:
    del ds
    del dstopo
except:
    pass

# Display the data
plt.tight_layout()
# ax.set_rasterized(True)
try:
    fig.savefig('{}_map.pdf'.format(outfile), format='PDF',dpi=180)
except:
    pass

plot=False
if plot:
    plt.show()
