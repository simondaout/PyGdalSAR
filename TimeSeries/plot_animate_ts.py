#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Mathieu Volat 
############################################


"""\
plot_geots.py
-------------
Plot annimated time series maps 

Usage: plot_annimate_ts.py --cube=<path> [--vmin=<value>] [--vmax=<value>] \
        [--wrap=<values>] [--cpt=<values>] [--crop=<values>] [--output=<path>] 

Options:
-h --help           Show this screen.
--cube              time series displacement cube 
--vmax VALUE        Max colorscale 
--vmin VALUE        Min colorscale
--wrap  VALUE       Wrapped phase between value [default: no]
--cpt               Indicate colorscale
--crop VALUE        Crop option [default: 0,nlign,0,ncol]
--output            Optinional saving as mp4
"""


import sys
import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
import os

import docopt
arguments = docopt.docopt(__doc__)
if arguments["--cpt"] is  None:
    cmap=cm.rainbow
else:
    cmap=arguments["--cpt"]

# Open input dataset
ds = gdal.Open(arguments["--cube"])
if not ds:
    exit(1)
md = ds.GetMetadata()
print(md)

# Last band (should) give us a min/max displacement
band = ds.GetRasterBand(ds.RasterCount)
m = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
if arguments["--crop"] ==  None:
    crop = [0,ds.RasterYSize,0,ds.RasterXSize]
else:
    crop = map(float,arguments["--crop"].replace(',',' ').split())
ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])
vmax = np.nanpercentile(m[ibeg:iend,jbeg:jend], 95)
vmin = np.nanpercentile(m[ibeg:iend,jbeg:jend], 5)

if arguments["--vmax"] is not  None:
    vmax = float(arguments["--vmax"])

if arguments["--vmin"] is not  None:
    vmin = float(arguments["--vmin"])
    
if arguments["--wrap"] is not None:
    vmax=float(arguments["--wrap"])
    vmin=-vmax

# Create figure
fig = plt.figure(sys.argv[1])

# Figure update function, will be used at init too
def f(i):
    global ds
    i = i % ds.RasterCount + 1
    b = ds.GetRasterBand(i)
    plt.title(md["Band_%d"%i])
    if arguments["--wrap"] is not None:
      los = np.mod(b.ReadAsArray()[ibeg:iend,jbeg:jend]+float(arguments["--wrap"]),2*float(arguments["--wrap"]))-float(arguments["--wrap"])
    else:
      los = b.ReadAsArray()[ibeg:iend,jbeg:jend]
    return los

# Initialize
i = 0
im = plt.imshow(f(i), cmap=cmap, vmin=vmin, vmax=vmax, animated=True)

# Animation update function
def updatefig(*args):
    global i
    i = i+1
    im.set_array(f(i))
    return im,

# Animate
ani = animation.FuncAnimation(fig, updatefig, interval=200, blit=True)

if arguments["--output"] is not None:
    base = os.path.splitext(arguments["--output"])[0]
    #Writer = animation.writers['ffmpeg']
    #writer = Writer('fps=15',metadata=dict('Me'),bitrate=1800)
    ani.save(base+'.mp4')

# Display
plt.colorbar()
plt.show()

