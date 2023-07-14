#!/usr/bin/env python3

'''
nsb_atmos_corr_histo.py
-------------
Histogram plotting of atmospheric correction parameters saved in IFG = ATMOS + RAMP estimation.

.. Author: Nicholas Dodds, University of Oxford. February, 2020.

Usage:
    nsb_gacos2rdr.py --atmos_corrs_grads=<path>
    
Options:
    --atmos_corrs_grads= path    Path to output file from correction application. Needs correlations and gradients for each IFG.

'''

print()
print()
print('Authors: Nicholas Dodds, Simon DAOUT')
print('Please cite:')
print('Dodds, N., Daout, S., Walker, R. T., Begenjev, G., Bezmenov, Y., Mirzin, R., & Parsons, B. (2022). Interseismic deformation and strain-partitioning along the Main Köpetdag Fault, Turkmenistan, with Sentinel-1 InSAR time-series. Geophysical Journal International, 230(3), 1612-1629.')
print()
print()


import os
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
import docopt

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import matplotlib.ticker as plticker
from matplotlib.colorbar import Colorbar
import matplotlib.patches as patches
import matplotlib.ticker as ticker

###

# read input parameters and arguments
arguments = docopt.docopt(__doc__)

if arguments["--atmos_corrs_grads"] != None:
    atmos_file = arguments["--atmos_corrs_grads"]

# read in data from text file
atmos_data = np.loadtxt(atmos_file, delimiter=',')
atmos_grads = np.asarray([i[2] for i in atmos_data])
atmos_corrs = np.asarray([i[3] for i in atmos_data])

# Calculate mean and standard deviation
grads_mean = np.mean(atmos_grads)
grads_std = np.std(atmos_grads)
corrs_mean = np.mean(atmos_corrs)
corrs_std = np.std(atmos_corrs)

###
# Plotting
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16,4))

##
sns.distplot(atmos_grads, norm_hist=True, hist=True, color="dodgerblue", ax=ax1)
ax1.set_xlabel('Gradient: $\phi_{ATMOS}$ vs. $(\phi_{IFG} - RAMP)$')
ax1.set_ylabel('Norm. PDF')
ax1.set_ylim(bottom=0)

f = ax1.lines[0]
xf = f.get_xydata()[:,0]
yf = f.get_xydata()[:,1]

ax1.fill_between(xf, yf, color="dodgerblue", alpha=0.5, where=(xf>(grads_mean-1*grads_std)) & (xf<(grads_mean+1*grads_std)))
ax1.axvline(xf[np.argmax(yf)], color="dodgerblue",linestyle='dashed')

##
sns.distplot(atmos_corrs, norm_hist=True, hist=True, color="lightcoral", ax=ax2)
ax2.set_xlabel('Correlation: $\phi_{IFG}$ vs. $(\phi_{ATMOS} - RAMP)$')
ax2.set_ylabel('Norm. PDF')
ax2.set_ylim(bottom=0)
ax2.set_xlim(left=0)

f = ax2.lines[0]
xf = f.get_xydata()[:,0]
yf = f.get_xydata()[:,1]

ax2.fill_between(xf, yf, color="lightcoral", alpha=0.5, where=(xf>(corrs_mean-1*corrs_std)) & (xf<(corrs_mean+1*corrs_std)))
ax2.axvline(xf[np.argmax(yf)], color="lightcoral",linestyle='dashed')

## 
sns.scatterplot(atmos_grads, atmos_corrs, color='black', ax=ax3)
ax3.set_xlabel('Gradient')
ax3.set_ylabel('Correlation')
ax3.set_ylim(bottom=0)

fig.savefig('phase-atmos_histo.png', format='PNG', bbox_inches='tight')
plt.clf()
