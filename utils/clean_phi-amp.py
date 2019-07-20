#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Simon Daout
############################################

"""\
clean_phi-amp.py
-------------
Clean amplitude and phase maps obtained with the  seasonal decomposition
and wrapped phase between 0 to 2pi

Usage: clean_phi-amp.py [--ampfile=<path>] [--phifile=<path>] [--linfile=<path>] [--topofile=<path>] \
[--sigampfile=<path>] [--sigphifile=<path>] [--lectfile=<path>] [--threshold_amp=<value>] [--perc_sig=<value>] \
[--outampfile=<path>] [--outphifile=<path>] 


Options:
-h --help             Show this screen.
--ampfile=<file>      Amplitude map file [default: ampwt_coeff.r4]
--phifile=<file>      Phase map file [default: phiwt_coeff.r4]
--linfile=<file>      Linear map file [default: lin_coeff.r4]
--topofile=<file>     Linear map file [default: dem]
--sigampfile=<file>   Uncertainty amplitude map file [default: ampwt_sigcoeff.r4]
--sigphifile=<file>   Uncertainty phase map file [default: phiwt_sigcoeff.r4]
--lectfile=<file>     Path of the lect.in file for r4 format [default: lect_ts.in]
--theshold_amp=<value>      Mask on minimum Amplitude for Phase [default: 1.5]
--perc_sig=<value>            Percentile uncertainty for map cleaning [default: 99.]
"""

import numpy as np
from numpy.lib.stride_tricks import as_strided
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from pylab import *
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
np.warnings.filterwarnings('ignore')
try:
    from nsbas import docopt
except:
    import docopt

arguments = docopt.docopt(__doc__)
if arguments["--ampfile"] ==  None:
    ampf = 'ampwt_coeff.r4'
else:
    ampf = arguments["--ampfile"]
if arguments["--phifile"] ==  None:
    phif = 'phiwt_coeff.r4'
else:
    phif = arguments["--phifile"]
if arguments["--linfile"] ==  None:
    linf = 'lin_coeff.r4'
else:
    linf = arguments["--linfile"]
if arguments["--sigampfile"] ==  None:
    ampsigf = 'ampwt_sigcoeff.r4'
else:
    ampsigf = arguments["--sigampfile"]
if arguments["--sigphifile"] ==  None:
    phisigf = 'phiwt_sigcoeff.r4'
else:
    phisigf = arguments["--sigphifile"]
if arguments["--lectfile"] ==  None:
    lectf = "lect_ts.in"
else:
    lectf = arguments["--lectfile"]
if arguments["--outampfile"] ==  None:
    ampoutf = 'ampwt_coeff_clean.r4'
else:
    ampoutf = arguments["--outampfile"]
if arguments["--outphifile"] ==  None:
    phioutf = 'phiwt_coeff_clean.r4'
else:
    phioutf = arguments["--outphifile"]

# disp. maps
amp,phi,lin=np.fromfile(ampf,dtype=np.float32),np.fromfile(phif,dtype=np.float32),np.fromfile(linf,dtype=np.float32)
# uncertainties maps
sigamp,sigphi=np.fromfile(ampsigf,dtype=np.float32),np.fromfile(phisigf,dtype=np.float32)
# convrt to map
ncols, nlines = map(int, open(lectf).readline().split(None, 2)[0:2])
amp_map, phi_map, lin_map = amp.reshape(nlines,ncols),phi.reshape(nlines,ncols),lin.reshape(nlines,ncols)
sigamp_map, sigphi_map = sigamp.reshape(nlines,ncols),sigphi.reshape(nlines,ncols)

if arguments["--topofile"] ==  None:
    try:
        dem = np.fromfile('dem',dtype=np.float32)[:nlines*ncols].reshape(nlines,ncols)
    except:
        demf = None
        dem = False
else:
    demf = arguments["--topofile"]
    dem = np.fromfile(demf,dtype=np.float32)[:nlines*ncols].reshape(nlines,ncols)

# # clean
if arguments["--perc_sig"] ==  None:
    perc_sig = 99.
else:
    perc_sig = np.float(arguments["--perc_sig"])
phi_map[sigamp_map>np.nanpercentile(sigamp_map,perc_sig)] = np.float('NaN')
amp_map[sigamp_map>np.nanpercentile(sigamp_map,perc_sig)] = np.float('NaN')
lin_map[sigamp_map>np.nanpercentile(sigamp_map,perc_sig)] = np.float('NaN')

phi_map[sigphi_map>np.nanpercentile(sigphi_map,perc_sig)] = np.float('NaN')
amp_map[sigphi_map>np.nanpercentile(sigphi_map,perc_sig)] = np.float('NaN')
lin_map[sigphi_map>np.nanpercentile(sigphi_map,perc_sig)] = np.float('NaN')

if arguments["--threshold_amp"] ==  None:
    threshold_amp = 1.5
else:
    threshold_amp = np.float(arguments["--threshold_amp"])
phi_map[amp_map<threshold_amp] = np.float('NaN')
lin_map[amp_map<threshold_amp] = np.float('NaN')

# convert phi between 0 and 2pi
phi[phi<0] = phi[phi<0] + 2*np.pi

# plot
fig_ampphi = plt.figure(0,figsize=(12,6))

#save output
fid1 = open(ampoutf,'wb')
amp_map.flatten().astype('float32').tofile(fid1)
fid1.close()
fid2 = open(phioutf,'wb')
phi_map.flatten().astype('float32').tofile(fid2)
fid2.close()

# initiate figure amp and phase
# rad2mm=-4.456

ax = fig_ampphi.add_subplot(1,3,1)
cmap = cm.rainbow
if dem is not False:
    cax = ax.imshow(dem,cmap=cm.Greys,zorder=1)
# no point to take negatif amplitude
im = ax.imshow(np.ma.array(amp_map, mask=np.isnan(amp_map)),cmap=cmap,vmin=threshold_amp,\
    vmax=np.nanpercentile(amp_map,98), alpha=0.8, zorder=2)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
ax.set_title('Amplitude (mm)',fontsize=12)

ax = fig_ampphi.add_subplot(1,3,2)
cmap = cm.tab20b
if dem is not False:
    cax = ax.imshow(dem,cmap=cm.Greys,zorder=1)
im = ax.imshow(np.ma.array(phi_map, mask=np.isnan(phi_map)),cmap=cmap,vmin=0,vmax=2*np.pi,alpha=0.8,zorder=2)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
ax.set_title('Timing (rad)',fontsize=12)

ax = fig_ampphi.add_subplot(1,3,3)
cpt = '/home/cometsoft/PyGdalSAR/contrib/python/colormaps/roma.txt'
cmap = LinearSegmentedColormap.from_list(cpt.split("/")[-1].split('.')[0], np.loadtxt(cpt))
if dem is not False:
    cax = ax.imshow(dem,cmap=cm.Greys,zorder=1)
vmax = np.nanpercentile(lin_map,95)
im = ax.imshow(np.ma.array(lin_map, mask=np.isnan(lin_map)),vmin=-vmax,vmax=vmax, cmap=cmap, alpha=0.8,zorder=2)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
ax.set_title('Velocity (mm/yr)',fontsize=12)

fig_ampphi.tight_layout()
plt.show()