#!/usr/bin/env python2.7
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
[--outampfile=<path>] [--outphifile=<path>] [--outlinfile=<path>] [--slopefile=<path>] [--slopelosfile=<path>] [--aspectfile=<path>]  [--plotcorr=<yes/no>] \
[--maxamp=<value>] [--minelev=<value>] [--rad2mm=<value>] [--outfig=<path>] [--crop=<values>]


Options:
-h --help               Show this screen.
--ampfile=<file>        Amplitude map file [default: ampwt_coeff.r4]
--phifile=<file>        Phase map file [default: phiwt_coeff.r4]
--linfile=<file>        Linear map file [default: lin_coeff.r4]
--topofile=<file>       Linear map file [default: dem]
--sigampfile=<file>     Uncertainty amplitude map file [default: ampwt_sigcoeff.r4]
--sigphifile=<file>     Uncertainty phase map file [default: phiwt_sigcoeff.r4]
--lectfile=<file>       Path of the lect.in file for r4 format [default: lect_ts.in]
--threshold_amp=<value>  Mask on minimum Amplitude for Phase (rad) [default: 1.5]
--perc_sig=<value>      Percentile uncertainty for map cleaning [default: 100]
--slopelosfile=<file>   SLope in the LOS direction [default: None]
--slopefile=<file>      SLope file [default: None]
--aspectfile=<file>     Aspect file [default: None]
--plotcorr=<yes/no>     Plot correlation plots [default: no]
--maxamp=<value>        Maximum Amplitude limit [default: 3.]
--minelev=<value>       Minimum Elevation limit [default: 0.]
--rad2mm                Scaling value between input data (rad) and desired output [default: -4.4563]
--outfig                Name for the output figure [default: clean_phi-amp]
--crop=<values>       Crop data [default: 0,nlines,0,ncol]
"""

from __future__ import print_function
import numpy as np
from numpy.lib.stride_tricks import as_strided
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from pylab import *
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import scipy.ndimage

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
if arguments["--outlinfile"] ==  None:
    linoutf = 'lin_coeff_clean.r4'
else:
    linoutf = arguments["--outlinfile"]
if arguments["--topofile"] == None:
    demf = None
else:
    demf = arguments["--topofile"]
if arguments["--plotcorr"] == None:
    plotcorr = 'no'
else:
    plotcorr = arguments["--plotcorr"]
if arguments["--maxamp"] ==  None:
    maxamp = 3
else:
    maxamp = np.float(arguments["--maxamp"])
if arguments["--minelev"] ==  None:
    minelev = 0
else:
    minelev = np.float(arguments["--minelev"])
if arguments["--rad2mm"] ==  None:
    rad2mm = -4.4563
else:
    rad2mm = float(arguments["--rad2mm"]) 
if arguments["--threshold_amp"] ==  None:
    threshold_amp = 1.5
else:
    threshold_amp = np.float(arguments["--threshold_amp"])

if arguments["--outfig"] ==  None:
    outfig = 'clean_phi-amp'
else:
    outfig = arguments["--outfig"]

ncols, nlines = map(int, open(lectf).readline().split(None, 2)[0:2])
if arguments["--crop"] ==  None:
    crop = [0,nlines,0,ncols]
else:
    crop = list(map(float,arguments["--crop"].replace(',',' ').split()))
ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])

amp,phi,lin=np.fromfile(ampf,dtype=np.float32)[:nlines*ncols],np.fromfile(phif,dtype=np.float32)[:nlines*ncols],np.fromfile(linf,dtype=np.float32)[:nlines*ncols]
sigamp,sigphi=np.fromfile(ampsigf,dtype=np.float32)[:nlines*ncols],np.fromfile(phisigf,dtype=np.float32)[:nlines*ncols]
amp_map, phi_map, lin_map = amp.reshape(nlines,ncols)[ibeg:iend,jbeg:jend],phi.reshape(nlines,ncols)[ibeg:iend,jbeg:jend],lin.reshape(nlines,ncols)[ibeg:iend,jbeg:jend]
sigamp_map, sigphi_map = sigamp.reshape(nlines,ncols)[ibeg:iend,jbeg:jend],sigphi.reshape(nlines,ncols)[ibeg:iend,jbeg:jend]
del amp,phi,lin

if demf ==  None:
    try:
        demf = 'dem'
        dem_map = np.fromfile('dem',dtype=np.float32)[:nlines*ncols].reshape(nlines,ncols)[ibeg:iend,jbeg:jend]
    except:
        print('DEM file is not readible. Set elelvation to zeros.')
        dem_map = np.zeros((nlines,ncols))[ibeg:iend,jbeg:jend]
else:
    dem_map = np.fromfile(demf,dtype=np.float32)[:nlines*ncols].reshape(nlines,ncols)

if arguments["--slopefile"] is not None:
    slope_map = np.fromfile(arguments["--slopefile"],dtype=np.float32)[:nlines*ncols].reshape(nlines,ncols)[ibeg:iend,jbeg:jend]*100
else:
    toposmooth = scipy.ndimage.filters.gaussian_filter(dem_map,3.)
    Py, Px = np.gradient(toposmooth)
    slope_map = np.sqrt(Px**2+Py**2)

if arguments["--slopelosfile"] is not None:
    slopelos_map = np.fromfile(arguments["--slopelosfile"],dtype=np.float32)[:nlines*ncols].reshape(nlines,ncols)[ibeg:iend,jbeg:jend]*100
    slopelos_map[slope_map<3.] = np.float('NaN') # mask aspect for small slope
else:
    slopelos_map = np.zeros((nlines,ncols))[ibeg:iend,jbeg:jend]

if arguments["--aspectfile"] is not None:
    aspect_map = np.fromfile(arguments["--aspectfile"],dtype=np.float32)[:nlines*ncols].reshape(nlines,ncols)[ibeg:iend,jbeg:jend]
    aspect_map[slope_map<3] = np.float('NaN') # mask aspect for small slope
    aspect_map[aspect_map<0] = aspect_map[aspect_map<0] + 360
else:
    aspect_map = np.zeros((nlines,ncols))[ibeg:iend,jbeg:jend]

# # clean
if arguments["--perc_sig"] ==  None:
    perc_sig = 100
else:
    perc_sig = np.float(arguments["--perc_sig"])
    
phi_map[sigamp_map>np.nanpercentile(sigamp_map,perc_sig)] = np.float('NaN')
amp_map[sigamp_map>np.nanpercentile(sigamp_map,perc_sig)] = np.float('NaN')
lin_map[sigamp_map>np.nanpercentile(sigamp_map,perc_sig)] = np.float('NaN')

phi_map[sigphi_map>np.nanpercentile(sigphi_map,perc_sig)] = np.float('NaN')
amp_map[sigphi_map>np.nanpercentile(sigphi_map,perc_sig)] = np.float('NaN')
lin_map[sigphi_map>np.nanpercentile(sigphi_map,perc_sig)] = np.float('NaN')

# clean based on elev
phi_map[dem_map<minelev] = np.float('NaN')
amp_map[dem_map<minelev] = np.float('NaN')
lin_map[dem_map<minelev] = np.float('NaN')

phi_map[amp_map<threshold_amp] = np.float('NaN')
phi_map[amp_map>maxamp] = np.float('NaN')
# lin_map[amp_map<threshold_amp] = np.float('NaN')

# convert phi between 0 and 2pi
phi_map[phi_map<0] = phi_map[phi_map<0] + 2*np.pi

#save output
fid1 = open(ampoutf,'wb')
amp_map.flatten().astype('float32').tofile(fid1)
fid1.close()
fid2 = open(phioutf,'wb')
phi_map.flatten().astype('float32').tofile(fid2)
fid2.close()
fid3 = open(linoutf,'wb')
lin_map.flatten().astype('float32').tofile(fid3)
fid3.close()

# conversion
amp_map,sigamp_map,sigphi_map = amp_map*abs(rad2mm),sigamp_map*abs(rad2mm),sigphi_map*abs(rad2mm)
lin_map = lin_map*rad2mm
threshold_amp,maxamp = threshold_amp*abs(rad2mm),maxamp*abs(rad2mm)

# initiate figure amp and phase
fig_ampphi = plt.figure(0,figsize=(12,6))
ax = fig_ampphi.add_subplot(2,3,1)
cmap = cm.rainbow
# if arguments["--slopefile"] is not None:
cax = ax.imshow(slope_map,vmax=np.nanpercentile(slope_map,98),vmin=np.nanpercentile(slope_map,2),cmap=cm.Greys,zorder=1)
# no point to take negatif amplitude
im = ax.imshow(np.ma.array(amp_map, mask=np.isnan(amp_map)),cmap=cmap,vmin=threshold_amp,\
    vmax=np.nanpercentile(amp_map,98), alpha=0.8, zorder=2)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
ax.set_title('Amplitude (mm)',fontsize=12)

ax = fig_ampphi.add_subplot(2,3,2)
# cmap_phi = cm.tab20b
cmap_phi = cm.PiYG_r
# if arguments["--slopefile"] is not None:
cax = ax.imshow(slope_map,vmax=np.nanpercentile(slope_map,98),vmin=np.nanpercentile(slope_map,2),cmap=cm.Greys,zorder=1)
im = ax.imshow(np.ma.array(phi_map, mask=np.isnan(phi_map)),cmap=cmap_phi,vmin=0,vmax=2*np.pi,alpha=0.8,zorder=2)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
ax.set_title('Timing (rad)',fontsize=12)

ax = fig_ampphi.add_subplot(2,3,3)
# if arguments["--slopefile"] is not None:
cax = ax.imshow(slope_map,vmax=np.nanpercentile(slope_map,98),vmin=np.nanpercentile(slope_map,2),cmap=cm.Greys,zorder=1)
vmax = np.nanpercentile(lin_map,95)
vmin = np.nanpercentile(lin_map,5)
im = ax.imshow(np.ma.array(lin_map, mask=np.isnan(lin_map)),vmin=vmin,vmax=vmax, cmap=cmap, alpha=0.8,zorder=2)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
ax.set_title('Velocity (mm/yr)',fontsize=12)
fig_ampphi.tight_layout()

ax = fig_ampphi.add_subplot(2,3,4)
# if arguments["--slopefile"] is not None:
im = ax.imshow(slope_map,vmax=np.nanpercentile(slope_map,98),vmin=np.nanpercentile(slope_map,2),cmap=cm.Greys,zorder=1)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
ax.set_title('Slope',fontsize=12)

ax = fig_ampphi.add_subplot(2,3,5)
# if arguments["--slopefile"] is not None:
im = ax.imshow(slopelos_map,vmax=np.nanpercentile(slopelos_map,90),vmin=np.nanpercentile(slopelos_map,10),cmap=cm.Greys,zorder=1)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
ax.set_title('Slope LOS',fontsize=12)

ax = fig_ampphi.add_subplot(2,3,6)
# if arguments["--slopefile"] is not None:
im = ax.imshow(aspect_map,vmax=np.nanpercentile(aspect_map,98),vmin=np.nanpercentile(aspect_map,2),cmap=cm.hsv,zorder=1)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
ax.set_title('Aspect',fontsize=12)

fig_ampphi.tight_layout()
fig_ampphi.savefig('{}.pdf'.format(outfig), format='PDF',dpi=150)
# 

# remove Nan 
kk = np.nonzero(~np.isnan(phi_map))
phi,amp,lin,dem,slopelos,slope, aspect = phi_map[kk], amp_map[kk], lin_map[kk], dem_map[kk], slopelos_map[kk], slope_map[kk], aspect_map[kk]

# figure histo phi
fig_histo = plt.figure(1,figsize=(6,4))
ax = fig_histo.add_subplot(1,1,1)
# N is the count in each bin, bins is the lower-limit of the bin
N, bins, patches = ax.hist(phi,bins=100, alpha=.8, color='blue',density=False)
# color scale based on the x axis
fracs = bins / bins.max()
# we need to normalize the data to 0..1 for the full range of the colormap
# print(fracs.min(), fracs.max())
norm = colors.Normalize(vmin=fracs.min(), vmax=fracs.max())
xmin,xmax = 0,2*np.pi
# norm = colors.Normalize(vmin=0, vmax=1)
# Now, we'll loop through our objects and set the color of each accordingly
for thisfrac, thispatch in zip(fracs, patches):
    color = cmap_phi(norm(thisfrac))
    thispatch.set_facecolor(color)
# ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
# ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.0%}'.format(y))) 
ax.set_xlim([xmin,xmax])
ax.set_xticks(np.arange(xmin,xmax , (xmax-xmin)/12))
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ax.set_title('Histo timing',fontsize=6)
ax.legend(loc='best')
fig_histo.tight_layout()

fig_histo.savefig('{}_histo.pdf'.format(outfig), format='PDF',dpi=150)


if plotcorr == 'yes':

    # some cleaning for correlations
    kk = np.nonzero(
        np.logical_and(lin>np.nanpercentile(lin,1),
        np.logical_and(lin<np.nanpercentile(lin,99),
        np.logical_and(phi>4.,
        np.logical_and(slope<np.nanpercentile(slope,90),
        np.logical_and(dem<np.nanpercentile(dem,95),
        np.logical_and(slopelos>np.nanpercentile(slopelos,10), 
        np.logical_and(slopelos<np.nanpercentile(slopelos,90), 
        np.logical_and(amp>5, amp <  maxamp
        )))))))))
    phi,amp,lin,dem,slopelos,slope,aspect = phi[kk], amp[kk], lin[kk], dem[kk], slopelos[kk],slope[kk],aspect[kk]

    # compute correlations
    m = np.vstack([phi,amp,lin,dem,slopelos,slope])
    labels=['Phi. (rad)','Amp. (mm)','Vel. (mm/yr)','DEM (m)','Slope in Los','Slope']
    cov = (100*np.corrcoef(m)).astype(int)
    # print()
    # print(labels)
    # print(cov)
    fig = plt.figure(2,figsize=(7,5))
    ax = fig.add_subplot(1,1,1)
    cax = ax.pcolor(cov,vmax=50,vmin=-50,cmap = 'seismic')
    cbar = plt.colorbar(cax, orientation='vertical')
    cbar.set_label('Correlation coefficient')
    plt.xticks(arange(0.5,len(labels)+0.5),labels,rotation='vertical')
    plt.yticks(arange(0.5,len(labels)+0.5),labels)
    fig.savefig('cov_phiall.eps', format='EPS')

    import plot2Ddist
    scatterstyle={'color':'black', 'alpha':0.5, 's':3.}
    styleargs = {'color':'k', 'scatterstyle':scatterstyle}

    # plot2Ddist.plot2DdistsPairs(amp,[slopelos,dem,slope,aspect], mainlabel='Amp. (mm)', labels = [ 'Slope in Los (%)','DEM (m)','Slope (%)','Aspect (deg)'], \
    #     plotcontours=False,plotscatter=True,plotKDE=False,contourKDEthin=5000,thin=100,\
    #     out='jointPDFcorrampall.pdf',corr='no',**styleargs)
    # plot2Ddist.plot2DdistsPairs(phi,[amp,slopelos,dem,slope,aspect], mainlabel='Phi. (rad)', labels = ['Amp. (mm)','Slope in Los (%)','DEM (m)','Slope (%)','Aspect (deg)'], \
    #     plotcontours=False,plotscatter=True,plotKDE=False,scaleview=True,contourKDEthin=5000,thin=100,\
    #     out='jointPDFcorrphiall.pdf',corr='no',**styleargs)
    plot2Ddist.plot2DdistsPairs(lin,[amp,phi,slopelos,dem,slope,aspect], mainlabel='Vel. (mm/yr)', labels = [ 'Amp. (mm)','Phi. (rad)','Slope in Los (%)','DEM (m)','Slope (%)','Aspect (deg)'], \
        plotcontours=False,plotscatter=True,plotKDE=False,contourKDEthin=5000,thin=100,\
        out='jointPDFcorrlinall.pdf',corr='no', **styleargs)

plt.show()
