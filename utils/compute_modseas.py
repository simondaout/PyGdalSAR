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
compute_modseas.py
-------------
Compute the mod of the seasonality 

Usage: clean_phi-amp.py [--cube=<path> ] [--ampfile=<path>] [--phifile=<path>] [--linfile=<path>] [--demerrfile=<path>] \
[--ref_file=<path>] [--slopefile=<path>] [--outampfile=<path>] [--outphifile=<path>] [--crop=<values>] [--lectfile=<path>] \
[--sigampfile=<path>] [--sigphifile=<path>] [--minamp=<value>] [--maxamp=<value>] [--perc_sig=<value>] [--name=<value>]


Options:
-h --help             Show this screen.
--cube=<file>         Path to time series [default: depl_cumule_flat]
--ampfile=<file>      Amplitude map file [default: ampwt_coeff.r4]
--phifile=<file>      phi map file [default: phiwt_coeff.r4]
--linfile=<file>      Linear map file [default: lin_coeff.r4]
--demerrfile=<file>   Path to the dem error file [default: corrdem_coeff.r4]
--ref_file=<file>      Path to the reference file [default: ref_coeff.r4]
--slopefile=<file>     Linear map file [default: LOSslope.r4]
--lectfile=<file>     Path of the lect.in file for r4 format [default: lect_ts.in]
--crop=<values>       Crop data [default: 0,nlines,0,ncol]
--sigampfile=<file>   Uncertainty amplitude map file [default: ampwt_sigcoeff.r4]
--sigphifile=<file>   Uncertainty phi map file [default: phiwt_sigcoeff.r4]
--minamp=<value>      Mask on minimum Amplitude [default: 1.]
[--maxamp=<value>]     Maximum Amplitude limit [default: 2.]
--perc_sig=<value>    Percentile uncertainty for map cleaning [default: 99.]
--name=<value>        Output file name
"""

import numpy as np
from numpy.lib.stride_tricks import as_strided

import scipy.stats as stat
import scipy.optimize as opt
import scipy.linalg as lst
from scipy import ndimage
import gdal, sys

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

np.warnings.filterwarnings('ignore')
try:
    from nsbas import docopt
except:
    import docopt

arguments = docopt.docopt(__doc__)
if arguments["--ampfile"] ==  None:
    arguments["--ampfile"] = 'ampwt_coeff.r4'
if arguments["--phifile"] ==  None:
    arguments["--phifile"] = 'phiwt_coeff.r4'
if arguments["--ref_file"] ==  None:
    arguments["--ref_file"] = 'ref_coeff.r4'
if arguments["--linfile"] ==  None:
    arguments["--linfile"] = 'lin_coeff.r4'
if arguments["--lectfile"] ==  None:
    arguments["--lectfile"] = "lect_ts.in"
if arguments["--slopefile"] ==  None:
    arguments["--slopefile"] = "LOSslope.r4"
if arguments["--sigampfile"] ==  None:
    arguments["--sigampfile"] = 'ampwt_sigcoeff.r4'
if arguments["--sigphifile"] ==  None:
    arguments["--sigphifile"] = 'phiwt_sigcoeff.r4'
if arguments["--cube"] ==  None:
    arguments["--cube"] = 'depl_cumule_flat'
if arguments["--minamp"] ==  None:
    arguments["--minamp"] = 1.5
if arguments["--maxamp"] ==  None:
    arguments["--maxamp"] = 3

fimages='images_retenues'
imref = 0
rad2mm = 1

##################################

# lect cube
nb,idates,dates,base=np.loadtxt(fimages, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
ds = gdal.Open(arguments["--cube"])
if not ds:
  print ('.hdr file time series cube {0}, not found, open {1}'.format(arguments["--cube"],arguments["--lectfile"]))
  ncols, nlines, N = list(map(int, open(arguments["--lectfile"]).readline().split(None, 3)[0:3]))
else:
  hdr = arguments["--cube"] + '.hdr'
  print ('Read ', hdr)
  ncols, nlines, N = ds.RasterXSize, ds.RasterYSize, ds.RasterCount()
if arguments["--crop"] ==  None:
    crop = [0,nlines,0,ncols]
else:
    crop = list(map(float,arguments["--crop"].replace(',',' ').split()))
ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])
# print(ibeg,iend,jbeg,jend)

# Open maps
amp_map=np.fromfile(arguments["--ampfile"],dtype=np.float32)[:nlines*ncols].reshape(nlines,ncols)[ibeg:iend,jbeg:jend]
phi_map = np.fromfile(arguments["--phifile"],dtype=np.float32)[:nlines*ncols].reshape(nlines,ncols)[ibeg:iend,jbeg:jend]
lin_map = np.fromfile(arguments["--linfile"],dtype=np.float32)[:nlines*ncols].reshape(nlines,ncols)[ibeg:iend,jbeg:jend]
sigamp_map=np.fromfile(arguments["--sigampfile"],dtype=np.float32)[:nlines*ncols].reshape(nlines,ncols)[ibeg:iend,jbeg:jend]
sigphi_map = np.fromfile(arguments["--sigphifile"],dtype=np.float32)[:nlines*ncols].reshape(nlines,ncols)[ibeg:iend,jbeg:jend]
slope_map = np.fromfile(arguments["--slopefile"],dtype=np.float32)[:nlines*ncols].reshape(nlines,ncols)[ibeg:iend,jbeg:jend]
ref_map = np.fromfile(arguments["--ref_file"],dtype=np.float32)[:nlines*ncols].reshape(nlines,ncols)[ibeg:iend,jbeg:jend]
if arguments["--demerrfile"] is not None:
    bperp_map = np.fromfile(arguments["--demerrfile"],dtype=np.float32)[:nlines*ncols].reshape(nlines,ncols)[ibeg:iend,jbeg:jend]
else:
    bperp_map = np.zeros((nlines,ncols))[ibeg:iend,jbeg:jend]

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


phi_map[amp_map<np.float(arguments["--minamp"])] = np.float('NaN')
lin_map[amp_map<np.float(arguments["--minamp"])] = np.float('NaN')

# convert phi between 0 and 2pi
phi_map[phi_map<0] = phi_map[phi_map<0] + 2*np.pi

# Initialisation
dmod=[]
flat_disp_pos = []
flat_disp_neg = []
time_pos,time_neg = [],[]
amp_pos,amp_neg = [],[]

##################################

## function invers seasonality
def invers_seas(x,y,std):
    G=np.zeros((len(x),3))
    G[:,0]=np.cos(2*np.pi*(x))
    G[:,1]=np.sin(2*np.pi*(x))
    G[:,2]=np.ones(len(x))
    
    x0 = lst.lstsq(G,y)[0]
    _func = lambda x: np.sum(((np.dot(G,x)-y)/std)**2)
    _fprime = lambda x: 2*np.dot(G.T/std, (np.dot(G,x)-y)/std)
    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=20000,full_output=True,iprint=0)[0] 

    Cd = np.diag(std_neg**2,k=0)
    Cov = np.linalg.inv(Cd)
    a,b = pars[0],pars[1]
    phi=np.arctan2(b,a)
    amp = np.sqrt(a**2+b**2)
    cova = np.linalg.inv(np.dot(G.T,np.dot(Cov,G)))
    siga,sigb,sigc = np.sqrt(np.diag(cova))
    sigamp = np.sqrt(siga**2+sigb**2)
    sigphi = (a*siga+b*sigb)/(a**2+b**2)

    return pars, amp, phi, sigamp, sigphi

def seasonal(time,a,b,c):
    return a*np.cos(2*np.pi*time) + b*np.sin(2*np.pi*time) + c

##################################

# plot
fig_ampphi = plt.figure(0,figsize=(14,4))

ax = fig_ampphi.add_subplot(1,3,1)
im = ax.imshow(slope_map,cmap=cm.Greys,zorder=1)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
ax.set_title('Slope in LOS',fontsize=12)

ax = fig_ampphi.add_subplot(1,3,2)
cmap = cm.rainbow
cax = ax.imshow(slope_map,cmap=cm.Greys,zorder=1)
# no point to take negatif amplitude
im = ax.imshow(np.ma.array(amp_map, mask=np.isnan(amp_map)),cmap=cmap,vmin=arguments["--minamp"],\
    vmax=np.nanpercentile(amp_map,98), alpha=0.2, zorder=2)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
ax.set_title('Amplitude (mm)',fontsize=12)

ax = fig_ampphi.add_subplot(1,3,3)
cmap_phi = cm.tab20b
cax = ax.imshow(slope_map,cmap=cm.Greys,zorder=1)
im = ax.imshow(np.ma.array(phi_map, mask=np.isnan(phi_map)),cmap=cmap_phi,vmin=0,vmax=2*np.pi,alpha=0.8,zorder=2)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
ax.set_title('Timing (rad)',fontsize=12)
fig_ampphi.tight_layout()
# plt.show()
# sys.exit(0)

# Open cube
cube = np.fromfile(arguments["--cube"],dtype=np.float32)[:nlines*ncols*N]
maps = cube.reshape((nlines,ncols,N))[ibeg:iend,jbeg:jend]

dmodt = np.fmod(dates,1)
for t in range((N)):
    dmod.append(dmodt[t])

for i in range(iend-ibeg):
    for j in range(jend-jbeg):
        if (~np.isnan(amp_map[i,j]) and ~np.isnan(phi_map[i,j]) and amp_map[i,j]>np.float(arguments["--maxamp"])):
            temp_pos = (maps[i,j,:] - lin_map[i,j]*(dates[:]-dates[imref]) - ref_map[i,j] \
                - bperp_map[i,j]*(base[:]-base[imref]))/amp_map[i,j]
            for t in range((N)):
                flat_disp_pos.append(temp_pos[t])
                time_pos.append(dmodt[t])
                amp_pos.append(amp_map[i,j])
        elif (~np.isnan(amp_map[i,j]) and ~np.isnan(phi_map[i,j]) and amp_map[i,j]<np.float(arguments["--maxamp"])):
            temp_neg = (maps[i,j,:] - lin_map[i,j]*(dates[:]-dates[imref]) - ref_map[i,j] \
                - bperp_map[i,j]*(base[:]-base[imref]))/amp_map[i,j]
            for t in range((N)):
                flat_disp_neg.append(temp_neg[t])
                time_neg.append(dmodt[t])
                amp_neg.append(amp_map[i,j])

# wrap everything
time_pos = np.array(time_pos).flatten()
amp_pos = np.array(amp_pos).flatten()
time_neg = np.array(time_neg).flatten()
amp_neg = np.array(amp_neg).flatten()
flat_disp_pos = np.array(flat_disp_pos).flatten()
flat_disp_neg = np.array(flat_disp_neg).flatten()
dmod=np.unique(np.array(dmod).flatten())

mean_pos,std_pos=[],[]
mean_neg,std_neg=[],[]
for d in dmod:
    uu = np.flatnonzero(time_pos==d)
    mean_pos.append(np.nanmedian(flat_disp_pos[uu]))
    std_pos.append(np.nanstd(flat_disp_pos[uu]))
for d in dmod:
    uu = np.flatnonzero(time_neg==d)
    mean_neg.append(np.nanmedian(flat_disp_neg[uu]))
    std_neg.append(np.nanstd(flat_disp_neg[uu]))
mean_pos,std_pos = np.array(mean_pos),np.array(std_pos)
mean_neg,std_neg = np.array(mean_neg),np.array(std_neg)

# plot slope postive
fig=plt.figure(1,figsize=(14,5))
ax=fig.add_subplot(1,2,1)
ax.plot(dmod,mean_pos,'o',c='blue',ms=6.,label='Thesh. Amp: {}'.format(arguments["--minamp"])) 
ax.errorbar(dmod,mean_pos,yerr=std_pos,fmt='none',ecolor='blue',alpha=0.1)

try:
    pars, amp, phi, sigamp, sigphi = invers_seas(dmod,mean_pos,std_pos) 
    t = np.arange(1,100)/100.
    ax.plot(t,seasonal(t,pars[0],pars[1],pars[2]),'red',\
        lw=2,label='{:0.1f}+-{:0.1f} * cos(wt - {:0.1f}+-{:0.1f})'.format(amp,sigamp,phi,sigphi))
except:
    pass


ax.set_ylim([-2,2])
ax.set_xlim([0,1])
ax.set_xticks(np.arange(0,1, 1./12))
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
plt.legend(loc='best')
ax.set_title('Amplitudes > {}'.format(arguments["--maxamp"]))

# plot slope negative
ax2=fig.add_subplot(1,2,2)
ax2.plot(dmod,mean_neg,'o',c='blue',ms=6.,label='Thesh. Amp: {}'.format(arguments["--minamp"])) 
ax2.errorbar(dmod,mean_neg,yerr=std_neg,fmt='none',ecolor='blue',alpha=0.1)

try:
    pars, amp, phi, sigamp, sigphi = invers_seas(dmod,mean_neg,std_neg) 
    t = np.arange(1,100)/100.
    ax2.plot(t,seasonal(t,pars[0],pars[1],pars[2]),'red',\
        lw=2,label='{:0.1f}+-{:0.1f} * cos(wt - {:0.1f}+-{:0.1f})'.format(amp,sigamp,phi,sigphi))
except:
    pass


ax2.set_ylim([-2,2])
ax2.set_xlim([0,1])
ax2.set_xticks(np.arange(0,1, 1./12))
ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
plt.legend(loc='best')
ax2.set_title(' {} < Amplitudes < {}'.format(arguments["--minamp"],arguments["--maxamp"]))

fig.tight_layout()
fig.savefig('{}-modseas.pdf'.format(arguments["--name"]), format='PDF',dpi=80)
plt.show()
