#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################

############################################
# Author        : Simon DAOUT (Oxford)
############################################

"""\
invert_ramp_topo_unw.py
-------------
Estimates atmospheric phase/elevation correlations or/and azimuthal and range ramps polynomial coefficients
on unwrapped interferograms (ROI_PAC/GAMMA/GTIFF). Temporal inversion of all coeffient with a strong weight for
small temporal baselines and lage cover interferograms. Reconstruction of the empirical phase correction.

usage: invert_ramp_topo_unw.py --int_list=<path> [--ref_zone=<jstart,jend,istart,iend> ] [--int_path=<path>] [--out_path=<path>] \
[--prefix=<value>] [--suffix=<value>] [--rlook=<value>]  [--ref=<path>] [--format=<value>] \
[--flat=<0/1/2/3/4/5/6>] [--topofile=<path>] [--ivar=<0/1>] [--nfit=<0/1>] [--tsinv=<yes/no>]\
[--estim=yes/no] [--mask=<path>] [--threshold_mask=<value>]  \
[--cohpixel=<yes/no>] [--threshold_coh=<value>] [--ibeg_mask=<value>] [--iend_mask=<value>] \
[--perc=<value>] [--perc_topo=<value>] [--min_topo=<value>] [--max_topo=<value>] [--perc_slope=<value>] [--samp=<value>] \
[--plot=<yes/no>] [--suffix_output=<value>]\
[--crop_emp=<values>] [--nproc=<nb_cores>] [--cpt=<path>] 

--int_list PATH       Text file containing list of interferograms dates in two colums, $data1 $date2
--int_path PATh       Relative path to input interferograms directory
--ref_zone=<lin_start,lin_end,col_start,col_end> Starting and ending lines and col numbers where phase is set to zero [default: None] 
--prefix VALUE        Prefix name $prefix$date1-$date2$suffix_$rlookrlks.unw [default: '']
--suffix value        Suffix name $prefix$date1-$date2$suffix_$rlookrlks.unw [default: '']
--rlook value         look int. $prefix$date1-$date2$suffix_$rlookrlks.unw [default: 0]
--ref PATH            Path to a reference image to define format and size. Necessary if topofile is not given [default:None]     
--flat PATH           Remove a spatial ramp.If short acquisition, short is automatically set to 3.
0: ref frame [default], 1: range ramp ax+b , 2: azimutal ramp ay+b, 
3: ax+by+c, 4: ax+by+cxy+d 5: ax**2+bx+d, 6: ay**2+by+c
--topofile  PATH      Path to the radar_look.hgt file. If not None add phase/elevation relationship in the relation [default:None]
--nfit VALUE          fit degree in azimuth or in elevation
--ivar VALUE          define phase/elevation relationship: ivar=0 function of elevation, ivar=1 crossed function of azimuth and elevation
--nproc=<nb_cores>    Use <nb_cores> local cores to create delay maps [Default: 4]


if ivar=0 and nfit=0, add linear elev. term (z) to ramps estimation defined by the flat argument such as:
0: ref frame+ez [default], 1: range ramp ax+b+ez , 2: azimutal ramp ay+b+ez, 
3: ax+by+c+ez, 4: ax+by+cxy+d+ez, 5: ax**2+bx+d+ez, 6: ay**2+by+c+ez

if ivar=0 and nfit=1, add quadratic elev. term (z) to ramps estimation defined by the flat argument such as:
0: ref frame+ez+fz**2 [default], 1: range ramp ax+b+ez+fz**2 , 2: azimutal ramp ay+b+ez+fz**2, 
3: ax+by+c+ez+fz**2, 4: ax+by+cxy+d+ez+fz**2, 5: ax**2+bx+d+ez+fz**2, 6: ay**2+by+c+ez+fz**2

if ivar=1 and nfit=0, add cross function of elev. (z) and azimuth to ramps estimation defined by the flat argument such as:
0: ref frame+ez+fz*az [default], 1: range ramp ax+b+ez+fz*az , 2: azimutal ramp ay+b+ez+fz*az, 
3: ax+by+c+ez+fz*az, 4: ax+by+cxy+d+ez+fz*az, 5: ax**2+bx+d+ez+fz*az, 6: ay**2+by+c+ez+fz*az

if ivar=1 and nfit=1, add quadratic cross function of elev. (z) and azimuth to ramps estimation defined by the flat argument such as:
0: ref frame+ez+fz*az+g(z*az)**2 [default], 1: range ramp ax+b+ez+fz*az+g*(z*az)**2 , 2: azimutal ramp ay+b+ez+fz*az+g*(z*az)**2, 
3: ax+by+c+ez+fz*az+g*(z*az)**2, 4: ax+by+cxy+d+ez+fz*az+g*(z*az)**2, 5: ax**2+bx+d+ez+fz*az+g*(z*az)**2, 6: ay**2+by+c+ez+fz*az+g*(z*az)**2

--tsinv yes/no        If yes, invert corrected phase into time series [default:no]
--estim yes/no        If yes, do the estimation, otherwise read input files corection_matrix and liste_coeff_ramps.txt [default:yes]
--mask PATH           Mask in .r4 format. Keep only values > threshold_mask. [default:None]
--threshold_mask      Thresclean_r4.pyhold on mask: take only values > threshold_mask [default: -1]
--cohpixel  yes/no    If Yes, use amplitude interferogram to weight and mask pixels (e.g Coherence, Colinearity, Amp Filter) [default: no]
--format VALUE        Format input files: ROI_PAC, GAMMA, GTIFF [default: ROI_PAC]
--threshold_coh VALUE Threshold on cohpixel file [default:0]
--ibeg_mask VALUE     Start line number for the mask [default: None]
--iend_mask VALUE     Stop line number for the mask [default: None]  
--perc VALUE          Percentile of hidden LOS pixel for the estimation and clean outliers [default:98.]
--perc_topo VALUE     Percentile of hidden elevation pixel for the estimation and clean outliers [default:98.]
--min_topo VALUE      Minimum elevation for the estimation. Overwrite perc_topo argument [default: perc_topo]
--max_topo VALUE      Maxmimum elevation for the estimation. Overwrite perc_topo argument [default: perc_topo]
--perc_slope VALUE    Percentile of hidden slope-elevation pixel for the estimation and clean outliers [default:92.]
--samp=<value>        Undersampling for empirical estimation [default: 2]
--plot yes/no         If yes, plot figures for each ints [default: no]
--suffix_output value Suffix output file name $prefix$date1-$date2$suffix$suffix_output [default:_corrunw]
--cpt==<path>         Indicate colormap for plots [default: cm.rainbow]
--crop_emp VALUE    Define a region of interest for the spatial estimatiom (ramp+phase/topo) [default: 0,nlines,0,ncol]
"""


# gdal
import gdal
gdal.UseExceptions()

# system
from os import path, environ, getcwd
import os, sys

# plot
import matplotlib
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap

from datetime import datetime
import time
import logging

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided

# scipy
import scipy
import scipy.optimize as opt
import scipy.linalg as lst
import scipy.ndimage

try:
    from nsbas import docopt
except:
    import docopt
import shutil

from contextlib import contextmanager
from functools import wraps, partial
import multiprocessing

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning) 

##################################################################################
###  Extras functions and context maganers
##################################################################################

def makedirs(name):
    if os.path.exists(name):
        return
    os.makedirs(name)

def date2dec(dates):
    dates  = np.atleast_1d(dates)
    times = []
    for date in dates:
        x = datetime.strptime('{}'.format(date),'%Y%m%d')
        #x = datetime.strptime('{}'.format(date),'%Y-%m-%dT%H:%M:%S.%fZ')
        dec = float(x.strftime('%j'))/365.1
        year = float(x.strftime('%Y'))
        times.append(year + dec)
    return times

# Timer for all the functions
class ContextDecorator(object):
    def __call__(self, f):
        @wraps(f)
        def decorated(*args, **kwds):
            with self:
                try:
                    return f(*args, **kwds)
                except (KeyboardInterrupt, SystemExit):
                    raise
                except:
                    Exception('{0} Failed !'.format(f))
                    raise
        return decorated

class TimeIt(ContextDecorator):
    def __enter__(self):
        self.start = datetime.now()
        print('Starting time process: {0}'.format(self.start))
    def __exit__(self, type, value, traceback):
        print('Time process: {0}s'.format((datetime.now() - self.start).total_seconds()))


def checkinfile(file):
    if path.exists(file) is False:
        logger.critical("File: {0} not found, Exit !".format(file))
        print("File: {0} not found in {1}, Exit !".format(file,getcwd()))

# create generator for pool
@contextmanager
def poolcontext(*arg, **kargs):
    pool = multiprocessing.Pool(*arg, **kargs)
    yield pool
    pool.terminate()
    pool.join()

#####################################################################################
# FUNCTIONS
#####################################################################################

def estim_ramp(los,los_clean,topo_clean,az,rg,order,rms,nfit,ivar,los_ref,rg_ref,az_ref,topo_ref):
    """
    Empircal estmation function on flatten los vector
    """

    # initialise full vector 
    sol = np.zeros((13))
    #0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2
    # initialize correction
    corr = np.zeros((mlines,mcols))
 
    if radar is None:
        data = los_clean
        losbins = los_clean
        topobins = topo_clean
        rgbins, azbins = rg, az
        losstd = rms

    else:
        # lets try to digitize to improve the fit
        # digitize data in bins, compute median and std
        bins = np.arange(minelev,maxelev,abs(maxelev-minelev)/500.)
        inds = np.digitize(topo_clean,bins)
        topobins = []
        losbins = []
        losstd = []
        azbins, rgbins = [], []
        los_clean2, topo_clean2, az_clean2, rg_clean2, rms_clean2 = [], [], [], [], []
        for j in range(len(bins)-1):
                uu = np.flatnonzero(inds == j)
                if len(uu)>200:
                    topobins.append(bins[j] + (bins[j+1] - bins[j])/2.)

                    # do a clean within the bin
                    indice = np.flatnonzero(np.logical_and(los_clean[uu]>np.percentile(\
                        los_clean[uu],100-perc),los_clean[uu]<np.percentile(los_clean[uu],perc)))

                    losstd.append(np.nanstd(los_clean[uu][indice]))
                    losbins.append(np.nanmedian(los_clean[uu][indice]))
                    azbins.append(np.nanmedian(az[uu][indice]))
                    rgbins.append(np.nanmedian(rg[uu][indice]))

                    # remove outliers from data
                    los_clean2.append(los_clean[uu][indice])
                    az_clean2.append(az[uu][indice])
                    rg_clean2.append(rg[uu][indice])
                    topo_clean2.append(topo_clean[uu][indice])
                    # new rms is clean rms time the standard deviation of the los within the bin
                    rms_clean2.append(rms[uu][indice]*np.nanstd(los_clean[uu][indice]))

        # digitize data and model
        losbins = np.array(losbins)
        losstd = np.array(losstd)
        topobins = np.array(topobins)
        rgbins, azbins = np.array(rgbins),np.array(azbins)
        
        # data and model cleaned a second time by sliding median
        los_clean = np.concatenate(los_clean2)
        az, rg = np.concatenate(az_clean2), np.concatenate(rg_clean2)
        topo_clean = np.concatenate(topo_clean2)
        rms = np.concatenate(rms_clean2)
        del los_clean2, az_clean2, rg_clean2, topo_clean2

        if order == 0 and ivar == 0:
            data = losbins
            rms = losstd

        else:
            data = los_clean

    if order==0:  
    #0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2

        if radar is None:
            pars = los_ref
            sol[8] = los_ref
            logger.info('Remove ref frame within the ref area %f'%(los_ref))

            # build total G matrix
            G=np.zeros((len(los),1))
            G[:,0] = 1

        else:
            if ivar==0 and nfit==0:
                G=np.zeros((len(data),2))
                G[:,0] = 1
                G[:,1] = topobins

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0

                sol[8] = pars[0]; sol[9] = pars[1]
                logger.info('Remove ref frame %f + %f z'%(pars[0],pars[1]))

                # build total G matrix
                G=np.zeros((len(los),2))
                G[:,0] = 1
                G[:,1] = elev_map.flatten()

            elif ivar==0 and nfit==1:
                G=np.zeros((len(data),3))
                G[:,0] = 1
                G[:,1] = topobins
                G[:,2] = topobins**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                
                sol[8] = pars[0]; sol[9] = pars[1]; sol[10] = pars[2]
                logger.info('Remove ref frame %f + %f z + %f z**2'%(pars[0],pars[1],pars[2]))

                # build total G matrix
                G=np.zeros((len(los),3))
                G[:,0] = 1
                G[:,1] = elev_map.flatten()
                G[:,2] = elev_map.flatten()**2
            
            elif ivar==1 and nfit==0:
                G=np.zeros((len(data),3))
                G[:,0] = 1
                G[:,1] = topo_clean
                G[:,2] = az*topo_clean

                # ramp inversion 
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0

                sol[8] = pars[0]; sol[9] = pars[1]; sol[11] = pars[2]
                logger.info('Remove ref frame %f + %f z + %f az*z'%(pars[0],pars[1],pars[2]))

                # build total G matrix
                G=np.zeros((len(los),3))
                G[:,0] = 1
                G[:,1] = elev_map.flatten()
                G[:,2] = elev_map.flatten()
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,2] *= i - ibeg

            elif ivar==1 and nfit==1:
                G=np.zeros((len(data),5))
                G[:,0] = 1
                G[:,1] = topo_clean
                G[:,2] = topo_clean**2
                G[:,3] = az*topo_clean
                G[:,4] = (az*topo_clean)**2

                # ramp inversion 
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[8] = pars[0]; sol[9] = pars[1]; sol[10] = pars[2]; sol[11] = pars[3]; sol[12] = pars[4]
                logger.info('Remove ref frame %f + %f z + %f z**2 + %f az*z + %f (az*z)**2'%(pars[0],pars[1],pars[2],pars[3],pars[4]))

                # build total G matrix
                G=np.zeros((len(los),4))
                G[:,0] = 1
                G[:,1] = elev_map.flatten()
                G[:,2] = elev_map.flatten()**2
                G[:,3] = elev_map.flatten()
                G[:,4] = elev_map.flatten()**2
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,3] *= i - ibeg
                    G[i*mcols:(i+1)*mcols,4] *= (i-ibeg)**2


    elif order==1: # Remove a range ramp ay+b for each maps (y = col)

        if radar is None:
            G=np.zeros((len(data),2))
            G[:,0] = rg
            G[:,1] = 1

            # ramp inversion
            x0 = lst.lstsq(G,data)[0]
            try:
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
            except:
                pars = x0
            sol[2] = pars[0]; sol[8] = pars[1]
            logger.info('Remove ramp %f r + %f'%(pars[0],pars[1]))

            # build total G matrix
            G=np.zeros((len(los),2))
            for i in range(mlines):
                G[i*mcols:(i+1)*mcols,0] = np.arange((mcols))  - jbeg
            G[:,1] = 1


        else:
            if ivar==0 & nfit==0:
                G=np.zeros((len(data),3))
                G[:,0] = rg
                G[:,1] = 1
                G[:,2] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[8] = pars[1]; sol[9] = pars[2]
                logger.info('Remove ramp %f r + %f + %f z '%(pars[0],pars[1],pars[2]))

                # build total G matrix
                G=np.zeros((len(los),3))
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) - jbeg
                G[:,1] = 1
                G[:,2] = elev_map.flatten()

            if ivar==0 & nfit==1:
                G=np.zeros((len(data),4))
                G[:,0] = rg
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[8] = pars[1]; sol[9] = pars[2]; sol[10] = pars[3]
                logger.info('Remove ramp %f r + %f + %f z + %f z**2'%(pars[0],pars[1],pars[2],pars[3]))

                # build total G matrix
                G=np.zeros((len(los),4))
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) - jbeg 
                G[:,1] = 1
                G[:,2] = elev_map.flatten()
                G[:,3] = elev_map.flatten()**2
            
            elif ivar==1 and nfit==0:
                G=np.zeros((len(data),4))
                G[:,0] = y
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean*az

                # ramp inversion
                #y**3 y**2 y x**3 x**2 x xy**2 xy cst z z**2 yz yz**2 
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[8] = pars[1]; sol[9] = pars[2]; sol[11] = pars[3]
                logger.info('Remove ramp %f r + %f + %f z + %f z*az'%(pars[0],pars[1],pars[2],pars[3]))

                # build total G matrix
                G=np.zeros((len(los),4))
                G[:,1] = 1
                G[:,2] = elev_map.flatten()
                G[:,3] = elev_map.flatten()
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) - jbeg 
                    G[i*mcols:(i+1)*mcols,3] *= i - ibeg

            elif ivar==1 and nfit==1:
                G=np.zeros((len(data),6))
                G[:,0] = rg
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean**2
                G[:,4] = topo_clean*az
                G[:,5] = (topo_clean*az)**2

                # ramp inversion
                #y**3 y**2 y x**3 x**2 x xy**2 xy cst z z**2 yz yz**2 
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[8] = pars[1]; sol[9] = pars[2]; sol[10] = pars[3]; sol[11] = pars[4]; sol[12] = pars[5]
                logger.info('Remove ramp %f r + %f + %f z + %f z**2 + %f z*az'%(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5]))

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,1] = 1
                G[:,2] = elev_map.flatten()
                G[:,3] = elev_map.flatten()**2
                G[:,4] = elev_map.flatten()
                G[:,5] = elev_map.flatten()**2
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) - jbeg
                    G[i*mcols:(i+1)*mcols,4] *= i - ibeg
                    G[i*mcols:(i+1)*mcols,5] *= (i - ibeg)**2 

        
    elif order==2: # Remove an azimutal ramp ax+b for each maps (x is row)
    #y**3 y**2 y x**3 x**2 x xy**2 xy cst z z**2 z*az az*z**2

        if radar is None:
            G=np.zeros((len(data),2))
            G[:,0] = az
            G[:,1] = 1

            # ramp inversion
            x0 = lst.lstsq(G,data)[0]
            try:
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
            except:
                pars = x0
            sol[5] = pars[0]; sol[8] = pars[1]
            logger.info('Remove ramp %f az + %f'%(pars[0],pars[1]))

            # build total G matrix
            G=np.zeros((len(los),2))
            for i in range(mlines):
                G[i*mcols:(i+1)*mcols,0] = i  
            G[:,1] = 1

        else:
            if ivar==0 & nfit==0:
                G=np.zeros((len(data),3))
                G[:,0] = az
                G[:,1] = 1
                G[:,2] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[5] = pars[0]; sol[8] = pars[1]; sol[9] = pars[2]
                logger.info('Remove ramp %f az + %f + %f z'%(pars[0],pars[1],pars[2]))

                # build total G matrix
                G=np.zeros((len(los),3))
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = i - ibeg 
                G[:,1] = 1
                G[:,2] = elev_map.flatten()

            if ivar==0 & nfit==1:
                G=np.zeros((len(data),4))
                G[:,0] = az
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[5] = pars[0]; sol[8] = pars[1]; sol[9] = pars[2]; sol[10] = pars[-1]
                logger.info('Remove ramp %f az + %f + %f z + %f z**2'%(pars[0],pars[1],pars[2],pars[3]))

                # build total G matrix
                G=np.zeros((len(los),4))
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = i - ibeg
                G[:,1] = 1
                G[:,2] = elev_map.flatten()
                G[:,3] = elev_map.flatten()**2
            
            elif ivar==1 and nfit==0:
                G=np.zeros((len(data),4))
                G[:,0] = az
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean*az

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[5] = pars[0]; sol[8] = pars[1]; sol[9] = pars[2]; sol[11] = pars[3]
                logger.info('Remove ramp %f az + %f + %f z + %f z*az'%(pars[0],pars[1],pars[2],pars[3]))

                # build total G matrix
                G=np.zeros((len(los),4))
                G[:,1] = 1
                G[:,2] = elev_map.flatten()
                G[:,3] = elev_map.flatten()
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = i - ibeg
                    G[i*mcols:(i+1)*mcols,3] *= i - ibeg

            elif ivar==1 and nfit==1:
                G=np.zeros((len(data),5))
                G[:,0] = az
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean*az
                G[:,4] = (topo_clean*az)**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[5] = pars[0]; sol[8] = pars[1]; sol[9] = pars[2]; sol[11] = pars[3]; sol[12] = pars[4]
                logger.info('Remove ramp %f az + %f + %f z + %f z*az + %f (z*az)**2'%(pars[0],pars[1],pars[2],pars[3],pars[4]))

                # build total G matrix
                G=np.zeros((len(los),5))
                G[:,1] = 1
                G[:,2] = elev_map.flatten()
                G[:,3] = elev_map.flatten()
                G[:,3] = elev_map.flatten()**2
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = i - ibeg
                    G[i*mcols:(i+1)*mcols,3] *= i - ibeg
                    G[i*mcols:(i+1)*mcols,4] *= (i-ibeg)**2

    elif order==3: # Remove a ramp ay+bx+c for each maps
    #y**3 y**2 y x**3 x**2 x xy**2 xy cst z z**2 z*az az*z**2

        if radar is None:
            G=np.zeros((len(data),3))
            G[:,0] = rg
            G[:,1] = az
            G[:,2] = 1

            # ramp inversion
            x0 = lst.lstsq(G,data)[0]
            try:
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
            except:
                pars = x0
            sol[2] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]
            logger.info('Remove ramp %f r  + %f az + %f'%(pars[0],pars[1],pars[2]))

            # build total G matrix
            G=np.zeros((len(los),3))
            for i in range(mlines):
                G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) - jbeg 
                G[i*mcols:(i+1)*mcols,1] = i - ibeg
            G[:,2] = 1

        else:
            if ivar==0 and nfit==0:
                G=np.zeros((len(data),4))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = 1
                G[:,3] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:    
                    pars = x0
                sol[2] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]
                logger.info('Remove ramp %f r  + %f az + %f + %f z '%(pars[0],pars[1],pars[2],pars[3]))

                # build total G matrix
                G=np.zeros((len(los),4))
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) - jbeg
                    G[i*mcols:(i+1)*mcols,1] = i - ibeg
                G[:,2] = 1
                G[:,3] = elev_map.flatten()

            if ivar==0 and nfit==1:
                G=np.zeros((len(data),5))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = 1
                G[:,3] = topo_clean
                G[:,4] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]; sol[10] = pars[4]
                logger.info('Remove ramp %f r  + %f az + %f + %f z + %f z**2'%(pars[0],pars[1],pars[2],pars[3],pars[4]))

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) - jbeg 
                    G[i*mcols:(i+1)*mcols,1] = i - ibeg   
                G[:,2] = 1
                G[:,3] = elev_map.flatten()
                G[:,4] = elev_map.flatten()**2
            
            elif ivar==1 and nfit==0:
                # y**3 y**2 y x**3 x**2 x xy**2 xy cst z z**2 z*az az*z**2
                G=np.zeros((len(data),))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = 1
                G[:,3] = topo_clean
                G[:,4] = topo_clean*az

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]; sol[11] = pars[4]
                logger.info('Remove ramp %f r  + %f az + %f + %f z + %f z*az'%(pars[0],pars[1],pars[2],pars[3],pars[4]))

                # build total G matrix
                G=np.zeros((len(los),5))
                G[:,2] = 1
                G[:,3] = elev_map.flatten()
                G[:,4] = elev_map.flatten()
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) - jbeg 
                    G[i*mcols:(i+1)*mcols,1] = i - ibeg
                    G[i*mcols:(i+1)*mcols,4] *= i - ibeg
            
            elif ivar==1 and nfit==1:
                G=np.zeros((len(data),7))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = 1
                G[:,3] = topo_clean
                G[:,4] = topo_clean**2
                G[:,5] = topo_clean*az
                G[:,6] = (topo_clean*az)**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                #0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2    
                sol[2] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]; sol[10] = pars[4]; sol[11] = pars[5]; sol[12] = pars[6]
                logger.info('Remove ramp %f r  + %f az + %f + %f z + %f z**2 + %f z*az + %f (z*az)**2'%(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5],pars[6]))

                # build total G matrix
                G=np.zeros((len(los),7))
                G[:,2] = 1
                G[:,3] = elev_map.flatten()
                G[:,4] = elev_map.flatten()**2
                G[:,5] = elev_map.flatten()
                G[:,6] = elev_map.flatten()**2
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) - jbeg
                    G[i*mcols:(i+1)*mcols,1] = i - ibeg
                    G[i*mcols:(i+1)*mcols,5] *= i - ibeg
                    G[i*mcols:(i+1)*mcols,6] *= (i - ibeg)**2


    elif order==4:
    #y**3 y**2 y x**3 x**2 x xy**2 xy cst z z**2 z*az az*z**2

        if radar is None:
            G=np.zeros((len(data),4))
            G[:,0] = rg
            G[:,1] = az
            G[:,2] = rg*az
            G[:,3] = 1

            # ramp inversion
            x0 = lst.lstsq(G,data)[0]
            try:
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
            except:
                pars = x0
            sol[2] = pars[0]; sol[5] = pars[1]; sol[7] = pars[2]; sol[8] = pars[3]
            logger.info('Remove ramp %f r + %f az  + %f r*az + %f'%(pars[0],pars[1],pars[2],pars[3]))

            # build total G matrix
            G=np.zeros((len(los),4))
            for i in range(mlines):
                G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) - jbeg
                G[i*mcols:(i+1)*mcols,1] = i - ibeg 
                G[i*mcols:(i+1)*mcols,2] = (i - ibeg) * (np.arange((mcols))-jbeg)    
            G[:,3] = 1

        else:
            if ivar==0 and nfit==0:
                G=np.zeros((len(data),5))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = rg*az
                G[:,3] = 1
                G[:,4] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[5] = pars[1]; sol[7] = pars[2]; sol[8] = pars[3]; sol[9] = pars[4]
                logger.info('Remove ramp %f r, %f az  + %f r*az + %f + %f z'%(pars[0],pars[1],pars[2],pars[3],pars[4]))

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) - jbeg
                    G[i*mcols:(i+1)*mcols,1] = i - ibeg
                    G[i*mcols:(i+1)*mcols,2] = (i - ibeg) * (np.arange((mcols))-jbeg)
                G[:,3] = 1
                G[:,4] = elev_map.flatten()

            if ivar==0 and nfit==1:
                G=np.zeros((len(data),6))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = rg*az
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[5] = pars[1]; sol[7] = pars[2]; sol[8] = pars[3]; sol[9] = pars[4]; sol[10] = pars[-1]
                logger.info('Remove ramp %f r, %f az  + %f r*az + %f + %f z+ %f z**2'%(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5]))

                # build total G matrix
                G=np.zeros((len(los),6))
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) - jbeg
                    G[i*mcols:(i+1)*mcols,1] = i - ibeg
                    G[i*mcols:(i+1)*mcols,2] = (i - ibeg) * (np.arange((mcols)) - jbeg)
                G[:,3] = 1
                G[:,4] = elev_map.flatten()
                G[:,5] = elev_map.flatten()**2

            elif ivar==1 and nfit==0:
                G=np.zeros((len(data),6))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = rg*az
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean*az

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[5] = pars[1]; sol[7] = pars[2]; sol[8] = pars[3]; sol[9] = pars[4]; sol[11] = pars[5]
                logger.info('Remove ramp %f r, %f az  + %f r*az + %f + %f z + %f z*az'%(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5]))

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,3] = 1
                G[:,4] = elev_map.flatten()
                G[:,5] = elev_map.flatten()
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols))
                    G[i*mcols:(i+1)*mcols,1] = i
                    G[i*mcols:(i+1)*mcols,2] = (i) * (np.arange((mcols)))
                    G[i*mcols:(i+1)*mcols,5] *= i

            elif ivar==1 and nfit==1:
                G=np.zeros((len(data),8))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = rg*az
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean**2
                G[:,6] = topo_clean*az
                G[:,7] = (topo_clean*az)**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[5] = pars[1]; sol[7] = pars[2]; sol[8] = pars[3]; sol[9] = pars[4]; sol[10] = pars[5]; sol[11] = pars[6]; sol[12] = pars[7]
                logger.info('Remove ramp %f r, %f az  + %f r*az + %f + %f z + %f z**2 + %f z*az + %f (z*az)**2'%(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5],pars[6],pars[7]))

                # build total G matrix
                G=np.zeros((len(los),8))
                G[:,3] = 1
                G[:,4] = elev_map.flatten()
                G[:,5] = elev_map.flatten()**2
                G[:,6] = elev_map.flatten()
                G[:,7] = elev_map.flatten()**2
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) - jbeg
                    G[i*mcols:(i+1)*mcols,1] = i - ibeg
                    G[i*mcols:(i+1)*mcols,2] = (i - ibeg) * (np.arange((mcols)) - jbeg)
                    G[i*mcols:(i+1)*mcols,6] *= i - ibeg
                    G[i*mcols:(i+1)*mcols,7] *= (i - ibeg)**2

    elif order==5:
    #0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2

        if radar is None:
            G=np.zeros((len(data),4))
            G[:,0] = rg**2
            G[:,1] = rg
            G[:,2] = az
            G[:,3] = 1

            # ramp inversion
            x0 = lst.lstsq(G,data)[0]
            try:
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
            except:
                pars = x0
            sol[1] = pars[0]; sol[2] = pars[1]; sol[5] = pars[2]; sol[8] = pars[3]
            logger.info('Remove ramp %f r**2 %f r + %f az + %f'%(pars[0],pars[1],pars[2],pars[3]))

            # build total G matrix
            G=np.zeros((len(los),4))
            for i in range(mlines):
                G[i*mcols:(i+1)*mcols,0] = (np.arange((mcols))-jbeg)**2
                G[i*mcols:(i+1)*mcols,1] = np.arange((mcols)) - jbeg 
                G[i*mcols:(i+1)*mcols,2] = i - ibeg 
            G[:,3] = 1

        else:
            if ivar==0 and nfit==0:
                G=np.zeros((len(data),5))
                G[:,0] = rg**2
                G[:,1] = rg
                G[:,2] = az
                G[:,3] = 1
                G[:,4] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[1] = pars[0]; sol[2] = pars[1]; sol[5] = pars[2]; sol[8] = pars[3]; sol[9] = pars[4]
                logger.info('Remove ramp %f r**2, %f r + %f az + %f + %f z'%(pars[0],pars[1],pars[2],pars[3],pars[4]))

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = (np.arange((mcols))-jbeg)**2
                    G[i*mcols:(i+1)*mcols,1] = np.arange((mcols)) - jbeg
                    G[i*mcols:(i+1)*mcols,2] = i - ibeg
                G[:,3] = 1
                G[:,4] = elev_map.flatten()

            elif ivar==0 and nfit==1:
                G=np.zeros((len(data),6))
                G[:,0] = rg**2
                G[:,1] = rg
                G[:,2] = az
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[1] = pars[0]; sol[2] = pars[1]; sol[5] = pars[2]; sol[8] = pars[3]; sol[9] = pars[4]; sol[10] = pars[-1]
                logger.info('Remove ramp %f r**2, %f r  + %f az + %f + %f z + %f z**2'%(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5]))

                # build total G matrix
                G=np.zeros((len(los),6))
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = (np.arange((mcols))-jbeg)**2
                    G[i*mcols:(i+1)*mcols,1] = np.arange((mcols)) - jbeg 
                    G[i*mcols:(i+1)*mcols,2] *= i - ibeg
                G[:,3] = 1
                G[:,4] = elev_map.flatten()
                G[:,5] = elev_map.flatten()**2
            
            elif ivar==1 and nfit==0:
                G=np.zeros((len(data),6))
                G[:,0] = rg**2
                G[:,1] = rg
                G[:,2] = az
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean*az

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[1] = pars[0]; sol[2] = pars[1]; sol[5] = pars[2]; sol[8] = pars[3]; sol[9] = pars[4]; sol[11] = pars[5]
                logger.info('Remove ramp %f r**2, %f r   + %f az + %f + %f z + %f z*az'%(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5]))

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,3] = 1
                G[:,4] = elev_map.flatten()
                G[:,5] = elev_map.flatten()
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = (np.arange((mcols)) - jbeg)**2
                    G[i*mcols:(i+1)*mcols,1] = np.arange((mcols)) - jbeg
                    G[i*mcols:(i+1)*mcols,2] = i - ibeg
                    G[i*mcols:(i+1)*mcols,5] *= i - ibeg

            elif ivar==1 and nfit==1:
                G=np.zeros((len(data),8))
                G[:,0] = rg**2
                G[:,1] = rg
                G[:,2] = az
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean**2
                G[:,6] = topo_clean*az
                G[:,7] = (topo_clean*az)**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[1] = pars[0]; sol[2] = pars[1]; sol[5] = pars[2]; sol[8] = pars[3]; sol[9] = pars[4]; sol[10] = pars[5]; sol[11] = pars[6]; sol[12] = pars[7]
                logger.info('Remove ramp %f r**2, %f r  + %f az + %f + %f z + %f z**2 + %f z*az + %f (z*az)**2'%(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5],pars[6],pars[7]))

                # build total G matrix
                G=np.zeros((len(los),8))
                G[:,3] = 1
                G[:,4] = elev_map.flatten()
                G[:,5] = elev_map.flatten()**2
                G[:,6] = elev_map.flatten()
                G[:,7] = elev_map.flatten()**2
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = (np.arange((mcols)) - jbeg)**2
                    G[i*mcols:(i+1)*mcols,1] = np.arange((mcols)) - jbeg
                    G[i*mcols:(i+1)*mcols,2] = i - ibeg
                    G[i*mcols:(i+1)*mcols,6] *= i - ibeg
                    G[i*mcols:(i+1)*mcols,7] *= (i-ibeg)**2

    elif order==6:
    #0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2

        if radar is None:
            G=np.zeros((len(data),3))
            G[:,0] = az**2
            G[:,1] = az
            G[:,2] = 1

            # ramp inversion
            x0 = lst.lstsq(G,data)[0]
            try:
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
            except:
                pars = x0
            sol[4] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]
            logger.info('Remove ramp %f az**2 %f az  + %f'%(pars[0],pars[1],pars[2]))

            # build total G matrix
            G=np.zeros((len(los),3))
            for i in range(mlines):
                G[i*mcols:(i+1)*mcols,0] = (i)**2
                G[i*mcols:(i+1)*mcols,1] = (i) 
            G[:,2] = 1

        else:
            if ivar==0 and nfit==0:
                G=np.zeros((len(data),4))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = 1
                G[:,3] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[4] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]
                logger.info('Remove ramp %f az**2, %f az  + %f + %f z'%(pars[0],pars[1],pars[2],pars[3]))

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = (i-ibeg)**2
                    G[i*mcols:(i+1)*mcols,1] = i - ibeg 
                G[:,2] = 1
                G[:,3] = elev_map.flatten()

            elif ivar==0 and nfit==1:
                G=np.zeros((len(data),5))
                G[:,0] = az**2
                G[:,1] = az
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[4] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]; sol[10] = pars[4]
                logger.info('Remove ramp %f az**2, %f az  + %f + %f z + %f z**2'%(pars[0],pars[1],pars[2],pars[3],pars[4]))

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = (i-ibeg)**2
                    G[i*mcols:(i+1)*mcols,1] = i - ibeg
                G[:,2] = 1
                G[:,3] = elev_map.flatten()
                G[:,4] = elev_map.flatten()**2
            
            elif ivar==1 and nfit==0:
                G=np.zeros((len(data),5))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = 1
                G[:,3] = topo_clean
                G[:,4] = topo_clean*az

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                sol[4] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]; sol[11] = pars[4];
                logger.info('Remove ramp %f az**2, %f az + %f + %f z + %f z*az'%(pars[0],pars[1],pars[2],pars[3],pars[4]))

                # build total G matrix
                G=np.zeros((len(los),5))
                G[:,2] = 1
                G[:,3] = elev_map.flatten()
                G[:,4] = elev_map.flatten()
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = (i-ibeg)**2
                    G[i*mcols:(i+1)*mcols,1] = i - ibeg
                    G[:,4] *= i 

            elif ivar==1 and nfit==1:
                G=np.zeros((len(data),7))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = 1
                G[:,3] = topo_clean
                G[:,4] = topo_clean**2
                G[:,5] = topo_clean*az
                G[:,6] = (topo_clean*az)**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1e-9)[0]
                except:
                    pars = x0
                #0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2
                sol[4] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]; sol[10] = pars[4]; sol[11] = pars[5]; sol[12] = pars[6]
                logger.info('Remove ramp %f az**2, %f az + %f + %f z + %f z**2 + %f z*az + %f (z*az)**2 '%(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5],pars[6]))

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,2] = 1
                G[:,3] = elev_map.flatten()
                G[:,4] = elev_map.flatten()**2
                G[:,5] = elev_map.flatten()
                G[:,6] = elev_map.flatten()**2
                for i in range(mlines):
                    G[i*mcols:(i+1)*mcols,0] = (i-ibeg)**2
                    G[i*mcols:(i+1)*mcols,1] = i -ibeg
                    G[:,5] *= i -ibeg
                    G[:,6] *= (i-ibeg)**2 


    corr = np.dot(G,pars).reshape(mlines,mcols)
    res = los - np.dot(G,pars).flatten()
    var = np.nanstd(res)
    # plt.imshow(los.reshape(mlines,mcols))
    # plt.show()
    # plt.imshow(corr)
    # plt.show()

    return sol, corr, var, rgbins, azbins, topobins, losbins, losstd

def empirical_cor(kk):
    """
    Function that preapare and run empirical estimaton for each interferogram kk
    """

    date1, date2 = date_1[kk], date_2[kk]
    idate = str(date1) + '-' + str(date2) 

    if sformat == 'ROI_PAC':
        folder =  'int_'+ str(date1) + '_' + str(date2) + '/'
        rscfile=int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.unw.rsc'
        infile=int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.unw'

        checkinfile(infile)
        checkinfile(rscfile)

        ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])
        # Get the band that have the data we want
        ds_band1 = ds.GetRasterBand(1)
        ds_band2 = ds.GetRasterBand(2)

        los_map = np.zeros((mlines,mcols))
        los_map[:ds.RasterYSize,:ds.RasterXSize] = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:mlines,:mcols]
        # los_map[los_map==0] = float('NaN')
        lines, cols = ds.RasterYSize, ds.RasterXSize


    elif sformat == 'GTIFF':
        infile = int_path + prefix + str(date1) + '_' + str(date2) + suffix + rlook + '.geo.unw.tif'

        checkinfile(infile)

        ds = gdal.Open(infile, gdal.GA_ReadOnly)
        # Get the band that have the data we want
        ds_band2 = ds.GetRasterBand(1)
        lines, cols = ds.RasterYSize, ds.RasterXSize

        los_map = np.zeros((mlines,mcols))
        los_map[:ds.RasterYSize,:ds.RasterXSize] = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:mlines,:mcols]
        # los_map[los_map==0] = float('NaN')

    elif sformat == 'GAMMA':
        # scfile=prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.unw.par'
        # par_file = ref 
        lines,cols = gm.readpar(int_path)
        infile= int_path + prefix + str(date1) + '_' + str(date2) + suffix + rlook + '.unw'
        checkinfile(infile)
        los_map = gm.readgamma(infile,int_path)

    logger.info('lines:{0}, cols:{1}, IFG:{2}'.format(lines, cols, idate))

    # load coherence or whatever
    spacial_mask = np.ones((mlines,mcols))*float('NaN')
    
    rms_map = np.ones((mlines,mcols))
    if rmsf=='yes':
        try:
            if sformat == 'ROI_PAC':
                rms_map[:ds.RasterYSize,:ds.RasterXSize] = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:mlines,:mcols]
                k = np.nonzero(np.logical_or(rms_map==0.0, rms_map==9999))
                rms_map[k] = float('NaN')
            elif sformat == 'GTIFF':
                rmsfile=  int_path + str(date1) + '_' + str(date2) + '.geo.cc.tif'
                checkinfile(rmsfile)
                ds = gdal.Open(rmsfile, gdal.GA_ReadOnly)
                ds_band1 = ds.GetRasterBand(1)
                lines, cols = ds.RasterYSize, ds.RasterXSize
                rms_map = np.zeros((mlines,mcols))
                rms_map[:ds.RasterYSize,:ds.RasterXSize] = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:mlines,:mcols]
                k = np.nonzero(np.logical_or(rms_map==0.0, rms_map==9999))
                rms_map[k] = float('NaN')
            elif sformat == 'GAMMA':
                rmsfile=  int_path + str(date1) + '_' + str(date2) + '.filt.cc'
                rms_map = gm.readgamma(rmsfile,int_path)
                # plt.imshow(rms_map,vmax=1,vmin=0)
                # plt.show()
        except:
           logger.warning('Coherence file cannot be read')

    # time.sleep(1.)
    # clean for estimation
    _los_map = np.copy(los_map)
    _los_map[los_map==0] = float('NaN')
    maxlos,minlos = np.nanpercentile(_los_map,perc),np.nanpercentile(_los_map,(100-perc))

    ## CRITICAL STEP ####
    # select points for estimation only: minmax elev, los not NaN, rms<rmsthreshold ....
    index = np.nonzero(
    np.logical_and(elev_map<maxelev,
    np.logical_and(elev_map>minelev,
    np.logical_and(slope_map>minslope,     
    np.logical_and(los_map!=0, 
    np.logical_and(los_map>minlos,
    np.logical_and(los_map<maxlos,
    np.logical_and(rms_map>threshold_rms, 
    np.logical_and(pix_az>ibeg,
    np.logical_and(pix_az<iend,
    np.logical_and(pix_rg>jbeg,
    np.logical_and(pix_rg<jend,
    np.logical_and(mask>threshold_mask,
    np.logical_and(~np.isnan(los_map),
    np.logical_or(pix_az<ibeg_mask,pix_az>iend_mask)
    ))))))))))))))
    
    indexref = np.nonzero(
    np.logical_and(elev_map<maxelev,
    np.logical_and(elev_map>minelev,    
    np.logical_and(los_map!=0, 
    np.logical_and(los_map>minlos,
    np.logical_and(los_map<maxlos,
    np.logical_and(rms_map>threshold_rms, 
    np.logical_and(pix_az>lin_start,
    np.logical_and(pix_az<lin_end,
    np.logical_and(pix_rg>col_start,
    np.logical_and(pix_rg<col_end,
    np.logical_and(mask>threshold_mask,
    np.logical_and(~np.isnan(los_map),
    np.logical_or(pix_az<ibeg_mask,pix_az>iend_mask)
    )))))))))))))

    spacial_mask[index] = np.copy(los_map[index])

    # extract range and azimuth coordinates
    temp = np.array(index).T
    az = temp[:,0]; rg = temp[:,1]

    # clean maps
    los_temp = np.matrix.copy(los_map)
    elev_temp = np.matrix.copy(elev_map)
    los_clean = los_temp[index].flatten()
    los_ref = los_temp[indexref].flatten()
    rms_ref = rms_map[indexref].flatten()
    cst = np.nansum(los_ref*rms_ref) / np.nansum(rms_ref)
    logger.info('Estimation of a constant within lines {0}-{1} and cols {2}-{3}'.format(lin_start,lin_end,col_start,col_end))
    logger.info('Average phase within ref area: {0}:'.format(cst))
    rg_ref, az_ref, topo_ref = np.nanmean(pix_rg[indexref]),np.nanmean(pix_az[indexref]),np.nanmean(elev_temp[indexref])


    # logger.info('Average rg: {0}, az:{1}, topo:{2}, within ref area'.format(np.int(rg_ref), np.int(az_ref), np.int(topo_ref)))
    # sys.exit()
    elev_clean = elev_temp[index].flatten()
    rms_clean = rms_map[index].flatten()
    del los_temp, elev_temp

    # Take care to not do high polynomial estimations for short int.
    # find the begining of the image
    itemp = ibeg
    for row in range(ibeg,iend,10):
      if np.isnan(np.nanmean(_los_map[row:row+10,:])):
          itemp = row  
      else:
          break
    del _los_map

    # 0: ref frame [default], 1: range ramp ax+b , 2: azimutal ramp ay+b, 
    # 3: ax+by+c, 4: ax+by+cxy+d 5: ax**2+bx+d, 6: ay**2+by+c
    if flat>5 and iend-itemp < .6*(iend-ibeg):
      logger.warning('Int. too short in comparison to ref, set flat to 5')
      temp_flat=5
    elif flat>5 and iend-itemp < .9*mcols:
      logger.warning('Lenght int. inferior to width, set flat to 5 and nfit to 0')
      temp_flat=5
    else:
      temp_flat=flat

    if ivar>0 and iend-itemp < .6*(iend-ibeg):
      logger.warning('Int. too short in comparison to ref, set ivar and nfit to 0')
      nfit_temp=0
      ivar_temp=0
    else:
      nfit_temp=nfit
      ivar_temp=ivar

    #try:
    sol, corr, var, rgbins, azbins, topobins, losbins, losstd = estim_ramp(los_map.flatten(),
        los_clean[::samp],elev_clean[::samp],az[::samp],rg[::samp],
        temp_flat,rms_clean[::samp],nfit_temp,ivar_temp,cst, rg_ref, az_ref, topo_ref)
    #except:
    #    sol = np.zeros((13))
    #    corr = np.zeros((mlines,mcols))
    #    var = 1
    #    rgbins, azbins = az[::samp],rg[::samp]
    #    topobins,losbins,losstd = elev_clean[::samp],los_clean[::samp],rms_clean[::samp]

    logger.info('RMS: {0} '.format(var))

    # 0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2
    func = sol[0]*rg**3 + sol[1]*rg**2 + sol[2]*rg + sol[3]*az**3 + sol[4]*az**2 \
    + sol[5]*az + sol[6]*(rg*az)**2 + sol[7]*rg*az + sol[11]*az*elev_clean + \
    sol[12]*((az*elev_clean)**2)

    if radar is not None: 
       # plot phase/elevation

       funcbins = sol[0]*rgbins**3 + sol[1]*rgbins**2 + sol[2]*rgbins + sol[3]*azbins**3 + sol[4]*azbins**2 \
       + sol[5]*azbins + sol[6]*(rgbins*azbins)**2 + sol[7]*rgbins*azbins + sol[11]*azbins*topobins + \
       sol[12]*((azbins*topobins)**2)

       fig2 = plt.figure(2,figsize=(9,4))
       ax = fig2.add_subplot(1,1,1)
       z = np.linspace(np.min(elev_clean), np.max(elev_clean), 100)
       # the sliding rms is only taken into account when flat =0 
       if ivar_temp == 0 and temp_flat == 0:
           ax.fill_between(topobins,losbins-losstd-funcbins,losbins+losstd-funcbins,color='red',alpha=0.1,label='sliding RMS')
       ax.scatter(elev_clean[::samp],los_clean[::samp] - func[::samp], s=0.005, alpha=0.05,rasterized=True)
       ax.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
       if nfit==0:
            ax.plot(z,sol[8]+sol[9]*z,'-r',lw =3.,label='{0:.3f}*z + {1:.3f}'.format(sol[9],sol[8])) 
       else:
            ax.plot(z,sol[8]+sol[9]*z+sol[10]*z**2, '-r', lw =3.,label='{0:.3f}*z**2 + {1:.3f}*z + {2:.3f}'.format(sol[10],sol[9],sol[8]))

       ax.set_xlabel('Elevation (m)')
       ax.set_ylabel('LOS (rad)')
       plt.legend(loc='best')
       if sformat == 'ROI_PAC':
          fig2.savefig( int_path + folder + idate+'phase-topo.png', format='PNG')
       else:
          fig2.savefig(out_path + idate+'phase-topo.png', format='PNG')

    _los_map = np.copy(los_map)
    _los_map[los_map==0] = float('NaN')
    vmax = np.nanpercentile(_los_map,98)
    vmin = np.nanpercentile(_los_map,2)

    fig = plt.figure(3,figsize=(11,4))

    ax = fig.add_subplot(2,3,1)
    cax = ax.imshow(los_map,cmap=cmap,vmax=vmax,vmin=vmin,interpolation=None)
    ax.set_title('LOS')
    plt.setp( ax.get_xticklabels(), visible=None)
    plt.setp( ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    ax = fig.add_subplot(2,3,2)
    cax = ax.imshow(spacial_mask,cmap=cmap,vmax=vmax,vmin=vmin,interpolation=None)
    ax.set_title('LOS ESTIMATION')
    plt.setp( ax.get_xticklabels(), visible=None)
    plt.setp( ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)
    
    ax = fig.add_subplot(2,3,3)
    cax = ax.imshow(rms_map,cmap=cmap,interpolation=None)
    ax.set_title('RMS')
    plt.setp( ax.get_xticklabels(), visible=None)
    plt.setp( ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    # colormap correction must be same than data!!
    ax = fig.add_subplot(2,3,4)
    cax = ax.imshow(corr,cmap=cmap,vmax=vmax,vmin=vmin,interpolation=None)
    ax.set_title('RAMP+TOPO')
    plt.setp( ax.get_xticklabels(), visible=None)
    plt.setp( ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    # for plot we can clean
    k = np.nonzero(np.logical_or(los_map==0.,abs(los_map)>999.))
    corr[k] = 0.

    flatlos = los_map - corr
    _los_map = np.copy(flatlos)
    _los_map[los_map==0] = float('NaN')
    vmax = np.nanpercentile(_los_map,98)
    vmin = np.nanpercentile(_los_map,2)

    ax = fig.add_subplot(2,3,5)
    hax = ax.imshow(rms_map, cm.Greys,vmax=1,vmin=0.)
    cax = ax.imshow(flatlos,cmap=cmap,vmax=vmax,vmin=-vmax,alpha=1.,interpolation=None)
    ax.set_title('CORR LOS')
    plt.setp( ax.get_xticklabels(), visible=None)
    plt.setp( ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)
    fig.tight_layout()

    if sformat == 'ROI_PAC':
        fig.savefig(int_path + folder + idate +'corrections.png', format='PNG')
    else:
        fig.savefig(out_path + idate +'corrections.png', format='PNG')

    if plot=='yes':
        plt.show()

    plt.close('all')
    del corr, flatlos

    del los_clean, rms_clean
    del elev_clean
    del az, rg
    try:
        del ds
    except:
        pass

    return iend-itemp, sol, var

def apply_cor(kk, sp, sp_inv):
    """
    Fonction that apply empirical estimatons for each interferograms kk
    """

    date1, date2 = date_1[kk], date_2[kk]
    folder = 'int_'+ str(date1) + '_' + str(date2) + '/'
    idate = str(date1) + '-' + str(date2) 

    if sformat == 'ROI_PAC':
        
        rscfile=int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix +  rlook + '.unw.rsc'
        infile=int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix +  rlook + '.unw'
        outfile = int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix +  suffout +  rlook + '.unw'  
        outrsc = int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix +  suffout +  rlook + '.unw.rsc' 
        # print(infile)
        # print(outfile)
        
        ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])
        # Get the band that have the data we want
        ds_band1 = ds.GetRasterBand(1)
        ds_band2 = ds.GetRasterBand(2)
        los_map = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        rms_map = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        lines, cols = ds.RasterYSize, ds.RasterXSize

    elif sformat == 'GTIFF':

        infile = int_path + prefix + str(date1) + '_' + str(date2) + suffix + rlook + '.geo.unw.tif'
        outfile = out_path + prefix + str(date1) + '_' + str(date2) + suffix + suffout +  rlook + '.geo.unw.tif' 
        ds = gdal.Open(infile, gdal.GA_ReadOnly)
        ds_band2 = ds.GetRasterBand(1)
        los_map = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        rms_map = np.ones((ds.RasterYSize,ds.RasterXSize))
        lines, cols = ds.RasterYSize, ds.RasterXSize

        if rmsf == 'yes':
          rmsfile=  int_path + str(date1) + '_' + str(date2) + 'geo.cc.tif'
          rmsfile=  int_path + str(date1) + '_' + str(date2) + '.geo.cc.tif'
          ds = gdal.Open(rmsfile, gdal.GA_ReadOnly)
          ds_band1 = ds.GetRasterBand(1)
          rms_map = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        else:
          rms_map = np.ones((lines,cols))

    elif sformat == 'GAMMA':

        infile = int_path + prefix + str(date1) + '_' + str(date2) + suffix +  rlook + '.unw'
        outfile = out_path + prefix + str(date1) + '_' + str(date2) + suffix +  suffout + rlook + '.unw' 
        # par_file = ref 
        lines,cols = gm.readpar(int_path)
        los_map = gm.readgamma(infile,int_path)
        
        if rmsf == 'yes':
          rmsfile=  int_path + str(date1) + '_' + str(date2) + '.filt.cc'
          rms_map = gm.readgamma(rmsfile,int_path)
        else:
          rms_map = np.ones((lines,cols))

    logger.info('mlines:{}, mcols:{}, int:{}:'.format(lines, cols, idate))

    # compute correction
    rg = np.tile(np.arange(cols), (lines,1))
    az = np.tile(np.arange(lines), (cols,1)).T
    z = np.zeros((lines, cols))
    z = elev_map[:lines,:cols]

    # 0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2

    # apply correction
    sol = sp_inv[kk,3:]
    corr_inv = sol[0]*rg**3 + sol[1]*rg**2 + sol[2]*rg + sol[3]*az**3 + sol[4]*az**2 \
    + sol[5]*az + sol[6]*(rg*az)**2 + sol[7]*rg*az + sol[8] + sol[9]*z + sol[10]*(z**2) \
    + sol[11]*az*z + sol[12]*((az*z)**2)

    sol = sp[kk,3:]
    corr = sol[0]*rg**3 + sol[1]*rg**2 + sol[2]*rg + sol[3]*az**3 + sol[4]*az**2 \
    + sol[5]*az + sol[6]*(rg*az)**2 + sol[7]*rg*az + sol[8] + sol[9]*z + sol[10]*(z**2) \
    + sol[11]*az*z + sol[12]*((az*z)**2)
            
    # apply corr
    flatlos = los_map - corr_inv
    
    # recompute ref frame
    zone = as_strided(flatlos[lin_start:lin_end,col_start:col_end])
    amp = as_strided(rms_map[lin_start:lin_end,col_start:col_end])
    minlos, maxlos = np.nanpercentile(zone[abs(zone)>1.e-6],5.), np.nanpercentile(zone[abs(zone)>1.e-6],95.)
    index = np.nonzero(
        np.logical_and(zone>minlos,
        np.logical_and(zone<maxlos,
        amp>threshold_rms,  
    )))
    # give small weights to small uncoherent vectors
    cst = np.nansum(zone[index]*amp[index]) / np.nansum(amp[index])
    if np.isnan(cst):
        pass
    else:
        flatlos = flatlos - cst
        corr_inv = corr_inv + cst
        corr = corr + cst
    logger.info('Remove reference frame: {}'.format(cst))

    # reset to 0 areas where no data (might change after time series inversion?)
    flatlos[np.isnan(los_map)],rms_map[np.isnan(los_map)] = 0.0, 0.0
    flatlos[np.isnan(flatlos)],rms_map[np.isnan(flatlos)] = 0.0, 0.0
    rms_map[np.isnan(rms_map)],flatlos[np.isnan(rms_map)] = 0.0, 0.0
    flatlos[los_map==0], rms_map[los_map==0] = 0.0, 0.0

    if sformat == 'ROI_PAC':
        dst_ds = driver.Create(outfile, cols, lines, 2, gdal.GDT_Float32)
        dst_band1 = dst_ds.GetRasterBand(1)
        dst_band2 = dst_ds.GetRasterBand(2)
        dst_band1.WriteArray(rms_map,0,0)
        dst_band2.WriteArray(flatlos,0,0)
        shutil.copy(rscfile,outrsc)
        dst_band1.FlushCache()
        dst_band2.FlushCache()
        del dst_ds, ds

    elif sformat == 'GTIFF':
        dst_ds = driver.Create(outfile, cols, lines, 1, gdal.GDT_Float32)
        dst_band2 = dst_ds.GetRasterBand(1)
        dst_band2.WriteArray(flatlos,0,0)
        dst_ds.SetGeoTransform(gt)
        dst_ds.SetProjection(proj)
        dst_band2.FlushCache()
        del dst_ds, ds

    elif sformat == 'GAMMA':
        fid = open(outfile, 'wb')
        flatlos.flatten().astype('>f4').tofile(fid)
        fid.close()

    fig = plt.figure(5,figsize=(9,6))

    _los_map = np.copy(los_map)
    _los_map[los_map==0] = float('NaN')
    vmax, vmin = np.nanpercentile(_los_map,99), np.nanpercentile(_los_map,1)
    
    ax = fig.add_subplot(1,4,1)
    cax = ax.imshow(los_map,cmap=cmap,vmax=vmax,vmin=vmin,alpha=1,interpolation=None)
    ax.set_title(str(date1) + '_' + str(date2))
    plt.setp( ax.get_xticklabels(), visible=None)
    plt.setp( ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    ax = fig.add_subplot(1,4,2)
    cax = ax.imshow(corr,cmap=cmap,vmax=vmax,vmin=vmin,alpha=1,interpolation=None)
    ax.set_title('Best-fit Ramp')
    plt.setp( ax.get_xticklabels(), visible=None)
    plt.setp( ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    ax = fig.add_subplot(1,4,3)
    cax = ax.imshow(corr_inv,cmap=cmap,vmax=vmax,vmin=vmin,alpha=1,interpolation=None)
    ax.set_title('Reconstructed')
    plt.setp( ax.get_xticklabels(), visible=None)
    plt.setp( ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    _los_map = np.copy(flatlos)
    _los_map[_los_map==0] = float('NaN')
    #flatlos[los_map==0] = float('NaN')
    vmax, vmin = np.nanpercentile(_los_map,99), np.nanpercentile(_los_map,1)

    ax = fig.add_subplot(1,4,4)
    cax = ax.imshow(_los_map,cmap=cmap,vmax=vmax,vmin=vmin,alpha=1,interpolation=None)
    ax.set_title('Flattened LOS')
    plt.setp( ax.get_xticklabels(), visible=None)
    plt.setp( ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)
    fig.tight_layout()

    if sformat == 'ROI_PAC':
        fig.savefig( int_path + folder + idate + '_reconstruc_corrections.png', format='PNG')
    else:
        fig.savefig( out_path + idate + '_reconstruc_corrections.png', format='PNG')

    if plot=='yes':
        plt.show()

    plt.close('all')
    del los_map, rms_map

#####################################################################################
# INIT LOG
#####################################################################################

# logging.basicConfig(level=logging.INFO,\
logging.basicConfig(level=logging.INFO,\
        format='line %(lineno)s -- %(levelname)s -- %(message)s')
logger = logging.getLogger('invert_ramp_topo_unw.log')

#####################################################################################
# READ INPUT PARAM
#####################################################################################

# read arguments
arguments = docopt.docopt(__doc__)

int_list=arguments["--int_list"]

if arguments["--int_path"] == None:
    int_path='./'
else:
    int_path=arguments["--int_path"] 

if arguments["--out_path"] == None:
    out_path=np.copy(int_path)
else:
    out_path=arguments["--out_path"] + '/'
    makedirs(out_path)

if arguments["--prefix"] == None:
    prefix = ''
else:
    prefix=arguments["--prefix"]
if arguments["--suffix"] == None:
    suffix = ''
else:
    suffix=arguments["--suffix"]
if arguments["--rlook"] == None:
    rlook = ''
else:
    rlook = '_' + arguments["--rlook"] + 'rlks'

if arguments["--flat"] == None:
    flat = 0
elif int(arguments["--flat"]) <  7:
    flat = int(arguments["--flat"])
else:
    logger.warning('flat > 6, set flat to 0')
    flat = 0

if arguments["--ivar"] == None:
    ivar = 0
elif int(arguments["--ivar"]) <  2:
    ivar = int(arguments["--ivar"])
    if arguments["--topofile"] is None:
        logger.critical('No topographic file is given. Empirical phase/topo will not be performed')

else:
    logger.warning('ivar > 1, set ivar to 0')
    ivar = 0
    if arguments["--topofile"] is None:
        logger.critical('No topographic file is given. Empirical phase/topo will not be performed')


if arguments["--nfit"] == None:
    nfit = 0
elif int(arguments["--nfit"]) <  2:
    nfit = int(arguments["--nfit"])
    if arguments["--topofile"] is None:
        logger.critical('No topographic file is given. Empirical phase/topo will not be performed')

else:
    logger.warning('nfit > 1, set nfit to 0')
    nfit = 0
    if arguments["--topofile"] is None:
        logger.critical('No topographic file is given. Empirical phase/topo will not be performed')


if arguments["--mask"] ==  None or not os.path.exists(arguments["--mask"]):
    maskfile = None
else:
    maskfile = arguments["--mask"]
if arguments["--threshold_mask"] ==  None:
    threshold_mask = -1
else:
    threshold_mask = float(arguments["--threshold_mask"])
if arguments["--threshold_coh"] ==  None:
    threshold_rms = -1
else:
    threshold_rms = float(arguments["--threshold_coh"])
if arguments["--cohpixel"] ==  None:
    rmsf = 'no'
    threshold_rms = -1
else:
    rmsf = arguments["--cohpixel"]
if arguments["--topofile"] ==  None or not os.path.exists(arguments["--topofile"]):
   radar = None
else:
   radar = arguments["--topofile"]

if arguments["--ref"] ==  None :
   ref = None
else:
   ref = arguments["--ref"]
if ref == None and radar == None:
    logger.critical('Argument error: Need to give ref or topographic file')
    sys.exit()

if arguments["--tsinv"] ==  None:
    tsinv = 'no'
else:
    tsinv = arguments["--tsinv"]

if arguments["--format"] ==  None:
    sformat = 'ROI_PAC'
else:
    sformat = arguments["--format"]

if arguments["--estim"] ==  None:
    estim = 'yes'
else:
    estim = arguments["--estim"]

if arguments["--ibeg_mask"] ==  None:
    ibeg_mask = np.inf
else:
    ibeg_mask = int(arguments["--ibeg_mask"])
if arguments["--iend_mask"] ==  None:
    iend_mask = -np.inf
else:
    iend_mask = int(arguments["--iend_mask"])
if arguments["--perc_topo"] ==  None:
    perc_topo = 99.
else:
    perc_topo = float(arguments["--perc_topo"])
if arguments["--perc"] ==  None:
    perc = 99.
else:
    perc = float(arguments["--perc"])
if arguments["--perc_slope"] ==  None:
    perc_slope = 92.
else:
    perc_slope = float(arguments["--perc_slope"])

if arguments["--nproc"] ==  None:
    nproc = 4
else:
    nproc = int(arguments["--nproc"])

if arguments["--plot"] ==  'yes':
    plot = 'yes'
    logger.warning('plot is yes. Set nproc to 1')
    nproc = 1
    if environ["TERM"].startswith("screen"):
        matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt
else:
    plot = 'no'
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt

if arguments["--suffix_output"] ==  None:
    suffout = '_corrunw'
else:
    suffout = arguments["--suffix_output"]

if arguments["--samp"] == None:
    samp = 2
else:
    samp = int(arguments["--samp"])

if arguments["--cpt"] is  None:
    #cmap=cm.rainbow
    ### Load colormaps
    try:
      cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
      cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt")).reversed()
    except:
      cmap=cm.rainbow
else:  
    try:
        cmap = LinearSegmentedColormap.from_list(arguments["--cpt"].split("/")[-1].split('.')[0], np.loadtxt(arguments["--cpt"]))
    except:
        cmap=arguments["--cpt"]

# print(nproc, plot)
# sys.exit()

#####################################################################################
# INITIALISE 
#####################################################################################

# print()
# read int
date_1,date_2=np.loadtxt(int_list,comments="#",unpack=True,usecols=(0,1),dtype='i,i')
Nifg=len(date_1)
logger.info("number of interferogram: {}".format(Nifg))

# list dates
im = []; bt = []
for date1,date2 in zip(date_1,date_2):
    if date1 not in im: im.append(date1)
    if date2 not in im: im.append(date2)
nmax=len(im)
logger.info("number of image: {} ".format(nmax))
imd = date2dec(im)
cst = np.copy(imd[0])
# compute temporal baseline for TS inv
for i in range((nmax)):
    bt.append(imd[i]-cst)

# # Now, write list_pair
# wf = open("list_dates", "w")
# for i in range((nmax)):
#     wf.write("%i %.6f %.6f\n" % (im[i], imd[i], bt[i]))
# wf.close()
# sys.exit()

# load ref to define mlines, mcols and format
if ref is not None:
    ds_extension = os.path.splitext(ref)[1]
    if sformat == 'GTIFF':
        geotiff = arguments["--ref"]
        georef = gdal.Open(ref)
        gt = georef.GetGeoTransform()
        proj = georef.GetProjection()
        driver = gdal.GetDriverByName('GTiff')
        ds = gdal.Open(ref, gdal.GA_ReadOnly)
        mlines,mcols = ds.RasterYSize, ds.RasterXSize
    elif sformat == 'ROI_PAC':
        driver = gdal.GetDriverByName("roi_pac")
        ds = gdal.OpenEx(ref, allowed_drivers=["ROI_PAC"])
        mlines,mcols = ds.RasterYSize, ds.RasterXSize
    elif sformat == 'GAMMA':
        from parsers import gamma as gm
        # par_file = ref 
        mlines,mcols = gm.readpar(int_path)

# laod elevation map
if radar is not None:
    extension = os.path.splitext(radar)[1]
    if extension == '.r4':  
      fid = open(radar,'r')
      elevi = np.fromfile(fid,dtype=np.float32)
      elev_map = elevi.reshape(mlines,mcols)
      fid.close()

    else:
        if sformat == 'GTIFF':
            geotiff = arguments["--geotiff"]
            georef = gdal.Open(radar)
            gt = georef.GetGeoTransform()
            proj = georef.GetProjection()
            driver = gdal.GetDriverByName('GTiff')
            ds = gdal.Open(radar, gdal.GA_ReadOnly)
            ds_band2 = ds.GetRasterBand(2)
            mlines,mcols = ds.RasterYSize, ds.RasterXSize
            elev_map = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
            del ds

        elif sformat == 'ROI_PAC':
            driver = gdal.GetDriverByName("roi_pac")
            ds = gdal.OpenEx(radar, allowed_drivers=["ROI_PAC"])
            ds_band2 = ds.GetRasterBand(2)
            mlines,mcols = ds.RasterYSize, ds.RasterXSize
            elev_map = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
            del ds
        
        elif sformat == 'GAMMA':
            from parsers import gamma as gm
            # par_file = ref 
            mlines,mcols = gm.readpar()
            elev_map = gm.readgamma(radar)
        
        if arguments["--max_topo"] !=  None:
            maxelev = float(arguments["--max_topo"])
        else:
            maxelev = np.nanpercentile(elev_map,perc_topo)
        if arguments["--min_topo"] !=  None:
            minelev = float(arguments["--min_topo"])
        else:
            minelev = np.nanpercentile(elev_map,100-perc_topo)
        logger.info('Max-Min topography for empirical estimation: {0:.1f}-{1:.1f}'.format(maxelev,minelev))

    # compute slope
    toposmooth = scipy.ndimage.filters.gaussian_filter(elev_map,3.)
    Py, Px = np.gradient(toposmooth)
    slope_map = np.sqrt(Px**2+Py**2)
    minslope = np.nanpercentile(slope_map,100-perc_slope)
    logger.info('Min relief for empirical estimation: {0}'.format(minslope))
    
    fig = plt.figure(0,figsize=(12,8))

    ax = fig.add_subplot(1,2,1)
    cax = ax.imshow(toposmooth, cm.RdBu_r, vmin=minelev, vmax=maxelev)
    plt.setp( ax.get_xticklabels(), visible=False)
    ax.set_title('Smoothed DEM',fontsize=6)

    ax = fig.add_subplot(1,2,2)
    cax = ax.imshow(slope_map, cm.RdBu_r, vmin=minslope, vmax=np.nanpercentile(slope_map,perc_slope))
    plt.setp( ax.get_xticklabels(), visible=False)
    ax.set_title('Mask Slope bellow: {0:.2f}'.format(minslope),fontsize=8)

    if plot == 'yes':
        plt.show()
 
else:
    maxelev,minelev = 1.,-1
    elev_map = np.zeros((mlines,mcols))
    slope_map = np.zeros((mlines,mcols))
    minslope = -1

# open mask file
if maskfile is not None:
    fid = open(maskfile,'r')
    maski = np.fromfile(fid,dtype=np.float32)[:mcols*mlines]
    mask = maski.reshape((mlines,mcols))
    k = np.nonzero(mask<threshold_mask)
    spacial_mask = np.copy(mask)
    spacial_mask[k] = float('NaN')

    if plot=='yes':

      fig = plt.figure(0,figsize=(5,4))
      ax = fig.add_subplot(1,1,1)
      cax = ax.imshow(spacial_mask,cmap=cmap)
      ax.set_title('Mask')
      plt.setp( ax.get_xticklabels(), visible=None)
      divider = make_axes_locatable(ax)
      c = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(cax, cax=c)
      plt.show()
      fid.close()
else:
    mask = np.zeros((mlines,mcols))
    threshold_mask = -1

if arguments["--crop_emp"] ==  None:
    crop_emp = [0,mlines,0,mcols]
else:
    crop_emp = list(map(float,arguments["--crop_emp"].replace(',',' ').split()))
    logger.warning('Crop empirical estimation between lines {}-{} and cols:{}'.format(int(crop_emp[0]),int(crop_emp[1]),int(crop_emp[2]),int(crop_emp[3])))
ibeg,iend,jbeg,jend = int(crop_emp[0]),int(crop_emp[1]),int(crop_emp[2]),int(crop_emp[3])

if arguments["--ref_zone"] == None:
    lin_start, lin_end, col_start, col_end = 0,mlines,0,mcols
else:
    ref = list(map(int,arguments["--ref_zone"].replace(',',' ').split()))
    try:
        lin_start,lin_end, col_start, col_end = ref[0], ref[1], ref[2], ref[3]
    except:
        lin_start,lin_end = ref[0], ref[1]
        col_start, col_end = 0, mcols

#####################################################################################
# MAIN
#####################################################################################

# extract range and azimuth coordinates from ref or radar file
pix_az, pix_rg = np.indices((mlines,mcols))

# initialise full vector correction 
# 16 values: date1 dates2  mlines y**3 y**2 y x**3 x**2 x xy**2 xy cst z z**2 z*az az*z**2 
spint = np.zeros((Nifg,16))
M = 13 # harcoding of the number of cols....

# fill dates
spint[:,0],spint[:,1] = date_1, date_2 

# initilise correction cube
rmsint = np.zeros((Nifg,3))
rmsint[:,0],rmsint[:,1] = date_1,date_2 

date1_err, date2_err = [], []

if estim=='yes':

    print() 
    #########################################
    print('#################################')
    print('Empirical estimations')
    print('#################################')
    #########################################
    print()

    output = []
    # go 
    with TimeIt():
        # for kk in range(Nifg):
        work = range(Nifg)
        with poolcontext(processes=nproc) as pool:
            results = pool.map(empirical_cor, work)
        output.append(results)
        # for work in range(Nifg):
        #     empirical_cor(work)

        for kk in range(Nifg):
            lenght, sol, rms = output[0][kk]
            # save size int to use as weight in the temporal inversion
            spint[kk,2] = lenght
            # fill correction matrix
            spint[kk,3:] = sol
            rmsint[kk,2] = rms

    # print(spint)

    # save spint
    if sformat == 'ROI_PAC':
      np.savetxt('list_coeff_ramps{}.txt'.format(suffout), spint , header='#date1   |   dates2   |   Lenght   |   y**3   |   y**2   |   y\
       |   **3   |   x**2   |   x   |   xy**2   |   xy   |   cst   |   z   |   z**2   |   z*az   |   z**2*az', fmt=('%i','%i','%.10f','%.10f','%.10f','%.10f','%.10f',\
        '%.10f','%.10f','%.10f','%.10f','%.10f','%.10f','%.10f','%.10f','%.10f'))

    else:
      np.savetxt(out_path+'list_coeff_ramps{}.txt'.format(suffout), spint , header='#date1   |   dates2   |   Lenght   |   y**3   |   y**2   |   y\
       |   **3   |   x**2   |   x   |   xy**2   |   xy   |   cst   |   z   |   z**2   |   z*az   |   z**2*az', fmt=('%i','%i','%.10f','%.10f','%.10f','%.10f','%.10f',\
        '%.10f','%.10f','%.10f','%.10f','%.10f','%.10f','%.10f','%.10f','%.10f'))

    # save rms
    if sformat == 'ROI_PAC':
        np.savetxt('rms{}.txt'.format(suffout), rmsint, header='# date1   |   dates2   |   RMS', fmt=('%i','%i','%.8f'))
    else:
        np.savetxt(out_path+'rms{}.txt'.format(suffout), rmsint, header='# date1   |   dates2   |   RMS', fmt=('%i','%i','%.8f'))

#####################################################################################

print()
logger.info('read input files list_coeff_ramps.txt')

# load spint
if sformat == 'ROI_PAC':
    date_1,date_2,length,a,b,c,d,e,f,g,h,i,j,k,l,m=np.loadtxt('list_coeff_ramps{}.txt'.format(suffout),comments="#",unpack=True,dtype='i,i,f,f,f,f,f,f,f,f,f,f,f,f,f,f')
else:
    date_1,date_2,length,a,b,c,d,e,f,g,h,i,j,k,l,m=np.loadtxt(out_path+'list_coeff_ramps{}.txt'.format(suffout),comments="#",unpack=True,dtype='i,i,f,f,f,f,f,f,f,f,f,f,f,f,f,f')
spint = np.vstack([date_1,date_2,length,a,b,c,d,e,f,g,h,i,j,k,l,m]).T
rec_spint = np.copy(spint)

# load rms
if sformat == 'ROI_PAC':
    rmsint = np.loadtxt('rms{}.txt'.format(suffout),comments="#")
else:
    rmsint = np.loadtxt(out_path+'rms{}.txt'.format(suffout),comments="#")

####################################################################################

# init inv solution
spint_inv = np.zeros((np.shape(spint)))

if tsinv=='yes':

    print()
    #########################################
    print('#################################')
    print('Temporal inversion of all coefficients')
    print('#################################')
    #########################################
    print()

    G_=np.zeros((Nifg,nmax))
    deltat = np.zeros((Nifg))
    for k in range((Nifg)):
      for n in range((nmax)):
        if (date_1[k]==im[n]): 
          G_[k,n]=-1
          t1 = bt[n]
        elif (date_2[k]==im[n]):
          G_[k,n]=1
          t2 = bt[n]
      deltat[k] = abs(t2 -t1)

    # 1) create weight based on temporal baseline: give stronger weight to short temporal baselines 
    #, where we dont expect def.
    w1 = np.exp(-(deltat/2.))

    # 2) create a weight based on the size of the int: give a stronger weight to long interferograms
    w2 = length/mlines

    # 3) rms weight
    w3 =  np.exp(-rmsint[:,2]/np.nanpercentile(rmsint[:,2],80)) + 0.01
  
    if sformat == 'ROI_PAC':
        wf = open('list_weigths{}.txt'.format(suffout),"w")
    else:
        wf = open(out_path+'list_weigths{}.txt'.format(suffout),"w")
    print('Weights for temporal inversion:') 
    for i in range(len(w3)):
        print(int(rmsint[i,0]), int(rmsint[i,1]) , rmsint[i,2], w1[i], w2[i], w3[i])
        wf.write("%i %i %f %f %f %f\n" % (int(rmsint[i,0]), int(rmsint[i,1]) , rmsint[i,2], w1[i], w2[i], w3[i]))
    wf.close() 
    print('Weights saved in: list_weights.txt') 
    
    # compute summ of weights
    sig_ = 1./w1 + 1./w2 + 1./w3 

    for j in range(3,np.shape(spint_inv)[1]):

        d = np.zeros(((Nifg+1)))
        sig = np.ones(((Nifg+1)))
        G = np.zeros(((Nifg+1),nmax))
    
        d[:Nifg] = as_strided(spint[:,j])
        G[:Nifg,:nmax] = G_ 
        G[-1,0] = 1 # ini phi first image to 0
        sig[:Nifg] = sig_

        try:
            x0 = lst.lstsq(G,d)[0]
            _func = lambda x: np.sum(((np.dot(G,x)-d)/sig)**2)
            _fprime = lambda x: 2*np.dot(G.T/sig, (np.dot(G,x)-d)/sig)
            pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=20000,full_output=True,iprint=0)[0]

            # reconstruct corr for selected int
            spint_inv[:,j] = np.dot(G,pars)[:Nifg]

        except:
            pass

    spint_inv[:,:3] = spint[:,:3]
    # save spint_inv
    if sformat == 'ROI_PAC':
        np.savetxt('list_coeff_ramps{}_inv.txt'.format(suffout), spint_inv , header='#date1   |   dates2   |   Lenght   |   y**3   |   y**2   |   y\
       |   **3   |   x**2   |   x   |   xy**2   |   xy   |   cst   |   z   |   z**2   |   z*az   |   z**2*az', fmt=('%i','%i','%.10f','%.10f','%.10f','%.10f','%.10f',\
        '%.10f','%.10f','%.10f','%.10f','%.10f','%.10f','%.12f','%.12f','%.12f')) 
    else:
        np.savetxt(out_path+'list_coeff_ramps{}_inv.txt'.format(suffout), spint_inv , header='#date1   |   dates2   |   Lenght   |   y**3   |   y**2   |   y\
       |   **3   |   x**2   |   x   |   xy**2   |   xy   |   cst   |   z   |   z**2   |   z*az   |   z**2*az', fmt=('%i','%i','%.10f','%.10f','%.10f','%.10f','%.10f',\
        '%.10f','%.10f','%.10f','%.10f','%.10f','%.10f','%.12f','%.12f','%.12f'))
####################################################################################

print() 
#########################################
print('#################################')
print('APPLY CORRECTION AND SAVE NEW INT.')
print('#################################')
#########################################
print()

# go 
with TimeIt():
    work = range(Nifg)
    # for j in range(Nifg):
    #     apply_cor(j, spint, spint_inv)
    with poolcontext(processes=nproc) as pool:
        pool.map(partial(apply_cor, sp=spint, sp_inv=spint_inv), work)
