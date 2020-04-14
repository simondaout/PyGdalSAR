#!/usr/bin/env python3

"""
nsb_unw_ERAcorr_stats.py
-------------------
Reads in ERA corrections in radar slant geometry. Following this apply
correction to IFGs, using IFG=ERA+RAMP. Enforce closure in the
ramp corrections, and save these resulting IFGs for timeseries 
analysis. In addition perform statistical analysis of corrections.


.. Author:
    Nicholas Dodds, University of Oxford; November, 2019.
.. Last Modified:
    26th January 2020

Usage:
    nsb_unw_ERAcorr_stats.py --int_dir=<path> --int_prefix=<string> --int_suffix=<string> --era_dir=<path> --int_list=<path> 
    --rdr_rsc=<path> --rdr_dem=<path> --ref_zone=<xstart,xend,ystart,yend> [--int_looks=<value>] [--nproc=<nb_cores>] [--plot=<yes/no>] [estim=<yes/no>]

Options:
    -h --help                   Show this screen
    --int_dir=PATH              Folder containing unw IFGs in radar.
    --int_prefix=STRING         Prefix string on unw IFGs.
    --int_suffix=STRING         Suffix string on unw IFGs.
    --era_dir=PATH              Folder containing era corrections in radar.
    --int_list=PATH             List of IFGs for TS processing/correction.
    --rdr_rsc=PATH              Path to radar param file
    --rdr_dem=PATH              Path to DEM in radar geometry.
    --ref_zone=VALUES           Area where phase is set to zero for double difference subtraction (Y0,Y1,X0,X1).
    --int_looks=VALUE           Looks on unwrapped IFG (default = 2rlks)
    --nproc=NB_CORES            Number of cores to use for parallelizing (default=1)
    --plot=YES/NO               If yes, plot figures for each ints [default: no]
    --estim=YES/NO              Estimate IFG = ATMOS + RAMP, must run first time (default=yes)

"""

import os, sys, logging, glob
import docopt
import numpy as np
import matplotlib.pyplot as plt 
from os import path, environ
import seaborn as sns
from copy import deepcopy

from scipy import ndimage
import scipy.constants as const
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.stats import pearsonr

import scipy.optimize as opt
import scipy.linalg as lst

from shutil import copyfile
from multiprocessing import Process, Pool
import subprocess

import matplotlib as mpl
# from matplotlib import pyplot as plt
import matplotlib.cm as cm
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import matplotlib.ticker as plticker
from matplotlib.colorbar import Colorbar
import matplotlib.patches as patches
import matplotlib.ticker as ticker

import gdal
gdal.UseExceptions()

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning) 

### Load colormaps
cm_locs = '/home/comethome/jdd/ScientificColourMaps5/by_platform/python/'
cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt")).reversed()

### Functions

def rdr_rsc(inname, full=False, verbose=False):
        '''Reading a ROI-PAC style radar rsc file.

        Args:
                * inname (str): Path to the RSC file.

        Returns:
                * nx  (np.int)   : Number of range bins.
                * ny  (np.int)   : Number of azimuth bins.

                .. note:: 
                        Currently set up to work with dem.rsc file from ROI-PAC.
            '''

        if verbose:
                logger.info("PROGRESS: READING %s RSC FILE" %inname)

        rpacdict = {}
        infile = open(inname,'r')
        line = infile.readline()
        while line:
                llist = line.split()
                if len(llist)>0 :
                        rpacdict[llist[0]] = llist[1]
                line = infile.readline()
        infile.close()

        nx = np.int(rpacdict['WIDTH'])
        ny = np.int(rpacdict['FILE_LENGTH'])

        if full:
                return nx,ny,rpacdict
        else:
                return nx,ny

def ramp_fit_int_era(rms_map, i_data, e_data):
    '''
    Returns ramp fit that ties IFG data to ERA:
    IFG - ERA = ax + by + c
    a = pars[0]; b = pars[1]; c = pars[2]
    '''
    # Clean NaNs from data, and pixels with ~0 colinearity
    index = np.nonzero(
        np.logical_and(~np.isnan(i_data),
        np.logical_and(~np.isnan(e_data), rms_map>1.e-6)
        ))
    
    temp = np.array(index).T
    x_clean = temp[:,0]
    y_clean = temp[:,1]
    del temp
    
    d = i_data[index] - e_data[index]
    
    # Clean the rms_map too
    rms_clean = rms_map[index].flatten()
    sigmad = rms_clean
    
    # build inverse matrix
    G=np.zeros((len(d),3))
    G[:,0] = x_clean
    G[:,1] = y_clean
    G[:,2] = 1
    
    # calculate lsq as a prior to gradient optimisation below
    x0 = lst.lstsq(G, d)[0]
    
    # gradient optimization method, find ramp that best fits data (ifg - era)
    _func = lambda x: np.sum(((np.dot(G,x)-d)/sigmad)**2)
    _fprime = lambda x: 2*np.dot(G.T/sigmad, (np.dot(G,x)-d)/sigmad)
    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=500,full_output=True,iprint=0)[0]
    
    a = pars[0]; b = pars[1]; c = pars[2]
    print ('Remove ramp  %f az  + %f r + %f.'%(a, b, c))
    
    return pars

def sliding_median(x_data, y_data):
    # Generate bins and bin data
    total_bins=int(round(y_data.shape[0]/5000))
    bins = np.linspace(x_data.min(), x_data.max(), total_bins)
    delta = bins[1] - bins[0]
    idx  = np.digitize(x_data, bins)
    # Calculate sliding median
    running_median = np.array([np.median(y_data[idx==k]) for k in range(total_bins)])
    # Calculate associated sliding standard deviation
    running_std = np.array([y_data[idx==k].std() for k in range(total_bins)])
    # Calculate bin centres for plotting
    binned_xdata = np.array(bins-delta/2)
    # return bin centres, median, and std
    return binned_xdata, running_median, running_std

def correct_int_era_ramp(m_date, s_date):
    print('Estimating IFG - ERA = RAMP for: {}_{}'.format(m_date, s_date))
    
    # Reading in IFG data
    # Add in: prefix, suffix, looks
    int_fid = int_prefix+'_' + str(m_date) + '-' + str(s_date) +'_'+int_suffix+'_'+int_looks+'.unw'
    int_did = 'int_' + str(m_date)+'_' + str(s_date)
    int_data = np.fromfile(int_dir + '/' + int_did + '/' + int_fid, np.float32).reshape(rdr_ny,rdr_nx*2)
    
    # split into bands...
    i_data = int_data[:,rdr_nx:]
    c_data = int_data[:,:rdr_nx]
    
    # Read in ERA data for master and slave respectively
    infile = era_dir+'/'+str(m_date)+'_mdel_'+int_looks+'.unw'
    ds = gdal.Open(infile, gdal.GA_ReadOnly)
    ds_band2 = ds.GetRasterBand(2)
    e_mdata = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
    del ds

    infile = era_dir+'/'+str(s_date)+'_mdel_'+int_looks+'.unw'
    ds = gdal.Open(infile, gdal.GA_ReadOnly)
    ds_band2 = ds.GetRasterBand(2)
    e_sdata = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
    del ds

    # Time difference ERA data (already in slant and radians)
    e_data = np.subtract(e_sdata, e_mdata)
    
    # Now minimise IFG = ERA + RAMP, use the coh/colinearity for weighting inversion
    ramp_pars = ramp_fit_int_era(c_data, i_data, e_data)
    
    return ramp_pars

def correct_int_era_recons(m_date, s_date, ramp_pars, ramp_pars_recons, maps='yes', phase_era='yes'):
    '''
    Compare ramp and reconstructed ramp, perform correction and save.
    '''
    print('Correcting IFG with IFG = ERA + RAMP: {}_{}'.format(m_date, s_date))

    # Read in ERA data for master and slave respectively
    infile = era_dir+'/'+str(m_date)+'_mdel_'+int_looks+'.unw'
    ds = gdal.Open(infile, gdal.GA_ReadOnly)
    ds_band2 = ds.GetRasterBand(2)
    e_mdata = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
    del ds
    
    infile = era_dir+'/'+str(s_date)+'_mdel_'+int_looks+'.unw'
    ds = gdal.Open(infile, gdal.GA_ReadOnly)
    ds_band2 = ds.GetRasterBand(2)
    e_sdata = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
    del ds

    # Reading in IFG data
    # Add in: prefix, suffix, looks
    int_fid = int_prefix+'_' + str(m_date) + '-' + str(s_date) +'_'+int_suffix+'_'+int_looks+'.unw'
    int_did = 'int_' + str(m_date)+'_' + str(s_date)

    ds = gdal.Open(int_dir + '/' + int_did + '/' + int_fid, gdal.GA_ReadOnly)
    ds_band1 = ds.GetRasterBand(1)
    ds_band2 = ds.GetRasterBand(2)
    c_data= ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
    i_data = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
    
    # Time difference ERA data (already in slant and radians)
    e_data = np.subtract(e_sdata, e_mdata)
    
    # rebuild total G matrix (for mapping)
    G=np.zeros((len(i_data.flatten()),3))
    for i in range(rdr_ny):
        G[i*rdr_nx:(i+1)*rdr_nx,0] = i  
        G[i*rdr_nx:(i+1)*rdr_nx,1] = np.arange(rdr_nx) 
    G[:,2] = 1
    
    # compute residual ramp for ramp and recons ramp
    ramp_res = np.dot(G, ramp_pars).reshape(rdr_ny, rdr_nx)
    ramp_res_recons = np.dot(G, ramp_pars_recons).reshape(rdr_ny, rdr_nx)
    
    # correct IFG
    ie_corr = i_data - e_data - ramp_res
    ie_corr_recons = i_data - e_data - ramp_res_recons
    
    # Slice out reference zone
    ie_data_refz = ie_corr[ref_zone[0]:ref_zone[1],ref_zone[2]:ref_zone[3]]
    ie_data_refz_recons = ie_corr_recons[ref_zone[0]:ref_zone[1],ref_zone[2]:ref_zone[3]]
    
    # Do the same for 
    c_data_refz = c_data[ref_zone[0]:ref_zone[1],ref_zone[2]:ref_zone[3]]
    indexref = np.nonzero(np.logical_and(~np.isnan(ie_data_refz), c_data_refz>1.e-6))
    
    # Use the colinearity to weight mean in reference zone constant
    ie_data_refz_cst = np.nansum(ie_data_refz[indexref]*c_data_refz[indexref])/ np.nansum(c_data_refz[indexref])
    ie_data_refz_recons_cst = np.nansum(ie_data_refz_recons[indexref]*c_data_refz[indexref])/ np.nansum(c_data_refz[indexref])
    
    # Apply reference zone constant shift to both IFG and ERA
    i_data = i_data -  ie_data_refz_cst
    e_data = e_data - ie_data_refz_cst
    ie_corr = ie_corr -  ie_data_refz_cst
    
    i_data_recons = i_data -  ie_data_refz_recons_cst
    e_data_recons = e_data - ie_data_refz_recons_cst
    ie_corr_recons = ie_corr_recons -  ie_data_refz_recons_cst
    
    # Save output in the output era int directory
    if not os.path.exists(int_edir+'/'+int_did):
        os.mkdir(int_edir+'/'+int_did)
    
    if maps=='yes':
        # Plot phase maps of ifg, era, ramp, ifg_corr
        maps_phase_era_ramp_corr_recons(int_did, m_date, s_date, i_data, e_data + ramp_res, ie_corr, ramp_res-ie_data_refz_cst, ramp_res_recons-ie_data_refz_recons_cst, ie_corr_recons, c_data)
    
    if phase_era=='yes':
        # Plot IFG phase vs ERA+RAMP (return correlation and gradient)
        ie_gradient, ie_correlation = phase_vs_era(int_did, m_date, s_date, i_data_recons-ramp_res_recons, e_data_recons, c_data, dem_data)
    
    # Set areas outside = 0 using c_data (amplitude)
    ie_corr_recons[c_data==0] = 0
    
    ie_corr_fid = int_edir+'/'+int_did+'/'+'eracorr_recons_'+int_fid

    dst_ds = driver.Create(ie_corr_fid, ds.RasterXSize, ds.RasterYSize, 2, gdal.GDT_Float32)
    dst_band1 = dst_ds.GetRasterBand(1)
    dst_band2 = dst_ds.GetRasterBand(2)
    dst_band1.WriteArray(c_data,0,0)
    dst_band2.WriteArray(ie_corr_recons,0,0)
    dst_band1.FlushCache()
    dst_band2.FlushCache()
    del dst_ds, ds
    
    # Copy across par files too...
    if not os.path.exists(ie_corr_fid+'.rsc'):
        copyfile(int_dir+'/'+int_did+'/'+int_fid+'.rsc', ie_corr_fid+'.rsc')
    
    # Caculate standard deviation for inversion weighting...
    e_std = np.nanstd(e_data + ramp_res_recons)
    
    return m_date, s_date, ie_gradient, ie_correlation, e_std

def maps_phase_era_ramp_corr_recons(int_did, m_date, s_date, map1, map2, map3, map4, map5, map6, map_mask):
    print('Plotting IFG - ERA - RAMP = IFGcorr phase maps: {}_{}'.format(m_date, s_date))
    
    map1 = deepcopy(map1)
    map2 = deepcopy(map2)
    map3 = deepcopy(map3)
    map4 = deepcopy(map4)
    map5 = deepcopy(map5)
    map6 = deepcopy(map6)
    
    ### Set to NaN outside areas for plotting...
    map1[map_mask==0] = np.nan
    map2[map_mask==0] = np.nan
    map3[map_mask==0] = np.nan
    map4[map_mask==0] = np.nan
    map5[map_mask==0] = np.nan
    map6[map_mask==0] = np.nan
    
    # Initiate figure
    fig = plt.figure(2,figsize=(18,18))
    
    # Set all figures on same colorbar
    vmax = np.nanmax([np.nanpercentile(map1, 97.5),
                        np.nanpercentile(map2, 97.5), 
                        np.nanpercentile(map3, 97.5),
                        np.nanpercentile(map4, 97.5),
                        np.nanpercentile(map5, 97.5),
                        np.nanpercentile(map6, 97.5),
                        ])
    vmin = np.nanmin([np.nanpercentile(map1, 2.5),
                        np.nanpercentile(map2, 2.5), 
                        np.nanpercentile(map3, 2.5),
                        np.nanpercentile(map4, 2.5),
                        np.nanpercentile(map5, 2.5),
                        np.nanpercentile(map6, 2.5),
                        ])
    
    ## Symmetrical colorbar
    #if vmax < 0:
        #vmax = -vmin
    #elif vmin > 0:
        #vmin = -vmax
    #elif abs(vmax)>abs(vmin):
        #vmin = -vmax
    #else:
        #vmax = -vmin
    
    # Plotting the raw unw IFG data
    ax1 = fig.add_subplot(2,3,1)
    cax1 = ax1.imshow(map1, cmap=cmap, vmax=vmax, vmin=vmin, interpolation='None', extent=None)
    # Add the patch to the Axes
    rect = patches.Rectangle((ref_zone[2],ref_zone[0]),(ref_zone[3]-ref_zone[2]),(ref_zone[1]-ref_zone[0]), linewidth=1, edgecolor='grey', facecolor='none')
    ax1.add_patch(rect)
    setp(ax1.get_xticklabels(), visible=True)
    setp(ax1.get_yticklabels(), visible=True)
    divider = make_axes_locatable(ax1)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax1, cax=c)
    ax1.set_title('$\phi_{IFG}$: '+str(m_date)+'_'+str(s_date))
    
    # Plotting the raw unw IFG data
    ax2 = fig.add_subplot(2,3,2)
    cax2 = ax2.imshow(map2, cmap=cmap, vmax=vmax, vmin=vmin, interpolation='None', extent=None)
    # Add the patch to the Axes
    rect = patches.Rectangle((ref_zone[2],ref_zone[0]),(ref_zone[3]-ref_zone[2]),(ref_zone[1]-ref_zone[0]), linewidth=1, edgecolor='grey', facecolor='none')
    ax2.add_patch(rect)
    setp(ax2.get_xticklabels(), visible=True)
    setp(ax2.get_yticklabels(), visible=False)
    divider = make_axes_locatable(ax2)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax2, cax=c)
    ax2.set_title('$\phi_{ERA} + \phi_{RAMP}$: '+str(m_date)+'_'+str(s_date))
    
    # Plotting the raw unw IFG data
    ax3 = fig.add_subplot(2,3,3)
    cax3 = ax3.imshow(map3, cmap=cmap, vmax=vmax, vmin=vmin, interpolation='None', extent=None)
    # Add the patch to the Axes
    rect = patches.Rectangle((ref_zone[2],ref_zone[0]),(ref_zone[3]-ref_zone[2]),(ref_zone[1]-ref_zone[0]), linewidth=1, edgecolor='grey', facecolor='none')
    ax3.add_patch(rect)
    setp(ax3.get_xticklabels(), visible=True)
    setp(ax3.get_yticklabels(), visible=False)
    divider = make_axes_locatable(ax3)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax3, cax=c)
    ax3.set_title('$\phi_{IFG} - \phi_{ERA} - \phi_{RAMP}$: '+str(m_date)+'_'+str(s_date))
    
    # Plotting the raw unw IFG data
    ax1 = fig.add_subplot(2,3,4)
    cax1 = ax1.imshow(map4, cmap=cmap, vmax=vmax, vmin=vmin, interpolation='None', extent=None)
    # Add the patch to the Axes
    rect = patches.Rectangle((ref_zone[2],ref_zone[0]),(ref_zone[3]-ref_zone[2]),(ref_zone[1]-ref_zone[0]), linewidth=1, edgecolor='grey', facecolor='none')
    ax1.add_patch(rect)
    setp(ax1.get_xticklabels(), visible=True)
    setp(ax1.get_yticklabels(), visible=True)
    divider = make_axes_locatable(ax1)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax1, cax=c)
    ax1.set_title('$\phi_{RAMP}$: '+str(m_date)+'_'+str(s_date))
    
    # Plotting the raw unw IFG data
    ax2 = fig.add_subplot(2,3,5)
    cax2 = ax2.imshow(map5, cmap=cmap, vmax=vmax, vmin=vmin, interpolation='None', extent=None)
    # Add the patch to the Axes
    rect = patches.Rectangle((ref_zone[2],ref_zone[0]),(ref_zone[3]-ref_zone[2]),(ref_zone[1]-ref_zone[0]), linewidth=1, edgecolor='grey', facecolor='none')
    ax2.add_patch(rect)
    setp(ax2.get_xticklabels(), visible=True)
    setp(ax2.get_yticklabels(), visible=False)
    divider = make_axes_locatable(ax2)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax2, cax=c)
    ax2.set_title('$\phi_{RAMPRECONS}$: '+str(m_date)+'_'+str(s_date))
    
    # Plotting the raw unw IFG data
    ax3 = fig.add_subplot(2,3,6)
    cax3 = ax3.imshow(map6, cmap=cmap, vmax=vmax, vmin=vmin, interpolation='None', extent=None)
    # Add the patch to the Axes
    rect = patches.Rectangle((ref_zone[2],ref_zone[0]),(ref_zone[3]-ref_zone[2]),(ref_zone[1]-ref_zone[0]), linewidth=1, edgecolor='grey', facecolor='none')
    ax3.add_patch(rect)
    setp(ax3.get_xticklabels(), visible=True)
    setp(ax3.get_yticklabels(), visible=False)
    divider = make_axes_locatable(ax3)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax3, cax=c)
    ax3.set_title('$\phi_{IFG} - \phi_{ERA} - \phi_{RAMPRECONS}$:: '+str(m_date)+'_'+str(s_date))
    
    fig.savefig(int_edir+'/'+int_did+'/'+'eracorr_recons_map_'+str(m_date)+'_'+str(s_date)+'.png', format='PNG', bbox_inches='tight')
    if plot == 'yes':
        print('Plotting correction...')
        plt.show()
    
    plt.clf()
    
def phase_vs_era(int_did, m_date, s_date, i_data, e_data, c_data, dem_data):
    print('Generating IFG phase vs ERA phase estimation: {}_{}'.format(m_date, s_date))
    
    i_data = deepcopy(i_data)
    e_data = deepcopy(e_data)
    dem_data = deepcopy(dem_data)
    
    # 
    i_data[c_data==0] = np.nan
    e_data[c_data==0] = np.nan
    dem_data[c_data==0] = np.nan
    
    #########################
    ## Plotting IFG against ERA corr
    data_interval = 10
    mask = ~np.isnan(i_data) & ~np.isnan(e_data)
    
    x_mask = i_data[mask]
    y_mask = e_data[mask]
    z_mask = dem_data[mask]
    
    x = x_mask[::data_interval].flatten()
    y = y_mask[::data_interval].flatten()
    z = z_mask[::data_interval].flatten()
    
    if x.size == 0 or y.size ==0:
        print('INPUT ARRAY ZERO:', m_date, s_date)
    
    # Calculate the linear trend, pearson coefficient
    
    # clean outliers
    index = np.nonzero(
        np.logical_and(x>np.nanpercentile(x, 1),
        np.logical_and(x<np.nanpercentile(x, 99),
        np.logical_and(y>np.nanpercentile(y, 1), y<np.nanpercentile(y, 99)
        ))))
    
    if x[index].shape[0] == 0:
        print('MASK ERROR:', m_date, s_date)
    
    # linear regression of data
    slope, intercept, rvalue, pvalue, stderr = linregress(x[index], y[index])
    
    # curvefit to sliding median with std weighting
    binned_xdata, running_median, running_std = sliding_median(x[index], y[index])
    
    def linear_f(x, a, b):
        return a*x + b
    
    # watch out for first element of sliding median equal to NaN
    popt, pcov = curve_fit(linear_f, binned_xdata[1:], running_median[1:], sigma=running_std[1:], absolute_sigma=True)
    
    # First, create the figure (size = x,y)
    fig = plt.figure(1, figsize=(10, 10))
    
    # Now, create the gridspec structure, as required
    gs = gridspec.GridSpec(3,2, height_ratios=[0.05, .8, .2], width_ratios=[.8, .2])
    # 3 rows, 4 columns, each with the required size ratios. 
    gs.update(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.02, hspace=0.02)
    
    # First, the scatter plot
    # --------------------------------------------------------
    ax1 = plt.subplot(gs[1,0]) # place it where it should be.
    # --------------------------------------------------------
    cmax = np.nanpercentile(z, 97.5)
    cmin = np.nanpercentile(z, 2.5)
    
    # The plot itself
    plt1 = ax1.scatter(x, y, c = z, 
                    marker = 's', s = 5, edgecolor = 'none', alpha = 0.05,
                    cmap = cmap, vmin = cmin , vmax = cmax, rasterized=True)
    ax1.plot(binned_xdata, running_median, color='grey', lw=2)
    ax1.fill_between(binned_xdata, running_median-running_std, running_median+running_std, color='grey',alpha=0.3)
    ax1.plot(x, linear_f(x,*popt), color='black', lw=2)
    
    # Define the limits, labels, ticks as required
    ax1.grid(True)
    ax1.set_axisbelow(True)
    xlim_buffer = (np.nanpercentile(x, 99.85) - np.nanpercentile(x, 0.15))*0.1
    ylim_buffer = (np.nanpercentile(y, 99.85) - np.nanpercentile(y, 0.15))*0.1
    ax1.set_xlim(np.nanpercentile(x, 0.15)-xlim_buffer, np.nanpercentile(x, 99.85)+xlim_buffer)
    ax1.set_ylim(np.nanpercentile(y, 0.15)-ylim_buffer, np.nanpercentile(y, 99.85)+ylim_buffer)
    ax1.xaxis.set_major_locator(plticker.MultipleLocator(base=5))
    ax1.yaxis.set_major_locator(plticker.MultipleLocator(base=5))
    ax1.set_ylabel('$\phi_{ERA}$: '+str(m_date)+'_'+str(s_date))
    
    # scatter labels to fine for legend so...
    labels = ['Gradient: {:.2f}'.format(popt[0]), 'Correlation: {:.2f}'.format(rvalue)]
    label_colours = [(0,0,0,0),(0,0,0,0)]
    recs = []
    for i in range(0,len(label_colours)):
        recs.append(patches.Rectangle((0,0),1,1,fc=label_colours[i]))
    plt.legend(recs,labels,loc=2)
    
    # and let us not forget the colorbar  above !
    # --------------------------------------------------------
    cbax = plt.subplot(gs[0,0]) # Place it where it should be.
    # --------------------------------------------------------

    cb = Colorbar(ax = cbax, mappable = plt1, orientation = 'horizontal', ticklocation = 'top')
    cb.set_label(r'Elevation (m)', labelpad=10)
    cb.solids.set(alpha=1)
    
    # And now the histogram
    # --------------------------------------------------------
    ax1v = plt.subplot(gs[1,1])
    # --------------------------------------------------------

    # Plot the data
    binwidth = 0.1
    bins = np.arange(np.nanpercentile(y, 0.15),np.nanpercentile(y, 99.85)+binwidth,binwidth)
    ax1v.hist(y, bins=bins, orientation='horizontal', color='grey')
    ax1v.set_xlabel('No. of Pixels')

    # Define the limits, labels, ticks as required
    ax1v.set_yticks(ax1.get_yticks()) # Ensures we have the same ticks as the scatter plot !
    ax1v.set_yticklabels([])
    ax1v.set_ylim(ax1.get_ylim())
    #ax1v.xaxis.set_major_locator(plticker.MultipleLocator(base=1000))
    ax1v.grid(True)
    ax1v.set_axisbelow(True)
    
    # And now another histogram
    # --------------------------------------------------------
    ax1h = plt.subplot(gs[2,0])
    # --------------------------------------------------------
    # Plot the data
    bins = np.arange(np.nanpercentile(x, 0.15),np.nanpercentile(x, 99.85)+binwidth,binwidth)
    ax1h.hist(x, bins=bins, orientation='vertical', color='grey')
    ax1h.set_xlabel('$\phi_{IFG} - \phi_{RAMP}$: '+str(m_date)+'_'+str(s_date))
    ax1h.set_ylabel(r'No. of pixels')

    # Define the limits, labels, ticks as required
    ax1h.set_xticks(ax1.get_xticks()) # Ensures we have the same ticks as the scatter plot !
    ax1h.set_xlim(ax1.get_xlim())
    #ax1h.yaxis.set_major_locator(plticker.MultipleLocator(base=1000))
    ax1h.grid(True)
    ax1h.set_axisbelow(True)
        
    ## Save figure
    fig.savefig(int_edir+'/'+int_did+'/'+'phase-era'+'_'+str(m_date)+'_'+str(s_date)+'.png', format='PNG', dpi=300, bbox_inches='tight')
    if plot == 'yes':
        print('Plotting phase vs. atmos...')
        plt.show()
        
    plt.clf()
    
    # Note slope = gradient, rvalue = correlation coefficient
    return popt[0], rvalue


#####################################################################################
# INIT LOG
#####################################################################################

# logging.basicConfig(level=logging.INFO,\
logging.basicConfig(level=logging.INFO,\
        format='line %(lineno)s -- %(levelname)s -- %(message)s')
logger = logging.getLogger('nsb_unw_ERAcorr_stats.log')


#####################################################################################
# READ INPUT PARAM
#####################################################################################

# read input parameters and arguments
arguments = docopt.docopt(__doc__)

if arguments["--int_dir"] != None:
    int_dir = arguments["--int_dir"]
    
if arguments["--int_prefix"] != None:
    int_prefix = arguments["--int_prefix"]

if arguments["--int_suffix"] != None:
    int_suffix = arguments["--int_suffix"]
    
if arguments["--rdr_rsc"] != None:
    rdr_rsc_file = arguments["--rdr_rsc"]

if arguments["--era_dir"] != None:
    era_dir = arguments["--era_dir"]

if arguments["--int_list"] != None:
    int_list = arguments["--int_list"]

if arguments["--rdr_dem"] != None:
    rdr_dem = arguments["--rdr_dem"]
    
if arguments["--ref_zone"] != None:
    ref_zone = list(map(int,arguments["--ref_zone"].replace(',',' ').split()))

if arguments["--int_looks"] == None:
    int_looks = '2rlks'
else:
    int_looks = arguments["--int_looks"]+'rlks'
    
if arguments["--nproc"] == None:
    nproc = 1
else:
    nproc = int(arguments["--nproc"])
    
if arguments["--estim"] == 'no':
    estim== 'no'
else:
    estim = 'yes'

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


#####################################################################################
# INITIALISE 
#####################################################################################

if __name__ == '__main__':
    print('If using multi-processing and crashing, try (bash): unset DISPLAY')

    # Use int_list to obtain IFG date list for ts processing
    m_dates, s_dates = np.loadtxt(int_list, comments='#', usecols=(0,1), unpack=True, dtype='i,i')
    N_ifgs = m_dates.shape[0]
    
    img_list = np.sort(np.unique(np.concatenate((m_dates, s_dates), axis=0)))
    M_imgs = img_list.shape[0]

    # read in radar geometry parameters
    rdr_nx, rdr_ny, rdr_dict = rdr_rsc(rdr_rsc_file, full=True)
    rdr_wl = np.float(rdr_dict['WAVELENGTH'])
    rdr_rad2mm = (rdr_wl/(4*const.pi))*1000

    # Read in DEM (assume nsb format)
    # dem = np.fromfile(rdr_dem, np.float32).reshape(rdr_ny,rdr_nx)
    # dem_data = dem[:,rdr_nx:]
    driver = gdal.GetDriverByName("roi_pac")
    ds = gdal.Open(rdr_dem, gdal.GA_ReadOnly)
    ds_band2 = ds.GetRasterBand(2)
    rdr_ny, rdr_nx = ds.RasterYSize, ds.RasterXSize
    dem_data = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)

    # Saving corrected IFGs in a new dir
    int_edir = 'INTERFERO_era'
    if not os.path.exists(int_edir):
        os.mkdir(int_edir)
    
    # Paralellize phase-era ramp estimation
    ifg_pairs = list(zip(m_dates, s_dates))
    
    if estim == 'yes':
        print('Estimating IFG = ERA + RAMP:')
    
        # Run this with a pool of agents (nproc)
        with Pool(processes=nproc) as pool:
            ramps_results = pool.starmap(correct_int_era_ramp, ifg_pairs)
        
        # Save ramp estimation results: m_date, s_date, ramp_pars
        ramp_results_array = np.concatenate((m_dates[:,None], s_dates[:,None], ramps_results), axis=1) 
        np.savetxt(int_edir+'/phase-era_ramp_estim.txt', ramp_results_array, fmt='%i, %i, %.7f, %.7f, %.7f', delimiter=' ')
    else:
        print("Skipping ramp estimation (estim='no")
        
    # Read in results from txt file saved above...
    ramps_results_txt = np.loadtxt(int_edir+'/phase-era_ramp_estim.txt', delimiter=', ')
        
    # Read in results (avoid errors from chaging interf_pair...)
    ramp_pars = np.zeros((N_ifgs, 3))
    for i in range(N_ifgs):
        idxs = np.logical_and(ramps_results_txt[:,0]==m_dates[i], ramps_results_txt[:,1]==s_dates[i])
        ramp_pars[i,0] = ramps_results_txt[idxs,2]
        ramp_pars[i,1] = ramps_results_txt[idxs,3]
        ramp_pars[i,2] = ramps_results_txt[idxs,4]
    
    ### Enforcing closure
    # Need to invert each ramp parameter independently
    d_a = ramp_pars[:,0]
    d_b = ramp_pars[:,1]
    d_c = ramp_pars[:,2]
    
    # Create a desgin matrix for inverting ramps
    G_cls = np.zeros((N_ifgs, M_imgs))
    for i in range(N_ifgs):
        m_idx = np.where(img_list==m_dates[i])[0][0]
        s_idx = np.where(img_list==s_dates[i])[0][0]
        G_cls[i,m_idx] = -1
        G_cls[i,s_idx] = 1
    
    # Calculate lsq inversion for ramp parameter a
    m_a = lst.lstsq(G_cls, d_a)[0]
    m_b = lst.lstsq(G_cls, d_b)[0]
    m_c = lst.lstsq(G_cls, d_c)[0]
    
    # Now reconstruct d_a' with inverted parameters
    d_a_recons = np.dot(G_cls, m_a)
    d_b_recons = np.dot(G_cls, m_b)
    d_c_recons = np.dot(G_cls, m_c)
    
    # Create array for storing reconstruct params m_date, s_date, ax, + by, + c
    ramp_pars_recons = np.zeros((N_ifgs, 3))
    ramp_pars_recons[:,0] = d_a_recons
    ramp_pars_recons[:,1] = d_b_recons
    ramp_pars_recons[:,2] = d_c_recons

    # Paralellize phase-era correction and estimation
    ifg_pairs_ramp_pars = list(zip(m_dates, s_dates, ramp_pars, ramp_pars_recons))
    
    # Save results as text file: m_date, s_date, ramp, ramp_recons
    ramp_recons_results_array = np.concatenate((m_dates[:,None], s_dates[:,None], ramp_pars_recons), axis=1)
    np.savetxt(int_edir+'/phase-era_ramp_recons.txt', ramp_recons_results_array, fmt='%i, %i, %.7f, %.7f, %.7f', delimiter=' ')
    
    print('Correcting reconstr unw IFGs with ERA:')
    
    # Run this with a pool of agents (nproc)
    with Pool(processes=nproc) as pool:
        results = pool.starmap(correct_int_era_recons, ifg_pairs_ramp_pars)

    print('Saving results for histogram of IFG vs. ERA (correlation and gradient)')
    
    # Save results as text file
    results_array = np.asarray(results)
    np.savetxt(int_edir+'/phase-era_stats_estim.txt', results_array, fmt='%i, %i, %.3f, %.3f, %.3f', delimiter=' ')
    
    print('Calculating inversion weighting file using correlation and atmos standard deviation.')
    
    # Calculate inversion weightings 
    sigma_corr = [i[3] for i in results]
    sigma_std = [i[4] for i in results]
    sigma_std_max = np.max(sigma_std)
    sigma_weight = np.zeros(N_ifgs)
    
    for i in range (N_ifgs):
        #sigma_weight[i] = (sigma_corr[i] + np.exp(-sigma_std[i]))/2
        sigma_weight[i] = (sigma_corr[i] + (1 - (sigma_std[i]/sigma_std_max)) )/2
        
    # Save to text file: m_date, s_date, corr, std, weighting
    sigma_corr_std_weights = np.concatenate((m_dates[:,None], s_dates[:,None], np.asarray(sigma_corr)[:,None], np.asarray(sigma_std)[:,None], np.asarray(sigma_weight)[:,None]), axis=1)
    np.savetxt(int_edir+'/phase-era_sigma_corr_weights.txt', sigma_corr_std_weights, fmt='%i, %i, %.3f, %.3f, %.3f', delimiter=' ')
