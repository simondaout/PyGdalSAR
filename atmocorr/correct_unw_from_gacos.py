#!/usr/bin/env python3

"""
ND_correct_unw_from_gacos.py
-------------------
Reads in GACOS corrections in radar geometry (gacos2rdr.py),
and converts ZTD to LOS using LiCSAR inc file. Then apply the
correction to unwrapped IFGs output from LiCSAR.

.. Author:
    Nicholas Dodds, University of Oxford; July, 2019.
.. Last Modified:
    5th August 2019 (created)

Usage:
    ND_correct_unw_from_gacos.py --unw_dir=<path> --gacos_dir=<path> --int_list=<path> --mli_par=<path> --inc_file=<path> --ref_zone=<xstart,xend,ystart,yend> --dem_file=<path> [--plot=<yes/no>] [--nproc=<value>]

Options:
    -h --help                Show this screen
    --unw_dir= PATH          Folder containing LiCSAR unw IFGs in radar.
    --gacos_dir= PATH        Folder containing gacos corrections in radar.
    --int_list= PATH         List of IFGs for TS processing/correction.
    --mli_par= PATH          MLI file path to extract radar dimensions.
    --inc_file= PATH         Incidence angle file for S1 TOPS in radar.
    --ref_zone= VALUES       Area where phase is set to zero for double difference subtraction (X0,X1,Y0,Y1).
    --dem_file= PATH         DEM for comparison to GACOS.
    --plot= YES/NO           Plot figures as script runs, useful for testing (default = no).
    --nproc= VALUE           Number of cores to use for parallelizing (default =1)

"""
import os
import docopt
import glob
import numpy as np
import scipy.constants as const
from multiprocessing import Process, Pool

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as patches

###
# Functions

def pygrep(infile, rstring):
    for line in open(infile):
        if rstring in line:
            return line.split("\n")[0]

### Load colormaps
cm_locs = '/home/comethome/jdd/ScientificColourMaps5/by_platform/python/'
cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt")).reversed()

###
# Read input arguments

arguments = docopt.docopt(__doc__)

if arguments["--unw_dir"] != None:
    unw_dir = arguments["--unw_dir"]
else:
    raise Exception('Input Missing: No input unwrapped IFG directory.')

if arguments["--gacos_dir"] != None:
    gacos_dir = arguments["--gacos_dir"]
else:
    raise Exception('Input Missing: No input GACOS directory.')

if arguments["--int_list"] != None:
    int_list = arguments["--int_list"]
else:
    raise Exception('Input Missing: No input int_list file given.')

if arguments["--mli_par"] != None:
    mli_par = arguments["--mli_par"]
else:
    raise Exception('Input Missing: No input MLI par file given.')

if arguments["--inc_file"] != None:
    inc_file = arguments["--inc_file"]
else:
    raise Exception('Input Missing: No input incidence file given.')

if arguments["--ref_zone"] != None:
    ref_zone = list(map(int,arguments["--ref_zone"].replace(',',' ').split()))
else:
    raise Exception('Input Missing: Need a reference zone to be given.')

if arguments["--dem_file"] != None:
    dem_file = arguments["--dem_file"]
else:
    raise Exception('Input Missing: No input incidence file given.')

if arguments["--plot"] == 'yes':
    figs = arguments["--plot"]
else:
    figs = 'no'

if arguments["--nproc"] != None:
    n_proc = int(arguments["--nproc"])
else:
    n_proc = 1

###
# Start main function

# Use int_list to obtain IFG date list for ts processing
m_dates, s_dates = np.loadtxt(int_list, comments='#', usecols=(0,1), unpack=True, dtype='i,i')

# Extract radar dimensions from mli file
if os.path.isfile(mli_par) == True:
    mli_wid = int(pygrep(mli_par, "range_samples").split(" ")[-1])
    mli_len = int(pygrep(mli_par, "azimuth_lines").split(" ")[-1])
    rad_freq = float(pygrep(mli_par, "radar_frequency").split(" ")[-3])
else:
    raise Exception('MLI par (input) file does not exist.')

# Read in incidence file for S1 TOPS from gamma (in radians)
inc = np.fromfile(inc_file, dtype='>f4').reshape(mli_len,mli_wid)
cos_los = np.cos(inc)
inc_d = np.degrees(inc)

# rad2mm conversion factor
rad_lambda = const.c/rad_freq
rad2mm = (rad_lambda/(4*const.pi))*1000

# Saving in new dir
unw_gacos = 'unw_gacos'
if not os.path.exists(unw_gacos):
    os.mkdir(unw_gacos)

fig = plt.figure(1,figsize=(14,8))

# Plotting the incidence file
ax = fig.add_subplot(1,2,1)
vmin = np.nanpercentile(inc_d, 0.3)
vmax = np.nanpercentile(inc_d, 99.7)
cax = ax.imshow(inc_d, cmap=cmap, vmax=vmax, vmin=vmin, interpolation='bicubic', extent=None)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
ax.set_title('TOPS incidence (degrees)')

# Plotting the dem file for comparison
dem = np.fromfile(dem_file, dtype='>f4').reshape(mli_len,mli_wid)
ax = fig.add_subplot(1,2,2)
vmax = np.nanpercentile(dem, 99.7)
vmin = np.nanpercentile(dem, 0.3)
cax = ax.imshow(dem, cmap=cmap, vmax=vmax, vmin=vmin, interpolation='bicubic', extent=None)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
ax.set_title('DEM elevation (m)')

fig.savefig(unw_gacos+'/TOPS_inc_DEM.pdf', format='PDF', bbox_inches='tight')
if figs == 'yes':
    plt.show()
plt.close()

def correct_unw_gacos(m_date, s_date):
    print(m_date, s_date)
    
    # Reading in IFG data
    i_data = np.fromfile(unw_dir+'/'+str(m_date)+'_'+str(s_date)+'.unw', dtype='>f4').reshape(mli_len,mli_wid)
    
    # Read in GACOS data for master and slave respectively
    g_mdata = np.fromfile(gacos_dir+'/'+str(m_date)+'_crop.ztd.rdr.unw', dtype='>f4').reshape(mli_len,mli_wid)
    g_sdata = np.fromfile(gacos_dir+'/'+str(s_date)+'_crop.ztd.rdr.unw', dtype='>f4').reshape(mli_len,mli_wid)
    
    # Time differential of GACOS data
    g_ztd_data = np.subtract(g_mdata, g_sdata)
    # Convert GACOS ZTD to LOS
    g_los_data = np.divide(g_ztd_data, cos_los)
    # Convert GACOS LOS (in m) to radians
    g_rads = g_los_data*1000/rad2mm
    
    # Need to reference the IFG and GACOS to same pixel for subtraaction.
    i_data_c = np.mean(i_data[ref_zone[0]:ref_zone[1],ref_zone[2]:ref_zone[3]])
    i_data = i_data - i_data_c
    g_rads_c = np.mean(g_rads[ref_zone[0]:ref_zone[1],ref_zone[2]:ref_zone[3]])
    g_rads = g_rads - g_rads_c
    
    # Save just the gacos correction for reference also
    gid = gacos_dir+'/'+str(m_date)+'_'+str(s_date)+'.gcorr.unw'
    g_rads.astype('>f4').tofile(gid)
    
    # Subtract the gacos correction from the IFG
    ig_corr = np.add(i_data, g_rads)
    
    # Mask areas based on coherence file for saving...
    c_data = np.fromfile(unw_dir+'/'+str(m_date)+'_'+str(s_date)+'.filt.cc', dtype='>f4').reshape(mli_len,mli_wid)
    i_data[c_data==0] = np.nan
    g_rads[c_data==0] = np.nan
    ig_corr[c_data==0] = np.nan
    
    # Print out average phase in ref zone
    print('Mean phase in ref. zone after correction: ', "%.2f" % np.mean(ig_corr[ref_zone[0]:ref_zone[1],ref_zone[2]:ref_zone[3]]) )
    
    # Save output in a new directory
    fid = unw_gacos+'/'+str(m_date)+'_'+str(s_date)+'.unw'
    ig_corr.astype('>f4').tofile(fid)
    
    fig = plt.figure(2,figsize=(18,8))
    
    vmax = np.nanmax([np.nanpercentile(i_data, 99.7), np.nanpercentile(g_rads, 99.7), np.nanpercentile(ig_corr, 99.7)])
    vmin = np.nanmin([np.nanpercentile(i_data, 0.3), np.nanpercentile(g_rads, 0.3), np.nanpercentile(ig_corr, 0.3)])
    
    # Plotting the raw unw IFG data
    ax = fig.add_subplot(1,3,1)
    cax = ax.imshow(i_data, cmap=cmap, vmax=vmax, vmin=vmin, interpolation='bicubic', extent=None)
    # Add the patch to the Axes
    rect = patches.Rectangle((ref_zone[2],ref_zone[0]),(ref_zone[3]-ref_zone[2]),(ref_zone[1]-ref_zone[0]), linewidth=1, edgecolor='grey', facecolor='none')
    ax.add_patch(rect)
    setp(ax.get_xticklabels(), visible=True)
    setp(ax.get_yticklabels(), visible=True)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)
    ax.set_title('IFG '+str(m_date)+'_'+str(s_date)+' - $\phi_{ref}$')
    
    # Plotting the raw unw IFG data
    ax = fig.add_subplot(1,3,2)
    cax = ax.imshow(g_rads, cmap=cmap, vmax=vmax, vmin=vmin, interpolation='bicubic', extent=None)
    # Add the patch to the Axes
    rect = patches.Rectangle((ref_zone[2],ref_zone[0]),(ref_zone[3]-ref_zone[2]),(ref_zone[1]-ref_zone[0]), linewidth=1, edgecolor='grey', facecolor='none')
    ax.add_patch(rect)
    setp(ax.get_xticklabels(), visible=True)
    setp(ax.get_yticklabels(), visible=False)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)
    ax.set_title('GACOS LOS - $\phi_{ref}$')
    
    # Plotting the raw unw IFG data
    ax = fig.add_subplot(1,3,3)
    cax = ax.imshow(ig_corr, cmap=cmap, vmax=vmax, vmin=vmin, interpolation='bicubic', extent=None)
    setp(ax.get_xticklabels(), visible=True)
    setp(ax.get_yticklabels(), visible=False)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)
    ax.set_title('IFG corr. GACOS (rad)')
    
    fig.savefig(unw_gacos+'/gacos_corr_unw_'+str(m_date)+'_'+str(s_date)+'.pdf', format='PDF', bbox_inches='tight')
    if figs == 'yes':
        plt.show()
    plt.close()
    
    return

### parallelizing data analysis
# paralellize phase-gacos estimation
if __name__ == '__main__':
    # Define the dataset (list of master and slave pairs)
    ifg_pairs = list(zip(m_dates, s_dates))

    # Run this with a pool of agents (nproc)
    with Pool(processes=n_proc) as pool:
        results = pool.starmap(correct_unw_gacos, ifg_pairs)

    
'''
# Looping through IFGs
i=0
for m_date in m_dates:
    # Reading in IFG data
    print(m_date, s_dates[i])
    s_date = s_dates[i]
    i_data = np.fromfile(unw_dir+'/'+str(m_date)+'_'+str(s_date)+'.unw', dtype='>f4').reshape(mli_len,mli_wid)
    
    # Read in GACOS data for master and slave respectively
    g_mdata = np.fromfile(gacos_dir+'/'+str(m_date)+'_crop.ztd.rdr.unw', dtype='>f4').reshape(mli_len,mli_wid)
    g_sdata = np.fromfile(gacos_dir+'/'+str(s_date)+'_crop.ztd.rdr.unw', dtype='>f4').reshape(mli_len,mli_wid)
    
    # Time differential of GACOS data
    g_ztd_data = np.subtract(g_mdata, g_sdata)
    
    # Convert GACOS ZTD to LOS
    cos_los = np.cos(inc)
    
    g_los_data = np.divide(g_ztd_data, cos_los)
    # Convert GACOS LOS (in m) to radians
    g_rads = g_los_data*1000/rad2mm
    
    # Need to reference the IFG and GACOS to same pixel for subtraaction.
    i_data_c = np.mean(i_data[ref_zone[0]:ref_zone[1],ref_zone[2]:ref_zone[3]])
    i_data = i_data - i_data_c
    g_rads_c = np.mean(g_rads[ref_zone[0]:ref_zone[1],ref_zone[2]:ref_zone[3]])
    g_rads = g_rads - g_rads_c
    
    # Save just the gacos correction for reference also
    gid = gacos_dir+'/'+str(m_date)+'_'+str(s_date)+'.gcorr.unw'
    g_rads.astype('>f4').tofile(gid)
    
    # Subtract the gacos correction from the IFG
    ig_corr = np.add(i_data, g_rads)
    
    # Mask areas based on coherence file for saving...
    c_data = np.fromfile(unw_dir+'/'+str(m_date)+'_'+str(s_date)+'.filt.cc', dtype='>f4').reshape(mli_len,mli_wid)
    i_data[c_data==0] = np.nan
    g_rads[c_data==0] = np.nan
    ig_corr[c_data==0] = np.nan
    
    # Print out average phase in ref zone
    print('Mean phase in ref. zone after correction: ', "%.4f" % np.mean(ig_corr[ref_zone[0]:ref_zone[1],ref_zone[2]:ref_zone[3]]) )
    
    # Save output in a new directory
    fid = unw_gacos+'/'+str(m_date)+'_'+str(s_date)+'.unw'
    ig_corr.astype('>f4').tofile(fid)
    
    fig = plt.figure(2,figsize=(18,8))
    
    vmax = np.nanmax([np.nanpercentile(i_data, 99.7), np.nanpercentile(g_rads, 99.7), np.nanpercentile(ig_corr, 99.7)])
    vmin = np.nanmin([np.nanpercentile(i_data, 0.3), np.nanpercentile(g_rads, 0.3), np.nanpercentile(ig_corr, 0.3)])
    
    # Plotting the raw unw IFG data
    ax = fig.add_subplot(1,3,1)
    cax = ax.imshow(i_data, cmap=cmap, vmax=vmax, vmin=vmin, interpolation='bicubic', extent=None)
    # Add the patch to the Axes
    rect = patches.Rectangle((ref_zone[2],ref_zone[0]),(ref_zone[3]-ref_zone[2]),(ref_zone[1]-ref_zone[0]), linewidth=1, edgecolor='grey', facecolor='none')
    ax.add_patch(rect)
    setp(ax.get_xticklabels(), visible=True)
    setp(ax.get_yticklabels(), visible=True)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)
    ax.set_title('IFG '+str(m_date)+'_'+str(s_date)+' - $\phi_{ref}$')
    
    # Plotting the raw unw IFG data
    ax = fig.add_subplot(1,3,2)
    cax = ax.imshow(g_rads, cmap=cmap, vmax=vmax, vmin=vmin, interpolation='bicubic', extent=None)
    # Add the patch to the Axes
    rect = patches.Rectangle((ref_zone[2],ref_zone[0]),(ref_zone[3]-ref_zone[2]),(ref_zone[1]-ref_zone[0]), linewidth=1, edgecolor='grey', facecolor='none')
    ax.add_patch(rect)
    setp(ax.get_xticklabels(), visible=True)
    setp(ax.get_yticklabels(), visible=False)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)
    ax.set_title('GACOS LOS - $\phi_{ref}$')
    
    # Plotting the raw unw IFG data
    ax = fig.add_subplot(1,3,3)
    cax = ax.imshow(ig_corr, cmap=cmap, vmax=vmax, vmin=vmin, interpolation='bicubic', extent=None)
    setp(ax.get_xticklabels(), visible=True)
    setp(ax.get_yticklabels(), visible=False)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)
    ax.set_title('IFG corr. GACOS (rad)')
    
    fig.savefig(unw_gacos+'/gacos_corr_unw_'+str(m_date)+'_'+str(s_date)+'.pdf', format='PDF', bbox_inches='tight')
    if figs == 'yes':
        plt.show()
    plt.close()
    
    i+=1

# Program end.
'''
