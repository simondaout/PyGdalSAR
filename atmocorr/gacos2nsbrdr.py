#!/usr/bin/env python3

"""
nsb_gacos2rdr.py
-------------
Reads in GACOS corrections and converts to radar geometry,
for application to NSBAS generated IFGS.

.. Author:
    Nicholas Dodds (Oxford University). Contact for information.
        Created: November, 2019.
    JD Dianala (Oxford University)
        
Usage:
    nsb_gacos2rdr.py --geo_rsc=<path> --rdr_rsc=<path> --lt=<path> --inc_file=<path> --sim_dem_rdr=<path> --gacos_dir=<path> --int_list=<path> [--nproc=<nb_cores>]

Options:
    --geo_rsc=path              path to geoc param file
    --rdr_rsc=path              path to radar param file
    --lt=path                   path to lookup table (rdr to geoc)
    --inc_file=PATH             path to incidence angle file for S1 TOPS in radar.
    --sim_dem_rdr=path          path to sim dem in rdr
    --gacos_dir=path            path to folder containing gacos ztds
    --int_list=path             path to interf pair list
    --nproc=values              number of processor (default: 1)

"""

import os, sys, logging, glob
import docopt
import numpy as np
import matplotlib.pyplot as plt 
from os import path
from scipy import ndimage
import scipy.constants as const

from shutil import copyfile
from multiprocessing import Process, Pool
import subprocess

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as patches

### Load colormaps
cm_locs = '/home/comethome/jdd/ScientificColourMaps5/by_platform/python/'
cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt")).reversed()

def geo_rsc(inname, full=False, verbose=False):
        '''Reading a ROI-PAC style geocoded rsc file.

        Args:
                * inname (str): Path to the RSC file.

        Returns:
                * lon (np.array) : Array of min and max lon values.
                * lat (np.array) : Array of min and max lat values.
                * nx  (np.int)   : Number of lon bins.
                * ny  (np.int)   : Number of lat bins.

                .. note:: 
                        Currently set up to work with dem.rsc file from ROI-PAC.'''

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
        lat=np.zeros((2,1))
        lon=np.zeros((2,1))
        lat[1] = np.float(rpacdict['Y_FIRST'])
        lon[0] = np.float(rpacdict['X_FIRST'])
        if(lon[0] < 0):
                lon[0] = lon[0] + 360.0

        dx = np.float(rpacdict['X_STEP'])
        dy = np.float(rpacdict['Y_STEP'])

        lat[0] = lat[1] + dy*ny
        lon[1] = lon[0] + dx*nx

        if full:
                return lon,lat,nx,ny,rpacdict
        else:
                return lon,lat,nx,ny

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

def gacos2geo(int_date):
    # file path for gacos ztd file for date
    gacos_ztd_file = gacos_dir + '/' + str(int_date) + '.ztd'

    # if radar geom gacos file does not already exist...
    if os.path.isfile(gacos_dir+'/'+str(int_date)+'.ztd.geo') != True:
        
        # read in gacos correction parameters from ztd.rsc file header (using geo_rsc function copied from GIAnT). 
        gacos_lon, gacos_lat, gacos_nx, gacos_ny, gacos_dict = geo_rsc(gacos_ztd_file+'.rsc', full=True, verbose=False)
        gacos_dx = np.float(gacos_dict['X_STEP'])
        gacos_dy = np.float(gacos_dict['Y_STEP'])

        # Make sure that GACOS coverage >= SAR geo coverage
        if ( geo_lat[0] < gacos_lat[0] or geo_lat[1] > gacos_lat[1] or \
                geo_lon[0] < gacos_lon[0] or geo_lon[1] > gacos_lon[1] ):
            #raise Exception('GACOS correction covers a smaller area than your SAR (re-request GACOS for larger polygon).')
            print('GACOS correction covers a smaller area than your SAR (re-request GACOS for larger polygon).')

        # read in gacos correction ztd and reshape
        gacos_ztd = np.fromfile(gacos_ztd_file, dtype=np.float32).reshape(gacos_ny, gacos_nx)
        
        # resample factor so the GACOS data matches SAR data resolution
        yzoom = gacos_dy / geo_dy
        xzoom = gacos_dx / geo_dx
        
        # resampling
        import warnings
        warnings.filterwarnings('ignore', '.*output shape of zoom.*')
        gacos_ztd_res = ndimage.zoom(gacos_ztd, [yzoom,xzoom], order=1)
        
        [gacos_res_ny, gacos_res_nx] = gacos_ztd_res.shape

        # crop GACOS ztd to SAR area, then reshape GACOS ztd
        cutmaxlat_i = int(np.round( ((geo_lat[1] - gacos_lat[1])) / geo_dy) - 1)
        cutminlat_i = cutmaxlat_i + geo_ny
        cutminlon_i = int(np.round( ((geo_lon[0] - gacos_lon[0])) / geo_dx) - 1)
        cutmaxlon_i = cutminlon_i + geo_nx
        gacos_ztd_geo = gacos_ztd_res[cutmaxlat_i:cutminlat_i, cutminlon_i:cutmaxlon_i]
        
        # save output...
        print('Saving GACOS correction in NSBAS geocoded geometry to: {}'.format(gacos_dir+'/'+str(int_date)+'.ztd.geo'))
        gacos_ztd_geo.astype(np.float32).tofile(gacos_dir+'/'+str(int_date)+'.ztd.geo')
        
        # copy across the lookup table rsc file to be used as rsc file for the  gacos.ztd.geo
        if os.path.isfile(gacos_dir+'/'+str(int_date)+'.ztd.geo.rsc') != True:
            copyfile(lt_file+'.rsc', gacos_dir+'/'+str(int_date)+'.ztd.geo.rsc')
        
    elif os.path.isfile(gacos_dir+'/'+str(int_date)+'.ztd.geo') == True:
        print('GACOS geo file for {} already exists. Must delete first to re-generate.'.format(int_date))
        
def gacosgeo2rdr(int_date):
    # in file, ztd cropped and resampled to geo
    gacos_geo_file = gacos_dir + '/' + str(int_date) + '.ztd.geo'
    
    # out file for saving ztd in rdr
    gacos_rdr_file = gacos_dir + '/' + str(int_date) + '.ztd.rdr'
    
    if os.path.exists(gacos_rdr_file) != True:
        print('Generating GACOS correction in NSBAS radar geometry: {}'.format(gacos_dir+'/'+str(int_date)+'.ztd.rdr'))
        # calling roipac perl script for geocoding
        subprocess.call("geo2radar.pl %s %s %s %s" %(lt_file, gacos_geo_file, sim_dem_rdr, gacos_rdr_file), shell=True)
        if os.stat(gacos_rdr_file).st_size == 0:
            #raise Exception('Warning: Output radar geometry GACOS correction is size zero.')
            print('Warning: Output radar geometry GACOS correction is size zero.')
            print(gacos_rdr_file)
            sys.exit()
    else:
        print('GACOS rdr file for {} already exists. Must delete *.rdr *geo.hst first to re-generate.'.format(int_date))

def gacosztd2los(int_date):
    print('Converting GACOS in ztd, to los geometry and radians.')
    
    # In file from gacosgeo2rdr
    gacos_ztd_file = gacos_dir + '/' + str(int_date) + '.ztd.rdr'
    # Out file in los 
    gacos_los_file = gacos_dir + '/' + str(int_date) + '.rdr'
    
    if os.path.exists(gacos_los_file) != True:
        # Import gacos ztd file
        gacos_ztd_data = np.fromfile(gacos_ztd_file, np.float32).reshape(rdr_ny,rdr_nx)
        # Convert GACOS ZTD to LOS
        gacos_los_data = np.divide(gacos_ztd_data, np.cos(inc_rad))
        # Convert GACOS LOS (in m) to radians
        gacos_los_data_rads = gacos_los_data*1000/rdr_rad2mm
        # Revert zero inc areas to zero
        gacos_los_data_rads[inc_rad==0] = 0
        # Save output
        gacos_los_data_rads.astype(np.float32).tofile(gacos_los_file)
    else:
        print('GACOS rdr file for {} already exists. Must delete *.rdr *geo.hst first to re-generate.'.format(int_date))

def gacosrdr2png(int_date, plot='yes'):
    # out file of gacosgeo2rdr
    gacos_rdr_file = gacos_dir + '/' + str(int_date) + '.rdr'
    
    if plot=='yes' and os.path.exists(gacos_rdr_file+'.png')!=True:
        print('Converting GACOS in radar geometry to PNG format: {}'.format(gacos_rdr_file+'.png'))
        # plot png of rdr gacos corrections for checking
        gacos_rdr_data = np.fromfile(gacos_rdr_file, np.float32).reshape((rdr_ny,rdr_nx))
        
        # plotting
        fig = plt.figure(1,figsize=(6,8))
        ax = fig.add_subplot(1,1,1)
        vmin = np.nanpercentile(gacos_rdr_data, 0.5)
        vmax = np.nanpercentile(gacos_rdr_data, 99.5)
        cax = ax.imshow(gacos_rdr_data, cmap=cmap, vmax=vmax, vmin=vmin, interpolation='bicubic', extent=None)
        divider = make_axes_locatable(ax)
        c = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax, cax=c)
        ax.set_title('GACOS in radar: {}'.format(int_date))
        
        # saving
        fig.savefig(gacos_rdr_file+'.png', format='PNG', bbox_inches='tight')
        fig.clf()
        
#####################################################################################
# READ INPUT PARAM
#####################################################################################

# read input parameters and arguments
arguments = docopt.docopt(__doc__)

if arguments["--geo_rsc"] != None:
    geo_rsc_file = arguments["--geo_rsc"]

if arguments["--rdr_rsc"] != None:
    rdr_rsc_file = arguments["--rdr_rsc"]
    
if arguments["--lt"] != None:
    lt_file = arguments["--lt"]

if arguments["--inc_file"] != None:
    inc_file = arguments["--inc_file"]

if arguments["--sim_dem_rdr"] != None:
    sim_dem_rdr = arguments["--sim_dem_rdr"]
    
if arguments["--gacos_dir"] != None:
    gacos_dir = arguments["--gacos_dir"]

if arguments["--int_list"] != None:
    int_list = arguments["--int_list"]

if arguments["--nproc"] == None:
    nproc = 1
else:
    nproc = int(arguments["--nproc"])
    
#####################################################################################
# INITIALISE 
######################################################################################

if __name__ == '__main__':
    
    print('NSBAS/ROI_PAC must be loaded prior to using this script.')
    print('If using multi-processing and crashing, try (bash): unset DISPLAY')
    print('Some issues with multiprocessing (parallel writes) with geo2radar.pl; try setting --nproc=1')

    # read in radar geometry parameters
    rdr_nx, rdr_ny, rdr_dict = rdr_rsc(rdr_rsc_file, full=True)
    rdr_wl = np.float(rdr_dict['WAVELENGTH'])
    rdr_rad2mm = (rdr_wl/(4*const.pi))*1000
    
    # incidence file for S1 TOPS (nsbas, in degrees)
    inc_head = np.fromfile(inc_file, np.float32)[:rdr_ny*rdr_nx*2].reshape(rdr_ny, rdr_nx*2)
    inc = inc_head[:,:rdr_nx]
    inc_rad = np.deg2rad(inc)

    # read in geocoded geometry parameters
    geo_lon, geo_lat, geo_nx, geo_ny, geo_dict = geo_rsc(geo_rsc_file, full=True, verbose=False)
    geo_dx = np.float(geo_dict['X_STEP'])
    geo_dy = np.float(geo_dict['Y_STEP'])

    # extract interferogram date list from int pairs
    date1_list = np.loadtxt(int_list, unpack=True, dtype='i')[0]
    date2_list = np.loadtxt(int_list, unpack=True, dtype='i')[1]
    int_all_list = np.concatenate((date1_list, date2_list), axis=0)
    int_date_list = list(set(int_all_list))
    int_date_list.sort()    
    
    # check all int dates have a gacos correction
    for int_date in int_date_list:
        if os.path.exists(gacos_dir + '/' + str(int_date) + '.ztd') != True:
            raise Exception('No GACOS correction in {} for {} in {}.'.format(gacos_dir,int_date,int_list))

    # Run gacos2geo with a pool of agents (nproc)
    with Pool(processes=nproc) as pool:
        results = pool.starmap(gacos2geo, list(zip(int_date_list)))
        
    # Run gacosgeo2rdr with a pool of agents (nproc)
    with Pool(processes=nproc) as pool:
        results = pool.starmap(gacosgeo2rdr, list(zip(int_date_list)))
        
    # Run gacosztd2los with a pool of agents (nproc)
    with Pool(processes=nproc) as pool:
        results = pool.starmap(gacosztd2los, list(zip(int_date_list)))
        
    # Run gacosrdr2png with a pool of agents (nproc)
    with Pool(processes=nproc) as pool:
        results = pool.starmap(gacosrdr2png, list(zip(int_date_list)))
