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
    nsb_gacos2rdr.py --rdr_rsc=<path> --lt=<path> --inc_file=<path> --sim_dem_rdr=<path> --gacos_dir=<path> --int_list=<path> [--nproc=<nb_cores>] [-f]

Options:
    --rdr_rsc=path              path to radar param file
    --lt=path                   path to lookup table (rdr to geoc)
    --inc_file=PATH             path to incidence angle file for S1 TOPS in radar.
    --sim_dem_rdr=path          path to sim dem in rdr
    --gacos_dir=path            path to folder containing gacos ztds
    --int_list=path             path to interf pair list
    --nproc=values              number of processor (default: 1)
    -f                          Force mode. Overwrite output files
"""

import os, sys, logging, glob
import docopt
import numpy as np
import matplotlib.pyplot as plt 
from os import path
from scipy import ndimage
import scipy.constants as const
from osgeo import osr

import shutil, filecmp
from multiprocessing import Process, Pool
import subprocess,gdal

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as patches

### Load colormaps
cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt")).reversed()

def rm(file):
    if path.exists(file):
        os.remove(file)

def copyrsc(src,dest):
    if path.exists(dest):
        if filecmp.cmp(src,dest):
            pass
        else:
            shutil.copy(src,dest)
    else:
            shutil.copy(src,dest)

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

def convert(int_date):
    gacos_ztd_file = gacos_dir + '/' + str(int_date) + '.ztd'
    rsc = gacos_dir + '/' + str(int_date) + '.ztd.rsc'
    temp = gacos_dir + '/' + str(int_date) + '_temp.tif'
    gacos_tif_file = gacos_dir + '/' + str(int_date) + '.tif'
    if force:
        rm(gacos_tif_file);rm(temp)

    # read .rsc
    nncol = np.int(open(rsc).readlines()[0].split(None, 2)[1])
    nnlines = np.int(open(rsc).readlines()[1].split(None, 2)[1])
    xfirst = np.float(open(rsc).readlines()[6].split(None, 2)[1])
    yfirst = np.float(open(rsc).readlines()[7].split(None, 2)[1])
    xstep = np.float(open(rsc).readlines()[8].split(None, 2)[1])
    ystep = np.float(open(rsc).readlines()[9].split(None, 2)[1])
    geotransform = (xfirst, xstep, 0, yfirst, 0, ystep)
    logger.info('Read .rsc and set Geotransform: {}'.format(geotransform))
    ztd_ = np.fromfile(gacos_ztd_file, dtype='float32').reshape(nnlines,nncol)

     # transform gacos format to samething humanily readeable
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(temp, nncol, nnlines, 1, gdal.GDT_Float32)
    band = ds.GetRasterBand(1)
    band.WriteArray(ztd_)
    ds.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    band.FlushCache()
    del ds,band

    logger.info("gdalwarp -overwrite -tr "+str(geo_dx)+" "+str(geo_dy)+" -te "+str(geo_lon[0])+" "+str(geo_lat[0])+" "+str(geo_lon[1])+" "+str(geo_lat[1])+" -r bilinear "+temp+" "+gacos_tif_file+" -of GTiff")
    r = subprocess.call("gdalwarp -overwrite -tr "+str(geo_dx)+" "+str(geo_dy)+" -te "+str(geo_lon[0])+" "+str(geo_lat[0])+" "+str(geo_lon[1])+" "+str(geo_lat[1])+" -r bilinear "+temp+" "+gacos_tif_file+" -of GTiff", shell=True,stdout=sys.stdout, stderr=subprocess.STDOUT,env=os.environ)
    if r != 0:
        logger.critical(r)

    outfile = gacos_dir+'/'+str(int_date)+'.ztd.geo'
    if force:
        rm(outfile); rm(outfile+'.rsc')

    # if radar geom gacos file does not already exist...
    if os.path.isfile(outfile) != True:
        
        # read in gacos correction metadata 
        ds = gdal.Open(gacos_tif_file, gdal.GA_ReadOnly)
        ds_geo = ds.GetGeoTransform()
        gacos_dx, gacos_dy = ds_geo[1], ds_geo[5]
        gacos_nx, gacos_ny = ds.RasterXSize, ds.RasterYSize
        gacos_lon = ds_geo[0], ds_geo[0]+ ds_geo[1]*ds.RasterXSize
        gacos_lat = ds_geo[3]+ds_geo[5]*ds.RasterYSize, ds_geo[3] 

        # Make sure that GACOS coverage >= SAR geo coverage
        if ( geo_lat[0] < gacos_lat[0] or geo_lat[1] > gacos_lat[1] or \
                geo_lon[0] < gacos_lon[0] or geo_lon[1] > gacos_lon[1] ):
            logger.warning('GACOS correction covers a smaller area than your SAR (re-request GACOS for larger polygon).')

        # read in gacos correction ztd and reshape
        gacos_ztd = ds.GetRasterBand(1).ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)

        # save output...
        logger.info('Saving GACOS correction in real4 ROI_PAC format to: {}'.format(gacos_dir+'/'+str(int_date)+'.ztd.geo'))
        gacos_ztd.astype(np.float32).tofile(outfile)
        # create rsc file
        f = open(outfile+'.rsc', "w")
        f.write("""\
WIDTH                 %d
FILE_LENGTH           %d
X_FIRST               %f
Y_FIRST               %f
X_STEP                %.8f
Y_STEP                %.8f
        """ % (ds.RasterXSize, ds.RasterYSize, ds_geo[0], ds_geo[3], ds_geo[1], ds_geo[5]))
        f.close()

    elif os.path.isfile(gacos_dir+'/'+str(int_date)+'.ztd.geo') == True:
        print('GACOS geo file for {} already exists. Use force mode.'.format(int_date))
    
    # del temporary files
    rm(temp); rm(gacos_tif_file)
        
def gacosgeo2rdr(int_date):
    # in file, ztd cropped and resampled to geo
    gacos_geo_file = gacos_dir + '/' + str(int_date) + '.ztd.geo'
    
    # out file for saving ztd in rdr
    gacos_rdr_file = gacos_dir + '/' + str(int_date) + '.ztd.rdr'
   
    if force:
        rm(gacos_rdr_file)

    if os.path.exists(gacos_rdr_file) != True:
        logger.info('Generating GACOS correction in NSBAS radar geometry: {}'.format(gacos_dir+'/'+str(int_date)+'.ztd.rdr'))
        # calling roipac perl script for geocoding
        logger.info("geo2radar.pl %s %s %s %s" %(lt_file, gacos_geo_file, sim_dem_rdr, gacos_rdr_file))
        subprocess.call("geo2radar.pl %s %s %s %s" %(lt_file, gacos_geo_file, sim_dem_rdr, gacos_rdr_file), shell=True)
        copyrsc(sim_dem_rdr+'.rsc',gacos_rdr_file+'.rsc')
        if os.stat(gacos_rdr_file).st_size == 0:
            #raise Exception('Warning: Output radar geometry GACOS correction is size zero.')
            logger.critical('Output radar geometry GACOS correction is size zero.')
            print(gacos_rdr_file)
            sys.exit()
    else:
        logger.info('GACOS rdr file for {} already exists. Use force mode.'.format(int_date))

def gacosztd2los(int_date):
    logger.info('Converting GACOS in ztd, to los geometry and radians.')
    
    # In file from gacosgeo2rdr
    gacos_ztd_file = gacos_dir + '/' + str(int_date) + '.ztd.rdr'
    # Out file in los 
    gacos_los_file = gacos_dir + '/' + str(int_date) + '.rdr'
    
    if force:
        rm(gacos_los_file)

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
        copyrsc(sim_dem_rdr+'.rsc',gacos_los_file+'.rsc')
    else:
        print('GACOS rdr file for {} already exists. Use force mode.'.format(int_date))
    
    # clean ztd.rdr files
    rm(gacos_ztd_file)

def gacosrdr2png(int_date, plot='yes'):
    # out file of gacosgeo2rdr
    gacos_rdr_file = gacos_dir + '/' + str(int_date) + '.rdr'
   
    if force:
        rm(gacos_rdr_file+'.png')
    if plot=='yes' and os.path.exists(gacos_rdr_file+'.png')!=True:
        logger.info('Converting GACOS in radar geometry to PNG format: {}'.format(gacos_rdr_file+'.png'))
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
if arguments["-f"]:
    force = True
else:
    force = False

#####################################################################################
# INITIALISE 
######################################################################################

if __name__ == '__main__':
   
    logging.basicConfig(level=logging.INFO,\
        format='line %(lineno)s -- %(levelname)s -- %(message)s')
    logger = logging.getLogger('gacos2nsbrdr.log')

    logger.info('NSBAS/ROI_PAC must be loaded prior to using this script.')
    logger.info('If using multi-processing and crashing, try (bash): unset DISPLAY')
    logger.info('Some issues with multiprocessing (parallel writes) with geo2radar.pl; try setting --nproc=1')

    # read in radar geometry parameters
    rdr_nx, rdr_ny, rdr_dict = rdr_rsc(rdr_rsc_file, full=True)
    rdr_wl = np.float(rdr_dict['WAVELENGTH'])
    rdr_rad2mm = (rdr_wl/(4*const.pi))*1000
    
    # incidence file for S1 TOPS (nsbas, in degrees)
    ds = gdal.Open(inc_file, gdal.GA_ReadOnly)
    inc = ds.GetRasterBand(2).ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:rdr_ny,:rdr_nx]
    inc_rad = np.deg2rad(inc)

    # read in geocoded geometry parameters
    ds = gdal.Open(lt_file, gdal.GA_ReadOnly)
    ds_geo = ds.GetGeoTransform()
    geo_dx, geo_dy = ds.GetGeoTransform()[1], ds.GetGeoTransform()[5]
    geo_nx, geo_ny = ds.RasterXSize, ds.RasterYSize
    geo_lon = ds_geo[0], ds_geo[0]+ ds_geo[1]*ds.RasterXSize
    geo_lat = ds_geo[3]+ds_geo[5]*ds.RasterYSize, ds_geo[3] 
    
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
    
    # Convert gacos model to tiff format of the  same size than geomap.trans
    with Pool(processes=nproc) as pool:
        results = pool.starmap(convert, list(zip(int_date_list)))

    # Run gacos2geo with a pool of agents (nproc)
    #with Pool(processes=nproc) as pool:
    #    results = pool.starmap(gacos2geo, list(zip(int_date_list)))
        
    # Run gacosgeo2rdr with a pool of agents (nproc)
    with Pool(processes=nproc) as pool:
        results = pool.starmap(gacosgeo2rdr, list(zip(int_date_list)))
        
    # Run gacosztd2los with a pool of agents (nproc)
    with Pool(processes=nproc) as pool:
        results = pool.starmap(gacosztd2los, list(zip(int_date_list)))
        
    # Run gacosrdr2png with a pool of agents (nproc)
    #with Pool(processes=nproc) as pool:
    #results = pool.starmap(gacosrdr2png, list(zip(int_date_list)))
