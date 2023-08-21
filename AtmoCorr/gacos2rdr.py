#!/usr/bin/env python3

"""\
gacos2rdr.py
-------------
Reads in GACOS corrections and prepares them for PyGDALSAR processing,
taken from ProcessStack.py (and ProcessStackGACOS.py, from JD) from GIAnT.

.. Author:
    JD Dianala (Oxford University). Contact for information.
    Nicholas Dodds
        Created: April, 2019.
.. Last modified:
    2019 April 28 (JD) - addded check for coverage, terminates script if GACOS area is smaller than SAR area; added docopt parser.

Usage:
    gacos2rdr.py --master_mlipar=<value> --baseline_file=<path> --gacos_dir=<path> --geo_dir=<path>\
     [--looks=<value>] [--giant_format=<yes/no>] [--nproc=<nb_cores>]

--master_mlipar= value      master par file, needed to get dimensions of interferogram
--baseline_file= path       baseline.rsc file, columns [aq_date, b_perp, b_temp, zeros, aq_date_decimal]
--gacos_dir= path           dir path containing gacos corrections (includes .ztd and .ztd.rsc files for all aqs)
--geo_dir= path             geo folder from licsar/gamma output path, needed for EQA.dem_par and look up tables
--looks= value              number of looks used in timeseries processing
--giant_format= yes/no      use yes if wanting to output corrections in GIAnT (little endian) format.
--nproc=<values>            number of processor (default: 1)

"""

print()
print()
print('Authors: Nicholas Dodds, Simon DAOUT')
print('Please cite:')
print('Dodds, N., Daout, S., Walker, R. T., Begenjev, G., Bezmenov, Y., Mirzin, R., & Parsons, B. (2022). Interseismic deformation and strain-partitioning along the Main Köpetdag Fault, Turkmenistan, with Sentinel-1 InSAR time-series. Geophysical Journal International, 230(3), 1612-1629.')
print()
print()


import os, sys, logging, glob
import numpy as np
from scipy import ndimage, misc
import matplotlib.pyplot as plt
import docopt
from contextlib import contextmanager
from functools import wraps, partial
import multiprocessing, subprocess
import glob
from datetime import datetime

##################
###  Main function
##################

def gacos2rdr(k):
    date = str(gacos_dates[k])
    gacos_ztd = gacos_dir+'/'+date+'.ztd'
    
    # read in gacos correction parameters from ztd.rsc file header (using geo_rsc function copied from GIAnT). 
    gaclons, gaclats, gacnx, gacny, gdict = geo_rsc(gacos_ztd, full=True, verbose=False)
    gacdx = np.float(gdict['X_STEP'])
    gacdy = np.float(gdict['Y_STEP'])
    
    # Make sure that GACOS coverage >= SAR geo coverage
    if ( geolats[0] < gaclats[0] or geolats[1] > gaclats[1] or \
         geolons[0] < gaclons[0] or geolons[1] > gaclons[1] ):
        
        logger.info("SAR bounds: %f - %f Long., %f - %f Lat" % (geolons[0], geolons[1], geolats[0],  geolats[1]))
        logger.info("GACOS bounds: %f - %f Long, %f - %f Lat" % (gaclons[0], gaclons[1], gaclats[0],  gaclats[1]))
        raise Exception('GACOS correction covers a smaller area than your SAR. This script is unlikely to work (rerequest GACOS for larger polygon)!')

    # read in gacos correction ztd and reshape
    gztd = np.fromfile(gacos_ztd, dtype=np.float32).reshape(gacny, gacnx)
    
    # if radar geom gacos file does not exist, create it.
    if os.path.isfile(gacos_dir+'/'+date+'_crop.rdr') == False :
        
        # resample factor so the GACOS data matches SAR data resolution
        yzoom = gacdy / geody
        xzoom = gacdx / geodx
        
        # resampling
        gztdres = ndimage.zoom(gztd, [yzoom,xzoom], order=1) ## Resampling through bilinear interpolation
        
        [gny_n,gnx_n] = gztdres.shape

        # crop GACOS ztd to SAR area, then reshape GACOS ztd
        cutmaxlat_i = int(np.round( ((geolats[1] - gaclats[1])) / geody) - 1)
        cutminlat_i = cutmaxlat_i + geolen
        cutminlon_i = int(np.round( ((geolons[0] - gaclons[0])) / geodx) - 1)
        cutmaxlon_i = cutminlon_i + geowid
        gztdres2 = gztdres[cutmaxlat_i:cutminlat_i, cutminlon_i:cutmaxlon_i]
        
        # Save into a binary file to transform to radar geometry with GAMMA
        gztdres2.astype('float32').byteswap().tofile(gacos_dir+'/'+date+'_crop.ztd.unw')
        
        # Geocode (result is in big endian)
        run('geocode %s %s %i %s %i %i 0 0 ' % (ltf, gacos_dir+'/'+date+'_crop.ztd.unw', int(geowid), gacos_dir+'/'+date+'_crop.ztd.rdr.unw', rmli_wid, rmli_len) )
        
        if giant_fmt == 'yes':
            ## Swap bytes back to little endian
            run('swap_bytes %s %s 4' % (gacos_dir+'/'+date+'_crop.ztd.rdr.unw', gacos_dir+'/'+date+'_crop.ztd.rdr'))
        
        # multi-looking data production
        if looks != 1:
            # Multi-look
            run('multi_look_MLI %s %s %s %s %i %i - - - ' % ( gacos_dir+'/'+date+'_crop.ztd.rdr.unw', rmli_par, gacos_dir+'/'+date+'_crop.ztd.rdr.unw.ml'+str(looks), gacos_dir+'/'+date+'_crop.ztd.rdr.unw.ml'+str(looks)+'.par', looks, looks))
            
            if giant_fmt == 'yes':
                ## Swap bytes back to little endian
                run('swap_bytes %s %s 4' % (gacos_dir+'/'+date+'_crop.ztd.rdr.unw.ml'+str(looks), gacos_dir+'/'+date+'_crop.ztd.rdr.ml'+str(looks)))

####################
###  Extra Functions 
####################
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
        logger.info('Starting time process: {0}'.format(self.start))
    def __exit__(self, type, value, traceback):
        logger.info('Time process: {0}s'.format((datetime.now() - self.start).total_seconds()))

def checkinfile(file):
    if os.path.exists(file) is False:
        logger.critical("File: {0} not found, Exit !".format(file))
        logger.info("File: {0} not found in {1}, Exit !".format(file,os.getcwd()))

# create generator for pool
@contextmanager
def poolcontext(*arg, **kargs):
    pool = multiprocessing.Pool(*arg, **kargs)
    yield pool
    pool.terminate()
    pool.join()

def run(cmd):
    """
    Runs a shell command, and print it before running.

    Arguments:
        cmd: string to be passed to a shell

    Both stdout and stderr of the shell in which the command is run are those
    of the parent process.
    """
    
    logger.info(cmd)
    subprocess.call(cmd, shell=True, stdout=sys.stdout, stderr=subprocess.STDOUT,
        env=os.environ)
    
    return

def geo_rsc(inname,full=False,verbose=False):
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
        infile = open(inname+'.rsc','r')
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

def pygrep(infile, rstring):
    for line in open(infile):
        if rstring in line:
            return line.split("\n")[0]


#####################################################################################
# READ INPUT PARAM
#####################################################################################

# read input parameters and arguments
arguments = docopt.docopt(__doc__)

if arguments["--master_mlipar"] == None:
    raise Exception('Input Missing: No --master_mlipar file given, look for .mli.par file of master date')
else:
    rmli_par = arguments["--master_mlipar"]

if arguments["--baseline_file"] == None:
    raise Exception('Input Missing: No --baseline_file flag given, suggest looking for baseline.rsc.')
else:
    baseline = arguments["--baseline_file"]

if arguments["--gacos_dir"] == None:
    raise Exception('Input Missing: No --gacos_dir flag given, input dir containing gacos corrections.')
else:
    gacos_dir= arguments["--gacos_dir"]

if arguments["--geo_dir"] == None:
    raise Exception('Input Missing: No --geo_dir flag given, geo folder from licsar/gamma output.')
else:
    geo_dir= arguments["--geo_dir"]
    dempar=geo_dir+'/EQA.dem_par' # geocoded dem par
    ltf = glob.glob(geo_dir+'/*.lt_fine')[0]
    
if arguments["--looks"] == None:
    looks = 1
else:
    looks = int(arguments["--looks"])
    
if arguments["--giant_format"] != 'yes':
    giant_fmt = 'no'
else:
    giant_fmt = 'yes'

if arguments["--nproc"] == None:
    nproc = 1
else:
    nproc = int(arguments["--nproc"])

#####################################################################################
# INITIALISE 
#####################################################################################

# logging.basicConfig(level=logging.INFO,\
logging.basicConfig(level=logging.INFO,\
        format='%(asctime)s -- %(levelname)s -- %(message)s')
logger = logging.getLogger('gacos2rdr.log')

logger.info("Must load module gamma prior to use! Cancel and restart if gamma has not been loaded.")


logger.debug('Read rslc mli input (radar geom, not multilooked)')
checkinfile(rmli_par)
rmli_wid = int(pygrep(rmli_par, "range_samples").split(" ")[-1])
rmli_len = int(pygrep(rmli_par, "azimuth_lines").split(" ")[-1])

logger.debug('eqa dem par input (geocoded, not multilooked)')
checkinfile(dempar)
geolon = float(pygrep(dempar, "corner_lon").split(" ")[-4])
geolat = float(pygrep(dempar, "corner_lat").split(" ")[-4])
geodx = float(pygrep(dempar, "post_lon").split(" ")[-4])
geody = float(pygrep(dempar, "post_lat").split(" ")[-4])
geowid = int(pygrep(dempar, "width").split(" ")[-1])
geolen = int(pygrep(dempar, "nlines").split(" ")[-1])
logger.debug('Get coverage (min and max lon and lat) of SAR DEM')
geolon2 = geolon + (geowid * geodx)
geolat2 = geolat + (geolen * geody)
geolons = [min(geolon,geolon2),max(geolon,geolon2)]
geolats = [min(geolat,geolat2),max(geolat,geolat2)]

logger.debug('Exract list of gacos dates')
gacos_dates = np.array(glob.glob(gacos_dir+'/*.ztd'))
for i in range(0, len(gacos_dates)):
    g_date = gacos_dates[i]
    gacos_dates[i] = os.path.splitext(os.path.basename(g_date))[0]

logger.info("Extracting date list from baseline.rsc")
base_dates = np.loadtxt(baseline, unpack=True, dtype='i')[0]

# compare if a gacos correction exists for all dates in baseline
for base_date in base_dates:
    if (np.isin(base_date, gacos_dates)) == False:
        logger.critical('{0} in {1} does not have a corresponding GACOS correction in {2}'.format(base_date,baseline,gacos_dir))
        sys.exit()

Nim = len(gacos_dates)
logger.info('Number of images: {}'.format(Nim))     

#####################################################################################
# MAIN
#####################################################################################

logger.info("Generate radar geometry gacos correction for every date")
with TimeIt():
    work = range(Nim)
    with poolcontext(processes=nproc) as pool:
        results = pool.map(gacos2rdr, work)
