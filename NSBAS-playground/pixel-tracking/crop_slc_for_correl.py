#!/usr/bin/env python3
# -*- coding:utf-8 -*-

#__projet__ = "masque"
#__nom_fichier__ = "fonctions"
#__author__ = "Laure MANCEAU"
#__date__ = "avril 2022"

"""\
prep_slc_for_correl.py.py
----------------------
usage:
  prep_slc_for_correl.py.py [-v] [-f] [--base_file=<path>] [--crop=<jstart,jend,ibeg,iend>] 
options:
  --base_file=<path>                Baseline file containing the dates list in the first column [Default: baseline.rsc]
  --crop=<jstart,jend,ibeg,iend>    Crop windows between colunms jstart,jend and lines  istart,iend 
  -v                                Verbose mode. Show more information about the processing
  -h --help                         Show this screen.
"""

import numpy as np
from math import *
from osgeo import gdal
import os
import logging
import getopt
from sys import argv,exit
import docopt

def checkinfile(file):
    if os.path.exists(file) is False:
        logger.critical("File: {0} not found, Exit!".format(file))
        print("File: {0} not found in {1}, Exit!".format(file,os.os.getcwd()))

class Cd(object):
    def __init__(self,dirname):
        self.dirname = dirname
    def __enter__(self):
        self.curdir = os.getcwd()
        logger.debug('Enter {0}'.format(self.dirname))
        os.chdir(self.dirname)
    def __exit__(self, type, value, traceback):
        logger.debug('Exit {0}'.format(self.dirname))
        logger.debug('Enter {0}'.format(self.curdir))
        os.chdir(self.curdir)

def mask(filename, crop):

    basename = os.path.splitext(filename)[0]; outfile = basename + 'crop.tiff'

    ds = gdal.OpenEx(filename, allowed_drivers=["ROI_PAC"])
    ds_band = ds.GetRasterBand(1) #extraction d'une bande
    proj = ds.GetProjectionRef()
    gt = ds.GetGeoTransform()
    logger.debug("Driver: {}".format(ds.GetDriver().ShortName))
    logger.debug("Size: ncols:{}, nlines:{}, n bands:{}".format(ds.RasterXSize,ds.RasterYSize,ds.RasterCount))

    # read complex
    data = ds_band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
    amp = np.absolute(data)
    phi = np.angle(data)

    # to do: crop options
    if crop == None:
        pass
    else:
        jbeg,jend,ibeg,iend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])
        amp2=amp[ibeg:iend,jbeg:jend]

    # take the log of the amplitude
    amp2 = np.log(amp2)

    # take the size of the new amp
    (y, x)=np.shape(amp2)

    # création du nouveau fichier
    drv = gdal.GetDriverByName('GTiff')
    # Carefull not to write in the input file, give different name for output ds and bands
    dst_ds = drv.Create(outfile, x, y, 1, gdal.GDT_UInt16)
    dst_band = dst_ds.GetRasterBand(1)
    dst_band.WriteArray(amp2)
    dst_ds.SetGeoTransform(gt)
    dst_ds.SetProjection(proj)
    dst_band.FlushCache()
    del dst_ds, ds
    
###########
#   MAIN  # 
###########
 gdal.GDT_UInt16
# Parse arguments: read arguments with dococpt
arguments = docopt.docopt(__doc__)

# init logger 
if arguments["-v"]:
    logging.basicConfig(level=logging.DEBUG,\
        format='%(asctime)s -- %(levelname)s -- %(message)s')
else:
    logging.basicConfig(level=logging.INFO,\
        format='%(asctime)s -- %(levelname)s -- %(message)s')
logger = logging.getLogger('filtflatunw_log.log')

if arguments["--base_file"] == None:
    base_file = 'baseline.rsc'
else:
    base_file = arguments["--base_file"]
#dates,bid=np.loadtxt(base_file,comments="#",unpack=True,dtype='i,f')
dates=np.loadtxt(base_file,usecols=0,unpack=True,dtype='int')

if arguments["--crop"]==None:
    crop = None
else:   
    crop = list(map(float,arguments["--crop"].replace(',',' ').split()))

for date in dates:
    date = str(date)
    with Cd(date):
      filename = date + '_coreg.slc'; checkinfile(filename)
      rsc = filename + '.rsc'; checkinfile(rsc)
      
      logger.info(print('Run {} ....'.format(date)))
      mask(filename, crop)

