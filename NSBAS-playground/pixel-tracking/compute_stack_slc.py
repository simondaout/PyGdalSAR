#!/usr/bin/env python3
# -*- coding:utf-8 -*-

############################################
# Author        : Simon DAOUT (CRPG)
############################################

"""\
compute_stack_slc.py
----------------------
usage:
  compute_stack_slc.py [-v] [-f] [--base_file=<path>] [--outfile=<path>] 
options:
  --base_file=<path>                Baseline file containing the dates list in the first column [Default: baseline.rsc]
  --outfile=<path>                  Name of the output file [Default: stack_amp_da.unw]
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
import shutil

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

def extract_log_amp(filename):
    
    logger.info(print('Extract Amplitude {} ....'.format(date)))
    ds = gdal.OpenEx(filename, allowed_drivers=["ROI_PAC"])
    ds_band = ds.GetRasterBand(1) #extraction d'une bande
    logger.debug("Driver: {}".format(ds.GetDriver().ShortName))
    logger.debug("Size: ncols:{}, nlines:{}, n bands:{}".format(ds.RasterXSize,ds.RasterYSize,ds.RasterCount))

    # read complex
    data = ds_band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
    amp = np.absolute(data)

    # take the log of the amplitude
    amp[amp>0] = np.log(amp[amp>0])

    return amp, ds.RasterXSize, ds.RasterYSize

###########
#   MAIN  # 
###########

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
dates,bid=np.loadtxt(base_file,comments="#",unpack=True,dtype='i,f')

if arguments["--outfile"] == None:
    outfile = 'stack_amp_da.unw'
else:
    outfile = arguments["--outfile"]

logger.info('Extrack maxumun length from all files')
ncols,nlines = 0, 0
for date in dates:
    date = str(date)
    with Cd(date):
      filename = date + '_coreg.slc'; checkinfile(filename)
      ds = gdal.OpenEx(filename, allowed_drivers=["ROI_PAC"])
      ds_band = ds.GetRasterBand(1) #extraction d'une bande
      if ds.RasterYSize > nlines :
        nlines = ds.RasterYSize
        save_date = date
      if ds.RasterXSize > ncols :
        ncols = ds.RasterXSize

logger.info('Compute Mean Amplitude and STD Amplitude')
stack,sigma,weight = np.zeros((nlines,ncols)),np.zeros((nlines,ncols)),np.zeros((nlines,ncols))
for date in dates:
    date = str(date)
    with Cd(date):
      filename = date + '_coreg.slc'; checkinfile(filename)
      
      amp = np.zeros((nlines,ncols))
      a,w,l = extract_log_amp(filename)
      amp[:l,:w] = a
      
      stack = stack + amp
      sigma = sigma + amp**2
      
      w = np.zeros((nlines,ncols))
      index = np.nonzero(amp)
      w[index] = 1
      weight = weight + w

# compute stack and sigma
stack[weight>0] = stack[weight>0]/weight[weight>0]
sigma[weight>0] = np.sqrt(sigma[weight>0]/weight[weight>0] - (stack[weight>0]/weight[weight>0])**2)
da = np.zeros((nlines,ncols))
da[sigma>0] = 1./sigma[sigma>0]
#print(np.nanmean(stack), np.mean(da))

logger.info('Save output file {}'.format(outfile))
ncols,nlines = np.shape(stack)
drv = gdal.GetDriverByName("roi_pac")
dst_ds = drv.Create(outfile, nlines, ncols, 2, gdal.GDT_Float32)
dst_band1 = dst_ds.GetRasterBand(1)
dst_band2 = dst_ds.GetRasterBand(2)
dst_band1.WriteArray(stack,0,0)
dst_band2.WriteArray(da,0,0)
del dst_ds
# copy rsc date with maximum length
shutil.copy(save_date + '/' + save_date + '_coreg.slc.rsc', outfile + '.rsc')



