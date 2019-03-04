#!/usr/bin/env python2
# -*- coding: utf-8 -*-
############################################

############################################
# Author        : Simon DAOUT (Oxford)
############################################


# gdal
import gdal, shutil
gdal.UseExceptions()
# system
from os import path, environ, system
# plot
import matplotlib
if environ["TERM"].startswith("screen"):
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from datetime import datetime
# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided
# scipy
import scipy.optimize as opt
import scipy.linalg as lst
from __future__ import print_function
import logging


##################################################################################
###  INITIALISE
##################################################################################

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('filtcorunw')

# figures
nfigure=0

# int and dates
date_1,date_2=np.loadtxt(int_list,comments="#",unpack=True,usecols=(0,1),dtype='i,i')
kmax=len(date_1)
print("number of interferogram: ",kmax)
im = []; bt = []
for date1,date2 in zip(date_1,date_2):
    if date1 not in im: im.append(date1)
    if date2 not in im: im.append(date2)
nmax=len(im)
print("number of image: ",nmax)
imd = date2dec(im)
cst = np.copy(imd[0])
for i in xrange((nmax)):
    bt.append(imd[i]-cst)


##################################################################################
###  DEF FUNCTIONS
##################################################################################

def look(kk, look, prefix, suffix):
	date1, date2 = date_1[kk], date_2[kk]
	corfile = str(date1) + '-' + str(date2) + '_' + rlook + 'rlks.cor'
	intfile = str(prefix) + str(date1) + '-' + str(date2) + str(suffix) + '_' +  rlook + 'rlks.int'
	
	# lets assume that int are always created at 2looks...
	rlook = int(look - 2)

	logger.debug('Look file {0} in {1} look'.format(intfile,rlook))
	try:
		os.system("look.pl "+str(intfile)+" "+str(rlook))
		os.system("look.pl "+str(corfile)+" "+str(rlook))
	except:
		logger.warning('Cant look file {0} in {1} look'.format(intfile,rlook))

def filter(file, ftype, outfile, corfile, window):
	logger.debug('Filter filter {0} with {1} filter type'.format(file,ftype))
	try:
		os.system("nsb_SWfilter.pl "+str(file)+" "+str(outfile)+" "+str(corfile)+" "+str(window)+" "+str(0.05)+" "+str(ftype))
	except:
		logger.warning('Cant filter {0} with {1} filter type'.format(file,ftype))

def replaceAmp(kk, rlook, prefix, suffix):
	date1, date2 = date_1[kk], date_2[kk]
	logger.debug('Replace Amplitude by coherence on int. {0}-{1}'.format(date1,date2))

	corfile = str(date1) + '-' + str(date2) + '_' + rlook + 'rlks.cor'
	ds_cor = gdal.Open(corfile, gdal.GA_ReadOnly)
	intfile = str(prefix) + str(date1) + '-' + str(date2) + str(suffix) + '_' +  rlook + 'rlks.int'
	rscfile = str(prefix) + str(date1) + '-' + str(date2) + str(suffix) + '_' +  rlook + 'rlks.int.rsc'
	ds_int = gdal.Open(infile, gdal.GA_ReadOnly)
	driver = ds.GetDriver()

	if os.path.exists(corfile):
		outfile = 'coh_' + str(prefix) + str(date1) + '-' + str(date2) + str(suffix) + '_' +  rlook + 'rlks.int'
		
		cor_band = ds_cor.GetRasterBand(2)
		cor = cor_band.ReadAsArray(0, 0,
                   ds_cor.RasterXSize, ds_cor.RasterYSize,
                   ds_cor.RasterXSize, ds_cor.RasterYSize)

		phase_band = ds_int.GetRasterBand(1)
		phi = np.angle(phase_band.ReadAsArray(0, 0,
               ds.RasterXSize, ds.RasterYSize,
               ds.RasterXSize, ds.RasterYSize))

		newphi = np.complex(cor,phi)
		dst_ds = driver.Create(outfile, ds_int.RasterXSize, ds_int.RasterYSize, ds_int.RasterCount, phi_band.DataType)
		dst_band1 = dst_ds.GetRasterBand(1)
		dst_band1.WriteArray(newphi)
		shutil.copy(rscfile,outrsc)
		
		dst_band1.FlushCache()
		phase_band.FlushCache()
		cor_band.FlushCache()
		del dst_ds, phi, newphi, cor

	else:
		logger.warning('Coherence file does not exit...')

	del ds_int, ds_cor 

def flattenrange(kk, rlook, prefix, suffix, ftype, fit, threshold):
	date1, date2 = date_1[kk], date_2[kk]
	logger.debug('Flatten range on int. {0}-{1}'.format(date1,date2))

	prefixfilt = 'filt' + str(ftype) + '_' + str(prefix)
	suffixrange = suffix + '_flatr'

	corfile = str(date1) + '-' + str(date2) + '_' + rlook + 'rlks.cor'
	intfile = str(prefix) + str(date1) + '-' + str(date2) + str(suffix) + '_' +  rlook + 'rlks.int'
	rscfile = str(prefix) + str(date1) + '-' + str(date2) + str(suffix) + '_' +  rlook + 'rlks.int.rsc'
	filtint = str(prefixfilt) + str(date1) + '-' + str(date2) + str(suffix) + '_' +  rlook + 'rlks.int'

	if os.path.exists(filtint) == False:
		logger.debug('{0} doesnot exist'.format(filtint))
		filter(intfile, ftype, filtint, corfile, window)

	outfile = str(prefix) + str(date1) + '-' + str(date2) + str(suffixrange) + '_' +  rlook + 'rlks.int'
	outrsc = str(prefix) + str(date1) + '-' + str(date2) + str(suffixrange) + '_' +  rlook + 'rlks.int.rsc'
	filtout = str(prefixfilt) + str(date1) + '-' + str(date2) + str(suffixrange) + '_' +  rlook + 'rlks.int'

	r = subprocess.call("flatten_range "+str(infile)+" "+str(filtint)+" "+str(outfile)+" "+str(filtout)+" "+str(fit)+" "+str(threshold), shell=True)
	if r != 0:
		logger.warning("Flatten range failed for int. {0}-{1}".format(date1,date2))
	else:
		shutil.copy(rscfile,outrsc)	

def flattentopo(kk, prefix, suffix):
	date1, date2 = date_1[kk], date_2[kk]
	logger.debug('Flatten topo on int. {0}-{1}'.format(date1,date2))

	prefixfilt = 'filt' + str(ftype) + '_' + str(prefix)


def unwrapping(kk):
	date1, date2 = date_1[kk], date_2[kk]
	return 

def addback(kk):
	date1, date2 = date_1[kk], date_2[kk]
	return 


##################################################################################
###  READ IMPUT PARAMETERS
##################################################################################

# # input parameters 
home='/home/cometraid14/daouts/work/tibet/qinghai/processing/Sentinel/iw1'
int_list='interf_pair.rsc'
master='20160608'
look_int=int(2)
look_unw=int(4)
look=int(look_int - look_unw)
fit_range = -1
threshold_coh = 0.3
ftype='SWc'
SWamplim=0.1
SWwindowsize=8
thresholdfiltSW=0.25
thresholdcol=0.04
seedx=268
seedy=1766
prefix = '' 
suffix = '_sd'


##################################################################################
###  WORK
##################################################################################

if look_int > 2:

	work = [(kk, look_int, prefix, suffix) for kk in xrange(kmax)]
	# pool = multiprocessing.Pool(nproc)
	# pool.map(look, work)
	# pool.close()
	for w in work:
		look(w)

if replace_amp == 'yes':

	work = [(kk, look_int, prefix, suffix) for kk in xrange(kmax)]
	for w in work:
		replaceAmp(w)

if flat_range == 'yes': 
	prefixflat = 'coh_' + preffix

	work = [(kk, look_int, prefixflat, suffix, ftype, fit_range, threshold_coh) for kk in xrange(kmax)]
	# pool = multiprocessing.Pool(nproc)
	# pool.map(flattenrange, work)
	# pool.close()
	for w in work:
		flattenrange(w)

if flat_topo == 'yes':
	pass

if do_unw == 'yes':
	pass

if add_back == 'yes':
	pass

