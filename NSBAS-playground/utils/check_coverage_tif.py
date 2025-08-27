#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package
# written in Python-Gdal
#
############################################
# Author        : Simon DAOUT (Oxford) / Hugo Watine (CRPG)
############################################

"""\
check_coverage_tiff.py
-------------
Check coverage unwrapped ifgs, can be used for Flatsim ifgs

Usage: check_coverage_tiff.py --int_list=<file> --int_path=<path> [--prefix=<value>] [--suffix=<value>] [--rlook=<value>]

Options:
-h --help           Show this screen.
--int_list=<file>   Text file containing list of interferograms dates in two colums, $data1 $date2
--int_path=<path>   Path to input interferograms directory
--prefix=<value>    Prefix name $prefix$date1_$date2$suffix_$rlookrlks.tiff [default: '']
--suffix=<value>    Suffix name $prefix$date1_$date2$suffix_$rlookrlks.tiff [default: '']
--rlook=<value>     look int. $prefix$date1_$date2$suffix_$rlookrlks.tiff [default: 0]
"""

from osgeo import gdal
import numpy as np
import docopt
import warnings, sys
from os import getcwd, path
from osgeo import gdalconst


warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# read arguments
arguments = docopt.docopt(__doc__)
int_list=arguments["--int_list"]
if arguments["--int_path"] == None:
    int_path='./'
else:
    int_path=arguments["--int_path"] + '/'
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

def check(kk):
    date1, date2 = date_1[kk], date_2[kk]
    idate = str(date1) + '_' + str(date2)
    folder =  'int_'+ str(date1) + '_' + str(date2) + '/'
    infile=int_path + folder + prefix + str(date1) + '_' + str(date2) + suffix + rlook + '.tiff'
    name = prefix + str(date1) + '_' + str(date2) + suffix + rlook + '.tiff'

    if path.exists(infile) is not False:
        print("Open: {0} in {1}".format(name,folder))
        ds = gdal.Open(infile, gdalconst.GA_ReadOnly)
        los = ds.GetRasterBand(1).ReadAsArray()
        del ds
    else:
        print("File: {0} not found in {1}".format(name,folder))
        los = np.zeros((2,2))
    los[np.isnan(los)] = 0.
    return los, name

print('----------------------------------')
print('Check interferograms list:', int_list)
date_1,date_2=np.loadtxt(int_list,comments="#",unpack=True,usecols=(0,1),dtype='i,i')
Nifg=len(date_1)
print("number of interferogram: {}".format(Nifg))
print()

ListInterfero = "interf_pair_coverage.txt"
print('Save empty interferogram list:', ListInterfero)
wf = open(ListInterfero, 'w')

for j in range(Nifg):
  try:
    los,name = check(j)
    size = np.shape(los)[0]* np.shape(los)[1]
    unw = np.count_nonzero(los==0)
    cov = 1 - (unw / size)
  except:
    print("Cannot open {0}: !!!!".format(name))
    cov = 0.0
  print("Unwrapping coverage {0}: {1} ".format(name,cov))
  wf.write("%i  %i %f\n" % (date_1[j], date_2[j], cov))
wf.close()

print('----------------------------------')
print()
