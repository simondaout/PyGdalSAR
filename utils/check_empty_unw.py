#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package
# written in Python-Gdal
#
############################################
# Author        : Simon DAOUT (Oxford)
############################################

"""\
check_empty_unw.py
-------------
Check for empty unwrapped ifgs

Usage: check_empty_unw.py --int_list=<file> --int_path=<path> [--prefix=<value>] [--suffix=<value>] [--rlook=<value>]

Options:
-h --help           Show this screen.
--int_list=<file>   Text file containing list of interferograms dates in two colums, $data1 $date2
--int_path=<path>   Path to input interferograms directory
--prefix=<value>    Prefix name $prefix$date1-$date2$suffix_$rlookrlks.unw [default: '']
--suffix=<value>    Suffix name $prefix$date1-$date2$suffix_$rlookrlks.unw [default: '']
--rlook=<value>     look int. $prefix$date1-$date2$suffix_$rlookrlks.unw [default: 0]
"""

import gdal
import numpy as np
import docopt
import warnings, sys
from os import getcwd, path

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
    idate = str(date1) + '-' + str(date2)
    folder =  'int_'+ str(date1) + '_' + str(date2) + '/'
    rscfile=int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.unw.rsc'
    infile=int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.unw'
    name = prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.unw'

    if path.exists(infile) is not False:
        print("Open: {0} in {1}".format(name,folder))
        ds = gdal.Open(infile, gdal.GA_ReadOnly)
        los = ds.GetRasterBand(2).ReadAsArray()
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

ListInterfero = "interf_pair_empty.rsc"
print('Save empty interferogram list:', ListInterfero)
wf = open(ListInterfero, 'w')

for j in range(Nifg):
    los,name = check(j)
    occurences = np.count_nonzero(los != 0)
    #if np.sum(los) == 0:
    if occurences < 500:
        print("File {0} empty!".format(name))
        wf.write("%i  %i\n" % (date_1[j], date_2[j]))
wf.close()

print('----------------------------------')
print()





