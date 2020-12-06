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
check_event.py
-------------
Check for event date in YYYYMMDD format in a list of ifgs

Usage: check_event.py --int_list=<file> --date=<value> --out_list=<file> 

Options:
-h --help           Show this screen.
--int_list=<file>   Text file containing list of interferograms dates in two colums, $data1 $date2
--date=<value>      Event date in YYYYMMDD format
--out_list=<file>    Output list name containing event 
"""

import gdal
import numpy as np
import docopt
import warnings, sys
from os import getcwd, path
from datetime import datetime

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# read arguments
arguments = docopt.docopt(__doc__)
int_list=arguments["--int_list"]
date=arguments["--date"]

def date2dec(dates):
    dates  = np.atleast_1d(dates)
    times = []
    for date in dates:
        x = datetime.strptime('{}'.format(date),'%Y%m%d')
        #x = datetime.strptime('{}'.format(date),'%Y-%m-%dT%H:%M:%S.%fZ')
        dec = float(x.strftime('%j'))/365.1
        year = float(x.strftime('%Y'))
        # print date,dec,year
        times.append(year + dec)
    return times

# convert event date to dec
event=date2dec(date)

print('----------------------------------')
print('Check interferograms list:', int_list)
date_1,date_2=np.loadtxt(int_list,comments="#",unpack=True,usecols=(0,1),dtype='i,i')
Nifg=len(date_1)
print("number of interferogram: {}".format(Nifg))
print()

ListInterfero = arguments["--out_list"]
print('Save empty interferogram list:', ListInterfero)
wf = open(ListInterfero, 'w')
for j in range(Nifg):
        date1 = date2dec(date_1[j])
        date2 = date2dec(date_2[j])
        if (event>date1) and (event<date2):
            print(date_1[j], date_2[j])
            wf.write("%i  %i\n" % (date_1[j], date_2[j]))
wf.close()

print('----------------------------------')
print()





