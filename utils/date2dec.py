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
date2dec.py
-------------
convert dates in YYYYMMDD format to numerical format

Usage: date2dec.py [--dates=<values>] [--datefile=<path>] [--outfile=<path>]

Options:
    -h --help           Show this screen
    --dates values      list of values separated by a coma (eg. 20130924,20131012)
    --datefile path     file containing dates in one column
"""

from datetime import datetime
import numpy as np
import docopt

arguments = docopt.docopt(__doc__)

if arguments["--dates"] is not  None:
  dates = map(int,arguments["--dates"].replace(',',' ').split())
elif arguments["--datefile"] is not None:
  infile=file(arguments["--datefile"],'r')
  dates=np.loadtxt(infile,comments="#",dtype='str')
else:
  print('No input dates')
  sys.exit()

print(dates)
# sys.exit(0)

def date2dec(dates):
    dates  = np.atleast_1d(dates)
    times = []
    for date in dates:
        x = datetime.strptime('{}'.format(date),'%Y%m%d')
        #x = datetime.strptime('{}'.format(date),'%Y-%m-%dT%H:%M:%S.%fZ')
        dec = float(x.strftime('%j'))/365.1
        year = float(x.strftime('%Y'))
        # print(date,dec,year)
        times.append(year + dec)
    return times

times = date2dec(dates)
print(times)

if arguments["--outfile"] is not None:
    np.savetxt(arguments["--outfile"], np.vstack(np.array(times).T),fmt=('%f'))

