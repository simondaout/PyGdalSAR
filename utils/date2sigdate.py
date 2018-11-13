#!/usr/bin/env python2
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
date2sigdate.py
-------------

Usage: date2sigdate.py --input=<path> --output=<path> --nsigma=<value>
date2sigdate.py -h | --help

Options:
-h --help           Show this screen.
--input PATH        Infile containing liste of dates.
--output PATH       Outfile containing liste of dates bounding at + and- nsigma the orignal dates  
--nsigma VALUE      define number of date    

"""
# docopt (command line parser)
import docopt

import numpy as np
import math
from datetime import datetime
from datetime import timedelta

# read arguments
arguments = docopt.docopt(__doc__)
infile = arguments["--input"]
outfile = arguments["--output"]
nsigma = abs(int(arguments["--nsigma"]))

#Data loading
print "dates list in: ",infile
source1=file(infile,'r')
dates=np.loadtxt(source1,comments="#",dtype='str')
imax=len(dates)
print "number of dates: ",imax

sigdates = [] 
for date in dates:
    x = datetime.strptime('{}'.format(date),'%Y%m%d')
    print 
    for dx in xrange(-nsigma,1+nsigma):
        print (x + timedelta(days=dx)).strftime('%Y%m%d')
        sigdates.append(int((x + timedelta(days=dx)).strftime('%Y%m%d')))

print 'Saving in the output file', outfile
np.savetxt(outfile, np.vstack(np.array(sigdates).T),fmt=('%i'))
#print np.array(sigdates).T
