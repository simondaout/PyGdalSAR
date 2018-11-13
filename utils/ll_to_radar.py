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
ll_to_radar.py
-------------
transform a list of geographic coordinates into radar corrdinates 

Usage: ll_to_radar.py --infile=<path> --outfile=<path> [--latfile=<path>] [--lonfile=<value>] [--precision=<value>] [--lectfile=<path>]

Options:
-h --help           Show this screen.
--infile PATH       path to the text file containing list of geographic coordinates
--outfile PATH      output file with radar coordinates
--latfile PATH      path to lon_erai_4rlks.r4  [default: ./lon_erai_4rlks.r4]
--lonfile PATH      path to lat_erai_4rlks.r4  [default: ./lat_erai_4rlks.r4]  
--lectfile PATH     Path of the lect.in file [default: lect.in]
--precision VALUE   Desired precision (number of decimal ) for matching [default:3]
"""

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided
from decimal import Decimal

import docopt
arguments = docopt.docopt(__doc__)
infile = arguments["--infile"]
outfile = arguments["--outfile"]
if arguments["--latfile"] ==  None:
   latfile = 'lat_erai_4rlks.r4'
else:
   latfile = arguments["--latfile"]
if arguments["--lonfile"] ==  None:
   lonfile = 'lat_erai_4rlks.r4'
else:
   lonfile = arguments["--lonfile"]
if arguments["--lectfile"] ==  None:
   lecfile = "lect.in"
else:
   lecfile = arguments["--lectfile"]
if arguments["--precision"] ==  None:
   prec = 3
else:
   prec = int(arguments["--precision"])

ncol, nlign = map(int, open(lecfile).readline().split(None, 2)[0:2])
fid = open(latfile, 'r')
lat = np.fromfile(fid,dtype=np.float32)
for i in xrange(len(lat)):
    lat[i] = round(lat[i],prec)
fid = open(lonfile, 'r')
lon = np.fromfile(fid,dtype=np.float32)
for i in xrange(len(lon)):
    lon[i] = round(lon[i],prec)

# for i in xrange(len(lon)):
#      if lat[i] > 37.45 and lat[i] < 37.48 and lon[i] > 95.63 and lon[i] < 95.65 :
#          print lat[i], lon[i]

lon = lon.reshape((nlign,ncol))
lat = lat.reshape((nlign,ncol))

list_lat,list_lon=np.loadtxt(infile, comments='#', usecols=(0,1), unpack=True,dtype='f,f')
list_lat,list_lon = np.atleast_1d(list_lat),np.atleast_1d(list_lon)
# print list_lat,list_lon
# Decimal(list_lat[0])
# sys.exit()

ligns,cols = [], []
epsi = 10**(-prec)
for llat,llon in zip(list_lat,list_lon):
    kk = np.nonzero(np.logical_and(np.logical_and(lat>llat-epsi,lat<llat+epsi), \
        np.logical_and(lon>llon-epsi,lon<llon+epsi)))
    print kk
    ligns.append(kk[0][0])
    cols.append(kk[1][0])

print ligns
ligns=np.array(ligns)
cols=np.array(cols)
print ligns

print 'Saving in the output file', outfile
np.savetxt(outfile, np.vstack([ligns,cols]).T, header='ligns cols', fmt='%d')