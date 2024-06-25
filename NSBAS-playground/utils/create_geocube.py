#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT 
################################################################################

"""\
create_geocube.py 
-------------
Create a NCDF cube file for the whoel time series

Usage: create_geocube.py [--lectfile=<path>] [--images=<path>] 

Options:
-h --help           Show this screen.
--images PATH       Path to image_retuenues file [default: images_retenues]
"""

print()
print()
print('Author: Simon Daout')
print()
print()

# gdal
from osgeo import gdal
gdal.UseExceptions()
import subprocess

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided

import docopt
arguments = docopt.docopt(__doc__)

if arguments["--images"] ==  None:
   fimages = "images_retenues"
else:
   fimages = arguments["--images"]

# load images_retenues file
nb,idates,dates,base=np.loadtxt(fimages, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
# nb images
N=len(dates)
print('Number images: ', N)

outfile = 'geo_TimeSeries_LOSmm_nan.tiff'
r = subprocess.call("rm -f "+outfile, shell=True)

for l in range((N)):
    drv = gdal.GetDriverByName('GTiff')
    infile = 'geo_'+str(idates[l])+'_LOSmm_nan.tiff'
    print('Read file:', infile) 
    # read data
    ds = gdal.Open(infile)
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()

    if l == 0:
        dst_ds = drv.Create(outfile, ds.RasterXSize, ds.RasterYSize, N, gdal.GDT_Float32)

    # los in band2 for unw file
    ds_band2 = ds.GetRasterBand(1)
    los = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)

    # write data
    dst_band1 = dst_ds.GetRasterBand(l+1)
    dst_band1.WriteArray(los,0,0)
    dst_ds.SetGeoTransform(gt)
    dst_ds.SetProjection(proj)
    
    ds_band2.FlushCache()
    del los

del dst_ds




