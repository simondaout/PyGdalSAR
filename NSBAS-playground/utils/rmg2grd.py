#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################

############################################
# Author        : Simon DAOUT
############################################

"""\
rmg2grd.py
-------------
Convert raster file in rmg format to NetCDF format

Usage: rmg2grd.py <infile>

Options:
-h --help           Show this screen.
--infile PATH       Input complex file
"""

import os, math, sys
import numpy as np
from osgeo import gdal
gdal.UseExceptions()
# docopt (command line parser)
import docopt

def checkinfile(file):
  if os.path.exists(file) is False:
    print("File: {0} not found, Exit!".format(file))
    print("File: {0} not found in {1}, Exit!".format(file,os.getcwd()))
    sys.exit()

def CreateOutfile(outfile, ds ,band):
  driver = gdal.GetDriverByName("NetCDF")
  dst_ds = driver.Create(outfile, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Float32)
  dst_band1 = dst_ds.GetRasterBand(1)
  gt = ds.GetGeoTransform()
  proj = ds.GetProjection()
  dst_band1.WriteArray(band,0,0)
  dst_ds.SetGeoTransform(gt)
  dst_ds.SetProjection(proj)
  dst_band1.FlushCache()
  del dst_ds 

# read arguments
arguments = docopt.docopt(__doc__)
infile = arguments["<infile>"]; checkinfile(infile)
extension = os.path.splitext(infile)[1]
basename = os.path.splitext(infile)[0]
rsc = arguments["<infile>"] + '.rsc' 

try:
    import nsb_utilproc
    wavelength = nsb_utilproc.get_value_rscfile(rsc, "WAVELENGTH")
except:
    print('Cannot read WAVELENGTH in rsc file. Set WAVELENGTH to 0.05546576')
    wavelength = 0.05546576
phs_factor = float(wavelength) / (4*math.pi)
print('RAD to m:', phs_factor)

# Open phiset (image)
ds = gdal.OpenEx(infile, nOpenFlags=gdal.GA_ReadOnly, allowed_drivers=["ROI_PAC"])
ds_band1 = ds.GetRasterBand(1)
ds_band2 = ds.GetRasterBand(2)
# Attributes
print("> Driver:   ", ds.GetDriver().ShortName)
print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
print("> Datatype: ", gdal.GetDataTypeName(ds_band1.DataType))

# read and convert NaN
pha = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
kk = np.nonzero(np.logical_or(pha>9990, pha==0.))
pha[kk] = float('NaN')
mag = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
ds_band1.FlushCache; ds_band2.FlushCache

outfile = basename + '.mag.grd'; CreateOutfile(outfile, ds, mag)
outfile = basename + '.phs.grd'; CreateOutfile(outfile, ds, pha)
outfile = basename + '.los.grd'; CreateOutfile(outfile, ds, pha*phs_factor)

del ds

