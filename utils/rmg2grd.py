#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package
# written in Python-Gdal
#
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

import os, math
import numpy as np
from osgeo import gdal
gdal.UseExceptions()
# docopt (command line parser)
import docopt

# read arguments
arguments = docopt.docopt(__doc__)
infile = arguments["<infile>"]
extension = os.path.splitext(infile)[1]
basename = os.path.splitext(infile)[0]
rsc = arguments["<infile>"] + '.rsc' 

try:
    import nsb_utilproc
    wavelength = nsb_utilproc.get_value_rscfile(rsc, "WAVELENGTH")
    phs_factor = wavelength / (4*math.pi)
except:
    phs_factor = 1

# Open phiset (image)
ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])
ds_band1 = ds.GetRasterBand(1)
ds_band2 = ds.GetRasterBand(2)
# Attributes
print("> Driver:   ", ds.GetDriver().ShortName)
print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
print("> Datatype: ", gdal.GetDataTypeName(ds_band1.DataType))
gt = ds.GetGeoTransform()
proj = ds.GetProjection()

# read and convert NaN
pha = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
kk = np.flatnonzero(np.logical_or(pha>9990, pha==0.))
pha[kk] = float('NaN')
mag = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)

outfile = basename + '.mag.grd'
driver = gdal.GetDriverByName("NetCDF")
dst_ds = driver.Create(outfile, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Float32)
dst_band1 = dst_ds.GetRasterBand(1)
dst_band1.WriteArray(mag,0,0)
dst_ds.SetGeoTransform(gt)
dst_ds.SetProjection(proj)
dst_band1.FlushCache()
del dst_ds 

outfile = basename + '.phs.grd'
driver = gdal.GetDriverByName("NetCDF")
dst_ds = driver.Create(outfile, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Float32)
dst_band1 = dst_ds.GetRasterBand(1)
dst_band1.WriteArray(pha,0,0)
dst_ds.SetGeoTransform(gt)
dst_ds.SetProjection(proj)
dst_band1.FlushCache()
del dst_ds

outfile = basename + '_nan_mmyr.grd'
driver = gdal.GetDriverByName("NetCDF")
dst_ds = driver.Create(outfile, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Float32)
dst_band1 = dst_ds.GetRasterBand(1)
dst_band1.WriteArray(pha*phs_factor,0,0)
dst_ds.SetGeoTransform(gt)
dst_ds.SetProjection(proj)
dst_band1.FlushCache()
del dst_ds, ds


