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
cpx2rmg.py
-------------
Convert raster file from complex to rmg format

Usage: cpx2rmg.py <infile>  

Options:
-h --help           Show this screen.
--infile PATH       Input complex file
"""

import os
import numpy as np
from osgeo import gdal
# docopt (command line parser)
import docopt

# read arguments
arguments = docopt.docopt(__doc__)
infile = arguments["<infile>"]
extension = os.path.splitext(infile)[1]
basename = os.path.splitext(infile)[0]
outfile = infile + '.tiff'

gdal.UseExceptions()
# Open phiset (image)
ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])
ds_band1 = ds.GetRasterBand(1)
# Attributes
print("> Driver:   ", ds.GetDriver().ShortName)
print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
print("> Datatype: ", gdal.GetDataTypeName(ds_band1.DataType))
gt = ds.GetGeoTransform()
proj = ds.GetProjection()
pha = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
amp = np.absolute(pha)
phi = np.angle(pha)

# create new GDAL imaga
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create(outfile, ds.RasterXSize, ds.RasterYSize, 2, gdal.GDT_Float32)
dst_band1 = dst_ds.GetRasterBand(1)
dst_band2 = dst_ds.GetRasterBand(2)
dst_band1.WriteArray(amp,0,0)
dst_band2.WriteArray(phi,0,0)
dst_ds.SetGeoTransform(gt)
dst_ds.SetProjection(proj)
dst_band1.FlushCache()
dst_band2.FlushCache()
del dst_ds, ds



