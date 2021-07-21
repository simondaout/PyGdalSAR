#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""\
netcdf2kite.py
-------------
Reading geocoded los maps in netcdf format  

Usage: netcdf2kite.py --infile=<path> [--lookfile=<path>] [--headingfile=<path>] \
[--headingm=<value>] [--lookm=<value>] [--id=<value>] 
Options:
-h --help			      Show this screen
--infile PATH 	    Netcdf LOS map (.grd or .unw)
--lookfile PATH     Netcdf Incidence map (incidence in .grd format or incidence+heading in .unw format)
--headingfile PATH  Average Heading angle (in .grd format)
--headingm VALUE			Average Heading angle
--lookm VALUE			  Average Incidence angle
--id value          Add id to kite scene
"""

import sys,os
import numpy as np

from osgeo import gdal
from kite import Scene

# docopt (command line parser http://docopt.org)
from docopt import docopt

arguments = docopt(__doc__, argv=None, help=True, version=None, options_first=False)
infile = arguments["--infile"] 

if arguments["--lookfile"] == None:
   lookfile = None
   if arguments["--lookm"] ==  None:
     lookm = 20.
   else:
     lookm = float(arguments["--lookm"])
   print('Set look angle to average value:', lookm)
else:
   lookfile = arguments["--lookfile"]

if arguments["--headingfile"] == None:
   headingfile = None
   if arguments["--headingm"] ==  None:
     headingm = -76
     
   else:
     headingm = float(arguments["--headingm"])
   print('Set heading angle to average value:', headingm)
else:
   headingfile = arguments["--headingfile"]


### Open infile
ds_name = os.path.splitext(infile)[0]
ds_extension = os.path.splitext(infile)[1]
gdal.UseExceptions()
ds = gdal.Open(infile,gdal.GA_ReadOnly)
print("> Driver:   ", ds.GetDriver().ShortName)
print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
ds_geo = ds.GetGeoTransform()
print('Coordinates:', ds_geo)
print 


if ds_extension == ".unw":
  phs_band = ds.GetRasterBand(2)
  los = np.flipud(phs_band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize))
  # los = phs_band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
else:
  phs_band = ds.GetRasterBand(1)
  los = np.flipud(phs_band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize))
  # los = phs_band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)

# also need to flip for ascending data?
print("> Datatype: ", gdal.GetDataTypeName(phs_band.DataType))

### Open incidence file
if lookfile is not None:
  ds2 = gdal.Open(lookfile,gdal.GA_ReadOnly)
  ds2_extension = os.path.splitext(lookfile)[1]
  print("> Driver:   ", ds2.GetDriver().ShortName)
  print("> Size:     ", ds2.RasterXSize,'x',ds2.RasterYSize,'x',ds2.RasterCount)
  
  if ds2_extension == ".unw":
    # ROI_PAC: look angle in band1, heading in band 2 
    look = np.flipud(ds2.GetRasterBand(1).ReadAsArray\
    (0, 0, ds.RasterXSize, ds.RasterYSize))
    # look = ds2.GetRasterBand(1).ReadAsArray\
    # (0, 0, ds2.RasterXSize, ds2.RasterYSize)
    heading = np.flipud(ds2.GetRasterBand(2).ReadAsArray\
    (0, 0, ds2.RasterXSize, ds2.RasterYSize))
    # heading = ds2.GetRasterBand(2).ReadAsArray\
    # (0, 0, ds2.RasterXSize, ds2.RasterYSize)
    print('Found heading angle in .unw file')
  
  else:   
     look = np.flipud(ds2.GetRasterBand(1).ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize))
     # look = ds2.GetRasterBand(1).ReadAsArray(0, 0, ds2.RasterXSize, ds2.RasterYSize)
     ### Open heading file
     if headingfile is not None:
       ds3 = gdal.Open(headingfile,gdal.GA_ReadOnly)
       ds3_extension = os.path.splitext(headingfile)[1]
       print("> Driver:   ", ds3.GetDriver().ShortName)
       print("> Size:     ", ds3.RasterXSize,'x',ds3.RasterYSize,'x',ds3.RasterCount)
       # heading = ds3.GetRasterBand(1).ReadAsArray(0, 0, ds3.RasterXSize, ds3.RasterYSize)
       heading = np.flipud(ds3.GetRasterBand(1).ReadAsArray(0, 0, ds3.RasterXSize, ds3.RasterYSize))
     else:
       heading = np.ones((ds.RasterYSize,ds.RasterXSize))*headingm
       print('Set heading angle to default Envisat value {}'.format(headingm))
       
else:  
  look = np.ones((ds.RasterYSize,ds.RasterXSize))*lookm
  print('Set look angle to default Envisat value {}'.format(lookm))
  if headingfile is not None:
     ds3 = gdal.Open(headingfile,gdal.GA_ReadOnly)
     ds3_extension = os.path.splitext(headingfile)[1]
     print("> Driver:   ", ds3.GetDriver().ShortName)
     print("> Size:     ", ds3.RasterXSize,'x',ds3.RasterYSize,'x',ds3.RasterCount)
     # heading = ds3.GetRasterBand(1).ReadAsArray(0, 0, ds3.RasterXSize, ds3.RasterYSize)
     heading = np.flipud(ds3.GetRasterBand(1).ReadAsArray(0, 0, ds3.RasterXSize, ds3.RasterYSize))
  else:
     heading = np.ones((ds.RasterYSize,ds.RasterXSize))*headingm
     print('Set heading angle to default Envisat value {}'.format(headingm))

print('Av. Heading:', np.nanmean(heading))
print('Av Look:', np.nanmean(look))

theta = np.deg2rad(90.-look)
phi=np.ones((ds.RasterYSize,ds.RasterXSize))*np.deg2rad(-90-heading)

print('Av. theta:', np.rad2deg(np.nanmean(theta)))
print('Av phi:', np.rad2deg(np.nanmean(phi)))

# sys.exit()
#los[np.isnan(los)]=0.0
#theta[np.isnan(theta)]=0.0
#phi[np.isnan(phi)]=0.0

sc = Scene()
sc.displacement = -los
# !!! lower left corner !!! 
sc.frame.llLat = ds_geo[3] + ds_geo[5] * ds.RasterYSize
sc.frame.llLon = ds_geo[0] 

sc.frame.dN = -ds_geo[5]
sc.frame.dE = ds_geo[1]
sc.frame.spacing = 'degree'

sc.theta = theta
sc.phi = phi

sc.meta.scene_id = arguments["--id"]

#get ref point
print(sc.frame.llEutm, sc.frame.llNutm)
print(sc.frame.utm_zone)


sc.spool()
sys.exit()

# sc.save('kite_scene')

qt = sc.quadtree  # this initializes and constructs the tree
qt.epsilon = 0.000444
qt.nan_allowed = .5
qt.tile_size_min = 2500
qt.tile_size_max = 4500
# sc.spool()

qt.export('{}_tree.csv'.format(ds_name))

