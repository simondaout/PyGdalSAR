#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT 
################################################################################

"""\
geocode_cube.py 
-------------
Geocode cube of cumulative deplacements: create date.unw, geo_date.unw, and geo_date.tiff for each dates
!!! Need geocode.pl from ROI_PAC

Usage: geocode_cube.py --cube=<path> --geomaptrans=<path> --amp=<path> \
[--lectfile=<path>] [--rscfile=<path>] 

Options:
-h --help           Show this screen.
--infile PATH       path to time series (depl_cumule)
--geomaptrans PATH  path to geomaptrans
--amp PATH          path to amplitude file  
--lectfile PATH     Path of the lect.in file [default: lect.in]
--rscfile PATH      Path to a rsc file [default: radar_4rlks.hgt.rsc]
"""

print()
print()
print('Please cite:')
print ('Daout, S., Doin, M. P., Peltzer, G., Socquet, A., & Lasserre, C. (2017). Large‐scale InSAR monitoring of permafrost freeze‐thaw cycles on the Tibetan Plateau. Geophysical Research Letters, 44(2), 901-909.')
print('Daout, S., Sudhaus, H., Kausch, T., Steinberg, A., & Dini, B. (2019). Interseismic and postseismic shallow creep of the North Qaidam Thrust faults detected with a multitemporal InSAR analysis. Journal of Geophysical Research: Solid Earth, 124(7), 7259-7279.')
print()
print()

# gdal
from osgeo import gdal
gdal.UseExceptions()

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided

import subprocess, shutil, sys, os

import docopt
arguments = docopt.docopt(__doc__)
infile = arguments["--cube"]
geomapf = arguments["--geomaptrans"]
ampf = arguments["--amp"]
if arguments["--lectfile"] ==  None:
   lecfile = "lect.in"
else:
   lecfile = arguments["--lectfile"]
if arguments["--rscfile"] ==  None:
   rscf = "radar_4rlks.hgt.rsc"
else:
   rscf = arguments["--rscfile"]

# load images_retenues file
fimages='images_retenues'
nb,idates,dates,base=np.loadtxt(fimages, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
# nb images
N=len(dates)
print 'Number images: ', N

# read lect.in 
ncol, nlign = map(int, open(lecfile).readline().split(None, 2)[0:2])

# lect cube
cubei = np.fromfile(infile,dtype=np.float32)
cube = as_strided(cubei[:nlign*ncol*N])
maps = cube.reshape((nlign,ncol,N))

del cubei, cube 

# open amp file
amp = np.zeros((nlign,ncol))
fid = open(ampf, 'r')
amp[:nlign,:ncol] = np.fromfile(fid,dtype=np.float32)[:nlign*ncol].reshape((nlign,ncol))

for l in xrange((N)):
# for l in xrange(3):
    data = as_strided(maps[:,:,l])
    drv = gdal.GetDriverByName("roi_pac")
    
    outfile=str(idates[l])+'_'+str(l)+'.unw'
    outrsc=str(idates[l])+'_'+str(l)+'.unw.rsc'
    geooutfile='geo_'+outfile
    tiffoutfile='geo_'+str(idates[l])+'_'+str(l)+'.tiff'

    print "+rm -f "+outfile+" "+outrsc+" "+tiffoutfile+" "+geooutfile
    r = subprocess.call("rm -f "+outfile+" "+outrsc+" "+tiffoutfile+" "+geooutfile, shell=True)

    dst_ds = drv.Create(outfile, ncol, nlign, 2, gdal.GDT_Float32)
    dst_band1 = dst_ds.GetRasterBand(1)
    dst_band2 = dst_ds.GetRasterBand(2)
    dst_band1.WriteArray(amp,0,0)
    dst_band2.WriteArray(data,0,0)
    # shutil.copy2(rscf,outrsc)
    del dst_ds, drv, dst_band1, dst_band2

    print "+cp "+rscf+" "+outrsc
    r = subprocess.call("cp "+rscf+" "+outrsc, shell=True)

del amp

for l in xrange((N)):
# for l in xrange(3):
    outfile=str(idates[l])+'_'+str(l)+'.unw'
    outrsc=str(idates[l])+'_'+str(l)+'.unw.rsc'
    geooutfile='geo_'+outfile
    tiffoutfile='geo_'+str(idates[l])+'_'+str(l)+'.tiff'

    print "+geocode.pl "+geomapf+" "+outfile+" "+geooutfile
    r = subprocess.call("geocode.pl "+geomapf+" "+outfile+" "+geooutfile,
                        shell=True)
    if r != 0:
        raise Exception("geocode.pl failed")

    print "+gdal_translate -ot Float32 -b 2 -co COMPRESS=DEFLATE -co COMPRESS=PREDICTOR "+geooutfile+" "+tiffoutfile
    r = subprocess.call("gdal_translate -ot Float32 -b 2 -co COMPRESS=DEFLATE -co COMPRESS=PREDICTOR "+geooutfile+" "+tiffoutfile,
                        shell=True)
    if r != 0:
        raise Exception("gdal_translate failed")

    ds = gdal.Open(geooutfile, gdal.GA_ReadOnly)
    ds_band2 = ds.GetRasterBand(2)
    los = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
    print 'Nlign:{}, Ncol:{}, geodate:{}:'.format(ds.RasterYSize, ds.RasterXSize, l)
    nlign,ncol=ds.RasterYSize, ds.RasterXSize
    # if l==0:
    #     geomaps=np.zeros((nlign,ncol,N))
    
    # geomaps[:,:,l] = los
    
    del ds, ds_band2

fid = open('lect_geo.in','w')
np.savetxt(fid, (nlign,ncol,N),fmt='%6i',newline='\t')
fid.close()



# fid = open('geo_depl_cumule', 'wb')
# geomaps.ravel().astype('float32').tofile(fid)

# import matplotlib as mpl
# from matplotlib import pyplot as plt
# import matplotlib.cm as cm
# from pylab import *
# fig = plt.figure(0,figsize=(6,4))
# ax = fig.add_subplot(1,1,1)
# cax = ax.imshow(geomaps[:,:,2],cmap=cm.jet,vmax=10,vmin=-10)
# plt.show()
