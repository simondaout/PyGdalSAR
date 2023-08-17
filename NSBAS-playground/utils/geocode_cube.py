#!/usr/bin/env python3
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
[--lectfile=<path>] [--rscfile=<path>] [--smooth=<value>] [--imref=<value>] 

Options:
-h --help           Show this screen.
--infile PATH       path to time series (depl_cumule)
--geomaptrans PATH  path to geomaptrans
--amp PATH          path to amplitude file  
--lectfile PATH     Path of the lect.in file [default: lect.in]
--rscfile PATH      Path to a rsc file [default: radar_4rlks.hgt.rsc]
--smooth VALUE      If >0 smooth time series in time with the value the windows lenght [default: 0]
--imref VALUE       Reference image number [default: 1]
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
from matplotlib import pyplot as plt

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
if arguments["--smooth"] ==  None:
   smooth = 0
else:
   smooth = int(arguments["--smooth"])
if arguments["--imref"] !=  None:
    if int(arguments["--imref"]) < 1:
        print('--imref must be between 1 and Nimages')
    else:
        imref = int(arguments["--imref"]) - 1

# load images_retenues file
fimages='images_retenues'
nb,idates,dates,base=np.loadtxt(fimages, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
# nb images
N=len(dates)
print('Number images: ', N)

# read lect.in 
ncol, nlign = map(int, open(lecfile).readline().split(None, 2)[0:2])

# lect cube
cubei = np.fromfile(infile,dtype=np.float32)
cube = as_strided(cubei[:nlign*ncol*N])
kk = np.flatnonzero(np.logical_or(cube==9990, cube==9999))
cube[kk] = float('NaN')
cube[cube==0] = float('NaN')
maps = cube.reshape((nlign,ncol,N))

del cubei, cube 

# open amp file
amp = np.zeros((nlign,ncol))
fid = open(ampf, 'r')
amp[:nlign,:ncol] = np.fromfile(fid,dtype=np.float32)[:nlign*ncol].reshape((nlign,ncol))

if arguments["--imref"] !=  None:
    cst = np.copy(maps[:,:,imref])
    for l in range((N)):
        maps[:,:,l] = maps[:,:,l] - cst
        if l != imref:
            index = np.nonzero(np.logical_and(maps[:,:,l]==0.0,maps[:,:,l]>999.))
            maps[:,:,l][index] = float('NaN')

for l in range((N)):
#for l in range(10):
    
    if smooth > 0 :
        beg = np.max([0, l - smooth]) 
        end = np.min([l+smooth, N])
        data = np.nanmean(as_strided(maps[:,:,beg:end]),axis=2)
        #plt.imshow(data,vmin=-2,vmax=2)
        #plt.show()
    else: 
        data = as_strided(maps[:,:,l])
    
    drv = gdal.GetDriverByName("roi_pac")
    outfile=str(idates[l])+'.unw'
    outrsc=str(idates[l])+'.unw.rsc'
    geooutfile='geo_'+outfile
    tiffoutfile='geo_'+str(idates[l])+'.tiff'

    print("+rm -f "+outfile+" "+outrsc+" "+tiffoutfile+" "+geooutfile)
    r = subprocess.call("rm -f "+outfile+" "+outrsc+" "+tiffoutfile+" "+geooutfile, shell=True)

    dst_ds = drv.Create(outfile, ncol, nlign, 2, gdal.GDT_Float32)
    dst_band1 = dst_ds.GetRasterBand(1)
    dst_band2 = dst_ds.GetRasterBand(2)
    dst_band1.WriteArray(amp,0,0)
    dst_band2.WriteArray(data,0,0)
    # shutil.copy2(rscf,outrsc)
    del dst_ds, drv, dst_band1, dst_band2

    print("+cp "+rscf+" "+outrsc)
    r = subprocess.call("cp "+rscf+" "+outrsc, shell=True)

del amp

for l in range((N)):
#for l in range(10):
    outfile=str(idates[l])+'.unw'
    outrsc=str(idates[l])+'.unw.rsc'
    geooutfile='geo_'+outfile
    geolosoutfile='geo_'+str(idates[l])+'.los.grd'
    geomagoutfile='geo_'+str(idates[l])+'.mag.grd'
    geophsoutfile='geo_'+str(idates[l])+'.phs.grd'
    mmoutfile='geo_'+str(idates[l])+'_LOSmm_nan.grd'

    print("+geocode.pl "+geomapf+" "+outfile+" "+geooutfile)
    r = subprocess.call("geocode.pl "+geomapf+" "+outfile+" "+geooutfile,
                        shell=True)
    if r != 0:
        raise Exception("geocode.pl failed")

    print("+rmg2grd.py "+geooutfile)
    r = subprocess.call("rmg2grd.py "+geooutfile,
                        shell=True)
    if r != 0:
        raise Exception("rmg2grd.py failed")
   
    print("+rm -f "+geomagoutfile+" "+geophsoutfile)
    r = subprocess.call("rm -f "+geomagoutfile+" "+geophsoutfile,shell=True)
 
    print("+gmt grdmath "+geolosoutfile+" 0 NAN -1000 MUL = "+mmoutfile)
    r = subprocess.call("gmt grdmath "+geolosoutfile+" 0 NAN -1000 MUL = "+mmoutfile,
                        shell=True)
    if r != 0:
        raise Exception("grdmath failed")
    
    #print("+gdal_translate -ot Float32 -b 2 -co COMPRESS=DEFLATE"+geooutfile+" "+tiffoutfile)
    #r = subprocess.call("gdal_translate -ot Float32 -b 2 -co COMPRESS=DEFLATE "+geooutfile+" "+tiffoutfile,
    #                    shell=True)
    #if r != 0:
    #    raise Exception("gdal_translate failed")

    del ds, ds_band2
