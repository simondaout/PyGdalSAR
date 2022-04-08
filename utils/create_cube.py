#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from os import environ
import numpy as np
from osgeo import gdal
import sys
import matplotlib.pyplot as plt
import subprocess

def run(cmd):
    """
    Runs a shell command, and print it before running.
    """
    print(cmd)
    r = subprocess.call(cmd, shell=True, stdout=sys.stdout, stderr=subprocess.STDOUT,
        env=environ)
    if r != 0:
        logger.critical(r)
    return

# read lect.in 
lectfile='lect.in'
ncol, nlign = map(int, open(lectfile).readline().split(None, 2)[0:2])
# read dates
fimages='list_images.txt'
nb,idates,dates,base=np.loadtxt(fimages, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
# nb images
N=len(dates)
print('Size expected time series cube: ({0} {1} {2}) '.format(nlign, ncol, N))

imref = 0
inlook = 4
outlook = 4
rlook = int(outlook/inlook)
ref = [1788, 1100]

maps = np.zeros((nlign,ncol,N))
i=0
for d in idates:
    infile = '{}_mdel_{}rlks.unw'.format(d,inlook)
    run("look.pl "+str(infile)+" "+str(rlook))
    infile = '{}_mdel_{}rlks.unw'.format(d,outlook)
    ds = gdal.OpenEx(infile,allowed_drivers=["ROI_PAC"])
    band = ds.GetRasterBand(2)
    print("> Driver:   ", ds.GetDriver().ShortName)
    print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
    print("> Datatype: ", gdal.GetDataTypeName(band.DataType))
    print()
    phi = band.ReadAsArray() 
    maps[:,:,i] = phi[:nlign,:ncol]
    plt.imshow(maps[:,:,i])
    i+=1

# ref to imref
cst = np.copy(maps[:,:,imref])
for l in range((N)):
	# ref in time     
	maps[:,:,l] = maps[:,:,l] - cst
	# ref in space
	# maps[:,:,l] = maps[:,:,l] - maps[ref[0],ref[1],l]

fid = open('cube_era5', 'wb')
maps.flatten().astype('float32').tofile(fid)
plt.show()
