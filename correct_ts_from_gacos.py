#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
@author: Simon Daout
"""

import gdal
from osgeo import osr
import sys,os
import numpy as np
from scipy.ndimage import zoom
import scipy.optimize as opt
import scipy.linalg as lst
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from numpy.lib.stride_tricks import as_strided
import subprocess
np.warnings.filterwarnings('ignore')

# read lect.in: size maps
ncol, nlign = map(int, open('lect.in').readline().split(None, 2)[0:2])
# Load list images file made of 4 columns containing for each images 1) number 2) date in YYYYMMDD format 
# 3) numerical date 4) perpendicular baseline
nb,idates,dates,base=np.loadtxt('list_images.txt', comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
N = len(dates)

# Choose plot option
plot = True
cmap = cm.gist_rainbow_r
cmap.set_bad('white')

# convert m to your format
gacos2data = 1000.

# Choose referencing between 'pixel' or 'ramp'
# refx,refy = 180,1658
ref = 'ramp'

# Reprojection option
proj = True
EPSG = 32645
# refzone = [col_beg,col_end,line_beg,line_end]
# refzone =[800,1500,300,1400]

# load cube of displacements
cubei = np.fromfile('depl_cumule',dtype=np.float32)
cube = as_strided(cubei[:nlign*ncol*N])
print 'Number of line in the cube: ', cube.shape
maps = cube.reshape((nlign,ncol,N))


gacos = np.zeros((nlign,ncol,N))
for i in xrange((N)):

    print 'Read ',idates[i]
    infile = 'GACOS/{}.ztd'.format(int(idates[i]))
    rsc = 'GACOS/{}.ztd.rsc'.format(int(idates[i]))
    temp1 = 'GACOS/{}_temp1.tif'.format(int(idates[i]))
    outtif = 'GACOS/{}_gacos.tif'.format(int(idates[i]))

    # read .rsc
    nncol = np.int(open(rsc).readlines()[0].split(None, 2)[1])
    nnlign = np.int(open(rsc).readlines()[1].split(None, 2)[1])
    xfirst = np.float(open(rsc).readlines()[6].split(None, 2)[1])
    yfirst = np.float(open(rsc).readlines()[7].split(None, 2)[1])
    xstep = np.float(open(rsc).readlines()[8].split(None, 2)[1])
    ystep = np.float(open(rsc).readlines()[9].split(None, 2)[1])
    geotransform = (xfirst, xstep, 0, yfirst, 0, ystep)
    print 'Read .rsc and set Geotransform', geotransform
    ztd_ = np.fromfile(infile, dtype='float32').reshape(nnlign,nncol)
    
    # transform gacos format to samething humanily readeable
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(temp1, nncol, nnlign, 1, gdal.GDT_Float32)
    band = ds.GetRasterBand(1)
    band.WriteArray(ztd_)
    ds.SetGeoTransform(geotransform)
    srs = osr.SpatialReference() 
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt()) 
    band.FlushCache()
    del ds,band

    # re-preoject and resample 
    if proj:
        print "gdalwarp -overwrite -s_srs EPSG:4326 -t_srs EPSG:"+str(EPSG)+" -ts "+str(ncol)+" "+str(nlign)+" -r average "+temp1+" "+outtif+" -of GTiff"
        r = subprocess.call("gdalwarp -overwrite -s_srs EPSG:4326 -t_srs EPSG:"+str(EPSG)+" -ts "+str(ncol)+" "+str(nlign)+" -r average "+temp1+" "+outtif+" -of GTiff", shell=True)
        if r != 0:
            raise Exception("gdalwarp failed")
    else:
        print "gdalwarp -overwrite -ts "+str(nlign)+" "+str(ncol)+" -r average "+temp1+" "+outtif+" -of GTiff"
        r = subprocess.call("gdalwarp -overwrite -ts "+str(nlign)+" "+str(ncol)+" -r average "+temp1+" "+outtif+" -of GTiff", shell=True)
        if r != 0:
            raise Exception("gdalwarp failed")

    # del temporary files
    cmd = 'rm -f {}'.format(temp1)
    os.system(cmd)

    # open gacos
    ds = gdal.Open('GACOS/{}_gacos.tif'.format(int(idates[i])), gdal.GA_ReadOnly)
    # print ncol, nlign
    # print ds.RasterXSize, ds.RasterYSize
    band = ds.GetRasterBand(1)
    gacos[:,:,i] = band.ReadAsArray()*gacos2data
    del ds,band

# Ref to the first image
cst = np.copy(gacos[:,:,0])
for l in xrange((N)):
    d = as_strided(gacos[:,:,l])
    gacos[:,:,l] = gacos[:,:,l] - cst

# Plot
fig = plt.figure(0,figsize=(14,10))
fig.subplots_adjust(wspace=0.001)
vmax = np.nanpercentile(gacos[:,:,-1],99)
vmin = np.nanpercentile(gacos[:,:,-1],1)
for l in xrange((N)):
    d = as_strided(gacos[:,:,l])
    ax = fig.add_subplot(4,int(N/4)+1,l+1)
    cax = ax.imshow(d,cmap=cmap,vmax=vmax,vmin=vmin)
    ax.set_title(idates[l],fontsize=6)
    plt.setp( ax.get_xticklabels(), visible=False)
    plt.setp( ax.get_yticklabels(), visible=False)

plt.suptitle('GACOS models')
fig.colorbar(cax, orientation='vertical',aspect=10)
fig.savefig('gocos.eps', format='EPS',dpi=150)
if plot:
    plt.show()

# Apply correction
maps_flat = np.zeros((nlign,ncol,N))
for i in xrange(1,N):   

    data = as_strided(maps[:,:,i])
    data_flat = as_strided(maps_flat[:,:,i])
    model = as_strided(gacos[:,:,i])

    if ref == 'ramp':
        pix_lin, pix_col = np.indices((nlign,ncol))
        try:
            col_beg,col_end,line_beg,line_end = refzone[0],refzone[1],refzone[2],refzone[3]
        except:
            col_beg,col_end,line_beg,line_end = 0 , ncol, 0., nlign

        # find proportionality between data and model
        index = np.nonzero(
            np.logical_and(~np.isnan(data),
            np.logical_and(pix_lin>line_beg, 
            np.logical_and(pix_lin<line_end,
            np.logical_and(pix_col>col_beg, pix_col<col_end
            )))))

        temp = np.array(index).T
        x = temp[:,0]; y = temp[:,1]
        los_clean = data[index].flatten()
        model_clean = model[index].flatten()

        G=np.zeros((len(los_clean),2))
        G[:,0] = 1
        G[:,1] = model_clean
        x0 = lst.lstsq(G,los_clean)[0]
        _func = lambda x: np.sum(((np.dot(G,x)-los_clean))**2)
        pars = opt.least_squares(_func,x0,jac='3-point',loss='cauchy').x
        a = pars[0]
        b = pars[1]
        print 'ref frame %f + %f model for date: %i'%(a,b,idates[i])
        model = a + b*model

    if ref == 'pixel':
        model = model - model[refy,refx]
        data = data - data[refy,refx]

    # correction
    model[model==0.] = 0. 
    data_flat[:,:] = data - model
    data_flat[np.isnan(data)] = np.float('NaN')

    if plot:
        vmax = np.nanpercentile(data,98)
        vmin = np.nanpercentile(data,2)
        # initiate figure depl
        fig = plt.figure(1,figsize=(10,5))
        ax = fig.add_subplot(1,3,1)
        cax = ax.imshow(data,cmap=cmap,vmax=vmax,vmin=vmin)
        ax.set_title('Data {}'.format(idates[i]),fontsize=6)

        # initiate figure depl
        ax = fig.add_subplot(1,3,2)
        cax = ax.imshow(data-data_flat,cmap=cmap,vmax=vmax,vmin=vmin)
        ax.set_title('Model {}'.format(idates[i]),fontsize=6)

        # initiate figure depl
        ax = fig.add_subplot(1,3,3)
        cax = ax.imshow(data_flat,cmap=cmap,vmax=vmax,vmin=vmin)
        ax.set_title('Correct Data {}'.format(idates[i]),fontsize=6)
        fig.colorbar(cax, orientation='vertical',aspect=10)
        plt.show()
        # sys.exit()

# save new cube
fid = open('depl_cumule_gacos', 'wb')
maps_flat.flatten().astype('float32').tofile(fid)













