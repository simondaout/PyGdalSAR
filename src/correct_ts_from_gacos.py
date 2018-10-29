#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT (ISTerre)
################################################################################

"""\
correct_ts_from_gacos.py
-------------
Correct Time Series data from Gacos atmospheric models. 1) convert .ztd files to .tif format, 2) re-project and re-resampel atmospheric models
3) correct data

Usage: correct_ts_from_gacos.py [--cube=<path>] [--path=<path>] [--list_images=<path>] [--imref=<value>] \
[--gacos2data=<value>] [--proj=<value>] [--ref=<values>] [--zone=<values>] [--plot=<yes/no>] [--load=<True/False>]

correct_ts_from_gacos.py -h | --help

Options:
-h --help           Show this screen.
--cube PATH         Path to the cube of TS displacements file to be corrected from atmo models [default: depl_cumul]
--path PATH         Path to .ztd data (default: ./GACOS/)
--list_images PATH  Path to list images file made of 4 columns containing for each images 1) number 2) date in YYYYMMDD format 3) numerical date 4) perpendicular baseline [default: list_images.txt] 
--imref VALUE       Reference image number [default: 1]
--ref  VALUES       Column and line number for referencing. If None then estimate the best linear relationship between model and data [default: None]. 
--zone VALUES       Crop option for ramp estimation [default: 0,ncol,0,nlign]
--proj VALUE        EPSG for projection GACOS map [default: 4326]
--gacos2data        Scaling value between gacos data (m) and desired output (e.g data in mm) [default: 1000.]
--plot              Display results [default: yes]   
--load              If False, do not load data again and directly read cube_gacos [default: True]   
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

from nsbas import docopt

# read arguments
arguments = docopt.docopt(__doc__)

if arguments["--path"] ==  None:
   path = "./GACOS/"
else:
   path = arguments["--path"]
if arguments["--cube"] ==  None:
    cubef = "depl_cumule"
else:
    cubef = arguments["--cube"]
if arguments["--list_images"] ==  None:
    listim = "list_images.txt"
else:
    listim = arguments["--list_images"]
if arguments["--imref"] ==  None:
    imref = 0
else:
    imref = int(arguments["--imref"]) - 1
if arguments["--plot"] ==  None:
    plot = 'yes'
else:
    plot = arguments["--plot"]
if arguments["--ref"] is None:
    ref = 'ramp'
else:
    ref = 'pixel'
    ref_col,ref_line = map(int,arguments["--ref"].replace(',',' ').split())
if arguments["--proj"] ==  None:
    proj = False
    EPSG = 4326
else:
    proj = True
    EPSG = int(arguments["--proj"])
if arguments["--gacos2data"] ==  None:
    gacos2data = 1000.
else:
    gacos2data = float(arguments["--gacos2data"])

if arguments["--load"] ==  None:
    load = True
else:
    load = arguments["--load"] 


# read lect.in: size maps
ncol, nlign = map(int, open('lect.in').readline().split(None, 2)[0:2])
if arguments["--zone"] ==  None:
    refzone = [0,ncol,0,nlign]
else:
    #refzone = [col_beg,col_end,line_beg,line_end]
    refzone = map(float,arguments["--clean"].replace(',',' ').split())

nb,idates,dates,base=np.loadtxt(listim, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
N = len(dates)

# Choose plot option
cmap = cm.gist_rainbow_r
cmap.set_bad('white')

# Choose referencing between 'pixel' or 'ramp'
# ref_col,ref_line = 910,1180 502540
# ref_col,ref_line = 271, 1691 503540

# Reprojection option
# EPSG = 32645

# load cube of displacements
cubei = np.fromfile('depl_cumule',dtype=np.float32)
cube = as_strided(cubei[:nlign*ncol*N])
print 'Number of line in the cube: ', cube.shape
maps = cube.reshape((nlign,ncol,N))

if load == True:
    gacos = np.zeros((nlign,ncol,N))
    for i in xrange((N)):

        print 'Read ',idates[i]
        infile = path+'{}.ztd'.format(int(idates[i]))
        rsc = path+'{}.ztd.rsc'.format(int(idates[i]))
        temp1 = path+'{}_temp1.tif'.format(int(idates[i]))
        outtif = path+'{}_gacos.tif'.format(int(idates[i]))

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
        ds = gdal.Open(path+'{}_gacos.tif'.format(int(idates[i])), gdal.GA_ReadOnly)
        # print ncol, nlign
        # print ds.RasterXSize, ds.RasterYSize
        band = ds.GetRasterBand(1)
        gacos[:,:,i] = band.ReadAsArray()*gacos2data
        del ds,band

    # Ref atmo models to the reference image
    cst = np.copy(gacos[:,:,imref])
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
    if plot == 'yes':
        plt.show()

    # save new cube
    fid = open('cube_gacos', 'wb')
    gacos.flatten().astype('float32').tofile(fid)

# # load gacos cube
gcubei = np.fromfile('cube_gacos',dtype=np.float32)
# gcube = as_strided(gcubei[:nlign*ncol*N])
gacos = gcubei.reshape((nlign,ncol,N))

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

    elif ref == 'pixel':
        model = model - np.nanmean(model[ref_line-2:ref_line+2,ref_col-2:ref_col+2])
        data = data - np.nanmean(data[ref_line-2:ref_line+2,ref_col-2:ref_col+2])

    # correction
    model[model==0.] = 0. 
    data_flat[:,:] = data - model
    data_flat[np.isnan(data)] = np.float('NaN')

    if plot == 'yes':
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













