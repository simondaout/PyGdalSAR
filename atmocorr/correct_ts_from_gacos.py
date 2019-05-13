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

Usage: 
    correct_ts_from_gacos.py [--cube=<path>] [--path=<path>] [--list_images=<path>] [--imref=<value>] [--crop=<values>] [--gacos2data=<value>] [--proj=<value>] [--plot=<yes|no>] [--load=<path>] [--rmspixel=<path>] [--threshold_rms=<value>] [--refstart=<values>] [--refend=<values>]  
    correct_ts_from_gacos.py [--cube=<path>] [--path=<path>] [--list_images=<path>] [--imref=<value>] [--crop=<values>] [--gacos2data=<value>] [--proj=<value>] [--plot=<yes|no>] [--load=<path>] [--rmspixel=<path>] [--threshold_rms=<value>] [--refstart=<values>] [--refend=<values>]  [--ramp=<cst|lin>] [--zone=<values>] [--topofile=<path>] 
    correct_ts_from_gacos.py [--cube=<path>] [--path=<path>] [--list_images=<path>] [--imref=<value>] [--crop=<values>] [--gacos2data=<value>] [--proj=<value>] [--plot=<yes|no>] [--load=<path>] [--rmspixel=<path>] [--threshold_rms=<value>] [--refstart=<values>] [--refend=<values>]  [--fitmodel=<yes|no>] [--zone=<values>] [--topofile=<path>] 

correct_ts_from_gacos.py -h | --help

Options:
-h --help             Show this screen.
--cube=<file>         Path to the cube of TS displacements file to be corrected from atmo models [default: depl_cumul]
--path=<dir>          Path to .ztd data [default: ./GACOS/]
--list_images=<file>  Path to list images file made of 5 columns containing for each images 1) number 2) Doppler freq (not read) 3) date in YYYYMMDD format 4) numerical date 5) perpendicular baseline [default: list_images.txt]
--imref=<value>       Reference image number [default: 1]
--proj=<value>        EPSG for projection GACOS map [default: 4326]
--crop=<values>       Crop GACOS data to data extent in the output projection, eg. --crop=xmin,ymin,xmax,ymax 
--refstart=<value>    Stating line number of the area where phase is set to zero [default: 0] 
--refend=<value>      Ending line number of the area where phase is set to zero [default: length]
--ramp=<cst|lin|quad> Estimate a constant, linear or quadratic ramp in x and y in addition to the los/gacos relationship [default: cst].
--zone=<value>        Crop option for empirical estimation (eg: 0,ncol,0,nlign)
--topo=<file>         Use DEM file to mask very low or high altitude values in the empirical estimation 
--gacos2data=<value>  Scaling value between zenithal gacos data (m) and desired output (e.g data in mm positive toward the sat. and LOS angle 23°: -1000*cos(23)=-920.504) [default: -920.504]
--rmspixel=<file>     Path to the RMS map that gives an error for each pixel (e.g RMSpixel, output of invers_pixel)
--threshold_rms=<value>  Thresold on rmspixel for ramp estimation [default: 2.]   
--plot=<yes|no>      Display results [default: yes]
--load=<file>        If a file is given, load directly GACOS cube 
--fitmodel=<yes|no>  If yes, then estimate the proportionlality between gacos and data in addition to a polynomial ramp
"""

import gdal
from osgeo import osr
import sys,os
import numpy as np
import scipy.optimize as opt
import scipy.linalg as lst
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy.lib.stride_tricks import as_strided
import subprocess
np.warnings.filterwarnings('ignore')

try:
    from nsbas import docopt
except:
    import docopt
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
if arguments["--proj"] ==  None:
    proj = False
    EPSG = 4326
else:
    proj = True
    EPSG = int(arguments["--proj"])
if arguments["--gacos2data"] ==  None:
    gacos2data = -920.504 
else:
    gacos2data = float(arguments["--gacos2data"])
if arguments["--ramp"] ==  None:
    ramp = 'cst'
else:
    ramp = arguments["--ramp"]
if arguments["--load"] ==  None:
    load = 'yes'
    loadf = 'cube_gacos'
else:
    load = 'no'
    loadf = arguments["--load"]
if arguments["--topofile"] ==  None:
   radar = None
else:
   radar = arguments["--topofile"]
# read lect.in: size maps
ncol, nlign = map(int, open('lect.in').readline().split(None, 2)[0:2])
if arguments["--zone"] ==  None:
    refzone = [0,ncol,0,nlign]
    col_beg,col_end,line_beg,line_end = 0 , ncol, 0., nlign
else:
    refzone = map(float,arguments["--clean"].replace(',',' ').split())

if arguments["--refstart"] == None:
    refstart = 0
else:
    refstart = int(arguments["--refstart"])
if arguments["--refend"] == None:
    refend = nlign
else:
    refend = int(arguments["--refend"])
   col_beg,col_end,line_beg,line_end = refzone[0],refzone[1],refzone[2],refzone[3]

if arguments["--rmspixel"] ==  None:
    rms = np.ones((nlign,ncol))
else:
    rmsf = arguments["--rmspixel"]
    rms = np.fromfile(rmsf,dtype=np.float32).reshape((nlign,ncol))
    plt.imshow(rms)
    #plt.show()

if arguments["--threshold_rms"] ==  None:
    threshold_rms = 2.
else:
    threshold_rms = float(arguments["--threshold_rms"])

if arguments["--fitmodel"] ==  None:
    fitmodel = 'no'
else:
    fitmodel = arguments["--fitmodel"]

if arguments["--crop"] ==  None:
    crop = False
else:
    crop = map(float,arguments["--crop"].replace(',',' ').split())
nb,idates,dates,base=np.loadtxt(listim, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
N = len(dates)

# Choose plot option
cmap = cm.gist_rainbow_r
# cmap = cm.jet
cmap.set_bad('white')

if radar is not None:
    extension = os.path.splitext(radar)[1]
    if extension == ".tif":
      ds = gdal.Open(radar, gdal.GA_ReadOnly)
      band = ds.GetRasterBand(1)
      elev = band.ReadAsArray()
      del ds
    else:
      fid = open(radar,'r')
      elevi = np.fromfile(fid,dtype=np.float32)
      elevi = elevi[:nlign*ncol]
      elev = elevi.reshape((nlign,ncol))
      fid.close()
else:
    elev = np.ones((nlign,ncol))

# load cube of displacements
cubei = np.fromfile(cubef,dtype=np.float32)
cube = as_strided(cubei[:nlign*ncol*N])
print 'Number of line in the cube: ', cube.shape
maps = cube.reshape((nlign,ncol,N))

# ini
nfigure = 0

if load == 'yes':
    gacos = np.zeros((nlign,ncol,N))
    for i in xrange((N)):

        print 'Read ',idates[i], i
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

        # crop, re-preoject and resample
        if proj and crop is not False:
            print "gdalwarp -overwrite -s_srs EPSG:4326 -t_srs EPSG:"+str(EPSG)+\
            " -ts "+str(ncol)+" "+str(nlign)+" -r average -te "+str(crop[0])+" "+str(crop[1])+" "+str(crop[2])+" "+str(crop[3])+" "\
            +temp1+" "+outtif+" -of GTiff"
            r = subprocess.call("gdalwarp -overwrite -s_srs EPSG:4326 -t_srs EPSG:"+str(EPSG)+\
                " -te "+str(crop[0])+" "+str(crop[1])+" "+str(crop[2])+" "+str(crop[3])+" "\
                +" -ts "+str(ncol)+" "+str(nlign)+" -r average "+temp1+" "+outtif+" -of GTiff", shell=True)
            if r != 0:
                raise Exception("gdalwarp failed")

        elif crop is not False and proj is False:
            print "gdalwarp -overwrite -s_srs EPSG:4326-ts "+str(ncol)+" "+str(nlign)+" -r average \
            -te "+str(crop[0])+" "+str(crop[1])+" "+str(crop[2])+" "+str(crop[3])+" "+temp1+" "+outtif+" -of GTiff"
            r = subprocess.call("gdalwarp -overwrite -s_srs EPSG:4326 -ts "+str(ncol)+" "+str(nlign)+" -r average \
                -te "+str(crop[0])+" "+str(crop[1])+" "+str(crop[2])+" "+str(crop[3])+" "+temp1+" "+outtif+" -of GTiff", shell=True)
            if r != 0:
                raise Exception("gdalwarp failed")

        elif proj and crop is False:
            print "gdalwarp -overwrite -s_srs EPSG:4326 -t_srs EPSG:"+str(EPSG)+" -ts "+str(ncol)+" "+str(nlign)+" -r average "+temp1+" "+outtif+" -of GTiff"
            r = subprocess.call("gdalwarp -overwrite -s_srs EPSG:4326 -t_srs EPSG:"+str(EPSG)+" -ts "+str(ncol)+" "+str(nlign)+" -r average "+temp1+" "+outtif+" -of GTiff", shell=True)
            if r != 0:
                raise Exception("gdalwarp failed")

        else:
            print "gdalwarp -overwrite -ts "+str(ncol)+" "+str(nlign)+" -r average "+temp1+" "+outtif+" -of GTiff"
            r = subprocess.call("gdalwarp -overwrite -ts "+str(ncol)+" "+str(nlign)+" -r average "+temp1+" "+outtif+" -of GTiff", shell=True)
            if r != 0:
                raise Exception("gdalwarp failed")

        # del temporary files
        cmd = 'rm -f {}'.format(temp1)
        os.system(cmd)

        # open gacos
        ds = gdal.Open(path+'{}_gacos.tif'.format(int(idates[i])), gdal.GA_ReadOnly)
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
    nfigure += 1
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

    loadf = 'cube_gacos'
    # save new cube
    fid = open(loadf, 'wb')
    gacos.flatten().astype('float32').tofile(fid)
    fid.close()

# # load gacos cube
gcubei = np.fromfile(loadf,dtype=np.float32)
gcube = as_strided(gcubei[:nlign*ncol*N])
# why did you mask gacos2data here? put gacos2data=1 if you already did the conversion before
gacos = gcube.reshape((nlign,ncol,N))*gacos2data

# Apply correction
maps_flat = np.zeros((nlign,ncol,N))
model_flat = np.zeros((nlign,ncol,N))

# initialise variance
var = np.zeros((N))

for l in xrange((N)):

  data = as_strided(maps[:,:,l]) 
  data_flat = as_strided(maps_flat[:,:,l])
  model = as_strided(gacos[:,:,l])
  
  if l == imref:
    # check if data and model are set to zeros
    print "REF DATE:", np.nanmean(data-model)
  
  else:

    _los_map = np.copy(data)
    _los_map[np.logical_or(data==0,data>=9990)] = np.float('NaN')
    losmin,losmax = np.nanpercentile(_los_map,2.),np.nanpercentile(_los_map,98.)
    gacosmin,gacosmax = np.nanpercentile(model,2),np.nanpercentile(model,98)

    if radar is not None:
        maxtopo,mintopo = np.nanpercentile(elev,98), np.nanpercentile(elev,2)
    else:
        maxtopo,mintopo = 1e8, -1e8

    # index = np.nonzero(data>2)
    # data[index] = np.float('NaN')
    # plt.imshow(data)
    # plt.show()
    # sys.exit()

    funct = 0.
    pix_lin, pix_col = np.indices((nlign,ncol))
   
    # find proportionality between data and model
    index = np.nonzero(
        np.logical_and(~np.isnan(data),
        np.logical_and(pix_lin>line_beg,
        np.logical_and(pix_lin<line_end,
        np.logical_and(data<losmax,
        np.logical_and(data>losmin,
        np.logical_and(rms<threshold_rms,
        np.logical_and(rms>1.e-6,
        np.logical_and(model<gacosmax,
        np.logical_and(model>gacosmin,
        np.logical_and(elev<maxtopo,
        np.logical_and(elev>mintopo,
        np.logical_and(data!=0.0,
        np.logical_and(model!=0.0,
        np.logical_and(pix_col>col_beg, pix_col<col_end
        )))))))))))))))

    indexref = np.nonzero(
        np.logical_and(~np.isnan(data),
        np.logical_and(pix_lin>refstart,
        np.logical_and(pix_lin<refend,
        np.logical_and(data<losmax,
        np.logical_and(data>losmin,
        np.logical_and(rms<threshold_rms,
        np.logical_and(rms>1.e-6,
        np.logical_and(model<gacosmax,
        np.logical_and(model>gacosmin,
        np.logical_and(elev<maxtopo,
        np.logical_and(elev>mintopo,
        np.logical_and(data!=0.0,
        np.logical_and(model!=0.0,
        np.logical_and(pix_col>col_beg, pix_col<col_end
        )))))))))))))))

    temp = np.array(index).T
    x = temp[:,0]; y = temp[:,1]
    los_clean = data[index].flatten()
    model_clean = model[index].flatten()
   
    # compute average phase in the ref area 
    # we want los ref area to be zero
    los_ref = data[indexref].flatten() - model[indexref].flatten()
    rms_ref = rms[indexref].flatten()
    amp_ref = 1./rms_ref
    amp_ref = amp_ref/np.nanpercentile(amp_ref,99)
    print 'Ref area set to zero:', refstart,refend
    # weigth average of the phase
    rg_ref, az_ref,model_ref = np.nanmean(pix_col[indexref]),np.nanmean(pix_lin[indexref]),np.nanmean(model[indexref])
    cst = np.nansum(los_ref*amp_ref) / np.nansum(amp_ref)
    print 'Average phase within ref area:', cst
    logger.info('Average rg: {0}, az:{1}, topo:{2}, within ref area'.format(np.int(rg_ref), np.int(az_ref), np.int(model_ref)))

    # digitize data in bins, compute median and std
    bins = np.arange(gacosmin,gacosmax,abs(gacosmax-gacosmin)/500.)
    inds = np.digitize(model_clean,bins)
    modelbins = []
    losbins = []
    losstd = []
    xbins, ybins = [], []
    for j in range(len(bins)-1):
            uu = np.flatnonzero(inds == j)
            if len(uu)>500:
                modelbins.append(bins[j] + (bins[j+1] - bins[j])/2.)

                indice = np.flatnonzero(np.logical_and(los_clean[uu]>np.percentile(\
                    los_clean[uu],10.),los_clean[uu]<np.percentile(los_clean[uu],90.)))

                losstd.append(np.std(los_clean[uu][indice]))
                losbins.append(np.median(los_clean[uu][indice]))
                xbins.append(np.median(x[uu][indice]))
                ybins.append(np.median(y[uu][indice]))

    losbins = np.array(losbins)
    losstd = np.array(losstd)
    modelbins = np.array(modelbins)
    xbins, ybins = np.array(xbins),np.array(ybins)

    if (ramp == 'cst' and  fitmodel=='no'):
        
        a = cst
        print 'Remove cst: %f for date: %i'%(a,idates[l])

        remove_ramp = np.ones((nlign,ncol)).reshape((nlign,ncol))*cst
        remove_ramp[model==0.] = 0.
        remove_ramp[np.isnan(data)] = np.float('NaN')
        
        # data = gacos + cst
        remove = remove_ramp + model
       
        # Compute ramp for data = f(model) 
        functbins = a
        funct = a
        # set coef gacos to 1
        f = 1

    elif (ramp == 'lin' and  fitmodel=='no'):
        # here we want to minimize data-gacos = ramp
        # invers both digitized data and ref frame together
        d = np.hstack([losbins-modelbins,cst])
        # give a strong weigth to ref frame
        sigmad = np.hstack([losstd,1e-3])       

        G=np.zeros((len(d),3))
        G[:-1,0] = xbins
        G[:-1,1] = ybins
        G[:,2] = 1
        G[-1,0] = az_ref
        G[-1,1] = rg_ref

        x0 = lst.lstsq(G,d)[0]
        # print x0
        _func = lambda x: np.sum(((np.dot(G,x)-d)/sigmad)**2)
        _fprime = lambda x: 2*np.dot(G.T/sigmad, (np.dot(G,x)-d)/sigmad)
        pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
        a = pars[0]; b = pars[1]; c = pars[2]
        print 'Remove ramp  %f az  + %f r + %f for date: %i'%(a,b,c,idates[l])
            
        # Compute ramp for data = f(model)
        functbins = a*xbins + b*ybins + c
        funct = a*x + b*y + c
        # set coef gacos to 1
        f = 1

        # build total G matrix
        G=np.zeros((len(data.flatten()),3))
        for i in xrange(nlign):
            G[i*ncol:(i+1)*ncol,0] = i - line_beg 
            G[i*ncol:(i+1)*ncol,1] = np.arange(ncol) - col_beg
        G[:,2] = 1

        # compute ramp
        remove_ramp = np.dot(G,pars).reshape(nlign,ncol)
        remove_ramp[model==0.] = 0.
        remove_ramp[np.isnan(data)] = np.float('NaN')
        
        # data = gacos + (a*rg + b*az + cst) 
        remove = model + remove_ramp

    elif fitmodel=='yes':
        # invers both digitized data and ref frame together
        # here we invert data = a*gacos + ramp
        d = np.hstack([losbins,cst])
        # give a strong weigth to ref frame
        sigmad = np.hstack([losstd,1e-3])       
            
        G=np.zeros((len(d),6))
        G[:-1,0] = xbins**2
        G[:-1,1] = xbins
        G[:-1,2] = ybins**2
        G[:-1,3] = ybins
        G[:,4] = 1
        G[:-1,5] = modelbins
        G[-1,0] = az_ref**2
        G[-1,1] = az_ref
        G[-1,2] = rg_ref**2
        G[-1,3] = rg_ref
        G[-1,5] = modelref

        x0 = lst.lstsq(G,d)[0]
        # print x0
        _func = lambda x: np.sum(((np.dot(G,x)-d)/sigmad)**2)
        _fprime = lambda x: 2*np.dot(G.T/sigmad, (np.dot(G,x)-d)/sigmad)
        pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
        a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
        print 'Remove ramp %f az**2, %f az  + %f r**2 + %f r + %f + %f model for date: %i'%(a,b,c,d,e,f,idates[l])

        #Compute ramp for data = f(model)
        funct = a*x**2 + b*x + c*y**2 + d*y + e
        functbins = a*xbins**2 + b*xbins + c*ybins**2 + d*ybins + e

        # build total G matrix
        G=np.zeros((len(data.flatten()),6))
        for i in xrange(nlign):
                G[i*ncol:(i+1)*ncol,0] = (i - line_beg)**2
                G[i*ncol:(i+1)*ncol,1] = i - line_beg
                G[i*ncol:(i+1)*ncol,2] = np.arange((ncol - col_beg))**2
                G[i*ncol:(i+1)*ncol,3] = np.arange(ncol - col_beg)
        G[:,4] = 1
        G[:,5] = model.flatten()

        # data = a*gacos + ramp
        remove = np.dot(G,pars).reshape(nlign,ncol)
        remove[model==0.] = 0.
        remove[np.isnan(data)] = np.float('NaN')

        # build total G matrix
        G=np.zeros((len(data.flatten()),5))
        for i in xrange(nlign):
            G[i*ncol:(i+1)*ncol,0] = (i - line_beg)**2
            G[i*ncol:(i+1)*ncol,1] = i - line_beg
            G[i*ncol:(i+1)*ncol,2] = np.arange((ncol - col_beg))**2
            G[i*ncol:(i+1)*ncol,3] = np.arange(ncol - col_beg)
        G[:,4] = 1

        # ramp only
        remove_ramp = np.dot(G,pars[:-1]).reshape(nlign,ncol)
        remove_ramp[model==0.] = 0.
        remove_ramp[np.isnan(data)] = np.float('NaN')
    
    # correction
    # data_flat = data - (gacos + ramp)
    data_flat[:,:] = data - remove
    data_flat[np.isnan(data)] = np.float('NaN')
    data_flat[data_flat>999.]= np.float('NaN')
    
    # model = gacos + ramp
    model_flat[:,:,l] = gacos[:,:,l] + remove_ramp
 
    # Refer data again (just to check)
    rms_ref = rms[indexref].flatten()
    amp_ref = 1./rms_ref
    amp_ref = amp_ref/np.nanmax(amp_ref)
    cst = np.nansum(data_flat[indexref]*amp_ref) / np.nansum(amp_ref)
    print 'Average phase within ref area, iter=2:', cst 
    data_flat = data_flat - cst

    # compute variance flatten data
    var[l] = np.sqrt(np.nanmean(data_flat**2))
    print 'Var: ', var[l]

    if plot == 'yes':
        # initiate figure depl
        fig = plt.figure(nfigure,figsize=(14,7))
        nfigure += 1
        ax = fig.add_subplot(3,2,1)
        im = ax.imshow(data,cmap=cmap,vmax=losmax,vmin=losmin)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
        ax.set_title('Data {}'.format(idates[l]),fontsize=6)

        # initiate figure depl
        ax = fig.add_subplot(3,2,2)
        im = ax.imshow(model,cmap=cmap,vmax=losmax,vmin=losmin)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
        ax.set_title('Model {}'.format(idates[l]),fontsize=6)
        
        # initiate figure depl
        ax = fig.add_subplot(3,2,3)
        im = ax.imshow(model_flat[:,:,l],cmap=cmap)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
        ax.set_title('Flatten Model {}'.format(idates[l]),fontsize=6)

        # initiate figure depl
        ax = fig.add_subplot(3,2,4)
        im = ax.imshow(data_flat,cmap=cmap,vmax=losmax,vmin=losmin)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
        ax.set_title('Correct Data {}'.format(idates[l]),fontsize=6)

        ax = fig.add_subplot(3,2,5)
        g = np.linspace(np.nanmax(model_clean),np.nanmin(model_clean),100)
        ax.scatter(model_clean,los_clean - funct, s=0.005, alpha=0.1, rasterized=True)
        ax.plot(modelbins,losbins - functbins,'-r', lw =.5)
        ax.plot(g,f*g,'-r', lw =4.)
        ax.set_ylim([(los_clean-funct).min(),(los_clean-funct).max()])
        ax.set_xlim([model_clean.min(),model_clean.max()])
        ax.set_xlabel('GACOS ZTD')
        ax.set_ylabel('LOS delay')
        ax.set_title('Data/Model')

        fig.tight_layout()
        fig.savefig('{}-gacos-cor.eps'.format(idates[l]), format='EPS',dpi=150)
        plt.show()
        # sys.exit()

# save new cube
fid = open('depl_cumule_gacos', 'wb')
maps_flat.flatten().astype('float32').tofile(fid)
fid.close()
# save gacos ref spatialy
fid = open('gacos_ref', 'wb')
model_flat.flatten().astype('float32').tofile(fid)
fid.close()

# save rms
np.savetxt('rms_gacos.txt', var, header='# date   |   RMS', fmt=('%.8f'))
fig_rms = plt.figure(100,figsize=(14,7))
plt.plot(dates, var)
fig_rms.savefig('rms.eps', format='EPS',dpi=150)

