#!/usr/bin/env python2
# -*- coding: utf-8 -*-
############################################


############################################
# Author        : Simon DAOUT (Oxford)
############################################

"""\
invert_ramp_topo_unw.py
-------------
Removes atmospheric phase/elevation correlations or/and azimuthal and range ramps polynomial coeffice
ints on unwrapped interferograms (2 bands file). Reconstruction of the empirical phase correction by time series inversion.

usage: invert_ramp_topo_unw.py --int_list=<path> [--int_path=<path>] [--out_path=<path>] \
[--prefix=<value>] [--suffix=<value>] [--rlook=<value>]  [--ref=<path>] [--format=<value>] \
[--flat=<0/1/2/3/4/5/6>] [--topofile=<path>] [--ivar=<0/1>] [--nfit=<0/1>] [--tsinv=<yes/no>]\
[--estim=yes/no] [--mask=<path>] [--threshold_mask=<value>]  \
[--cohpixel=<yes/no>] [--threshold_coh=<value>] \
[--ibeg_mask=<value>] [--iend_mask=<value>] [--perc=<value>] \
[--plot=<yes/no>] [--suffix_output=<value>]\
[<ibeg>] [<iend>] [<jbeg>] [<jend>] 

--int_list PATH       Text file containing list of interferograms dates in two colums, $data1 $date2
--int_path PATh       Relative path to input interferograms directory
--int_path PATh       Relative path to output interferograms directory
--prefix VALUE        Prefix name $prefix$date1-$date2$suffix_$rlookrlks.unw [default: '']
--suffix value        Suffix name $prefix$date1-$date2$suffix_$rlookrlks.unw [default: '']
--rlook value         look int. $prefix$date1-$date2$suffix_$rlookrlks.unw [default: 0]
--ref PATH            Path to reference image to define format and size. Necessary if topofile is None [default:None]     
--flat PATH           Remove a spatial ramp.If short acquisition, short is automatically set to 3.
0: ref frame [default], 1: range ramp ax+b , 2: azimutal ramp ay+b, 
3: ax+by+c, 4: ax+by+cxy+d 5: ax**2+bx+d, 6: ay**2+by+c
--topofile  PATH      Path to the radar_look.hgt file. If not None add phase/elevation relationship in the relation [default:None]
--nfit VALUE          fit degree in azimuth or in elevation
--ivar VALUE          define phase/elevation relationship: ivar=0 function of elevation, ivar=1 crossed function of azimuth and elevation

if ivar=0 and nfit=0, add linear elev. term (z) to ramps estimation defined by the flat argument such as:
0: ref frame+ez [default], 1: range ramp ax+b+ez , 2: azimutal ramp ay+b+ez, 
3: ax+by+c+ez, 4: ax+by+cxy+d+ez, 5: ax**2+bx+d+ez, 6: ay**2+by+c+ez

if ivar=0 and nfit=1, add quadratic elev. term (z) to ramps estimation defined by the flat argument such as:
0: ref frame+ez+fz**2 [default], 1: range ramp ax+b+ez+fz**2 , 2: azimutal ramp ay+b+ez+fz**2, 
3: ax+by+c+ez+fz**2, 4: ax+by+cxy+d+ez+fz**2, 5: ax**2+bx+d+ez+fz**2, 6: ay**2+by+c+ez+fz**2

if ivar=1 and nfit=0, add cross function of elev. (z) and azimuth to ramps estimation defined by the flat argument such as:
0: ref frame+ez+fz*az [default], 1: range ramp ax+b+ez+fz*az , 2: azimutal ramp ay+b+ez+fz*az, 
3: ax+by+c+ez+fz*az, 4: ax+by+cxy+d+ez+fz*az, 5: ax**2+bx+d+ez+fz*az, 6: ay**2+by+c+ez+fz*az

if ivar=1 and nfit=1, add quadratic cross function of elev. (z) and azimuth to ramps estimation defined by the flat argument such as:
0: ref frame+ez+fz*az+g(z*az)**2 [default], 1: range ramp ax+b+ez+fz*az+g*(z*az)**2 , 2: azimutal ramp ay+b+ez+fz*az+g*(z*az)**2, 
3: ax+by+c+ez+fz*az+g*(z*az)**2, 4: ax+by+cxy+d+ez+fz*az+g*(z*az)**2, 5: ax**2+bx+d+ez+fz*az+g*(z*az)**2, 6: ay**2+by+c+ez+fz*az+g*(z*az)**2

--tsinv yes/no        If yes, invert corrected phase into time series [default:no]
--estim yes/no        If yes, do the estimation, otherwise read input files corection_matrix and liste_coeff_ramps.txt [default:yes]
--mask PATH           Mask in .r4 format. Keep only values > threshold_mask. [default:None]
--threshold_mask      Threshold on mask: take only values > threshold_mask [default: -1]
--cohpixel  yes/no    If Yes, use amplitude interferogram to weight and mask pixels (e.g Coherence, Colinearity, Amp Filter) [default: no]
--format VALUE        Format input files: ROI_PAC, GAMMA, GTIFF [default: ROI_PAC]
--threshold_coh VALUE Threshold on cohpixel file [default:0]
--ibeg_mask VALUE     Start line number mask [default: None]
--iend_mask VALUE     Strop line number mask [default: None]  
--perc VALUE          Percentile of hidden LOS pixel for the estimation and clean outliers [default:98.]
--plot yes/no         If yes, plot figures for each ints [default: no]
--suffix_output value Suffix output file name $prefix$date1-$date2$suffix$suffix_output [default:_corrunw]
--ibeg VALUE          Line number bounding the estimation zone [default: 0]
--iend VALUE          Line number bounding the estimation zone [default: mlines]
--jbeg VALUE          Column numbers bounding the estimation zone [default: 0]
--jend VALUE          Column numbers bounding the estimation zone [default: mcols] 
"""


# gdal
import gdal
gdal.UseExceptions()

# system
from os import path, environ
import os

# plot
import matplotlib
if environ["TERM"].startswith("screen"):
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from datetime import datetime
import time

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided

# scipy
import scipy
import scipy.optimize as opt
import scipy.linalg as lst

import docopt
import gamma as gm
import shutil

# import shutil

# A more predictable makedirs
def makedirs(name):
    if os.path.exists(name):
        return
    os.makedirs(name)

def date2dec(dates):
    dates  = np.atleast_1d(dates)
    times = []
    for date in dates:
        x = datetime.strptime('{}'.format(date),'%Y%m%d')
        #x = datetime.strptime('{}'.format(date),'%Y-%m-%dT%H:%M:%S.%fZ')
        dec = float(x.strftime('%j'))/365.1
        year = float(x.strftime('%Y'))
        # print date,dec,year
        times.append(year + dec)
    return times


# read arguments
arguments = docopt.docopt(__doc__)

int_list=arguments["--int_list"]

if arguments["--int_path"] == None:
    int_path='./'
else:
    int_path=arguments["--int_path"] + '/'

if arguments["--out_path"] == None:
    out_path=np.copy(int_path)
else:
    out_path=arguments["--out_path"] + '/'

makedirs(out_path)

if arguments["--prefix"] == None:
    prefix = ''
else:
    prefix=arguments["--prefix"]
if arguments["--suffix"] == None:
    suffix = ''
else:
    suffix=arguments["--suffix"]
if arguments["--rlook"] == None:
    rlook = ''
else:
    rlook = '_' + arguments["--rlook"] + 'rlks'

if arguments["--flat"] == None:
    flat = 0
elif int(arguments["--flat"]) <  7:
    flat = int(arguments["--flat"])
else:
    print 'Carefull: flat > 6, set flat to 0'
    flat = 0

if arguments["--ivar"] == None:
    ivar = 0
elif int(arguments["--ivar"]) <  2:
    ivar = int(arguments["--ivar"])
else:
    print 'Carefull: ivar > 1, set ivar to 0'
    ivar = 0

if arguments["--nfit"] == None:
    nfit = 0
elif int(arguments["--nfit"]) <  2:
    nfit = int(arguments["--nfit"])
else:
    print 'Error: nfit > 1, set nfit to 0'
    nfit = 0

if arguments["--mask"] ==  None or not os.path.exists(arguments["--mask"]):
    maskfile = None
else:
    maskfile = arguments["--mask"]
if arguments["--threshold_mask"] ==  None:
    threshold_mask = -1
else:
    threshold_mask = float(arguments["--threshold_mask"])
if arguments["--cohpixel"] ==  None:
    rmsf = 'no'
else:
    rmsf = arguments["--cohpixel"]
if arguments["--threshold_coh"] ==  None:
    threshold_rms = 0.
else:
    threshold_rms = float(arguments["--threshold_coh"])
if arguments["--topofile"] ==  None or not os.path.exists(arguments["--topofile"]):
   radar = None
else:
   radar = arguments["--topofile"]

if arguments["--ref"] ==  None :
   ref = None
else:
   ref = arguments["--ref"]
if ref == None and radar == None:
    print 'Argument error: Need to give ref or topographic file'
    sys.exit()

if arguments["--tsinv"] ==  None:
    tsinv = 'no'
else:
    tsinv = arguments["--tsinv"]

if arguments["--format"] ==  None:
    sformat = 'ROI_PAC'
else:
    sformat = arguments["--format"]

if arguments["--estim"] ==  None:
    estim = 'yes'
else:
    estim = arguments["--estim"]

if arguments["--ibeg_mask"] ==  None:
    ibeg_mask = np.inf
else:
    ibeg_mask = int(arguments["--ibeg_mask"])
if arguments["--iend_mask"] ==  None:
    iend_mask = -np.inf
else:
    iend_mask = int(arguments["--iend_mask"])
if arguments["--perc"] ==  None:
    perc = 98.
else:
    perc = float(arguments["--perc"])
if arguments["--plot"] ==  None:
    plot = 'no'
else:
    plot = str(arguments["--plot"])
if arguments["--suffix_output"] ==  None:
    suffout = '_corrunw'
else:
    suffout = arguments["--suffix_output"]


#####################################################################################

print
# read int
date_1,date_2=np.loadtxt(int_list,comments="#",unpack=True,dtype='i,i')
kmax=len(date_1)
print "number of interferogram: ",kmax

# list dates
im = []; bt = []
for date1,date2 in zip(date_1,date_2):
    if date1 not in im: im.append(date1)
    if date2 not in im: im.append(date2)
nmax=len(im)
print "number of image: ",nmax
imd = date2dec(im)
cst = np.copy(imd[0])
for i in xrange((nmax)):
    bt.append(imd[i]-cst)

# # Now, write list_pair
# wf = open("list_dates", "w")
# for i in xrange((nmax)):
#     wf.write("%i %.6f %.6f\n" % (im[i], imd[i], bt[i]))
# wf.close()
# sys.exit()

# initialise number of figure
nfigure=0
# fig, ax = plt.subplots(1)
# x = [date2num(datetime.strptime('{}'.format(d),'%Y%m%d')) for d in im]
# ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
# ax.plot(x,bp,"ro",label='acquistions')
# ax.plot(x,,"green")
# fig.autofmt_xdate()
# ax.set_xlabel('Time (Year/month/day)')
# ax.set_ylabel('Perpendicular Baseline')
# plt.legend(loc=2)
# plt.show()

# load ref to define mlines, mcols and format
if ref is not None:
    ds_extension = os.path.splitext(ref)[1]
    if sformat == 'GTIFF':
        geotiff = arguments["--geotiff"]
        georef = gdal.Open(ref)
        gt = georef.GetGeoTransform()
        proj = georef.GetProjection()
        driver = gdal.GetDriverByName('GTiff')
        ds = gdal.Open(ref, gdal.GA_ReadOnly)
        mlines,mcols = ds.RasterYSize, ds.RasterXSize
    elif sformat == 'ROI_PAC':
        driver = gdal.GetDriverByName("roi_pac")
        ds = gdal.Open(ref, gdal.GA_ReadOnly)
        mlines,mcols = ds.RasterYSize, ds.RasterXSize
    elif sformat == 'GAMMA':
        # par_file = ref 
        mlines,mcols = gm.readpar(int_path)

# laod elevation map
if radar is not None:

    if sformat == 'GTIFF':
        geotiff = arguments["--geotiff"]
        georef = gdal.Open(radar)
        gt = georef.GetGeoTransform()
        proj = georef.GetProjection()
        driver = gdal.GetDriverByName('GTiff')
        ds = gdal.Open(radar, gdal.GA_ReadOnly)
        ds_band2 = ds.GetRasterBand(2)
        mlines,mcols = ds.RasterYSize, ds.RasterXSize
        elev_map = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        del ds

    elif sformat == 'ROI_PAC':
        driver = gdal.GetDriverByName("roi_pac")
        ds = gdal.Open(radar, gdal.GA_ReadOnly)
        ds_band2 = ds.GetRasterBand(2)
        mlines,mcols = ds.RasterYSize, ds.RasterXSize
        elev_map = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        del ds
    
    elif sformat == 'GAMMA':
        # par_file = ref 
        mlines,mcols = gm.readpar()
        elev_map = gm.readgamma(radar)

    maxelev,minelev = np.nanpercentile(elev_map,99),np.nanpercentile(elev_map,1)
    
else:
    maxelev,minelev = 1.,-1
    elev_map = np.zeros((mlines,mcols))


# open mask file
if maskfile is not None:
    fid = open(maskfile,'r')
    maski = np.fromfile(fid,dtype=np.float32)[:mcols*mlines]
    mask = maski.reshape((mlines,mcols))
    k = np.nonzero(mask<threshold_mask)
    spacial_mask = np.copy(mask)
    spacial_mask[k] = float('NaN')

    if plot=='yes':

      nfigure=nfigure+1
      fig = plt.figure(nfigure,figsize=(5,4))
      ax = fig.add_subplot(1,1,1)
      cax = ax.imshow(spacial_mask,cmap=cm.jet)
      ax.set_title('Mask')
      setp( ax.get_xticklabels(), visible=None)
      fig.colorbar(cax, orientation='vertical',aspect=10)
      plt.show()
      fid.close()

else:
    mask = np.zeros((mlines,mcols))
    threshold_mask = -1

# define estimation area
if arguments["<ibeg>"] ==  None:
  ibeg = 0
else:
  ibeg = int(arguments["<ibeg>"])
if arguments["<iend>"] ==  None:
  iend = mlines
else:
  iend = int(arguments["<iend>"])
if arguments["<jbeg>"] ==  None:
  jbeg = 0
else:
  jbeg = int(arguments["<jbeg>"])
if arguments["<jend>"] ==  None:
  jend = mcols
else:
  jend = int(arguments["<jend>"])

# extract range and azimuth coordinates from ref or radar file
pix_az, pix_rg = np.indices((mlines,mcols))
# print np.shape(pix_az)
# print np.shape(mask)
# sys.exit()

def estim_ramp(los,los_clean,topo_clean,x,y,order,rms,nfit,ivar):


    # initialise full vector 
    sol = np.zeros((13))
    #0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2

    # initialize correction
    corr = np.zeros((mlines,mcols))

    if order==0:  
    #0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2

        if radar is None:
            G=np.zeros((len(los_clean),1))
            G[:,0] = 1

            # ramp inversion
            x0 = lst.lstsq(G,los_clean)[0]
            try:
                _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=50,full_output=True,iprint=0)[0]
            except:
                pars = x0
            sol[8] = pars[0]
            print 'Remove ref frame %f'%(pars[0])

            # build total G matrix
            G=np.zeros((len(los),1))
            G[:,0] = 1

        else:
            if ivar==0 & nfit==0:
                G=np.zeros((len(los_clean),2))
                G[:,0] = 1
                G[:,1] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=50,full_output=True,iprint=0)[0]
                except:
                    pars = x0

                sol[8] = pars[0]; sol[9] = pars[1]
                print 'Remove ref frame %f + %f z'%(pars[0],pars[1])

                # build total G matrix
                G=np.zeros((len(los),2))
                G[:,0] = 1
                G[:,1] = elev_map.flatten()

            elif ivar==0 & nfit==1:
                G=np.zeros((len(los_clean),3))
                G[:,0] = 1
                G[:,1] = topo_clean
                G[:,2] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                
                sol[8] = pars[0]; sol[9] = pars[1]; sol[10] = pars[2]
                print 'Remove ref frame %f + %f z + %f z**2'%(pars[0],pars[1],pars[2])

                # build total G matrix
                G=np.zeros((len(los),3))
                G[:,0] = 1
                G[:,1] = elev_map.flatten()
                G[:,2] = elev_map.flatten()**2
            
            elif ivar==1 and nfit==0:
                G=np.zeros((len(los_clean),3))
                G[:,0] = 1
                G[:,1] = topo_clean
                G[:,2] = x*topo_clean

                # ramp inversion 
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0

                sol[8] = pars[0]; sol[9] = pars[1]; sol[11] = pars[2]
                print 'Remove ref frame %f + %f z + %f az*z'%(pars[0],pars[1],pars[2])

                # build total G matrix
                G=np.zeros((len(los),3))
                G[:,0] = 1
                G[:,1] = elev_map.flatten()
                G[:,2] = elev_map.flatten()
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,2] *= i

            elif ivar==1 and nfit==1:
                G=np.zeros((len(los_clean),4))
                G[:,0] = 1
                G[:,1] = topo_clean
                G[:,2] = x*topo_clean
                G[:,3] = (x*topo_clean)**2

                # ramp inversion 
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[8] = pars[0]; sol[9] = pars[1]; sol[11] = pars[2];sol[12] = pars[3]
                print 'Remove ref frame %f + %f z + %f az*z + %f (az*z)**2'%(pars[0],pars[1],pars[2],pars[3])

                # build total G matrix
                G=np.zeros((len(los),4))
                G[:,0] = 1
                G[:,1] = elev_map.flatten()
                G[:,2] = elev_map.flatten()
                G[:,3] = elev_map.flatten()**2
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,2] *= i
                    G[i*mcols:(i+1)*mcols,3] *= i**2


    elif order==1: # Remove a range ramp ay+b for each maps (y = col)

        if radar is None:
            G=np.zeros((len(los_clean),2))
            G[:,0] = y
            G[:,1] = 1

            # ramp inversion
            x0 = lst.lstsq(G,los_clean)[0]
            try:
                _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
            except:
                pars = x0
            sol[2] = pars[0]; sol[8] = pars[1]
            print 'Remove ramp %f r + %f'%(pars[0],pars[1])

            # build total G matrix
            G=np.zeros((len(los),2))
            for i in xrange(mlines):
                G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) 
            G[:,1] = 1


        else:
            if ivar==0 & nfit==0:
                G=np.zeros((len(los_clean),3))
                G[:,0] = y
                G[:,1] = 1
                G[:,2] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[8] = pars[1]; sol[9] = pars[2]
                print 'Remove ramp %f r + %f + %f z '%(pars[0],pars[1],pars[2])

                # build total G matrix
                G=np.zeros((len(los),3))
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) 
                G[:,1] = 1
                G[:,2] = elev_map.flatten()

            if ivar==0 & nfit==1:
                G=np.zeros((len(los_clean),4))
                G[:,0] = y
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[8] = pars[1]; sol[9] = pars[2]; sol[10] = pars[3]
                print 'Remove ramp %f r + %f + %f z + %f z**2'%(pars[0],pars[1],pars[2],pars[3])

                # build total G matrix
                G=np.zeros((len(los),4))
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) 
                G[:,1] = 1
                G[:,2] = elev_map.flatten()
                G[:,3] = elev_map.flatten()**2
            
            elif ivar==1 and nfit==0:
                G=np.zeros((len(los_clean),4))
                G[:,0] = y
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean*x

                # ramp inversion
                #y**3 y**2 y x**3 x**2 x xy**2 xy cst z z**2 yz yz**2 
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[8] = pars[1]; sol[9] = pars[2]; sol[11] = pars[3]
                print 'Remove ramp %f r + %f + %f z + %f z*az'%(pars[0],pars[1],pars[2],pars[3])

                # build total G matrix
                G=np.zeros((len(los),4))
                G[:,1] = 1
                G[:,2] = elev_map.flatten()
                G[:,3] = elev_map.flatten()
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) 
                    G[i*mcols:(i+1)*mcols,3] *= i

            elif ivar==1 and nfit==1:
                G=np.zeros((len(los_clean),4))
                G[:,0] = y
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean*x
                G[:,4] = (topo_clean*x)**2

                # ramp inversion
                #y**3 y**2 y x**3 x**2 x xy**2 xy cst z z**2 yz yz**2 
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[8] = pars[1]; sol[9] = pars[2]; sol[11] = pars[3]; sol[12] = pars[4]
                print 'Remove ramp %f r + %f + %f z + %f z*az'%(pars[0],pars[1],pars[2],pars[3],pars[4])

                # build total G matrix
                G=np.zeros((len(los),5))
                G[:,1] = 1
                G[:,2] = elev_map.flatten()
                G[:,3] = elev_map.flatten()
                G[:,4] = elev_map.flatten()**2
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) 
                    G[i*mcols:(i+1)*mcols,3] *= i
                    G[i*mcols:(i+1)*mcols,4] *= i**2

        
    elif order==2: # Remove an azimutal ramp ax+b for each maps (x is lign)
    #y**3 y**2 y x**3 x**2 x xy**2 xy cst z z**2 z*az az*z**2

        if radar is None:
            G=np.zeros((len(los_clean),2))
            G[:,0] = x
            G[:,1] = 1

            # ramp inversion
            x0 = lst.lstsq(G,los_clean)[0]
            try:
                _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
            except:
                pars = x0
            sol[5] = pars[0]; sol[8] = pars[1]
            print 'Remove ramp %f az + %f'%(pars[0],pars[1])

            # build total G matrix
            G=np.zeros((len(los),2))
            for i in xrange(mlines):
                G[i*mcols:(i+1)*mcols,0] = i  
            G[:,1] = 1

        else:
            if ivar==0 & nfit==0:
                G=np.zeros((len(los_clean),3))
                G[:,0] = x
                G[:,1] = 1
                G[:,2] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[5] = pars[0]; sol[8] = pars[1]; sol[9] = pars[2]
                print 'Remove ramp %f az + %f + %f z'%(pars[0],pars[1],pars[2])

                # build total G matrix
                G=np.zeros((len(los),3))
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = i 
                G[:,1] = 1
                G[:,2] = elev_map.flatten()

            if ivar==0 & nfit==1:
                G=np.zeros((len(los_clean),4))
                G[:,0] = x
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[5] = pars[0]; sol[8] = pars[1]; sol[9] = pars[2]; sol[10] = pars[-1]
                print 'Remove ramp %f az + %f + %f z + %f z**2'%(pars[0],pars[1],pars[2],pars[3])

                # build total G matrix
                G=np.zeros((len(los),4))
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = i 
                G[:,1] = 1
                G[:,2] = elev_map.flatten()
                G[:,3] = elev_map.flatten()**2
            
            elif ivar==1 and nfit==0:
                G=np.zeros((len(los_clean),4))
                G[:,0] = x
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean*x

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[5] = pars[0]; sol[8] = pars[1]; sol[9] = pars[2]; sol[11] = pars[3]
                print 'Remove ramp %f az + %f + %f z + %f z*az'%(pars[0],pars[1],pars[2],pars[3])

                # build total G matrix
                G=np.zeros((len(los),4))
                G[:,1] = 1
                G[:,2] = elev_map.flatten()
                G[:,3] = elev_map.flatten()
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = i
                    G[i*mcols:(i+1)*mcols,3] *= i

            elif ivar==1 and nfit==1:
                G=np.zeros((len(los_clean),5))
                G[:,0] = x
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean*x
                G[:,4] = (topo_clean*x)**2

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[5] = pars[0]; sol[8] = pars[1]; sol[9] = pars[2]; sol[11] = pars[3]; sol[12] = pars[4]
                print 'Remove ramp %f az + %f + %f z + %f z*az + %f (z*az)**2'%(pars[0],pars[1],pars[2],pars[3],pars[4])

                # build total G matrix
                G=np.zeros((len(los),5))
                G[:,1] = 1
                G[:,2] = elev_map.flatten()
                G[:,3] = elev_map.flatten()
                G[:,3] = elev_map.flatten()**2
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = i
                    G[i*mcols:(i+1)*mcols,3] *= i
                    G[i*mcols:(i+1)*mcols,4] *= i**2

    elif order==3: # Remove a ramp ay+bx+c for each maps
    #y**3 y**2 y x**3 x**2 x xy**2 xy cst z z**2 z*az az*z**2

        if radar is None:
            G=np.zeros((len(los_clean),3))
            G[:,0] = y
            G[:,1] = x
            G[:,2] = 1

            # ramp inversion
            x0 = lst.lstsq(G,los_clean)[0]
            try:
                _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
            except:
                pars = x0
            sol[2] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]
            print 'Remove ramp %f r  + %f az + %f'%(pars[0],pars[1],pars[2])

            # build total G matrix
            G=np.zeros((len(los),3))
            for i in xrange(mlines):
                G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) 
                G[i*mcols:(i+1)*mcols,1] = i 
            G[:,2] = 1


        else:
            if ivar==0 and nfit==0:
                G=np.zeros((len(los_clean),4))
                G[:,0] = y
                G[:,1] = x
                G[:,2] = 1
                G[:,3] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:    
                    pars = x0
                sol[2] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]
                print 'Remove ramp %f r  + %f az + %f + %f z '%(pars[0],pars[1],pars[2],pars[3])

                # build total G matrix
                G=np.zeros((len(los),4))
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) 
                    G[i*mcols:(i+1)*mcols,1] = i    
                G[:,2] = 1
                G[:,3] = elev_map.flatten()

            if ivar==0 and nfit==1:
                G=np.zeros((len(los_clean),5))
                G[:,0] = y
                G[:,1] = x
                G[:,2] = 1
                G[:,3] = topo_clean
                G[:,4] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]; sol[9] = pars[-2]; sol[10] = pars[-1]
                print 'Remove ramp %f r  + %f az + %f + %f z + %f z**2'%(pars[0],pars[1],pars[2],pars[3],pars[4])

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) 
                    G[i*mcols:(i+1)*mcols,1] = i    
                G[:,2] = 1
                G[:,3] = elev_map.flatten()
                G[:,4] = elev_map.flatten()**2
            
            elif ivar==1 and nfit==0:
                G=np.zeros((len(los_clean),5))
                G[:,0] = y
                G[:,1] = x
                G[:,2] = 1
                G[:,3] = topo_clean
                G[:,4] = topo_clean*x

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]; sol[11] = pars[4]
                print 'Remove ramp %f r  + %f az + %f + %f z + %f z*az'%(pars[0],pars[1],pars[2],pars[3],pars[4])

                # build total G matrix
                G=np.zeros((len(los),5))
                G[:,2] = 1
                G[:,3] = elev_map.flatten()
                G[:,4] = elev_map.flatten()
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) 
                    G[i*mcols:(i+1)*mcols,1] = i
                    G[i*mcols:(i+1)*mcols,4] *= i   
            
            elif ivar==1 and nfit==1:
                G=np.zeros((len(los_clean),6))
                G[:,0] = y
                G[:,1] = x
                G[:,2] = 1
                G[:,3] = topo_clean
                G[:,4] = topo_clean*x
                G[:,5] = (topo_clean*x)**2

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                #0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2    
                sol[2] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]; sol[11] = pars[4]; sol[12] = pars[5]
                print 'Remove ramp %f r  + %f az + %f + %f z + %f z*az + %f (z*az)**2'%(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5])

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,2] = 1
                G[:,3] = elev_map.flatten()
                G[:,4] = elev_map.flatten()
                G[:,5] = elev_map.flatten()**2
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols)) 
                    G[i*mcols:(i+1)*mcols,1] = i
                    G[i*mcols:(i+1)*mcols,4] *= i   
                    G[i*mcols:(i+1)*mcols,5] *= i**2


    elif order==4:
    #y**3 y**2 y x**3 x**2 x xy**2 xy cst z z**2 z*az az*z**2

        if radar is None:
            G=np.zeros((len(los_clean),4))
            G[:,0] = y
            G[:,1] = x
            G[:,2] = y*x
            G[:,3] = 1

            # ramp inversion
            x0 = lst.lstsq(G,los_clean)[0]
            try:
                _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
            except:
                pars = x0
            sol[2] = pars[0]; sol[5] = pars[1]; sol[7] = pars[2]; sol[8] = pars[3]
            print 'Remove ramp %f r %f az  + %f r*az + %f'%(pars[0],pars[1],pars[2],pars[3])

            # build total G matrix
            G=np.zeros((len(los),4))
            for i in xrange(mlines):
                G[i*mcols:(i+1)*mcols,0] = np.arange((mcols))
                G[i*mcols:(i+1)*mcols,1] = i 
                G[i*mcols:(i+1)*mcols,2] = (i) * (np.arange((mcols)))    
            G[:,3] = 1

        else:
            if ivar==0 and nfit==0:
                G=np.zeros((len(los_clean),5))
                G[:,0] = y
                G[:,1] = x
                G[:,2] = y*x
                G[:,3] = 1
                G[:,4] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[5] = pars[1]; sol[7] = pars[2] = pars[3]; sol[9] = pars[4]
                print 'Remove ramp %f r, %f az  + %f r*az + %f + %f z'%(pars[0],pars[1],pars[2],pars[3],pars[4])

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols))
                    G[i*mcols:(i+1)*mcols,1] = i
                    G[i*mcols:(i+1)*mcols,2] = (i) * (np.arange((mcols)))
                G[:,3] = 1
                G[:,4] = elev_map.flatten()

            if ivar==0 and nfit==1:
                G=np.zeros((len(los_clean),6))
                G[:,0] = y
                G[:,1] = x
                G[:,2] = y*x
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[5] = pars[1]; sol[7] = pars[2]; sol[8] = pars[3]; sol[9] = pars[4]; sol[10] = pars[-1]
                print 'Remove ramp %f r, %f az  + %f r*az + %f + %f z+ %f z**2'%(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5])

                # build total G matrix
                G=np.zeros((len(los),6))
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols))
                    G[i*mcols:(i+1)*mcols,1] = i
                    G[i*mcols:(i+1)*mcols,2] = (i) * (np.arange((mcols)))
                G[:,3] = 1
                G[:,4] = elev_map.flatten()
                G[:,5] = elev_map.flatten()**2

            elif ivar==1 and nfit==0:
                G=np.zeros((len(los_clean),6))
                G[:,0] = y
                G[:,1] = x
                G[:,2] = y*x
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean*x

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[5] = pars[1]; sol[7] = pars[2]; sol[8] = pars[3]; sol[9] = pars[4]; sol[11] = pars[5]
                print 'Remove ramp %f r, %f az  + %f r*az + %f + %f z + %f z*az'%(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5])

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,3] = 1
                G[:,4] = elev_map.flatten()
                G[:,5] = elev_map.flatten()
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols))
                    G[i*mcols:(i+1)*mcols,1] = i
                    G[i*mcols:(i+1)*mcols,2] = (i) * (np.arange((mcols)))
                    G[i*mcols:(i+1)*mcols,5] *= i

            elif ivar==1 and nfit==1:
                G=np.zeros((len(los_clean),7))
                G[:,0] = y
                G[:,1] = x
                G[:,2] = y*x
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean*x
                G[:,6] = (topo_clean*x)**2

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[2] = pars[0]; sol[5] = pars[1]; sol[7] = pars[2]; sol[8] = pars[3]; sol[9] = pars[4]; sol[11] = pars[5]; sol[12] = pars[6]
                print 'Remove ramp %f r, %f az  + %f r*az + %f + %f z + %f z*az+ %f (z*az)**2'%(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5],pars[6])

                # build total G matrix
                G=np.zeros((len(los),7))
                G[:,3] = 1
                G[:,4] = elev_map.flatten()
                G[:,5] = elev_map.flatten()
                G[:,6] = elev_map.flatten()**2
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = np.arange((mcols))
                    G[i*mcols:(i+1)*mcols,1] = i
                    G[i*mcols:(i+1)*mcols,2] = (i) * (np.arange((mcols)))
                    G[i*mcols:(i+1)*mcols,5] *= i
                    G[i*mcols:(i+1)*mcols,6] *= i**2

    elif order==5:
    #0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2

        if radar is None:
            G=np.zeros((len(los_clean),3))
            G[:,0] = y**2
            G[:,1] = y
            G[:,2] = 1

            # ramp inversion
            x0 = lst.lstsq(G,los_clean)[0]
            try:
                _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
            except:
                pars = x0
            sol[1] = pars[0]; sol[2] = pars[1]; sol[8] = pars[2]
            print 'Remove ramp %f r**2 %f r  + %f'%(pars[0],pars[1],pars[2])

            # build total G matrix
            G=np.zeros((len(los),3))
            for i in xrange(mlines):
                G[i*mcols:(i+1)*mcols,0] = (np.arange((mcols)))**2
                G[i*mcols:(i+1)*mcols,1] = np.arange((mcols)) 
            G[:,2] = 1

        else:
            if ivar==0 and nfit==0:
                G=np.zeros((len(los_clean),4))
                G[:,0] = y**2
                G[:,1] = y
                G[:,2] = 1
                G[:,3] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[1] = pars[0]; sol[2] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]
                print 'Remove ramp %f r**2, %f r  + %f + %f z'%(pars[0],pars[1],pars[2],pars[3])

                # build total G matrix
                G=np.zeros((len(los),4))
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = (np.arange((mcols)))**2
                    G[i*mcols:(i+1)*mcols,1] = np.arange((mcols)) 
                G[:,2] = 1
                G[:,3] = elev_map.flatten()

            elif ivar==0 and nfit==1:
                G=np.zeros((len(los_clean),5))
                G[:,0] = y**2
                G[:,1] = y
                G[:,2] = 1
                G[:,3] = topo_clean
                G[:,4] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[1] = pars[0]; sol[2] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]; sol[10] = pars[-1]
                print 'Remove ramp %f r**2, %f r  + %f + %f z + %f z**2'%(pars[0],pars[1],pars[2],pars[3],pars[4])

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = (np.arange((mcols)))**2
                    G[i*mcols:(i+1)*mcols,1] = np.arange((mcols)) 
                G[:,2] = 1
                G[:,3] = elev_map.flatten()
                G[:,4] = elev_map.flatten()**2
            
            elif ivar==1 and nfit==0:
                G=np.zeros((len(los_clean),5))
                G[:,0] = y**2
                G[:,1] = y
                G[:,2] = 1
                G[:,3] = topo_clean
                G[:,4] = topo_clean*x

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[1] = pars[0]; sol[2] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]; sol[11] = pars[4]
                print 'Remove ramp %f r**2, %f r  + %f + %f z + %f z*az'%(pars[0],pars[1],pars[2],pars[3],pars[4])

                # build total G matrix
                G=np.zeros((len(los),5))
                G[:,2] = 1
                G[:,3] = elev_map.flatten()
                G[:,4] = elev_map.flatten()
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = (np.arange((mcols)))**2
                    G[i*mcols:(i+1)*mcols,1] = np.arange((mcols))
                    G[i*mcols:(i+1)*mcols,4] *= i

            elif ivar==1 and nfit==1:
                G=np.zeros((len(los_clean),6))
                G[:,0] = y**2
                G[:,1] = y
                G[:,2] = 1
                G[:,3] = topo_clean
                G[:,4] = topo_clean*x
                G[:,5] = (topo_clean*x)**2

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[1] = pars[0]; sol[2] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]; sol[11] = pars[4]; sol[12] = pars[5]
                print 'Remove ramp %f r**2, %f r  + %f + %f z + %f z*az + %f (z*az)**2'%(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5])

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,2] = 1
                G[:,3] = elev_map.flatten()
                G[:,4] = elev_map.flatten()
                G[:,5] = elev_map.flatten()**2
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = (np.arange((mcols)))**2
                    G[i*mcols:(i+1)*mcols,1] = np.arange((mcols))
                    G[i*mcols:(i+1)*mcols,4] *= i
                    G[i*mcols:(i+1)*mcols,5] *= i**2

    elif order==6:
    #0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2

        if radar is None:
            G=np.zeros((len(los_clean),3))
            G[:,0] = x**2
            G[:,1] = x
            G[:,2] = 1

            # ramp inversion
            x0 = lst.lstsq(G,los_clean)[0]
            try:
                _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
            except:
                pars = x0
            sol[4] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]
            print 'Remove ramp %f az**2 %f az  + %f'%(pars[0],pars[1],pars[2])

            # build total G matrix
            G=np.zeros((len(los),3))
            for i in xrange(mlines):
                G[i*mcols:(i+1)*mcols,0] = (i)**2
                G[i*mcols:(i+1)*mcols,1] = (i) 
            G[:,2] = 1

        else:
            if ivar==0 and nfit==0:
                G=np.zeros((len(los_clean),4))
                G[:,0] = x**2
                G[:,1] = x
                G[:,3] = 1
                G[:,4] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[4] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]
                print 'Remove ramp %f az**2, %f az  + %f + %f z'%(pars[0],pars[1],pars[2],pars[3])

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = (i)**2
                    G[i*mcols:(i+1)*mcols,1] = i 
                G[:,2] = 1
                G[:,3] = elev_map.flatten()

            elif ivar==0 and nfit==1:
                G=np.zeros((len(los_clean),5))
                G[:,0] = x**2
                G[:,1] = x
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[4] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]; sol[10] = pars[4]
                print 'Remove ramp %f az**2, %f az  + %f + %f z + %f z**2'%(pars[0],pars[1],pars[2],pars[3],pars[4])

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = (i)**2
                    G[i*mcols:(i+1)*mcols,1] = i 
                G[:,2] = 1
                G[:,3] = elev_map.flatten()
                G[:,4] = elev_map.flatten()**2
            
            elif ivar==1 and nfit==0:
                G=np.zeros((len(los_clean),5))
                G[:,0] = x**2
                G[:,1] = x
                G[:,2] = 1
                G[:,3] = topo_clean
                G[:,4] = topo_clean*x

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                sol[4] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]; sol[11] = pars[4];
                print 'Remove ramp %f az**2, %f az + %f + %f z + %f z*az'%(pars[0],pars[1],pars[2],pars[3],pars[4])

                # build total G matrix
                G=np.zeros((len(los),5))
                G[:,2] = 1
                G[:,3] = elev_map.flatten()
                G[:,4] = elev_map.flatten()
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = (i)**2
                    G[i*mcols:(i+1)*mcols,1] = i 
                    G[:,4] *= i 

            elif ivar==1 and nfit==1:
                G=np.zeros((len(los_clean),6))
                G[:,0] = x**2
                G[:,1] = x
                G[:,2] = 1
                G[:,3] = topo_clean
                G[:,4] = topo_clean*x
                G[:,5] = (topo_clean*x)**2

                # ramp inversion
                x0 = lst.lstsq(G,los_clean)[0]
                try:
                    _func = lambda x: np.sum(((np.dot(G,x)-los_clean)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-los_clean)/rms)
                    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                except:
                    pars = x0
                #0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2
                sol[4] = pars[0]; sol[5] = pars[1]; sol[8] = pars[2]; sol[9] = pars[3]; sol[11] = pars[4]; sol[12] = pars[5]
                print 'Remove ramp %f az**2, %f az + %f + %f z + %f z*az + %f (z*az)**2 '%(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5])

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,2] = 1
                G[:,3] = elev_map.flatten()
                G[:,4] = elev_map.flatten()
                G[:,5] = elev_map.flatten()**2
                for i in xrange(mlines):
                    G[i*mcols:(i+1)*mcols,0] = (i)**2
                    G[i*mcols:(i+1)*mcols,1] = i 
                    G[:,4] *= i 
                    G[:,5] *= i**2 


    corr = np.dot(G,pars).reshape(mlines,mcols)
    res = los - np.dot(G,pars)
    rms = np.sqrt(np.nanmean(res**2))

    # plt.imshow(los.reshape(mlines,mcols))
    # plt.show()
    # plt.imshow(corr)
    # plt.show()

    return sol, corr, rms


#####################################################################################

# initialise full vector correction 
# 16 values: date1 dates2  mlines y**3 y**2 y x**3 x**2 x xy**2 xy cst z z**2 z*az az*z**2 
spint = np.zeros((kmax,16))
M = 13

# fill dates
spint[:,0],spint[:,1] = date_1,date_2 

# initilise correction cube
# corr_maps = np.zeros((mlines,mcols,kmax))
rms = np.zeros((kmax,3))
rms[:,0],rms[:,1] = date_1,date_2 

if estim=='yes':

    print 
    #########################################
    print '#################################'
    print 'Empirical estimations'
    print '#################################'
    #########################################
    print

    for kk in xrange((kmax)):
        date1, date2 = date_1[kk], date_2[kk]
        idate = str(date1) + '-' + str(date2) 

        if sformat == 'ROI_PAC':
            folder =  'int_'+ str(date1) + '_' + str(date2) 
            rscfile=int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.unw.rsc'
            infile=int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.unw'

            ds = gdal.Open(infile, gdal.GA_ReadOnly)
            # Get the band that have the data we want
            ds_band1 = ds.GetRasterBand(1)
            ds_band2 = ds.GetRasterBand(2)

            los_map = np.zeros((mlines,mcols))
            los_map[:ds.RasterYSize,:ds.RasterXSize] = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:mlines,:mcols]
            # los_map[los_map==0] = np.float('NaN')
            lines, cols = ds.RasterYSize, ds.RasterXSize


        elif sformat == 'GTIFF':
            infile = int_path + prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.tiff'

            ds = gdal.Open(infile, gdal.GA_ReadOnly)
            # Get the band that have the data we want
            ds_band2 = ds.GetRasterBand(1)
            lines, cols = ds.RasterYSize, ds.RasterXSize

            los_map[:ds.RasterYSize,:ds.RasterXSize] = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:mlines,:mcols]
            # los_map[los_map==0] = np.float('NaN')

        elif sformat == 'GAMMA':
            # scfile=prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.unw.par'
            # par_file = ref 
            lines,cols = gm.readpar(int_path)
            infile= int_path + prefix + str(date1) + '_' + str(date2) + suffix + rlook + '.unw'
            los_map = gm.readgamma(infile,int_path)

        print 
        print 'lines:{}, cols:{}, int:{}:'.format(lines, cols, idate)

        # load coherence or whatever
        spacial_mask = np.ones((mlines,mcols))*np.float('NaN')

        rms_map = np.ones((mlines,mcols))
        try:
            if rmsf=='yes':
                rms_map[:lines,:cols] = ds_band1.ReadAsArray(0, 0, lines, cols)[:mlines,:mcols]
                # rmsi = np.fromfile(folder + 'cor',dtype=np.float32)
                # _rms_map = rmsi.reshape(len(rmsi)/1420,1420)
                # if len(rmsi)/1420 < mlines:
                #     rms_map[:len(rmsi)/1420,:1420] = _rms_map
                # else:
                #     rms_map = _rms_map[:mlines,:mcols]
                k = np.nonzero(np.logical_or(rms_map==0.0, rms_map==9999))
                rms_map[k] = float('NaN')
            else:
                threshold_rms = -1
        except:
            threshold_rms = -1

        # time.sleep(1.)
        # clean for estimation
        _los_map = np.copy(los_map)
        _los_map[los_map==0] = np.float('NaN')
        maxlos,minlos = np.nanpercentile(_los_map,perc),np.nanpercentile(_los_map,(100-perc))
        # print maxlos,minlos


        # print np.shape(los_map), np.shape(elev_map), np.shape(rms_map), np.shape(pix_az)
        ## CRITICAL STEP ####
        # select points for estimation only: minmax elev, los not NaN, rms<rmsthreshold ....
        index = np.nonzero(
        np.logical_and(elev_map<maxelev,
        np.logical_and(elev_map>minelev,    
        np.logical_and(los_map!=0, 
        np.logical_and(los_map>minlos,
        np.logical_and(los_map<maxlos,
        np.logical_and(rms_map>threshold_rms, 
        np.logical_and(pix_az>ibeg,
        np.logical_and(pix_az<iend,
        np.logical_and(pix_rg>jbeg,
        np.logical_and(pix_rg<jend,
        np.logical_and(mask>threshold_mask,
        np.logical_and(~np.isnan(los_map),
        np.logical_or(pix_az<ibeg_mask,pix_az>iend_mask)
        )
        )
        )
        )
        )
        )
        )
        )
        )
        )
        )
        )
        )

        spacial_mask[index] = np.copy(los_map[index])

        # extract range and azimuth coordinates
        temp = np.array(index).T
        az = temp[:,0]; rg = temp[:,1]
        # print az
        # print rg

        # clean maps
        los_temp = np.matrix.copy(los_map)
        elev_temp = np.matrix.copy(elev_map)
        los_clean = los_temp[index].flatten()
        elev_clean = elev_temp[index].flatten()
        rms_clean = rms_map[index].flatten()
        del los_temp, elev_temp

        # Take care to not do high polynomial estimations for short int.
        # find the begining of the image
        itemp = ibeg
        for lign in xrange(ibeg,iend,10):
          if np.isnan(np.nanmean(_los_map[lign:lign+10,:])):
              itemp = lign  
          else:
              break
        del _los_map

        # print itemp
        # 0: ref frame [default], 1: range ramp ax+b , 2: azimutal ramp ay+b, 
        # 3: ax+by+c, 4: ax+by+cxy+d 5: ax**2+bx+d, 6: ay**2+by+c
        if flat>5 and iend-itemp < .6*(iend-ibeg):
          print
          print 'Int. too short in comparison to ref, set flat to 5'
          temp_flat=5
        elif flat>5 and iend-itemp < .9*mcols:
          print
          print 'Lenght int. inferior to width, set flat to 5 and nfit to 0'
          temp_flat=5
        else:
          temp_flat=flat

        if ivar>0 and iend-itemp < .6*(iend-ibeg):
          print
          print 'Int. too short in comparison to ref, set ivar and nfit to 0'
          nfit_temp=0
          ivar_temp=0
        else:
          nfit_temp=nfit
          ivar_temp=ivar

        # save size int to use as weight in the temporal inversion
        spint[kk,2] = iend-itemp

        # hard-coding subsample 
        samp = 1
        sol, corr, rms[kk,2] = estim_ramp(los_map.flatten(),
        los_clean[::samp],elev_clean[::samp],az[::samp],rg[::samp],
        temp_flat,rms_clean[::samp],nfit_temp,ivar_temp)

        print 'RMS: ',rms[kk,2]

        # clean corr : non car il faut extrapoler pour l'inversion ens erie temp
        # k = np.nonzero(np.logical_or(los_map==0.,abs(los_map)>999.))
        # corr[k] = 0.
        # corr_maps[:,:,kk] = np.copy(corr)

        # print sol
        # 0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2
        func = sol[0]*rg**3 + sol[1]*rg**2 + sol[2]*rg + sol[3]*az**3 + sol[4]*az**2 \
        + sol[5]*az + sol[6]*(rg*az)**2 + sol[7]*rg*az + sol[11]*az*elev_clean + \
        sol[12]*((az*elev_clean)**2)

        if radar is not None: 
           # plot phase/elevation
           nfigure=nfigure+1
           fig2 = plt.figure(nfigure,figsize=(9,4))
           ax = fig2.add_subplot(1,1,1)
           z = np.linspace(np.min(elev_clean), np.max(elev_clean), 100)
           ax.scatter(elev_clean,los_clean - func, s=0.005, alpha=0.05,rasterized=True)

           # if nfit==0:
           ax.plot(z,sol[8]+sol[9]*z, '-r', lw =3.,label='{1:6f}*z + {1:.3f}'.format(sol[9],sol[8])) 
           # else:
           #      ax.plot(z,sol[8]+sol[9]*z+sol[10]*z**2, '-r', lw =3.,label='{1:6f}*z**2 + {1:6f}*z + {1:.3f}'.format(sol[10],sol[9],sol[8]))

           ax.set_xlabel('Elevation (m)')
           ax.set_ylabel('LOS (rad)')
           plt.legend(loc='best')
           fig2.savefig(out_path + idate+'phase-topo.png', format='PNG')

        #corected map
        vmax = np.max(np.array([abs(maxlos),abs(minlos)]))

        nfigure=nfigure+1
        fig = plt.figure(nfigure,figsize=(11,4))

        ax = fig.add_subplot(1,4,1)
        hax = ax.imshow(rms_map, cm.Greys,vmax=1,vmin=0.)
        cax = ax.imshow(los_map,cmap=cm.gist_rainbow,vmax=vmax,vmin=-vmax,interpolation=None,alpha=1.)
        ax.set_title('LOS')
        setp( ax.get_xticklabels(), visible=None)
        fig.colorbar(cax, orientation='vertical',aspect=10)

        ax = fig.add_subplot(1,4,2)
        cax = ax.imshow(spacial_mask,cmap=cm.gist_rainbow,vmax=vmax,vmin=-vmax)
        ax.set_title('LOS ESTIMATION')
        setp( ax.get_xticklabels(), visible=None)
        fig.colorbar(cax, orientation='vertical',aspect=10)

        ax = fig.add_subplot(1,4,3)
        cax = ax.imshow(corr,cmap=cm.gist_rainbow,vmax=vmax,vmin=-vmax)
        ax.set_title('RAMP+TOPO')
        setp( ax.get_xticklabels(), visible=None)
        fig.colorbar(cax, orientation='vertical',aspect=10)

        # for plot we can clean
        k = np.nonzero(np.logical_or(los_map==0.,abs(los_map)>999.))
        corr[k] = 0.

        ax = fig.add_subplot(1,4,4)
        hax = ax.imshow(rms_map, cm.Greys,vmax=1,vmin=0.)
        cax = ax.imshow(los_map - corr,cmap=cm.gist_rainbow,vmax=vmax,vmin=-vmax,alpha=1.,interpolation=None)
        ax.set_title('CORR LOS')
        setp( ax.get_xticklabels(), visible=None)
        fig.colorbar(cax, orientation='vertical',aspect=10)
        fig.tight_layout()

        if sformat == 'ROI_PAC':
            fig.savefig(int_path + folder + idate +'corrections.png', format='PNG')
        else:
            fig.savefig(out_path + idate +'corrections.png', format='PNG')

        if plot=='yes':
            plt.show()

        # fill correction matrix
        spint[kk,3:] = sol

        plt.close('all')
        del corr, los_map, rms_map
        del los_clean, rms_clean
        del elev_clean
        del az, rg
        try:
            del ds
        except:
            pass


    # save spint 
    np.savetxt('liste_coeff_ramps.txt', spint , header='#date1   |   dates2   |   Lenght   |   y**3   |   y**2   |   y\
       |   **3   |   x**2   |   x   |   xy**2   |   xy   |   cst   |   z   |   z**2   |   z*az   |   z**2*az', fmt=('%i','%i','%.8f','%.8f','%.8f','%.8f','%.8f',\
        '%.8f','%.8f','%.8f','%.8f','%3.8f','%.8f','%.8f','%.8f','%.8f'))

    # # save correction matrix
    # fid = open('corection_matrix', 'wb') 
    # corr_maps.flatten().astype('float32').tofile(fid)
    # fid.close()

    # save rms
    np.savetxt('rms_unwcor.txt', rms, header='#date1   |   dates2   |   RMS', fmt=('%i','%i','%.8f'))

#####################################################################################

print 
print 'read input files liste_coeff_ramps.txt and corection_matrix'
print 

# load spint
date_1,date_2,length,a,b,c,d,e,f,g,h,i,j,k,l,m=np.loadtxt('liste_coeff_ramps.txt',comments="#",unpack=True,dtype='i,i,f,f,f,f,f,f,f,f,f,f,f,f,f,f')
spint = np.vstack([date_1,date_2,length,a,b,c,d,e,f,g,h,i,j,k,l,m]).T
rec_spint = np.copy(spint)

# # # load correction matrix
# corr_maps = np.fromfile("corection_matrix",dtype=np.float32).reshape((mlines,mcols,kmax))

####################################################################################

if tsinv=='yes':

    print
    #########################################
    print '#################################'
    print 'Temporal inversion of all coefficients'
    print '#################################'
    #########################################
    print

    spint_inv = np.zeros((np.shape(spint)))

    G_=np.zeros((kmax,nmax))
    deltat = np.zeros((kmax))
    for k in xrange((kmax)):
      for n in xrange((nmax)):
        if (date_1[k]==im[n]): 
          G_[k,n]=-1
          t1 = bt[n]
        elif (date_2[k]==im[n]):
          G_[k,n]=1
          t2 = bt[n]
      deltat[k] = abs(t2 -t1)

    # 1) create weight based on temporal baseline: give stronger weight to short temporal baselines 
    #, where we dont expect def.
    w1 = np.exp(-(deltat/2.))
    # print deltat
    # print w1
    # print 

    # 2) create a weight based on the size of the int: give a stronger weight to long interferograms
    # print length
    w2 = length/mlines
    # print w2
    # sys.exit()
    sig_ = 1./w1 + 1./w2

    for j in xrange(3,shape(spint_inv)[1]):

        d = np.zeros(((kmax+1)))
        sig = np.ones(((kmax+1)))
        G = np.zeros(((kmax+1),nmax))
    
        d[:kmax] = as_strided(spint[:,j])
        G[:kmax,:nmax] = G_ 
        G[-1,0] = 1 # ini phi first image to 0
        sig[:kmax] = sig_

        try:
            x0 = lst.lstsq(G,d)[0]
            _func = lambda x: np.sum(((np.dot(G,x)-d)/sig)**2)
            _fprime = lambda x: 2*np.dot(G.T/sig, (np.dot(G,x)-d)/sig)
            pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=500,full_output=True,iprint=0)[0]

            # reconstruct corr for selected int
            spint_inv[:,j] = np.dot(G,pars)[:kmax]
            # print 
            # print spint_inv[:,j+3] - spint[:,j+3]
            # print

        except:
            pass

    spint_inv[:,:3] = spint[:,:3]
    # I dont think the cst should be inverted ??
    # spint_inv[:,8] = spint[:,8]


#####################################################################################
####

print 
#########################################
print '#################################'
print 'APPLY CORRECTION AND SAVE NEW INT.'
print '#################################'
#########################################
print


# apply correction
for kk in xrange((kmax)):
    date1, date2 = date_1[kk], date_2[kk]
    folder = 'int_'+ str(date1) + '_' + str(date2) 
    idate = str(date1) + '-' + str(date2) 

    if sformat == 'ROI_PAC':
        
        rscfile=int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix +  rlook + '.unw.rsc'
        infile=int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix +  rlook + '.unw'
        outfile = int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix +  suffout +  rlook + '.unw'  
        outrsc = int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix +  suffout +  rlook + '.unw.rsc' 
        
        ds = gdal.Open(infile, gdal.GA_ReadOnly)
        # Get the band that have the data we want
        ds_band1 = ds.GetRasterBand(1)
        ds_band2 = ds.GetRasterBand(2)
        los_map = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        rms_map = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        lines, cols = ds.RasterYSize, ds.RasterXSize

    elif sformat == 'GTIFF':

        infile = int_path + prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.tiff'
        outfile = out_path + prefix + str(date1) + '-' + str(date2) + suffix + suffout +  rlook + '.tiff' 
        ds = gdal.Open(infile, gdal.GA_ReadOnly)
        ds_band2 = ds.GetRasterBand(1)
        los_map = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        rms_map = np.zeros((ds.RasterYSize,ds.RasterXSize))
        lines, cols = ds.RasterYSize, ds.RasterXSize

    elif sformat == 'GAMMA':

        infile = int_path + prefix + str(date1) + '_' + str(date2) + suffix +  rlook + '.unw'
        outfile = out_path + prefix + str(date1) + '_' + str(date2) + suffix +  suffout + rlook + '.unw' 
        # par_file = ref 
        lines,cols = gm.readpar(int_path)
        los_map = gm.readgamma(infile,int_path)
        rms_map = np.zeros((lines,cols))

    # rms_map = np.zeros((mlines,mcols))
    # print 'Apply correction and clean coh < 0.03....'
    # los_map[rms_map<0.03] = 0.0
    # rms_map[rms_map<0.03] = 0.0

    print 
    print 'mlines:{}, mcols:{}, int:{}:'.format(lines, cols, idate)

    # compute correction
    rg = np.tile(np.arange(cols), (lines,1))
    az = np.tile(np.arange(lines), (cols,1)).T
    z = np.zeros((lines, cols))
    z = elev_map[:lines,:cols]

    # 0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2
    if tsinv=='yes':

            sol = spint_inv[kk,3:]
            corr_inv = sol[0]*rg**3 + sol[1]*rg**2 + sol[2]*rg + sol[3]*az**3 + sol[4]*az**2 \
            + sol[5]*az + sol[6]*(rg*az)**2 + sol[7]*rg*az + sol[8] + sol[9]*z + sol[10]*(z**2) \
            + sol[11]*az*z + sol[12]*((az*z)**2)

            sol = spint[kk,3:]
            corr = sol[0]*rg**3 + sol[1]*rg**2 + sol[2]*rg + sol[3]*az**3 + sol[4]*az**2 \
            + sol[5]*az + sol[6]*(rg*az)**2 + sol[7]*rg*az + sol[8] + sol[9]*z + sol[10]*(z**2) \
            + sol[11]*az*z + sol[12]*((az*z)**2)
            
    
    else:   
            
            sol = spint[kk,3:]
            corr_inv = sol[0]*rg**3 + sol[1]*rg**2 + sol[2]*rg + sol[3]*az**3 + sol[4]*az**2 \
            + sol[5]*az + sol[6]*(rg*az)**2 + sol[7]*rg*az + sol[8] + sol[9]*z + sol[10]*(z**2) \
            + sol[11]*az*z + sol[12]*((az*z)**2)

            corr = corr_inv

    # reset to 0 areas where no data (might change after time series inversion?)
    flatlos = los_map - corr_inv
    flatlos[los_map==0], rms_map[los_map==0] = 0.0, 0.0
    flatlos[isnan(flatlos)],rms_map[isnan(flatlos)] = 0.0, 0.0
    flatlos[isnan(los_map)],rms_map[isnan(los_map)] = 0.0, 0.0
    rms_map[isnan(rms_map)],flatlos[isnan(rms_map)] = 0.0, 0.0
    
    # print(suffout, rlook)
    # print(outfile)
    # print(outrsc)

    if sformat == 'ROI_PAC':
        dst_ds = driver.Create(outfile, cols, lines, 2, gdal.GDT_Float32)
        dst_band1 = dst_ds.GetRasterBand(1)
        dst_band2 = dst_ds.GetRasterBand(2)
        dst_band1.WriteArray(rms_map,0,0)
        dst_band2.WriteArray(flatlos,0,0)
        shutil.copy(rscfile,outrsc)
        dst_band1.FlushCache()
        dst_band2.FlushCache()

    elif sformat == 'GTIFF':
        dst_ds = driver.Create(outfile, cols, lines, 1, gdal.GDT_Float32)
        dst_band2 = dst_ds.GetRasterBand(1)
        dst_band2.WriteArray(flatlos,0,0)
        dst_ds.SetGeoTransform(gt)
        dst_ds.SetProjection(proj)
        dst_band2.FlushCache()

    elif sformat == 'GAMMA':
        fid = open(outfile, 'wb')
        flatlos.flatten().astype('>f4').tofile(fid)

    nfigure=nfigure+1
    fig = plt.figure(nfigure,figsize=(9,4))

    _los_map = np.copy(flatlos)
    _los_map[los_map==0] = np.float('NaN')
    maxlos,minlos=np.nanpercentile(_los_map,perc),np.nanpercentile(_los_map,(100-perc))
    vmax = np.max(np.array([abs(maxlos),abs(minlos)]))
    # vmax = np.percentile(los_map, 98)

    ax = fig.add_subplot(1,4,1)
    cax = ax.imshow(los_map,cmap=cm.gist_rainbow,vmax=vmax,vmin=-vmax,alpha=0.7,interpolation=None)
    ax.set_title('LOS')
    setp( ax.get_xticklabels(), visible=None)

    ax = fig.add_subplot(1,4,2)
    cax = ax.imshow(corr,cmap=cm.gist_rainbow,vmax=vmax,vmin=-vmax,alpha=0.7,interpolation=None)
    ax.set_title('RAMP+TOPO ORIG')
    setp( ax.get_xticklabels(), visible=None)

    ax = fig.add_subplot(1,4,3)
    cax = ax.imshow(corr_inv,cmap=cm.gist_rainbow,vmax=vmax,vmin=-vmax,alpha=0.7,interpolation=None)
    ax.set_title('RAMP+TOPO RECONST.')
    setp( ax.get_xticklabels(), visible=None)

    ax = fig.add_subplot(1,4,4)
    cax = ax.imshow(flatlos,cmap=cm.gist_rainbow,vmax=vmax,vmin=-vmax,alpha=0.7,interpolation=None)
    ax.set_title('CORR LOS RECONST.')
    setp( ax.get_xticklabels(), visible=None)
    fig.colorbar(cax, orientation='vertical',aspect=10)
    fig.tight_layout()

    if sformat == 'ROI_PAC':
        fig.savefig( int_path + folder + idate + '_reconstruc_corrections.png', format='PNG')
    else:
        fig.savefig( out_path + idate + '_reconstruc_corrections.png', format='PNG')

    if plot=='yes':
        plt.show()

    plt.close('all')

    try:
        del dst_ds, ds, drv
    except:
        pass
    del los_map, rms_map

