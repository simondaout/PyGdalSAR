#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT (Oxford)
################################################################################

"""\
correct_ifg_from_gacos.py
-------------
Correct unwrapped IFGs from Gacos atmospheric models. 

Usage: 
    correct_ifg_from_gacos.py --int_list=<path> [--int_path=<path>] [--gacos_path=<path>]  [--gacos2los=<value>] [--prefix=<value>] \
    [--suffix=<value>] [--rlook=<value>] [--plot=<yes|no>] [--cohpixel=<yes/no>] [--threshold_coh=<value>] \
    [--refstart=<values>] [--refend=<values>] [--format=<value>]  \
    [--ramp=<cst|lin>] [--crop=<values>] [--fitmodel=<yes|no>] [--perc=<value>]  \
    [--suffix_output=<value>] [--nproc=<nb_cores>]

correct_ifg_from_gacos.py -h | --help

Options:
-h --help               Show this screen.
--int_list=<path>       Text file containing list of interferograms dates in two colums, $los_map1 $date2
--int_path=<path>       Relative path to input interferograms directory
--gacos_path=<dir>      Path to GACOS atmospheric models [default: ./GACOS/]
--prefix=<value>        Prefix name $prefix$date1-$date2$suffix_$rlookrlks.unw [default: '']
--suffix=<value>        Suffix name $prefix$date1-$date2$suffix_$rlookrlks.unw [default: '']
--suffix_output=<value> Suffix output file name $prefix$date1-$date2$suffix$suffix_output [default:_gacos]
--rlook=<value>         look int. $prefix$date1-$date2$suffix_$rlookrlks.unw [default: 0]
--refstart=<value>      Stating line number of the area where phase is set to zero [default: 0] 
--refend=<value>        Ending line number of the area where phase is set to zero [default: length]
--ramp=<cst|lin>        Estimate a constant or linear ramp in x and y [default: cst].
--crop=<value>          Crop option for empirical estimation (eg: 0,ncol,0,nlign)
--gacos2los=<value>     Scaling value between gacos los_map (m) and desired output (e.g los_map in mm positive toward the sat. and LOS angle 23°: -1000*cos(23)=-920.504)
--cohpixel=<yes/no>     If Yes, use amplitude interferogram to weight and mask pixels (e.g Coherence, Colinearity, Amp Filter) [default: no]
--threshold_coh=<value> Thresold on rmspixel for ramp estimation [default: 2.]   
--plot=<yes|no>         Display results [default: yes]
--format VALUE          Format input files: ROI_PAC, GAMMA, GTIFF [default: GAMMA]
--perc=<value>          Percentile of hidden LOS pixel for the estimation and clean outliers [default:98.]
--fitmodel=<yes|no>     If yes, then estimate the proportionlality between gacos and los_map in addition to a polynomial ramp
--nproc=<values>            number of processor (default: 1)
"""

import gdal
gdal.UseExceptions()

import sys,os,logging
import numpy as np
import scipy.optimize as opt
import scipy.linalg as lst
from numpy.lib.stride_tricks import as_strided
from os import environ

import matplotlib
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from datetime import datetime
import time
from mpl_toolkits.axes_grid1 import make_axes_locatable

import docopt
import shutil

from contextlib import contextmanager
from functools import wraps, partial
import multiprocessing

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning) 

#####################################################################################
# READ INPUT PARAM
#####################################################################################

# read arguments
arguments = docopt.docopt(__doc__)

int_list=arguments["--int_list"]

if arguments["--int_path"] == None:
    int_path='.'
else:
    int_path=arguments["--int_path"] + '/'

if arguments["--gacos_path"] ==  None:
   gacos_path = "./GACOS/"
else:
   gacos_path = arguments["--gacos_path"] + '/'

if arguments["--refstart"] == None:
    refstart = 0
else:
    refstart = int(arguments["--refstart"])

if arguments["--refend"] == None:
    refend = 200
else:
    refend = int(arguments["--refend"])

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

if arguments["--nproc"] == None:
    nproc = 1
else:
    nproc = int(arguments["--nproc"])

if arguments["--plot"] ==  'yes':
    plot = 'yes'
    logger.warning('plot is yes. Set nproc to 1')
    nproc = 1
    if environ["TERM"].startswith("screen"):
        matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt
else:
    plot = 'no'
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt

if arguments["--gacos2los"] ==  None:
    gacos2los = 1
else:
    gacos2los = float(arguments["--gacos2los"])

if arguments["--ramp"] ==  None:
    ramp = 'cst'
else:
    ramp = arguments["--ramp"]

if arguments["--refstart"] == None:
    refstart = 0
else:
    refstart = int(arguments["--refstart"])
if arguments["--refend"] == None:
    refend = 200
else:
    refend = int(arguments["--refend"])

if arguments["--format"] ==  None:
    sformat = 'GAMMA'
else:
    sformat = arguments["--format"]

if arguments["--perc"] ==  None:
    perc = 98.
else:
    perc = float(arguments["--perc"])

if arguments["--cohpixel"] ==  None:
    rmsf = 'no'
else:
    rmsf = arguments["--cohpixel"]

if arguments["--threshold_coh"] ==  None:
    threshold_rms = 2.
else:
    threshold_rms = float(arguments["--threshold_coh"])

if arguments["--fitmodel"] ==  None:
    fitmodel = 'no'
else:
    fitmodel = arguments["--fitmodel"]

if arguments["--crop"] ==  None:
    crop = False
else:
    crop = map(float,arguments["--crop"].replace(',',' ').split())

if arguments["--suffix_output"] ==  None:
    suffout = '_gacos'
else:
    suffout = arguments["--suffix_output"]

# Choose plot option
cmap = cm.gist_rainbow_r
# cmap = cm.jet
cmap.set_bad('white')

##################################################################################
###  Extras functions and context maganers
##################################################################################

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
        times.append(year + dec)
    return times

# Timer for all the functions
class ContextDecorator(object):
    def __call__(self, f):
        @wraps(f)
        def decorated(*args, **kwds):
            with self:
                try:
                    return f(*args, **kwds)
                except (KeyboardInterrupt, SystemExit):
                    raise
                except:
                    Exception('{0} Failed !'.format(f))
                    raise
        return decorated

class TimeIt(ContextDecorator):
    def __enter__(self):
        self.start = datetime.now()
        print('Starting time process: {0}'.format(self.start))
    def __exit__(self, type, value, traceback):
        print('Time process: {0}s'.format((datetime.now() - self.start).total_seconds()))

def checkinfile(file):
    if os.path.exists(file) is False:
        logger.critical("File: {0} not found, Exit !".format(file))
        print("File: {0} not found in {1}, Exit !".format(file,os.getcwd()))

# create generator for pool
@contextmanager
def poolcontext(*arg, **kargs):
    pool = multiprocessing.Pool(*arg, **kargs)
    yield pool
    pool.terminate()
    pool.join()

#####################################################################################
# FUNCTIONS
#####################################################################################


def gacos2ifg(kk):
    """
    Function that compute modeled ifgs
    """ 

    date1, date2 = date_1[kk], date_2[kk]
    idate = str(date1) + '-' + str(date2) 

    if sformat == 'GAMMA':
        import gamma as gm
        lines,cols = gm.readpar(gacos_path)
        infile1 = gacos_path +  str(date1) + '_crop.ztd.unw'
        infile2 = gacos_path +  str(date2) + '_crop.ztd.unw'
        checkinfile(infile1); checkinfile(infile2)
        gacos1 = gm.readgamma(infile1,int_path)
        gacos2 = gm.readgamma(infile2,int_path)

    if arguments["--zone"] ==  None:
        refzone = [0,cols,0,lines]
    else:
        refzone = map(float,arguments["--clean"].replace(',',' ').split())
    col_beg,col_end,line_beg,line_end = refzone[0],refzone[1],refzone[2],refzone[3]

    # extract range and azimuth coordinates from ref or radar file
    pix_az, pix_rg = np.indices((lines,cols))

    _los_map = np.copy(los_map)
    _los_map[np.logical_or(gacos2==0,gacos2>=9990)] = np.float('NaN')
    _los_map[np.logical_or(gacos1==0,gacos1>=9990)] = np.float('NaN')

    indexref = np.nonzero(
        np.logical_and(~np.isnan(gacos2),
        np.logical_and(pix_lin>refstart,
        np.logical_and(pix_lin<refend,
        np.logical_and(los_map!=0.0,
        np.logical_and(pix_col>col_beg, pix_col<col_end
        ))))))
   
    # compute average phase in the ref area 
    # we want los ref area to be zero
    cst1 = np.nanmean(gacos1[indexref].flatten())
    cst2 = np.nanmean(gacos2[indexref].flatten())

    # remove ref frame
    gacos2 = gacos2 - cst2
    gacos1 = gacos1 - cst1

    # compute differential los
    gacosm = (gacos2 - gacos1) * gacos2los
    # gacosm = (gacos2 - gacos1) * 4 * math.pi / wavelength

    # open gacos corrections
    gacosf = gacos_path + str(date1) + '-' + str(date2) + '_gacosmdel' +'.r4'
    fid = open(gacosf,'wb')
    gacosm.flatten().astype('float32').tofile(fid)
    fid.close()

    del gacosm, gacos2, gacos1

def correct_ifg(kk):
    """
    Function that corrects each ifg from gacos
    """

    date1, date2 = date_1[kk], date_2[kk]
    idate = str(date1) + '-' + str(date2) 

    if sformat == 'GAMMA':
        # scfile=prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.unw.par'
        # par_file = ref 
        lines,cols = gm.readpar(int_path)
        infile= int_path + prefix + str(date1) + '_' + str(date2) + suffix + rlook + '.unw'
        checkinfile(infile)
        los_map = gm.readgamma(infile,int_path)

    logger.info('lines:{0}, cols:{1}, IFG:{2}:'.format(lines, cols, idate))

    # open gacos corrections
    gacosf = gacos_path + str(date1) + '-' + str(date2) + '_gacosmdel' +'.r4'
    model = np.fromfile(gacosf,dtype=np.float32)

    # load coherence or whatever
    rms_map = np.ones((lines,cols))

    if rmsf=='yes':
        try:
            if sformat == 'ROI_PAC':
                rms_map[:ds.RasterYSize,:ds.RasterXSize] = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:lines,:cols]
                k = np.nonzero(np.logical_or(rms_map==0.0, rms_map==9999))
                rms_map[k] = float('NaN')
            elif sformat == 'GAMMA':
                rmsfile=  int_path + str(date1) + '_' + str(date2) + '.filt.cc'
                rms_map = gm.readgamma(rmsfile,int_path)

        except:
            logger.warning('Coherence file cannot be read')

    if arguments["--zone"] ==  None:
        refzone = [0,cols,0,lines]
    else:
        refzone = map(float,arguments["--clean"].replace(',',' ').split())
    col_beg,col_end,line_beg,line_end = refzone[0],refzone[1],refzone[2],refzone[3]

    # extract range and azimuth coordinates from ref or radar file
    pix_az, pix_rg = np.indices((lines,cols))

    _los_map = np.copy(los_map)
    _los_map[np.logical_or(los_map==0,los_map>=9990)] = np.float('NaN')
    losmin,losmax = np.nanpercentile(_los_map,2.),np.nanpercentile(_los_map,98.)
    gacosmin,gacosmax = np.nanpercentile(model,2),np.nanpercentile(model,98)

    funct = 0.
    pix_lin, pix_col = np.indices((lines,cols))

    # find proportionality between los_map and model
    index = np.nonzero(
        np.logical_and(~np.isnan(los_map),
        np.logical_and(pix_lin>line_beg,
        np.logical_and(pix_lin<line_end,
        np.logical_and(los_map<losmax,
        np.logical_and(los_map>losmin,
        np.logical_and(rms_map<threshold_rms,
        np.logical_and(rms_map>1.e-6,
        np.logical_and(model<gacosmax,
        np.logical_and(model>gacosmin,
        np.logical_and(los_map!=0.0,
        np.logical_and(model!=0.0,
        np.logical_and(pix_col>col_beg, pix_col<col_end
        )))))))))))))

    indexref = np.nonzero(
        np.logical_and(~np.isnan(los_map),
        np.logical_and(pix_lin>refstart,
        np.logical_and(pix_lin<refend,
        np.logical_and(los_map<losmax,
        np.logical_and(los_map>losmin,
        np.logical_and(rms_map<threshold_rms,
        np.logical_and(rms_map>1.e-6,
        np.logical_and(model<gacosmax,
        np.logical_and(model>gacosmin,
        np.logical_and(los_map!=0.0,
        np.logical_and(model!=0.0,
        np.logical_and(pix_col>col_beg, pix_col<col_end
        )))))))))))))

    # set los to zero in the ref area too
    # weights pixels by rms
    los_ref = los_map[indexref].flatten() 
    rms_ref = rms_map[indexref].flatten()
    amp_ref = 1./rms_ref
    amp_ref = amp_ref/np.nanpercentile(amp_ref,99)
    cst = np.nansum(los_ref*amp_ref) / np.nansum(amp_ref)
    los_map  = los_map - cst
    
    # extract los for empirical estimation
    temp = np.array(index).T
    x = temp[:,0]; y = temp[:,1]
    los_clean = los_map[index].flatten()
    model_clean = model[index].flatten()

    # compute average residual los - gacos in the ref area 
    los_ref = los_map[indexref].flatten() - model[indexref].flatten()
    rms_ref = rms_map[indexref].flatten()
    amp_ref = 1./rms_ref
    amp_ref = amp_ref/np.nanpercentile(amp_ref,99)
    # weigth average of the phase
    cst = np.nansum(los_ref*amp_ref) / np.nansum(amp_ref)

    # digitize los_map in bins, compute median and std
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
        remove_ramp[np.isnan(los_map)] = np.float('NaN')
        
        # los_map = gacos + cst
        remove = remove_ramp + model
       
        # Compute ramp for los_map = f(model) 
        functbins = a
        funct = a
        # set coef gacos to 1
        f = 1

    elif (ramp == 'lin' and  fitmodel=='no'):
        # here we want to minimize los_map-gacos = ramp
        # invers both digitized los_map and ref frame together
        d = losbins-modelbins
        sigmad = losstd

        G=np.zeros((len(d),3))
        G[:,0] = xbins
        G[:,1] = ybins
        G[:,2] = 1

        x0 = lst.lstsq(G,d)[0]
        # print x0
        _func = lambda x: np.sum(((np.dot(G,x)-d)/sigmad)**2)
        _fprime = lambda x: 2*np.dot(G.T/sigmad, (np.dot(G,x)-d)/sigmad)
        pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
        a = pars[0]; b = pars[1]; c = pars[2]
        print 'Remove ramp  %f az  + %f r + %f for date: %i'%(a,b,c,idates[l])
            
        # Compute ramp for los_map = f(model)
        functbins = a*xbins + b*ybins + c
        funct = a*x + b*y + c
        # set coef gacos to 1
        f = 1

        # build total G matrix
        G=np.zeros((len(los_map.flatten()),3))
        for i in xrange(nlign):
            G[i*ncol:(i+1)*ncol,0] = i - line_beg 
            G[i*ncol:(i+1)*ncol,1] = np.arange(ncol) - col_beg
        G[:,2] = 1

        # compute ramp
        remove_ramp = np.dot(G,pars).reshape(nlign,ncol)
        remove_ramp[model==0.] = 0.
        remove_ramp[np.isnan(los_map)] = np.float('NaN')
        
        # los_map = gacos + (a*rg + b*az + cst) 
        remove = model + remove_ramp

    elif fitmodel=='yes':
        # invers both digitized los_map and ref frame together
        # here we invert los_map = a*gacos + ramp
        d = losbins
        sigmad = losstd   
            
        G=np.zeros((len(d),6))
        G[:,0] = xbins**2
        G[:,1] = xbins
        G[:,2] = ybins**2
        G[:,3] = ybins
        G[:,4] = 1
        G[:,5] = modelbins

        x0 = lst.lstsq(G,d)[0]
        # print x0
        _func = lambda x: np.sum(((np.dot(G,x)-d)/sigmad)**2)
        _fprime = lambda x: 2*np.dot(G.T/sigmad, (np.dot(G,x)-d)/sigmad)
        pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
        a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
        print 'Remove ramp %f az**2, %f az  + %f r**2 + %f r + %f + %f model for date: %i'%(a,b,c,d,e,f,idates[l])

        #Compute ramp for los_map = f(model)
        funct = a*x**2 + b*x + c*y**2 + d*y + e
        functbins = a*xbins**2 + b*xbins + c*ybins**2 + d*ybins + e

        # build total G matrix
        G=np.zeros((len(los_map.flatten()),6))
        for i in xrange(nlign):
                G[i*ncol:(i+1)*ncol,0] = (i - line_beg)**2
                G[i*ncol:(i+1)*ncol,1] = i - line_beg
                G[i*ncol:(i+1)*ncol,2] = np.arange((ncol - col_beg))**2
                G[i*ncol:(i+1)*ncol,3] = np.arange(ncol - col_beg)
        G[:,4] = 1
        G[:,5] = model.flatten()

        # los_map = a*gacos + ramp
        remove = np.dot(G,pars).reshape(nlign,ncol)
        remove[model==0.] = 0.
        remove[np.isnan(los_map)] = np.float('NaN')

        # build total G matrix
        G=np.zeros((len(los_map.flatten()),5))
        for i in xrange(nlign):
            G[i*ncol:(i+1)*ncol,0] = (i - line_beg)**2
            G[i*ncol:(i+1)*ncol,1] = i - line_beg
            G[i*ncol:(i+1)*ncol,2] = np.arange((ncol - col_beg))**2
            G[i*ncol:(i+1)*ncol,3] = np.arange(ncol - col_beg)
        G[:,4] = 1

        # ramp only
        remove_ramp = np.dot(G,pars[:-1]).reshape(nlign,ncol)
        remove_ramp[model==0.] = 0.
        remove_ramp[np.isnan(los_map)] = np.float('NaN')

    # correction
    # los_map_flat = los_map - (gacos + ramp)
    los_map_flat[:,:] = los_map - remove
    los_map_flat[np.isnan(los_map)] = np.float('NaN')
    los_map_flat[los_map_flat>999.]= np.float('NaN')

    # model = gacos + ramp
    model_flat[:,:,l] = gacos[:,:,l] + remove_ramp

    # Refer los_map again (just to check)
    rms_ref = rms[indexref].flatten()
    amp_ref = 1./rms_ref
    amp_ref = amp_ref/np.nanmax(amp_ref)
    cst = np.nansum(los_map_flat[indexref]*amp_ref) / np.nansum(amp_ref)
    print 'Ref. area set to zero:', refstart,refend
    print 'Average phase within ref. area:', cst 
    los_map_flat = los_map_flat - cst
    los_map_flat[np.isnan(los_map)] = np.float('NaN')
    los_map_flat[los_map_flat>999.]= np.float('NaN')

    # compute variance flatten los_map
    var[l] = np.sqrt(np.nanmean(los_map_flat**2))
    print 'Var: ', var[l]

    # initiate figure depl
    fig = plt.figure(nfigure,figsize=(14,7))
    nfigure += 1
    ax = fig.add_subplot(3,2,1)
    im = ax.imshow(los_map,cmap=cmap,vmax=losmax,vmin=losmin)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.set_title('los_map {}'.format(idates[l]),fontsize=6)

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
    im = ax.imshow(los_map_flat,cmap=cmap,vmax=losmax,vmin=losmin)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.set_title('Correct los_map {}'.format(idates[l]),fontsize=6)

    ax = fig.add_subplot(3,2,5)
    g = np.linspace(np.nanmax(model_clean),np.nanmin(model_clean),100)
    ax.scatter(model_clean,los_clean - funct, s=0.005, alpha=0.1, rasterized=True)
    ax.plot(modelbins,losbins - functbins,'-r', lw =.5)
    ax.plot(g,f*g,'-r', lw =4.)
    ax.set_ylim([(los_clean-funct).min(),(los_clean-funct).max()])
    ax.set_xlim([model_clean.min(),model_clean.max()])
    ax.set_xlabel('GACOS ZTD')
    ax.set_ylabel('LOS delay')
    ax.set_title('los_map/Model')

    fig.tight_layout()
    fig.savefig('{}-gacos-cor.png'.format(idates[l]), format='PNG',dpi=150)

    if plot == 'yes':
        plt.show()
        # sys.exit()

    elif sformat == 'GAMMA':
        outfile = out_path + prefix + str(date1) + '_' + str(date2) + suffix +  suffout + rlook + '.unw' 
        fid = open(outfile, 'wb')
        los_map_flat.flatten().astype('>f4').tofile(fid)

    if plot=='yes':
        plt.show()

    plt.close('all')
    try:
        del dst_ds, ds, drv
    except:
        pass
    del los_map, rms_map

#####################################################################################
# INITIALISE 
#####################################################################################

# logging.basicConfig(level=logging.INFO,\
logging.basicConfig(level=logging.INFO,\
        format='%(asctime)s -- %(levelname)s -- %(message)s')
logger = logging.getLogger('correct_ifg_from_gacos.log')

# read int
date_1,date_2=np.loadtxt(int_list,comments="#",unpack=True,usecols=(0,1),dtype='i,i')
Nifg=len(date_1)
logger.info("number of interferogram: {}".format(Nifg))

# list dates
im = []; bt = []
for date1,date2 in zip(date_1,date_2):
    if date1 not in im: im.append(date1)
    if date2 not in im: im.append(date2)
nmax=len(im)
logger.info("number of image: {} ".format(nmax))
imd = date2dec(im)
cst = np.copy(imd[0])
# compute temporal baseline for TS inv
for i in range((nmax)):
    bt.append(imd[i]-cst)

#####################################################################################
# MAIN
#####################################################################################

# compute model IFGs
with TimeIt():
    # for kk in range(Nifg):
    work = range(Nifg)
    with poolcontext(processes=nproc) as pool:
        results = pool.map(gacos2ifg, work)

# apply corrections
with TimeIt():
    # for kk in range(Nifg):
    work = range(Nifg)
    with poolcontext(processes=nproc) as pool:
        results = pool.map(empirical_cor, work)
