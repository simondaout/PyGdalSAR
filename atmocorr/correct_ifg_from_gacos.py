#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT (Oxford)
#                 Nicholas Dodds
################################################################################

"""\
correct_ifg_from_gacos.py
-------------
Correct unwrapped IFGs from Gacos atmospheric models. 
Usage: 
    correct_ifg_from_gacos.py --int_list=<path> --inc_file=<path>  --ref_zone=<values> \
[--int_path=<path>] [--gacos_path=<path>]  [--prefix=<value>] \
    [--suffix=<value>] [--rlook=<value>] [--plot=<yes|no>] [--cohpixel=<yes/no>] [--threshold_coh=<value>] \
    [--format=<value>] [--mli_par=<path>] [--ramp=<cst|lin>] [--perc=<value>]\
    [--suffix_output=<value>] [--nproc=<nb_cores>] 

correct_ifg_from_gacos.py -h | --help
Options:
-h --help               Show this screen.
--int_list=<path>       Text file containing list of interferograms dates in two colums, $los_map1 $date2
--int_path=<path>       Relative path to input interferograms directory
--inc_file=<path>	Path to incidence file 
--gacos_path=<dir>      Path to GACOS atmospheric models [default: ./GACOS/]
--prefix=<value>        Prefix name $prefix$date1-$date2$suffix_$rlookrlks.unw [default: '']
--suffix=<value>        Suffix name $prefix$date1-$date2$suffix_$rlookrlks.unw [default: '']
--suffix_output=<value> Suffix output file name $prefix$date1-$date2$suffix$suffix_output [default:_gacos]
--rlook=<value>         look int. $prefix$date1-$date2$suffix_$rlookrlks.unw [default: 0]
--ref_zone=<lin_start,lin_end,col_start,refcol_end> Starting and ending lines and col numbers where phase is set to zero on the corrected data
--ramp=<cst|lin>        Estimate a constant or linear ramp in x and y [default: cst].
--cohpixel=<yes/no>     If Yes, use amplitude interferogram to weight and mask pixels (e.g Coherence, Colinearity, Amp Filter) [default: no]
--threshold_coh=<value> Thresold on rmspixel for ramp estimation [default: 2.]   
--plot=<yes|no>         Display results [default: yes]
--format VALUE          Format input files: ROI_PAC, GAMMA, GTIFF [default: GAMMA]
--mli_par= PATH         MLI file path to extract radar dimensions.
--perc=<value>          Percentile of hidden LOS pixel for the estimation and clean outliers [default:98.]
--nproc=<values>            number of processor (default: 1)
"""

import gdal
gdal.UseExceptions()

import sys,os,logging
import numpy as np
import scipy.optimize as opt
import scipy.linalg as lst
import scipy.constants as const
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

if arguments["--int_list"] != None:
  int_list=arguments["--int_list"]
else:
  raise Exception('Input Missing: No input unwrapped IFG list.')
if arguments["--inc_file"] != None:
  inc_file=arguments["--inc_file"]
else:
  raise Exception('Input Missing: No input Incidence angle file.')
if arguments["--ref_zone"] != None:
  ref = list(map(int,arguments["--ref_zone"].replace(',',' ').split()))
  refline_beg,refline_end, refcol_beg, refcol_end = ref[0], ref[1], ref[2], ref[3]
else:
  raise Exception('Input Missing: No reference zone defined.')

if arguments["--int_path"] == None:
    int_path='.'
else:
    int_path=arguments["--int_path"] + '/'
if arguments["--gacos_path"] ==  None:
   gacos_path = "./gacos/"
else:
   gacos_path = arguments["--gacos_path"] + '/'
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
    print('Plot is yes: Set nproc to 1')
    nproc = 1
    if environ["TERM"].startswith("screen"):
        matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt
else:
    plot = 'no'
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt
if arguments["--ramp"] ==  None:
    ramp = 'cst'
else:
    ramp = arguments["--ramp"]
if arguments["--format"] ==  None:
    sformat = 'GAMMA'
else:
    sformat = arguments["--format"]
if sformat == 'GAMMA':
  if arguments["--mli_par"] != None:
    mli_par = arguments["--mli_par"]
  else:
    raise Exception('Input Missing: No input MLI par file given with GAMMA format.')

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
if arguments["--suffix_output"] ==  None:
    suffout = '_gacos'
else:
    suffout = arguments["--suffix_output"]


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

def pygrep(infile, rstring):
    for line in open(infile):
        if rstring in line:
            return line.split("\n")[0]

#####################################################################################
# FUNCTIONS
#####################################################################################


def correct_unw_gacos(kk):
    """
    Function that compute differential gacos model and correct ifgs 
    """ 

    date1, date2 = date_1[kk], date_2[kk]
    idate = str(date1) + '-' + str(date2) 

    if sformat == 'GAMMA':
        infile1 = gacos_path +  str(date1) + '_crop.ztd.rdr.unw'
        infile2 = gacos_path +  str(date2) + '_crop.ztd.rdr.unw'
        checkinfile(infile1); checkinfile(infile2)
        gacos1 = np.fromfile(infile1, dtype='>f4').reshape(lines,cols)
        gacos2 = np.fromfile(infile2, dtype='>f4').reshape(lines,cols)
        infile= int_path + prefix + str(date1) + '_' + str(date2) + suffix + rlook + '.unw'
        checkinfile(infile)
        los_map = np.fromfile(infile,dtype='>f4').reshape(lines,cols)

    logger.info('lines:{0}, cols:{1}, IFG:{2}:'.format(lines, cols, idate))
    
    # Time differential of GACOS data
    g_ztd_data = np.subtract(gacos1, gacos2) 
    # Convert GACOS ZTD to LOS
    g_los_data = np.divide(g_ztd_data, cos_los)
    # Convert GACOS LOS (in m) to radians
    model = g_los_data / rad2m
    
    # save gacos corrections
    if sformat == 'GAMMA':
      gacosf = unw_gacos + '/' + str(date1) + '-' + str(date2) + '_gacosmdel' + '.unw'
      fid = open(gacosf,'wb')
      model.flatten().astype('>f4').tofile(fid)
      fid.close()

    # load coherence or whatever
    rms_map = np.ones((lines,cols))

    if rmsf=='yes':
      if sformat == 'GAMMA':
        rmsfile=  int_path + str(date1) + '_' + str(date2) + '.filt.cc'
        rms_map = np.fromfile(rmsfile,dtype='>f4').reshape(lines,cols)

    # extract range and azimuth coordinates from ref or radar file
    pix_az, pix_rg = np.indices((lines,cols))

    _los_map = np.copy(los_map)
    _los_map[np.logical_or(los_map==0,los_map>=9990)] = np.float('NaN')
    losmin,losmax = np.nanpercentile(_los_map,100-perc),np.nanpercentile(_los_map,perc)
    gacosmin,gacosmax = np.nanpercentile(model,100-perc),np.nanpercentile(model,perc)
    del _los_map

    funct = 0.
    pix_lin, pix_col = np.indices((lines,cols))

    # find proportionality between los_map and model
    index = np.nonzero(
        np.logical_and(~np.isnan(los_map),
        np.logical_and(los_map<losmax,
        np.logical_and(los_map>losmin,
        np.logical_and(rms_map<threshold_rms,
        np.logical_and(rms_map>1.e-6,
        np.logical_and(model<gacosmax,
        np.logical_and(model>gacosmin,
        np.logical_and(los_map!=0.0,
        model!=0.0,
        )))))))))

    indexref = np.nonzero(
        np.logical_and(~np.isnan(los_map),
        np.logical_and(pix_lin>refline_beg,
        np.logical_and(pix_lin<refline_end,
        np.logical_and(los_map<losmax,
        np.logical_and(los_map>losmin,
        np.logical_and(rms_map<threshold_rms,
        np.logical_and(rms_map>1.e-6,
        np.logical_and(model<gacosmax,
        np.logical_and(model>gacosmin,
        np.logical_and(los_map!=0.0,
        np.logical_and(model!=0.0,
        np.logical_and(pix_col>refcol_beg, pix_col<refcol_end
        )))))))))))))
    
    # extract los for empirical estimation
    temp = np.array(index).T
    x_clean = temp[:,0]; y_clean = temp[:,1]
    los_clean = los_map[index].flatten()
    model_clean = model[index].flatten()
    rms_clean = rms_map[index].flatten()
    del temp

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
    los_clean2, model_clean2, x_clean2, y_clean2, rms_clean2 = [], [], [], [], []
    for j in range(len(bins)-1):
            uu = np.flatnonzero(inds == j)
            if len(uu)>500:
                modelbins.append(bins[j] + (bins[j+1] - bins[j])/2.)
                
                # do a clean within the sliding median
                indice = np.flatnonzero(np.logical_and(los_clean[uu]>np.percentile(\
                    los_clean[uu],10.),los_clean[uu]<np.percentile(los_clean[uu],90.)))

                losstd.append(np.std(los_clean[uu][indice]))
                losbins.append(np.median(los_clean[uu][indice]))
                xbins.append(np.median(x_clean[uu][indice]))
                ybins.append(np.median(y_clean[uu][indice]))

                # remove outliers from data
                los_clean2.append(los_clean[uu][indice])
                x_clean2.append(x_clean[uu][indice])
                y_clean2.append(y_clean[uu][indice])
                model_clean2.append(model_clean[uu][indice])
                # new rms is clean rms time the standard deviation of the los within the bin
                rms_clean2.append(rms_clean[uu][indice]*np.nanstd(los_clean[uu][indice]))

    losbins = np.array(losbins)
    losstd = np.array(losstd)
    modelbins = np.array(modelbins)
    xbins, ybins = np.array(xbins),np.array(ybins)

    los_clean = np.concatenate(los_clean2)
    x_clean, y_clean = np.concatenate(x_clean2), np.concatenate(y_clean2)
    model_clean = np.concatenate(model_clean2)
    rms_clean = np.concatenate(rms_clean2)
    del los_clean2, x_clean2, y_clean2, model_clean2, bins, inds
    
    if (ramp == 'cst' ):
        # adjust ifg and gacos model with a simple constante estimated within the ref area
        a = cst
        logger.info('Remove cst: %f for date: %s'%(a,idate))
       
        remove_ramp = np.ones((lines,cols))*cst
        remove_ramp[model==0.] = 0.
        remove_ramp[np.isnan(los_map)] = np.float('NaN')
        
        # Compute ramp for los_map = f(model) 
        functbins = a
        funct = a

    elif (ramp == 'lin' ):
        # here we want to minimize los_map-gacos = ramp
        # invers both clean data and ref frame together

        d = d = los_clean - model_clean
        sigmad = rms_clean

        G=np.zeros((len(d),3))
        G[:,0] = x_clean
        G[:,1] = y_clean
        G[:,2] = 1

        x0 = lst.lstsq(G,d)[0]
        # print x0
        _func = lambda x: np.sum(((np.dot(G,x)-d)/sigmad)**2)
        _fprime = lambda x: 2*np.dot(G.T/sigmad, (np.dot(G,x)-d)/sigmad)
        pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=500,full_output=True,iprint=0)[0]
        a = pars[0]; b = pars[1]; c = pars[2]
        print ('Remove ramp  %f az  + %f r + %f for date: %i'%(a,b,c,idate))
            
        # Compute ramp for los_map = f(model)
        functbins = a*xbins + b*ybins + c
        funct = a*x_clean + b*y_clean + c

        # build total G matrix
        G=np.zeros((len(los_map.flatten()),3))
        for i in xrange(nlign):
            G[i*ncol:(i+1)*ncol,0] = i  
            G[i*ncol:(i+1)*ncol,1] = np.arange(ncol) 
        G[:,2] = 1

        # compute ramp
        remove_ramp = np.dot(G,pars).reshape(lines,cols)
        remove_ramp[model==0.] = 0.
        remove_ramp[np.isnan(los_map)] = np.float('NaN')

    # Aply correction

    # los_map_flat = los_map - (gacos + ramp)
    los_map_flat = los_map - remove_ramp
    los_map_flat[np.isnan(los_map)] = np.float('NaN')
    los_map_flat[los_map_flat>999.]= np.float('NaN')
    # model = gacos + ramp
    model_flat = model + remove_ramp

    # Refer los_map again (just to check)
    rms_ref = rms_map[indexref].flatten()
    amp_ref = 1./rms_ref
    amp_ref = amp_ref/np.nanmax(amp_ref)
    cst = np.nansum(los_map_flat[indexref]*amp_ref) / np.nansum(amp_ref)
    logger.info('Ref area set to zero: lines: {0}-{1}, cols: {2}-{3}'.format(refline_beg,refline_end,refcol_beg,refcol_end))
    logger.info('Average phase within ref area: {0}'.format(cst))
    los_map_flat = los_map_flat - cst
    los_map_flat[np.isnan(los_map)] = np.float('NaN')
    los_map_flat[los_map_flat>999.]= np.float('NaN')

    # compute variance flatten los_map
    var = np.nanstd(los_map_flat)
    logger.info('Var: {0} '.format(var))    

    # initiate figure depl
    fig = plt.figure(0,figsize=(14,7))
    ax = fig.add_subplot(2,2,1)
    im = ax.imshow(los_map,cmap=cmap,vmax=losmax,vmin=losmin)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.set_title('los_map {}'.format(idate),fontsize=6)

    ax = fig.add_subplot(2,2,2)
    im = ax.imshow(model,cmap=cmap)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.set_title('Model {}'.format(idate),fontsize=6)
    
    ax = fig.add_subplot(2,2,3)
    im = ax.imshow(model_flat,cmap=cmap,vmax=losmax,vmin=losmin)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.set_title('Flatten Model {}'.format(idate),fontsize=6)

    ax = fig.add_subplot(2,2,4)
    im = ax.imshow(los_map_flat,cmap=cmap,vmax=losmax,vmin=losmin)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.set_title('Correct los_map {}'.format(idate),fontsize=6)
    fig.tight_layout()
    fig.savefig(unw_gacos+'/{}-data-model-maps.png'.format(idate), format='PNG',dpi=150)

    fig = plt.figure(2,figsize=(7,5))
    ax = fig.add_subplot(1,1,1)
    g = np.linspace(np.nanmax(model_clean),np.nanmin(model_clean),100)
    ax.scatter(model_clean,los_clean - funct, s=0.005, alpha=0.1, rasterized=True)
    ax.plot(modelbins,losbins - functbins,'-r', lw =.5)
    ax.set_xlim([model_clean.min(),model_clean.max()])
    ax.set_xlabel('GACOS ZTD')
    ax.set_ylabel('LOS delay - RAMP')
    ax.set_title('Data/Model')
    plt.legend(loc='best')
    fig.tight_layout()
    fig.savefig(unw_gacos+'/{}-gacos-cor.png'.format(idate), format='PNG',dpi=150)

    if plot == 'yes':
        plt.show()
        # sys.exit()

    plt.close('all')

    if sformat == 'GAMMA':
        outfile = unw_gacos + prefix + str(date1) + '_' + str(date2) + suffix +  suffout + rlook + '.unw' 
        fid = open(outfile, 'wb')
        los_map_flat.flatten().astype('>f4').tofile(fid)

    del los_map, rms_map, los_map_flat, model_flat

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

if sformat == 'GAMMA':
    # Extract radar dimensions from mli file
    checkinfile(mli_par)
    if os.path.isfile(mli_par) == True:
      cols = int(pygrep(mli_par, "range_samples").split(" ")[-1])
      lines = int(pygrep(mli_par, "azimuth_lines").split(" ")[-1])
      rad_freq = float(pygrep(mli_par, "radar_frequency").split(" ")[-3])
    else:
      logger.critical('MLI par (input) file does not exist.')
    
    # Read incidence file
    inc = np.fromfile(inc_file, dtype='>f4').reshape(lines,cols)
    inc_d = np.deg2rad(inc)
    cos_los = np.cos(inc_d)
else:
    logger.critical('Other formats not implemented yet... Exit!')

# gacos2los conversion factor
rad_lambda = const.c/rad_freq
# rad2mm = 56mm / (4*math.pi) = 4.4563384065730 
rad2m = rad_lambda/(4*const.pi)

# Saving in new dir
unw_gacos = 'unw_gacos/'
if not os.path.exists(unw_gacos):
    os.mkdir(unw_gacos)

# Choose plot option
cmap = cm.gist_rainbow_r
# cmap = cm.jet
cmap.set_bad('white')

fig = plt.figure(1,figsize=(8,4))
# Plotting the incidence file
ax = fig.add_subplot(1,1,1)
cax = ax.imshow(inc, cmap=cmap,  interpolation='bicubic', extent=None)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)
ax.set_title('TOPS incidence (degrees)')
fig.savefig(unw_gacos+'/TOPS_inc_DEM.png', format='PNG', bbox_inches='tight')

#####################################################################################
# MAIN
#####################################################################################

with TimeIt():
    # for kk in range(Nifg):
    work = range(Nifg)
    with poolcontext(processes=nproc) as pool:
        results = pool.map(correct_unw_gacos, work)

