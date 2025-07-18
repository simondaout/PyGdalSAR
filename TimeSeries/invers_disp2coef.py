#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Simon DAOUT (CRPG, Nancy)
############################################

"""\
invers_disp2coef.py
-------------
Spatial and temporal inversions of the cumulative time series delay maps based on an iteration procedure. 
Requieres depl_cumule (BIP format) and images_retenues, output of invers_pixel.
At each iteration, (1) estimation of spatial ramps and/or phase/topography relationship, (2) linear decomposition
in time based on a library of temporal functions (linear, heaviside, logarithm, seasonal, slow-slip...) 
or given in as vector text file with the vector argument, (3) estimation of the misfit that will be then used as 
weight for the next iteration. 

Usage: invers_disp2coef.py  [--cube=<path>] [--lectfile=<path>] [--list_images=<path>] [--aps=<path>] [--rms=<path>]\
[--rmspixel=<path>] [--threshold_rms=<value>] [--ref_zone=<jstart,jend,istart,iend>] [--niter=<value>]  [--spatialiter=<yes/no>] \
[--linear=<yes/no>] [--steps=<value,value>] [--postseismic=<value,value>] [--seasonal=<yes/no>] [--seasonal_increase=<yes/no>] [--slowslip=<value,value>] \
[--semianual=<yes/no>] [--bianual=<yes/no>]  [--bperp=<yes/no>] [--vector=<path>] \
[--flat=<0/1/2/3/4/5/6/7/8/9>] [--nfit=<0/1>] [--ivar=<0/1>] \
[--sampling=<value>] [--emp_sampling=<value>] [--imref=<value>]  [--cond=<value>] [--ineq=<yes/no>]  \
[--mask=<path>] [--rampmask=<yes/no>] [--threshold_mask=<value>] [--scale_mask=<value>] [--tempmask=<yes/no>]\
[--topofile=<path>] [--aspect=<path>] [--perc_topo=<value>] [--perc_los=<value>] \
[--crop=<value,value,value,value>] [--crop_emp=<value,value,value,value>] [--fulloutput=<yes/no>] [--geotiff=<path>] [--plot=<yes/no>] \
[--dateslim=<values_min,value_max>]  [--nproc=<nb_cores>] [--ndatasets=<nb_data_sets>] 

-h --help               Show this screen
--cube=<path>           Path to time series displacements cube file [default: no]
--lectfile=<path>       Path to the lect.in file. Simple text file containing width and length and number of images of the time series cube (output of invers_pixel). By default the program will try to find an .hdr file. [default: lect.in].
--list_images=<path>    Path to list images file. text file made of 5 columns containing for each images 1) number 2) Doppler freq (not read) 3) date in YYYYMMDD format 4) numerical date 5) perpendicular baseline [default: images_retenues].
--aps=<path>            Path to the APS file. Text file with two columns (dates, APS) giving an input APS to each dates. By default, no weigthing if no empirical estimation or misfit from the spatial estimation used as input uncertianties [default: None].
--rms=<path>            Path to the RMS file. Text file with two columns (dates, RMS) giving an input RMS to each dates. By default, no weigthing if no empirical estimation or misfit from the spatial estimation used as input uncertianties [default: None].
--rmspixel=<path>       Path to the RMS map. Map in r4 or tiff format that gives an error for each pixel (e.g RMSpixel, output of invers_pixel) [default: None].
--threshold_rms=<value> Threshold on rmspixel argument for empitical spatial estimations [default: 1.]
--ref_zone=<lin_start,lin_end,col_start,col_end> Starting and ending lines and col numbers where phase is set to zero [default: None]  
--niter=<value>         Number of iterations. At the first iteration, image uncertainties is given by aps file or misfit spatial iteration, while for the next itarations, uncertainties are equals to the global RMS previous temporal decomposition [default: 0].
--linear=<yes/no>       Add a linear function in the inversion [default:yes]
--steps=<value,value>     Add heaviside functions to the inversion, indicate steps time (e.g 2004.,2006.) [default: None]
--postseismic=<value,value>   Add logarithmic transients to each steps step. Indicate characteristic time of the log function, must be a serie of values of the same lenght than steps (e.g 1.,1.). To not associate postseismic function to a give steps step, put None (e.g None,1.) [default: None].
--slowslip=<value,value>      Add slow-slip function in the inversion . As defined by Larson et al., 2004. Indicate median and characteristic time of the events (e.g. 2004.,1,2006,0.5) [default: None].
--vector=<path>         Path to the vector text files containing a value for each dates [default: None]
--seasonal=<yes/no>       If yes, add seasonal terms in the decomposition [default: no]
--semianual=<yes/no>      If yes, add semianual terms in the decomposition [default: no]
--seasonal_increase PATH         If yes, add seasonal terms function of time in the inversion
--bianual=<yes/no>       If yes, add bianual terms in the decomposition [default: no]
--bperp=<yes/no>           If yes, add term proportional to the perpendicular baseline in the inversion [default: no]
--ivar=<0/1>            Define the phase/elevation relationship: ivar=0 function of elevation, ivar=1 crossed function of azimuth and elevation [default: 0]
--nfit=<0/1>            Fit degree in azimuth or in elevation (0:linear (default), 1: quadratic) [default: 0]
--flat=<0/1/2/3/4/5/6/7/8/9>             Remove a spatial ramp at each iteration [default: 0].
--spatialiter=<yes/no>   If 'yes' iterate the spatial estimations at each iterations (defined by niter arguments) on the maps minus the temporal terms (ie. linear, steps...) [default: yes]
--sampling=<value>      Downsampling factor for temporal decomposition [default: 1]
--emp_sampling=<value>      Downsampling factor for empirical estimations [default: 1]
--imref=<value>         Reference image number [default: 1]
--mask=<path>           Path to mask file in r4 or tif format for the empirical spatial estimations. Keep only values > threshold_mask for ramp estimation [default: no].
--rampmask=<yes/no>     Remove a quadratic ramp in range and linear ramp in azimuth on the mask [default: no].
--threshold_mask=<value> Threshold on mask: take only < values (use scale factor for convenience) [default: 1].
--scale_mask=<value>     Scale factor to apply on mask [default: 1]
--tempmask=<yes/no>       If yes, also use the mask for the temporal decomposition [default: no]
--topofile=<path>         Path to topographic file in r4 or tif format. If not None, add a phase-elevation relationship in the saptial estimation [default: None].
--aspect=<path>            Path to aspect file in r4 or tif format: take into account the slope orientation in  the phase/topo relationship .
--perc_los=<value>        Percentile of hidden LOS pixel for the spatial estimations to clean outliers [default:99.]
--perc_topo=<value>       Percentile of topography ranges for the spatial estimations to remove some very low valleys or peaks [default:99.]
--cond=<value>            Condition value for optimization: Singular value smaller than cond are considered zero [default: 1e-3]
--ineq=<yes/no>           If yes, sequential least-square optimisation. If no, SVD inversion with mask on eigenvalues smaller than --cond value. If postseimsic functions, add ineguality constraints in the inversion. Use least square results without post-seismic functions as a first guess to iterate the inversion. Then, force postseismic to be the same sign and inferior than steps steps of the first guess [default: yes].
--fulloutput=<yes/no>      If yes produce maps of models, residuals, ramps, as well as flatten cube without seasonal and linear term [default: no]
--geotiff=<path>           Path to Geotiff to save outputs in tif format. If None save output are saved as .r4 files 
--plot=<yes/no>         Display plots [default: no]
--dateslim=<value,value>     Datemin,Datemax time series 
--crop=<value,value,value,value>            Define a region of interest for the temporal decomposition 
--crop_emp=<value,value,value,value>    Define a region of interest for the spatial estimatiom (ramp+phase/topo) 
--nproc=<nb_cores>        Use <nb_cores> local cores to create delay maps [Default: 4]
--ndatasets=<nb_data_sets> Number of data sets [Default:1] 
"""

# 0: ref frame [default], 1: range ramp ax+b , 2: azimutal ramp ay+b, 3: ax+by+c,
# 4: ax+by+cxy+d 5: ax**2+bx+cy+d, 6: ay**2+by+cx+d, 7: ay**2+by+cx**2+dx+e,
# 8: ay**2+by+cx**3+dx**2+ex+f, 9: ax+by+cxy**2+dxy+e

print()
print('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #')
print('#'                                                                 '#')
print('#         Linear Inversion of InSAR time series displacements       #')
print('#         with a decomposition in time                              #')
print('#'                                                                 '#')
print('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #')
print()
print()
print('Author: Simon Daout')
print('Please cite:') 
print ('Daout, S., Doin, M. P., Peltzer, G., Socquet, A., & Lasserre, C. (2017). Large‐scale InSAR monitoring of permafrost freeze‐thaw cycles on the Tibetan Plateau. Geophysical Research Letters, 44(2), 901-909.')
print('Daout, S., Sudhaus, H., Kausch, T., Steinberg, A., & Dini, B. (2019). Interseismic and postseismic shallow creep of the North Qaidam Thrust faults detected with a multitemporal InSAR analysis. Journal of Geophysical Research: Solid Earth, 124(7), 7259-7279.')
print()
print()

import numpy as np
from numpy.lib.stride_tricks import as_strided
import scipy as sp
import scipy.optimize as opt
import numpy.linalg as lst
from osgeo import gdal, osr
import math, sys, getopt, shutil
from os import path, environ, getcwd
import os
import matplotlib
import matplotlib.cm as cm
import matplotlib.dates as mdates
from datetime import datetime as dt

try:
    from nsbas import docopt
except:
    import docopt

from contextlib import contextmanager
from functools import wraps, partial
# import multiprocessing
import logging, time

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning) 

logging.basicConfig(level=logging.INFO,\
        format='line %(lineno)s -- %(levelname)s -- %(message)s')
logger = logging.getLogger('invers_disp2coef.log')
start_time = time.time()

################################
# Create lib of wavelet functions
################################

class pattern:
    def __init__(self,name,reduction,date):
        self.name=name
        self.reduction=reduction
        self.date=date

    def info(self):
        print(self.name, self.date)

### BASIS FUNCTIONS: function of time

def Heaviside(t):
        h=np.zeros((len(t)))
        h[t>=0]=1.0
        return h

def Box(t):
        return Heaviside(t+0.5)-Heaviside(t-0.5)

class steps(pattern):
      def __init__(self,name,reduction,date):
          pattern.__init__(self,name,reduction,date)
          self.to=date

      def g(self,t):
        return Heaviside(t-self.to)

class postseismic(pattern):
      def __init__(self,name,reduction,date,tcar=1):
          pattern.__init__(self,name,reduction,date)
          self.to=date
          self.tcar=tcar

      def g(self,t):
        t=(t-self.to)/self.tcar
        t[t<=0] = 0
        g = np.log10(1+t)
        return g

class reference(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
    def g(self,t):
        return np.ones((t.size))

class linear(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
        self.to=date

    def g(self,t):
        func=(t-self.to)
        return func

class sint(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
        self.to=date

    def g(self,t):
        func=np.zeros(t.size)
        for i in range(t.size):
            func[i]=(t[i]-self.to)*math.sin(2*math.pi*(t[i]-self.to))
        return func

class cost(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
        self.to=date

    def g(self,t):
        func=np.zeros(t.size)
        for i in range(t.size):
            func[i]=(t[i]-self.to)*math.cos(2*math.pi*(t[i]-self.to))
        return func

class sin2var(pattern):
     def __init__(self,name,reduction,date):
         pattern.__init__(self,name,reduction,date)
         self.to=date

     def g(self,t):
         func=np.zeros(t.size)
         for i in range(t.size):
             func[i]=math.sin(4*math.pi*(t[i]-self.to))
         return func

class cos2var(pattern):
     def __init__(self,name,reduction,date):
         pattern.__init__(self,name,reduction,date)
         self.to=date

     def g(self,t):
         func=np.zeros(t.size)
         for i in range(t.size):
             func[i]=math.cos(4*math.pi*(t[i]-self.to))
         return func

class sin5var(pattern):
     def __init__(self,name,reduction,date):
         pattern.__init__(self,name,reduction,date)
         self.to=date

     def g(self,t):
         func=np.zeros(t.size)
         for i in range(t.size):
             func[i]=math.sin(math.pi*(t[i]-self.to))
         return func

class cos5var(pattern):
     def __init__(self,name,reduction,date):
         pattern.__init__(self,name,reduction,date)
         self.to=date

     def g(self,t):
         func=np.zeros(t.size)
         for i in range(t.size):
             func[i]=math.cos(math.pi*(t[i]-self.to))
         return func

class sinvar(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
        self.to=date

    def g(self,t):
        func=np.zeros(t.size)
        for i in range(t.size):
            func[i]=math.sin(2*math.pi*(t[i]-self.to))
        return func

class cosvar(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
        self.to=date

    def g(self,t):
        func=np.zeros(t.size)
        for i in range(t.size):
            func[i]=math.cos(2*math.pi*(t[i]-self.to))
        return func

class slowslip(pattern):
      def __init__(self,name,reduction,date,tcar=1):
          pattern.__init__(self,name,reduction,date)
          self.to=date
          self.tcar=tcar

      def g(self,t):
          t=(t-self.to)/self.tcar
          funct = 0.5*(np.tanh(t)-1) + 1
          return funct

### KERNEL FUNCTIONS: not function of time
class corrBperp(pattern):
    def __init__(self,name,reduction,bp0,bp):
        self.name = name
        self.reduction = reduction
        self.bpo=bp0
        self.bp=bp

    def info(self):
        print(self.name)

    def g(self,index):
        func = (self.bp-self.bpo)
        return func[index]

class vector(pattern):
    def __init__(self,name,reduction,vect):
        self.name = name
        self.reduction = reduction
        self.func=vect

    def info(self):
        print(self.name)

    def g(self,index):
        return self.func[index]

def date2dec(dates):
    dates  = np.atleast_1d(dates)
    times = []
    for date in dates:
        x = dt.strptime('{}'.format(date),'%Y%m%d')
        dec = float(x.strftime('%j'))/365.1
        year = float(x.strftime('%Y'))
        times.append(year + dec)
    return times


##################################################################################
###  Extras functions and context maganers
##################################################################################

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
        self.start = dt.now()
        logger.info('Starting time process: {0}'.format(self.start))
    def __exit__(self, type, value, traceback):
        logger.info('Time process: {0}s'.format((dt.now() - self.start).total_seconds()))


def checkinfile(file):
    if path.exists(file) is False:
        logger.critical("File: {0} not found, Exit !".format(file))
        logger.info("File: {0} not found in {1}, Exit !".format(file,getcwd()))
        sys.exit()

def plot_displacement_maps(maps, idates, nfigure=0, cmap='RdBu', plot='yes', filename='maps.eps', title='Time series maps', fig_dpi=150):
    """
    Affiche et sauvegarde une série de cartes de déplacement (type cube t,x,y).

    """
    N = maps.shape[2]
    ncols = int(N / 4) + 1
    vmax = np.nanpercentile(maps, 99.)
    vmin = np.nanpercentile(maps, 1.)

    fig = plt.figure(nfigure, figsize=(14, 10))
    fig.subplots_adjust(wspace=0.001)

    for l in range(N):
        d = as_strided(maps[:, :, l])
        ax = fig.add_subplot(4, ncols, l + 1)
        cax = ax.imshow(d, cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
        ax.set_title(str(idates[l]), fontsize=6)
        ax.set_xticks([])
        ax.set_yticks([])

    plt.suptitle('{}'.format(title), fontsize=12)
    fig.colorbar(cax, orientation='vertical', aspect=10)
    fig.subplots_adjust(hspace=.001,wspace=0.001)
    fig.savefig(filename, format='EPS', dpi=fig_dpi)

    if plot == 'yes':
        plt.show()
    
    plt.close(fig)
    del fig


################################
# Initialization
################################

# read arguments
arguments = docopt.docopt(__doc__)
if arguments["--lectfile"] ==  None:
    arguments["--lectfile"] = "lect.in"
if arguments["--list_images"] ==  None:
    try:
        checkinfile("images_retenues")
        arguments["--list_images"] = "images_retenues"
    except:
        checkinfile("list_images.txt")
        arguments["--list_images"] = "list_images.txt"
if arguments["--cube"] ==  None:
    arguments["--cube"] = "depl_cumule"
if arguments["--linear"] ==  None:
    arguments["--linear"] = 'yes'
if arguments["--seasonal"] ==  None:
    arguments["--seasonal"] = 'no'
if arguments["--seasonal_increase"] ==  None:
    arguments["--seasonal_increase"] = 'no'
if arguments["--semianual"] ==  None:
    arguments["--semianual"] = 'no'
if arguments["--bianual"] ==  None:
    arguments["--bianual"] = 'no'
if arguments["--bperp"] ==  None:
    arguments["--bperp"] = 'no'
if arguments["--niter"] ==  None:
    arguments["--niter"] = 1
if arguments["--spatialiter"] ==  None:
    arguments["--spatialiter"] = 'yes'
if arguments["--flat"] == None:
    flat = 0
elif int(arguments["--flat"]) <  10:
    flat = int(arguments["--flat"])
else:
    flat = 0
if arguments["--sampling"] ==  None:
    arguments["--sampling"] = 1
if arguments["--emp_sampling"] ==  None:
    arguments["--emp_sampling"] = 1
if arguments["--mask"] ==  None:
    arguments["--mask"] = None
if arguments["--rampmask"] ==  None:
    arguments["--rampmask"] = 'no'
if arguments["--threshold_mask"] ==  None:
    arguments["--threshold_mask"] = 0.
if arguments["--threshold_rms"] ==  None:
    arguments["--threshold_rms"]  = 1.
if arguments["--tempmask"] ==  None:
    arguments["--tempmask"] = 'no'
if arguments["--scale_mask"] ==  None:
    arguments["--scale_mask"] = 1
if arguments["--topofile"] ==  None:
   arguments["--topofile"] = None
if arguments["--cond"] ==  None:
    arguments["--cond"] = 1e-3
if arguments["--rmspixel"] ==  None:
    arguments["--rmspixel"] = None
if arguments["--ineq"] ==  None:
    arguments["--ineq"] = 'no'
if arguments["--fulloutput"] ==  None:
    arguments["--fulloutput"] = 'no'
if arguments["--geotiff"] is not None:
    logger.warning('Load geotiff: {}'.format(arguments["--geotiff"]))
    georef = gdal.Open(arguments["--geotiff"])
    gt = georef.GetGeoTransform()
    proj = georef.GetProjection()
    driver = gdal.GetDriverByName('GTiff')
    logger.warning('Set geotiff projection: {}'.format(proj))
if arguments["--ivar"] == None:
    ivar = 0
elif int(arguments["--ivar"]) <  2:
    ivar = int(arguments["--ivar"])
    if arguments["--topofile"] is None:
        logger.critical('No topographic file is given. Empirical phase/topo will not be performed')
else:
    logger.warning('Error: ivar > 1, set ivar to 0')
    if arguments["--topofile"] is None:
        logger.critical('No topographic file is given. Empirical phase/topo will not be performed')
    ivar = 0
if arguments["--nfit"] == None:
    nfit = 0
elif int(arguments["--nfit"]) <  2:
    nfit = int(arguments["--nfit"])
    if arguments["--topofile"] is None:
        logger.critical('No topographic file is given. Empirical phase/topo will not be performed')
else:
    logger.warning('Error: nfit > 1, set nfit to 0')
    nfit = 0
    if arguments["--topofile"] is None:
        logger.critical('No topographic file is given. Empirical phase/topo will not be performed')
if arguments["--perc_topo"] ==  None:
    arguments["--perc_topo"] = 90.
if arguments["--perc_los"] ==  None:
    arguments["--perc_los"] = 98.
if arguments["--nproc"] ==  None:
    nproc = 5
else:
    nproc = int(arguments["--nproc"])
if arguments["--ndatasets"] ==  None:
    ndata = 1
else:
    ndata = int(arguments["--ndatasets"])
if arguments["--plot"] ==  'yes':
    plot = 'yes'
    logger.warning('plot is yes. Set nproc to 1')
    nproc = 1
    if environ["TERM"].startswith("screen"):
        matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt
    from pylab import date2num
else:
    plot = 'no'
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt
    from pylab import date2num
if arguments["--imref"] ==  None:
    imref = 0
elif int(arguments["--imref"]) < 1:
    logger.warning('--imref must be between 1 and Nimages')
else:
    imref = int(arguments["--imref"]) - 1

#####################################################################################
# INITIALISATION
#####################################################################################

# cm
cmap = cm.jet
cmap.set_bad('white')

logger.debug('Load list of dates file: {}'.format(arguments["--list_images"]))
checkinfile(arguments["--list_images"])

# Initialisation des variables
idates = None
dates = None
bases = []
baserefs = []

# Chargement du fichier images_retenues avec une structure dynamique pour ndata > 2
if ndata >= 1:
    # Charger les colonnes pour idates et dates (colonnes 1 et 3)
    # Générer dynamiquement les colonnes pour chaque base (à partir de la colonne 5)
    columns = [1, 3] + list(range(5, 5 + ndata))

    # Charger les données en utilisant `usecols=columns`
    data = np.loadtxt(arguments["--list_images"], comments='#', usecols=columns, unpack=True)

    # Extraire les données de dates
    idates = data[0].astype(int)
    dates = data[1]

    # Extraire les bases et initialiser les références
    for i in range(ndata):
        base = data[2 + i]
        bases.append(base)

        # Choisir la première date comme référence (index `imref = 0`)
        baserefs.append(base[0])  # baseref correspondant à chaque base

# Nombre total d'images (dates)
N = len(dates)

# datemin est la reference temporelle pour les fonctions de base
# dmin, dmax sont les dates limites des figures
if arguments["--dateslim"] is not  None:
    dmin,dmax = arguments["--dateslim"].replace(',',' ').split()
    datemin = int(date2dec(dmin)[0])
    datemax = int(date2dec(dmax)[0]) + 1
    dmin = str(datemin) + '0101'
    dmax = str(datemax) + '0101'
else:
    datemin, datemax = int(np.min(dates)), int(np.max(dates)) + 1
    dmin = str(int(np.min(dates))) + '0101'
    dmax = str(int(np.max(dates))+1) + '0101'

# clean dates
indexd = np.flatnonzero(np.logical_and(dates<=datemax,dates>=datemin))
idates,dates = idates[indexd],dates[indexd]
# Mettre à jour les bases en filtrant également selon indexd
bases = [bp[indexd] for bp in bases]

# lect cube
checkinfile(arguments["--cube"])
ds = gdal.Open(arguments["--cube"])
if not ds:
  logger.info('.hdr file time series cube {0}, not found, open {1}'.format(arguments["--cube"],arguments["--lectfile"]))
  ncol, nlines = list(map(int, open(arguments["--lectfile"]).readline().split(None, 2)[0:2]))
else:
  ncol, nlines = ds.RasterXSize, ds.RasterYSize
  N = ds.RasterCount

logger.debug('Read reference zones: {}'.format(arguments["--ref_zone"]))
if arguments["--ref_zone"] == None:
    lin_start, lin_end, col_start, col_end = None,None,None,None
else:
    ref = list(map(int,arguments["--ref_zone"].replace(',',' ').split()))
    try:
        lin_start,lin_end, col_start, col_end = ref[0], ref[1], ref[2], ref[3]
    except:
        lin_start,lin_end = ref[0], ref[1]
        col_start, col_end = 0, ncol

logger.debug('Read crop zones time decomposition: {}, and empirical estimations'.\
  format(arguments["--crop"],arguments["--crop_emp"]))
if arguments["--crop"] ==  None:
    crop = [0,nlines,0,ncol]
else:
    crop = list(map(float,arguments["--crop"].replace(',',' ').split()))
    logger.warning('Crop time series data between lines {}-{} and cols {}-{}'.format(int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])))
ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])

# extract time series
try:
    maps_temp = np.zeros((nlines, ncol, N), dtype=np.float32)
    for band_index in range(1, N + 1):
        band = ds.GetRasterBand(band_index)
        maps_temp[:, :, band_index - 1] = band.ReadAsArray()
    logger.info('Load time series cube: {0}, with length: {1}'.format(arguments["--cube"], len(maps_temp.flatten())))
    kk = np.nonzero(np.logical_or(maps_temp==9990, maps_temp==9999))
    maps_temp[kk] = float('NaN')
except:
    cubei = np.fromfile(arguments["--cube"],dtype=np.float32)
    cube = as_strided(cubei[:nlines*ncol*N])
    logger.info('Load time series cube: {0}, with length: {1}'.format(arguments["--cube"], len(cube)))
    kk = np.flatnonzero(np.logical_or(cube==9990, cube==9999))
    cube[kk] = float('NaN')
    maps_temp = as_strided(cube.reshape((nlines,ncol,N)))
    del cube, cubei

# set at NaN zero values for all dates
# kk = np.nonzero(maps_temp[:,:,-1]==0)
cst = np.copy(maps_temp[:,:,imref])
cst[np.isnan(cst)] = 0.0
for l in range((N)):
    maps_temp[:,:,l] = maps_temp[:,:,l] - cst
    if l != imref:
        index = np.nonzero(maps_temp[:,:,l]==0.0)
        maps_temp[:,:,l][index] = float('NaN')

N=len(dates)
maps = np.copy(maps_temp[ibeg:iend,jbeg:jend,indexd])
logger.info('Number images between {0} and {1}: {2}'.format(dmin,dmax,N))
logger.info('Reshape cube: {}'.format(maps.shape))
new_lines, new_cols = maps.shape[0], maps.shape[1]

if arguments["--crop_emp"] ==  None:
    crop_emp = [0,new_lines,0,new_cols]
else:
    crop_emp = list(map(float,arguments["--crop_emp"].replace(',',' ').split()))
    logger.warning('Crop empirical estimation between lines {}-{} and cols {}-{}'.format(int(crop_emp[0]),int(crop_emp[1]),int(crop_emp[2]),int(crop_emp[3])))
ibeg_emp,iend_emp,jbeg_emp,jend_emp = int(crop_emp[0]),int(crop_emp[1]),int(crop_emp[2]),int(crop_emp[3])

# clean
del maps_temp

nfigure=0
# open mask file
if arguments["--mask"] is not None:
    extension = os.path.splitext(arguments["--mask"])[1]
    checkinfile(arguments["--mask"])
    if extension == ".tif":
      ds = gdal.Open(arguments["--mask"], gdal.GA_ReadOnly)
      band = ds.GetRasterBand(1)
      mask = band.ReadAsArray()[ibeg:iend,jbeg:jend]*float(arguments["--scale_mask"])
      del ds
    else:
      fid = open(arguments["--mask"],'r')
      mask = np.fromfile(fid,dtype=np.float32).reshape(nlines, ncol)[ibeg:iend,jbeg:jend]*float(arguments["--scale_mask"])
      fid.close()
    maski = mask.flatten()
else:
    mask_flat = np.ones((new_lines,new_cols))
    mask = np.ones((new_lines,new_cols))
    maski = mask.flatten()

# open elevation map
if arguments["--topofile"] is not None:
    extension = os.path.splitext(arguments["--topofile"])[1]
    checkinfile(arguments["--topofile"])
    if extension == ".tif":
      ds = gdal.Open(arguments["--topofile"], gdal.GA_ReadOnly)
      band = ds.GetRasterBand(1)
      elev = band.ReadAsArray()[ibeg:iend,jbeg:jend]
      del ds
    else:
      fid = open(arguments["--topofile"],'r')
      elev = np.fromfile(fid,dtype=np.float32).reshape(nlines, ncol)[ibeg:iend,jbeg:jend]
      fid.close()
    elev[np.isnan(maps[:,:,-1])] = float('NaN')
    kk = np.nonzero(abs(elev)>9999.)
    elev[kk] = float('NaN')
    elevi = elev.flatten()

    # define max min topo for empirical relationship
    maxtopo,mintopo = np.nanpercentile(elev,float(arguments["--perc_topo"])),np.nanpercentile(elev,100-float(arguments["--perc_topo"]))
    logger.info('Max-Min topography for empirical estimation: {0:.1f}-{1:.1f}'.format(maxtopo,mintopo))

else:
    elev = np.ones((new_lines,new_cols),dtype=np.float32)
    elevi = elev.flatten()
    maxtopo,mintopo = 2, 0 

if arguments["--aspect"] is not None:
    extension = os.path.splitext(arguments["--aspect"])[1]
    checkinfile(arguments["--aspect"])
    if extension == ".tif":
      ds = gdal.Open(arguments["--aspect"], gdal.GA_ReadOnly)
      band = ds.GetRasterBand(1)
      aspect = band.ReadAsArray()[ibeg:iend,jbeg:jend]
      del ds
    else:
      fid = open(arguments["--aspect"],'r')
      aspect = np.fromfile(fid,dtype=np.float32).reshape(nlines, ncol)[ibeg:iend,jbeg:jend]
      fid.close()
    aspect[np.isnan(maps[:,:,-1])] = float('NaN')
    kk = np.nonzero(abs(aspect>9999.))
    aspect[kk] = float('NaN')
    aspecti = aspect.flatten()
else:
    aspect = np.ones((new_lines,new_cols))
    aspecti = aspect.flatten()

if arguments["--rmspixel"] is not None:
    extension = os.path.splitext(arguments["--rmspixel"])[1]
    checkinfile(arguments["--rmspixel"])
    if extension == ".tif":
        ds = gdal.Open(arguments["--rmspixel"], gdal.GA_ReadOnly)
        band = ds.GetRasterBand(1)
        rmsmap = band.ReadAsArray()[ibeg:iend,jbeg:jend]
        del ds
    else:
        rmsmap = np.fromfile(arguments["--rmspixel"],dtype=np.float32).reshape((nlines,ncol))[ibeg:iend,jbeg:jend]

    kk = np.nonzero(np.logical_or(rmsmap==0.0, rmsmap>999.))
    rmsmap[kk] = float('NaN')
    kk = np.nonzero(rmsmap>float(arguments["--threshold_rms"]))
    spacial_mask = np.copy(rmsmap)
    spacial_mask[kk] = float('NaN')
    fig = plt.figure(nfigure,figsize=(9,4))
    nfigure = nfigure + 1
    ax = fig.add_subplot(1,1,1)
    cax = ax.imshow(spacial_mask,cmap=cmap,interpolation='nearest')
    ax.set_title('Mask on spatial estimation based on RMSpixel')
    plt.setp( ax.get_xticklabels(), visible=False)
    fig.colorbar(cax, orientation='vertical',aspect=10)
    del spacial_mask
else:
    rmsmap = np.ones((new_lines,new_cols))
    spacial_mask = np.ones((new_lines,new_cols))
    arguments["--threshold_rms"] = 2.

if arguments["--bperp"]=='yes':
  for idx, bp in enumerate(bases):
    # Création d'une figure pour chaque bp dans bases
    fig = plt.figure(nfigure, figsize=(10, 4))
    nfigure += 1

    # Graphique 1: Baseline perpendiculaire en fonction du temps
    ax1 = fig.add_subplot(1, 2, 1)

    # selection des baselines non nulles
    idx = np.flatnonzero(bp)
    dates_temp = dates[idx]
    idates_temp = idates[idx]
    bp_temp = bp[idx]
    N_temp = len(bp_temp)

    # Conversion des dates au format numérique pour matplotlib
    x = [mdates.date2num(dt.strptime(f'{d}', '%Y%m%d')) for d in idates_temp]

    # Configuration de l'affichage des dates
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
    ax1.plot(x, bp_temp, "ro", label=f'Baseline history of the {N_temp} images')
    ax1.plot(x, bp_temp, "green")
    fig.autofmt_xdate()

    ax1.set_xlabel('Time (Year/month/day)')
    ax1.set_ylabel('Perpendicular Baseline')
    ax1.legend(loc='best')

    # Graphique 2: Baseline en fonction de la saisonnalité (modulo 1 année)
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.plot(np.mod(dates_temp, 1), bp_temp, "ro", label=f'Baseline seasonality of the {N_temp} images')
    ax2.legend(loc='best')

    # Sauvegarde de la figure et des données
    fig.savefig(f'baseline_{idx}.eps', format='EPS', dpi=150)
    np.savetxt(f'bp_t_{idx}.in', np.vstack([dates_temp, bp_temp]).T, fmt='%.6f')

if arguments["--mask"] is not None:
    los_temp = as_strided(mask[ibeg_emp:iend_emp,jbeg_emp:jend_emp]).flatten()

    if arguments["--rampmask"]=='yes':
        logger.info('Flatten mask...')
        temp = [(i,j) for i in range(iend_emp-ibeg_emp) for j in range(jend_emp-jbeg_emp) \
        if np.logical_and((math.isnan(los_temp[i*(jend_emp-jbeg_emp)+j]) is False), \
            (los_temp[i*(jend_emp-jbeg_emp)+j]>float(arguments["--threshold_mask"])))]

        temp2 = np.array(temp)
        x = temp2[:,0]; y = temp2[:,1]
        los_clean = los_temp[x*(new_cols)+y]
        G=np.zeros((len(los_clean),4))
        G[:,0], G[:,1], G[:,2], G[:,3] = y**2, y, x, 1
        # ramp inversion
        pars = np.dot(np.dot(np.linalg.inv(np.dot(G.T,G)),G.T),los_clean)
        a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
        logger.info('Remove ramp mask %f x**2 %f x  + %f y + %f for : %s'%(a,b,c,d,arguments["--mask"]))

        # remove 0 values
        kk = np.flatnonzero(np.logical_or(maski==0, maski==9999))
        #kk = np.flatnonzero(los==9999)
        maski[kk] = float('NaN')
        G=np.zeros((len(maski),4))
        for i in range(nlines):
            G[i*ncol:(i+1)*ncol,0] = (np.arange((ncol)) - jbeg_emp)**2
            G[i*ncol:(i+1)*ncol,1] = np.arange((ncol)) - jbeg_emp
            G[i*ncol:(i+1)*ncol,2] = i - ibeg_emp
        G[:,3] = 1
        mask_flat = (maski - np.dot(G,pars)).reshape(new_lines,new_cols)
        mask_flat = mask_flat - np.nanmean(mask_flat)

    else:

        # remove 0 values
        kk = np.flatnonzero(np.logical_or(np.logical_or(maski==0, maski==9999),np.isnan(los_temp)))
        #kk = np.flatnonzero(los==9999)
        maski[kk] = float('NaN')
        mask_flat = maski.reshape(new_lines,new_cols)

    del maski

    # check seuil
    kk = np.flatnonzero(mask_flat>float(arguments["--threshold_mask"]))
    mask_flat_clean=np.copy(mask_flat.flatten())
    mask_flat_clean[kk]=float('NaN')
    mask_flat_clean = mask_flat_clean.reshape(new_lines,new_cols)

    # mask maps if necessary for temporal inversion
    if arguments["--tempmask"]=='yes':
        kk = np.nonzero(np.logical_or(mask_flat<float(arguments["--threshold_mask"]),
          np.isnan(mask_flat)))
        for l in range((N)):
            # clean only selected area
            d = as_strided(maps[ibeg_emp:iend_emp,jbeg_emp:jend_emp,l])
            d[kk] = float('NaN')

    # plots
    nfigure+=1
    fig = plt.figure(nfigure,figsize=(7,6))
    vmax = np.abs([np.nanmedian(mask_flat) + np.nanstd(mask_flat),\
        np.nanmedian(mask_flat) - np.nanstd(mask_flat)]).max()
    vmin = -vmax

    ax = fig.add_subplot(1,3,1)
    cax = ax.imshow(mask,cmap=cmap,vmax=vmax,vmin=vmin,interpolation='nearest')
    ax.set_title('Original Mask')
    plt.setp( ax.get_xticklabels(), visible=False)

    ax = fig.add_subplot(1,3,2)
    cax = ax.imshow(mask_flat,cmap=cmap,vmax=vmax,vmin=vmin,interpolation='nearest')
    ax.set_title('Flat Mask')
    plt.setp( ax.get_xticklabels(), visible=False)
    #cbar = fig.colorbar(cax, orientation='vertical',aspect=10)

    ax = fig.add_subplot(1,3,3)
    cax = ax.imshow(mask_flat_clean,cmap=cmap,vmax=vmax,vmin=vmin,interpolation='nearest')
    ax.set_title('Final Mask')
    plt.setp( ax.get_xticklabels(), visible=False)
    #cbar = fig.colorbar(cax, orientation='vertical',aspect=10)
    fig.savefig('mask.eps', format='EPS',dpi=150)
    del mask_flat_clean

# plot diplacements maps
nfigure+=1
plot_displacement_maps(maps, idates, nfigure=nfigure, cmap=cmap, plot=plot, title='Time series maps', filename='maps.eps')

#######################################################
# Save new lect.in file
#######################################################

fid = open('lect_ts.in','w')
np.savetxt(fid, (new_cols,new_lines,N),fmt='%6i',newline='\t')
fid.close()

#######################################################
# Create functions of decomposition
######################################################

if arguments["--steps"] ==  None:
    cos = []
else:
    cos = list(map(float,arguments["--steps"].replace(',',' ').split()))

if arguments["--postseismic"] ==  None:
    pos = []
else:
    pos = list(map(float,arguments["--postseismic"].replace('None','-1').replace(',',' ').split()))

if len(pos)>0 and len(cos) != len(pos):
    raise Exception("coseimic and postseismic lists are not the same size")

if arguments["--slowslip"] == None:
    sse, sse_time, sse_car = [], [], []
else:
    try:
        sse = list(map(float,arguments["--slowslip"].replace(',',' ').split()))
        sse_time = sse[::2]
        sse_car = sse[1::2]
    except:
        sse_time, sse_car = [], []

if arguments["--vector"] != None:
    vectf = arguments["--vector"].replace(',',' ').split()
else:
    vectf = []
    vect = None

basis=[
    reference(name='reference',date=datemin,reduction='ref'),
    ]
index = len(basis)

# initialise iteration with linear alone
iteration=False

if arguments["--linear"]=='yes':
    indexinter=index
    basis.append(linear(name='linear',reduction='lin',date=datemin))
    index = index + 1

if arguments["--seasonal"] =='yes':
    indexseas = index
    basis.append(cosvar(name='seas. var (cos)',reduction='cos',date=datemin))
    basis.append(sinvar(name='seas. var (sin)',reduction='sin',date=datemin))
    index = index + 2

if arguments["--seasonal_increase"] =='yes':
   # 2
   indexseast = index
   basis.append(cost(name='increased seas. var (cos)',reduction='cost',date=datemin))
   basis.append(sint(name='increased seas. var (sin)',reduction='sint',date=datemin))
   index = index + 2

if arguments["--semianual"]=='yes':
     indexsemi = index
     basis.append(cos2var(name='semi-anual var (cos)',reduction='cosw2t',date=datemin))
     basis.append(sin2var(name='semi-anual var (sin)',reduction='sinw2t',date=datemin))
     index = index + 2

if arguments["--bianual"]=='yes':
     indexbi = index
     basis.append(cos5var(name='bi-anual var (cos)',reduction='cos5wt',date=datemin))
     basis.append(sin5var(name='bi-anual var (sin)',reduction='sin5wt',date=datemin))
     index = index + 2

indexco = np.zeros(len(cos))
for i in range(len(cos)):
    basis.append(steps(name='steps {}'.format(i),reduction='step{}'.format(i),date=cos[i])),
    indexco[i] = index
    index = index + 1
    iteration=True

indexpo,indexpofull = [],[]
for i in range(len(pos)):
  if pos[i] > 0. :
    basis.append(postseismic(name='postseismic {}'.format(i),reduction='log{}'.format(i),date=cos[i],tcar=pos[i])),
    indexpo.append(int(index))
    indexpofull.append(int(index))
    index = index + 1
    arguments["--ineq"] = 'yes'
  else:
    indexpofull.append(0)
indexpo = np.array(indexpo)
indexpofull = np.array(indexpofull)

indexsse = np.zeros(len(sse_time))
for i in range(len(sse_time)):
    basis.append(slowslip(name='sse {}'.format(i),reduction='sse{}'.format(i),date=sse_time[i],tcar=sse_car[i])),
    indexsse[i] = int(index)
    index = index + 1
    iteration=True

kernels = []
if arguments["--bperp"]=='yes':
    indexbperp = []
    i = 0
    for bp,bref in zip(bases,baserefs):
        i = i + 1
        kernels.append(corrBperp(name='bperp correction', reduction=f'corrBperp{i}', bp0=bref, bp=bp))
        indexbperp.append(index)
        index = index + 1

if arguments["--vector"] != None:
    fig = plt.figure(nfigure,figsize=(6,4))
    nfigure = nfigure + 1
    indexvect = np.zeros(len(vectf))
    for i in range(len(vectf)):
      ax = fig.add_subplot(len(vectf),1,i+1)
      v = np.loadtxt(vectf[i], comments='#', unpack = False, dtype='f')
      kernels.append(vector(name=vectf[i],reduction='vector_{}'.format(i),vect=v))
      ax.plot(v,label='Vector')
      plt.legend(loc='best')
      indexvect[i] = index
      index = index + 1    
    indexvect = indexvect.astype(int)
    if plot=='yes':
      plt.show()
    # sys.exit()

indexpo = indexpo.astype(int)
indexco = indexco.astype(int)
indexsse = indexsse.astype(int)

eguality = False
if arguments["--seasonal_increase"] == 'yes' and arguments["--seasonal"] == 'yes':
    eguality = True
    arguments["--ineq"] = 'yes'

print()
Mbasis=len(basis)
logger.info('Number of basis functions: {}'.format(Mbasis))
Mker=len(kernels)
logger.info('Number of kernel functions: {}'.format(Mker))
M = Mbasis + Mker

print('Basis functions, Time:')
for i in range((Mbasis)):
    basis[i].info()
if Mker > 1:
  print('Kernels functions, Time:')
for i in range((Mker)):
    kernels[i].info()

# initialize matrix model to NaN
for l in range((Mbasis)):
    basis[l].m = np.ones((new_lines,new_cols))*float('NaN')
    basis[l].sigmam = np.ones((new_lines,new_cols))*float('NaN')
for l in range((Mker)):
    kernels[l].m = np.ones((new_lines,new_cols))*float('NaN')
    kernels[l].sigmam = np.ones((new_lines,new_cols))*float('NaN')

# initialize aps
if  arguments["--aps"] is None:
    in_aps = np.ones((N)) # no weigthing for the first itertion
else:
    fimages =  arguments["--aps"]
    try:
        in_aps = np.loadtxt(fimages, unpack=True, comments='#', usecols=(2), dtype='f')
    except:
        logger.warning('APS file is in decrepicated format, requiered two columns text file`')
        in_aps = np.loadtxt(fimages, comments='#', dtype='f')
    logger.info('Input APS: {}'.format(in_aps))
    logger.info('Set very low values to the 2 percentile to avoid overweighting...')
    min_aps= np.nanpercentile(in_aps,2)
    index = np.flatnonzero(in_aps<min_aps)
    in_aps[index] = min_aps
    in_aps = in_aps[indexd]

# initialize rms
if  arguments["--rms"] is None:
    in_rms = np.ones((N)) # no weigthing for the first itertion
else:
    fimages =  arguments["--rms"]
    try:
        in_rms = np.loadtxt(fimages, unpack=True, comments='#', usecols=(2), dtype='f')
    except:
        logger.warning('RMS file is in decrepicated format, requiered two columns text file`')
        in_rms = np.loadtxt(fimages, comments='#', dtype='f')
    logger.info('Input RMS: {}'.format(in_rms))
    logger.info('Set very low values to the 2 percentile to avoid overweighting...')
    min_rms= np.nanpercentile(in_rms,2)
    index = np.flatnonzero(in_rms < min_rms)
    in_rms[index] = min_rms
    in_rms = in_rms[indexd]

## initialize input uncertainties
in_sigma = in_rms * in_aps
logger.info('Input uncertainties: {}'.format(in_rms))

## inversion procedure
def consInvert(A,b,sigmad,ineq='yes',cond=1.0e-3, iter=100,acc=1e-6, eguality=False):
    '''Solves the constrained inversion problem.

    Minimize:

    ||Ax-b||^2

    Subject to:
    mmin < m < mmax
    '''

    if A.shape[0] != len(b):
        raise ValueError('Incompatible dimensions for A and b')

    if ineq == 'no':
        try:
          Cd = np.diag(sigmad**2, k = 0)
          fsoln = np.dot(np.linalg.inv(np.dot(np.dot(A.T,np.linalg.inv(Cd)),A)),np.dot(np.dot(A.T,np.linalg.inv(Cd)),b))
        except:
          fsoln = lst.lstsq(A,b,rcond=None)[0]  
    else:
        if len(indexpo>0):
          # invert first without post-seismic
          Ain = np.delete(A,indexpo,1)
          mtemp = lst.lstsq(Ain,b,rcond=cond)[0]
    
          # rebuild full vector
          for z in range(len(indexpo)):
            mtemp = np.insert(mtemp,indexpo[z],0)
          minit = np.copy(mtemp)
          # # initialize bounds
          mmin,mmax = -np.ones(len(minit))*np.inf, np.ones(len(minit))*np.inf

          # We here define bounds for postseismic to be the same sign than steps
          # and steps inferior or egual to the coseimic initial
          for i in range(len(indexco)):
            if (pos[i] > 0.) and (minit[int(indexco[i])]>0.):
                mmin[int(indexpofull[i])], mmax[int(indexpofull[i])] = 0, np.inf
                mmin[int(indexco[i])], mmax[int(indexco[i])] = 0, minit[int(indexco[i])]
            if (pos[i] > 0.) and (minit[int(indexco[i])]<0.):
                mmin[int(indexpofull[i])], mmax[int(indexpofull[i])] = -np.inf , 0
                mmin[int(indexco[i])], mmax[int(indexco[i])] = minit[int(indexco[i])], 0

        else:
          minit = lst.lstsq(A,b,rcond=None)[0]
          mmin,mmax = -np.ones(len(minit))*np.inf, np.ones(len(minit))*np.inf

        bounds=list(zip(mmin,mmax))
        def eq_cond(x, *args):
           return (x[indexseast+1]/x[indexseast]) - (x[indexseas+1]/x[indexseas])

        ####Objective function and derivative
        _func = lambda x: np.sum(((np.dot(A,x)-b)/sigmad)**2)
        _fprime = lambda x: 2*np.dot(A.T/sigmad, (np.dot(A,x)-b)/sigmad)
        if eguality:
            res = opt.fmin_slsqp(_func,minit,bounds=bounds,fprime=_fprime,eqcons=[eq_cond], \
                iter=iter,full_output=True,iprint=0,acc=acc)
        else:
            res = opt.fmin_slsqp(_func,minit,bounds=bounds,fprime=_fprime, \
                iter=iter,full_output=True,iprint=0,acc=acc)
        fsoln = res[0]
        #print('Optimization:', fsoln)

    try:
       varx = np.linalg.inv(np.dot(A.T,A))
       res2 = np.sum(pow((b-np.dot(A,fsoln)),2))
       scale = 1./(A.shape[0]-A.shape[1])
       sigmam = np.sqrt(scale*res2*np.diag(varx))
    except:
       sigmam = np.ones((A.shape[1]))*float('NaN')
    return fsoln,sigmam

def linear_inv(G, data, sigma):
      'Iterative linear inversion'

      x0 = lst.lstsq(G,data)[0]
      _func = lambda x: np.sum(((np.dot(G,x)-data)/sigma)**2)
      _fprime = lambda x: 2*np.dot(G.T/sigma, (np.dot(G,x)-data)/sigma)
      pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1.e-9)[0]

      return pars

def estim_ramp(los,los_clean,topo_clean,az,rg,order,sigma,nfit,ivar,l,ax_dphi):
      'Ramp/Topo estimation and correction. Estimation is performed on sliding median'

      # initialize topo
      topo = np.zeros((new_lines,new_cols))
      ramp = np.zeros((new_lines,new_cols))      

      # y: range, x: azimuth
      if arguments["--topofile"] is None:
          topobins = topo_clean
          rgbins, azbins = rg, az
          data = np.copy(los_clean)

      else:
          # lets try to digitize to improve the fit
          # digitize data in bins, compute median and std
          bins = np.arange(mintopo,maxtopo,abs(maxtopo-mintopo)/500.)
          inds = np.digitize(topo_clean,bins)
          topobins = []
          losbins = []
          losstd = []
          azbins, rgbins = [], []
          los_clean2, topo_clean2, az_clean2, rg_clean2, sigma_clean2 = [], [], [], [], []
          for j in range(len(bins)-1):
                  uu = np.flatnonzero(inds == j)
                  if len(uu)>200:
                      topobins.append(bins[j] + (bins[j+1] - bins[j])/2.)

                      # do a small clean within the bin
                      indice = np.flatnonzero(np.logical_and(los_clean[uu]>np.percentile(\
                          los_clean[uu],100-float(arguments["--perc_los"])),los_clean[uu]<np.percentile(los_clean[uu],float(arguments["--perc_los"]))))

                      losstd.append(np.nanstd(los_clean[uu][indice]))
                      losbins.append(np.nanmedian(los_clean[uu][indice]))
                      azbins.append(np.nanmedian(az[uu][indice]))
                      rgbins.append(np.nanmedian(rg[uu][indice]))

                      # remove outliers from data
                      los_clean2.append(los_clean[uu][indice])
                      az_clean2.append(az[uu][indice])
                      rg_clean2.append(rg[uu][indice])
                      topo_clean2.append(topo_clean[uu][indice])
                      # new sigma is clean sigma time the standard deviation of the los within the bin
                      sigma_clean2.append(sigma[uu][indice]*np.nanstd(los_clean[uu][indice]))

          topobins = np.array(topobins)
          rgbins, azbins = np.array(rgbins),np.array(azbins)

          # data and model cleaned a second time by sliding median
          los_clean = np.concatenate(los_clean2)
          az, rg = np.concatenate(az_clean2), np.concatenate(rg_clean2)
          topo_clean = np.concatenate(topo_clean2)
          sigma = np.concatenate(sigma_clean2)
          del los_clean2, topo_clean2, az_clean2, rg_clean2, sigma_clean2 
          
          if order == 0 and ivar == 0:
              data = np.array(losbins)
              sigma = np.array(losstd)

          else:
              data = np.array(los_clean)

          if len(data) < 10:
            logger.critical('Too small area for empirical phase/topo relationship. Re-defined crop values Exit!')
            sys.exit()
    
      if order==0:

        if arguments["--topofile"] is None:

            a = 0.
            ramp = np.zeros((new_lines,new_cols))
            res = los 

        else:

            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),2))
                G[:,0] = 1
                G[:,1] = topobins

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]
                logger.info('Remove ref frame %f + %f z for date: %i'%(a,b,idates[l]))

                # plot phase/elev
                funct = a
                funcbins = a
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,b*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),2))
                G[:,0] = 1
                G[:,1] = elevi

                res = los - np.dot(G,pars)
                topo = np.dot(G,pars).reshape(new_lines,new_cols)


            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),3))
                G[:,0] = 1
                G[:,1] = topobins
                G[:,2] = topobins**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c=pars[2]
                print ('Remove ref frame %f + %f z + %f z**2 for date: %i'%(a,b,c,idates[l]))

                # plot phase/elev
                funct = a
                funcbins = a
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,b*x+c*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),3))
                G[:,0] = 1
                G[:,1] = elevi
                G[:,2] = elevi**2

                res = los - np.dot(G,pars)
                topo = np.dot(G,pars).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),3))
                G[:,0] = 1
                G[:,1] = topo_clean
                G[:,2] = az*topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]
                print ('Remove ref frame %f + %f z + %f az*z for date: %i'%(a,b,c,idates[l]))

                # plot phase/elev
                funct = a + c*topo_clean*az
                funcbins = a + c*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,b*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),3))
                G[:,0] = 1
                G[:,1] = elevi
                G[:,2] = elevi
                for i in range(nlines):
                    G[i*new_cols:(i+1)*new_cols,2] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)

                topo = np.dot(G,pars).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),4))
                G[:,0] = 1
                G[:,1] = az*topo_clean
                G[:,2] = topo_clean
                G[:,3] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
                print ('Remove ref frame %f + %f az*z + %f z + %f z**2 for date: %i'%(a,b,c,d,idates[l]))

                # plot phase/elev
                funct = a + b*topo_clean*az
                funcbins = a + b*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,c*x+d*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                G[:,0] = 1
                G[:,1] = elevi
                G[:,2] = elevi
                G[:,3] = elevi**2
                for i in range(nlines):
                    G[i*new_cols:(i+1)*new_cols,1] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                topo = np.dot(G,pars).reshape(new_lines,new_cols)

      elif order==1: # Remove a range ramp ay+b for each maps (y = col)

        if arguments["--topofile"] is None:
            G=np.zeros((len(data),2))
            G[:,0] = rg
            G[:,1] = 1

            # ramp inversion
            x0 = lst.lstsq(G,data)[0]
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]
            print ('Remove ramp %f r + %f for date: %i'%(a,b,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),2))
            for i in range(nlines):
                G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
            G[:,1] = 1

            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),3))
                G[:,0] = rg
                G[:,1] = 1
                G[:,2] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]
                print ('Remove ramp %f r + %f + %f z for date: %i'%(a,b,c,idates[l]))

                # plot phase/elev
                funct = a*rg + b
                funcbins = a*rgbins + b
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,c*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),3))
                for i in range(nlines):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbeg_emp
                G[:,1] = 1
                G[:,2] = elevi

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),4))
                G[:,0] = rg
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d=pars[3]
                print ('Remove ramp %f r + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,idates[l]))

                # plot phase/elev
                funct = a*rg+ b
                funcbins = a*rgbins+ b
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,c*x+d*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                for i in range(nlines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                G[:,1] = 1
                G[:,2] = elevi
                G[:,3] = elevi**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),4))
                G[:,0] = rg
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
                print ('Remove ramp %f r + %f + %f z + %f z*az for date: %i'%(a,b,c,d,idates[l]))

                # plot phase/elev
                funct = a*rg+ b + d*topo_clean*az
                funcbins = a*rgbins+ b + d*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,c*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                G[:,1] = 1
                G[:,2] = elevi
                G[:,3] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,3] *= (i - ibeg_emp)


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),5))
                G[:,0] = rg
                G[:,1] = 1
                G[:,2] = topo_clean*az
                G[:,3] = topo_clean
                G[:,4] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
                print ('Remove ramp %f r + %f +  %f z*az + %f z + %f z**2 for date: %i'%(a,b,c,d,e,idates[l]))

                # plot phase/elev
                funct = a*rg+ b + c*topo_clean*az
                funcbins = a*rgbins+ b + c*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,d*x+e*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                G[:,1] = 1
                G[:,2] = elevi
                G[:,3] = elevi
                G[:,4] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)


      elif order==2: # Remove a azimutal ramp ax+b for each maps (x is lign)
        if arguments["--topofile"] is None:
            G=np.zeros((len(data),2))
            G[:,0] = az
            G[:,1] = 1

            # ramp inversion
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]
            print ('Remove ramp %f az + %f for date: %i'%(a,b,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),2))
            for i in range(new_lines):
                G[i*new_cols:(i+1)*new_cols,0] =(i - ibeg_emp)
            G[:,1] = 1


            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),3))
                G[:,0] = az
                G[:,1] = 1
                G[:,2] = topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]
                print ('Remove ramp %f az + %f + %f z for date: %i'%(a,b,c,idates[l]))

                # plot phase/elev
                funct = a*az + b
                funcbins = a*azbins + b
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,c*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),3))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] =(i - ibeg_emp)
                G[:,1] = 1
                G[:,2] = elevi


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),4))
                G[:,0] = az
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
                print ('Remove ramp %f az + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,idates[l]))

                # plot phase/elev
                funct = a*az + b
                funcbins = a*azbins + b
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,c*x + d*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] =(i - ibeg_emp)
                G[:,1] = 1
                G[:,2] = elevi
                G[:,3] = elevi**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),4))
                G[:,0] = az
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
                print ('Remove ramp %f az + %f + %f z + %f z*az for date: %i'%(a,b,c,d,idates[l]))

                # plot phase/elev
                funct = a*az + b + d*topo_clean*az
                funcbins = a*azbins + b + d*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,c*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                G[:,1] = 1
                G[:,2] = elevi
                G[:,3] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,3] *= (i - ibeg_emp)


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),5))
                G[:,0] = az
                G[:,1] = 1
                G[:,2] = topo_clean*az
                G[:,3] = topo_clean
                G[:,4] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
                print ('Remove ramp %f az + %f + %f z*az + %f z + %f z**2 for date: %i'%(a,b,c,d,e,idates[l]))

                # plot phase/elev
                funct = a*az + b + c*topo_clean*az
                funcbins = a*azbins + b + c*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,d*x+e*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                G[:,1] = 1
                G[:,2] = elevi
                G[:,3] = elevi
                G[:,4] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)

      elif order==3: # Remove a ramp ay+bx+c for each maps
        if arguments["--topofile"] is None:
            G=np.zeros((len(data),3))
            G[:,0] = rg
            G[:,1] = az
            G[:,2] = 1

            # ramp inversion
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]; c = pars[2]
            print ('Remove ramp %f r  + %f az + %f for date: %i'%(a,b,c,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),3))
            for i in range(new_lines):
                G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                G[i*new_cols:(i+1)*new_cols,1] = (i - ibeg_emp)
            G[:,2] = 1

            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),4))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = 1
                G[:,3] = topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
                print ('Remove ramp %f r  + %f az + %f + %f z for date: %i'%(a,b,c,d,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c
                funcbins = a*rgbins+ b*azbins + c
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,d*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] =(i - ibeg_emp)
                G[:,2] = 1
                G[:,3] = elevi


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),5))
                G[:-1,0] = rg
                G[:-1,1] = az
                G[:,2] = 1
                G[:-1,3] = topo_clean
                G[:-1,4] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
                print ('Remove ramp %f r  + %f az + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,e,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c
                funcbins = a*rgbins+ b*azbins + c
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,d*x+e*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] =(i - ibeg_emp)
                G[:,2] = 1
                G[:,3] = elevi
                G[:,4] = elevi**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),5))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = 1
                G[:,3] = topo_clean
                G[:,4] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e=pars[4]
                print ('Remove ramp %f r  + %f az + %f + %f z +  %f z*az for date: %i'%(a,b,c,d,e,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c + e*topo_clean*az
                funcbins = a*rgbins+ b*azbins + c + e*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,d*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                G[:,2] = 1
                G[:,3] = elevi
                G[:,4] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] =(i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,4] *= (i - ibeg_emp)


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),6))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = 1
                G[:,3] = topo_clean*az
                G[:,4] = topo_clean
                G[:,5] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e=pars[4]; f=pars[5]
                print ('Remove ramp %f r  + %f az + %f +  %f z*az + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c + d*topo_clean*az
                funcbins = a*rgbins+ b*azbins + c + d*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.1, alpha=0.01, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x+f*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,2] = 1
                G[:,3] = elevi
                G[:,4] = elevi
                G[:,5] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] =(i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,3] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)

      elif order==4:
        if arguments["--topofile"] is None:
            G=np.zeros((len(data),4))
            G[:,0] = rg
            G[:,1] = az
            G[:,2] = rg*az
            G[:,3] = 1

            # ramp inversion
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
            print ('Remove ramp %f r %f az  + %f r*az + %f for date: %i'%(a,b,c,d,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),4))
            for i in range(new_lines):
                G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                G[i*new_cols:(i+1)*new_cols,2] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
            G[:,3] = 1

            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),5))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = rg*az
                G[:,3] = 1
                G[:,4] = topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]

                print ('Remove ramp %f r, %f az  + %f r*az + %f + %f z for date: %i'%(a,b,c,d,e,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c*az*rg+ d
                funcbins = a*rgbins+ b*azbins + c*azbins*rgbins+ d
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
                G[:,3] = 1
                G[:,4] = elevi


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),6))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = rg*az
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]

                print ('Remove ramp %f r, %f az  + %f r*az + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c*az*rg+ d
                funcbins = a*rgbins+ b*azbins + c*azbins*rgbins+ d
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x+f*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),6))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = rg*az
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]

                print ('Remove ramp %f r, %f az  + %f r*az + %f + %f z + %f az*z for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c*az*rg+ d + f*topo_clean*az
                funcbins = a*rgbins+ b*azbins + c*azbins*rgbins+ d + f*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
                    G[i*new_cols:(i+1)*new_cols,5] *=  (i - ibeg_emp)


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),6))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = rg*az
                G[:,3] = 1
                G[:,4] = topo_clean*az
                G[:,5] = topo_clean
                G[:,6] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]

                print ('Remove ramp %f r, %f az  + %f r*az + %f + + %f az*z +  %f z + %f z**2  for date: %i'%(a,b,c,d,e,f,g,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c*az*rg+ d + e*topo_clean*az
                funcbins = a*rgbins+ b*azbins + c*azbins*rgbins+ d + e*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,f*x+g*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi
                G[:,6] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
                    G[i*new_cols:(i+1)*new_cols,4] *=  (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)

      elif order==5:

        if arguments["--topofile"] is None:
            G=np.zeros((len(data),4))
            G[:,0] = rg**2
            G[:,1] = rg
            G[:,2] = az
            G[:,3] = 1

            # ramp inversion
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
            print ('Remove ramp %f r**2 + %f r  + %f az + %f for date: %i'%(a,b,c,d,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),4))
            for i in range(new_lines):
                G[i*new_cols:(i+1)*new_cols,0] = (np.arange((new_cols)) - jbeg_emp)**2
                G[i*new_cols:(i+1)*new_cols,1] = np.arange((new_cols)) - jbeg_emp
                G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
            G[:,3] = 1


            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:

            if (ivar==0 and nfit==0):

                G=np.zeros((len(data),5))
                G[:,0] = rg**2
                G[:,1] = rg
                G[:,2] = az
                G[:,3] = 1
                G[:,4] = topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
                print ('Remove ramp %f r**2 + %f r  + %f az + %f + %f z for date: %i'%(a,b,c,d,e,idates[l]))

                # plot phase/elev
                funct = a*rg**2 + b*rg+ c*az + d
                funcbins = a*rgbins**2 + b*rgbins+ c*azbins + d
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = np.arange((new_cols)) -  jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                G[:,3] = 1
                G[:,4] = elevi


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),6))
                G[:,0] = rg**2
                G[:,1] = rg
                G[:,2] = az
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
                print ('Remove ramp %f r**2 + %f r  + %f az + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funct = a*rg**2 + b*rg+ c*az + d
                funcbins = a*rgbins**2 + b*rgbins+ c*azbins + d
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x+f*x**2,'-r', lw =4.)


                # build total G matrix
                G=np.zeros((len(los),6))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = np.arange((new_cols)) -  jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi**2

                res = los - np.dot(G,pars)

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==0):

                G=np.zeros((len(data),6))
                G[:,0] = rg**2
                G[:,1] = rg
                G[:,2] = az
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
                print ('Remove ramp %f r**2 + %f r  + %f az + %f + %f z + %f z*az for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funct = a*rg**2 + b*rg+ c*az + d + f*topo_clean*az
                funcbins = a*rgbins**2 + b*rgbins+ c*azbins + d + f*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = np.arange((new_cols)) -  jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,5] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)


            elif (ivar==1 and nfit==1):

                G=np.zeros((len(data),7))
                G[:,0] = rg**2
                G[:,1] = rg
                G[:,2] = az
                G[:,3] = 1
                G[:,4] = topo_clean*az
                G[:,5] = topo_clean
                G[:,6] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]
                print ('Remove ramp %f r**2 + %f r  + %f az + %f + + %f z*az + %f z +%f z**2 for date: %i'%(a,b,c,d,e,f,g,idates[l]))

                # plot phase/elev
                funct = a*rg**2 + b*rg+ c*az + d + e*topo_clean*az
                funcbins = a*rgbins**2 + b*rgbins+ c*azbins + d + e*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,f*x+g*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi
                G[:,6] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = np.arange((new_cols)) -  jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,4] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)

            else:
                pass

      elif order==6:
        if arguments["--topofile"] is None:
            G=np.zeros((len(data),4))
            G[:,0] = az**2
            G[:,1] = az
            G[:,2] = rg
            G[:,3] = 1

            # ramp inversion
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
            print ('Remove ramp %f az**2 + %f az  + %f r + %f for date: %i'%(a,b,c,d,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),4))
            for i in range(new_lines):
                G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                G[i*new_cols:(i+1)*new_cols,1] = (i - ibeg_emp)
                G[i*new_cols:(i+1)*new_cols,2] = np.arange((new_cols)) - jbeg_emp
            G[:,3] = 1


            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:
            if (ivar==0 and nfit==0) :
                G=np.zeros((len(data),5))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = rg
                G[:,3] = 1
                G[:,4] = topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
                print ('Remove ramp %f az**2 + %f az  + %f r + %f + %f z for date: %i'%(a,b,c,d,e,idates[l]))

                # plot phase/elev
                funct = a*az**2 + b*az + c*rg+ d
                funcbins = a*azbins**2 + b*azbins + c*rgbins+ d
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = i -  ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = np.arange((new_cols)) - jbeg_emp
                G[:,3] = 1
                G[:,4] = elevi


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),6))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = rg
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
                print ('Remove ramp %f az**2 + %f az  + %f r + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funct = a*az**2 + b*az + c*rg+ d
                funcbins = a*azbins**2 + b*azbins + c*rgbins+ d
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x+f*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = i -  ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = np.arange((new_cols)) - jbeg_emp
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),6))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = rg
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
                print ('Remove ramp %f az**2 + %f az  + %f r + %f + %f z + %f z*az for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funct = a*az**2 + b*az + c*rg+ d + f*topo_clean*az
                funcbins = a*azbins**2 + b*azbins + c*rgbins+ d + f*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = i -  ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,5] *= (i - ibeg_emp)


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),8))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = rg
                G[:,3] = 1
                G[:,4] = topo_clean*az
                G[:,5] = topo_clean
                G[:,6] = topo_clean**2
                G[:,7] = (topo_clean*az)**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g=pars[6]; h = pars[7]
                print ('Remove ramp %f az**2 + %f az  + %f r + %f + %f z*az + %f z + %f z**2 + %f (z*az)**2 for date: %i'%(a,b,c,d,e,f,g,h,idates[l]))

                # plot phase/elev
                funct = a*az**2 + b*az + c*rg+ d + e*topo_clean*az + h*(topo_clean*az)**2
                funcbins = a*azbins**2 + b*azbins + c*rgbins+ d + e*topobins*azbins + h*(topobins*azbins)**2
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,f*x+g*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),8))
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi
                G[:,6] = elevi**2
                G[:,7] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = i -  ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,4] *= (i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,7] *= (i - ibeg_emp)**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)


      elif order==7:
        if arguments["--topofile"] is None:
            G=np.zeros((len(data),5))
            G[:,0] = az**2
            G[:,1] = az
            G[:,2] = rg**2
            G[:,3] = rg
            G[:,4] = 1

            # ramp inversion
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
            print ('Remove ramp %f az**2 + %f az  + %f r**2 + %f r + %f for date: %i'%(a,b,c,d,e,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),5))
            for i in range(new_lines):
                G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                G[i*new_cols:(i+1)*new_cols,2] = (np.arange((new_cols)) - jbeg_emp)**2
                G[i*new_cols:(i+1)*new_cols,3] = np.arange((new_cols)) - jbeg_emp
            G[:,4] = 1


            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:
            if (ivar==0 and nfit ==0):
                G=np.zeros((len(data),6))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = rg**2
                G[:,3] = rg
                G[:,4] = 1
                G[:,5] = topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
                print ('Remove ramp %f az**2 + %f az  + %f r**2 + %f r + %f + %f z for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funct = a*az**2 + b*az + c*rg**2 + d*az+ e
                funcbins = a*azbins**2 + b*azbins + c*rgbins**2 + d*rgbins+ e
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,f*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,3] = np.arange((new_cols)) - jbeg_emp
                G[:,4] = 1
                G[:,5] = elevi


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            if (ivar==0 and nfit ==1):
                G=np.zeros((len(data),7))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = rg**2
                G[:,3] = rg
                G[:,4] = 1
                G[:,5] = topo_clean
                G[:,6] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]
                print ('Remove ramp %f az**2 + %f az  + %f r**2 + %f r + %f + %f z + %f z**2  for date: %i'%(a,b,c,d,e,f,g,idates[l]))

                # plot phase/elev
                funct = a*az**2 + b*az + c*rg**2 + d*rg+ e
                funcbins = a*azbins**2 + b*azbins + c*rgbins**2 + d*rgbins+ e
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,f*x+g*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,3] = np.arange((new_cols)) - jbeg_emp
                G[:,4] = 1
                G[:,5] = elevi
                G[:,6] = elevi**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit ==0):
                G=np.zeros((len(data),7))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = rg**2
                G[:,3] = rg
                G[:,4] = 1
                G[:,5] = topo_clean
                G[:,6] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g=pars[6]
                print ('Remove ramp %f az**2 + %f az  + %f r**2 + %f r + %f + %f z + %f az*z for date: %i'%(a,b,c,d,e,f,g,idates[l]))

                # plot phase/elev
                funct = a*az**2 + b*az + c*rg**2 + d*rg+ e + g*topo_clean*az
                funcbins = a*azbins**2 + b*azbins + c*rgbins**2 + d*rgbins+ e + g*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,f*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                G[:,4] = 1
                G[:,5] = elevi
                G[:,6] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,3] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,6] *= (i - ibeg_emp)


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),8))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = rg**2
                G[:,3] = rg
                G[:,4] = 1
                G[:,5] = topo_clean*az
                G[:,6] = topo_clean
                G[:,7] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g=pars[6]; h=pars[7]
                print ('Remove ramp %f az**2 + %f az  + %f r**2 + %f r + %f +  %f az*z + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,g,h,idates[l]))

                # plot phase/elev
                funct = a*az**2 + b*az + c*rg**2 + d*rg+ e + f*topo_clean*az
                funcbins = a*azbins**2 + b*azbins + c*rgbins**2 + d*rgbins+ e + f*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,g*x+h*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),8))
                G[:,4] = 1
                G[:,5] = elevi
                G[:,6] = elevi
                G[:,7] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,3] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,5] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)

      elif order==8:
        if arguments["--topofile"] is None:
            G=np.zeros((len(data),6))
            G[:,0] = az**3
            G[:,1] = az**2
            G[:,2] = az
            G[:,3] = rg**2
            G[:,4] = rg
            G[:,5] = 1

            # ramp inversion
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
            print ('Remove ramp %f az**3 + %f az**2  + %f az + %f r**2 + %f r + %f for date: %i'%(a,b,c,d,e,f,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),6))
            for i in range(new_lines):
                G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**3
                G[i*new_cols:(i+1)*new_cols,1] = (i - ibeg_emp)**2
                G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                G[i*new_cols:(i+1)*new_cols,3] = (np.arange((new_cols)) - jbeg_emp)**2
                G[i*new_cols:(i+1)*new_cols,4] = (np.arange((new_cols)) - jbeg_emp)
            G[:,5] = 1


            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),7))
                G[:,0] = az**3
                G[:,1] = az**2
                G[:,2] = az
                G[:,3] = rg**2
                G[:,4] = rg
                G[:,5] = 1
                G[:,6] = topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]
                print ('Remove ramp %f az**3 + %f az**2  + %f az + %f r**2 + %f r + %f + %f z for date: %i'%(a,b,c,d,e,f,g,idates[l]))

                # plot phase/elev
                funct = a*az**3 + b*az**2 + c*az + d*rg**2 + e*rg+ f
                funcbins = a*azbins**3 + b*azbins**2 + c*azbins + d*rgbins**2 + e*rgbins+ f
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,g*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**3
                    G[i*new_cols:(i+1)*new_cols,1] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,3] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,4] = np.arange((new_cols)) - jbeg_emp
                G[:,5] = 1
                G[:,6] = elevi


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            if (ivar==0 and nfit==1):
                G=np.zeros((len(data),8))
                G[:,0] = az**3
                G[:,1] = az**2
                G[:,2] = az
                G[:,3] = rg**2
                G[:,4] = rg
                G[:,5] = 1
                G[:,6] = topo_clean
                G[:,7] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]; h = pars[7]
                print ('Remove ramp %f az**3 + %f az**2  + %f az + %f r**2 + %f r + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,g,h,idates[l]))

                # plot phase/elev
                funct = a*az**3 + b*az**2 + c*az + d*rg**2 + e*rg+ f
                funcbins = a*azbins**3 + b*azbins**2 + c*azbins + d*rgbins**2 + e*rgbins+ f
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,g*x+h*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),8))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**3
                    G[i*new_cols:(i+1)*new_cols,1] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,3] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,4] = np.arange((new_cols)) - jbeg_emp
                G[:,5] = 1
                G[:,6] = elevi
                G[:,7] = elevi**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)


            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),8))
                G[:,0] = az**3
                G[:,1] = az**2
                G[:,2] = az
                G[:,3] = rg**2
                G[:,4] = rg
                G[:,5] = 1
                G[:,6] = topo_clean
                G[:,7] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]; h=pars[7]
                print ('Remove ramp %f az**3 + %f az**2  + %f az + %f r**2 + %f r + %f + %f z + %f z*az for date: %i'%(a,b,c,d,e,f,g,h,idates[l]))

                # plot phase/elev
                funct = a*az**3 + b*az**2 + c*az + d*rg**2 + e*rg + f + h*topo_clean*az
                funcbins = a*azbins**3 + b*azbins**2 + c*azbins + d*rgbins**2 + e*rgbins+ f + h*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,g*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),8))
                G[:,5] = 1
                G[:,6] = elevi
                G[:,7] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**3
                    G[i*new_cols:(i+1)*new_cols,1] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,3] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,4] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,7] *= (i - ibeg_emp)


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),10))
                G[:,0] = az**3
                G[:,1] = az**2
                G[:,2] = az
                G[:,3] = rg**2
                G[:,4] = rg
                G[:,5] = 1
                G[:,6] = topo_clean*az
                G[:,7] = topo_clean
                G[:,8] = topo_clean**2
                G[:,9] = (topo_clean*az)**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]; h=pars[7]; i=pars[8]; k=pars[9]
                print ('Remove ramp %f az**3 + %f az**2  + %f az + %f r**2 + %f r + %f z*az + %f + %f z + %f z**2 + %f (z*az)**2 for date: %i'%(a,b,c,d,e,f,g,h,i,k,idates[l]))

                # plot phase/elev
                funct = a*az**3 + b*az**2 + c*az + d*rg**2 + e*rg + f + g*topo_clean*az + k*(topo_clean*az)**2
                funcbins = a*azbins**3 + b*azbins**2 + c*azbins + d*rgbins**2 + e*rgbins+ f + g*topobins*azbins + k*(topobins*azbins)**2
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,h*x+i*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),10))
                G[:,5] = 1
                G[:,6] = elevi
                G[:,7] = elevi
                G[:,8] = elevi**2
                G[:,9] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**3
                    G[i*new_cols:(i+1)*new_cols,1] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,3] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,4] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,6] *= (i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,9] *= (i - ibeg_emp)**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)

      elif order==9:
        if arguments["--topofile"] is None:
            G=np.zeros((len(data),5))
            G[:,0] = rg
            G[:,1] = az
            G[:,2] = (rg*az)**2
            G[:,3] = rg*az
            G[:,4] = 1

            # ramp inversion
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
            print ('Remove ramp %f r + %f az  + %f r*az**2 + %f r*az + %f for date: %i'%(a,b,c,d,e,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),5))
            for i in range(new_lines):
                G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                G[i*new_cols:(i+1)*new_cols,2] = ((i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp))**2
                G[i*new_cols:(i+1)*new_cols,3] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
            G[:,4] = 1

            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),6))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = (rg*az)**2
                G[:,3] = rg*az
                G[:,4] = 1
                G[:,5] = topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]

                print ('Remove ramp %f r + %f az  + %f (r*az)**2 + %f r*az + %f + %f z for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funcbins = a*rgbins+ b*azbins + c*(azbins*rgbins)**2 + d*azbins*rgbins+ e
                funct = a*rg+ b*az + c*(az*rg)**2 + d*az*rg+ e
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = ((i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp))**2
                    G[i*new_cols:(i+1)*new_cols,3] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
                G[:,4] = 1
                G[:,5] = elevi

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            if (ivar==0 and nfit==1):
                G=np.zeros((len(data),7))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = (rg*az)**2
                G[:,3] = rg*az
                G[:,4] = 1
                G[:,5] = topo_clean
                G[:,6] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]

                print ('Remove ramp %f r + %f az  + %f (r*az)**2 + %f r*az + %f + %f z + %f z**2  for date: %i'%(a,b,c,d,e,f,g,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c*(rg*az)**2 + d*rg*az+ e
                funcbins = a*rgbins+ b*azbins + c*(azbins*rgbins)**2 + d*azbins*rgbins+ e
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,f*x+g*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = ((i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp))**2
                    G[i*new_cols:(i+1)*new_cols,3] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
                G[:,4] = 1
                G[:,5] = elevi
                G[:,6] = elevi**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),7))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = (rg*az)**2
                G[:,3] = rg*az
                G[:,4] = 1
                G[:,5] = topo_clean
                G[:,6] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5] ; g = pars[6]

                print ('Remove ramp %f r + %f az  + %f (r*az)**2 + %f r*az + %f + %f z + %f az*z for date: %i'%(a,b,c,d,e,f,g,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c*(rg*az)**2 + d*rg*az+ e + g*topo_clean*az
                funcbins = a*rgbins+ b*azbins + c*(azbins*rgbins)**2 + d*azbins*rgbins+ e + g*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,f*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                G[:,4] = 1
                G[:,5] = elevi
                G[:,6] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = ((i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp))**2
                    G[i*new_cols:(i+1)*new_cols,3] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
                    G[i*new_cols:(i+1)*new_cols,6] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),8))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = (rg*az)**2
                G[:,3] = rg*az
                G[:,4] = 1
                G[:,5] = topo_clean*az
                G[:,6] = topo_clean
                G[:,7] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5] ; g = pars[6]; h=pars[7]

                print ('Remove ramp %f r + %f az  + %f (r*az)**2 + %f r*az + %f + %f az*z + %f z + %f z**2  for date: %i'%(a,b,c,d,e,f,g,h,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c*(rg*az)**2 + d*rg*az+ e + f*topo_clean*az
                funcbins = a*rgbins+ b*azbins + c*(azbins*rgbins)**2 + d*azbins*rgbins+ e + f*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,g*x+h*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),8))
                G[:,4] = 1
                G[:,5] = elevi
                G[:,6] = elevi
                G[:,7] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = ((i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp))**2
                    G[i*new_cols:(i+1)*new_cols,3] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
                    G[i*new_cols:(i+1)*new_cols,5] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)

      rms = np.sqrt(np.nanmean(res**2))
      logger.info('RMS dates %i: %f'%(idates[l], rms))
      flata = los.reshape(new_lines,new_cols) - ramp - topo
      
      try:
         del G; del los
      except:
         pass
      return ramp, flata, topo, rms
 

def empirical_cor(l, ibeg_emp, iend_emp, mintopo, maxtopo, topofile, perc_los, threshold_rms, threshold_mask, emp_sampling, lin_start, lin_end, col_start, col_end): 
  """
  Function that preapare and run empirical estimaton for each interferogram kk
  """

  global fig_dphi
  global maps, models, elev_map, aspect_map, rms_map

  # first clean los
  map_temp = as_strided(maps[:,:,l]) - as_strided(models[:,:,l])

  # no estimation on the ref image set to zero 
  if np.nansum(maps[:,:,l]) != 0:

    maxlos,minlos=np.nanpercentile(map_temp[ibeg_emp:iend_emp,jbeg_emp:jend_emp],float(perc_los)),np.nanpercentile(map_temp[ibeg_emp:iend_emp,jbeg_emp:jend_emp],100-float(perc_los))
    logger.debug('Set Max-Min LOS for empirical estimation: {0}-{1}'.format(maxlos,minlos))
    kk = np.nonzero(np.logical_or(map_temp==0.,np.logical_or((map_temp>maxlos),(map_temp<minlos))))
    map_temp[kk] = float('NaN')

    itemp = ibeg_emp
    for lign in range(ibeg_emp,iend_emp,10):
        if np.isnan(np.nanmean(maps[lign:lign+10,:,l])):
            itemp = lign
        else:
            break
    logger.debug('Begining of the image: {}'.format(itemp))

    if topofile is not None:
        ax_dphi = fig_dphi.add_subplot(4,int(N/4)+1,l+1)
    else:
        ax_dphi = None

    logger.debug('Threshold RMS: {}'.format(float(threshold_rms)))

    # selection pixels
    index = np.nonzero(np.logical_and(elev<maxtopo,
        np.logical_and(elev>mintopo,
            np.logical_and(mask_flat>float(threshold_mask),
            np.logical_and(~np.isnan(map_temp),
                np.logical_and(~np.isnan(rmsmap),
                np.logical_and(~np.isnan(elev),
                np.logical_and(rmsmap<float(threshold_rms),
                np.logical_and(rmsmap>1.e-6,
                np.logical_and(~np.isnan(map_temp),
                np.logical_and(pix_az>ibeg_emp,
                np.logical_and(pix_az<iend_emp,
                np.logical_and(pix_rg>jbeg_emp,
                np.logical_and(pix_rg<jend_emp, 
                    aspect>0.,
                    ))))))))
                ))))))

    # extract coordinates for estimation
    temp = np.array(index).T
    x = temp[:,0]; y = temp[:,1]
    # clean maps
    los_clean = map_temp[index].flatten()
    topo_clean = elev[index].flatten()
    rms_clean = rmsmap[index].flatten()
    
    logger.debug('Number of points for empirical estimation: {}'.format(len(los_clean)))
    if len(los_clean) < 1:
      logger.critical('No points left for empirical estimation. Exit!')
      logger.critical('threshold RMS: {0}, threshold Mask: {1}, Min-Max LOS: {2}-{3}, Min-Max topo: {4}-{5}, lines: {6}-{7}, \
        cols: {8}- {9}'.format(float(threshold_rms),float(threshold_mask),minlos,maxlos,mintopo,maxtopo,ibeg_emp,iend_emp,jbeg_emp,jend_emp))
      sys.exit()

    # print itemp, iend_emp
    if flat>5 and iend_emp-itemp < .6*(iend_emp-ibeg_emp):
        logger.warning('Image too short in comparison to master, set flat to 5')
        temp_flat=5
    else:
        temp_flat=flat

    if ivar>0 and iend_emp-itemp < .6*(iend_emp-ibeg_emp):
      logger.warning('Image too short in comparison to master, set ivar to 0')
      ivar_temp=0
      nfit_temp=0
    else:
      ivar_temp=ivar
      nfit_temp=nfit

    # call ramp estim
    los = as_strided(maps[:,:,l]).flatten()
    samp = int(emp_sampling)

    map_ramp, map_flata, map_topo, rmsi = estim_ramp(los,los_clean[::samp],topo_clean[::samp],x[::samp],\
      y[::samp],temp_flat,rms_clean[::samp],nfit_temp, ivar_temp, l, ax_dphi)

    if (lin_start is not None) and (lin_end is not None):
        indexref = np.nonzero(np.logical_and(elev<maxtopo,
        np.logical_and(elev>mintopo,
            np.logical_and(mask_flat>float(threshold_mask),
            np.logical_and(~np.isnan(map_temp),
                np.logical_and(~np.isnan(rmsmap),
                np.logical_and(~np.isnan(elev),
                np.logical_and(rmsmap<float(threshold_rms),
                np.logical_and(rmsmap>1.e-6,
                np.logical_and(~np.isnan(map_temp),
                np.logical_and(pix_az>lin_start,
                np.logical_and(pix_az<lin_end,
                np.logical_and(pix_rg>col_start,
                np.logical_and(pix_rg<col_end, 
                    aspect>0.,
                ))))))))))))
                ))
        
        if len(indexref[0]) == 0:
             logger.warning('Ref zone is empty! Re-define --ref_zone argument. Exit!')
             sys.exit()

        ## Set data minus temporal model to zero in the ref area
        zone = as_strided(map_flata[:,:] - models[:,:,l])
        los_ref2 = zone[indexref].flatten()
        rms_ref = rmsmap[indexref].flatten()
        amp_ref = 1./rms_ref
        amp_ref = amp_ref/np.nanmax(amp_ref)
        # weigth avera of the phase
        cst = np.nansum(los_ref2*amp_ref) / np.nansum(amp_ref)
        logger.info('Re-estimation of a constant within lines {0}-{1} and cols {2}-{3}'.format(lin_start,lin_end,col_start,col_end))
        logger.info('Average phase within ref area: {0}:'.format(cst))
        if np.isnan(cst):
          cst = 0.
        map_ramp, map_flata = map_ramp + cst, map_flata - cst
        del zone
      
  else:
    map_flata = np.copy(maps[:,:,l])
    map_ramp, map_topo  = np.zeros(np.shape(map_flata)), np.zeros(np.shape(map_flata))
    rmsi = 1

  # set ramp to NaN to have ramp of the size of the images
  kk = np.nonzero(np.isnan(map_flata))
  ramp = as_strided(map_ramp)
  ramp[kk] = float('NaN')
  topo = as_strided(map_topo)
  topo[kk] = float('NaN')
  del ramp, topo, map_temp
 
  if arguments["--topofile"] is not None: 
    return map_flata, map_topo, rmsi 
  else:
    return map_flata, rmsi 

def temporal_decomp(pix, disp, sigma, cond, ineq, eguality):
    j = pix  % (new_cols)
    i = int(pix/(new_cols))

    # Initialisation
    mdisp=np.ones((N), dtype=np.float32)*float('NaN')
    mlin=np.ones((N), dtype=np.float32)*float('NaN')
    mseas=np.ones((N), dtype=np.float32)*float('NaN')
    mvect=np.ones((N), dtype=np.float32)*float('NaN')
    k = np.flatnonzero(~np.isnan(disp)) # invers of isnan
    # do not take into account NaN data
    kk = len(k)
    tabx = dates[k]
    taby = disp[k]

    # Inisilize m to zero
    m = np.ones((M), dtype=np.float32)*float('NaN')
    sigmam = np.ones((M), dtype=np.float32)*float('NaN')

    if kk > N/6:
        G=np.zeros((kk,M), dtype=np.float32)
        # Build G family of function k1(t),k2(t),...,kn(t): #
        #                                                   #
        #           |k1(0) .. kM(0)|                        #
        # Gfamily = |k1(1) .. kM(1)|                        #
        #           |..    ..  ..  |                        #
        #           |k1(N) .. kM(N)|                        #
        #                                                   #
        for l in range((Mbasis)):
            G[:,l]=basis[l].g(tabx)
        for l in range((Mker)):
            G[:,Mbasis+l]=kernels[l].g(k)
        
        # inversion
        m,sigmam = consInvert(G,taby,sigma[k],cond=cond,ineq=ineq,eguality=eguality)

        # forward model in original order
        mdisp[k] = np.dot(G,m)

    return m, sigmam, mdisp 

# initialization
models = np.zeros((new_lines,new_cols,N), dtype=np.float32)

# prepare flatten maps
if arguments["--topofile"] is not None:
    maps_topo = np.zeros((new_lines,new_cols,N),dtype=np.float32)

for ii in range(int(arguments["--niter"])):
    print()
    print('---------------')
    print('iteration: {}'.format(ii+1))
    print('---------------')

    #############################
    # SPATIAL ITERATION N  ######
    #############################

    rms = np.zeros((N),dtype=np.float32)
    pix_az, pix_rg = np.indices((new_lines,new_cols))
    # if radar file just initialise figure
    if arguments["--topofile"] is not None:
      nfigure +=1
      fig_dphi = plt.figure(nfigure,figsize=(14,10))
    
    # if iteration = 0 or spatialiter==yes, then spatial estimation
    if (ii==0) or (arguments["--spatialiter"]=='yes') :

      print() 
      #########################################
      print('---------------')
      print('Empirical estimations')
      print('---------------')
      #########################################
      print()
    
      for l in range((N)):
          if arguments["--topofile"] is not None:
            maps[:,:,l], maps_topo[:,:,l], rms[l] = empirical_cor(l, ibeg_emp, iend_emp, mintopo, maxtopo, arguments["--topofile"], arguments["--perc_los"], arguments["--threshold_rms"], arguments["--threshold_mask"], arguments["--emp_sampling"], lin_start, lin_end, col_start, col_end)
          else:
            maps[:,:,l], rms[l] = empirical_cor(l, ibeg_emp, iend_emp, mintopo, maxtopo, arguments["--topofile"], arguments["--perc_los"], arguments["--threshold_rms"], arguments["--threshold_mask"], arguments["--emp_sampling"], lin_start, lin_end, col_start, col_end)

      if N < 30:
        # plot corrected ts
        nfigure +=1
        plot_displacement_maps(maps, idates, nfigure=nfigure, cmap=cmap, plot=plot, title='Corrected time series maps from empirical estimations', filename='maps_flat.eps')

        if arguments["--topofile"] is not None:
            fig_dphi.savefig('phase-topo.eps', format='EPS',dpi=150)
            nfigure +=1
            plot_displacement_maps(maps_topo, idates, nfigure=nfigure, cmap=cmap, plot=plot, title='Time series RAMPS+TROPO', filename='maps_models_tropo.eps') 
            del maps_topo

    if plot=='yes':
        plt.show()
    plt.close('all')

    # save rms
    if (arguments["--aps"] is None and ii==0):
        # aps from rms
        logger.info('Use RMS empirical estimations as input APS for time decomposition')
        in_aps = np.copy(rms)
        logger.info('Set very low values to the 2 percentile to avoid overweighting...')
        # scale between 0 and 1 
        maxaps = np.nanmax(in_aps)
        in_aps = in_aps/maxaps
        in_aps[imref] = 0.5
        min_aps= np.nanpercentile(in_aps,2)
        index = np.flatnonzero(in_aps<min_aps)
        in_aps[index] = min_aps
        in_sigma = in_aps * in_rms
        np.savetxt('rms_empcor.txt', in_aps.T)
        del rms

    #######################################################
    # Save new cubes
    #######################################################

    # create new cube
    if flat>0:
        logger.info('Save flatten time series cube: {}'.format('disp_cumule_flat'))
        fid = open('disp_cumule_flat', 'wb')
        maps[:,:,:].astype('float32').tofile(fid)
        in_hdr = arguments["--cube"] + '.hdr'
        arguments["--cube"] = 'disp_cumule_flat' # update depl_cumul
        out_hdr = 'disp_cumule_flat.hdr'
        try:
            shutil.copy(in_hdr, out_hdr)
        except:
            pass
        fid.close()
        del fid

    ########################
    # TEMPORAL ITERATION N #
    ########################

    print() 
    #########################################
    print('---------------')
    print('Time Decomposition')
    print('---------------')
    #########################################
    print()

    logger.info('Input uncertainties: {}'.format(in_sigma))

    # reiinitialize maps models
    models = np.zeros((new_lines,new_cols,N),dtype=np.float32)

    with TimeIt():
        for pix in range(0,(new_lines)*(new_cols),int(arguments["--sampling"])):
              j = pix  % (new_cols)
              i = int(pix/(new_cols))
              if ((i % 20) == 0) and (j==0):
                  logger.info('Processing line: {} --- {} seconds ---'.format(i,time.time() - start_time))
              disp = as_strided(maps[i,j,:])
              m, sigmam, models[i,j,:]  = temporal_decomp(pix, disp, in_sigma, arguments['--cond'], arguments['--ineq'], eguality)
            
              # save m
              for l in range((Mbasis)):
                  basis[l].m[i,j] = m[l]
                  basis[l].sigmam[i,j] = sigmam[l]

              for l in range((Mker)):
                  kernels[l].m[i,j] = m[Mbasis+l]
                  kernels[l].sigmam[i,j] = sigmam[Mbasis+l]

              del m, sigmam

    # compute RMSE
    # remove outiliers
    index = np.logical_or(models>9999., models<-9999)
    models[index] = 0.
    squared_diff = (np.nan_to_num(maps,nan=0) - np.nan_to_num(models, nan=0))**2
    res = np.sqrt(np.nanmean(squared_diff, axis=(0,1))**2)  

    # remove low res to avoid over-fitting in next iter
    min_res= np.nanpercentile(res,2)
    index = np.flatnonzero(res < min_res)
    res[index] = min_res

    print('Dates      Residuals  ')
    for l in range(N):
        print (idates[l], res[l])
    np.savetxt('aps_{}.txt'.format(ii), res.T, fmt=('%.6f'))
    # set apsf is yes for next iteration
    arguments["--aps"] == 'yes'
    # update aps for next iterations taking into account in_aps and the residues of the last iteration
    in_sigma = res * in_aps * in_rms

#######################################################
# Save functions in binary file
#######################################################

if arguments["--geotiff"] is not None:
    for l in range((Mbasis)):
        outname = '{}_coeff.tif'.format(basis[l].reduction)
        logger.info('Save: {}'.format(outname))
        ds = driver.Create(outname, new_cols, new_lines, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(basis[l].m)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

        outname = '{}_sigcoeff.tif'.format(basis[l].reduction)
        logger.info('Save: {}'.format(outname))
        ds = driver.Create(outname, new_cols, new_lines, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(basis[l].sigmam)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

    for l in range((Mker)):
        outname = '{}_coeff.tif'.format(kernels[l].reduction)
        logger.info('Save: {}'.format(outname))
        ds = driver.Create(outname, new_cols, new_lines, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(kernels[l].m)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

        outname = '{}_sigcoeff.tif'.format(kernels[l].reduction)
        logger.info('Save: {}'.format(outname))
        ds = driver.Create(outname, new_cols, new_lines, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(kernels[l].sigmam)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

else:
    for l in range((Mbasis)):
        outname = '{}_coeff.r4'.format(basis[l].reduction)
        logger.info('Save: {}'.format(outname))
        fid = open(outname, 'wb')
        basis[l].m.astype('float32').tofile(fid)
        fid.close()
        outname = '{}_sigcoeff.r4'.format(basis[l].reduction)
        logger.info('Save: {}'.format(outname))
        fid = open(outname, 'wb')
        basis[l].sigmam.astype('float32').tofile(fid)
        fid.close()
    for l in range((Mker)):
        outname = '{}_coeff.r4'.format(kernels[l].reduction)
        logger.info('Save: {}'.format(outname))
        fid = open('{}_coeff.r4'.format(kernels[l].reduction), 'wb')
        kernels[l].m.astype('float32').tofile(fid)
        fid.close()
        outname = '{}_sigcoeff.r4'.format(kernels[l].reduction)
        logger.info('Save: {}'.format(outname))
        fid = open('{}_sigcoeff.r4'.format(kernels[l].reduction), 'wb')
        kernels[l].sigmam.astype('float32').tofile(fid)
        fid.close()

#######################################################
# Save models
#######################################################

if arguments["--fulloutput"]=='yes':
        logger.info('Save  time series cube model: {}'.format('disp_cumule_models'))
        fid = open('depl_cumule_models', 'wb')
        maps[:,:,:].astype('float32').tofile(fid)
        in_hdr = arguments["--cube"] + '.hdr'
        out_hdr = 'depl_cumule_models.hdr'
        try:
            shutil.copy(in_hdr, out_hdr)
        except:
            pass
        fid.close()
        del fid

#######################################################
# Plot models and residuals
#######################################################

if N < 30:
  nfigure +=1
  figclr = plt.figure(nfigure)
  # plot color map
  ax = figclr.add_subplot(1,1,1)
  vmax = np.nanpercentile(maps[:,:,:],99.)
  vmin = np.nanpercentile(maps[:,:,:],1.)  
  cax = ax.imshow(maps[:,:,-1],cmap=cmap,vmax=vmax,vmin=vmin)
  plt.setp( ax.get_xticklabels(), visible=False)
  cbar = figclr.colorbar(cax, orientation='horizontal',aspect=5)
  figclr.savefig('colorscale.eps', format='EPS',dpi=150) 

  nfigure +=1
  plot_displacement_maps(models, idates, nfigure=nfigure, cmap=cmap, plot=plot, title='Time series models', filename='maps-time-models.eps')
  nfigure +=1
  plot_displacement_maps(maps-models, idates, nfigure=nfigure, cmap=cmap, plot=plot, title='Time series residuals', filename='maps-residuals.eps')

# clean memory
del maps, models

#######################################################
# Compute Amplitude and phase seasonal
#######################################################

if arguments["--seasonal"]  == 'yes':
    cosine = as_strided(basis[indexseas].m)
    sine = as_strided(basis[indexseas+1].m)
    amp = np.sqrt(cosine**2+sine**2)
    phi = np.arctan2(sine,cosine)

    sigcosine = as_strided(basis[indexseas].sigmam)
    sigsine = as_strided(basis[indexseas+1].sigmam)
    sigamp = np.sqrt(sigcosine**2+sigsine**2)
    sigphi = (sigcosine*abs(sine)+sigsine*abs(cosine))/(sigcosine**2+sigsine**2)

    if arguments["--geotiff"] is not None:
        logger.info('Save: {}'.format(outname))
        ds = driver.Create('ampwt_coeff.tif', new_cols, new_lines, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(amp)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

        logger.info('Save: {}'.format('ampwt_sigcoeff.tif'))
        ds = driver.Create('ampwt_sigcoeff.tif', new_cols, new_lines, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(sigamp)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

    else:
        logger.info('Save: {}'.format('ampwt_coeff.r4'))
        fid = open('ampwt_coeff.r4', 'wb')
        amp.astype('float32').tofile(fid)
        fid.close()

        logger.info('Save: {}'.format('ampwt_sigcoeff.r4'))
        fid = open('ampwt_sigcoeff.r4', 'wb')
        sigamp.astype('float32').tofile(fid)
        fid.close()

    if arguments["--geotiff"] is not None:
        logger.info('Save: {}'.format('phiwt_coeff.tif'))
        ds = driver.Create('phiwt_coeff.tif', new_cols, new_lines, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(phi)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

        logger.info('Save: {}'.format('phiwt_sigcoeff.tif'))
        ds = driver.Create('phiwt_sigcoeff.tif', new_cols, new_lines, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(sigphi)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

    else:
        logger.info('Save: {}'.format('phiwt_coeff.r4'))
        fid = open('phiwt_coeff.r4', 'wb')
        phi.astype('float32').tofile(fid)
        fid.close()

        logger.info('Save: {}'.format('phiwt_sigcoeff.r4'))
        fid = open('phiwt_sigcoeff.r4', 'wb')
        sigphi.astype('float32').tofile(fid)
        fid.close()

if arguments["--seasonal_increase"]  == 'yes':
    cosine = as_strided(basis[indexseast].m)
    sine = as_strided(basis[indexseast+1].m)
    amp = np.sqrt(cosine**2+sine**2)
    phi = np.arctan2(sine,cosine)

    sigcosine = as_strided(basis[indexseast].sigmam)
    sigsine = as_strided(basis[indexseast+1].sigmam)
    sigamp = np.sqrt(sigcosine**2+sigsine**2)
    sigphi = (sigcosine*abs(sine)+sigsine*abs(cosine))/(sigcosine**2+sigsine**2)

    if arguments["--geotiff"] is not None:
        logger.info('Save: {}'.format(outname))
        ds = driver.Create('ampwt_increase_coeff.tif', new_cols, new_lines, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(amp)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

        logger.info('Save: {}'.format('ampwt_increase_sigcoeff.tif'))
        ds = driver.Create('ampwt_increase_sigcoeff.tif', new_cols, new_lines, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(sigamp)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

    else:
        logger.info('Save: {}'.format('ampwt_increase_coeff.r4'))
        fid = open('ampwt_increase_coeff.r4', 'wb')
        amp.astype('float32').tofile(fid)
        fid.close()

        logger.info('Save: {}'.format('ampwt_increase_sigcoeff.r4'))
        fid = open('ampwt_increase_sigcoeff.r4', 'wb')
        sigamp.astype('float32').tofile(fid)
        fid.close()

    if arguments["--geotiff"] is not None:
        logger.info('Save: {}'.format('phiwt_increase_coeff.tif'))
        ds = driver.Create('phiwt_increase_coeff.tif', new_cols, new_lines, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(phi)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

        logger.info('Save: {}'.format('phiwt_increase_sigcoeff.tif'))
        ds = driver.Create('phiwt_increase_sigcoeff.tif', new_cols, new_lines, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(sigphi)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

    else:
        logger.info('Save: {}'.format('phiwt_increase_coeff.r4'))
        fid = open('phiwt_increase_coeff.r4', 'wb')
        phi.astype('float32').tofile(fid)
        fid.close()

        logger.info('Save: {}'.format('phiwt_increase_sigcoeff.r4'))
        fid = open('phiwt_increase_sigcoeff.r4', 'wb')
        sigphi.astype('float32').tofile(fid)
        fid.close()

if arguments["--semianual"] == 'yes':
    cosine = as_strided(basis[indexsemi].m)
    sine = as_strided(basis[indexsemi+1].m)
    amp = np.sqrt(cosine**2+sine**2)
    phi = np.arctan2(sine,cosine)

    sigcosine = as_strided(basis[indexseas].sigmam)
    sigsine = as_strided(basis[indexseas+1].sigmam)
    sigamp = np.sqrt(sigcosine**2+sigsine**2)
    sigphi = (sigcosine*abs(sine)+sigsine*abs(cosine))/(sigcosine**2+sigsine**2)

    if arguments["--geotiff"] is not None:
        logger.info('Save: {}'.format('amp_simiwt_coeff.tif'))
        ds = driver.Create('amp_simiwt_coeff.tif', new_cols, new_lines, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(amp)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

        logger.info('Save: {}'.format('amp_simiwt_sigcoeff.tif'))
        ds = driver.Create('amp_simiwt_sigcoeff.tif', new_cols, new_lines, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(sigamp)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

    else:
        logger.info('Save: {}'.format('amp_simiwt_coeff.r4'))
        fid = open('amp_simiwt_coeff.r4', 'wb')
        amp.astype('float32').tofile(fid)
        fid.close()

        logger.info('Save: {}'.format('amp_simiwt_sigcoeff.r4'))
        fid = open('amp_simiwt_sigcoeff.r4', 'wb')
        sigamp.astype('float32').tofile(fid)
        fid.close()

    if arguments["--geotiff"] is not None:
        logger.info('Save: {}'.format('phi_simiwt_coeff.tif'))
        ds = driver.Create('phi_simiwt_coeff.tif', new_cols, new_lines, 1, gdal.GDT_Float32)
        band.WriteArray(phi)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

        logger.info('Save: {}'.format('phi_simiwt_sigcoeff.tif'))
        ds = driver.Create('phi_simiwt_sigcoeff.tif', new_cols, new_lines, 1, gdal.GDT_Float32)
        band.WriteArray(sigphi)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

    else:
        logger.info('Save: {}'.format('phi_simiwt_coeff.r4'))
        fid = open('phi_simiwt_coeff.r4', 'wb')
        phi.astype('float32').tofile(fid)
        fid.close()

        logger.info('Save: {}'.format('phi_simiwt_sigcoeff.r4'))
        fid = open('phi_simiwt_sigcoeff.r4', 'wb')
        sigphi.astype('float32').tofile(fid)
        fid.close()

if arguments["--bianual"] == 'yes':
    cosine = as_strided(basis[indexbi].m)
    sine = as_strided(basis[indexbi+1].m)
    amp = np.sqrt(cosine**2+sine**2)
    phi = np.arctan2(sine,cosine)

    sigcosine = as_strided(basis[indexbi].sigmam)
    sigsine = as_strided(basis[indexbi+1].sigmam)
    sigamp = np.sqrt(sigcosine**2+sigsine**2)
    sigphi = (sigcosine*abs(sine)+sigsine*abs(cosine))/(sigcosine**2+sigsine**2)

    if arguments["--geotiff"] is not None:
        logger.info('Save: {}'.format('amp_biwt_coeff.tif'))
        ds = driver.Create('ampw.5t_coeff.tif', new_cols, new_lines, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(amp)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

        logger.info('Save: {}'.format('amp_biwt_sigcoeff.tif'))
        ds = driver.Create('ampw.5t_sigcoeff.tif', new_cols, new_lines, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(sigamp)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

    else:
        logger.info('Save: {}'.format('amp_biwt_coeff.r4'))
        fid = open('ampw.5t_coeff.r4', 'wb')
        amp.astype('float32').tofile(fid)
        fid.close()

        logger.info('Save: {}'.format('amp_biwt_sigcoeff.r4'))
        fid = open('ampw.5t_sigcoeff.r4', 'wb')
        sigamp.astype('float32').tofile(fid)
        fid.close()

    if arguments["--geotiff"] is not None:
        logger.info('Save: {}'.format('phi_biwt_coeff.tif'))
        ds = driver.Create('phiw.5t_coeff.tif', new_cols, new_lines, 1, gdal.GDT_Float32)
        band.WriteArray(phi)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

        logger.info('Save: {}'.format('phi_biwt_sigcoeff.tif'))
        ds = driver.Create('phiw.5t_sigcoeff.tif', new_cols, new_lines, 1, gdal.GDT_Float32)
        band.WriteArray(sigphi)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

    else:
        logger.info('Save: {}'.format('phi_biwt_coeff.r4'))
        fid = open('phiw.5t_coeff.r4', 'wb')
        phi.astype('float32').tofile(fid)
        fid.close()

        logger.info('Save: {}'.format('phi_biwt_sigcoeff.r4'))
        fid = open('phiw.5t_sigcoeff.r4', 'wb')
        sigphi.astype('float32').tofile(fid)
        fid.close()

#######################################################
# Plot
#######################################################

nfigure +=1
fig = plt.figure(nfigure,figsize=(14,10))

for l in range(Mbasis):
    vmax = np.abs([np.nanpercentile(basis[l].m,98.),np.nanpercentile(basis[l].m,2.)]).max()
    vmin = -vmax

    ax = fig.add_subplot(1,M,l+1)
    cax = ax.imshow(basis[l].m,cmap=cmap,vmax=vmax,vmin=vmin)
    ax.set_title(basis[l].reduction)
    # add colorbar
    cbar = fig.colorbar(cax, orientation='vertical',shrink=0.2)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)

for l in range(Mker):
    vmax = np.abs([np.nanpercentile(kernels[l].m,98.),np.nanpercentile(kernels[l].m,2.)]).max()
    vmin = -vmax

    ax = fig.add_subplot(1,M,Mbasis+l+1)
    cax = ax.imshow(kernels[l].m,cmap=cmap,vmax=vmax,vmin=vmin)
    ax.set_title(kernels[l].reduction)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    cbar = fig.colorbar(cax, orientation='vertical',shrink=0.2)

plt.suptitle('Time series decomposition')

nfigure += 1
fig.tight_layout()
fig.savefig('inversion.eps', format='EPS',dpi=150)

if plot=='yes':
    plt.show()

