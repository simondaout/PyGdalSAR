#!/usr/bin/env python2
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Simon DAOUT (Oxford)
############################################

"""\
invers_disp2coef.py
-------------
Spatial and temporal inversions of the time series delay maps (used depl_cumule (BIP format) and images_retenues, output of invers_pixel), based on an iteration procedure.
At each iteration, (1) estimation of spatial ramps, (2) linear decomposition in time based on a library of temporal functions (linear, heaviside, logarithm, seasonal),
(3) estimation of RMS that will be then used as weight for the next iteration. Possibility to also to correct for a term proportional to the topography.

Usage: invers_disp2coef.py  [--cube=<path>] [--lectfile=<path>] [--list_images=<path>] [--aps=<path>] [--refstart=<value>] [--refend=<value>] [--interseismic=<yes/no>] [--threshold_rmsd=<value>] \
[--coseismic=<values>] [--postseismic=<values>]  [--seasonal=<yes/no>] [--slowslip=<values>] [--semianual=<yes/no>]  [--dem=<yes/no>] [--vector=<path>] \
[--flat=<0/1/2/3/4/5/6/7/8/9>] [--nfit=<0/1>] [--ivar=<0/1>] [--niter=<value>]  [--spatialiter=<yes/no>]  [--sampling=<value>] [--imref=<value>] [--mask=<path>] \
[--rampmask=<yes/no>] [--threshold_mask=<value>] [--scale_mask=<value>] [--topofile=<path>] [--aspect=<path>] [--perc_topo=<value>] [--perc_los=<value>] \
[--tempmask=<yes/no>] [--cond=<value>] [--ineq=<value>] [--rmspixel=<path>] [--threshold_rms=<path>] \
[--crop=<values>] [--fulloutput=<yes/no>] [--geotiff=<path>] [--plot=<yes/no>] \
[<ibeg>] [<iend>] [<jbeg>] [<jend>]

invers_disp2coef.py -h | --help

Options:
-h --help               Show this screen
--cube PATH             Path to displacement file [default: depl_cumul]
--lectfile PATH         Path to the lect.in file (output of invers_pixel) [default: lect.in]
--list_images PATH      Path to list images file made of 5 columns containing for each images 1) number 2) Doppler freq (not read) 3) date in YYYYMMDD format 4) numerical date 5) perpendicular baseline [default: images_retenues]
--refstart VALUE          Stating line number of the area where phase is set to zero [default: 0]
--refend VALUE            Ending line number of the area where phase is set to zero [default: 200]
--aps PATH              Path to the APS file giving an input error to each dates [default: No weigthing if no spatial estimation or misfit spatial estimation used as input uncertianties]
--rmspixel PATH         Path to the RMS map that gives an error for each pixel (e.g RMSpixel, output of invers_pixel) [default: None]
--threshold_rms VALUE   Threshold on rmsmap for spatial estimations [default: 1.]
--interseismic YES/NO   Add a linear function in the inversion
--threshold_rmsd VALUE  If interseismic = yes: first try inversion without coseismic and postseismic, if RMDS inversion > threshold_rmsd then add other basis functions [default: 1.]
--coseismic PATH        Add heaviside functions to the inversion, indicate coseismic time (e.g 2004.,2006.)
--postseismic PATH      Add logarithmic transients to each coseismic step, indicate characteristic time of the log function, must be a serie of values of the same lenght than coseismic (e.g 1.,1.). To not associate postseismic function to a give coseismic step, put None (e.g None,1.)
--slowslip   VALUE      Add slow-slip function in the inversion (as defined by Larson et al., 2004). Indicate median and characteristic time of the events (e.g. 2004.,1,2006,0.5), default: None
--vector PATH           Path to the vector text files containing a value for each dates [default: None]
--seasonal YES/NO       If yes, add seasonal terms in the inversion
--semianual YES/NO      If yes, add semianual terms in the inversion
--dem Yes/No            If yes, add term proportional to the perpendicular baseline in the inversion
--ivar VALUE            Define phase/elevation relationship: ivar=0 function of elevation, ivar=1 crossed function of azimuth and elevation
--nfit VALUE            Fit degree in azimuth or in elevation [0:linear (default), 1: quadratic]
--flat PATH             Remove a spatial ramp at each iteration.
0: ref frame [default], 1: range ramp ax+b , 2: azimutal ramp ay+b, 3: ax+by+c,
4: ax+by+cxy+d 5: ax**2+bx+cy+d, 6: ay**2+by+cx+d, 7: ay**2+by+cx**2+dx+e,
8: ay**2+by+cx**3+dx**2+ex+f, 9: ax+by+cxy**2+dxy+e
--niter VALUE           Number of iterations. At the first iteration, image uncertainties is given by aps file or misfit spatial iteration, while for the next itarations, uncertainties are equals to the global RMS of the previous iteration for each map [default: 1]
--spatialiter  YES/NO   If yes iterate the spatial estimations at each iterations (defined by niter) on the maps minus the temporal terms (ie. interseismic, coseismic...) [default: no]
--sampling VALUE        Downsampling factor [default: 1]
--imref VALUE           Reference image number [default: 1]
--mask PATH             Path to mask file in r4 or tif format for the spatial estimations (Keep only values > threshold_mask for ramp estimation).
--rampmask YES/NO       Remove a quadratic ramp in range a linear in azimuth on the mask before computing threshold [default: no]
--threshold_mask VALUE  Threshold on mask: take only > values (use scale factor for convenience) [default: 0]
--scale_mask  VALUE     Scale factor to apply on mask
--tempmask YES/NO       If yes, also use the spatial mask for the temporal inversion [default: no]
--topofile PATH         Path to topographic file in r4 or tif format. If not None, add a phase-elevation relationship in the saptial estimation.
--aspect PATH           Path to aspect file in r4 or tif format: take into account the slope orientation in the phase/topo relationship [default: None].
--perc_los VALUE        Percentile of hidden LOS pixel for the spatial estimations to clean outliers [default:98.]
--perc_topo VALUE       Percentile of topography ranges for the spatial estimations to remove some very low valleys or peaks [default:90.]
--crop VALUE            Define a region of interest for the temporal decomposition [default: 0,nlign,0,ncol]
--cond VALUE            Condition value for optimization: Singular value smaller than cond are considered zero [default: 1e-3]
--ineq VALUE            If yes, add ineguality constraints in the inversion: use least square result without post-seismic functions as a first guess to iterate the inversion. Force postseismic to be the same sign and inferior than coseismic steps of the first guess [default: no].
--fulloutput YES/NO     If yes produce maps of models, residuals, ramps, as well as flatten cube without seasonal and linear term [default: no]
--geotiff PATH          Path to Geotiff to save outputs in tif format. If None save output are saved as .r4 files [default: .r4]
--plot YES/NO           Display plots [default: no]
--ibeg VALUE            Line numbers bounding the ramp estimation zone [default: 0]
--iend VALUE            Line numbers bounding the ramp estimation zone [default: nlign]
--jbeg VALUE            Column numbers bounding the ramp estimation zone [default: 0]
--jend VALUE            Column numbers bounding the ramp estimation zone [default: ncol]
"""

print
print '# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #'
print '#'                                                                 '#'
print '#         Linear Inversion of InSAR time series displacements       #'
print '#         with a decomposition in time                              #'
print '#'                                                                 '#'
print '# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #'
print

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided
# import numpy.linalg as lst
# scipy
import scipy as sp
import scipy.optimize as opt
import scipy.linalg as lst
# from nsbas import gdal, osr
import gdal, osr
# basic
import math,sys,getopt
from os import path, environ
import os
# plot
import matplotlib
#matplotlib.use('TkAgg') # Must be before importing matplotlib.pyplot or pylab!
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.dates as mdates
from datetime import datetime
import datetime

# docopt (command line parser)
import docopt

np.warnings.filterwarnings('ignore')

################################
# Create lib of wavelet functions
################################

class pattern:
    def __init__(self,name,reduction,date):
        self.name=name
        self.reduction=reduction
        self.date=date

    def info(self):
        print self.name, self.date

### BASIS FUNCTIONS: function of time

def Heaviside(t):
        h=np.zeros((len(t)))
        h[t>=0]=1.0
        return h

def Box(t):
        return Heaviside(t+0.5)-Heaviside(t-0.5)

class coseismic(pattern):
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

class interseismic(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
        self.to=date

    def g(self,t):
        func=(t-self.to)
        return func

class sin2var(pattern):
     def __init__(self,name,reduction,date):
         pattern.__init__(self,name,reduction,date)
         self.to=date

     def g(self,t):
         func=np.zeros(t.size)
         for i in xrange(t.size):
             func[i]=math.sin(4*math.pi*(t[i]-self.to))
         return func

class cos2var(pattern):
     def __init__(self,name,reduction,date):
         pattern.__init__(self,name,reduction,date)
         self.to=date

     def g(self,t):
         func=np.zeros(t.size)
         for i in xrange(t.size):
             func[i]=math.cos(4*math.pi*(t[i]-self.to))
         return func

class sinvar(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
        self.to=date

    def g(self,t):
        func=np.zeros(t.size)
        for i in xrange(t.size):
            func[i]=math.sin(2*math.pi*(t[i]-self.to))
        return func

class cosvar(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
        self.to=date

    def g(self,t):
        func=np.zeros(t.size)
        for i in xrange(t.size):
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
class corrdem(pattern):
    def __init__(self,name,reduction,bp0,bp):
        self.name = name
        self.reduction = reduction
        self.bpo=bp0
        self.bp=bp

    def info(self):
        print self.name

    def g(self,index):
        func = (self.bp-self.bpo)
        return func[index]

class vector(pattern):
    def __init__(self,name,reduction,vect):
        self.name = name
        self.reduction = reduction
        self.func=vect

    def info(self):
        print self.name

    def g(self,index):
        return self.func[index]

################################
# Initialization
################################

# read arguments
arguments = docopt.docopt(__doc__)
if arguments["--lectfile"] ==  None:
    infile = "lect.in"
else:
    infile = arguments["--lectfile"]
if arguments["--list_images"] ==  None:
    listim = "images_retenues"
else:
    listim = arguments["--list_images"]
if arguments["--aps"] ==  None:
    apsf = 'no'
else:
    apsf = arguments["--aps"]
if arguments["--interseismic"] ==  None:
    inter = 'no'
else:
    inter = arguments["--interseismic"]
if arguments["--seasonal"] ==  None:
    seasonal = 'no'
else:
    seasonal = arguments["--seasonal"]
if arguments["--semianual"] ==  None:
    semianual = 'no'
else:
    semianual = arguments["--semianual"]
if arguments["--dem"] ==  None:
    dem = 'no'
else:
    dem = arguments["--dem"]
if arguments["--coseismic"] ==  None:
    cos = []
else:
    cos = map(float,arguments["--coseismic"].replace(',',' ').split())

if arguments["--postseismic"] ==  None:
    pos = []
else:
    pos = map(float,arguments["--postseismic"].replace('None','-1').replace(',',' ').split())

if len(pos)>0 and len(cos) != len(pos):
    raise Exception("coseimic and postseismic lists are not the same size")

if arguments["--slowslip"] == None:
    sse=[]
else:
    sse = map(float,arguments["--slowslip"].replace(',',' ').split())
sse_time = sse[::2]
sse_car = sse[1::2]

if arguments["--vector"] != None:
    vectf = arguments["--vector"].replace(',',' ').split()
else:
    vect = None

if arguments["--refstart"] == None:
    refstart = 0
else:
    refstart = int(arguments["--refstart"])

if arguments["--refend"] == None:
    refend = 200
else:
    refend = int(arguments["--refend"])

# read lect.in
ncol, nlign = map(int, open(infile).readline().split(None, 2)[0:2])

if arguments["--niter"] ==  None:
    niter = 1
else:
    niter = int(arguments["--niter"])
if arguments["--spatialiter"] ==  None:
    spatialiter = 'no'
else:
    spatialiter = arguments["--spatialiter"]
if arguments["--flat"] == None:
    flat = 0
elif int(arguments["--flat"]) <  10:
    flat = int(arguments["--flat"])
else:
    flat = 0
if arguments["--sampling"] ==  None:
    sampling = 1
else:
    sampling = int(arguments["--sampling"])
if arguments["--mask"] ==  None:
    maskfile = None
else:
    maskfile = arguments["--mask"]
if arguments["--rampmask"] ==  None:
    rampmask = 'no'
else:
    rampmask = arguments["--rampmask"]
if arguments["--threshold_mask"] ==  None:
    seuil = 0.
else:
    seuil = float(arguments["--threshold_mask"])
if arguments["--threshold_rms"] ==  None:
    seuil_rms = 1.
else:
    seuil_rms = float(arguments["--threshold_rms"])
if arguments["--tempmask"] ==  None:
    tempmask = 'no'
else:
    tempmask = arguments["--tempmask"]
if arguments["--scale_mask"] ==  None:
    scale = 1
else:
    scale = float(arguments["--scale_mask"])
if arguments["--topofile"] ==  None:
   radar = None
else:
   radar = arguments["--topofile"]

if arguments["--aspect"] ==  None:
   aspect = None
else:
   aspect = arguments["--aspect"]

if arguments["--imref"] ==  None:
    imref = 0
elif arguments["--imref"] < 1:
    print '--imref must be between 1 and Nimages'
else:
    imref = int(arguments["--imref"]) - 1

if arguments["--crop"] ==  None:
    crop = [0,nlign,0,ncol]
else:
    crop = map(float,arguments["--crop"].replace(',',' ').split())
ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])

if arguments["<ibeg>"] ==  None:
    ibegref = 0
else:
    ibegref = int(arguments["<ibeg>"])
if arguments["<iend>"] ==  None:
    iendref = nlign
else:
    iendref = int(arguments["<iend>"])
if arguments["<jbeg>"] ==  None:
    jbegref = 0
else:
    jbegref = int(arguments["<jbeg>"])
if arguments["<jend>"] ==  None:
    jendref = ncol
else:
    jendref = int(arguments["<jend>"])

if arguments["--cond"] ==  None:
    rcond = 1e-3
else:
    rcond = float(arguments["--cond"])
if arguments["--rmspixel"] ==  None:
    rmsf = None
else:
    rmsf = arguments["--rmspixel"]
if arguments["--ineq"] ==  None:
    ineq = 'no'
else:
    ineq = arguments["--ineq"]
if arguments["--threshold_rmsd"] ==  None:
    maxrmsd = 0.
else:
    maxrmsd = float(arguments["--threshold_rmsd"])
if arguments["--fulloutput"] ==  None:
    fulloutput = 'no'
else:
    fulloutput = arguments["--fulloutput"]

if arguments["--plot"] ==  None:
    plot = 'no'
else:
    plot = arguments["--plot"]

if arguments["--cube"] ==  None:
    cubef = "depl_cumule"
else:
    cubef = arguments["--cube"]

if arguments["--geotiff"] ==  None:
    geotiff = None
else:
    geotiff = arguments["--geotiff"]
    georef = gdal.Open(geotiff)
    gt = georef.GetGeoTransform()
    proj = georef.GetProjection()
    driver = gdal.GetDriverByName('GTiff')

if arguments["--ivar"] == None:
    ivar = 0
elif int(arguments["--ivar"]) <  2:
    ivar = int(arguments["--ivar"])
else:
    print 'Error: ivar > 1, set ivar to 0'
    ivar = 0

if arguments["--nfit"] == None:
    nfit = 0
elif int(arguments["--nfit"]) <  2:
    nfit = int(arguments["--nfit"])
else:
    print 'Error: nfit > 1, set nfit to 0'
    nfit = 0

if arguments["--perc_topo"] ==  None:
    perc_topo = 90.
else:
    perc_topo = float(arguments["--perc_topo"])

if arguments["--perc_los"] ==  None:
    perc_los = 98.
else:
    perc_los = float(arguments["--perc_los"])


if len(cos) > 0:
    print
    print 'Define a maximal RMSD for adding coseismic and postseismic basis functions in the inversion'
    print 'Max RMSD:',  maxrmsd
    print

#######################################################

# cm
cmap = cm.jet
cmap.set_bad('white')

# load images_retenues file
nb,idates,dates,base=np.loadtxt(listim, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')

N=len(dates)
print 'Number images: ', N
datemin, datemax = np.int(np.nanmin(dates)), np.int(np.nanmax(dates))+1

# lect cube
cubei = np.fromfile(cubef,dtype=np.float32)

# extract
cube = as_strided(cubei[:nlign*ncol*N])
print 'Number of line in the cube: ', cube.shape

# !!! remove crazy values !!!
kk = np.flatnonzero(cube>9990)
cube[kk] = float('NaN')
#cube[kk] = 0

maps = cube.reshape((nlign,ncol,N))
print 'Reshape cube: ', maps.shape
# set at NaN zero values for all dates
# kk = np.nonzero(maps[:,:,-1]==0)
# ref displacements to ref date
cst = np.copy(maps[:,:,imref])
for l in xrange((N)):
    maps[:,:,l] = maps[:,:,l] - cst
    if l != imref:
        index = np.nonzero(maps[:,:,l]==0.0)
        maps[:,:,l][index] = np.float('NaN')

# fig = plt.figure(0)
# plt.imshow(cst,vmax=1,vmin=-1)
# fig = plt.figure(1)
# plt.imshow(maps[ibeg:iend,jbeg:jend,imref],vmax=1,vmin=-1)
# plt.show()
# sys.exit()

fig = plt.figure(12)
nfigure=0

# open mask file
mask = np.zeros((nlign,ncol))
if maskfile is not None:
    extension = os.path.splitext(maskfile)[1]
    if extension == ".tif":
      ds = gdal.Open(maskfile, gdal.GA_ReadOnly)
      band = ds.GetRasterBand(1)
      maski = band.ReadAsArray().flatten()*scale
      del ds
    else:
      fid = open(maskfile,'r')
      maski = np.fromfile(fid,dtype=np.float32)*scale
      fid.close()
    maski = maski[:nlign*ncol]
    mask = maski.reshape((nlign,ncol))
else:
    mask_flat = np.ones((nlign,ncol))

elev = np.zeros((nlign,ncol))
# open elevation map
if radar is not None:
    extension = os.path.splitext(radar)[1]
    if extension == ".tif":
      ds = gdal.Open(radar, gdal.GA_ReadOnly)
      band = ds.GetRasterBand(1)
      elevi = band.ReadAsArray().flatten()
      del ds
    else:
      fid = open(radar,'r')
      elevi = np.fromfile(fid,dtype=np.float32)
      fid.close()

    elevi = elevi[:nlign*ncol]
    # fig = plt.figure(10)
    # plt.imshow(elevi.reshape(nlign,ncol)[ibeg:iend,jbeg:jend])
    elev = elevi.reshape((nlign,ncol))
    elev[np.isnan(maps[:,:,-1])] = float('NaN')
    kk = np.nonzero(abs(elev)>9999.)
    elev[kk] = float('NaN')
    # fig = plt.figure(11)
    # plt.imshow(elev[ibeg:iend,jbeg:jend])
    # plt.show()
    # sys.exit()
else:
   elev = np.ones((nlign,ncol))

if aspect is not None:
    extension = os.path.splitext(aspect)[1]
    if extension == ".tif":
      ds = gdal.Open(aspect, gdal.GA_ReadOnly)
      band = ds.GetRasterBand(1)
      aspecti = band.ReadAsArray().flatten()
      # print ds.RasterYSize, ds.RasterXSize
      del ds
    else:
      fid = open(aspect,'r')
      aspecti = np.fromfile(fid,dtype=np.float32)
      fid.close()
    aspecti = aspecti[:nlign*ncol]
    slope = aspecti.reshape((nlign,ncol))
    slope[np.isnan(maps[:,:,-1])] = float('NaN')
    kk = np.nonzero(abs(slope>9999.))
    slope[kk] = float('NaN')
    # print slope[slope<0]
    # fig = plt.figure(11)
    # plt.imshow(elev[ibeg:iend,jbeg:jend])
    # plt.show()
    # sys.exit()
else:
    slope = np.ones((nlign,ncol))

if rmsf is not None:
    rmsmap = np.fromfile(rmsf,dtype=np.float32).reshape((nlign,ncol))
    rmsmap = rmsmap[:nlign,:ncol]
    kk = np.nonzero(np.logical_or(rmsmap==0.0, rmsmap>999.))
    rmsmap[kk] = float('NaN')
    kk = np.nonzero(rmsmap>seuil_rms)
    spacial_mask = np.copy(rmsmap)
    spacial_mask[kk] = float('NaN')
    fig = plt.figure(nfigure,figsize=(9,4))
    nfigure = nfigure + 1
    ax = fig.add_subplot(1,1,1)
    cax = ax.imshow(spacial_mask,cmap=cmap)
    ax.set_title('Mask on spatial estimation based on RMSpixel')
    setp( ax.get_xticklabels(), visible=False)
    fig.colorbar(cax, orientation='vertical',aspect=10)
    del spacial_mask
    # if plot=='yes':
    #    plt.show()
else:
    rmsmap = np.ones((nlign,ncol))

#fig, ax = plt.subplots(1)
#plt.imshow(rmsmap[:,:], vmax=1, vmin=0)
#plt.show()
#sys.exit()

# plot last date
#fig, ax = plt.subplots(1)
#plt.imshow(maps[:,:,-1], vmax=1, vmin=-1)
#plt.title('Last image')
#plt.show()
#sys.exit()

# plot bperp vs time
fig = plt.figure(nfigure,figsize=(10,4))
nfigure = nfigure + 1
ax = fig.add_subplot(1,2,1)
# convert idates to num
x = [date2num(datetime.datetime.strptime('{}'.format(d),'%Y%m%d')) for d in idates]
# format the ticks
ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
ax.plot(x,base,"ro",label='Baseline history of the {} images'.format(N))
ax.plot(x,base,"green")
# rotates and right aligns the x labels, and moves the bottom of the
# axes up to make room for them
fig.autofmt_xdate()
ax.set_xlabel('Time (Year/month/day)')
ax.set_ylabel('Perpendicular Baseline')
plt.legend(loc='best')

ax = fig.add_subplot(1,2,2)
ax.plot(np.mod(dates,1),base,"ro",label='Baseline seasonality of the {} images'.format(N))
plt.legend(loc='best')

fig.savefig('baseline.eps', format='EPS',dpi=150)
np.savetxt('bp_t.in', np.vstack([dates,base]).T, fmt='%.6f')

if vect is not None:
    v = np.loadtxt(vect, comments='#', unpack = False, dtype='f')
    fig = plt.figure(nfigure,figsize=(6,4))
    nfigure = nfigure + 1
    ax = fig.add_subplot(1,1,1)
    ax.plot(v,label='Vector')
    plt.legend(loc='best')
    if plot=='yes':
        plt.show()
    # sys.exit()


if maskfile is not None:
    print
    print 'Flatten mask...'
    print
    los_temp = as_strided(mask[ibegref:iendref,jbegref:jendref]).flatten()

    if rampmask=='yes':
        temp = [(i,j) for i in xrange(iendref-ibegref) for j in xrange(jendref-jbegref) \
        if np.logical_and((math.isnan(los_temp[i*(jendref-jbegref)+j]) is False), \
            (los_temp[i*(jendref-jbegref)+j]>seuil))]

        temp2 = np.array(temp)
        x = temp2[:,0]; y = temp2[:,1]
        los_clean = los_temp[x*(jend-jbeg)+y]
        G=np.zeros((len(los_clean),4))
        G[:,0], G[:,1], G[:,2], G[:,3] = y**2, y, x, 1
        # ramp inversion
        pars = np.dot(np.dot(np.linalg.inv(np.dot(G.T,G)),G.T),los_clean)
        a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
        print 'Remove ramp mask %f x**2 %f x  + %f y + %f for : %s'%(a,b,c,d,maskfile)

        # remove 0 values
        kk = np.flatnonzero(np.logical_or(maski==0, maski==9999))
        #kk = np.flatnonzero(los==9999)
        maski[kk] = float('NaN')

        G=np.zeros((len(maski),4))
        for i in xrange(nlign):
            G[i*ncol:(i+1)*ncol,0] = (np.arange((ncol)) - jbegref)**2
            G[i*ncol:(i+1)*ncol,1] = np.arange((ncol)) - jbegref
            G[i*ncol:(i+1)*ncol,2] = i - ibegref
        G[:,3] = 1
        mask_flat = (maski - np.dot(G,pars)).reshape(nlign,ncol)
        mask_flat = mask_flat - np.nanmean(mask_flat)

    else:
        # if rampmask is no yes then  mask_flat is mask and los_clean is selected area mask without NaN
        temp = [(i,j) for i in xrange(iendref-ibegref) for j in xrange(jendref-jbegref) \
        if math.isnan(los_temp[i*(jendref-jbegref)+j]) is False]

        temp2 = np.array(temp)
        x = temp2[:,0]; y = temp2[:,1]
        los_clean = los_temp[x*(jendref-jbegref)+y]

        # remove 0 values
        kk = np.flatnonzero(np.logical_or(maski==0, maski==9999))
        #kk = np.flatnonzero(los==9999)
        maski[kk] = float('NaN')
        mask_flat = maski.reshape(nlign,ncol)

    del maski

    # check seuil
    kk = np.flatnonzero(mask_flat<seuil)
    mask_flat_clean=np.copy(mask_flat.flatten())
    mask_flat_clean[kk]=float('NaN')
    mask_flat_clean = mask_flat_clean.reshape(nlign,ncol)

    # mask maps if necessary for temporal inversion
    if tempmask=='yes':
        kk = np.nonzero(mask_flat[ibegref:iendref,jbegref:jendref]<seuil)
        for l in xrange((N)):
            # clean only selected area
            d = as_strided(maps[ibegref:iendref,jbegref:jendref,l])
            d[kk] = np.float('NaN')

    # plots
    nfigure+=1
    fig = plt.figure(nfigure,figsize=(7,6))
    vmax = np.abs([np.nanmedian(mask_flat) + nanstd(mask_flat),\
        np.nanmedian(mask_flat) - nanstd(mask_flat)]).max()
    vmin = -vmax

    ax = fig.add_subplot(1,3,1)
    cax = ax.imshow(mask,cmap=cmap,vmax=vmax,vmin=vmin)
    ax.set_title('Original Mask')
    setp( ax.get_xticklabels(), visible=False)

    ax = fig.add_subplot(1,3,2)
    cax = ax.imshow(mask_flat,cmap=cmap,vmax=vmax,vmin=vmin)
    ax.set_title('Flat Mask')
    setp( ax.get_xticklabels(), visible=False)
    #cbar = fig.colorbar(cax, orientation='vertical',aspect=10)

    ax = fig.add_subplot(1,3,3)
    cax = ax.imshow(mask_flat_clean,cmap=cmap,vmax=vmax,vmin=vmin)
    ax.set_title('Final Mask')
    setp( ax.get_xticklabels(), visible=False)
    #cbar = fig.colorbar(cax, orientation='vertical',aspect=10)
    fig.savefig('mask.eps', format='EPS',dpi=150)
    del mask_flat_clean

    if plot=='yes':
        plt.show()
    # sys.exit()


# plot diplacements maps
nfigure+=1
fig = plt.figure(nfigure,figsize=(14,10))
fig.subplots_adjust(wspace=0.001)
vmax = np.nanpercentile(maps[:,:,-1],98.)
vmin = np.nanpercentile(maps[:,:,-1],2.)
# vmax = np.abs([np.nanmedian(maps[:,:,-1]) + 1.*np.nanstd(maps[:,:,-1]),\
#     np.nanmedian(maps[:,:,-1]) - 1.*np.nanstd(maps[:,:,-1])]).max()
# vmin = -vmax

for l in xrange((N)):
    d = as_strided(maps[ibeg:iend,jbeg:jend,l])
    #ax = fig.add_subplot(1,N,l+1)
    ax = fig.add_subplot(4,int(N/4)+1,l+1)
    #cax = ax.imshow(d,cmap=cmap,vmax=vmax,vmin=vmin)
    cax = ax.imshow(d,cmap=cmap,vmax=vmax,vmin=vmin)
    ax.set_title(idates[l],fontsize=6)
    setp( ax.get_xticklabels(), visible=False)
    setp( ax.get_yticklabels(), visible=False)

plt.suptitle('Time series maps')
fig.colorbar(cax, orientation='vertical',aspect=10)
fig.tight_layout()
fig.savefig('maps.eps', format='EPS',dpi=150)

# plt.show()
#sys.exit()

#######################################################
# Save new lect.in file
#######################################################

fid = open('lect_ts.in','w')
np.savetxt(fid, (jend-jbeg,iend-ibeg),fmt='%6i',newline='\t')
fid.close()

#######################################################
# Create functions of decomposition
######################################################

basis=[
    reference(name='reference',date=datemin,reduction='ref'),
    ]
index = len(basis)

# initialise iteration with interseismic alone
iteration=False

if inter=='yes':
    indexinter=index
    basis.append(interseismic(name='interseismic',reduction='lin',date=datemin))
    index = index + 1

if seasonal=='yes':
    indexseas = index
    basis.append(cosvar(name='seas. var (cos)',reduction='coswt',date=datemin))
    basis.append(sinvar(name='seas. var (sin)',reduction='sinwt',date=datemin))
    index = index + 2

if semianual=='yes':
     indexsemi = index
     basis.append(cos2var(name='semi-anual var (cos)',reduction='cosw2t',date=datemin))
     basis.append(sin2var(name='semi-anual var (sin)',reduction='sinw2t',date=datemin))
     index = index + 2

indexco = np.zeros(len(cos))
for i in xrange(len(cos)):
    basis.append(coseismic(name='coseismic {}'.format(i),reduction='cos{}'.format(i),date=cos[i])),
    indexco[i] = index
    index = index + 1
    iteration=True


indexpo,indexpofull = [],[]
for i in xrange(len(pos)):
  if pos[i] > 0. :
    basis.append(postseismic(name='postseismic {}'.format(i),reduction='post{}'.format(i),date=cos[i],tcar=pos[i])),
    indexpo.append(int(index))
    indexpofull.append(int(index))
    index = index + 1
  else:
    indexpofull.append(0)
indexpo = np.array(indexpo)
indexpofull = np.array(indexpofull)


indexsse = np.zeros(len(sse_time))
for i in xrange(len(sse_time)):
    basis.append(slowslip(name='sse {}'.format(i),reduction='sse{}'.format(i),date=sse_time[i],tcar=sse_car[i])),
    indexsse[i] = int(index)
    index = index + 1
    iteration=True

kernels=[]

if dem=='yes':
   kernels.append(corrdem(name='dem correction',reduction='corrdem',bp0=base[imref],bp=base))
   indexdem = index
   index = index + 1

if vect != None:
   indexvect = np.zeros(len(vectf))
   for i in xrange(len(vectf)):
     kernels.append(vector(name=vectf[i],reduction='vector_{}'.format(i),vect=v[i]))
     indexvect[i] = index
     index = index + 1

indexpo = indexpo.astype(int)
indexco = indexco.astype(int)
indexsse = indexsse.astype(int)

Mbasis=len(basis)
print 'Number of basis functions:', Mbasis
Mker=len(kernels)
print 'Number of kernel functions:', Mker
M = Mbasis + Mker

print
print 'Basis functions, Time:'
for i in xrange((Mbasis)):
    basis[i].info()
print 'Kernels functions, Time:'
for i in xrange((Mker)):
    kernels[i].info()

# initialize matrix model to NaN
for l in xrange((Mbasis)):
    basis[l].m = np.ones((iend-ibeg,jend-jbeg))*np.float('NaN')
    basis[l].sigmam = np.ones((iend-ibeg,jend-jbeg))*np.float('NaN')
for l in xrange((Mker)):
    kernels[l].m = np.ones((iend-ibeg,jend-jbeg))*np.float('NaN')
    kernels[l].sigmam = np.ones((iend-ibeg,jend-jbeg))*np.float('NaN')

# initialize qual
if apsf=='no':
    inaps=np.ones((N)) # no weigthing for the first itertion
else:
    fimages=apsf
    inaps=np.loadtxt(fimages, comments='#', dtype='f')
    print
    print 'Input uncertainties:', inaps
    print 'Scale input uncertainties between 0 and 1 and set very low values to the 2 percentile to avoid overweighting...'
    # maxinaps = np.nanmax(inaps)
    # inaps= inaps/maxinaps
    minaps= np.nanpercentile(inaps,2)
    index = flatnonzero(inaps<minaps)
    inaps[index] = minaps
    print 'Output uncertainties for first iteration:', inaps
    print

# SVD inversion with cut-off eigenvalues
def invSVD(A,b,cond):
    try:
        U,eignv,V = lst.svd(A, full_matrices=False)
        s = np.diag(eignv)
        index = np.nonzero(s<cond)
        inv = lst.inv(s)
        inv[index] = 0.
        fsoln = np.dot( V.T, np.dot( inv , np.dot(U.T, b) ))
    except:
        fsoln = lst.lstsq(A,b)[0]
        #fsoln = lst.lstsq(A,b,rcond=cond)[0]
    
    return fsoln

## inversion procedure 
def consInvert(A,b,sigmad,ineq='no',cond=1.0e-3, iter=2000,acc=1e-12):
    '''Solves the constrained inversion problem.

    Minimize:
    
    ||Ax-b||^2

    Subject to:
    mmin < m < mmax
    '''

    if A.shape[0] != len(b):
        raise ValueError('Incompatible dimensions for A and b')

    if ineq == 'no':
        
        fsoln = invSVD(A,b,cond)
        
    else:

        # prior solution without postseismic 
        Ain = np.delete(A,indexpo,1)
        mtemp = invSVD(Ain,b,cond)
        #print mtemp
        
        # rebuild full vector
        for z in xrange(len(indexpo)):
            mtemp = np.insert(mtemp,indexpo[z],0)
        minit = np.copy(mtemp)

        # # initialize bounds
        mmin,mmax = -np.ones(M)*np.inf, np.ones(M)*np.inf 

        # We here define bounds for postseismic to be the same sign than coseismic
        # and coseismic inferior or egual to the coseimic initial 
        for i in xrange(len(indexco)):
            if (pos[i] > 0.) and (minit[int(indexco[i])]>0.):
                mmin[int(indexpofull[i])], mmax[int(indexpofull[i])] = 0, np.inf 
                mmin[int(indexco[i])], mmax[int(indexco[i])] = 0, minit[int(indexco[i])] 
            if (pos[i] > 0.) and (minit[int(indexco[i])]<0.):
                mmin[int(indexpofull[i])], mmax[int(indexpofull[i])] = -np.inf , 0
                mmin[int(indexco[i])], mmax[int(indexco[i])] = minit[int(indexco[i])], 0
        
        # print mmin,mmax
        ####Objective function and derivative
        _func = lambda x: np.sum(((np.dot(A,x)-b)/sigmad)**2)
        _fprime = lambda x: 2*np.dot(A.T/sigmad, (np.dot(A,x)-b)/sigmad)
        
        bounds=zip(mmin,mmax)
        res = opt.fmin_slsqp(_func,minit,bounds=bounds,fprime=_fprime, \
            iter=iter,full_output=True,iprint=0,acc=acc)  
        fsoln = res[0]
  
    # tarantola:
    # Cm = (Gt.Cov.G)-1 --> si sigma=1 problems
    # sigma m **2 =  misfit**2 * diag([G.TG]-1)
    try:
       varx = np.linalg.inv(np.dot(A.T,A))
       res2 = np.sum(pow((b-np.dot(A,fsoln)),2))
       scale = 1./(A.shape[0]-A.shape[1])
       # scale = 1./A.shape[0]
       sigmam = np.sqrt(scale*res2*np.diag(varx))
    except:
       sigmam = np.ones((A.shape[1]))*float('NaN')

    return fsoln,sigmam

# initialization
maps_flata = np.copy(maps)
models = np.zeros((nlign,ncol,N))

# prepare flatten maps
maps_ramp = np.zeros((nlign,ncol,N))
maps_topo = np.zeros((nlign,ncol,N))
maps_noramps = np.zeros((nlign,ncol,N))
rms = np.zeros((N))

for ii in xrange(niter):
    print
    print '---------------'
    print 'iteration: ', ii
    print '---------------'

    #############################
    # SPATIAL ITERATION N  ######
    #############################

    print
    print 'Spatial correction..'
    print

    def estim_ramp(los,los_clean,topo_clean,x,y,order,rms,nfit,ivar,los_ref):

      # initialize topo
      topo = np.zeros((nlign,ncol))
      ramp = np.zeros((nlign,ncol))

      if radar is not None:
        # calc elevi as los
        elev_temp = np.matrix.copy(elevi)
      
      # create new data matrix with cst 
      data = np.hstack([los_clean,los_ref])
      rms = np.hstack([rms,1e-3])
    
      if order==0:

        if radar is None:
            
            a = cst
            print 'Remove ref frame %f  for date: %i'%(a,idates[l])
            ramp = np.ones((nlign,ncol))*a
            rms = np.sqrt(np.nanmean((los-a)**2))
            print 'RMS:', rms

        else:

            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),2))
                G[:,0] = 1
                G[:-1,1] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]
                print 'Remove ref frame %f + %f z for date: %i'%(a,b,idates[l])

                # plot phase/elev
                funct = a
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,b*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),2))
                G[:,0] = 1
                G[:,1] = elev_temp

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                topo = np.dot(G,pars).reshape(nlign,ncol)


            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),3))
                G[:,0] = 1
                G[:-1,1] = topo_clean
                G[:-1,2] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c=pars[2]
                print 'Remove ref frame %f + %f z + %f z**2 for date: %i'%(a,b,c,idates[l])

                # plot phase/elev
                funct = a
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,b*x+c*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),3))
                G[:,0] = 1
                G[:,1] = elev_temp
                G[:,2] = elev_temp**2

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                topo = np.dot(G,pars).reshape(nlign,ncol)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),3))
                G[:,0] = 1
                G[:-1,1] = topo_clean
                G[:-1,2] = x*topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]
                print 'Remove ref frame %f + %f z + %f az*z for date: %i'%(a,b,c,idates[l])

                # plot phase/elev
                funct = a + c*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,b*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),3))
                G[:,0] = 1
                G[:,1] = elev_temp
                G[:,2] = elev_temp
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,2] *= (i - ibegref)

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                topo = np.dot(G,pars).reshape(nlign,ncol)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),4))
                G[:,0] = 1
                G[:-1,1] = x*topo_clean
                G[:-1,2] = topo_clean
                G[:-1,3] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
                print 'Remove ref frame %f + %f az*z + %f z + %f z**2 for date: %i'%(a,b,c,d,idates[l])

                # plot phase/elev
                funct = a + b*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,c*x+d*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                G[:,0] = 1
                G[:,1] = elev_temp
                G[:,2] = elev_temp
                G[:,3] = elev_temp**2
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,1] *= (i - ibegref)

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                topo = np.dot(G,pars).reshape(nlign,ncol)

      elif order==1: # Remove a range ramp ay+b for each maps (y = col)

        if radar is None:
            G=np.zeros((len(data),2))
            G[:-1,0] = y
            G[:,1] = 1

            # ramp inversion
            x0 = lst.lstsq(G,data)[0]
            _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
            _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
            pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
            a = pars[0]; b = pars[1]
            print 'Remove ramp %f r + %f for date: %i'%(a,b,idates[l])

            # build total G matrix
            G=np.zeros((len(los),2))
            for i in xrange(nlign):
                G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
            G[:,1] = 1


            res = los - np.dot(G,pars)
            rms = np.sqrt(np.nanmean(res**2))
            print 'RMS:', rms

            ramp = np.dot(G,pars).reshape(nlign,ncol)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),3))
                G[:-1,0] = y
                G[:,1] = 1
                G[:-1,2] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]
                print 'Remove ramp %f r + %f + %f z for date: %i'%(a,b,c,idates[l])

                # plot phase/elev
                funct = a*y + b
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,c*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),3))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                G[:,1] = 1
                G[:,2] = elev_temp

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(nlign,ncol)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),4))
                G[:-1,0] = y
                G[:,1] = 1
                G[:-1,2] = topo_clean
                G[:-1,3] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d=pars[3]
                print 'Remove ramp %f r + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,idates[l])

                # plot phase/elev
                funct = a*y + b
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,c*x+d*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                G[:,1] = 1
                G[:,2] = elev_temp
                G[:,3] = elev_temp**2

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),4))
                G[:-1,0] = y
                G[:-1,1] = 1
                G[:-1,2] = topo_clean
                G[:-1,3] = topo_clean*x

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
                print 'Remove ramp %f r + %f + %f z + %f z*az for date: %i'%(a,b,c,d,idates[l])

                # plot phase/elev
                funct = a*y + b + d*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,c*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                G[:,1] = 1
                G[:,2] = elev_temp
                G[:,3] = elev_temp
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,3] *= (i - ibegref)


                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),5))
                G[:-1,0] = y
                G[:,1] = 1
                G[:-1,2] = topo_clean*x
                G[:-1,3] = topo_clean
                G[:-1,4] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
                print 'Remove ramp %f r + %f +  %f z*az + %f z + %f z**2 for date: %i'%(a,b,c,d,e,idates[l])

                # plot phase/elev
                funct = a*y + b + c*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,d*x+e*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                G[:,1] = 1
                G[:,2] = elev_temp
                G[:,3] = elev_temp
                G[:,4] = elev_temp**2
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,2] *= (i - ibegref)

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(nlign,ncol)


      elif order==2: # Remove a azimutal ramp ax+b for each maps (x is lign)
        if radar is None:
            G=np.zeros((len(data),2))
            G[:-1,0] = x
            G[:,1] = 1

            # ramp inversion
            x0 = lst.lstsq(G,data)[0]
            _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
            _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
            pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
            a = pars[0]; b = pars[1]
            print 'Remove ramp %f az + %f for date: %i'%(a,b,idates[l])

            # build total G matrix
            G=np.zeros((len(los),2))
            for i in xrange(nlign):
                G[i*ncol:(i+1)*ncol,0] =(i - ibegref)
            G[:,1] = 1


            res = los - np.dot(G,pars)
            rms = np.sqrt(np.nanmean(res**2))
            print 'RMS:', rms

            ramp = np.dot(G,pars).reshape(nlign,ncol)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),3))
                G[:-1,0] = x
                G[:,1] = 1
                G[:-1,2] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]
                print 'Remove ramp %f az + %f + %f z for date: %i'%(a,b,c,idates[l])

                # plot phase/elev
                funct = a*x + b
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,c*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),3))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] =(i - ibegref)
                G[:,1] = 1
                G[:,2] = elev_temp


                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(nlign,ncol)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),4))
                G[:-1,0] = x
                G[:,1] = 1
                G[:-1,2] = topo_clean
                G[:-1,3] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
                print 'Remove ramp %f az + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,idates[l])

                # plot phase/elev
                funct = a*x + b
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,c*x + d*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] =(i - ibegref)
                G[:,1] = 1
                G[:,2] = elev_temp
                G[:,3] = elev_temp**2

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),4))
                G[:-1,0] = x
                G[:,1] = 1
                G[:-1,2] = topo_clean
                G[:-1,3] = topo_clean*x

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
                print 'Remove ramp %f az + %f + %f z + %f z*az for date: %i'%(a,b,c,d,idates[l])

                # plot phase/elev
                funct = a*x + b + d*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,c*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                G[:,1] = 1
                G[:,2] = elev_temp
                G[:,3] = elev_temp
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = i - ibegref
                    G[i*ncol:(i+1)*ncol,3] *= (i - ibegref)


                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),5))
                G[:-1,0] = x
                G[:,1] = 1
                G[:-1,2] = topo_clean*x
                G[:-1,3] = topo_clean
                G[:-1,4] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
                print 'Remove ramp %f az + %f + %f z*az + %f z + %f z**2 for date: %i'%(a,b,c,d,e,idates[l])

                # plot phase/elev
                funct = a*x + b + c*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,d*x+e*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                G[:,1] = 1
                G[:,2] = elev_temp
                G[:,3] = elev_temp
                G[:,4] = elev_temp**2
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = i - ibegref
                    G[i*ncol:(i+1)*ncol,2] *= (i - ibegref)

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(nlign,ncol)

      elif order==3: # Remove a ramp ay+bx+c for each maps
        if radar is None:
            G=np.zeros((len(data),3))
            G[:-1,0] = y
            G[:-1,1] = x
            G[:,2] = 1

            # ramp inversion
            x0 = lst.lstsq(G,data)[0]
            # print x0
            # x0 = np.zeros((3))
            _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
            _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
            pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
            a = pars[0]; b = pars[1]; c = pars[2]
            print 'Remove ramp %f r  + %f az + %f for date: %i'%(a,b,c,idates[l])

            # build total G matrix
            G=np.zeros((len(los),3))
            for i in xrange(nlign):
                G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                G[i*ncol:(i+1)*ncol,1] = (i - ibegref)
            G[:,2] = 1

            res = los - np.dot(G,pars)
            rms = np.sqrt(np.nanmean(res**2))
            print 'RMS:', rms

            ramp = np.dot(G,pars).reshape(nlign,ncol)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),4))
                G[:-1,0] = y
                G[:-1,1] = x
                G[:,2] = 1
                G[:-1,3] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
                print 'Remove ramp %f r  + %f az + %f + %f z for date: %i'%(a,b,c,d,idates[l])

                # plot phase/elev
                funct = a*y + b*x + c
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,d*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,1] =(i - ibegref)
                G[:,2] = 1
                G[:,3] = elev_temp


                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(nlign,ncol)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),5))
                G[:-1,0] = y
                G[:-1,1] = x
                G[:,2] = 1
                G[:-1,3] = topo_clean
                G[:-1,4] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]

                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
                print 'Remove ramp %f r  + %f az + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,e,idates[l])

                # plot phase/elev
                funct = a*y + b*x + c
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,d*x+e*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,1] =(i - ibegref)
                G[:,2] = 1
                G[:,3] = elev_temp
                G[:,4] = elev_temp**2

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),5))
                G[:-1,0] = y
                G[:-1,1] = x
                G[:,2] = 1
                G[:-1,3] = topo_clean
                G[:-1,4] = topo_clean*x

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                # print x0
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                # print pars - x0
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                # print pars - x0
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e=pars[4]
                print 'Remove ramp %f r  + %f az + %f + %f z +  %f z*az for date: %i'%(a,b,c,d,e,idates[l])

                # plot phase/elev
                funct = a*y + b*x + c + e*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,d*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                G[:,2] = 1
                G[:,3] = elev_temp
                G[:,4] = elev_temp
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,1] =(i - ibegref)
                    G[i*ncol:(i+1)*ncol,4] *= (i - ibegref)


                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),6))
                G[:-1,0] = y
                G[:-1,1] = x
                G[:,2] = 1
                G[:-1,3] = topo_clean*x
                G[:-1,4] = topo_clean
                G[:-1,5] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e=pars[4]; f=pars[5]
                print 'Remove ramp %f r  + %f az + %f +  %f z*az + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,idates[l])

                # plot phase/elev
                funct = a*y + b*x + c + d*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.1, alpha=0.01, rasterized=True)
                ax.plot(x,e*x+f*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,2] = 1
                G[:,3] = elev_temp
                G[:,4] = elev_temp
                G[:,5] = elev_temp**2
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,1] =(i - ibegref)
                    G[i*ncol:(i+1)*ncol,3] *= (i - ibegref)

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(nlign,ncol)

      elif order==4:
        if radar is None:
            G=np.zeros((len(data),4))
            G[:-1,0] = y
            G[:-1,1] = x
            G[:-1,2] = y*x
            G[:,3] = 1

            # ramp inversion
            x0 = lst.lstsq(G,data)[0]
            _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
            _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
            pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
            print 'Remove ramp %f r %f az  + %f r*az + %f for date: %i'%(a,b,c,d,idates[l])

            # build total G matrix
            G=np.zeros((len(los),4))
            for i in xrange(nlign):
                G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                G[i*ncol:(i+1)*ncol,1] = i - ibegref
                G[i*ncol:(i+1)*ncol,2] = (i-ibegref) * (np.arange((ncol))-jbegref)
            G[:,3] = 1

            res = los - np.dot(G,pars)
            rms = np.sqrt(np.nanmean(res**2))
            print 'RMS:', rms

            ramp = np.dot(G,pars).reshape(nlign,ncol)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),5))
                G[:-1,0] = y
                G[:-1,1] = x
                G[:-1,2] = y*x
                G[:,3] = 1
                G[:-1,4] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]

                print 'Remove ramp %f r, %f az  + %f r*az + %f + %f z for date: %i'%(a,b,c,d,e,idates[l])

                # plot phase/elev
                funct = a*y + b*x + c*x*y + d
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,1] = i - ibegref
                    G[i*ncol:(i+1)*ncol,2] = (i-ibegref) * (np.arange((ncol))-jbegref)
                G[:,3] = 1
                G[:,4] = elev_temp


                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(nlign,ncol)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),6))
                G[:-1,0] = y
                G[:-1,1] = x
                G[:-1,2] = y*x
                G[:,3] = 1
                G[:-1,4] = topo_clean
                G[:-1,5] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]

                print 'Remove ramp %f r, %f az  + %f r*az + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,idates[l])

                # plot phase/elev
                funct = a*y + b*x + c*x*y + d
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,e*x+f*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,1] = i - ibegref
                    G[i*ncol:(i+1)*ncol,2] = (i-ibegref) * (np.arange((ncol))-jbegref)
                G[:,3] = 1
                G[:,4] = elev_temp
                G[:,5] = elev_temp**2

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),6))
                G[:-1,0] = y
                G[:-1,1] = x
                G[:-1,2] = y*x
                G[:,3] = 1
                G[:-1,4] = topo_clean
                G[:-1,5] = topo_clean*x

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]

                print 'Remove ramp %f r, %f az  + %f r*az + %f + %f z + %f az*z for date: %i'%(a,b,c,d,e,f,idates[l])

                # plot phase/elev
                funct = a*y + b*x + c*x*y + d + f*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,3] = 1
                G[:,4] = elev_temp
                G[:,5] = elev_temp
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,1] = i - ibegref
                    G[i*ncol:(i+1)*ncol,2] = (i-ibegref) * (np.arange((ncol))-jbegref)
                    G[i*ncol:(i+1)*ncol,5] *=  (i - ibegref)


                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),6))
                G[:-1,0] = y
                G[:-1,1] = x
                G[:-1,2] = y*x
                G[:,3] = 1
                G[:-1,4] = topo_clean*x
                G[:-1,5] = topo_clean
                G[:-1,6] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]

                print 'Remove ramp %f r, %f az  + %f r*az + %f + + %f az*z +  %f z + %f z**2  for date: %i'%(a,b,c,d,e,f,g,idates[l])

                # plot phase/elev
                funct = a*y + b*x + c*x*y + d + e*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,f*x+g*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                G[:,3] = 1
                G[:,4] = elev_temp
                G[:,5] = elev_temp
                G[:,6] = elev_temp**2
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,1] = i - ibegref
                    G[i*ncol:(i+1)*ncol,2] = (i-ibegref) * (np.arange((ncol))-jbegref)
                    G[i*ncol:(i+1)*ncol,4] *=  (i - ibegref)

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(nlign,ncol)

      elif order==5:

        if radar is None:
            G=np.zeros((len(data),4))
            G[:-1,0] = y**2
            G[:-1,1] = y
            G[:-1,2] = x
            G[:,3] = 1

            # ramp inversion
            x0 = lst.lstsq(G,data)[0]
            _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
            _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
            pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
            print 'Remove ramp %f r**2 %f r  + %f az + %f for date: %i'%(a,b,c,d,idates[l])

            # build total G matrix
            G=np.zeros((len(los),4))
            for i in xrange(nlign):
                G[i*ncol:(i+1)*ncol,0] = (np.arange((ncol)) - jbegref)**2
                G[i*ncol:(i+1)*ncol,1] = np.arange((ncol)) - jbegref
                G[i*ncol:(i+1)*ncol,2] =(i - ibegref)
            G[:,3] = 1


            res = los - np.dot(G,pars)
            rms = np.sqrt(np.nanmean(res**2))
            print 'RMS:', rms

            ramp = np.dot(G,pars).reshape(nlign,ncol)

        else:

            if (ivar==0 and nfit==0):

                G=np.zeros((len(data),5))
                G[:-1,0] = y**2
                G[:-1,1] = y
                G[:-1,2] = x
                G[:,3] = 1
                G[:-1,4] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
                print 'Remove ramp %f r**2, %f r  + %f az + %f + %f z for date: %i'%(a,b,c,d,e,idates[l])

                # plot phase/elev
                funct = a*y**2 + b*y + c*x + d
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = (np.arange((ncol)) - jbegref)**2
                    G[i*ncol:(i+1)*ncol,1] = np.arange((ncol)) -  jbegref
                    G[i*ncol:(i+1)*ncol,2] =(i - ibegref)
                G[:,3] = 1
                G[:,4] = elev_temp


                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(nlign,ncol)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),6))
                G[:-1,0] = y**2
                G[:-1,1] = y
                G[:-1,2] = x
                G[:,3] = 1
                G[:-1,4] = topo_clean
                G[:-1,5] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
                print 'Remove ramp %f r**2, %f r  + %f az + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,idates[l])

                # plot phase/elev
                funct = a*y**2 + b*y + c*x + d
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,e*x+f*x**2,'-r', lw =4.)


                # build total G matrix
                G=np.zeros((len(los),6))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = (np.arange((ncol)) - jbegref)**2
                    G[i*ncol:(i+1)*ncol,1] = np.arange((ncol)) -  jbegref
                    G[i*ncol:(i+1)*ncol,2] =(i - ibegref)
                G[:,3] = 1
                G[:,4] = elev_temp
                G[:,5] = elev_temp**2

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)

            elif (ivar==1 and nfit==0):

                G=np.zeros((len(data),6))
                G[:-1,0] = y**2
                G[:-1,1] = y
                G[:-1,2] = x
                G[:,3] = 1
                G[:-1,4] = topo_clean
                G[:-1,5] = topo_clean*x

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
                print 'Remove ramp %f r**2, %f r  + %f az + %f + %f z + %f z*az for date: %i'%(a,b,c,d,e,f,idates[l])

                # plot phase/elev
                funct = a*y**2 + b*y + c*x + d + f*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,3] = 1
                G[:,4] = elev_temp
                G[:,5] = elev_temp
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = (np.arange((ncol)) - jbegref)**2
                    G[i*ncol:(i+1)*ncol,1] = np.arange((ncol)) -  jbegref
                    G[i*ncol:(i+1)*ncol,2] =(i - ibegref)
                    G[i*ncol:(i+1)*ncol,5] *= (i - ibegref)

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)


            elif (ivar==1 and nfit==1):

                G=np.zeros((len(data),7))
                G[:-1,0] = y**2
                G[:-1,1] = y
                G[:-1,2] = x
                G[:,3] = 1
                G[:-1,4] = topo_clean*x
                G[:-1,5] = topo_clean
                G[:-1,6] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]
                print 'Remove ramp %f r**2, %f r  + %f az + %f + + %f z*az + %f z +%f z**2 for date: %i'%(a,b,c,d,e,f,g,idates[l])

                # plot phase/elev
                funct = a*y**2 + b*y + c*x + d + e*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,f*x+g*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                G[:,3] = 1
                G[:,4] = elev_temp
                G[:,5] = elev_temp
                G[:,6] = elev_temp**2
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = (np.arange((ncol)) - jbegref)**2
                    G[i*ncol:(i+1)*ncol,1] = np.arange((ncol)) -  jbegref
                    G[i*ncol:(i+1)*ncol,2] =(i - ibegref)
                    G[i*ncol:(i+1)*ncol,4] *= (i - ibegref)

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(nlign,ncol)

            else:
                pass

      elif order==6:
        if radar is None:
            G=np.zeros((len(data),4))
            G[:-1,0] = x**2
            G[:-1,1] = x
            G[:-1,2] = y
            G[:,3] = 1

            # ramp inversion
            x0 = lst.lstsq(G,data)[0]
            _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
            _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
            pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
            print 'Remove ramp %f az**2 %f az  + %f r + %f for date: %i'%(a,b,c,d,idates[l])

            # build total G matrix
            G=np.zeros((len(los),4))
            for i in xrange(nlign):
                G[i*ncol:(i+1)*ncol,0] = (i - ibegref)**2
                G[i*ncol:(i+1)*ncol,1] = (i - ibegref)
                G[i*ncol:(i+1)*ncol,2] = np.arange((ncol)) - jbegref
            G[:,3] = 1


            res = los - np.dot(G,pars)
            rms = np.sqrt(np.nanmean(res**2))
            print 'RMS:', rms

            ramp = np.dot(G,pars).reshape(nlign,ncol)

        else:
            if (ivar==0 and nfit==0) :
                G=np.zeros((len(data),5))
                G[:-1,0] = x**2
                G[:-1,1] = x
                G[:-1,2] = y
                G[:,3] = 1
                G[:-1,4] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
                print 'Remove ramp %f az**2, %f az  + %f r + %f + %f z for date: %i'%(a,b,c,d,e,idates[l])

                # plot phase/elev
                funct = a*x**2 + b*x + c*y + d
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = (i - ibegref)**2
                    G[i*ncol:(i+1)*ncol,1] = i -  ibegref
                    G[i*ncol:(i+1)*ncol,2] = np.arange((ncol)) - jbegref
                G[:,3] = 1
                G[:,4] = elev_temp


                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(nlign,ncol)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),6))
                G[:-1,0] = x**2
                G[:-1,1] = x
                G[:-1,2] = y
                G[:,3] = 1
                G[:-1,4] = topo_clean
                G[:-1,5] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
                print 'Remove ramp %f az**2, %f az  + %f r + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,idates[l])

                # plot phase/elev
                funct = a*x**2 + b*x + c*y + d
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,e*x+f*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = (i - ibegref)**2
                    G[i*ncol:(i+1)*ncol,1] = i -  ibegref
                    G[i*ncol:(i+1)*ncol,2] = np.arange((ncol)) - jbegref
                G[:,3] = 1
                G[:,4] = elev_temp
                G[:,5] = elev_temp**2

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),6))
                G[:-1,0] = x**2
                G[:-1,1] = x
                G[:-1,2] = y
                G[:,3] = 1
                G[:-1,4] = topo_clean
                G[:-1,5] = topo_clean*x

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
                print 'Remove ramp %f az**2, %f az  + %f r + %f + %f z + %f z*az for date: %i'%(a,b,c,d,e,f,idates[l])

                # plot phase/elev
                funct = a*x**2 + b*x + c*y + d + f*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,3] = 1
                G[:,4] = elev_temp
                G[:,5] = elev_temp
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = (i - ibegref)**2
                    G[i*ncol:(i+1)*ncol,1] = i -  ibegref
                    G[i*ncol:(i+1)*ncol,2] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,5] *= (i - ibegref)


                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),7))
                G[:-1,0] = x**2
                G[:-1,1] = x
                G[:-1,2] = y
                G[:,3] = 1
                G[:-1,4] = topo_clean*x
                G[:-1,5] = topo_clean
                G[:-1,6] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g=pars[6]
                print 'Remove ramp %f az**2, %f az  + %f r + %f + %f z*az + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,g,idates[l])

                # plot phase/elev
                funct = a*x**2 + b*x + c*y + d + e*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,f*x+g*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                G[:,3] = 1
                G[:,4] = elev_temp
                G[:,5] = elev_temp
                G[:,6] = elev_temp**2
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = (i - ibegref)**2
                    G[i*ncol:(i+1)*ncol,1] = i -  ibegref
                    G[i*ncol:(i+1)*ncol,2] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,4] *= (i - ibegref)

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(nlign,ncol)


      elif order==7:
        if radar is None:
            G=np.zeros((len(data),5))
            G[:-1,0] = x**2
            G[:-1,1] = x
            G[:-1,2] = y**2
            G[:-1,3] = y
            G[:,4] = 1

            # ramp inversion
            x0 = lst.lstsq(G,data)[0]
            _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
            _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
            pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
            print 'Remove ramp %f az**2 %f az  + %f r**2 + %f r + %f for date: %i'%(a,b,c,d,e,idates[l])

            # build total G matrix
            G=np.zeros((len(los),5))
            for i in xrange(nlign):
                G[i*ncol:(i+1)*ncol,0] = (i - ibegref)**2
                G[i*ncol:(i+1)*ncol,1] = i - ibegref
                G[i*ncol:(i+1)*ncol,2] = (np.arange((ncol)) - jbegref)**2
                G[i*ncol:(i+1)*ncol,3] = np.arange((ncol)) - jbegref
            G[:,4] = 1


            res = los - np.dot(G,pars)
            rms = np.sqrt(np.nanmean(res**2))
            print 'RMS:', rms

            ramp = np.dot(G,pars).reshape(nlign,ncol)

        else:
            if (ivar==0 and nfit ==0):
                G=np.zeros((len(data),6))
                G[:-1,0] = x**2
                G[:-1,1] = x
                G[:-1,2] = y**2
                G[:-1,3] = y
                G[:,4] = 1
                G[:-1,5] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
                print 'Remove ramp %f az**2, %f az  + %f r**2 + %f r + %f + %f z for date: %i'%(a,b,c,d,e,f,idates[l])

                # plot phase/elev
                funct = a*x**2 + b*x + c*y**2 + d*y + e
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,f*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = (i - ibegref)**2
                    G[i*ncol:(i+1)*ncol,1] = i - ibegref
                    G[i*ncol:(i+1)*ncol,2] = (np.arange((ncol)) - jbegref)**2
                    G[i*ncol:(i+1)*ncol,3] = np.arange((ncol)) - jbegref
                G[:,4] = 1
                G[:,5] = elev_temp


                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(nlign,ncol)

            if (ivar==0 and nfit ==1):
                G=np.zeros((len(data),7))
                G[:-1,0] = x**2
                G[:-1,1] = x
                G[:-1,2] = y**2
                G[:-1,3] = y
                G[:,4] = 1
                G[:-1,5] = topo_clean
                G[:-1,6] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]
                print 'Remove ramp %f az**2, %f az  + %f r**2 + %f r + %f + %f z + %f z**2  for date: %i'%(a,b,c,d,e,f,g,idates[l])

                # plot phase/elev
                funct = a*x**2 + b*x + c*y**2 + d*y + e
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,f*x+g*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = (i - ibegref)**2
                    G[i*ncol:(i+1)*ncol,1] = i - ibegref
                    G[i*ncol:(i+1)*ncol,2] = (np.arange((ncol)) - jbegref)**2
                    G[i*ncol:(i+1)*ncol,3] = np.arange((ncol)) - jbegref
                G[:,4] = 1
                G[:,5] = elev_temp
                G[:,6] = elev_temp**2

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)

            elif (ivar==1 and nfit ==0):
                G=np.zeros((len(data),7))
                G[:-1,0] = x**2
                G[:-1,1] = x
                G[:-1,2] = y**2
                G[:-1,3] = y
                G[:,4] = 1
                G[:-1,5] = topo_clean
                G[:-1,6] = topo_clean*x

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g=pars[6]
                print 'Remove ramp %f az**2, %f az  + %f r**2 + %f r + %f + %f z + %f az*z for date: %i'%(a,b,c,d,e,f,g,idates[l])

                # plot phase/elev
                funct = a*x**2 + b*x + c*y**2 + d*y + e + g*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,f*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                G[:,4] = 1
                G[:,5] = elev_temp
                G[:,6] = elev_temp
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = (i - ibegref)**2
                    G[i*ncol:(i+1)*ncol,1] = i - ibegref
                    G[i*ncol:(i+1)*ncol,2] = (np.arange((ncol)) - jbegref)**2
                    G[i*ncol:(i+1)*ncol,3] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,6] *= (i - ibegref)


                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),8))
                G[:-1,0] = x**2
                G[:-1,1] = x
                G[:-1,2] = y**2
                G[:-1,3] = y
                G[:,4] = 1
                G[:-1,5] = topo_clean*x
                G[:-1,6] = topo_clean
                G[:,7] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g=pars[6]; h=pars[7]
                print 'Remove ramp %f az**2, %f az  + %f r**2 + %f r + %f +  %f az*z + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,g,h,idates[l])

                # plot phase/elev
                funct = a*x**2 + b*x + c*y**2 + d*y + e + f*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,g*x+h*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),8))
                G[:,4] = 1
                G[:,5] = elev_temp
                G[:,6] = elev_temp
                G[:,7] = elev_temp**2
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = (i - ibegref)**2
                    G[i*ncol:(i+1)*ncol,1] = i - ibegref
                    G[i*ncol:(i+1)*ncol,2] = (np.arange((ncol)) - jbegref)**2
                    G[i*ncol:(i+1)*ncol,3] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,5] *= (i - ibegref)

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(nlign,ncol)

      elif order==8:
        if radar is None:
            G=np.zeros((len(data),6))
            G[:-1,0] = x**3
            G[:-1,1] = x**2
            G[:-1,2] = x
            G[:-1,3] = y**2
            G[:-1,4] = y
            G[:,5] = 1

            # ramp inversion
            x0 = lst.lstsq(G,data)[0]
            _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
            _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
            pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
            print 'Remove ramp %f az**3 %f az**2  + %f az + %f r**2 + %f r + %f for date: %i'%(a,b,c,d,e,f,idates[l])

            # build total G matrix
            G=np.zeros((len(los),6))
            for i in xrange(nlign):
                G[i*ncol:(i+1)*ncol,0] = (i - ibegref)**3
                G[i*ncol:(i+1)*ncol,1] = (i - ibegref)**2
                G[i*ncol:(i+1)*ncol,2] =(i - ibegref)
                G[i*ncol:(i+1)*ncol,3] = (np.arange((ncol)) - jbegref)**2
                G[i*ncol:(i+1)*ncol,4] = (np.arange((ncol)) - jbegref)
            G[:,5] = 1


            res = los - np.dot(G,pars)
            rms = np.sqrt(np.nanmean(res**2))
            print 'RMS:', rms

            ramp = np.dot(G,pars).reshape(nlign,ncol)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),7))
                G[:-1,0] = x**3
                G[:-1,1] = x**2
                G[:-1,2] = x
                G[:-1,3] = y**2
                G[:-1,4] = y
                G[:,5] = 1
                G[:-1,6] = topo_clean

                # ramp inversion1
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]
                print 'Remove ramp %f az**3, %f az**2  + %f az + %f r**2 + %f r + %f + %f z for date: %i'%(a,b,c,d,e,f,g,idates[l])

                # plot phase/elev
                funct = a*x**3 + b*x**2 + c*x + d*y**2 + e*y + f
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,g*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = (i - ibegref)**3
                    G[i*ncol:(i+1)*ncol,1] = (i - ibegref)**2
                    G[i*ncol:(i+1)*ncol,2] =(i - ibegref)
                    G[i*ncol:(i+1)*ncol,3] = (np.arange((ncol)) - jbegref)**2
                    G[i*ncol:(i+1)*ncol,4] = np.arange((ncol)) - jbegref
                G[:,5] = 1
                G[:,6] = elev_temp


                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(nlign,ncol)

            if (ivar==0 and nfit==1):
                G=np.zeros((len(data),8))
                G[:-1,0] = x**3
                G[:-1,1] = x**2
                G[:-1,2] = x
                G[:-1,3] = y**2
                G[:-1,4] = y
                G[:,5] = 1
                G[:-1,6] = topo_clean
                G[:-1,7] = topo_clean**2

                # ramp inversion1
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]; h = pars[7]
                print 'Remove ramp %f az**3, %f az**2  + %f az + %f r**2 + %f r + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,g,h,idates[l])

                # plot phase/elev
                funct = a*x**3 + b*x**2 + c*x + d*y**2 + e*y + f
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,g*x+h*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),8))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = (i - ibegref)**3
                    G[i*ncol:(i+1)*ncol,1] = (i - ibegref)**2
                    G[i*ncol:(i+1)*ncol,2] =(i - ibegref)
                    G[i*ncol:(i+1)*ncol,3] = (np.arange((ncol)) - jbegref)**2
                    G[i*ncol:(i+1)*ncol,4] = np.arange((ncol)) - jbegref
                G[:,5] = 1
                G[:,6] = elev_temp
                G[:,7] = elev_temp**2

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)


            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),8))
                G[:-1,0] = x**3
                G[:-1,1] = x**2
                G[:-1,2] = x
                G[:-1,3] = y**2
                G[:-1,4] = y
                G[:,5] = 1
                G[:-1,6] = topo_clean
                G[:-1,7] = topo_clean*x

                # ramp inversion1
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]; h=pars[7]
                print 'Remove ramp %f az**3, %f az**2  + %f az + %f r**2 + %f r + %f + %f z + %f z*az for date: %i'%(a,b,c,d,e,f,g,h,idates[l])

                # plot phase/elev
                funct = a*x**3 + b*x**2 + c*x + d*y**2 + e*y + f + h*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,g*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),8))
                G[:,5] = 1
                G[:,6] = elev_temp
                G[:,7] = elev_temp
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = (i - ibegref)**3
                    G[i*ncol:(i+1)*ncol,1] = (i - ibegref)**2
                    G[i*ncol:(i+1)*ncol,2] =(i - ibegref)
                    G[i*ncol:(i+1)*ncol,3] = (np.arange((ncol)) - jbegref)**2
                    G[i*ncol:(i+1)*ncol,4] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,7] *= (i - ibegref)


                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),9))
                G[:-1,0] = x**3
                G[:-1,1] = x**2
                G[:-1,2] = x
                G[:-1,3] = y**2
                G[:-1,4] = y
                G[:,5] = 1
                G[:-1,6] = topo_clean*x
                G[:-1,7] = topo_clean
                G[:-1,8] = topo_clean

                # ramp inversion1
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]; h=pars[7]; i=pars[8]
                print 'Remove ramp %f az**3, %f az**2  + %f az + %f r**2 + %f r + %f z*az + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,g,h,i,idates[l])

                # plot phase/elev
                funct = a*x**3 + b*x**2 + c*x + d*y**2 + e*y + f + g*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,h*x+i*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),9))
                G[:,5] = 1
                G[:,6] = elev_temp
                G[:,7] = elev_temp
                G[:,8] = elev_temp**2
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = (i - ibegref)**3
                    G[i*ncol:(i+1)*ncol,1] = (i - ibegref)**2
                    G[i*ncol:(i+1)*ncol,2] =(i - ibegref)
                    G[i*ncol:(i+1)*ncol,3] = (np.arange((ncol)) - jbegref)**2
                    G[i*ncol:(i+1)*ncol,4] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,6] *= (i - ibegref)

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(nlign,ncol)

      elif order==9:
        if radar is None:
            G=np.zeros((len(data),5))
            G[:-1,0] = y
            G[:-1,1] = x
            G[:-1,2] = (y*x)**2
            G[:-1,3] = y*x
            G[:,4] = 1

            # ramp inversion
            x0 = lst.lstsq(G,data)[0]
            _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
            _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
            pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
            print 'Remove ramp %f r %f az  + %f r*az**2 + %f r*az + %f for date: %i'%(a,b,c,d,e,idates[l])

            # build total G matrix
            G=np.zeros((len(los),5))
            for i in xrange(nlign):
                G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                G[i*ncol:(i+1)*ncol,1] = i - ibegref
                G[i*ncol:(i+1)*ncol,2] = ((i-ibegref) * (np.arange((ncol))-jbegref))**2
                G[i*ncol:(i+1)*ncol,3] = (i-ibegref) * (np.arange((ncol))-jbegref)
            G[:,4] = 1

            res = los - np.dot(G,pars)
            rms = np.sqrt(np.nanmean(res**2))
            print 'RMS:', rms

            ramp = np.dot(G,pars).reshape(nlign,ncol)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),6))
                G[:-1,0] = y
                G[:-1,1] = x
                G[:-1,2] = (y*x)**2
                G[:-1,3] = y*x
                G[:,4] = 1
                G[:-1,5] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]

                print 'Remove ramp %f r, %f az  + %f (r*az)**2 + %f r*az + %f + %f z for date: %i'%(a,b,c,d,e,f,idates[l])

                # plot phase/elev
                funct = a*y + b*x + c*(x*y)**2 + d*x*y + e
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,1] = i - ibegref
                    G[i*ncol:(i+1)*ncol,2] = ((i-ibegref) * (np.arange((ncol))-jbegref))**2
                    G[i*ncol:(i+1)*ncol,3] = (i-ibegref) * (np.arange((ncol))-jbegref)
                G[:,4] = 1
                G[:,5] = elev_temp

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(nlign,ncol)

            if (ivar==0 and nfit==1):
                G=np.zeros((len(data),7))
                G[:-1,0] = y
                G[:-1,1] = x
                G[:-1,2] = (y*x)**2
                G[:-1,3] = y*x
                G[:,4] = 1
                G[:-1,5] = topo_clean
                G[:-1,6] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]

                print 'Remove ramp %f r, %f az  + %f (r*az)**2 + %f r*az + %f + %f z + %f z**2  for date: %i'%(a,b,c,d,e,f,g,idates[l])

                # plot phase/elev
                funct = a*y + b*x + c*(x*y)**2 + d*x*y + e
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,f*x+g*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,1] = i - ibegref
                    G[i*ncol:(i+1)*ncol,2] = ((i-ibegref) * (np.arange((ncol))-jbegref))**2
                    G[i*ncol:(i+1)*ncol,3] = (i-ibegref) * (np.arange((ncol))-jbegref)
                G[:,4] = 1
                G[:,5] = elev_temp
                G[:,6] = elev_temp**2

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),7))
                G[:-1,0] = y
                G[:-1,1] = x
                G[:-1,2] = (y*x)**2
                G[:-1,3] = y*x
                G[:,4] = 1
                G[:-1,5] = topo_clean
                G[:-1,6] = topo_clean*x

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5] ; g = pars[6]

                print 'Remove ramp %f r, %f az  + %f (r*az)**2 + %f r*az + %f + %f z + %f az*z for date: %i'%(a,b,c,d,e,f,g,idates[l])

                # plot phase/elev
                funct = a*y + b*x + c*(x*y)**2 + d*x*y + e + g*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,f*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                G[:,4] = 1
                G[:,5] = elev_temp
                G[:,6] = elev_temp
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,1] = i - ibegref
                    G[i*ncol:(i+1)*ncol,2] = ((i-ibegref) * (np.arange((ncol))-jbegref))**2
                    G[i*ncol:(i+1)*ncol,3] = (i-ibegref) * (np.arange((ncol))-jbegref)
                    G[i*ncol:(i+1)*ncol,6] *= (i - ibegref)

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(nlign,ncol)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),8))
                G[:-1,0] = y
                G[:1,1] = x
                G[:-1,2] = (y*x)**2
                G[:-1,3] = y*x
                G[:,4] = 1
                G[:-1,5] = topo_clean*x
                G[:-1,6] = topo_clean
                G[:-1,7] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5] ; g = pars[6]; h=pars[7]

                print 'Remove ramp %f r, %f az  + %f (r*az)**2 + %f r*az + %f + %f az*z + %f z + %f z**2  for date: %i'%(a,b,c,d,e,f,g,h,idates[l])

                # plot phase/elev
                funct = a*y + b*x + c*(x*y)**2 + d*x*y + e + f*topo_clean*x
                x = np.linspace(np.nanmin(topo_clean), np.nanmax(topo_clean), 100)
                ax.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax.plot(x,g*x+h*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),8))
                G[:,4] = 1
                G[:,5] = elev_temp
                G[:,6] = elev_temp
                G[:,7] = elev_temp**2
                for i in xrange(nlign):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbegref
                    G[i*ncol:(i+1)*ncol,1] = i - ibegref
                    G[i*ncol:(i+1)*ncol,2] = ((i-ibegref) * (np.arange((ncol))-jbegref))**2
                    G[i*ncol:(i+1)*ncol,3] = (i-ibegref) * (np.arange((ncol))-jbegref)
                    G[i*ncol:(i+1)*ncol,5] *= (i - ibegref)

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                print 'RMS:', rms

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(nlign,ncol)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(nlign,ncol)

      # flata = (los - np.dot(G,pars)).reshape(nlign,ncol)
      flata = los.reshape(nlign,ncol) - ramp - topo
      noramps = los.reshape(nlign,ncol) - ramp

      if radar is not None:
        del elev_temp

      return ramp, flata, topo, rms, noramps
 
    pix_az, pix_rg = np.indices((nlign,ncol))
    # if radar file just initialise figure
    if radar is not None:
      nfigure +=1
      fig = plt.figure(nfigure,figsize=(14,10))
    
    # if iteration = 0 or spatialiter > 0, then spatial estimation
    if (ii==0) or (spatialiter=='yes') :

      # Loop over the dates
      for l in xrange((N)):
        print

        # no estimation on the ref image set to zero 
        if l is not imref:

          # first clean los
          maps_temp = np.matrix.copy(maps[:,:,l]) - np.matrix.copy(models[:,:,l])

          maxlos,minlos=np.nanpercentile(maps_temp[ibegref:iendref,jbegref:jendref],perc_los),np.nanpercentile(maps_temp[ibegref:iendref,jbegref:jendref],100-perc_los)
          kk = np.nonzero(np.logical_or(maps_temp==0.,np.logical_or((maps_temp>maxlos),(maps_temp<minlos))))
          maps_temp[kk] = np.float('NaN')

          #noise_level=np.nanpercentile(maps_temp,65) - np.nanpercentile(maps_temp,35)
          #print 'Accepted noise level in the ramp optimisation:', noise_level

          itemp = ibegref
          for lign in xrange(ibegref,iendref,10):
              # find the begining of the image
              if np.isnan(np.nanmean(maps[lign:lign+10,:,l])):
                  itemp = lign
              else:
                  break

          if radar is not None:
              maxtopo,mintopo = np.nanpercentile(elev,perc_topo),np.nanpercentile(elev,100-perc_topo)
              # initialize plot
              ax = fig.add_subplot(4,int(N/4)+1,l+1)
          else:
              maxtopo,mintopo = 2, 0

          if rmsf is None:
              seuil_rms = 2

          # selection pixels
          index = np.nonzero(np.logical_and(elev<maxtopo,
              np.logical_and(elev>mintopo,
                  np.logical_and(mask_flat>seuil,
                  np.logical_and(~np.isnan(maps_temp),
                      np.logical_and(~np.isnan(rmsmap),
                      np.logical_and(~np.isnan(elev),
                      np.logical_and(rmsmap<seuil_rms,
                      np.logical_and(rmsmap>1.e-6,
                      np.logical_and(~np.isnan(maps_temp),
                      np.logical_and(pix_az>ibeg,
                      np.logical_and(pix_az<iend,
                      np.logical_and(pix_rg>jbeg,
                      np.logical_and(pix_rg<jend, 
                          slope>0.,
                          ))))))))
                      ))))))

          indexref = np.nonzero(np.logical_and(elev<maxtopo,
              np.logical_and(elev>mintopo,
                  np.logical_and(mask_flat>seuil,
                  np.logical_and(~np.isnan(maps_temp),
                      np.logical_and(~np.isnan(rmsmap),
                      np.logical_and(~np.isnan(elev),
                      np.logical_and(rmsmap<seuil_rms,
                      np.logical_and(rmsmap>1.e-6,
                      np.logical_and(~np.isnan(maps_temp),
                      np.logical_and(pix_az>refstart,
                      np.logical_and(pix_az<refend,
                      np.logical_and(pix_rg>jbeg,
                      np.logical_and(pix_rg<jend, 
                          slope>0.,
                          ))))))))
                      ))))))


          # extract coordinates for estimation
          temp = np.array(index).T
          x = temp[:,0]; y = temp[:,1]

          # clean maps
          los_clean = maps_temp[index].flatten()
          topo_clean = elev[index].flatten()
          rms_clean = rmsmap[index].flatten()
          
          try:
            los_ref = maps_temp[indexref].flatten()
            rms_ref = rmsmap[indexref].flatten()
            amp_ref = 1./rms_ref
            amp_ref = amp_ref/np.nanmax(amp_ref)
            print 'Ref area set to zero:', refstart,refend
            cst = np.nansum(los_ref*amp_ref) / np.nansum(amp_ref)
            print 'Average phase within ref area:', cst
            if np.isnan(cst):
              cst = 0

          except:
            cst = 0.
          
          # print itemp, iendref
          #4: ax+by+cxy+d 5: ax**2+bx+cy+d, 6: ay**2+by+cx+d, 7: ay**2+by+cx**2+dx+e, 8: ay**2+by+cx**3+dx**2+ex+f
          if flat>5 and iendref-itemp < .6*(iendref-ibegref):
              print 'Image too short in comparison to master, set flat to 5'
              temp_flat=5
          # elif flat>5 and iendref-itemp < ncol:
          #     print 'Lenght image inferior to width, set flat to 5'
          #     temp_flat=5
          else:
              temp_flat=flat

          if ivar>0 and iendref-itemp < .6*(iendref-ibegref):
            print
            print 'Image too short in comparison to master, set ivar to 0'
            ivar_temp=0
            nfit_temp=0
          else:
            ivar_temp=ivar
            nfit_temp=nfit

          # call ramp estim
          los = as_strided(maps[:,:,l]).flatten()
          samp = 1

          # print los,los_clean[::samp],topo_clean[::samp],x[::samp],y[::samp],temp_flat,rms_clean[::samp]
          maps_ramp[:,:,l], maps_flata[:,:,l], maps_topo[:,:,l], rms[l], maps_noramps[:,:,l] = estim_ramp(los,los_clean[::samp],topo_clean[::samp],x[::samp],y[::samp],temp_flat,rms_clean[::samp],nfit_temp, ivar_temp,cst)

          # set ramp to NaN to have ramp of the size of the images
          kk = np.nonzero(np.isnan(maps_flata[:,:,l]))
          ramp = as_strided(maps_ramp[:,:,l])
          ramp[kk] = float('NaN')
          topo = as_strided(maps_topo[:,:,l])
          topo[kk] = float('NaN')
          
          ## Refer data again 
          zone = as_strided(maps_flata[:,:,l])
          los_ref2 = zone[indexref]
          
          # weigth avera of the phase
          cst = np.nansum(los_ref2*amp_ref) / np.nansum(amp_ref)
          print 'Average phase within ref area, iter=2:', cst
          if np.isnan(cst):
            cst = 0.
          maps_ramp[:,:,l], maps_flata[:,:,l], maps_noramps[:,:,l] = maps_ramp[:,:,l] + cst, maps_flata[:,:,l] - cst, maps_noramps[:,:,l] - cst 
    
          del zone
          del los_clean
          del rms_clean
          del topo_clean

      # plot corrected ts
      nfigure +=1
      figd = plt.figure(nfigure,figsize=(14,10))
      figd.subplots_adjust(hspace=0.001,wspace=0.001)
      for l in xrange((N)):
          axd = figd.add_subplot(4,int(N/4)+1,l+1)
          caxd = axd.imshow(maps_flata[ibeg:iend,jbeg:jend,l],cmap=cmap,vmax=vmax,vmin=vmin)
          axd.set_title(idates[l],fontsize=6)
          setp(axd.get_xticklabels(), visible=False)
          setp(axd.get_yticklabels(), visible=False)
      setp(axd.get_xticklabels(), visible=False)
      setp(axd.get_yticklabels(), visible=False)
      figd.colorbar(caxd, orientation='vertical',aspect=10)
      figd.suptitle('Corrected time series maps')
      fig.tight_layout()
      figd.savefig('maps_flat.eps', format='EPS',dpi=150)

      if radar is not None:
          fig.savefig('phase-topo.eps', format='EPS',dpi=150)
          nfigure +=1
          figtopo = plt.figure(nfigure,figsize=(14,10))
          figtopo.subplots_adjust(hspace=.001,wspace=0.001)
          for l in xrange((N)):
              axtopo = figtopo.add_subplot(4,int(N/4)+1,l+1)
              caxtopo = axtopo.imshow(maps_topo[ibeg:iend,jbeg:jend,l]+maps_ramp[ibeg:iend,jbeg:jend,l],cmap=cmap,vmax=vmax,vmin=vmin)
              axtopo.set_title(idates[l],fontsize=6)
              setp(axtopo.get_xticklabels(), visible=False)
              setp(axtopo.get_yticklabels(), visible=False)
              setp(axtopo.get_xticklabels(), visible=False)
              setp(axtopo.get_yticklabels(), visible=False)
          figtopo.colorbar(caxtopo, orientation='vertical',aspect=10)
          figtopo.suptitle('Time series RAMPS+TOPO')
          fig.tight_layout()
          figtopo.savefig('tropo.eps', format='EPS',dpi=150)
          

      else:
          # plot corrected ts
          nfigure +=1
          figref = plt.figure(nfigure,figsize=(14,10))
          figref.subplots_adjust(hspace=0.001,wspace=0.001)
          for l in xrange((N)):
              axref = figref.add_subplot(4,int(N/4)+1,l+1)
              caxref = axref.imshow(maps_ramp[ibeg:iend,jbeg:jend,l],cmap=cmap,vmax=vmax,vmin=vmin)
              axref.set_title(idates[l],fontsize=6)
              setp(axref.get_xticklabels(), visible=False)
              setp(axref.get_yticklabels(), visible=False)
          setp(axref.get_xticklabels(), visible=False)
          setp(axref.get_yticklabels(), visible=False)
          figref.suptitle('Time series RAMPS')
          figref.colorbar(caxref, orientation='vertical',aspect=10)
          fig.tight_layout()
          figref.savefig('maps_ramps.eps', format='EPS',dpi=150)
          

    if plot=='yes':
        plt.show()
    plt.close('all')

    # save rms
    if (apsf=='no' and ii==0):
        # aps from rms
        print
        print 'Use RMS empirical estimation as uncertainties for time decomposition'
        inaps = np.copy(rms)
        print 'Set very low values to the 2 percentile to avoid overweighting...'
        # scale between 0 and 1 for threshold_rmsd
        maxaps = np.nanmax(inaps)
        inaps = inaps/maxaps
        minaps= np.nanpercentile(inaps,2)
        index = flatnonzero(inaps<minaps)
        inaps[index] = minaps
        np.savetxt('rms_empcor.txt', inaps.T)
        del rms

    ########################
    # TEMPORAL ITERATION N #
    ########################

    print
    print 'Time decomposition..'
    print

    # initialize aps for each images to 1
    aps = np.ones((N))
    n_aps = np.ones((N)).astype(int)
    print inaps

    # reiinitialize maps models
    models = np.zeros((nlign,ncol,N))

    if seasonal=='yes' or semianual=='yes' or inter=='yes' or vect != None:
        models_trends = np.zeros((nlign,ncol,N))
        models_detrends = np.zeros((nlign,ncol,N))


    # ligns = [2014,2157,1840,1960,1951]
    # cols = [100,117,843,189,43]
    # for i,j in zip(ligns,cols):

    for i in xrange(ibeg,iend,sampling):
        for j in xrange(jbeg,jend,sampling):
            #print j

            # Initialisation
            mdisp=np.ones((N))*float('NaN')
            #!!!!!
            disp = as_strided(maps_flata[i,j,:])
            # disp = as_strided(maps[i,j,:])
            # print disp

            k = np.flatnonzero(~np.isnan(disp)) # invers of isnan
            # do not take into account NaN data
            kk = len(k)
            tabx = dates[k]
            taby = disp[k]
            bp = base[k]

            # Inisilize m to zero
            m = np.zeros((M))
            sigmam = np.ones((M))*float('NaN')

            if kk > N/6:
                G=np.zeros((kk,M))
                # Build G family of function k1(t),k2(t),...,kn(t): #
                #                                                   #
                #           |k1(0) .. kM(0)|                        #
                # Gfamily = |k1(1) .. kM(1)|                        #
                #           |..    ..  ..  |                        #
                #           |k1(N) .. kM(N)|                        #
                #                                                   #

                rmsd = maxrmsd + 1

                if inter=='yes' and iteration is True:
                    Glin=np.zeros((kk,2+Mker))
                    for l in xrange((2)):
                        Glin[:,l]=basis[l].g(tabx)
                    for l in xrange((Mker)):
                        Glin[:,2+l]=kernels[l].g(k)

                    mt,sigmamt = consInvert(Glin,taby,inaps[k],cond=rcond)

                    # compute rmsd
                    mdisp[k] = np.dot(Glin,mt)
                    # sum sur toutes les dates
                    # rmsd = np.sum(abs((disp[k] - mdisp[k])/inaps[k]))/kk  S
                    rmsd = np.sqrt(np.sum(pow((disp[k] - mdisp[k]),2))/kk)
                    # print i,j,rmsd,maxrmsd

                G=np.zeros((kk,M))
                for l in xrange((Mbasis)):
                    G[:,l]=basis[l].g(tabx)
                for l in xrange((Mker)):
                    G[:,Mbasis+l]=kernels[l].g(k)

                # if only ref + seasonal: ref + cos + sin
                #print rmsd
                if rmsd >= maxrmsd or inter!='yes':
                    mt,sigmamt = consInvert(G,taby,inaps[k],cond=rcond,ineq=ineq)

                # rebuild full vectors
                if Mker>0:
                    m[Mbasis:],sigmam[Mbasis:] = mt[-Mker:],sigmamt[-Mker:]
                    m[:mt.shape[0]-Mker],sigmam[:mt.shape[0]-Mker] = mt[:-Mker],sigmamt[:-Mker]
                else:
                    sigmam[:mt.shape[0]] = sigmamt
                    m[:mt.shape[0]] = mt

                # save m
                for l in xrange((Mbasis)):
                    basis[l].m[i-ibeg,j-jbeg] = m[l]
                    basis[l].sigmam[i-ibeg,j-jbeg] = sigmam[l]

                for l in xrange((Mker)):
                    kernels[l].m[i-ibeg,j-jbeg] = m[Mbasis+l]
                    kernels[l].sigmam[i-ibeg,j-jbeg] = sigmam[Mbasis+l]

                # forward model in original order
                mdisp[k] = np.dot(G,m)

                # compute aps for each dates
                # aps_tmp = pow((disp[k]-mdisp[k])/inaps[k],2)
                aps_tmp = abs((disp[k]-mdisp[k]))/inaps[k]

                # # remove NaN value for next iterations (but normally no NaN?)
                index = np.flatnonzero(np.logical_or(np.isnan(aps_tmp),aps_tmp==0))
                aps_tmp[index] = 1.0 # 1 is a bad misfit

                # save total aps of the map
                aps[k] = aps[k] + aps_tmp

                # count number of pixels per dates
                n_aps[k] = n_aps[k] + 1.0

                # save new aps for each maps
                # maps_aps[i,j,k] = aps_tmp

                # fill maps models
                models[i,j,:] = mdisp

                # Build seasonal and linear models
                if inter=='yes':
                    models_detrends[i,j,k] = models_detrends[i,j,k] + np.dot(G[:,indexinter],m[indexinter])

                if inter=='yes':
                    models_trends[i,j,k] = models_trends[i,j,k] + np.dot(G[:,indexinter],m[indexinter])
                if vect != None:
                    models_trends[i,j,k] = models_trends[i,j,k] + np.dot(G[:,indexvect],m[indexvect])

    # convert aps in rad
    aps = aps/n_aps
    # aps = np.sqrt(abs(aps/n_aps))
    minaps= np.nanpercentile(aps,2)
    index = flatnonzero(aps<minaps)
    aps[index] = minaps

    print
    print 'Dates      APS     # of points'
    for l in xrange(N):
        print idates[l], aps[l], n_aps[l]
    np.savetxt('aps_{}.txt'.format(ii), aps.T, fmt=('%.6f'))
    # set apsf is yes for iteration
    apsf=='yes'
    # update aps for next iterations
    inaps = np.copy(aps)

# del maps_aps

#######################################################
# Save new cubes
#######################################################

# create new cube
cube_flata = maps_flata[ibeg:iend,jbeg:jend,:].flatten()
cube_noramps = maps_noramps[ibeg:iend,jbeg:jend,:].flatten()

fid = open('depl_cumule_flat', 'wb')
cube_flata.flatten().astype('float32').tofile(fid)
fid.close()

if fulloutput=='yes':
    if (seasonal=='yes' or semianual=='yes') and (vect != None or inter=='yes'):
        fid = open('depl_cumule_dseas', 'wb')
        (maps_flata - models_trends).flatten().astype('float32').tofile(fid)
        fid.close()

    if inter=='yes':
        fid = open('depl_cumule_dtrend', 'wb')
        (maps_flata - models_detrends).flatten().astype('float32').tofile(fid)
        fid.close()

    if flat>0:
        fid = open('depl_cumule_noramps', 'wb')
        cube_noramps.flatten().astype('float32').tofile(fid)
        fid.close()

# # save APS
# print
# print 'Saving APS in liste_images_aps.txt'
# np.savetxt('liste_images_aps.txt', np.vstack([idates,dates,aps]).T,header='#dates #dates_dec #aps', fmt=('%i', '%.6f', '%.6f'))

# create MAPS directory to save .r4
if fulloutput=='yes':
    outdir = './MAPS/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

# plot displacements models and residuals
nfigure +=1
figres = plt.figure(nfigure,figsize=(14,10))
figres.subplots_adjust(hspace=.001,wspace=0.001)

nfigure +=1
fig = plt.figure(nfigure,figsize=(14,10))
fig.subplots_adjust(hspace=.001,wspace=0.01)

nfigure +=1
figall = plt.figure(nfigure,figsize=(20,9))
figall.subplots_adjust(hspace=0.00001,wspace=0.001)

nfigure +=1
figclr = plt.figure(nfigure)

# plot color map
ax = figclr.add_subplot(1,1,1)
cax = ax.imshow(maps[:,:,-1],cmap=cmap,vmax=vmax,vmin=vmin)
setp( ax.get_xticklabels(), visible=False)
cbar = figclr.colorbar(cax, orientation='horizontal',aspect=5)
figclr.savefig('colorscale.eps', format='EPS',dpi=150)

# vmax = np.abs([np.nanmedian(data) + 2*nanstd(data),np.nanmedian(data) - 2*nanstd(data)]).max()
# vmin = -vmax

for l in xrange((N)):
    data = as_strided(maps[ibeg:iend,jbeg:jend,l])
    if Mker>0:
        data_flat = as_strided(maps_flata[ibeg:iend,jbeg:jend,l])- as_strided(kernels[0].m[:,:]) - as_strided(basis[0].m[:,:])
        model = as_strided(models[ibeg:iend,jbeg:jend,l]) - as_strided(basis[0].m[:,:]) - as_strided(kernels[0].m[:,:])
    else:
        data_flat = as_strided(maps_flata[ibeg:iend,jbeg:jend,l]) - as_strided(basis[0].m[:,:])
        model = as_strided(models[ibeg:iend,jbeg:jend,l]) - as_strided(basis[0].m[:,:])

    res = data_flat - model
    ramp = as_strided(maps_ramp[ibeg:iend,jbeg:jend,l])
    tropo = as_strided(maps_topo[ibeg:iend,jbeg:jend,l])

    ax = fig.add_subplot(4,int(N/4)+1,l+1)
    axres = figres.add_subplot(4,int(N/4)+1,l+1)

    axall = figall.add_subplot(6,N,l+1)
    axall.imshow(data,cmap=cmap,vmax=vmax,vmin=vmin)
    axall.set_title(idates[l],fontsize=6)
    setp(axall.get_xticklabels(), visible=False)
    setp(axall.get_yticklabels(), visible=False)
    if l==0:
        axall.set_ylabel('DATA')
    axall = figall.add_subplot(6,N,l+1+N)
    axall.imshow(ramp,cmap=cmap,vmax=vmax,vmin=vmin)
    setp(axall.get_xticklabels(), visible=False)
    setp(axall.get_yticklabels(), visible=False)
    if l==0:
        axall.set_ylabel('RAMP')
    axall = figall.add_subplot(6,N,l+1+2*N)
    axall.imshow(tropo,cmap=cmap,vmax=vmax,vmin=vmin)
    setp(axall.get_xticklabels(), visible=False)
    setp(axall.get_yticklabels(), visible=False)
    if l==0:
        axall.set_ylabel('TROP0')
    axall = figall.add_subplot(6,N,l+1+3*N)
    axall.imshow(data_flat,cmap=cmap,vmax=vmax,vmin=vmin)
    setp(axall.get_xticklabels(), visible=False)
    setp(axall.get_yticklabels(), visible=False)
    if l==0:
        axall.set_ylabel('FLATTEN DATA')
    axall = figall.add_subplot(6,N,l+1+4*N)
    axall.imshow(model,cmap=cmap,vmax=vmax,vmin=vmin)
    setp(axall.get_xticklabels(), visible=False)
    setp(axall.get_yticklabels(), visible=False)
    if l==0:
        axall.set_ylabel('MODEL')
    axall = figall.add_subplot(6,N,l+1+5*N)
    axall.imshow(res,cmap=cmap,vmax=vmax,vmin=vmin)
    setp(axall.get_xticklabels(), visible=False)
    setp(axall.get_yticklabels(), visible=False)
    if l==0:
        axall.set_ylabel('RES')

    cax = ax.imshow(model,cmap=cmap,vmax=vmax,vmin=vmin)
    caxres = axres.imshow(res,cmap=cmap,vmax=vmax,vmin=vmin)

    ax.set_title(idates[l],fontsize=6)
    axres.set_title(idates[l],fontsize=6)

    setp(ax.get_xticklabels(), visible=False)
    setp(ax.get_yticklabels(), visible=False)

    setp(axres.get_xticklabels(), visible=False)
    setp(axres.get_yticklabels(), visible=False)

    fig.tight_layout()

    # ############
    # # SAVE .R4 #
    # ############

    # save flatten maps
    if fulloutput=='yes':

        if geotiff is not None:

            ds = driver.Create(outdir+'{}_flat.tif'.format(idates[l]), jend-jbeg, iend-ibeg, 1, gdal.GDT_Float32)
            band = ds.GetRasterBand(1)
            band.WriteArray(data_flat)
            ds.SetGeoTransform(gt)
            ds.SetProjection(proj)
            band.FlushCache()

            ds = driver.Create(outdir+'{}_ramp_tropo.tif'.format(idates[l]), jend-jbeg, iend-ibeg, 1, gdal.GDT_Float32)
            band = ds.GetRasterBand(1)
            band.WriteArray(ramp+tropo)
            ds.SetGeoTransform(gt)
            ds.SetProjection(proj)
            band.FlushCache()

            ds = driver.Create(outdir+'{}_model.tif'.format(idates[l]), jend-jbeg, iend-ibeg, 1, gdal.GDT_Float32)
            band = ds.GetRasterBand(1)
            band.WriteArray(model)
            ds.SetGeoTransform(gt)
            ds.SetProjection(proj)
            band.FlushCache()

            # ds = driver.Create(outdir+'{}_res.tif'.format(idates[l]), jend-jbeg, iend-ibeg, 1, gdal.GDT_Float32)
            # band = ds.GetRasterBand(1)
            # band.WriteArray(res)
            # ds.SetGeoTransform(gt)
            # ds.SetProjection(proj)
            # band.FlushCache()

        else:

            fid = open(outdir+'{}_flat.r4'.format(idates[l]), 'wb')
            data_flat.flatten().astype('float32').tofile(fid)
            fid.close()

            # save ramp maps
            fid = open(outdir+'{}_ramp_tropo.r4'.format(idates[l]), 'wb')
            (ramp+tropo).flatten().astype('float32').tofile(fid)
            fid.close()

            # save model maps
            fid = open(outdir+'{}_model.r4'.format(idates[l]), 'wb')
            model.flatten().astype('float32').tofile(fid)
            fid.close()

            # save residual maps
            fid = open(outdir+'{}_res.r4'.format(idates[l]), 'wb')
            res.flatten().astype('float32').tofile(fid)
            # fid.close()


fig.suptitle('Time series models')
figres.suptitle('Time series residuals')
figall.suptitle('Time series inversion')
fig.savefig('models.eps', format='EPS',dpi=150)
figres.savefig('residuals.eps', format='EPS',dpi=150)
figall.savefig('timeseries.eps', format='EPS',dpi=150)
if plot=='yes':
    plt.show()
plt.close('all')

#######################################################
# Save functions in binary file
#######################################################


if geotiff is not None:
    for l in xrange((Mbasis)):
        ds = driver.Create('{}_coeff.tif'.format(basis[l].reduction), jend-jbeg, iend-ibeg, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(basis[l].m)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

        ds = driver.Create('{}_sigcoeff.tif'.format(basis[l].reduction), jend-jbeg, iend-ibeg, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(basis[l].sigmam)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

    for l in xrange((Mker)):
        ds = driver.Create('{}_coeff.tif'.format(kernels[l].reduction), jend-jbeg, iend-ibeg, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(kernels[l].m)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

        ds = driver.Create('{}_sigcoeff.tif'.format(kernels[l].reduction), jend-jbeg, iend-ibeg, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(kernels[l].sigmam)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

else:
    for l in xrange((Mbasis)):
        fid = open('{}_coeff.r4'.format(basis[l].reduction), 'wb')
        basis[l].m.flatten().astype('float32').tofile(fid)
        fid.close()
        fid = open('{}_sigcoeff.r4'.format(basis[l].reduction), 'wb')
        basis[l].sigmam.flatten().astype('float32').tofile(fid)
        fid.close()
    for l in xrange((Mker)):
        fid = open('{}_coeff.r4'.format(kernels[l].reduction), 'wb')
        kernels[l].m.flatten().astype('float32').tofile(fid)
        fid.close()
        fid = open('{}_sigcoeff.r4'.format(kernels[l].reduction), 'wb')
        kernels[l].sigmam.flatten().astype('float32').tofile(fid)
        fid.close()


#######################################################
# Compute Amplitude and phase seasonal
#######################################################

if seasonal == 'yes':
    cosine = as_strided(basis[indexseas].m)
    sine = as_strided(basis[indexseas+1].m)
    amp = np.sqrt(cosine**2+sine**2)
    phi = np.arctan2(sine,cosine)

    sigcosine = as_strided(basis[indexseas].sigmam)
    sigsine = as_strided(basis[indexseas+1].sigmam)
    sigamp = np.sqrt(sigcosine**2+sigsine**2)
    sigphi = (sigcosine*abs(sine)+sigsine*abs(cosine))/(sigcosine**2+sigsine**2)

    if geotiff is not None:
        ds = driver.Create('ampwt_coeff.tif', jend-jbeg, iend-ibeg, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(amp)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

        ds = driver.Create('ampwt_sigcoeff.tif', jend-jbeg, iend-ibeg, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(sigamp)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

    else:
        fid = open('ampwt_coeff.r4', 'wb')
        amp.flatten().astype('float32').tofile(fid)
        fid.close()

        fid = open('ampwt_sigcoeff.r4', 'wb')
        sigamp.flatten().astype('float32').tofile(fid)
        fid.close()

    if geotiff is not None:
        ds = driver.Create('phiwt_coeff.tif', jend-jbeg, iend-ibeg, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(phi)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

        ds = driver.Create('phiwt_sigcoeff.tif', jend-jbeg, iend-ibeg, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(sigphi)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

    else:
        fid = open('phiwt_coeff.r4', 'wb')
        phi.flatten().astype('float32').tofile(fid)
        fid.close()

        fid = open('phiwt_sigcoeff.r4', 'wb')
        sigphi.flatten().astype('float32').tofile(fid)
        fid.close()

if semianual == 'yes':
    cosine = as_strided(basis[indexsemi].m)
    sine = as_strided(basis[indexsemi+1].m)
    amp = np.sqrt(cosine**2+sine**2)
    phi = np.arctan2(sine,cosine)

    sigcosine = as_strided(basis[indexseas].sigmam)
    sigsine = as_strided(basis[indexseas+1].sigmam)
    sigamp = np.sqrt(sigcosine**2+sigsine**2)
    sigphi = (sigcosine*abs(sine)+sigsine*abs(cosine))/(sigcosine**2+sigsine**2)

    if geotiff is not None:
        ds = driver.Create('ampw2t_coeff.tif', jend-jbeg, iend-ibeg, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(amp)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

        ds = driver.Create('ampw2t_sigcoeff.tif', jend-jbeg, iend-ibeg, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(sigamp)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

    else:
        fid = open('ampw2t_coeff.r4', 'wb')
        amp.flatten().astype('float32').tofile(fid)
        fid.close()

        fid = open('ampw2t_sigcoeff.r4', 'wb')
        sigamp.flatten().astype('float32').tofile(fid)
        fid.close()

    if geotiff is not None:
        ds = driver.Create('phiw2t_coeff.tif', jend-jbeg, iend-ibeg, 1, gdal.GDT_Float32)
        band.WriteArray(phi)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

        ds = driver.Create('phiw2t_sigcoeff.tif', jend-jbeg, iend-ibeg, 1, gdal.GDT_Float32)
        band.WriteArray(sigphi)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band.FlushCache()
        del ds

    else:
        fid = open('phiw2t_coeff.r4', 'wb')
        phi.flatten().astype('float32').tofile(fid)
        fid.close()

        fid = open('phiw2t_sigcoeff.r4', 'wb')
        sigphi.flatten().astype('float32').tofile(fid)
        fid.close()


#######################################################
# Plot
#######################################################

# plot ref term
vmax = np.abs([np.nanpercentile(basis[0].m,98.),np.nanpercentile(basis[0].m,2.)]).max()
vmin = -vmax

nfigure +=1
fig=plt.figure(nfigure,figsize=(14,12))

ax = fig.add_subplot(1,M,1)
cax = ax.imshow(basis[0].m,cmap=cmap,vmax=vmax,vmin=vmin)
cbar = fig.colorbar(cax, orientation='vertical',shrink=0.2)
setp(ax.get_xticklabels(), visible=False)
setp(ax.get_yticklabels(), visible=False)

# plot linear term
vmax = np.abs([np.nanpercentile(basis[1].m,98.),np.nanpercentile(basis[1].m,2.)]).max()
vmin = -vmax

ax = fig.add_subplot(1,M,2)
cax = ax.imshow(basis[1].m,cmap=cmap,vmax=vmax,vmin=vmin)
ax.set_title(basis[1].reduction)
cbar = fig.colorbar(cax, orientation='vertical',shrink=0.2)
setp(ax.get_xticklabels(), visible=False)
setp(ax.get_yticklabels(), visible=False)

# plot others
for l in range(2,Mbasis):
    vmax = np.abs([np.nanpercentile(basis[l].m,98.),np.nanpercentile(basis[l].m,2.)]).max()
    vmin = -vmax

    ax = fig.add_subplot(1,M,l+1)
    cax = ax.imshow(basis[l].m,cmap=cmap,vmax=vmax,vmin=vmin)
    ax.set_title(basis[l].reduction)
    # add colorbar
    cbar = fig.colorbar(cax, orientation='vertical',shrink=0.2)
    setp(ax.get_xticklabels(), visible=False)
    setp(ax.get_yticklabels(), visible=False)

for l in xrange(Mker):
    vmax = np.abs([np.nanpercentile(kernels[l].m,98.),np.nanpercentile(kernels[l].m,2.)]).max()
    vmin = -vmax

    ax = fig.add_subplot(1,M,Mbasis+l+1)
    cax = ax.imshow(kernels[l].m,cmap=cmap,vmax=vmax,vmin=vmin)
    ax.set_title(kernels[l].reduction)
    setp(ax.get_xticklabels(), visible=False)
    setp(ax.get_yticklabels(), visible=False)
    cbar = fig.colorbar(cax, orientation='vertical',shrink=0.2)

plt.suptitle('Time series decomposition')

nfigure += 1
fig.tight_layout()
fig.savefig('inversion.eps', format='EPS',dpi=150)

if plot=='yes':
    plt.show()
