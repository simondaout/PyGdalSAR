#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT 
################################################################################

"""\
invers_disp_pixel.py
-------------
Temporal decomposition of the time series delays of selected pixels (used depl_cumule (BIP format) and images_retenues, output of invers_pixel). 

Usage: invers_disp_pixel.py --cols=<values> --ligns=<values> [--cube=<path>] [--list_images=<path>] [--windowsize=<value>] [--windowrefsize=<value>]  [--lectfile=<path>] [--aps=<path>] \
[--linear=<value>] [--coseismic=<value>] [--postseismic=<value>] [--seasonal=<yes/no>] [--seasonal_increase=<yes/no>] [--vector=<path>] [--info=<path>]\
[--semianual=<yes/no>] [--bianual=<yes/no>] [--degreeday=<values>] [--dem=<yes/no>] [--imref=<value>] [--cond=<value>] [--slowslip=<value>] [--ineq=<value>] \
[--name=<value>] [--rad2mm=<value>] [--plot=<yes/no>] [<iref>] [<jref>] [--bounds=<value>] [--dateslim=<values>] [--plot_dateslim=<values>] [--color=<value>] [--fillstyle=<value>] 

invers_disp_pixel.py -h | --help
Options:

-h --help               Show this screen
--ncols VALUE           Pixel column numbers (eg. 200,400,450) 
--nligns VALUE          Pixel lines numbers  (eg. 1200,1200,3000) 
--cube PATH             Path to displacement file [default: depl_cumul_flat]
--list_images PATH      Path to list images file made of 5 columns containing for each images 1) number 2) Doppler freq (not read) 3) date in YYYYMMDD format 4) numerical date 5) perpendicular baseline [default: images_retenues]
--windowsize VALUE      Number of pixels around the pixel defining the window [default: 0]
--windowrefsize VALUE      Number of pixels around the referenced pixel defining the window [default: windowsize]
--lectfile PATH         Path to the lect.in file (output of invers_pixel) [default: lect.in]
--aps PATH              Path to the APS file giving the error associated to each dates [default: No weigthing]
--linear PATH     Add a linear function to the inversion
--coseismic PATH        Add heaviside functions to the inversion .Indicate coseismic time (e.g 2004.,2006.)
--postseismic PATH      Add logarithmic transients to each coseismic step. Indicate characteristic time of the log function, must be a serie of values of the same lenght than coseismic (e.g 1.,1.). To not associate postseismic function to a given coseismic step, put None (e.g None,1.) 
--slowslip   VALUE      Add slow-slip function in the inversion (as defined by Larson et al., 2004). Indicate median and characteristic time of the events (e.g. 2004.,1,2006,0.5) [default: None] 
--seasonal PATH         If yes, add seasonal terms in the inversion
--seasonal_increase PATH         If yes, add seasonal terms function of time in the inversion 
--degreeday Values        Add degree-day model (Stefan model) in the inversion. Indicate thawing and freezing onsets (e.g. 0.36,0.7) [default: None] 
--semianual PATH        If yes, add semianual  terms in the inversion
--bianual PATH          If yes, add bianual  terms in the inversion
--vector PATH           Path to the vector text files containing a value for each dates [default: None]
--info PATH             Path to extra file in r4 or tif format to plot is value on the selected pixel, e.g. aspect [default: None].
--dem PATH              If yes, add term proportional to the perpendicular baseline in the inversion
--imref VALUE           Reference image number [default: 1]
--cond VALUE            Condition value for optimization: Singular value smaller than cond are considered zero [default: 1.e-3]
--ineq VALUE            If yes, add ineguality constrained in the inversion: use least square result to iterate the inversion. Force postseismic to be the  same sign than coseismic [default: no].       
--name Value            Name output figures [default: None] 
--rad2mm                Scaling value between input data (rad) and desired output [default: -4.4563]
--plot                  Display results [default: yes]            
iref                  colum numbers of the reference pixel [default: None] 
jref                  lign number of the reference pixel [default: None]
--bounds                yMin,yMax time series plots 
--dateslim              Datemin,Datemax time series  
--plot_dateslim         Datemin,Datemax time series for plot only 
--color                 Colors time series points [default:blue]
--fillstyle             Fill Style time series points. Can be: none,full,top,bottom,right,left [default: none]
"""

print()
print()
print('Please cite:')
print ('Daout, S., Doin, M. P., Peltzer, G., Socquet, A., & Lasserre, C. (2017). Large‐scale InSAR monitoring of permafrost freeze‐thaw cycles on the Tibetan Plateau. Geophysical Research Letters, 44(2), 901-909.')
print('Daout, S., Sudhaus, H., Kausch, T., Steinberg, A., & Dini, B. (2019). Interseismic and postseismic shallow creep of the North Qaidam Thrust faults detected with a multitemporal InSAR analysis. Journal of Geophysical Research: Solid Earth, 124(7), 7259-7279.')
print()
print()

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided
import numpy.linalg as lst

# scipy
import scipy
import scipy.optimize as opt

# basic
import math,sys,getopt
from os import path, environ
import os

# gdal
from osgeo import gdal

# plot
import matplotlib
if environ["TERM"].startswith("screen"):
    matplotlib.use('Agg') 
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from datetime import datetime as datetimes
import datetime
import time

# docopt (command line parser)
import docopt


########################################################################
# Define basis functions
########################################################################


class pattern:
    def __init__(self,name,reduction,date):
        self.name=name
        self.reduction=reduction
        self.date=date
    
    def info(self):
        print(self.name, self.date)

def Heaviside(t):
        h=np.zeros((len(t)))
        h[t>=0]=1.0
        return h

class coseismic(pattern):
      def __init__(self,name,reduction,date):
          pattern.__init__(self,name,reduction,date)
          self.to=date

      def g(self,t):
        return (Heaviside(t-self.to))

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

class slowslip(pattern):
      def __init__(self,name,reduction,date,tcar=1):
          pattern.__init__(self,name,reduction,date)
          self.to=date
          self.tcar=tcar

      def g(self,t):
          t=(t-self.to)/self.tcar
          funct = 0.5*(np.tanh(t)-1) + 1
          return funct

class stefan(pattern):
    def __init__(self,name,reduction,date,tcar1,tcar2):
          pattern.__init__(self,name,reduction,date)
          self.to = date
          self.t1 = tcar1
          self.t2 = tcar2

    def g(self,t):
        # t = np.arange(2014.,2015.,0.001)
        day = np.mod(t-self.to,1)
        # ddt: degree-day of thawing
        ddt = day - self.t1; ddt[day<=self.t1] = self.t2 - self.t1; ddt[day>=self.t2] = self.t2 - self.t1
        # ddf: degree-day of freezing
        ddf = np.zeros(len(t))
        ddf[day>self.t2] = day[day>self.t2] - self.t2
        ddf[day<self.t1] = day[day<self.t1] + 1 - self.t2

        # lets assune same coef. of freezing and thawing
        alpha = np.sqrt((self.t2-self.t1) / (1 - self.t2 + self.t1))
        func = np.sqrt(ddt) - alpha*np.sqrt(ddf) 
        # plt.plot(t,func,'-b')
        # plt.show()
        # sys.exit()
        return func

### KERNEL FUNCTIONS: not function of time
class corrdem(pattern):
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
        self.func = vect

    def info(self):
        print(self.name)

    def g(self,index):
        return self.func[index]

def date2dec(dates):
    dates  = np.atleast_1d(dates)
    times = []
    for date in dates:
        x = datetimes.strptime('{}'.format(date),'%Y%m%d')
        dec = float(x.strftime('%j'))/365.1
        year = float(x.strftime('%Y'))
        times.append(year + dec)
    return times


########################################################################

# read arguments
arguments = docopt.docopt(__doc__)

# initialise iteration with interseismic alone
iteration=False

if arguments["--list_images"] ==  None:
    listim = "images_retenues"
else:
    listim = arguments["--list_images"]

if arguments["--ineq"] ==  None:
    ineq = 'yes' 
else:
    ineq = arguments["--ineq"] 
if arguments["--cond"] ==  None:
    rcond = 1e-3
else:
    rcond = float(arguments["--cond"]) 
if arguments["--imref"] ==  None:
    imref = 0
else:
    imref = int(arguments["--imref"]) - 1
if arguments["--cube"] ==  None:
    cubef = "depl_cumule"
else:
    cubef = arguments["--cube"]
if arguments["--lectfile"] ==  None:
    infile = "lect.in"
else:
    infile = arguments["--lectfile"]
if arguments["--windowsize"] ==  None:
    w = 0
else:
    w = int(arguments["--windowsize"])
if arguments["--windowrefsize"] ==  None:
    wref = np.copy(w)
else:
    wref = int(arguments["--windowrefsize"])
if arguments["--aps"] ==  None:
    apsf = None
else:
    apsf = arguments["--aps"]
if arguments["--linear"] ==  None:
    inter = 'no'
else:
    inter = arguments["--linear"]
if arguments["--seasonal"] ==  None:
    seasonal = 'no'
else:
    seasonal = arguments["--seasonal"]
if arguments["--seasonal_increase"] ==  None:
    seasonalt = 'no'
else:
    seasonalt = arguments["--seasonal_increase"]
if arguments["--semianual"] ==  None:
    semianual = 'no'
else:
    semianual = arguments["--semianual"]
if arguments["--bianual"] ==  None:
    bianual = 'no'
else:
    bianual = arguments["--bianual"]

if arguments["--dem"] ==  None:
    dem = 'no'
else:
    dem = arguments["--dem"]
if arguments["--coseismic"] ==  None:
    cos = []
else:
    cos = list(map(float,arguments["--coseismic"].replace(',',' ').split()))
if arguments["--postseismic"] ==  None:
    pos = np.zeros(len(cos))
else:
    pos = list(map(float,arguments["--postseismic"].replace('None','-1').replace(',',' ').split()))

if arguments["--slowslip"] == None:
    sse=[]
else:
    sse = list(map(float,arguments["--slowslip"].replace(',',' ').split())) 
sse_time = sse[::2]
sse_car = sse[1::2]  

if arguments["--degreeday"] ==  None:
    degreeday = 'no'
else:
    degreeday = 'yes'
    try:
        ddt,ddf = list(map(float,arguments["--degreeday"].replace(',',' ').split())) 
    except:
        print('degreeday argument must contain two float values corresponding to the thawing and freezing onsets')
        sys.exit(0)

if arguments["--vector"] != None:
    vectf = arguments["--vector"].replace(',',' ').split()
else:
    vectf = None

if arguments["--bounds"] is not  None:
    ylim = list(map(float,arguments["--bounds"].replace(',',' ').split()))

if arguments["<iref>"] ==  None:
    iref = None
else:
    iref = int(arguments["<iref>"])
if arguments["<jref>"] ==  None:
    jref = None
else:
    jref = int(arguments["<jref>"])

if arguments["--name"] ==  None:
    output = None
else:
    output = arguments["--name"]

if arguments["--rad2mm"] ==  None:
    rad2mm = -4.4563
else:
    rad2mm = float(arguments["--rad2mm"]) 

if arguments["--plot"] ==  None:
    plot = 'yes'
else:
    plot = arguments["--plot"]

if arguments["--info"] ==  None:
   infof = None
else:
   infof = arguments["--info"]

if arguments["--color"] ==  None:
   color = "blue"
else:
   color = arguments["--color"]
if arguments["--fillstyle"] ==  None:
   fillstyle = "none"
else:
   fillstyle = arguments["--fillstyle"]

if len(pos)>0 and len(cos) != len(pos):
    raise Exception("coseimic and postseismic lists are not the same size")

# markers = ['o','v','^','s','P','X','o','v','^','s','P','X']
markers = ['o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o']

ipix = list(map(int,arguments["--cols"].replace(',',' ').split()))
jpix = list(map(int,arguments["--ligns"].replace(',',' ').split()))
if len(jpix) != len(ipix):
   raise Exception("ncols and nligns lists are not the same size")
# number of pixels
Npix = len(ipix)

# read lect.in 
ncol, nlign = list(map(int, open(infile).readline().split(None, 2)[0:2]))

# bounds plots
istart,iend = np.min(ipix) - 100, np.max(ipix) + 100
jstart,jend = np.min(jpix) - 100, np.max(jpix) + 100
if istart < 0:
    istart = 0
if jstart < 0:
    jstart = 0
if iend > ncol:
    iend = ncol
if jend > nlign:
    jend = nlign

# load images_retenues file
nb,idates,dates,base=np.loadtxt(listim, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
N = len(dates)
# save bp0 before time crop
bp0=base[imref] 

if arguments["--dateslim"] is not  None:
    dmin,dmax = arguments["--dateslim"].replace(',',' ').split()
    datemin = date2dec(dmin)
    datemax = date2dec(dmax)
else:
    datemin, datemax = int(np.min(dates)), int(np.max(dates))+1
    dmax = str(datemax) + '0101'
    dmin = str(datemin) + '0101'


# clean dates
indexd = np.flatnonzero(np.logical_and(dates<datemax,dates>datemin))
nb,idates,dates,base = nb[indexd],idates[indexd],dates[indexd],base[indexd]

# lect cube
cubei = np.fromfile(cubef,dtype=np.float32)
cube = as_strided(cubei[:nlign*ncol*N])
print('Number of line in the cube: ', cube.shape)
kk = np.flatnonzero(np.logical_or(cube==9990, cube==9999))
cube[kk] = float('NaN')
# ATTENTION: here i convert rad to mm
cube = cube*rad2mm
maps = cube.reshape((nlign,ncol,N))

# set to NaN lakes and ocean, new image ref.
cst = np.copy(maps[:,:,imref])
for l in range((N)):
    d = as_strided(maps[:,:,l])
    maps[:,:,l] = maps[:,:,l] - cst

# new number of dates
N = len(dates)
maps = as_strided(maps[:,:,indexd])
print('Reshape cube: ', maps.shape)

# perp baseline term
base_moy = np.mean(base)

if apsf is not None:
    inaps=np.loadtxt(apsf, comments='#', unpack=True,dtype='f')*abs(rad2mm)
    print('Input uncertainties:', inaps)
    print('Set very low values to the 2 percentile to avoid overweighting...')
    # maxinaps = np.nanmax(inaps) 
    # inaps= inaps/maxinaps

    minaps= np.nanpercentile(inaps,2)
    index = flatnonzero(inaps<minaps)
    inaps[index] = minaps
    try:
        inaps = inaps[indexd]
    except:
        pass
    print('Output uncertainties for first iteration:', inaps)

if vectf is not None:
    v = np.zeros((len(vectf),N))
    for i in range(len(vectf)):
        v[i,:] = np.loadtxt(vectf[i], comments='#', unpack = False, dtype='f')[indexd]


if infof is not None:
    extension = os.path.splitext(infof)[1]
    if extension == ".tif":
      ds = gdal.Open(infof, gdal.GA_ReadOnly)
      band = ds.GetRasterBand(1)
      info = band.ReadAsArray()
      del ds
    else:
      fid = open(infof,'r')
      infoi = np.fromfile(fid,dtype=np.float32)
      info = infoi[:nlign*ncol].reshape((nlign,ncol))
      fid.close()

# plot pixels on map
fig = plt.figure(1,figsize=(12,8))
# if arguments["--bounds"] is not  None:
#     vmax,vmin = np.nanmax(ylim), np.nanmin(ylim)
# else:
vmax = np.nanpercentile(maps[:,:,-1],90)
vmin = np.nanpercentile(maps[:,:,-1],10)

ax = fig.add_subplot(1,2,1)
ax.imshow(maps[jstart:jend,istart:iend,-1], cmap=cm.rainbow, vmax=vmax, vmin=vmin, alpha=0.6)
for i in range(len(ipix)):
    ax.scatter(ipix[i]-istart,jpix[i]-jstart,marker=markers[i],color='black',s=10.)
if iref is not None and jref is not None:
    ax.scatter(iref-istart,jref-jstart,marker='x',color='red',s=20.)
for i in range((Npix)):
    ax.text(ipix[i]-istart,jpix[i]-jstart,i)
plt.suptitle('Black cross: pixels, red cross: reference point')

ax = fig.add_subplot(1,2,2)
im = ax.imshow(maps[:,:,-1], cmap=cm.rainbow, vmax=vmax, vmin=vmin, alpha=0.6)
for i in range(len(ipix)):
    ax.scatter(ipix[i],jpix[i],marker=markers[i],color='black',s=10.)
ax.scatter(iref,jref,marker='x',color='red',s=20.)
for i in range((Npix)):
    ax.text(ipix[i],jpix[i],i)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.suptitle('Black symbols: pixels, red cross: reference point')
plt.savefig('Map_{}.pdf'.format(output), format='PDF')
# print(iref, jref)
# plt.show()
# sys.exit()

# 0
basis=[
      reference(name='reference',date=datemin,reduction='ref'),
      ]
index = len(basis)

if inter=='yes':
      basis.append(interseismic(name='interseismic',reduction='lin',date=datemin))
      # 1
      indexinter=index
      index = index + 1

if degreeday=='yes':
   indexdd = index
   print(ddt, ddf, datemin)
   basis.append(stefan(name='degree-day',reduction='ddt',date=datemin,tcar1=ddt,tcar2=ddf))
   index = index + 1

if seasonal=='yes':
   # 2
   indexseas = index
   basis.append(cosvar(name='seas. var (cos)',reduction='coswt',date=datemin))
   basis.append(sinvar(name='seas. var (sin)',reduction='sinwt',date=datemin))
   index = index + 2

if seasonalt=='yes':
   # 2
   indexseast = index
   basis.append(cost(name='cost',reduction='cost',date=datemin))
   basis.append(sint(name='sint',reduction='sint',date=datemin))
   index = index + 2

if semianual=='yes':
   # 4
   indexsemi = index 
   basis.append(cos2var(name='semi. annual var (cos)',reduction='cos2wt',date=datemin))
   basis.append(sin2var(name='semi. annual var (sin)',reduction='sin2wt',date=datemin))
   index = index + 2

if bianual=='yes':
   # 4
   indexbi = index 
   basis.append(cos5var(name='bi-annual var (cos)',reduction='cos5wt',date=datemin))
   basis.append(sin5var(name='bi-annual var (sin)',reduction='sin5wt',date=datemin))
   index = index + 2

indexco = np.zeros(len(cos))
for i in range(len(cos)):
   # 6
   indexco[i] = int(index)
   basis.append(coseismic(name='coseismic {}'.format(i),reduction='cos{}'.format(i),date=cos[i])),
   index = index + 1
   iteration=True

indexsse = np.zeros(len(sse_time))
for i in range(len(sse_time)):
    basis.append(slowslip(name='sse {}'.format(i),reduction='sse{}'.format(i),date=sse_time[i],tcar=sse_car[i])),
    indexsse[i] = int(index)
    index = index + 1
    iteration=True

indexpo,indexpofull = [],[]
for i in range(len(pos)):
  if pos[i] > 0. :
    basis.append(postseismic(name='postseismic {}'.format(i),reduction='post{}'.format(i),date=cos[i],tcar=pos[i])),
    indexpo.append(int(index))
    indexpofull.append(int(index))
    index = index + 1
  else:
    indexpofull.append(0)
indexpo = np.array(indexpo)
indexpofull = np.array(indexpofull)

kernels=[]

if dem=='yes':
   kernels.append(corrdem(name='dem correction',reduction='corrdem',bp0=bp0,bp=base))
   indexdem = index
   index = index + 1

if vectf != None:
  indexvect = np.zeros(len(vectf))
  fig = plt.figure(100,figsize=(12,8))
  for i in range(len(vectf)):
    kernels.append(vector(name=vectf[i],reduction='vector_{}'.format(i),vect=v[i]))
    ax = fig.add_subplot(len(vectf),1,i+1)
    ax.plot(dates,kernels[-1].g(~np.isnan(dates)))
    indexvect[i] = index
    index = index + 1

indexpo = indexpo.astype(int)
indexco = indexco.astype(int)
indexsse = indexsse.astype(int)

# define size G matrix
print()
Mbasis=len(basis)
print('Number of basis functions:', Mbasis)
Mker=len(kernels)
print('Number of kernel functions:', Mker)
M = Mbasis + Mker
for i in range((Mbasis)):
    basis[i].info()

# SVD inversion with cut-off eigenvalues
def invSVD(A,b,cond):
    try:
        U,eignv,V = lst.svd(A, full_matrices=False)
        s = np.diag(eignv)
        print(s)
        index = np.nonzero(s<cond)
        inv = lst.inv(s)
        inv[index] = 0.
        fsoln = np.dot( V.T, np.dot( inv , np.dot(U.T, b) ))
    except:
        fsoln = lst.lstsq(A,b)[0]
        #fsoln = lst.lstsq(A,b,rcond=cond)[0]
    
    return fsoln

## inversion procedure 
#def consInvert(A,b,sigmad,ineq='yes',cond=1.0e-3, iter=2000,acc=1e-12, eguality=False):
def consInvert(A,b,sigmad,ineq='yes',cond=1.0e-3, iter=100,acc=1e-4, eguality=False):
    '''Solves the constrained inversion problem.

    Minimize:
    
    ||Ax-b||^2

    Subject to:
    mmin < m < mmax
    '''

    if A.shape[0] != len(b):
        raise ValueError('Incompatible dimensions for A and b')

    if ineq == 'no':
        print('ineq=no: SVD decomposition neglecting small eigenvectors inferior to {} (cond)'.format(cond))
        fsoln = invSVD(A,b,cond)
        print('SVD solution:', fsoln)

    else:
        print('ineq=yes: Iterative least-square decomposition. Prior obtained with SVD.')
        if len(indexpo>0):
          # invert first without post-seismic
          Ain = np.delete(A,indexpo,1)
          try:
              U,eignv,V = lst.svd(Ain, full_matrices=False)
              s = np.diag(eignv) 
              print('Eigenvalues:', eignv)
              index = np.nonzero(s<cond)
              inv = lst.inv(s)
              inv[index] = 0.
              mtemp = np.dot( V.T, np.dot( inv , np.dot(U.T, b) ))
          except:
              mtemp = lst.lstsq(Ain,b,rcond=cond)[0]
          print('SVD solution:', mtemp)

          # rebuild full vector
          for z in range(len(indexpo)):
            mtemp = np.insert(mtemp,indexpo[z],0)
          minit = np.copy(mtemp)
          # # initialize bounds
          mmin,mmax = -np.ones(len(minit))*np.inf, np.ones(len(minit))*np.inf 

          # We here define bounds for postseismic to be the same sign than coseismic
          # and coseismic inferior or egual to the coseimic initial 
          print('ineq=yes: Impose postseismic to be the same sign than coseismic')
          for i in range(len(indexco)):
            if (pos[i] > 0.) and (minit[int(indexco[i])]>0.):
                mmin[int(indexpofull[i])], mmax[int(indexpofull[i])] = 0, np.inf 
                mmin[int(indexco[i])], mmax[int(indexco[i])] = 0, minit[int(indexco[i])] 
            if (pos[i] > 0.) and (minit[int(indexco[i])]<0.):
                mmin[int(indexpofull[i])], mmax[int(indexpofull[i])] = -np.inf , 0
                mmin[int(indexco[i])], mmax[int(indexco[i])] = minit[int(indexco[i])], 0
          bounds=list(zip(mmin,mmax))

        else:
          minit=invSVD(A,b,cond)
          print('SVD solution:', minit)
          bounds=None
        
        def eq_cond(x, *args):
           return math.atan2(x[indexseast+1],x[[indexseast]]) - math.atan2(x[indexseas+1],x[[indexseas]])
       
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
        print('Optimization:', fsoln)

    # tarantola:
    # Cm = (Gt.Cov.G)-1 --> si sigma=1 problems
    # sigma m **2 =  misfit**2 * diag([G.TG]-1)
    try:
       varx = np.linalg.inv(np.dot(A.T,A))
       # res2 = np.sum(pow((b-np.dot(A,fsoln))/sigmad,2))
       res2 = np.sum(pow((b-np.dot(A,fsoln)),2))
       scale = 1./(A.shape[0]-A.shape[1])
       # scale = 1./A.shape[0]
       sigmam = np.sqrt(scale*res2*np.diag(varx))
    except:
       sigmam = np.ones((A.shape[1]))*float('NaN')
    print('model errors:', sigmam)

    return fsoln,sigmam

# plot diplacements maps
nfigure = 10
if Npix > 2:
    fig = plt.figure(nfigure+1,figsize=(12,8))
else:
    fig = plt.figure(nfigure+1,figsize=(10,4))
if inter=='yes':  
    if Npix > 2:  
        fig3 = plt.figure(nfigure+3,figsize=(12,8))
    else:
        fig3 = plt.figure(nfigure+3,figsize=(10,4))
if seasonal == 'yes' or semianual == 'yes' or bianual == 'yes' or seasonalt == 'yes':
    if Npix > 2:
        fig2 = plt.figure(nfigure+2,figsize=(12,8))
    else:
        fig2 = plt.figure(nfigure+2,figsize=(10,4))


for jj in range((Npix)):
    # initialize aps before iterations   
    if apsf is None:
        aps = np.ones((N))
    else:
        aps = inaps

    # for ll in range(niter):
    i, j = ipix[jj], jpix[jj]
    print()
    print('---------------------------------------')
    print('pixel:{} {}'.format(i,j))
    print('---------------------------------------')
    print() 

    x = [date2num(datetimes.strptime('{}'.format(d),'%Y%m%d')) for d in idates]
    if arguments["--dateslim"] is not  None:
        dmin,dmax = arguments["--dateslim"].replace(',',' ').split()
    elif arguments["--plot_dateslim"] is not  None:
        dmin,dmax = arguments["--plot_dateslim"].replace(',',' ').split()
    else:
        dmax = str(datemax) + '0101'
        dmin = str(datemin) + '0101'
    # dmin,dmax = 20030101,20110101
    # dmin,dmax = 20140101,20200101
    xmin = datetimes.strptime('{}'.format(dmin),'%Y%m%d')
    xmax = datetimes.strptime('{}'.format(dmax),'%Y%m%d')
    xlim=date2num(np.array([xmin,xmax]))

    # extract data
    wind = as_strided(maps[j-w:j+w+1,i-w:i+w+1,:])
    if infof is not None:
      infm = np.nanmedian(info[j-w:j+w+1,i-w:i+w+1])
    if iref is not None:
        windref = as_strided(maps[jref-wref:jref+wref+1,iref-wref:iref+wref+1,:])
        dispref = np.nanmedian(windref,axis=(0,1))
    else:
        dispref = np.zeros((N))

    disp = np.nanmedian(wind,axis=(0,1)) - dispref
    #aps = np.nanstd(wind,axis=(0,1))

    # inversion model
    mdisp= np.ones((N))*float('NaN')
    demerr = np.zeros((N))
    lin = np.zeros((N))
    k = np.flatnonzero(~np.isnan(disp))
    kk = len(k)
    tabx = dates[k]
    taby = disp[k] 
    bp = base[k]
    sigmad = aps  
    names = []
    # print('data uncertainties', sigmad)    

    # Inisilize 
    G=np.zeros((kk,M))
    m,sigmam = np.zeros((M)),np.zeros((M))
    # do only this if more than N/2 points left
    if kk > N/6:
        # inversion
        t = time.time()

        names = []
        for l in range((Mbasis)):
            G[:,l]=basis[l].g(tabx)
            names.append(basis[l].reduction)
        for l in range((Mker)):
            G[:,Mbasis+l]=kernels[l].g(k)
            names.append(kernels[l].reduction)

        print('basis functions:', names)
        eguality = False
        if seasonalt == 'yes' and seasonal == 'yes':
            eguality = True
        mt,sigmamt = consInvert(G,taby,sigmad[k],cond=rcond, ineq=ineq, eguality=eguality)

        # rebuild full vectors
        if Mker>0:
            m[Mbasis:],sigmam[Mbasis:] = mt[-Mker:],sigmamt[-Mker:]
            m[:mt.shape[0]-Mker],sigmam[:mt.shape[0]-Mker] = mt[:-Mker],sigmamt[:-Mker]
        else:
            m [:mt.shape[0]] = mt  

        # forward model in original order
        mdisp[k] = np.dot(G,m)
        if dem=='yes':
            demerr[k] =  np.dot(G[:,indexdem],m[indexdem])
        else:
            demerr = np.zeros((N))
        
        print()
        print('computation time:', time.time() - t)
        print()
        

    # plot at the end of iterations for each points
    ax = fig.add_subplot(Npix,1,jj+1)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
    
    if inter=='yes':
        ax3 = fig3.add_subplot(Npix,1,jj+1)
        ax3.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
    if seasonal == 'yes' or semianual == 'yes' or bianual == 'yes' or seasonalt == 'yes':
        ax2 = fig2.add_subplot(Npix,1,jj+1)
        ax2.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
    
    # initialize variables
    disp_seas = np.zeros((N))
    disp_seast = np.zeros((N))

    if inter=='yes':
        lin[k] = np.dot(G[:,indexinter],m[indexinter])

    # plot data and model minus dem error and seasonal terms
    if seasonal=='yes':
            G=np.zeros((kk,2))
            for l in range((2)):
                G[:,l]=basis[l+indexseas].g(tabx)
            disp_seas[k] = disp_seas[k] + np.dot(G[:,:],m[indexseas:indexseas+2])
            amp,phi = np.sqrt(m[indexseas]**2+m[indexseas+1]**2),np.arctan2(m[indexseas+1],m[indexseas])

            if phi<0: 
               phi = phi + 2*np.pi
            #print(phi)

    if seasonalt=='yes':
            G=np.zeros((kk,2))
            for l in range((2)):
                G[:,l]=basis[l+indexseast].g(tabx)
            disp_seas[k] = disp_seas[k] + np.dot(G[:,:],m[indexseast:indexseast+2])
            disp_seast[k] = disp_seast[k] + np.dot(G[:,:],m[indexseast:indexseast+2])

            amp_min = 0
            amp_max = np.sqrt(m[indexseast]**2+m[indexseast+1]**2)*(datemax-datemin)

            phi = np.arctan2(m[indexseast+1],m[indexseast])
            
            # convert between 0 and 2pi
            if phi<0: 
                phi = phi + 2*np.pi
            #print(phi)

            a_rate = (amp_max - amp_min)/(datemax-datemin)

            if phi<0: 
               phi = phi + 2*np.pi

    if seasonalt=='yes' and seasonal=='yes':
        a1_min, a2_min = np.sqrt(m[indexseas]**2+m[indexseas+1]**2),0
        a1_max, a2_max = np.sqrt(m[indexseas]**2+m[indexseas+1]**2),np.sqrt(m[indexseast]**2+m[indexseast+1]**2)*(datemax-datemin)

        phi1, phi2 = np.arctan2(m[indexseas+1],m[indexseas]), np.arctan2(m[indexseast+1],m[indexseast])
        
        if phi1<0: 
            phi1 = phi1 + 2*np.pi

        if phi2<0: 
            phi2 = phi2 + 2*np.pi

        amp_min = a1_min
        phi_min = phi1
        # convert between 0 and 2pi
        if phi_min<0: 
            phi_min = phi_min + 2*np.pi

        amp_max = np.sqrt( (a1_max*np.cos(phi1) + a2_max*np.cos(phi2))**2 +
            (a1_max*np.sin(phi1) + a2_max*np.sin(phi2))**2)

        phi_max = np.arctan2(a1_max*np.sin(phi1)+a2_max*np.sin(phi2),a1_max*np.cos(phi1)+a2_max*np.cos(phi2)) 

        if phi_max<0: 
            phi_max = phi_max + 2*np.pi

        a_rate = (amp_max - amp_min)/(datemax-datemin)
        
    if semianual=='yes':
            G=np.zeros((kk,2))
            for l in range((2)):
                G[:,l]=basis[l+indexsemi].g(tabx)
            disp_seas[k] = disp_seas[k] +  np.dot(G[:,:],m[indexsemi:indexsemi+2])

    if bianual=='yes':
            G=np.zeros((kk,2))
            for l in range((2)):
                G[:,l]=basis[l+indexbi].g(tabx)
            disp_seas[k] = disp_seas[k] +  np.dot(G[:,:],m[indexbi:indexbi+2])

    # plot data and model minus dem error
    if infof is not None:
      # print(infof, infm)
      ax.plot(x,disp-demerr,markers[jj],color=color,fillstyle=fillstyle,label='TS {}: lign:{}, column:{}, Info:{:.2f}'.format(jj,i,j,infm))
    else:
      ax.plot(x,disp-demerr,markers[jj],color=color,fillstyle=fillstyle,label='TS {}: lign:{}, column:{}'.format(jj,i,i,j,j))
    
    ax.errorbar(x,disp-demerr,yerr = sigmad, ecolor=color,fmt='none', alpha=0.5)
    
    # plot data and model minus dem error and linear term
    if inter=='yes':
        ax3.plot(x,disp-demerr-lin,markers[jj],color=color,fillstyle=fillstyle,markersize=4,label='detrended data')
        ax3.errorbar(x,disp-demerr-lin,yerr = sigmad, ecolor=color,fmt='none', alpha=0.5)
            
    if semianual=='yes' or seasonal=='yes' or bianual=='yes' or seasonalt == 'yes':
        ax2.plot(x,disp-disp_seas-demerr,markers[jj],color=color,fillstyle=fillstyle,label='data -seasonal')
        # ax2.plot(x,mdisp-disp_seas-demerr,'o',color='red',alpha=0.5,label='model -seasonal')
        ax2.errorbar(x,disp-disp_seas-demerr,yerr = sigmad, ecolor=color,fmt='none', alpha=0.3)

    # create synthetic time
    t = np.array([xmin + datetime.timedelta(days=d) for d in range(0, 2920)])
    tdec = np.array([float(date.strftime('%Y')) + float(date.strftime('%j'))/365.1 for date in t])
    mseas = np.zeros(len(tdec))
    mseast = np.zeros(len(tdec))

    G=np.zeros((len(tdec),M))
    for l in range((Mbasis)):
        G[:,l]=basis[l].g(tdec)
    for l in range((Mker)):
        G[:,Mbasis+l]=np.interp(tdec,tabx,kernels[l].g(k))
    model = np.dot(G,m)
    
    if inter=='yes':
        model_lin = np.dot(G[:,indexinter],m[indexinter])
    if dem=='yes':
        model_dem = np.dot(G[:,indexdem],m[indexdem])
    else:
        model_dem = np.zeros(len(model))
        
    if seasonal=='yes':
        G=np.zeros((len(tdec),2))
        for l in range((2)):
            G[:,l]=basis[l+indexseas].g(tdec)
        mseas = mseas + np.dot(G[:,:],m[indexseas:indexseas+2])
        
    if semianual=='yes':
        G=np.zeros((len(tdec),2))
        for l in range((2)):
            G[:,l]=basis[l+indexsemi].g(tdec)
        mseas = mseas + np.dot(G[:,:],m[indexsemi:indexsemi+2])

    if seasonalt=='yes':
        G=np.zeros((len(tdec),2))
        for l in range((2)):
            G[:,l]=basis[l+indexseast].g(tdec)
        mseas = mseas + np.dot(G[:,:],m[indexseast:indexseast+2])
        # mseast = mseast + np.dot(G[:,:],m[indexseast:indexseast+2])

    if bianual=='yes':
        G=np.zeros((len(tdec),2))
        for l in range((2)):
            G[:,l]=basis[l+indexbi].g(tdec)
        mseas = mseas + np.dot(G[:,:],m[indexbi:indexbi+2])

    # plot model
    if inter=='yes':
        if seasonal=='yes' and seasonalt=='no':
            ax.plot(t,model-model_dem,'-r',label='Rate: {:.2f}, Amp: {:.2f}, Phi: {:.2f}'.format(m[indexinter],amp,phi))
            ax3.plot(t,model-model_lin-model_dem,'-r')
        elif seasonal=='yes' and seasonalt=='yes':
            ax.plot(t,model-model_dem,'-r',label='Rate: {:.2f}, Amp_beg: {:.2f}, Phi1: {:.2f}, Amp_end: {:.2f}, Phi2: {:.2f}, Rate_Amp: {:.2f}'.format(m[indexinter],amp_min,phi1,amp_max,phi2,a_rate))
            ax3.plot(t,model-model_lin-model_dem,'-r')
        elif seasonal=='no' and seasonalt=='yes':
            ax.plot(t,model-model_dem,'-r',label='Rate: {:.2f}, Amp_beg: {:.2f}, Amp_end: {:.2f}, Phi: {:.2f}, Rate_Amp: {:.2f}'.format(m[indexinter],amp_min,amp_max,phi,a_rate))
            ax3.plot(t,model-model_lin-model_dem,'-r')
        else:
            ax.plot(t,model-model_dem,'-r',label='Rate: {:.2f}'.format(m[indexinter]))
            ax3.plot(t,model-model_lin-model_dem,'-r')
    else:
        ax.plot(t,model-model_dem,'-r')
           
    if seasonal=='yes' or semianual=='yes' or bianual=='yes' or seasonalt=='yes':
        ax2.plot(t,model-mseas-model_dem,'-r')
        # ax2.plot(t,model-mseast-model_dem,'-r')
        ax2.legend(loc='best',fontsize='x-small')
        ax2.set_xlim(xlim)
        # if arguments["--bounds"] is not  None:
        #     ax2.set_ylim(ylim)

    ax.legend(loc='best',fontsize='x-small')
    ax.set_xlim(xlim)
    # ax.set_xlim([2003,2011])
    if arguments["--bounds"] is not  None:
        ax.set_ylim(ylim)
    
    if inter=='yes':
        ax3.legend(loc='best',fontsize='x-small')
        ax3.set_xlim(xlim)
    #     if arguments["--bounds"] is not  None:
    #         ax3.set_ylim(ylim)

    fig.autofmt_xdate()
    ax.set_xlabel('Time (Year/month/day)')
    ax.set_ylabel('Displacements (mm)')
    fig.savefig('Disp_{}.pdf'.format(output), format='PDF')

    if inter=='yes':    
        fig3.autofmt_xdate()
        ax3.set_xlabel('Time (Year/month/day)')
        ax3.set_ylabel('Displacements (mm)')
        fig3.savefig('Disp_{}_detrended.pdf'.format(output), format='PDF')

    if seasonal == 'yes' or semianual =='yes' or bianual =='yes':
        fig2.autofmt_xdate()
        ax2.set_xlabel('Time (Year/month/day)')
        ax2.set_ylabel('Displacements (mm)')
        fig2.savefig('Disp_{}_deseas.pdf'.format(output), format='PDF')

    if seasonal == 'yes':
        cosine = as_strided(m[indexseas])
        sine = as_strided(m[indexseas+1])
        amp = np.sqrt(cosine**2+sine**2)
        phi = np.arctan2(sine,cosine)
        # print(amp,phi)

if plot == 'yes':
    plt.show()
sys.exit()
