#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT 
################################################################################


"""\
invers_vectors.py
-------------
Linear decomposition of temporal vectors

Usage: invers_vectors.py --vectors=<path> [--list_images=<path>] [--aps=<path>] [--interseismic=<value>] [--coseismic=<value>] [--postseismic=<value>] [--seasonal=<yes/no>] [--semianual=<yes/no>] [--slowslip=<value>] [--dem=<yes/no>] [--bounds=<value>] [--dateslim=<values>]

invers_vectors.py -h | --help
Options:
-h --help               Show this screen
--vectors PATH           Path to the vector text files containing a value for each dates [default: None]
--list_images PATH      Path to list images file made of 5 columns containing for each images 1) number 2) Doppler freq (not read) 3) date in YYYYMMDD format 4) numerical date 5) perpendicular baseline [default: images_retenues]
--lectfile PATH         Path to the lect.in file (output of invers_pixel) [default: lect.in]
--aps PATH              Path to the APS file giving the error associated to each dates [default: No weigthing]
--interseismic PATH     Add a linear function to the inversion
--coseismic PATH        Add heaviside functions to the inversion. Indicate coseismic time (e.g 2004.,2006.)
--postseismic PATH      Add logarithmic transients to each coseismic step, indicate characteristic time of the log function, must be a serie of values of the same lenght than coseismic (e.g 1.,1.). To not associate postseismic function to a give coseismic step, put None (e.g None,1.)
--slowslip   VALUE      Add slow-slip function in the inversion (as defined by Larson et al., 2004). Indicate median and characteristic time of the events (e.g. 2004.,1,2006,0.5), default: None
--seasonal PATH         If yes, add seasonal terms in the inversion
--semianual PATH        If yes, add semianual  terms in the inversion
--dem PATH              If yes, add term proportional to the perpendicular baseline in the inversion
--rad2mm                Scaling value between input data (rad) and desired output [default: -4.4563]
--bounds                Min,Max time series plots
--datelim               datemin,datemax time series plots
"""

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
from datetime import datetime
import datetime
import time
# docopt (command line parser)
import docopt


################################
# Create lib of basis functions
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
          self.tcar = tcar

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
        print(self.name)

    def g(self,index):
        func = (self.bp-self.bpo)
        return func[index]

################################
# Initialization
################################

# read arguments
arguments = docopt.docopt(__doc__)
vectf = arguments["--vectors"].replace(',',' ').split()

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

if arguments["--seasonal"] ==  None:
    seasonal = 'no'
else:
    seasonal = arguments["--seasonal"]
if arguments["--semianual"] ==  None:
    semianual = 'no'
else:
    semianual = arguments["--semianual"]
if arguments["--bounds"] is not  None:
    ylim = map(float,arguments["--bounds"].replace(',',' ').split())

Nv = len(vectf)
# load images_retenues file
nb,idates,dates,base=np.loadtxt(listim, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
datemin, datemax = np.int(np.min(dates)), np.int(np.max(dates))+1
N = len(dates)
# perp baseline term
base_moy = np.mean(base)

if apsf is not 'no':
    inaps=np.loadtxt(apsf, comments='#', unpack=True,dtype='f')*rad2mm
    print('Input uncertainties:', inaps)
    print('Set very low values to the 2 percentile to avoid overweighting...')
    # maxinaps = np.nanmax(inaps) 
    # inaps= inaps/maxinaps

    minaps= np.nanpercentile(inaps,2)
    index = flatnonzero(inaps<minaps)
    inaps[index] = minaps
    print('Output uncertainties for first iteration:', inaps)
else: 
    inaps = np.ones((N))

basis=[
      reference(name='reference',date=datemin,reduction='ref'),
      ]
index = len(basis)

if inter=='yes':
      basis.append(interseismic(name='interseismic',reduction='lin',date=datemin))
      # 1
      indexinter=index
      index = index + 1
if seasonal=='yes':
   # 2
   indexseas = index
   basis.append(cosvar(name='seas. var (cos)',reduction='coswt',date=datemin))
   basis.append(sinvar(name='seas. var (sin)',reduction='sinwt',date=datemin))
   index = index + 2

if semianual=='yes':
   # 4
   indexsemi = index
   basis.append(cos2var(name='semi. annual var (cos)',reduction='cos2wt',date=datemin))
   basis.append(sin2var(name='semi. annual var (sin)',reduction='sin2wt',date=datemin))
   index = index + 2


indexco = np.zeros(len(cos))
for i in range(len(cos)):
   # 6
   indexco[i] = int(index)
   basis.append(coseismic(name='coseismic {}'.format(i),reduction='cos{}'.format(i),date=cos[i])),
   index = index + 1

indexsse = np.zeros(len(sse_time))
for i in range(len(sse_time)):
    basis.append(slowslip(name='sse {}'.format(i),reduction='sse{}'.format(i),date=sse_time[i],tcar=sse_car[i])),
    indexsse[i] = int(index)
    index = index + 1

indexpo = []
for i in range(len(pos)):
   if pos[i] > 0. :
      basis.append(postseismic(name='postseismic {}'.format(i),reduction='post{}'.format(i),date=cos[i],tcar=pos[i])),
      indexpo.append(int(index))
      index = index + 1
indexpo = np.array(indexpo)

kernels=[]

if dem=='yes':
   kernels.append(corrdem(name='dem correction',reduction='corrdem',bp0=base[0],bp=base))
   indexdem = index
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
def invSVD(A,b,cond=0.1):
    try:
        U,eignv,V = lst.svd(A, full_matrices=False)
        s = np.diag(eignv)
        index = np.nonzero(s<cond)
        inv = lst.inv(s)
        inv[index] = 0.
        fsoln = np.dot( V.T, np.dot( inv , np.dot(U.T, b) ))
    except:
        fsoln = lst.lstsq(A,b,rcond=cond)[0]

    return fsoln

nfigure = 0
fig = plt.figure(nfigure,figsize=(10,4))
nfigure =+ 1 
for jj in range((Nv)):
    if apsf is None:
        aps = np.ones((N))
    else:
        aps = inaps

    x = [date2num(datetime.datetime.strptime('{}'.format(d),'%Y%m%d')) for d in idates]
    if arguments["--dateslim"] is not  None:
        dmin,dmax = arguments["--dateslim"].replace(',',' ').split()
    else:
        dmax = str(datemax) + '0101'
        dmin = str(datemin) + '0101'
    xmin = datetime.datetime.strptime('{}'.format(dmin),'%Y%m%d')
    xmax = datetime.datetime.strptime('{}'.format(dmax),'%Y%m%d')
    xlim=date2num(np.array([xmin,xmax]))

    basename = os.path.splitext(vectf[jj])[0]
    disp =  np.loadtxt(vectf[jj], comments='#', unpack = False, dtype='f')
    
    mdisp= np.ones((N))*float('NaN')
    demerr = np.zeros((N))
    
    k = np.flatnonzero(~np.isnan(disp))
    kk = len(k)
    tabx = dates[k]
    taby = disp[k]
    bp = base[k]
    sigmad = aps[k]
    print('data uncertainties', sigmad)
    
    m = np.zeros((M))
    G=np.zeros((kk,M))
    for l in range((Mbasis)):
        G[:,l]=basis[l].g(tabx)
    for l in range((Mker)):
        G[:,Mbasis+l]=kernels[l].g(k)

    m = invSVD(G,taby)
    mdisp[k] = np.dot(G,m)
    rmsd = np.sqrt(np.sum(pow((disp[k] - mdisp[k]),2))/kk)
    print()
    print('rmsd:', rmsd)
    print()

    if dem=='yes':
        demerr[k] =  np.dot(G[:,indexdem],m[indexdem])
    else:
        demerr = np.zeros((N))
    
    ax = fig.add_subplot(Nv,1,jj+1)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))

    ax.plot(x,disp-demerr,'o',label='TS:{}'.format(basename))
    if apsf is not 'no':
        ax.errorbar(x,disp-demerr,yerr = sigmad, ecolor='blue',fmt='none', alpha=0.5)

    # create synthetic time
    t = np.array([xmin + datetime.timedelta(days=d) for d in range(0, 2920)])
    tdec = np.array([float(date.strftime('%Y')) + float(date.strftime('%j'))/365.1 for date in t])

    G=np.zeros((len(tdec),M))
    for l in range((Mbasis)):
        G[:,l]=basis[l].g(tdec)
    for l in range((Mker)):
        G[:,Mbasis+l]=np.interp(tdec,tabx,kernels[l].g(k))
    model = np.dot(G,m)
    if dem=='yes':
        model_dem = np.dot(G[:,indexdem],m[indexdem])
    else:
        model_dem = np.zeros(len(model))
    
    # plot model
    ax.plot(t,model-model_dem,'-r')

    ax.legend(loc='best',fontsize='x-small')
    ax.set_xlim(xlim)
    if arguments["--bounds"] is not  None:
        ax.set_ylim(ylim)

    fig.autofmt_xdate()
    ax.set_xlabel('Time (Year/month/day)')
    ax.set_ylabel('Displacements (mm)')
    fig.savefig('{}.pdf'.format(basename), format='PDF')

plt.show()



















