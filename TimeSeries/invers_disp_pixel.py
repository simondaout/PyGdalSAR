#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT (CRPG-ENSG) 
################################################################################

"""\
invers_disp_pixel.py
-------------
Temporal decomposition of the time series delays of selected pixels (used depl_cumule (BIP format) and images_retenues, output of invers_pixel). 

Usage: invers_disp_pixel.py --cols=<values> --lines=<values> [--cube=<path>] [--list_images=<path>] [--windowsize=<value>] [--windowrefsize=<value>]  [--lectfile=<path>] [--aps=<path>] \
[--linear=<value>] [--steps=<value>] [--postseismic=<value>] [--seasonal=<yes/no>] [--seasonal_increase=<yes/no>] [--vector=<path>] [--info=<path>]\
[--semianual=<yes/no>] [--bianual=<yes/no>] [--degreeday=<values>] [--bperp=<yes/no>] [--imref=<value>] [--cond=<value>] [--slowslip=<value>] [--ineq=<value>] \
[--name=<value>] [--scale=<value>] [--plot=<yes/no>] [<iref>] [<jref>] [--bounds=<value>] [--dateslim=<values>] [--plot_dateslim=<values>] [--color=<value>] [--fillstyle=<value>] [--ndatasets=<nb_data_sets>] 

invers_disp_pixel.py -h | --help
Options:

-h --help               Show this screen
--ncols VALUE           Pixel column numbers (eg. 200,400,450) 
--nlines VALUE          Pixel lines numbers  (eg. 1200,1200,3000) 
--cube PATH             Path to displacement file [default: depl_cumul_flat]
--list_images PATH      Path to list images file made of 5 columns containing for each images 1) number 2) Doppler freq (not read) 3) date in YYYYMMDD format 4) numerical date 5) perpendicular baseline [default: images_retenues]
--windowsize VALUE      Number of pixels around the pixel defining the window [default: 0]
--windowrefsize VALUE      Number of pixels around the referenced pixel defining the window [default: windowsize]
--lectfile PATH         Path to the lect.in file (output of invers_pixel) [default: lect.in]
--aps PATH              Path to the APS file giving the error associated to each dates [default: No weigthing]
--linear PATH     Add a linear function to the inversion
--steps PATH        Add heaviside functions to the inversion .Indicate steps time (e.g 2004.,2006.)
--postseismic PATH      Add logarithmic transients to each steps step. Indicate characteristic time of the log function, must be a serie of values of the same lenght than steps (e.g 1.,1.). To not associate postseismic function to a given steps step, put None (e.g None,1.) 
--slowslip   VALUE      Add slow-slip function in the inversion (as defined by Larson et al., 2004). Indicate median and characteristic time of the events (e.g. 2004.,1,2006,0.5) [default: None] 
--seasonal PATH         If yes, add seasonal terms in the inversion
--seasonal_increase PATH         If yes, add seasonal terms function of time in the inversion 
--degreeday Values        Add degree-day model (Stefan model) in the inversion. Indicate thawing and freezing onsets (e.g. 0.36,0.7) [default: None] 
--semianual PATH        If yes, add semianual  terms in the inversion
--bianual PATH          If yes, add bianual  terms in the inversion
--vector PATH           Path to the vector text files containing a value for each dates [default: None]
--info PATH             Path to extra file in r4 or tif format to plot is value on the selected pixel, e.g. aspect [default: None].
--bperp PATH              If yes, add term proportional to the perpendicular baseline in the inversion
--imref VALUE           Reference image number [default: 1]
--cond VALUE            Condition value for optimization: Singular value smaller than cond are considered zero [default: 1.e-6]
--ineq VALUE            If yes, add inequality constrained in the inversion: use least square result to iterate the inversion. Force postseismic to be the  same sign than steps [default: no].       
--name Value            Name output figures [default: None] 
--scale                Scaling value between input data and desired output [default: 1]
--plot                  Display results [default: yes]            
iref                  colum numbers of the reference pixel [default: None] 
jref                  lign number of the reference pixel [default: None]
--bounds                yMin,yMax time series plots 
--dateslim              Datemin,Datemax time series  
--plot_dateslim         Datemin,Datemax time series for plot only 
--color                 Colors time series points [default:blue]
--fillstyle             Fill Style time series points. Can be: none,full,top,bottom,right,left [default: none]
--ndatasets=<nb_data_sets> Number of data sets [Default:1]
"""

print()
print()
print('Author: Simon Daout')
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
from datetime import datetime as dt
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

class steps(pattern):
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
        self.func = vect

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


########################################################################

# read arguments
arguments = docopt.docopt(__doc__)

if arguments["--list_images"] ==  None:
    arguments["--list_images"] = "images_retenues"
if arguments["--ineq"] ==  None:
    arguments["--ineq"] = 'no' 
if arguments["--cond"] ==  None:
    rcond = 1e-6
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
if arguments["--bperp"] ==  None:
    bperp = 'no'
else:
    bperp = arguments["--bperp"]
if arguments["--steps"] ==  None:
    cos = []
else:
    cos = list(map(float,arguments["--steps"].replace(',',' ').split()))
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
if arguments["--scale"] ==  None:
    scale = 1
else:
    scale = float(arguments["--scale"]) 
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
if arguments["--ndatasets"] ==  None:
    ndata = 1
else:
    ndata = int(arguments["--ndatasets"])

if len(pos)>0 and len(cos) != len(pos):
    raise Exception("coseimic and postseismic lists are not the same size")

markers = ['o','v','^','s','P','X','o','v','^','s','P','X']
#markers = ['o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o']

ipix = list(map(int,arguments["--cols"].replace(',',' ').split()))
jpix = list(map(int,arguments["--lines"].replace(',',' ').split()))
if len(jpix) != len(ipix):
   raise Exception("ncols and nlines lists are not the same size")
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

# Définir datemin comme référence temporelle pour les fonctions de base
# et dmin, dmax comme limites des figures
if arguments["--dateslim"] is not None:
    dmin, dmax = arguments["--dateslim"].replace(',', ' ').split()
    datemin = int(date2dec(dmin)[0])
    datemax = int(date2dec(dmax)[0]) + 1
    dmin = str(datemin) + '0101'
    dmax = str(datemax) + '0101'
else:
    datemin, datemax = int(np.min(dates)), int(np.max(dates)) + 1
    dmin = str(int(np.min(dates))) + '0101'
    dmax = str(int(np.max(dates)) + 1) + '0101'

if arguments["--plot_dateslim"] is not None:
    dmin, dmax = arguments["--plot_dateslim"].replace(',', ' ').split()

# Filtrer les dates selon les limites définies
indexd = np.flatnonzero(np.logical_and(dates <= datemax, dates >= datemin))
idates, dates = idates[indexd], dates[indexd]
# Mettre à jour les bases en filtrant également selon indexd
bases = [bp[indexd] for bp in bases]

# lect cube
cubei = np.fromfile(cubef,dtype=np.float32)
cube = as_strided(cubei[:nlign*ncol*N])
print('Number of line in the cube: ', cube.shape)
kk = np.flatnonzero(np.logical_or(cube==9990, cube==9999))
cube[kk] = float('NaN')
cube = cube*scale
maps = cube.reshape((nlign,ncol,N))

# set to NaN lakes and ocean, new image ref.
cst = np.copy(maps[:,:,imref])
for l in range((N)):
    maps[:,:,l] = maps[:,:,l] - cst
    if l != imref:
        index = np.nonzero(maps[:,:,l]==0.0)
        maps[:,:,l][index] = float('NaN')

# new number of dates
N = len(dates)
maps = as_strided(maps[:,:,indexd])
print('Reshape cube: ', maps.shape)

if apsf is not None:
    inaps=np.loadtxt(apsf, comments='#', unpack=True,dtype='f')*abs(scale)
    print('Input uncertainties:', inaps)
    print('Set very low values to the 2 percentile to avoid overweighting...')
    minaps= np.nanpercentile(inaps,2)
    index = np.flatnonzero(inaps<minaps)
    inaps[index] = minaps
    inaps = inaps[indexd]
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

nfigure = 0
# plot pixels on map
fig = plt.figure(nfigure,figsize=(12,8))
nfigure += 1
vmax = np.nanpercentile(maps[:,:,-1],90)
vmin = np.nanpercentile(maps[:,:,-1],10)

ax = fig.add_subplot(1,2,1)
ax.imshow(maps[jstart:jend,istart:iend,-1], cmap=cm.rainbow, vmax=vmax, vmin=vmin, alpha=0.6, interpolation='none')
for i in range(len(ipix)):
    ax.scatter(ipix[i]-istart,jpix[i]-jstart,marker=markers[i],color='black',s=10.)
if iref is not None and jref is not None:
    ax.scatter(iref-istart,jref-jstart,marker='x',color='red',s=20.)
for i in range((Npix)):
    ax.text(ipix[i]-istart,jpix[i]-jstart,i)
plt.suptitle('Black cross: pixels, red cross: reference point')

ax = fig.add_subplot(1,2,2)
im = ax.imshow(maps[:,:,-1], cmap=cm.rainbow, vmax=vmax, vmin=vmin, alpha=0.6, interpolation='none')
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

if arguments["--bperp"]=='yes':
  for idx, bp in enumerate(bases):
    # Création d'une figure pour chaque bp dans bases
    fig = plt.figure(nfigure, figsize=(10, 4))
    nfigure += 1
    
    # Graphique 1: Baseline perpendiculaire en fonction du temps
    ax1 = fig.add_subplot(1, 2, 1)
    
    # selection des baselines non nulles
    index = np.flatnonzero(bp)
    dates_temp = dates[index]
    idates_temp = idates[index]
    bp_temp = bp[index] 
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
    fig.savefig('baseline_{}.eps'.format(idx), format='EPS', dpi=150)
    np.savetxt('bp_t_{}.in'.format(idx), np.vstack([dates_temp, bp_temp]).T, fmt='%.6f')

# 0
basis=[
      reference(name='reference',date=datemin,reduction='ref'),
      ]
index = len(basis)

if inter=='yes':
      basis.append(linear(name='linear',reduction='lin',date=datemin))
      # 1
      indexinter=index
      index = index + 1

if degreeday=='yes':
   indexdd = index
   basis.append(stefan(name='degree-day',reduction='ddt',date=datemin,tcar1=ddt,tcar2=ddf))
   index = index + 1

if seasonal=='yes':
   # 2
   indexseas = index
   basis.append(cosvar(name='seas. var (cos)',reduction='cos',date=datemin))
   basis.append(sinvar(name='seas. var (sin)',reduction='sin',date=datemin))
   index = index + 2

if seasonalt=='yes':
   # 2
   indexseast = index
   basis.append(cost(name='increased seas. var (cos)',reduction='cost',date=datemin))
   basis.append(sint(name='increased seas. var (sin)',reduction='sint',date=datemin))
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
   basis.append(steps(name='steps {}'.format(i),reduction='step{}'.format(i),date=cos[i])),
   index = index + 1

indexsse = np.zeros(len(sse_time))
for i in range(len(sse_time)):
    basis.append(slowslip(name='sse {}'.format(i),reduction='sse{}'.format(i),date=sse_time[i],tcar=sse_car[i])),
    indexsse[i] = int(index)
    index = index + 1

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

kernels=[]
if arguments["--bperp"]=='yes':
    indexbperp = []
    i = 0
    for bp,bref in zip(bases,baserefs):
        i = i + 1
        kernels.append(corrBperp(name='bperp correction', reduction=f'corrBperp{i}', bp0=bref, bp=bp))
        indexbperp.append(index)
        index = index + 1

if vectf != None:
  indexvect = np.zeros(len(vectf))
  fig = plt.figure(nfigure,figsize=(12,8))
  nfigure += 1 
  for i in range(len(vectf)):
    kernels.append(vector(name=vectf[i],reduction='vector_{}'.format(i),vect=v[i]))
    ax = fig.add_subplot(len(vectf),1,i+1)
    ax.plot(dates,kernels[-1].g(~np.isnan(dates)))
    indexvect[i] = index
    index = index + 1

indexpo = indexpo.astype(int)
indexco = indexco.astype(int)
indexsse = indexsse.astype(int)

equality = False
if arguments["--seasonal_increase"] == 'yes' and arguments["--seasonal"] == 'yes':
    equality = True
    arguments["--ineq"] = 'yes'

# define size G matrix
print()
Mbasis=len(basis)
print('Number of basis functions:', Mbasis)
Mker=len(kernels)
print('Number of kernel functions:', Mker)
M = Mbasis + Mker
for i in range((Mbasis)):
    basis[i].info()

def consInvert(A, b, sigmad, ineq='yes', cond=1e-6, iter=60, acc=5e-4, equality=False):
    """
    Résout Ax ≈ b sous contraintes d'inégalité (et éventuellement d’égalité).

    Retourne :
        fsoln : vecteur de solution
        sigmam : incertitudes (diag de la covariance)
    """
    global indexpo, indexco, indexpofull, pos, indexseas, indexseast
    if A.shape[0] != len(b):
        raise ValueError('Dimensions incompatibles pour A et b')

    if ineq == 'no':
        W = np.diag(1.0 / sigmad)
        fsoln = np.linalg.lstsq(W @ A, W @ b, rcond=cond)[0]
    else:
        # Initialisation
        if indexpo is not None and len(indexpo) > 0:
            Ain = np.delete(A, indexpo, axis=1)
            Win = np.diag(1.0 / np.delete(sigmad, indexpo))
            mtemp = np.linalg.lstsq(Win @ Ain, Win @ b, rcond=cond)[0]

            # Réinsérer les post-sismiques
            for z in range(len(indexpo)):
                mtemp = np.insert(mtemp, indexpo[z], 0.0)
            minit = np.copy(mtemp)
        else:
            W = np.diag(1.0 / sigmad)
            minit = np.linalg.lstsq(W @ A, W @ b, rcond=cond)[0]

        # Définir les bornes
        n = len(minit)
        mmin = -np.ones(n) * np.inf
        mmax = np.ones(n) * np.inf

        if indexpo is not None and indexco is not None:
            for i in range(len(indexco)):
                ico = int(indexco[i])
                ipo = int(indexpofull[i])
                if pos[i] > 0. and minit[ico] > 0.:
                    mmin[ipo], mmax[ipo] = 0, np.inf
                    mmin[ico], mmax[ico] = 0, minit[ico]
                elif pos[i] > 0. and minit[ico] < 0.:
                    mmin[ipo], mmax[ipo] = -np.inf, 0
                    mmin[ico], mmax[ico] = minit[ico], 0

        bounds = list(zip(mmin, mmax))

        # Fonction à minimiser
        def _func(x):
            return np.sum(((A @ x - b) / sigmad) ** 2)

        def _fprime(x):
            return 2 * A.T @ ((A @ x - b) / sigmad**2)

        # Contraintes d'égalité
        if equality:
            def eq_cond(x):
                return (x[indexseast + 1] / x[indexseast]) - (x[indexseas + 1] / x[indexseas])
            res = opt.fmin_slsqp(_func, minit, bounds=bounds, fprime=_fprime,
                                 eqcons=[eq_cond], iter=iter, acc=acc,
                                 full_output=True, iprint=0)
        else:
            res = opt.fmin_slsqp(_func, minit, bounds=bounds, fprime=_fprime,
                                 iter=iter, acc=acc, full_output=True, iprint=0)

        fsoln, fx, its, imode, msg = res
        if imode != 0:
            logger.warning("SLSQP did not converge: {msg}")
            fsoln = minit  # ou np.full_like(minit, np.nan)

    # Calcul de l'incertitude
    try:
        varx = np.linalg.pinv(A.T @ A)
        res2 = np.sum((b - A @ fsoln) ** 2)
        scale = 1. / (A.shape[0] - A.shape[1])
        sigmam = np.sqrt(scale * res2 * np.diag(varx))
    except np.linalg.LinAlgError:
        sigmam = np.full(A.shape[1], np.nan)

    return fsoln, sigmam


# plot diplacements maps
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

    x = [date2num(dt.strptime('{}'.format(d),'%Y%m%d')) for d in idates]
    if arguments["--plot_dateslim"] is not  None:
        dmin,dmax = arguments["--plot_dateslim"].replace(',',' ').split()
    xmin = dt.strptime('{}'.format(dmin),'%Y%m%d')
    xmax = dt.strptime('{}'.format(dmax),'%Y%m%d')
    xlim=date2num(np.array([xmin,xmax]))

    # extract data
    if w>0:
      wind = np.nanmedian(as_strided(maps[j-w:j+w+1,i-w:i+w+1,:]),axis=(0,1))
      if infof is not None:
        infm = np.nanmedian(info[j-w:j+w+1,i-w:i+w+1])
      if iref is not None:
          windref = as_strided(maps[jref-wref:jref+wref+1,iref-wref:iref+wref+1,:])
          dispref = np.nanmedian(windref,axis=(0,1))
      else:
          dispref = np.zeros((N))
    else:
      wind = as_strided(maps[j,i,:]) 
      if infof is not None:    
        infm = info[j,i]
      if iref is not None:
        dispref = maps[jref,iref,:]
      else:
        dispref = np.zeros((N))

    disp = wind - dispref

    # inversion model
    mdisp= np.ones((N), dtype=np.float32)*float('NaN')
    bperperr = np.zeros((N))
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
    m = np.ones((M), dtype=np.float32)*float('NaN')
    sigmam = np.ones((M), dtype=np.float32)*float('NaN')
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
        #print(len(k), N)
        #print(G[:6,:6])
        #print(taby)
        #print(sigmad[k])
        m,sigmam = consInvert(G,taby,sigmad[k],cond=rcond, ineq=arguments["--ineq"], equality=equality)
        #sys.exit()   
     
        # forward model in original order
        mdisp[k] = np.dot(G,m)

        if bperp=='yes':
            bperperr[k] =  np.dot(G[:,indexbperp[0]:indexbperp[0]+ndata],m[indexbperp[0]:indexbperp[0]+ndata])
        else:
            bperperr = np.zeros((N))
        
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

    # plot data and model minus bperp error and seasonal terms
    if seasonal=='yes':
            G=np.zeros((kk,2))
            for l in range((2)):
                G[:,l]=basis[l+indexseas].g(tabx)
            disp_seas[k] = disp_seas[k] + np.dot(G[:,:],m[indexseas:indexseas+2])
            amp = np.sqrt(m[indexseas]**2+m[indexseas+1]**2)
            phi = np.arctan2(m[indexseas+1],m[indexseas])
            #phi = np.arctan(m[indexseas+1]/m[indexseas])

            if phi<0: 
               phi = phi + 2*np.pi

    if seasonalt=='yes':
            G=np.zeros((kk,2))
            for l in range((2)):
                G[:,l]=basis[l+indexseast].g(tabx)
            disp_seas[k] = disp_seas[k] + np.dot(G[:,:],m[indexseast:indexseast+2])
            disp_seast[k] = disp_seast[k] + np.dot(G[:,:],m[indexseast:indexseast+2])

            phit = np.arctan2(m[indexseast+1],m[indexseast])
            #phit = np.arctan(m[indexseast+1]/m[indexseast])
            ampt = np.sqrt(m[indexseast]**2+m[indexseast+1]**2)
 
            # convert between 0 and 2pi
            if phit<0: 
                phit = phit + 2*np.pi

            print(phi,phit)
            print(m[indexseas+1]/m[indexseas],m[indexseast+1]/m[indexseast] )
 
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

    # Tracé des données corrigées de l'erreur de baseline perpendiculaire (`bperp`) et du modèle
    for jj, bp in enumerate(bases):
        # Remplacer les valeurs non nulles de `bp` par 1, et les zéros par NaN pour les masquer
        bp_masked = np.where(bp != 0, 1.0, np.nan)
        dispp = disp * bp_masked  # Sélectionner les valeurs de déplacement en fonction de `bp`

         # Tracé avec ou sans information additionnelle `infof`
        if infof is not None:
            ax.plot(x, dispp - bperperr,
            markers[jj],
            color=color,
            fillstyle=fillstyle,
            label=f'TS {jj}: lines:{i}, column:{j}, Info:{infm:.2f}'
            )
        else:
            ax.plot(x, dispp - bperperr,
            markers[jj],
            color=color,
            fillstyle=fillstyle,
            label=f'TS {jj}: lines:{i}, column:{j}'
            )

         # Ajout des barres d'erreur
        ax.errorbar(x,dispp - bperperr,
        yerr=sigmad,
        ecolor=color,
        fmt='none',
        alpha=0.5
        )
  
    # plot data and model minus bperp error and linear term
    if inter=='yes':
        ax3.plot(x,disp-bperperr-lin,markers[jj],color=color,fillstyle=fillstyle,markersize=4,label='detrended data')
        ax3.errorbar(x,disp-bperperr-lin,yerr = sigmad, ecolor=color,fmt='none', alpha=0.5)
            
    if semianual=='yes' or seasonal=='yes' or bianual=='yes' or seasonalt == 'yes':
        ax2.plot(x,disp-disp_seas-bperperr,markers[jj],color=color,fillstyle=fillstyle,label='data -seasonal')
        # ax2.plot(x,mdisp-disp_seas-bperperr,'o',color='red',alpha=0.5,label='model -seasonal')
        ax2.errorbar(x,disp-disp_seas-bperperr,yerr = sigmad, ecolor=color,fmt='none', alpha=0.3)

    # create synthetic time
    #t = np.arange(xmin, xmax, 0.1)
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
   
    # we need to substrat bperp model to the whole model 
    if bperp=='yes':
        model_bperp = np.dot(G[:,indexbperp[0]:indexbperp[0]+ndata],m[indexbperp[0]:indexbperp[0]+ndata])
    else:
        model_bperp = np.zeros(len(model))
        
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

    if bianual=='yes':
        G=np.zeros((len(tdec),2))
        for l in range((2)):
            G[:,l]=basis[l+indexbi].g(tdec)
        mseas = mseas + np.dot(G[:,:],m[indexbi:indexbi+2])

    # plot model
    if inter=='yes':
        if seasonal=='yes' and seasonalt=='no':
            ax.plot(t,model-model_bperp,'-r',label='Rate: {:.2f}, Amp: {:.2f}, Phi: {:.2f}'.format(m[indexinter],amp,phi))
            ax3.plot(t,model-model_lin-model_bperp,'-r')
        elif seasonal=='yes' and seasonalt=='yes':
            ax.plot(t,model-model_bperp,'-r',label='Rate: {:.2f}, Ampt: {:.2f}, Phit: {:.2f}, Amp: {:.2f}, Phi: {:.2f}'.format(m[indexinter],ampt,phit,amp,phi))
            ax3.plot(t,model-model_lin-model_bperp,'-r')
        elif seasonal=='no' and seasonalt=='yes':
            ax.plot(t,model-model_bperp,'-r',label='Rate: {:.2f}, Ampt: {:.2f}, Phit: {:.2f}'.format(m[indexinter],ampt,phit))
            ax3.plot(t,model-model_lin-model_bperp,'-r')
        else:
            ax.plot(t,model-model_bperp,'-r',label='Rate: {:.2f}'.format(m[indexinter]))
            ax3.plot(t,model-model_lin-model_bperp,'-r')
    else:
        ax.plot(t,model-model_bperp,'-r')
           
    if seasonal=='yes' or semianual=='yes' or bianual=='yes' or seasonalt=='yes':
        ax2.plot(t,model-mseas-model_bperp,'-r')
        ax2.legend(loc='best',fontsize='x-small')
        ax2.set_xlim(xlim)

    ax.legend(loc='best',fontsize='x-small')
    ax.set_xlim(xlim)
    if arguments["--bounds"] is not  None:
        ax.set_ylim(ylim)
    
    if inter=='yes':
        ax3.legend(loc='best',fontsize='x-small')
        ax3.set_xlim(xlim)

    fig.autofmt_xdate()
    ax.set_xlabel('Time (Year/month/day)')
    ax.set_ylabel('Displacements')
    fig.savefig('Disp_{}.pdf'.format(output), format='PDF')

    if inter=='yes':    
        fig3.autofmt_xdate()
        ax3.set_xlabel('Time (Year/month/day)')
        ax3.set_ylabel('Displacements')
        fig3.savefig('Disp_{}_detrended.pdf'.format(output), format='PDF')

    if seasonal == 'yes' or semianual =='yes' or bianual =='yes':
        fig2.autofmt_xdate()
        ax2.set_xlabel('Time (Year/month/day)')
        ax2.set_ylabel('Displacements')
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
