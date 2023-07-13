#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT (Oxford)
################################################################################


"""\
correct_ts_from_model.py
-------------
Correct time series from coefficienf file obtained with invers_disp2coef.py
Input crep maps must be named *_coef_clean.r4  

Usage: correct_ts_from_model.py [--cube=<path>] [--list_images=<path>] [--lectfile=<path>] \
 [--slope=<path>] [--coseismic=<paths>] [--postseismic=<paths>] [--slowslip=<value>] \
 [--cos=<path>] [--sin=<path>] [--dem=<path>] [--aps=<path>] \
 [--rad2mm=<value>] [--imref=<value>] [--crop=<values>]\
[--dateslim=<values>] [--cols=<values>] [--ligns=<values>] [--bounds=<value>]


Options:
-h --help           Show this screen
--cube PATH         Path to displacement file [default: depl_cumul_flat]
--list_images PATH  Path to list images file made of 5 columns containing for each images 1) number 2) Doppler freq (not read) 3) date in YYYYMMDD format 4) numerical date 5) perpendicular baseline [default: images_retenues]
--lectfile PATH     Path of the lect.in file [default: lect.in]
--slope PATH        Path to velocity map in r4 or tif format (e.g lin_coeff.r4) [default: None]
--cos PATH          Path to cosinus map in r4 or tif format (e.g coswt_coeff.r4) [default: None]
--sin PATH          Path to sinus map in r4 or tif format (e.g sinwt_coeff.r4)  [default: None]
--coseismic         Time of the events (e.g 2013.2,2014.5) [default: None]
--postseismic       Characteristic time of the transient functions (e.g 1.,0.5) [default: None]
--slowslip   VALUE  Read slow-slip maps. Indicate median and characteristic time of the events
--dem PATH          Path to dem file error map in r4 or tif format (demerr_coeff.r4) [default: None]
--aps PATH          Path to aps file: 1 column file given aps for each dates (eg. aps.txt) [default: None]
--rad2mm            Scaling value between input data (rad) and desired output [default: -4.4563]
--imref VALUE       Reference image number [default: 1]
--cols VALUE        Pixels column numbers to display (eg. 200,400,450) 
--ligns VALUE       Pixels lign numbers pixel to display (eg. 1200,1200,3000)
--crop VALUE        Crop option [default: 0,nlign,0,ncol] 
--bounds            yMin,yMax time series plots
"""

print()
print()
print('Please cite:')
print ('Daout, S., Doin, M. P., Peltzer, G., Socquet, A., & Lasserre, C. (2017). Large‐scale InSAR monitoring of permafrost freeze‐thaw cycles on the Tibetan Plateau. Geophysical Research Letters, 44(2), 901-909.')
print('Daout, S., Sudhaus, H., Kausch, T., Steinberg, A., & Dini, B. (2019). Interseismic and postseismic shallow creep of the North Qaidam Thrust faults detected with a multitemporal InSAR analysis. Journal of Geophysical Research: Solid Earth, 124(7), 7259-7279.')
print()
print()

import numpy as np
from numpy.lib.stride_tricks import as_strided
import os
import matplotlib
if os.environ["TERM"].startswith("screen"):
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import *
import math
import sys
import matplotlib.dates as mdates
import datetime

###############################################################
# define basis functioGns for plot
def Heaviside(t):
    h=np.zeros((len(t)))
    h[t>=0]=1.0
    return h

def coseismic(time,to,step):
    return step*Heaviside(time-to)

def postseismic(time,to,tcar,trans):
       t=(time-to)/tcar
       t[t<=0] = 0
       g = np.log10(1+t)
       return trans*g

def linear(time,v):
    return v*(time-datemin)

def seasonal(time,a,b):
    return a*np.cos(2*np.pi*(time-datemin)) + b*np.sin(2*np.pi*(time-datemin))

def slowslip(time,to,tcar):
          t=(time-to)/tcar
          funct = 0.5*(np.tanh(t)-1) + 1
          return funct
###############################################################

# docopt (command line parser)
import docopt

# read arguments
arguments = docopt.docopt(__doc__)
if arguments["--list_images"] ==  None:
    listim = "images_retenues"
else:
    listim = arguments["--list_images"]
if arguments["--cube"] ==  None:
   cubef = "depl_cumule_flat"
else:
   cubef = arguments["--cube"]
if arguments["--lectfile"] ==  None:
   infile = "lect.in"
else:
   infile = arguments["--lectfile"]
if arguments["--slope"] ==  None:
   slopef = None
else:
   slopef = arguments["--slope"]

if arguments["--cos"] ==  None:
   cosf = None
else:
   cosf = arguments["--cos"]
if arguments["--sin"] ==  None:
   sinf = None
else:
   sinf = arguments["--sin"]
if sinf is None and cosf is not None or cosf is None and sinf is not None:
    raise Exception("cosinus and sinus files must be given")

if arguments["--coseismic"] ==  None:
   cotimes = []
else:
   cotimes = list(map(float,arguments["--coseismic"].replace(',',' ').split()))
if arguments["--postseismic"] ==  None:
   postimes = []
else:
   postimes = list(map(float,arguments["--postseismic"].replace('None','-1').replace(',',' ').split()))

if len(postimes)>0 and len(cotimes) != len(postimes):
    raise Exception("coseimic and postseismic lists are not the same size")
if arguments["--dem"] ==  None:
   demf = None
else:
   demf = arguments["--dem"]
if arguments["--rad2mm"] ==  None:
    rad2mm = -56/(4*np.pi)
else:
    rad2mm = float(arguments["--rad2mm"]) 
if arguments["--imref"] ==  None:
    imref = 0
elif arguments["--imref"] < 1:
    print('--imref must be between 1 and Nimages')
else:
    imref = int(arguments["--imref"]) - 1

if arguments["--bounds"] is not  None:
    ylim = list(map(float,arguments["--bounds"].replace(',',' ').split()))

if arguments["--slowslip"] == None:
    sse=[]
else:
    sse = list(map(float,arguments["--slowslip"].replace(',',' ').split())) 
sse_time = sse[::2]
sse_car = sse[1::2]  

###############################################################

# create a list of pixels
jpix = list(map(int,arguments["--cols"].replace(',',' ').split()))
ipix = list(map(int,arguments["--ligns"].replace(',',' ').split()))
if len(jpix) != len(ipix):
    raise Exception("ncols and nligns lists are not the same size")
# number of pixels
Npix = len(ipix)

# load images_retenues file
# load images_retenues file
nb,idates,tdec,base=np.loadtxt(listim, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
datemin, datemax = int(np.min(tdec)), int(np.max(tdec))+1
N = len(idates)
base = base - base[imref]
print('Number images: ', N)
# read lect.in 
ncol, nlign = list(map(int, open(infile).readline().split(None, 2)[0:2]))

# lect cube
#maps = np.fromfile(cubef,dtype=np.float32)[:nlign*ncol*N].reshape((nlign,ncol,N))
cube = np.fromfile(cubef,dtype=np.float32)[:nlign*ncol*N]
# clean
kk = np.flatnonzero(cube>9990)
cube[kk] = float('NaN')
maps = cube.reshape((nlign,ncol,N))
del cube

if arguments["--crop"] ==  None:
    crop = [0,nlign,0,ncol]
else:
    crop = list(map(float,arguments["--crop"].replace(',',' ').split()))
ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])
    
if slopef is not None:
  extension = os.path.splitext(slopef)[1]
  if extension == ".tif":
    ds = gdal.Open(slopef, gdal.GA_ReadOnly)
    slopemap = ds.GetRasterBand(1).ReadAsArray()
  else:
    slopemap = np.fromfile(slopef,dtype=np.float32).reshape((nlign,ncol))
else:
    slopemap = np.zeros((nlign,ncol))

try:
  ds = gdal.Open('ref_coeff.tif', gdal.GA_ReadOnly)
  band = ds.GetRasterBand(1)
  refmap =  band.ReadAsArray()
except:
  refmap = np.fromfile('ref_coeff.r4',dtype=np.float32).reshape((nlign,ncol))

if demf is not None:
  extension = os.path.splitext(demf)[1]
  if extension == ".tif":
    ds = gdal.Open(demf, gdal.GA_ReadOnly)
    demmap = ds.GetRasterBand(1).ReadAsArray()
  else:
    demmap = np.fromfile(demf,dtype=np.float32).reshape((nlign,ncol))   
else:
    demmap = np.zeros((nlign,ncol))

if cosf is not None:
  extension = os.path.splitext(cosf)[1]
  if extension == ".tif":
    ds = gdal.Open(cosf, gdal.GA_ReadOnly)
    cosmap = ds.GetRasterBand(1).ReadAsArray()
  else:
    cosmap = np.fromfile(cosf,dtype=np.float32).reshape((nlign,ncol))   

if sinf is not None:
  extension = os.path.splitext(sinf)[1]
  if extension == ".tif":
    ds = gdal.Open(sinf, gdal.GA_ReadOnly)
    sinmap = ds.GetRasterBand(1).ReadAsArray()
  else:
    sinmap = np.fromfile(sinf,dtype=np.float32).reshape((nlign,ncol))

M = len(cotimes)
coseismaps=np.zeros((M,nlign,ncol))
postmaps=np.zeros((M,nlign,ncol))
# if no postseismic functions defined then set all functions to 0
if len(postimes) == 0:
    postimes = np.ones((M))*-1

for i in range(M):
  cofile = 'cos{}_coeff.r4'.format(i)
  try:
    ds = gdal.Open(cofile, gdal.GA_ReadOnly)
    coseismaps[i,:,:] = ds.GetRasterBand(1).ReadAsArray()
  except:
    coseismaps[i,:,:] = np.fromfile(cofile,dtype=np.float32).reshape((nlign,ncol))

for i in range((M)):
  postfile = 'post{}_coeff.r4'.format(i)
  if postimes[i] > 0:
    try:
      ds = gdal.Open(postfile, gdal.GA_ReadOnly)
      postmaps[i,:,:] = ds.GetRasterBand(1).ReadAsArray()
    except:  
      postmaps[i,:,:] = np.fromfile(postfile,dtype=np.float32).reshape((nlign,ncol))
  else: 
      postmaps[i,:,:] = np.zeros((nlign,ncol))
  
sse_times = sse[::2]
sse_car = sse[1::2] 
L = len(sse_times)
ssemaps=np.zeros((M,nlign,ncol))
for i in range(L):
    try:
      ds = gdal.Open('sse{}_coeff.tif'.format(i), gdal.GA_ReadOnly)
      ssemaps[i,:,:] = ds.GetRasterBand(1).ReadAsArray()
    except:
      ssemaps[i,:,:] = np.fromfile('sse{}_coeff.r4'.format(i),dtype=np.float32).reshape((nlign,ncol))


###############################################################
# figure ts
fig = plt.figure(0,figsize=(10,9))
fig.subplots_adjust(wspace=0.001)

# extract model
k = 0
model = np.zeros((nlign,ncol,N))
maps_clean = np.zeros((nlign,ncol,N))
demcor = np.zeros((nlign,ncol,N))

# plt.imshow(refmap)
# plt.show()

# Ini Model to DEMCOR for the whole map
# Remove ref Frame for the whole map
for l in range((N)):
    demcor[:,:,l] = demmap*base[l]
    model[:,:,l] =  demcor[:,:,l] + refmap


for i in range(ibeg,iend): # lines
    for j in range(jbeg,jend): # cols

        lin, a, b = 0, 0, 0
        steps, trans, sse = np.zeros(len(cotimes)), np.zeros(len(postimes)),np.zeros(len(sse_times))
        if slopef is not None:
            lin = slopemap[i,j]
        if cosf is not None:
            a,b = cosmap[i,j], sinmap[i,j]
        for l in range(len(cotimes)):
            steps[l] = coseismaps[l,i,j]
            trans[l] = postmaps[l,i,j]
        for l in range(len(sse_times)):
            sse[l] = ssemaps[l,i,j]
        ###################################################
        model[i,j,:] = model[i,j,:] + linear(tdec,lin) + seasonal(tdec, a, b)
        for l in range((M)):
            model[i,j,:] = model[i,j,:] + coseismic(tdec, cotimes[l], steps[l]) + \
            postseismic(tdec,cotimes[l],postimes[l],trans[l])
        for l in range((L)):
            model[i,j,:] = model[i,j,:] + sse[l]*slowslip(tdec,sse_times[l],sse_car[l])
        ####################################################
        maps_clean[i,j,:] = maps[i,j,:] - model[i,j,:]

        ###############################################################
        # plot TS
        if (i in ipix and j==jpix[ipix.index(i)]):
            ax = fig.add_subplot(Npix,1,k+1)
            print('plot TS: {0}-{1}'.format(i,j))
            x = [date2num(datetime.datetime.strptime('{}'.format(d),'%Y%m%d')) for d in idates]

            dmax = str(datemax) + '0101'
            dmin = str(datemin) + '0101'
            xmin = datetime.datetime.strptime('{}'.format(dmin),'%Y%m%d') 
            xmax = datetime.datetime.strptime('{}'.format(dmax),'%Y%m%d')
            xlim=date2num(np.array([xmin,xmax]))
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
            # t = np.array([xmin + datetime.timedelta(days=d) for d in range(0, 2920)]) 
            # print(coseismic(tdec, cotimes[0], steps[0])*rad2mm)
            # print(postseismic(tdec,cotimes[0],postimes[0],trans[0])*rad2mm)
            # print()
            ax.plot(x,(maps[i,j,:]-demcor[i,j])*rad2mm,'o',label='line: {}, col: {}'.format(i,j))
            if np.std(model[i,j,:]) > 0:
                plt.plot(x,(model[i,j,:]-demcor[i,j])*rad2mm,'-r')
            plt.legend(loc='best', fontsize='x-small')
            ax.set_xlim(xlim)
            if arguments["--bounds"] is not  None:
                ax.set_ylim(ylim)
            k = k+1


# save fig
fig.autofmt_xdate()
plt.xlabel('Time (Year/month/day)')
plt.ylabel('Displacements')

# save
fid = open('depl_cumule_nomodel', 'wb')
maps_clean.flatten().astype('float32').tofile(fid)
fid.close()

###############################################################
# plot maps
fig = plt.figure(2,figsize=(10,9))
fig.subplots_adjust(wspace=0.001)

maps_clean = maps_clean*rad2mm
maps = maps*rad2mm
model = model*rad2mm
if arguments["--bounds"] is not  None:
      vmax,vmin = np.nanmax(ylim), np.nanmin(ylim)
else:
    vmax = np.nanpercentile(maps_clean[:,:,-1],90)
    vmin = np.nanpercentile(maps_clean[:,:,-1],10)

for l in range((N)):
    d = as_strided(maps_clean[ibeg:iend,jbeg:jend,l])
    ax = fig.add_subplot(4,int(N/4)+1,l+1)
    cmap = cm.jet
    cmap.set_bad('white')
    cax = ax.imshow(d,cmap=cm.rainbow,vmax=vmax,vmin=vmin)
    ax.set_title(idates[l],fontsize=6)
    setp( ax.get_xticklabels(), visible=False)
    setp( ax.get_yticklabels(), visible=False)

setp(ax.get_xticklabels(), visible=False)
setp(ax.get_yticklabels(), visible=False)
fig.tight_layout()
plt.suptitle('Corrected TS')
fig.colorbar(cax, orientation='vertical',aspect=10)
fig.savefig('maps_clean_models.eps', format='EPS',dpi=150)


# plot diplacements maps
fig = plt.figure(3,figsize=(10,9))
fig.subplots_adjust(wspace=0.001)

for l in range((N)):
    d = as_strided(maps[ibeg:iend,jbeg:jend,l])
    ax = fig.add_subplot(4,int(N/4)+1,l+1)
    cmap = cm.jet
    cmap.set_bad('white')
    cax = ax.imshow(d,cmap=cm.rainbow,vmax=vmax,vmin=vmin)
    ax.set_title(idates[l],fontsize=6)
    setp( ax.get_xticklabels(), visible=False)
    setp( ax.get_yticklabels(), visible=False)

setp(ax.get_xticklabels(), visible=False)
setp(ax.get_yticklabels(), visible=False)
fig.tight_layout()
plt.suptitle('TS')
fig.colorbar(cax, orientation='vertical',aspect=10)
plt.show()




