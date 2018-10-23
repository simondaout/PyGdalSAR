#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT 
################################################################################

"""\
lect_disp_pixel.py
-------------
Plot time series results obtained with invers_disp2coef.py for given pixels (reqiered image_retenues file)

Usage: lect_cube_pixel.py --cols=<values> --ligns=<values> [--cube=<path>] [--lectfile=<path>] \
[--ref=<path>] [--slope=<path>] [--coseismic=<paths>] [--postseismic=<paths>] [--slowslip=<value>] \
 [--cos=<path>] [--sin=<path>] [--dem=<path>] [--aps=<path>] \
 [--rad2mm=<value>] [--plot=<yes/no>] [--name=<value>] [--imref=<value>] \
 [<iref>] [<jref>] [--bounds=<value>] 


Options:
-h --help           Show this screen
--ncols VALUE       Pixels column numbers (eg. 200,400,450) 
--nligns VALUE      Pixels lign numbers pixel (eg. 1200,1200,3000) 
--cube PATH         Path to displacement file [default: depl_cumul_flat]
--lectfile PATH     Path of the lect.in file [default: lect.in]
--ref PATH          Path to the reference map in r4 format (e.g ref_coeff.r4) [default: ref_coeff.r4]
--slope PATH        Path to velocity map in r4 format (e.g lin_coeff.r4) [default: None]
--cos PATH          Path to cosinus map in r4 format (e.g coswt_coeff.r4) [default: None]
--sin PATH          Path to sinus map in r4 format (e.g sinwt_coeff.r4)  [default: None]
--coseismic         Time of the events (e.g 2013.2,2014.5) [default: None]
--postseismic       Characteristic time of the transient functions (e.g 1.,0.5) [default: None]
--slowslip   VALUE  Read slow-slip maps. Indicate median and characteristic time of the events
--ref PATH          Path to reference file (ref_coeff.r4) [default: None]
--dem PATH          Path to dem file error map in r4 format (demerr_coeff.r4) [default: None]
--aps PATH          Path to aps file: 1 column file given aps for each dates (eg. aps.txt) [default: None]
--rad2mm                Scaling value between input data (rad) and desired output [default: -4.4563]
--name Value            Name output figures [default: None] 
--plot                  Display results [default: yes]            
--iref              colum numbers of the reference pixel [default: None] 
--jref              lign number of the reference pixel [default: None]
--imref VALUE           Reference image number [default: 1]
--bounds                Min,Max time series plots 
"""

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided

# basic
import math,sys,getopt
from os import path, environ
import os

# plot
import matplotlib
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.dates as mdates
from datetime import datetime
import datetime

# docopt (command line parser)
from nsbas import docopt

# read arguments
arguments = docopt.docopt(__doc__)
if arguments["--cube"] ==  None:
   cubef = "depl_cumul"
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
   cotimes = map(float,arguments["--coseismic"].replace(',',' ').split())
if arguments["--postseismic"] ==  None:
   postimes = []
else:
   postimes = map(float,arguments["--postseismic"].replace('None','-1').replace(',',' ').split())

if len(postimes)>0 and len(cotimes) != len(postimes):
    raise Exception("coseimic and postseismic lists are not the same size")

if arguments["--ref"] ==  None:
   reff = 'ref_coeff.r4'
else:
   reff = arguments["--ref"]

if arguments["--dem"] ==  None:
   demf = None
else:
   demf = arguments["--dem"]
if arguments["--aps"] ==  None:
   apsf = None
else:
   apsf = arguments["--aps"]
if arguments["<iref>"] ==  None:
    iref = None
else:
    iref = int(arguments["<iref>"])
if arguments["<jref>"] ==  None:
    jref = None
else:
    jref = int(arguments["<jref>"])

if arguments["--rad2mm"] ==  None:
    rad2mm = -4.4563
else:
    rad2mm = float(arguments["--rad2mm"]) 

if arguments["--name"] ==  None:
    output = None
else:
    output = arguments["--name"]

if arguments["--plot"] ==  None:
    plot = 'yes'
else:
    plot = arguments["--plot"]

if arguments["--imref"] ==  None:
    imref = 0
elif arguments["--imref"] < 1:
    print '--imref must be between 1 and Nimages'
else:
    imref = int(arguments["--imref"]) - 1

if arguments["--bounds"] is not  None:
    ylim = map(float,arguments["--bounds"].replace(',',' ').split())

if arguments["--slowslip"] == None:
    sse=[]
else:
    sse = map(float,arguments["--slowslip"].replace(',',' ').split()) 
sse_time = sse[::2]
sse_car = sse[1::2]  

# create a list of pixels
ipix = map(int,arguments["--cols"].replace(',',' ').split())
jpix = map(int,arguments["--ligns"].replace(',',' ').split())
if len(jpix) != len(ipix):
    raise Exception("ncols and nligns lists are not the same size")
# number of pixels
Npix = len(ipix)

# read lect.in 
ncol, nlign = map(int, open(infile).readline().split(None, 2)[0:2])

#initialize figures
nfigure=0

# load images_retenues file
fimages='images_retenues'
nb,idates,dates,base=np.loadtxt(fimages, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
datemin, datemax = np.int(np.min(dates)), np.int(np.max(dates))+1
N = len(dates)
base = base - base[imref]

# load
if apsf is not None:
    inaps=np.loadtxt(apsf, comments='#', unpack=True,dtype='f')*rad2mm

# lect cube
cube = np.fromfile(cubef,dtype=np.float32)*rad2mm
maps = cube.reshape((nlign,ncol,N))
print 'Reshape cube: ', maps.shape
listplot = [maps[:,:,-1]] 
titles = ['Depl. Cumul.']

if slopef is not None:
  slopemap = np.fromfile(slopef,dtype=np.float32).reshape((nlign,ncol))*rad2mm
  listplot.append(slopemap)
  titles.append('Slope')

if reff is not None:
  refmap = np.fromfile(reff,dtype=np.float32).reshape((nlign,ncol))*rad2mm

if demf is not None:
  demmap = np.fromfile(demf,dtype=np.float32).reshape((nlign,ncol))*rad2mm

if cosf is not None:
    cosmap = np.fromfile(cosf,dtype=np.float32).reshape((nlign,ncol))*rad2mm
    listplot.append(cosmap)
    titles.append('Cosinus')
if sinf is not None:
    sinmap = np.fromfile(sinf,dtype=np.float32).reshape((nlign,ncol))*rad2mm
    listplot.append(slopemap)
    titles.append('Phase')


M = len(cotimes)
coseismaps=np.zeros((M,nlign,ncol))
postmaps=np.zeros((M,nlign,ncol))
# if no postseismic functions defined then set all functions to 0
if len(postimes) == 0:
    postimes = np.ones((M))*-1

for i in xrange(M):
    coseismaps[i,:,:] = np.fromfile('cos{}_coeff.r4'.format(i),dtype=np.float32).reshape((nlign,ncol))*rad2mm
    listplot.append(coseismaps[i,:,:])
    titles.append('Coseism.{}'.format(i))
for i in xrange((M)):
    if postimes[i] > 0:
        postmaps[i,:,:] = np.fromfile('post{}_coeff.r4'.format(i),dtype=np.float32).reshape((nlign,ncol))*rad2mm
    else: 
        postmaps[i,:,:] = np.zeros((nlign,ncol))
    listplot.append(postmaps[i,:,:])
    titles.append('Post.{}'.format(i))

sse_times = sse[::2]
sse_car = sse[1::2] 

L = len(sse_times)
ssemaps=np.zeros((M,nlign,ncol))
for i in xrange(L):
    ssemaps[i,:,:] = np.fromfile('sse{}_coeff.r4'.format(i),dtype=np.float32).reshape((nlign,ncol))*rad2mm
    listplot.append(ssemaps[i,:,:])
    titles.append('sse.{}'.format(i))

# plot pixels on map
fig = plt.figure(0,figsize=(10,8))
# fig.subplots_adjust(hspace=.001,wspace=0.001)

for l in xrange(len(listplot)):
  ax = fig.add_subplot(2,int(len(listplot)/2)+1,l+1)
  vmax = np.nanpercentile(listplot[l],98)
  vmin = np.nanpercentile(listplot[l],2)
  ax.imshow(listplot[l], vmax=vmax, vmin= vmin, alpha=0.6)
  ax.scatter(ipix,jpix,marker='x',color='black',s=15.)
  ax.scatter(iref,jref,marker='x',color='red',s=20.)
  for i in xrange((Npix)):
      ax.text(ipix[i],jpix[i],i)
  ax.set_title(titles[l],fontsize=6)
  setp(ax.get_xticklabels(), visible=False)
  setp(ax.get_yticklabels(), visible=False)
# plt.show()
# sys.exit()

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

# plot diplacements maps
fig = plt.figure(1,figsize=(12,8))

for k in xrange(len(ipix)):
    i, j = ipix[k], jpix[k]

    ax = fig.add_subplot(Npix,1,k+1)
    x = [date2num(datetime.datetime.strptime('{}'.format(d),'%Y%m%d')) for d in idates]
    dmin = str(datemin) + '0101'
    dmax = str(datemax) + '0101'
    xmin = datetime.datetime.strptime('{}'.format(dmin),'%Y%m%d') 
    xmax = datetime.datetime.strptime('{}'.format(dmax),'%Y%m%d')
    xlim=date2num(np.array([xmin,xmax]))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))

    demcor = np.zeros((N))
    alpha = 0
    if demf is not None:
      if iref is not None:
          alpha = demmap[j,i] - demmap[jref,iref]
      else:
          alpha = demmap[j,i]
      demcor = alpha*base

    # plot data
    disp = as_strided(maps[j,i,:])
    if iref is not None:
        dispref = as_strided(maps[jref,iref,:])
    else:
        dispref = np.zeros((N))
    disp = disp - dispref

    ax.plot(x,disp-demcor,'o',label='TS {}: lign: {}, column: {}'.format(k,i,j))
    ax.errorbar(x,disp-demcor,yerr = inaps, ecolor='blue',fmt='none', alpha=0.5)

    # extract model
    ref, lin, a, b = 0, 0, 0, 0
    steps, trans, sse = np.zeros(len(cotimes)), np.zeros(len(postimes)),np.zeros(len(sse_times))

    if slopef is not None:
        if iref is not None:
            lin = slopemap[j,i] - slopemap[jref,iref]
        else:
            lin = slopemap[j,i]
    if cosf is not None:
        if iref is not None:
            a,b = cosmap[j,i] - cosmap[jref,iref], sinmap[j,i] - sinmap[jref,iref]
        else:
            a,b = cosmap[j,i], sinmap[j,i]
    if reff is not None:
        if iref is not None:
            ref = refmap[j,i] - refmap[jref,iref]
        else:
            ref = refmap[j,i]

    for l in xrange(len(cotimes)):
        if iref is not None:
            steps[l] = coseismaps[l,j,i] - coseismaps[l,jref,iref] 
            trans[l] =postmaps[l,j,i] - postmaps[l,jref,iref] 
        else:
            steps[l] = coseismaps[l,j,i]
            trans[l] = postmaps[l,j,i]
   
    for l in xrange(len(sse_times)):
      if iref is not None:
          sse[l] = ssemaps[l,j,i] - ssemaps[l,jref,iref] 
      else:
          sse[l] = ssemaps[l,j,i]

    t = np.array([xmin + datetime.timedelta(days=d) for d in range(0, 2920)])
    tdec = np.array([float(date.strftime('%Y')) + float(date.strftime('%j'))/365.1 for date in t])
   
    print
    print i,j 
    print  ref, lin, steps, trans

    model = ref + linear(tdec,lin) + seasonal(tdec, a, b) 
    for l in xrange((M)):
        model = model + coseismic(tdec, cotimes[l], steps[l]) + \
        postseismic(tdec,cotimes[l],postimes[l],trans[l])
    for l in xrange((L)):
        model = model + sse[l]*slowslip(tdec,sse_times[l],sse_car[l])

    if np.std(model) > 0:
        plt.plot(t,model,'-r')

    plt.legend(loc='best', fontsize='x-small')
    ax.set_xlim(xlim)
    if arguments["--bounds"] is not  None:
        ax.set_ylim(ylim)

fig.autofmt_xdate()
plt.xlabel('Time (Year/month/day)')
plt.ylabel('Displacements')
fig.savefig('Model_{}.eps'.format(output), format='EPS', dpi=150)

if plot == 'yes':
    plt.show()
sys.exit()


