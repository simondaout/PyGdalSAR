#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Simon DAOUT (CPRG, Nancy)
############################################

"""\
plot_los_gps.py
-------------
Plot GPS time series in LOS

Usage: plot_los_gps.py --network=<path> --reduction=<path> [--wdir=<path>] [--extension=<value>] [--outputdir=<value>] 
plot_los_gps.py -h | --help

Options:
-h --help           Show this screen 
--network PATH      Text file of name lat, lon for each stations
--reduction PATH    Directory name containing the TS (each TS in a text file untitled by the GPS name: East, North, Down, sigma_East, sigma_North, sigma Down format)      
--wdir PATH         Path to the network file and reduction directory [default: ./]
--extension Value   Extention of the the GPS TS files [default: .neu]
--outputdir PATH    Output directory for pdf plots [default: ts_plots]
"""

print()
print()
print('Please cite:')
print('Daout, S., Steinberg, A., Isken, M. P., Heimann, S., & Sudhaus, H. (2020). Illuminating the spatio-temporal evolution of the 2008–2009 Qaidam earthquake sequence with the joint use of InSAR time series and teleseismic data. Remote Sensing, 12(17), 2850.')
print()
print()

import numpy as np
from numpy.lib.stride_tricks import as_strided
# scipy
import scipy
import scipy.optimize as opt
import scipy.linalg as lst
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import datetime, math, time
from os import path, makedirs
import sys
# docopt (command line parser)
import docopt

# read arguments
arguments = docopt.docopt(__doc__)
network = arguments["--network"]
reduction = arguments["--reduction"]
if arguments["--wdir"] ==  None:
    wdir = './'
else:
    wdir = arguments["--wdir"]
if arguments["--extension"] ==  None:
    ext = '.neu'
else:
    ext = arguments["--extension"]    
if arguments["--outputdir"] ==  None:
    out = 'ts_plots/'
else:
    out = arguments["--outputdir"] + '/'   

dim = 3 # cannot proj without vertical 

def makedir(name):
    if path.exists(name):
        return
    makedirs(name)

### Define GPS class
class point:
    def __init__(self,x,y):
        self.x=x
        self.y=y

class gpspoint(point):
    def __init__(self,x,y,name):
        point.__init__(self,x,y)
        self.name=name
        self.t=[]
        self.d=[] #self.east,self.north,self.up
        self.sigmad=[] #self.east,self.north, self.up

        # lenght of the time series
        self.Nt = 0

    def info(self):
        print('gps station: {}, self.east(km): {}, self.north(km): {}'.format(self.name,self.x,self.y))

def neu2los_roipac(vn, ve, vu, look, head):
    theta = np.deg2rad(90.-look)
    phi = np.deg2rad(-90-head)
    los = ve * np.cos(phi) * np.cos(theta) \
            + vn * np.sin(phi) * np.cos(theta) \
            + vu *np.sin(theta)
    #print([np.cos(phi) * np.cos(theta),np.sin(phi) * np.cos(theta),np.sin(theta)])
    #sys.exit()
    return los

class gpstimeseries:
    def __init__(self,network,reduction,dim,wdir,extension):

        self.network=network
        self.reduction=reduction
        self.dim=dim
        self.wdir=wdir
        self.extension=extension
    
        # inititialisation
        self.points=[]
        self.Npoints=0
        self.N=0
        self.x,self.y=[],[]
        self.tmin,self.tmax=10000000,0

    def load(self):
        fname=self.wdir + self.network
        if not path.isfile(fname):
            raise ValueError("invalid file name: " + fname)
        else:
            print('Load GPS time series: ', fname)
            print() 
        
        # name, self.east(km), self.north(km)
        name,x,y,look,heading=np.loadtxt(fname,comments='#',unpack=True,dtype='S4,f,f,f,f')

        self.name,self.x,self.y,self.look,self.heading=np.atleast_1d(name,x,y,look,heading)

        self.Npoints=len(self.name)
        self.N=0
        self.d = []
        
        print('Load time series... ')
        for i in range(self.Npoints):
            station=self.wdir+self.reduction+'/'+self.name[i].decode('utf-8')+self.extension 
            print(station)
            if not path.isfile(station):
                raise ValueError("invalid file name: " + station)
                pass
            else:
                self.points.append(gpspoint(self.x[i],self.y[i],self.name[i]))

            if 3==self.dim:
                dated,east,north,up,esigma,nsigma,upsigma=\
                np.loadtxt(station,comments='#',usecols=(0,1,2,3,4,5,6), unpack=True,\
                    dtype='f,f,f,f,f,f,f')
                
                self.points[i].dated,self.points[i].east,self.points[i].north,self.points[i].up,\
                self.points[i].esigma,self.points[i].nsigma,self.points[i].upsigma=\
                np.atleast_1d(dated,east,north,\
                    up,esigma,nsigma,upsigma)
                 
                self.points[i].los=neu2los_roipac(self.points[i].north, self.points[i].east, self.points[i].up, self.look[i], self.heading[i])
                self.points[i].lossigma=neu2los_roipac(self.points[i].nsigma, self.points[i].esigma, self.points[i].upsigma, self.look[i], self.heading[i])

                self.points[i].d=[self.points[i].east,self.points[i].north,self.points[i].up,self.points[i].los]
                self.points[i].sigmad=[self.points[i].esigma,self.points[i].nsigma,self.points[i].upsigma,self.points[i].lossigma]
                self.points[i].comp=['east','north','up','los']
                self.points[i].t=np.atleast_1d(self.points[i].dated)


            self.points[i].tmin,self.points[i].tmax=min(self.points[i].dated),max(self.points[i].dated)
            if self.points[i].tmin < self.tmin:
                self.tmin = self.points[i].tmin
            if self.points[i].tmax > self.tmax:
                self.tmax = self.points[i].tmax
            
            self.points[i].Nt=len(self.points[i].t)
            self.N += int(len(self.points[i].t)*self.dim)
            self.d.append(self.points[i].d)

        self.d = np.array(self.d).flatten()
        self.sigmad = np.ones(self.N)

    def info(self):
        print()
        print('GPS time series from network:',self.network)
        print('Number of stations:', self.Npoints)
        print('Lenght data vector:', self.N)


if arguments["--extension"] ==  None:
    ext = '.neu'
else:
    ext = arguments["--extension"]    

dim = 3 # cannot proj without vertical 

def makedir(name):
    if path.exists(name):
        return
    makedirs(name)

### Define GPS class
class point:
    def __init__(self,x,y):
        self.x=x
        self.y=y

class gpspoint(point):
    def __init__(self,x,y,name):
        point.__init__(self,x,y)
        self.name=name
        self.t=[]
        self.d=[] #self.east,self.north,self.up
        self.sigmad=[] #self.east,self.north, self.up

        # lenght of the time series
        self.Nt = 0

    def info(self):
        print('gps station: {}, self.east(km): {}, self.north(km): {}'.format(self.name,self.x,self.y))

def neu2los_roipac(vn, ve, vu, look, head):
    theta = np.deg2rad(90.-look)
    phi = np.deg2rad(-90-head)
    los = ve * np.cos(phi) * np.cos(theta) \
            + vn * np.sin(phi) * np.cos(theta) \
            + vu *np.sin(theta)
    #print([np.cos(phi) * np.cos(theta),np.sin(phi) * np.cos(theta),np.sin(theta)])
    #sys.exit()
    return los

class gpstimeseries:
    def __init__(self,network,reduction,dim,wdir,extension):

        self.network=network
        self.reduction=reduction
        self.dim=dim
        self.wdir=wdir
        self.extension=extension
    
        # inititialisation
        self.points=[]
        self.Npoints=0
        self.N=0
        self.x,self.y=[],[]
        self.tmin,self.tmax=10000000,0

    def load(self):
        fname=self.wdir + self.network
        if not path.isfile(fname):
            raise ValueError("invalid file name: " + fname)
        else:
            print('Load GPS time series: ', fname)
            print() 
        
        # name, self.east(km), self.north(km)
        name,x,y,look,heading=np.loadtxt(fname,comments='#',unpack=True,dtype='S4,f,f,f,f')

        self.name,self.x,self.y,self.look,self.heading=np.atleast_1d(name,x,y,look,heading)

        self.Npoints=len(self.name)
        self.N=0
        self.d = []
        
        print('Load time series... ')
        for i in range(self.Npoints):
            station=self.wdir+self.reduction+'/'+self.name[i].decode('utf-8')+self.extension 
            print(station)
            if not path.isfile(station):
                raise ValueError("invalid file name: " + station)
                pass
            else:
                self.points.append(gpspoint(self.x[i],self.y[i],self.name[i]))

            if 3==self.dim:
                dated,east,north,up,esigma,nsigma,upsigma=\
                np.loadtxt(station,comments='#',usecols=(0,1,2,3,4,5,6), unpack=True,\
                    dtype='f,f,f,f,f,f,f')
                
                self.points[i].dated,self.points[i].east,self.points[i].north,self.points[i].up,\
                self.points[i].esigma,self.points[i].nsigma,self.points[i].upsigma=\
                np.atleast_1d(dated,east*1000,north*1000,\
                    up*1000,esigma*1000,nsigma*1000,upsigma*1000)
                 
                self.points[i].los=neu2los_roipac(self.points[i].north, self.points[i].east, self.points[i].up, self.look[i], self.heading[i])
                self.points[i].lossigma=neu2los_roipac(self.points[i].nsigma, self.points[i].esigma, self.points[i].upsigma, self.look[i], self.heading[i])

                self.points[i].d=[self.points[i].east,self.points[i].north,self.points[i].up,self.points[i].los]
                self.points[i].sigmad=[self.points[i].esigma,self.points[i].nsigma,self.points[i].upsigma,self.points[i].lossigma]
                self.points[i].comp=['east','north','up','los']
                self.points[i].t=np.atleast_1d(self.points[i].dated)


            self.points[i].tmin,self.points[i].tmax=min(self.points[i].dated),max(self.points[i].dated)
            if self.points[i].tmin < self.tmin:
                self.tmin = self.points[i].tmin
            if self.points[i].tmax > self.tmax:
                self.tmax = self.points[i].tmax
            
            self.points[i].Nt=len(self.points[i].t)
            self.N += int(len(self.points[i].t)*self.dim)
            self.d.append(self.points[i].d)

        self.d = np.array(self.d).flatten()
        self.sigmad = np.ones(self.N)

    def info(self):
        print()
        print('GPS time series from network:',self.network)
        print('Number of stations:', self.Npoints)
        print('Lenght data vector:', self.N)


manifold = gpstimeseries(network=network,reduction=reduction,dim=dim,wdir=wdir,extension=ext)
manifold.load()
manifold.info()

#datemin, datemax = np.int(np.min(manifold.tmin)), np.int(np.max(manifold.tmax))+1 
datemin, datemax= 2015, 2022

# plot diplacements maps
nfigure = 0

colors = ['royalblue','darkgreen','darkorchid','firebrick']

ts_plots = path.join(wdir, out)
makedir(ts_plots)

for jj in range((manifold.Npoints)):
    fig = plt.figure(nfigure+jj,figsize=(12,8))
    name, i, j = manifold.name[jj], manifold.x[jj], manifold.y[jj]
    print('station:', name,i,j)
    fig.suptitle('{0} {1:.3f} {2:.3f}'.format(name.decode('utf-8'),i,j))
    pt = manifold.points[jj]
    
    # plot date
    i = 0
    for disp, sigma, name in zip(pt.d, pt.sigmad, pt.comp):
        ax = fig.add_subplot(len(pt.comp),1,i+1)
        ax.plot(pt.t, disp, 'o', color=colors[i], ms=2)
        ax.errorbar(pt.t,disp,yerr=sigma,ecolor=colors[i],fmt='none', alpha=0.1)
        ax.set_ylabel('{} (mm)'.format(name))
        ax.set_xlim([datemin,datemax])
        ax.set_ylim([np.nanpercentile(disp,.2),np.nanpercentile(disp,99.8)])
        #ax.legend(loc='best',fontsize='small')
        ax.grid(True)
        i = i+1

    ax.set_xlabel('Time (Year/month/day)')
    fig.savefig(out+'Disp_{}.pdf'.format(pt.name.decode('utf-8')), format='PDF')
    del fig

# plt.legend(loc='best')
#plt.show()
sys.exit()
