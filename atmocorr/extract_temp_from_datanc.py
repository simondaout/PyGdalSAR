#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Authors        : Mathieu Volat
#                  Simon DAOUT (Oxford)
############################################

import datetime, sys
import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pylab import date2num
import matplotlib.ticker as ticker
import scipy.optimize as opt
import scipy.linalg as lst


def seasonal(x,v,a,b,c):
    return v*(x-date0) + a*np.cos(2*np.pi*(x-date0)) + b*np.sin(2*np.pi*(x-date0)) + c

def invers(data,time,date0):
    # fit
    G = np.zeros((len(data),4))
    G[:,0] = time-date0
    G[:,1] = np.cos(2*np.pi*(time-date0))
    G[:,2] = np.sin(2*np.pi*(time-date0))
    G[:,3] = 1

    x0 = lst.lstsq(G,data)[0]
    _func = lambda x: np.sum(((np.dot(G,x)-data))**2)
    pars = opt.least_squares(_func,x0,jac='3-point',loss='cauchy').x
    # pars = np.dot(np.dot(np.linalg.inv(np.dot(G.T,G)),G.T),temp)
    v,a,b,c = pars[0],pars[1],pars[2],pars[3]
    fit = v*(time-date0) + a*np.cos(2*np.pi*(time-date0)) + b*np.sin(2*np.pi*(time-date0)) + c
    return v, a, b, c 

def plot(ldates,time,mode,temp,fit,amp,dphi,lat,lon,level,nfigure):
    # plot data
    fig = plt.figure(nfigure,figsize=(10,6))
    ax = fig.add_subplot(1,1,1)
    x = date2num(ldates)
    xmin,xmax=datetime.datetime(date0, 01, 01),datetime.datetime(2011, 01, 01)
    xlim=date2num(np.array([xmin,xmax]))
    # format the ticks
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
    ax.plot(x,temp, '.k',ms=1.0,color='blue',label='ERAI point {}, {}, {}'.format(lat,lon,level),alpha=0.8, rasterized=True)
    ax.plot(x,fit,'-r',label='v: {:0.2f} deg/yrs, dphi:{:0.2f} rad, Amp:{:0.2f} deg'.format(v,dphi,amp))
    #ax.plot(xlim,[273.15,273.15],'-b')
    # axes up to make room for them
    fig.autofmt_xdate()
    ax.set_xlim(xlim)
    # ax.set_ylim([-30,30])
    plt.xlabel('Time (Year/month/day)')
    plt.ylabel('Temperature (deg)')
    plt.legend(loc='best')
    fig.savefig('temperatures_{}_{}_{}.pdf'.format(lat,lon,level), format='PDF')
    nfigure=+1 
   
    # plot mod
    fig2 = plt.figure(nfigure,figsize=(7,5))
    ax2 = fig2.add_subplot(1,1,1)
    ax2.set_xlim([0,1])
    ax2.plot(mode,temp, '.k',ms=0.8,alpha=0.7,color='blue',label='ERAI point {}, {}, {}'.format(lat,lon, level))
    #ax2.plot([0,1],[273.15,273.15],'-b')
    ax2.set_xticks(np.arange(0, 1, 1./12))
    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    plt.xlabel('Time')
    plt.ylabel('Temperature (deg)')
    plt.legend(loc='best')
    fig2.savefig('Modtemperatures_{}_{}_{}.pdf'.format(lat,lon,level), format='PDF')
    plt.show()
    plt.close() 


# open .nc and select temp
ds = gdal.Open("NETCDF:"+sys.argv[1], gdal.GA_ReadOnly)
dsmd = ds.GetMetadata()
ds_geo = ds.GetGeoTransform()
print('Coordinates:', ds.GetGeoTransform())
print('Grid:', ds.RasterYSize,ds.RasterXSize, ds.RasterCount)
print()

pix_az, pix_rg = np.indices((ds.RasterYSize,ds.RasterXSize))
llat,llon = ds_geo[3]+ds_geo[5]*pix_az, ds_geo[0]+ds_geo[1]*pix_rg

latbeg,latend,lonbeg,lonend = 27 , 28.25, 89, 90.25
index = np.nonzero(
        np.logical_and(llat>latbeg,
        np.logical_and(llat<latend, 
        np.logical_and(llon>lonbeg,llon<lonend,
	))))
N = np.shape(index)[1]
lat,lon = llat[index],llon[index]

# convert dates
date_ref = datetime.datetime(1900, 01, 01)

# number of pressure levband = ds.GetRasterBand(b)el
levels_values = map(int,dsmd["NETCDF_DIM_level_VALUES"].replace('{','').replace('}','').replace(',',' ').split())
print('Pressure levels values:', levels_values)
M = len(levels_values)
temp = np.zeros((N,ds.RasterCount)) 
dates,time,mode,levels=[],[],[],[]

for b in range(1, ds.RasterCount+1):
    band = ds.GetRasterBand(b)
    bandmd = band.GetMetadata()

    netcdf_time = int(bandmd["NETCDF_DIM_time"])
    date = date_ref + datetime.timedelta(netcdf_time/24, netcdf_time%24)
    dates.append(date)
    # convert date to dec
    dec = float(date.strftime('%j'))/365.1
    year = float(date.strftime('%Y'))
    #print(date,dec,year)
    mode.append(dec)
    time.append(year+dec)
    
    # extract temp
    scale = float(bandmd["scale_factor"])
    offset = float(bandmd["add_offset"])
    data = band.ReadAsArray()[index]*scale + offset
   
    levels.append(bandmd["NETCDF_DIM_level"]) 
    #print(date, " -> ", np.mean(data))
    temp[:,b-1] = data-273.15    

time,dates,mode,levels = np.array(time), np.array(dates),np.array(mode),np.array(levels)

# initialise figure
nfigure = 0
date0 = 2007
doplot='yes'
   

for i in range(N):
# for i in range(2):
    for j in range(M):
        level = levels_values[j]
        print('lon: {}, lat: {}, level: {}'.format(lon[i],lat[i],level))
        data = temp[i,j::M]
        t = time[j::M]
        v, a, b, c = invers(data,t,date0)
        amp,dphi=np.sqrt(a**2+b**2),np.arctan2(b,a)
        print('Meandates T: {:0.3f} °'.format(np.mean(data)))
        print('Vitesse: {:0.3f} °/yr'.format(v))
        print('Amplitude: {:0.3f} °'.format(np.sqrt(a**2+b**2)))
        print('Phase Shift: {:0.3f} rad'.format(np.arctan2(b,a)))
        print()
        if doplot=='yes':
            fit = seasonal(t,v,a,b,c)
            plot(dates[j::M],t,mode[j::M],data,fit,amp,dphi,lat[i],lon[i],level,nfigure)
            nfigure+=2

        # save data
        fid=open('temperatures_{}_{}_{}.txt'.format(lat[i],lon[i],level),'wb')
        np.savetxt(fid,np.vstack([t,data]).T,fmt='%.6f',delimiter='\t',newline='\n')
        fid.flush()
        fid.close()

#plt.show()

