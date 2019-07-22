#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
import datetime, sys
import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pylab import date2num
import matplotlib.ticker as ticker

# open .nc and select temp
ds = gdal.Open("NETCDF:"+sys.argv[1]+":t2m", gdal.GA_ReadOnly)
print 'Coordinates:', ds.GetGeoTransform()
print 'Grid:', ds.RasterYSize,ds.RasterXSize

# convert dates
date_ref = datetime.datetime(1900, 01, 01)

temp = np.zeros((ds.RasterYSize,ds.RasterXSize,ds.RasterCount)) 
dates,time,mode=[],[],[]

for b in range(1, ds.RasterCount+1):
    band = ds.GetRasterBand(b)
    bandmd = band.GetMetadata()
    # print bandmd
    
    # extract time
    netcdf_time = int(bandmd["NETCDF_DIM_time"])
    # print netcdf_time
    # print netcdf_time/24, netcdf_time%24
    date = date_ref + datetime.timedelta(netcdf_time/24, netcdf_time%24)
    dates.append(date)
    # convert date to dec
    dec = float(date.strftime('%j'))/365.1
    year = float(date.strftime('%Y'))
    # print date,dec,year
    mode.append(dec)
    time.append(year+dec)
    
    # extract temp
    scale = float(bandmd["scale_factor"])
    offset = float(bandmd["add_offset"])
    data = band.ReadAsArray()*scale + offset
    
    #print date, " -> ", np.mean(data), "m"
    #temp.append(data[-1,-1]) 
    # K to C
    temp[:,:,b-1] = data-273.15    
    #sys.exit()

time,dates,mode = np.array(time), np.array(dates),np.array(mode)
#temp = temp[4,3,:] 

# Initialisation
nfigure = 0
phase,amplitude = np.zeros((ds.RasterYSize,ds.RasterXSize)),np.zeros((ds.RasterYSize,ds.RasterXSize))
datemin, datemax= np.int(np.min(time)), np.int(np.min(time))+1
doplot='yes'

def seasonal(t,v,a,b,c):
    return v*(t-datemin) + a*np.cos(2*np.pi*(time-datemin)) + b*np.sin(2*np.pi*(time-datemin)) + c

def invers(temp,time,datemin):
    # fit
    G = np.zeros((len(temp),4))
    G[:,0] = time-datemin
    G[:,1] = np.cos(2*np.pi*(time-datemin))
    G[:,2] = np.sin(2*np.pi*(time-datemin))
    G[:,3] = 1

    pars = np.dot(np.dot(np.linalg.inv(np.dot(G.T,G)),G.T),temp)
    v,a,b,c = pars[0],pars[1],pars[2],pars[3]
    fit = v*(time-datemin) + a*np.cos(2*np.pi*(time-datemin)) + b*np.sin(2*np.pi*(time-datemin)) + c
    return v, a, b, c 

def plot(dates,time,mode,temp,fit,amp,dphi,nfigure,datemin,datemax):
    # plot data
    fig = plt.figure(nfigure,figsize=(7,6))
    ax = fig.add_subplot(1,1,1)
    # convert idates to num
    x = date2num(dates)
    xmin,xmax=datetime.datetime(np.int(datemin), 01, 01),datetime.datetime(datemax, 01, 01)
    xlim=date2num(np.array([xmin,xmax]))
    # format the ticks
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
    ax.plot(x,temp, '-bo',ms=.5,lw=.1,color='blue',label='ERAI point {}, {}'.format(lat,lon))
    ax.plot(x,fit,'-r',lw=2.0,label='v: {:0.2f} deg/yrs, dphi:{:0.2f} rad, Amp:{:0.2f} deg'.format(v,dphi,amp))
    #ax.plot(xlim,[273.15,273.15],'-b')
    # axes up to make room for them
    fig.autofmt_xdate()
    ax.set_xlim(xlim)
    ax.set_ylim([-30,30])
    plt.xlabel('Time (Year/month/day)')
    plt.ylabel('Temperature (deg)')
    plt.legend(loc='best')
    fig.savefig('temperatures_{}_{}.pdf'.format(lat,lon), format='PDF')
    nfigure=+1 
    # plot mod
    fig2 = plt.figure(nfigure,figsize=(7,5))
    ax2 = fig2.add_subplot(1,1,1)
    ax2.set_xlim([0,1])
    ax2.plot(mode,temp, '.k',ms=0.5,color='blue',label='ERAI point {}, {}'.format(lat,lon))
    #ax2.plot([0,1],[273.15,273.15],'-b')
    ax2.set_xticks(np.arange(0, 1, 1./12))
    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    plt.xlabel('Time')
    plt.ylabel('Temperature (deg)')
    plt.legend(loc='best')
    fig2.savefig('Modtemperatures_{}_{}.pdf'.format(lat,lon), format='PDF')

for i in xrange(ds.RasterYSize):
    lat = ds.GetGeoTransform()[3]-0.125 + ds.GetGeoTransform()[5]*i
    for j in xrange(ds.RasterXSize):
        lon = ds.GetGeoTransform()[0]+0.125 + ds.GetGeoTransform()[1]*j
        print 'lon: {}, lat: {}'.format(lon,lat)
        t = temp[i,j,:]
        v, a, b, c = invers(t,time,datemin)
        amp,dphi=np.sqrt(a**2+b**2),np.arctan2(b,a)
        if dphi < 0:
            dphi = dphi + 2*np.pi
        phase[i,j], amplitude[i,j] = dphi,amp
        print('Mean T: {:0.3f} °'.format(np.mean(t)))
        print('Vitesse: {:0.3f} °/yr'.format(v))
        print('Amplitude: {:0.3f} °'.format(np.sqrt(a**2+b**2)))
        print('Phase Shift: {:0.3f} rad'.format(np.arctan2(b,a)))
        print
        if doplot=='yes':
            fit = seasonal(time,v,a,b,c)
            plot(dates,time,mode,t,fit,amp,dphi,nfigure,datemin,datemax)
            nfigure+=2
            # plt.show()
            # sys.exit()

        # save data
        fid=open('temperatures_{}_{}.txt'.format(lat,lon),'wb')
        np.savetxt(fid,np.vstack([time,t]).T,fmt='%.6f',delimiter='\t',newline='\n')

# create new GDAL NCDF
drv = gdal.GetDriverByName("NETCDF")
dst_ds = drv.Create('t2m_amp_phi.nc', ds.RasterXSize, ds.RasterYSize, 2, gdal.GDT_Float32)
dst_ds.SetGeoTransform(ds.GetGeoTransform())
# Sortir les 2 bandes
dst_bandA = dst_ds.GetRasterBand(1)
dst_bandP = dst_ds.GetRasterBand(2)
# Write 
dst_bandA.WriteArray(amplitude,0,0)
dst_bandP.WriteArray(phase,0,0)
# close image
del dst_ds

plt.show()

