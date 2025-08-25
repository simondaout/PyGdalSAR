#!/usr/bin/env python3

import datetime, sys
import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pylab import date2num
import matplotlib.ticker as ticker

# open .nc and select temp
ds = gdal.Open("NETCDF:"+sys.argv[1]+":t2m", gdal.GA_ReadOnly)
print('Coordinates:', ds.GetGeoTransform())
print('Grid:', ds.RasterYSize,ds.RasterXSize)

# convert dates
#date_ref = datetime.datetime(1900, 1, 1)

temp = np.zeros((ds.RasterYSize,ds.RasterXSize,ds.RasterCount)) 
dates,time,mode=[],[],[]

for b in range(1, ds.RasterCount+1):
    band = ds.GetRasterBand(b)
    bandmd = band.GetMetadata()
    # print(bandmd)
    
    # extract time
    print(bandmd.keys())
    netcdf_time = int(bandmd["NETCDF_DIM_valid_time"])
    #print(bandmd)
    date = datetime.datetime.utcfromtimestamp(netcdf_time)
    dates.append(date)
    # convert date to dec
    dec = float(date.strftime('%j'))/365.1
    year = float(date.strftime('%Y'))
    print(date,dec,year)
    mode.append(dec)
    time.append(year+dec)
    
    # extract temp
    #scale = float(bandmd["scale_factor"])
    #offset = float(bandmd["add_offset"])
    #data = band.ReadAsArray()*scale + offset
    data = band.ReadAsArray()
    
    #print(date, " -> ", np.mean(data), "m")
    #temp.append(data[-1,-1]) 
    # K to C
    temp[:,:,b-1] = data-273.15    
    #sys.exit()

time,dates,mode = np.array(time), np.array(dates),np.array(mode)
#temp = temp[4,3,:] 

# Initialisation
nfigure = 0
phase,amplitude = np.zeros((ds.RasterYSize,ds.RasterXSize)),np.zeros((ds.RasterYSize,ds.RasterXSize))
datemin, datemax= int(np.min(time)), int(np.max(time))+1
# print(datemin,datemax)
# sys.exit()
doplot='yes'

def seasonal(t,v,a,b,c):
    return v*t + a*np.cos(2*np.pi*t) + b*np.sin(2*np.pi*t) + c

def invers(temp,time):
    # fit
    G = np.zeros((len(temp),4))
    G[:,0] = time
    G[:,1] = np.cos(2*np.pi*(time))
    G[:,2] = np.sin(2*np.pi*(time))
    G[:,3] = 1

    pars = np.dot(np.dot(np.linalg.inv(np.dot(G.T,G)),G.T),temp)
    v,a,b,c = pars[0],pars[1],pars[2],pars[3]
    # fit = v*(time-datemin) + a*np.cos(2*np.pi*(time-datemin)) + b*np.sin(2*np.pi*(time-datemin)) + c
    return v, a, b, c 

def plot(dates,time,mode,temp):

    global nfigure, datemin, datemax

    v, a, b, c = invers(temp, time-datemin)
    amp,dphi=np.sqrt(a**2+b**2),np.arctan2(b,a)
    if dphi < 0:
        dphi = dphi + 2*np.pi
    phase[i,j], amplitude[i,j] = dphi,amp
    print('Mean T: {:0.3f} °'.format(np.mean(t)))
    print('Vitesse: {:0.3f} °/yr'.format(v))
    print('Amplitude: {:0.3f} °'.format(np.sqrt(a**2+b**2)))
    print('Phase Shift: {:0.3f} rad'.format(np.arctan2(b,a)))
    fit = seasonal(time-datemin,v,a,b,c)

    # plot data
    fig = plt.figure(nfigure,figsize=(12,5))
    ax = fig.add_subplot(1,1,1)
    # convert idates to num
    x = date2num(dates)
    xmin,xmax=datetime.datetime(int(datemin), 1, 1),datetime.datetime(datemax, 1, 1)
    xlim=date2num(np.array([xmin,xmax]))
    # format the ticks
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
    ax.plot(x,temp, '-bo',ms=.5,lw=.1,label='ERAI point {}, {}'.format(lat,lon),rasterized=True)
    ax.plot(x,fit,'-r',lw=2.0,label='MAAT: {0:.2f}, v: {1:.2f} deg/yrs, dphi:{2:.2f} rad, Amp:{3:.2f} deg'.format(np.mean(temp),v,dphi,amp))
    ax.plot(xlim,[0,0],'--',c='black')
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
    ax2.plot(mode,temp, '-bo',ms=.5,lw=.1,label='ERAI point {}, {}'.format(lat,lon),rasterized=True)
    x = np.arange(0,1,0.01)
    fit = seasonal(x,v,a,b,c)
    ax2.plot(x,fit,'-r',lw=2.0,label='MAAT: {0:.2f}'.format(np.mean(temp)))
    ax2.plot([0,1],[0,0],'--',c='black')

    idx = np.argwhere(np.diff(np.sign(fit))).flatten()
    for k in range(len(idx)):
        ax2.text(x[idx[k]],-25,'{:0.3f}'.format(x[idx[k]]))
        ax2.text(x[idx[k]],-28,(datetime.datetime(int(datemin), 1, 1) + datetime.timedelta(x[idx[k]]*365.1 - 1)).strftime('%Y/%m/%d'))
        ax2.plot([x[idx[k]],x[idx[k]]],[-30,30],'--',c='black')

    ax2.set_xticks(np.arange(0, 1, 1./12))
    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    plt.xlabel('Time')
    plt.ylabel('Temperature (deg)')
    plt.legend(loc='best')
    ax2.set_ylim([-30,30])
    fig2.savefig('Modtemperatures_{}_{}.pdf'.format(lat,lon), format='PDF')

for i in range(ds.RasterYSize):
    lat = ds.GetGeoTransform()[3]-0.125 + ds.GetGeoTransform()[5]*i
    for j in range(ds.RasterXSize):
        lon = ds.GetGeoTransform()[0]+0.125 + ds.GetGeoTransform()[1]*j
        print('lon: {}, lat: {}'.format(lon,lat))
        t = temp[i,j,:]
        if doplot=='yes':
            plot(dates,time,mode,t)
            nfigure+=2
            #plt.show()
            # plt.close()
            # sys.exit()

        # save data
        fid=open('temperatures_{}_{}.txt'.format(lat,lon),'wb')
        #np.savetxt(fid,np.vstack([time,t]).T,fmt='%.6f',delimiter='\t',newline='\n')
        np.savetxt(fid,np.vstack([dates,time,t]).T,fmt='%s',delimiter=',',newline='\n')

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

