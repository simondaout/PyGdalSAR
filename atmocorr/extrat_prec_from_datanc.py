#!/usr/bin/env python
# -*- coding: utf-8 -*-

import datetime, sys
import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
from pylab import date2num
import matplotlib.dates as mdates

list_dates = ['2020332','2020334','2020334','2020335'] 
time = []
new_date= []
precip = []

for date in list_dates:
    print("Read {}.nc file".format(date))
    # open .nc and read precipitations 
    ds = gdal.Open("NETCDF:"+date+".nc:precipitation", gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    bandmd = band.GetMetadata()
    print('Grid:', ds.RasterYSize,ds.RasterXSize)
    print('Metadata:', bandmd)
    data = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
   
    # extract lat,lon .nc
    ds_geo=ds.GetGeoTransform() 
    print('Coordinates:', ds.GetGeoTransform())
    pix_az, pix_rg = np.indices((ds.RasterYSize,ds.RasterXSize))
    lat,lon = ds_geo[3]+ds_geo[5]*pix_az, ds_geo[0]+ds_geo[1]*pix_rg
    
    # extrat time
    date = str(date)
    date_ref = datetime.datetime(int(date[:4]), 1, 1) 
    d = date_ref + datetime.timedelta(days=int(date[4:]))
    time.append(d)
    new_date.append(d.strftime("%d/%m/%y")) 
    
    # crop data within defined area and take mean
    lat_min,lat_max,lon_min,lon_max = 27., 28., 83., 84
    m = np.nanmean(data[np.logical_and(lon<lon_max,np.logical_and(lon>lon_min,np.logical_and(lat<lat_max,lat>lat_min)))])
    print(d, " -> ", m , bandmd["units"])
    precip.append(m)
    print()

# save data in a text file
list_dates, new_date, precip = np.array(list_dates),np.array(new_date), np.array(precip)
np.savetxt('precipitations.txt', np.vstack([new_date,precip]).T, header='times  |   precipitations', fmt="%s") 

# Plot
fig = plt.figure(0,figsize=(12,5))
ax = fig.add_subplot(1,1,1)
x = date2num(time)
ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
ax.plot(x,precip, '-o',ms=5.,lw=.5,c='blue',label='Mean Preciptation between Lats: {}-{} Lons: {}-{}'.format(lat_min,lat_max,lon_min,lon_max))
fig.autofmt_xdate()
plt.xlabel('Time (Year/month/day)')
plt.ylabel('Precipitation ({})'.format(bandmd["units"]))
plt.legend(loc='best')
plt.show()
