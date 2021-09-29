#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################################
#
# PyGdalSAR: An InSAR post-processing package
# written in Python-Gdal
#
############################################
# Author        :   Simon Daout (ISTerre)
############################################

import numpy as np
import scipy as sp
import scipy.optimize as opt
import scipy.linalg as lst
from osgeo import gdal
import pyproj
gdal.UseExceptions()
from sys import argv,exit,stdin,stdout
import getopt
from os import path

import os, sys ,re
import geopandas as gpd
import pandas as pd
import shapely.speedups
from shapely.geometry import Point, Polygon
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
shapely.speedups.enable()

############################
# USAGE: convertGNSS_to_LOS.py input.py
# EXAMPLE INPUT FILE: input.py
###########################

## define path to data
#wdir = '/Users/symeone/work/italy/'

# define track name
#track = 'Ascending'

## define gnss file
#gps_file = wdir + 'gps/02_COMB_99UncSc_Filt_CLASS_ABC_Rout_100_Tsig_3.txt'

## define format incidence and heading file
## not same convention in GAMMA or ROIPAC format
#iformat = 'ROIPAC'

## InSAR data T004
#heading_file= wdir + '/insar/T117_head_s90.tiff'
#inc_file=wdir + '/insar/T117_look_s90.tiff'
#av_heading = 76.71
#av_inc = 32

#plot optional shapefile
#shapefile = '/gps/italian-maps-shapefiles/italy-with-regions/reg2011_g.shp'

## define frame coordinates
#LAT_REF1=38
#LAT_REF2=38
#LAT_REF3=46
#LAT_REF4=46
#LON_REF1=9.5
#LON_REF2=17
#LON_REF3=9.5
#LON_REF4=17

##########################
def usage():
  print('convertGNSS_to_LOS.py infile.py [-v] [-h]')
  print('-v Verbose mode. Show more information about the processing')
  print('-h Show this screen')

def load_gps(file):
    gps_gdf = gpd.GeoDataFrame()
    gps_gdf['geometry'] = None
    index = 0
    fl = open(file, "r").readlines()
    for line in fl:
        if not line.startswith(('#')):
            # Name lon lat Ve Vn Vup Se Sn Sup C
            try:
              sta, lon, lat, Ve, dVe, Vn, dVn ,Vup, dVup = line.split()
            except:
              print("Fail reading GNSS data block. Please organise input file as follow:")
              print("sta, lon, lat, Ve, dVe, Vn, dVn ,Vup, dVup")
              print("Exit.")
            gps_gdf.loc[index, 'geometry'] = Point(float(lon), float(lat))
            gps_gdf.loc[index, 'station'] = sta[:4]
            gps_gdf.loc[index, 've'] = float(Ve)  # eastern velocity in mm/yr in fixed eurasia reference frame
            gps_gdf.loc[index, 'vn'] = float(Vn)  # northern velocity in mm/yr in fixed eurasia reference frame
            gps_gdf.loc[index, 'vu'] = float(Vup)
            gps_gdf.loc[index, 'se'] = float(dVe)  # sigma ve
            gps_gdf.loc[index, 'sn'] = float(dVn)  # sigma vn
            gps_gdf.loc[index, 'su'] = float(dVup)  # sigma vup
            index += 1
    return gps_gdf

def neu2los_gamma(vn, ve, vu, theta, phi):
    los = - ve * np.cos(phi * np.pi / 180) * np.sin(theta * np.pi / 180) \
          + vn * np.sin(phi * np.pi / 180) * np.sin(theta * np.pi / 180) \
          + vu * np.cos(theta * np.pi / 180)
    return los

def neu2los_roipac(vn, ve, vu, look, head):
    theta = np.deg2rad(90.-look)
    phi = np.deg2rad(-90-head)
    los = ve * np.cos(phi) * np.cos(theta) \
            + vn * np.sin(phi) * np.cos(theta) \
            + vu *np.sin(theta)
    #print([np.cos(phi) * np.cos(theta),np.sin(phi) * np.cos(theta),np.sin(theta)])
    #sys.exit()
    return los

def add_inc_head_los(gps_df, look, heading):
    # extract pixel values from inc and heading based on gps lon lat and insert into geodataframe
    gps_df.loc[:, 'look'] = [look.extract_pixel_value(point.x, point.y, 0)[0] for point in gps_df['geometry']]
    gps_df.loc[:, 'heading'] = [heading.extract_pixel_value(point.x, point.y, 0)[0] for point in gps_df['geometry']]

    # calculate los per gps point
    if iformat == 'ROIPAC':
            gps_df.loc[:, 'los'] = [neu2los_roipac(vn, ve, vu, look, heading) for vn, ve, vu, look, heading
                            in gps_df[['vn', 've', 'vu', 'look', 'heading']].to_numpy()]
            gps_df.loc[:, 'siglos'] = [neu2los_roipac(vn, ve, vu, look, heading) for vn, ve, vu, look, heading
                            in gps_df[['sn', 'se', 'su', 'look', 'heading']].to_numpy()]
    elif iformat == 'GAMMA':
            gps_df.loc[:, 'los'] = [neu2los_gamma(vn, ve, vu, look, heading) for vn, ve, vu, look, heading
                            in gps_df[['vn', 've', 'vu', 'look', 'heading']].to_numpy()]
            gps_df.loc[:, 'siglos'] = [neu2los_gamma(vn, ve, vu, look, heading) for vn, ve, vu, look, heading
                            in gps_df[['sn', 'se', 'su', 'look', 'heading']].to_numpy()]

def plot_gps_distribution(gps, poly):
    fig, ax = plt.subplots(1,figsize=(8, 9))
    x,y = poly.exterior.xy
    ax.plot(x, y, color='#6699cc', alpha=0.7,
    linewidth=3, solid_capstyle='round', zorder=2)
    gps.plot(ax=ax, facecolor='orange')
    ax.xaxis.set_visible(True)
    ax.yaxis.set_visible(True)
    #plt.show()
    fig.savefig(track+'_gps_frame.png', format='PNG', dpi=150, bbox_inches='tight')

def plot_gps_in_LOS(gps, poly):
    fig = plt.figure(2,figsize=(10,8))
    fig.subplots_adjust(wspace=0.001)
    ax = fig.add_subplot(1,1,1)

    x,y = poly.exterior.xy
    minx,maxx,miny,maxy = np.min(x),np.max(x),np.min(y),np.max(y)
    ax.axis((minx,maxx,miny,maxy))

    # plot contextfile
    try:
      import contextily as ctx    
      ctx.add_basemap(ax,crs=4326)        
    except:
      print('Contextily package not installed. Skip backgroup topography plot')

    # plot world
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    if shapefile: 
      bounds = gpd.read_file(shapefile)
      bounds = bounds.to_crs(world.crs)
      bounds.plot(ax=ax,facecolor='none', color='none',edgecolor='black',zorder=1)
    else:
      world.plot(ax=ax,facecolor='none',color='none', edgecolor='black',zorder=1)

    los = gps['los'].to_numpy()
    vmin=-4
    vmax=4
    #vmin = np.nanpercentile(los, 8)
    #vmax = np.nanpercentile(los, 92)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    try:
        from matplotlib.colors import LinearSegmentedColormap
        cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
        cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
        cmap = cmap.reversed()
    except:
        cmap=cm.rainbow   
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    facelos = m.to_rgba(los) 
    gps.plot(ax=ax, facecolor=facelos, zorder=3)
    
    ax.set_xticks(np.linspace(minx,maxx,3))
    ax.set_yticks(np.linspace(miny,maxy,3))
    
    ax.set_xlim([minx,maxx]) 
    ax.set_ylim([miny,maxy])  

    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(m, cax=c)
    cbar.set_label('LOS Velocity (mm/yr)', rotation=270)
    ax.title.set_text('{}'.format(track))
    fig.savefig(track+'_gps_LOS.png', format='PNG', dpi=150, bbox_inches='tight')
    fig.savefig(track+'_gps_LOS.pdf', format='PDF', dpi=150)

class OpenTif(object):
    """ a Class that stores the band array and metadata of a Gtiff file."""
    def __init__(self, filename, sigfile=None):
        self.ds = gdal.Open(filename)
        self.basename = os.path.splitext(os.path.basename(filename))[0]
        self.band = self.ds.GetRasterBand(1)
        self.data = self.band.ReadAsArray()
        self.xsize = self.ds.RasterXSize
        self.ysize = self.ds.RasterYSize
        self.left = self.ds.GetGeoTransform()[0]
        self.top = self.ds.GetGeoTransform()[3]
        self.xres = self.ds.GetGeoTransform()[1]
        self.yres = self.ds.GetGeoTransform()[5]
        self.right = self.left + self.xsize * self.xres
        self.bottom = self.top + self.ysize * self.yres
        self.projection = self.ds.GetProjection()
        pix_lin, pix_col = np.indices((self.ds.RasterYSize,self.ds.RasterXSize))
        self.lat,self.lon = self.top + self.yres*pix_lin, self.left+self.xres*pix_col

        # convert 0 and 255 to NaN
        self.data[self.data==0.] = float('NaN')
        self.data[self.data==255] = float('NaN')

    def extract_pixel_value(self, lon, lat, n):

        x = int((lon-self.left)/self.xres+0.5)
        y = int((lat - self.top) / self.yres + 0.5)

        pixel_values = self.data[y - n: y + n + 1, x - n: x + n + 1]

        index = np.nonzero(~np.isnan(pixel_values))
        if len(index[0]) > 0:
            m = np.nanmean(pixel_values)
            std = np.nanstd(pixel_values)
        else:
            m = float('NaN')
            std = float('NaN')

        if m == 0:  #if only NaN nanmean is 0
            m = float('NaN')
        if std == 0:
            std = float('NaN')

        return m, std

class CreateTif(object):
    """ a Class that stores the band array and metadata of a Gtiff file."""
    def __init__(self, data, poly):
        x,y = poly.exterior.xy
        self.left,self.top = x[2],y[2]
        self.xres = 0.0008333333333333537
        self.yres = 0.0008333333333333537
        self.right,self.bottom = x[0],y[0]
        self.xsize = np.int((self.right - self.left)/self.xres) + 1
        self.ysize = np.int((self.top - self.bottom)/self.yres) + 1
        pix_lin, pix_col = np.indices((self.ysize,self.xsize))
        self.lat,self.lon = self.top + self.yres*pix_lin, self.left+self.xres*pix_col

        # set average value to band
        self.data = np.ones((self.xsize,self.ysize))
        self.data.fill(data)

    def extract_pixel_value(self, lon, lat, n):

        x = int((lon-self.left)/self.xres+0.5)
        y = int((lat - self.top) / self.yres + 0.5)

        pixel_values = self.data[y - n: y + n + 1, x - n: x + n + 1]

        index = np.nonzero(~np.isnan(pixel_values))
        if len(index[0]) > 0:
            m = np.nanmean(pixel_values[index]) 
            std = np.nanstd(pixel_values)
        else:
            m = float('NaN')
            std = float('NaN')

        if m == 0:  #if only NaN nanmean is 0
            m = float('NaN')
        if std == 0:
            std = float('NaN')

        return m, std

if __name__ == "__main__":
    try:
        opts,args = getopt.getopt(argv[1:], "h", ["help"])
    except:
        print(str(err))
        print("for help use --help")
        exit()

    level = 'basic'
    for o in argv:
        if o in ("-h","--help"):
            usage()
            exit()
        if o in ("-v","--verbose"):
            level = 'debug'

    if len(argv)>1:
        try:
            fname=argv[1]
            print('Read input file {0}'.format(fname))
            try:
                sys.path.append(path.dirname(path.abspath(fname)))
                exec ("from "+path.basename(fname)+" import *")
            except:
                exec(open(path.abspath(fname)).read())

        except Exception as e:
            print('Problem in input file')
            print(e)
            print(network.__doc__)
            exit()

    # define poly
    try:
        epsi
    except NameError:
        epsi = 0
    try:
        LON_REF2;LON_REF1;LON_REF3; LON_REF4
    except NameError:
        print('LAT_REF* and LON_REF* not defined. Exit!')
        sys.exit()
    coords = [(LON_REF2-epsi, LAT_REF2-epsi), (LON_REF1-epsi,LAT_REF1-epsi), (LON_REF3-epsi,LAT_REF3-epsi), (LON_REF4-epsi,LAT_REF4-epsi)]
    poly = Polygon(coords)
    
    try:
        gps_file
    except NameError:
        print('gps_file not defined. Exit!')
        sys.exit()
    # load gnss
    all_gps = load_gps(gps_file)

    # gps within poly
    gps = all_gps[all_gps.within(poly)]
    try:
        print('Read incidence and heading file in {} format'.format(iformat))
    except:
        iformat='ROIPAC'
        print('Read incidence and heading file in {} format'.format(iformat))
    
    try:
        inc_file; heading_file 
        inc = OpenTif(inc_file)
        heading = OpenTif(heading_file)
    except NameError:
        print('inc_file or heading_file or insar_file not defined. Read av_inc and av_heading')
        try: 
           av_heading; av_inc 
           inc = CreateTif(av_inc,poly)      
           heading = CreateTif(av_heading,poly)      
        except:
           print('av_inc and av_heading not defined. Exit!')
           sys.exit()
    
    # extract track numbered
    try:
      track
      track = str(track)
    except NameError:
      print('track Name not defined. Set to T000')
      track = 'T000'

    # compute inc angle gps
    add_inc_head_los(gps, inc, heading)
    print(gps)

    # if heading, then increase window size
    for n in range(4,200,2):
        #print(n)
        for index,g in gps.iterrows():
            if np.isnan(g['heading']):

                gps.loc[index, 'look'] = inc.extract_pixel_value(g['geometry'].x, g['geometry'].y, n)[0]
                gps.loc[index, 'heading'] = heading.extract_pixel_value(g['geometry'].x, g['geometry'].y, n)[0]
                if iformat == 'ROIPAC':
                    gps.loc[index, 'los'] = neu2los_roipac(g['vn'], g['ve'], g['vu'], gps.loc[index,'look'], gps.loc[index,'heading'])
                    gps.loc[index, 'siglos'] = np.abs(neu2los_roipac(g['sn'], g['se'], g['su'], gps.loc[index,'look'], gps.loc[index,'heading']))
                elif iformat == 'GAMMA':
                    gps.loc[index, 'los'] = neu2los_gamma(g['vn'], g['ve'], g['vu'], gps.loc[index,'look'], gps.loc[index,'heading'])
                    gps.loc[index, 'siglos'] = np.abs(neu2los_gamma(g['sn'], g['se'], g['su'], gps.loc[index,'look'], gps.loc[index,'heading']))

    if gps.los.count() == 0:
        print('No GPS within polygon. You may increase epsi value. Exit!')
        sys.exit()

    # remove NaNs
    gps = gps.dropna()
    print(gps)

    # Plot frame and GPS
    #plot_gps_distribution(gps, poly)
    plot_gps_in_LOS(gps, poly)
    
    # save table
    gps['lon2'] = None; gps['lat2'] = None
    for index,g in gps.iterrows():
      gps.loc[index,'lon2'] = g['geometry'].x
      gps.loc[index,'lat2'] = g['geometry'].y
    gps.to_csv(path_or_buf=track+'_gps_table.txt', sep=' ', index=False, header=None)    
    gps.to_csv(path_or_buf=track+'_gps_table.txt', sep=' ', index=False, header='[geometry ,sta, Ve, dVe, Vn, dVn ,Vup, dVup, look, heading, los, siglos, lon, lat]')    
    #plt.show()
