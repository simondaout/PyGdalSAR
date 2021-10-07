#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################################
#
# PyGdalSAR: An InSAR post-processing package
# written in Python-Gdal
#
############################################
# Author        :   Qi Ou (Oxford)
############################################

import numpy as np
import scipy as sp
import scipy.optimize as opt
import scipy.linalg as lst
from osgeo import gdal
import pyproj
#gdal.UseExceptions()
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

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
shapely.speedups.enable()

##################################
# EXAMPLE INPUT FILE 
##################################

## define projection
#UTM = pyproj.Proj("+init=EPSG:32646")

## define path to data
#wdir = '/home/cometraid14/daouts/work/tibet/qinghai/'

## define gnss file
#gps_file = wdir + 'data/gps/table_liang_eurasia_3d_ll.dat'

## define format incidence and heading file
## not same convention in GAMMA or ROIPAC format
#iformat = 'ROIPAC'

## InSAR data T004
#insar_file = wdir + 'data/insar/T004_inter/T004_inter_LOSVelocity_nan_mmyr_s90.tiff'
#insar_sig_file = wdir + 'data/insar/T004_inter/T004_LOS_sigVelocity.tiff'
#heading_file= wdir + 'data/insar/T004_inter/T004_head_s90.tiff'
#inc_file=wdir + 'data/insar/T004_inter/T004_look_s90.tiff'

## define frame coordinates
#LAT_REF1=38.6042
#LAT_REF2=38.8860

##################################

def usage():
  print('tie_to_gps.py infile.py [-v] [-h]')
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
    fig, ax = plt.subplots(figsize=(3, 2))
    x,y = poly.exterior.xy
    ax.plot(x, y, color='#6699cc', alpha=0.7,
    linewidth=3, solid_capstyle='round', zorder=2)
    gps.plot(ax=ax, facecolor='orange')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    #ax.set_title(poly.loc[0, 'Name']) #[:4][:4] + '_GPS'
    #plt.show()
    fig.savefig(track+'gps_frame.png', format='PNG', dpi=150, bbox_inches='tight')

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
        
        if sigfile is not None:
            self.dst = gdal.Open(sigfile)
            self.bandt = self.dst.GetRasterBand(1)
            self.sigma = self.bandt.ReadAsArray()
            if self.dst.RasterXSize != self.xsize or self.dst.RasterYSize != self.ysize:
                print('Error: Sigma and Velocity file not the same size! Set sigma to 1.')
                self.sigma = np.ones((self.ysize,self.xsize))
        else:
            self.sigma = np.ones((self.ysize,self.xsize))

    def extract_pixel_value(self, lon, lat, n):
       
        x = int((lon-self.left)/self.xres+0.5)
        y = int((lat - self.top) / self.yres + 0.5)

        pixel_values = self.data[y - n: y + n + 1, x - n: x + n + 1]
        pixel_sigma = self.sigma[y - n: y + n + 1, x - n: x + n + 1]

        index = np.nonzero(~np.isnan(pixel_values))
        if len(index[0]) > 0:
            #m1 = np.nanmean(pixel_values)
            m = np.nansum(pixel_values[index] * pixel_sigma[index]) / np.nansum(pixel_sigma[index])
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
     # load insar
    # read format incidence heading angle
    # not same convention in GAMMA or ROIPAC
    try:
        print('Read incidence and heading file in {} format'.format(iformat))
    except:
        iformat='ROIPAC'
        print('Read incidence and heading file in {} format'.format(iformat))
   
    try:
        inc_file; heading_file; insar_file
    except NameError:
         print('inc_file or heading_file or insar_file not defined. Exit!')
         sys.exit()
    inc = OpenTif(inc_file)
    heading = OpenTif(heading_file)

    try:
        print('Load InSAR file:', insar_file)
        insar = OpenTif(insar_file,insar_sig_file)
        print('Load uncertainty file:', insar_sig_file)
    except:
        print('Load InSAR file:', insar_file)
        insar = OpenTif(insar_file)

    # extract track numbered
    m =re.search('\w+(?<=_)', insar_file)
    track = m.group(0)

    # compute inc angle gps
    add_inc_head_los(gps, inc, heading)
    #print(gps)
    
    
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

    # compute diff GPS - InSAR
    gps.loc[:, 'diff'] = [gps_los -  insar.extract_pixel_value(point.x, point.y, 2)[0] for point,gps_los in gps[['geometry','los']].to_numpy()]
    gps.loc[:, 'std'] = [insar.extract_pixel_value(point.x, point.y, 2)[1] for point in gps['geometry'].to_numpy()]
   

    # if GPS_LOS - LOS is NaN, then increase window size
    for n in range(4,200,2):
        for index,g in gps.iterrows():
            if np.isnan(g['diff']) or np.isnan(g['std']):
                gps.loc[index,'diff'] = g['los'] - insar.extract_pixel_value(g['geometry'].x, g['geometry'].y, n)[0]
                gps.loc[index,'std'] = insar.extract_pixel_value(g['geometry'].x, g['geometry'].y, n)[1]

    # remove NaNs
    gps = gps.dropna()
    print(gps)
    
    # Plot frame and GPS
    plot_gps_distribution(gps, poly)

    # data vector
    d = gps['diff'].to_numpy()
    sig = (gps['std'].to_numpy() + gps['siglos'].to_numpy()) / 2

    # convert GPS lat/lon to km
    x =  [ point.x for point in gps['geometry'] ]
    y =  [ point.y for point in gps['geometry'] ]
    east, north = UTM(x,y)
    #east, north = x, y

    # invers G matrix
    G = np.zeros((len(d),5))
    G[:,0] = np.array(north)**2
    G[:,1] = np.array(east)**2
    G[:,2] = north
    G[:,3] = east
    G[:,4] = 1
    
    # inversion
    x0 = lst.lstsq(G,d)[0]
    _func = lambda x: np.sum(((np.dot(G,x)-d)/sig)**2)
    _fprime = lambda x: 2*np.dot(G.T/sig, (np.dot(G,x)-d)/sig)
    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=20000,full_output=True,iprint=0)[0]
    #print(pars)
    
    # compute residual between data and model
    res = d - np.dot(G,pars).flatten()
    var = np.nanstd(res)
    print('Mean Residual:', np.mean(res))
    print('Variance:', var)
    
    # Build G matrix for all insar points
    east, north = UTM(insar.lon.flatten(), insar.lat.flatten())
    #east, north = insar.lon.flatten(), insar.lat.flatten()
    G = np.zeros((len(insar.data.flatten()),5))
    G[:,0] = np.array(north)**2
    G[:,1] = np.array(east)**2
    G[:,2] = north
    G[:,3] = east
    G[:,4] = 1

    ramp_array = np.dot(G, pars).reshape(insar.ysize, insar.xsize)
    model_array = ramp_array + insar.data
    
    # remove cst
    cst = np.nanmean(model_array)
    model_array = model_array - cst

    # Plot
    vmin = np.nanpercentile(insar.data, 2)
    vmax = np.nanpercentile(insar.data, 98)
    vmax_insar = np.max([np.abs(vmin),np.abs(vmax)])
    vmin_insar = -vmax_insar
    vmin_insar = -10
    vmax_insar = 10

    try:
        from matplotlib.colors import LinearSegmentedColormap
        cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
        cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
        cmap = cmap.reversed()
    except:
        cmap=cm.rainbow
    
    fig, (ax1,ax2,ax3) = plt.subplots(1, 3,  figsize=(10,4))

    cax = ax1.imshow(insar.data, cmap=cmap, vmin=vmin_insar, vmax=vmax_insar, interpolation='nearest')
    ax1.set_title('InSAR')
    plt.setp( ax1.get_xticklabels(), visible=None)
    plt.setp( ax1.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax1)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)
    c.set_label('mm/yr')

    cax = ax2.imshow(ramp_array, cmap=cmap, vmin=vmin_insar, vmax=vmax_insar,interpolation='nearest')
    ax2.set_title('Ramp')
    plt.setp( ax2.get_xticklabels(), visible=None)
    plt.setp( ax2.get_yticklabels(), visible=None)
    
    # Plot
    #vmin = np.nanpercentile(model_array, 2)
    #vmax = np.nanpercentile(model_array, 98)
    #vmax_insar = np.max([np.abs(vmin),np.abs(vmax)])
    #vmin_insar = -vmax_insar

    cax = ax3.imshow(model_array, cmap=cmap, vmin=vmin_insar, vmax=vmax_insar,interpolation='nearest')
    ax3.set_title('InSAR - Ramp')
    plt.setp( ax3.get_xticklabels(), visible=None)
    plt.setp( ax3.get_yticklabels(), visible=None)

    divider = make_axes_locatable(ax3)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)
    c.set_label('mm/yr')
    
    fig.savefig(track +'tie_to_gps.png', format='PNG', dpi=150, bbox_inches='tight')

    # save output file
    # Export merged data to tif format.
    driver = gdal.GetDriverByName("GTiff")
    outdata = driver.Create(track+'projected.tif', insar.xsize, insar.ysize, 1, gdal.GDT_Float32)
    outdata.SetGeoTransform([insar.left, insar.xres, 0, insar.top, 0, insar.yres])  ##sets same geotransform as input
    outdata.SetProjection(insar.projection)  ##sets same projection as input
    outdata.GetRasterBand(1).WriteArray(model_array)
    outdata.FlushCache()
    outdata.FlushCache()

    plt.show()
    

