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
import gdal, pyproj
gdal.UseExceptions()

import os,sys,re
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

# define projection
UTM = pyproj.Proj("+init=EPSG:32646")

# define gnss file
gps_file = '../gps/table_liang_eurasia_3d_ll.dat'

# InSAR data
insar_file = '../insar/T004_inter/T004_inter_LOSVelocity_nan_mmyr_s90.tiff'
heading_file='../insar/T004_inter/T004_head_s90.tiff'
inc_file='../insar/T004_inter/T004_look_s90.tiff'
# define frame coordinates
LAT_REF1=38.6042
LAT_REF2=38.8860
LAT_REF3=35.9663
LAT_REF4=36.2203
LON_REF1=98.6828
LON_REF2=95.7729
LON_REF3=98.0398
LON_REF4=95.2283
epsi=0.1
#epsi=0.

def load_gps(file):
    gps_gdf = gpd.GeoDataFrame()
    gps_gdf['geometry'] = None
    index = 0
    fl = open(file, "r").readlines()
    for line in fl:
        if not line.startswith(('#')):
            # lon lat Ve Vn Vup Se Sn Sup C Name
            lon, lat, Ve, Vn, Vup, dVe, dVn, dVup, Cen, sta = line.split()
            gps_gdf.loc[index, 'geometry'] = Point(float(lon), float(lat))
            gps_gdf.loc[index, 'station'] = sta[:4]
            gps_gdf.loc[index, 've'] = float(Ve)  # eastern velocity in mm/yr in fixed eurasia reference frame
            gps_gdf.loc[index, 'vn'] = float(Vn)  # northern velocity in mm/yr in fixed eurasia reference frame
            gps_gdf.loc[index, 'vu'] = float(Vup)
            gps_gdf.loc[index, 'se'] = float(dVe)  # sigma ve
            gps_gdf.loc[index, 'sn'] = float(dVn)  # sigma vn
            gps_gdf.loc[index, 'su'] = float(dVup)  # sigma vup
            gps_gdf.loc[index, 'cen'] = float(Cen)  # correlation between east and northern velocities
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
    return los

def add_inc_head_los(gps_df, look, heading):
    # extract pixel values from inc and heading based on gps lon lat and insert into geodataframe
    gps_df.loc[:, 'look'] = [look.extract_pixel_value(point.x, point.y, 0)[0] for point in gps_df['geometry']]
    gps_df.loc[:, 'heading'] = [heading.extract_pixel_value(point.x, point.y, 0)[0] for point in gps_df['geometry']]
    # calculate los per gps point
    gps_df.loc[:, 'los'] = [neu2los_roipac(vn, ve, vu, look, heading) for vn, ve, vu, look, heading
                            in gps_df[['vn', 've', 'vu', 'look', 'heading']].to_numpy()]

def plot_gps_distribution(gps, poly):
    fig, ax = plt.subplots(figsize=(3, 2))
    x,y = poly.exterior.xy
    ax.plot(x, y, color='#6699cc', alpha=0.7,
    linewidth=3, solid_capstyle='round', zorder=2)
    gps.plot(ax=ax, facecolor='orange')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    #ax.set_title(poly.loc[0, 'Name']) #[:4][:4] + '_GPS'
    plt.show()
    fig.savefig('../gps+frame.png', format='PNG', dpi=150, bbox_inches='tight')

class OpenTif(object):
    """ a Class that stores the band array and metadata of a Gtiff file."""
    def __init__(self, filename):
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

    def extract_pixel_value(self, lon, lat, n):
        x = int((lon-self.left)/self.xres+0.5)
        y = int((lat - self.top) / self.yres + 0.5)
        #print()
        #print(self.xsize,self.ysize)
        #print(x,y)
        #if y > self.ysize:
        #    y = self.ysize-1
        #if y < 0:
        #    y = 0
        #if x > self.xsize:
        #   x = self.xsize-1
        #if x < 0:
        #    x = 0
        pixel_values = self.data[y - n: y + n + 1, x - n: x + n + 1]
        return np.nanmean(pixel_values), np.nanstd(pixel_values)

if __name__ == "__main__":

    # define poly
    coords = [(LON_REF2-epsi, LAT_REF2-epsi), (LON_REF1-epsi,LAT_REF1-epsi), (LON_REF3-epsi,LAT_REF3-epsi), (LON_REF4-epsi,LAT_REF4-epsi)]
    poly = Polygon(coords)

    # load gnss
    all_gps = load_gps(gps_file)
    #print(all_gps)

    # gps within poly
    gps = all_gps[all_gps.within(poly)]
    #print(gps)

    # load insar
    inc = OpenTif(inc_file)
    heading = OpenTif(heading_file)
    insar = OpenTif(insar_file)
    
    # compute inc angle gps
    add_inc_head_los(gps, inc, heading)
    #print(gps)

    # Plot frame and GPS
    #plot_gps_distribution(gps, poly)

    # compute diff GPS - InSAR
    gps.loc[:, 'diff'] = [gps_los -  insar.extract_pixel_value(point.x, point.y, 2)[0] for point,gps_los in gps[['geometry','los']].to_numpy()]
    gps.loc[:, 'std'] = [gps_los -  insar.extract_pixel_value(point.x, point.y, 2)[1] for point,gps_los in gps[['geometry','los']].to_numpy()]
    #print(gps['diff'].to_numpy())
    print(gps)

    # data vector
    d = gps['diff'].to_numpy()
    sig = gps['std'].to_numpy()

    # convert GPS lat/lon to km
    x =  [ point.x for point in gps['geometry'] ]
    y =  [ point.y for point in gps['geometry'] ]
    east, north = UTM(x,y)

    # invers G matrix
    G = np.zeros((len(d),3))
    G[:,0] = north
    G[:,1] = east
    G[:,2] = 1
    
    # inversion
    x0 = lst.lstsq(G,d)[0]
    _func = lambda x: np.sum(((np.dot(G,x)-d)/sig)**2)
    _fprime = lambda x: 2*np.dot(G.T/sig, (np.dot(G,x)-d)/sig)
    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=20000,full_output=True,iprint=0)[0]
    print(pars)
    
    # compute residual between data and model
    res = d - np.dot(G,pars).flatten()
    var = np.nanstd(res)
    print('Mean Residual:', np.mean(res))
    print('Variance:', var)
    
    # Build G matrix for all insar points
    east, north = UTM(insar.lon.flatten(), insar.lat.flatten())
    G = np.zeros((len(insar.data.flatten()),3))
    G[:,0] = north
    G[:,1] = east
    G[:,2] = 1

    ramp_array = np.dot(G, pars).reshape(insar.ysize, insar.xsize)
    model_array = ramp_array + insar.data
    
    # Plot
    vmin_insar = np.nanpercentile([insar.data,model_array], 2)
    vmax_insar = np.nanpercentile([insar.data,model_array], 98)
    
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

    cax = ax2.imshow(ramp_array, cmap=cmap, vmin=vmin_insar, vmax=vmax_insar,interpolation='nearest')
    ax2.set_title('Ramp')
    plt.setp( ax2.get_xticklabels(), visible=None)
    plt.setp( ax2.get_yticklabels(), visible=None)

    cax = ax3.imshow(model_array, cmap=cmap, vmin=vmin_insar, vmax=vmax_insar,interpolation='nearest')
    ax3.set_title('InSAR - Ramp')
    plt.setp( ax3.get_xticklabels(), visible=None)
    plt.setp( ax3.get_yticklabels(), visible=None)

    divider = make_axes_locatable(ax3)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)
    c.set_label('mm/yr')
    
    plt.show()
    

