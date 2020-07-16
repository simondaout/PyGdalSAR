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
import os
import geopandas as gpd
import pandas as pd
import shapely.speedups
import re
from shapely.geometry import Point, Polygon
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
shapely.speedups.enable()

gps_file = '../gps/table_liang_eurasia_3d_ll.dat'
kml_file = '../insar/T004/T004_inter_LOSVelocity_nan_mmyr_s90.kml'

def load_gps(file):

    dated,east,north,down,esigma,nsigma,dsigma=\
                np.loadtxt(station,comments='#',usecols=(0,1,2,3,4,5,6), unpack=True,\
                    dtype='f,f,f,f,f,f,f')

    gps_gdf = gpd.GeoDataFrame()
    gps_gdf['geometry'] = None
    index = 0
    fl = open(file, "r").readlines()
    for line in fl:
        if not line.startswith(('*', "Table")):
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


if __name__ == "__main__":

    # load gnss
    # all_gps = load_gps(gps_file)
    # print(all_gps)

    # load frame kml into geopandas dataframe
    frames = gpd.read_file(kml_file, driver='KML')
    frames['track'] = frames['Name'].str[:4]
    tracks = frames.dissolve(by='track', aggfunc='sum')
    print(tracks)









