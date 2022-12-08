#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Simon DAOUT (Oxford)
############################################

# import sys,os
# from sys import argv,stdin,stdout
# from os import environ,path
# import logging
# import getopt

from sys import argv,exit,stdin,stdout
import getopt
import os, math, sys
from os import path
import logging

import numpy as np
import scipy.optimize as opt
import scipy.linalg as lst
from osgeo import gdal
gdal.UseExceptions()

# plot
import matplotlib
if os.environ["TERM"].startswith("screen"):
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

################################
# CLASSES
################################

class network:
    """
    Load InSAR displacements and LOS angle maps
    name: name of raster displacement map (convention positive towards satellite)
    reduction: reducted name 
    wdir: relative path input files
    lookf: name incidence angle file (angle between vertical and LOS)
    headf: name heading file (angle between North and LOS)
    scale: scale LOS displacement map
    format: format input files: GTiff, NetCDF, ROI_PAC (default: GTiff)
    scale_sig: scale sigma file
    bounds: optional bounds for plot [losmin,losmax]
    """
    
    def __init__(self,name,wdir,reduction,lookf,headf,sigmaf=None,scale=1,scale_sig=1000,format='GTiff',bounds=None, ref_zone=None):
        self.name = name
        self.wdir = wdir
        self.reduction = reduction
        self.path = wdir + name
        self.lookf = wdir + lookf
        self.headf = wdir + headf
        if sigmaf is not None:
          self.sigmaf = wdir + sigmaf
        else:
          self.sigmaf = None
        self.scale = scale
        self.scale_sig = scale_sig
        self.bounds = bounds
        self.ref_zone = ref_zone

    def load(self,rot):
        logger.info('Read track: {}'.format(self.name))
        ds = gdal.Open(self.path,gdal.GA_ReadOnly)
        band = ds.GetRasterBand(1)
        self.los = self.scale*band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        band.FlushCache()
        self.nlines,self.ncols = ds.RasterYSize, ds.RasterXSize
        logger.info('Number of pixel: {}'.format(len(self.los.flatten())))
        logger.info('Nlines: {}, Ncols: {}'.format(self.nlines,self.ncols))
        del ds, band

        # compute losmax,losmin
        if self.bounds == None:
            self.losmax = np.nanpercentile(self.los,98)
            self.losmin = np.nanpercentile(self.los,2)
        else:
            self.losmin = self.bounds[0]
            self.losmax = self.bounds[1]
        
        ds = gdal.Open(self.lookf,gdal.GA_ReadOnly)
        # param output files
        self.gt = ds.GetGeoTransform()
        self.projref = ds.GetProjectionRef()
        self.driver = gdal.GetDriverByName('GTiff')

        band = ds.GetRasterBand(1)
        self.look = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        # self.look[np.isnan(self.los)] = float('NaN') 
        band.FlushCache()
        del ds, band

        if self.sigmaf is not None:
            self.sigmaf = self.wdir + self.sigmaf
            ds = gdal.Open(self.sigmaf,gdal.GA_ReadOnly)
            band = ds.GetRasterBand(1)
            self.sigma = self.scale_sig*band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
            # self.sigma[np.isnan(self.los)] = float('NaN') 
            band.FlushCache()
            del ds, band
        else:
            self.sigma = np.ones((self.nlines,self.ncols))

        # compute siglosmax,siglosmin
        self.sig_losmax = np.nanpercentile(self.sigma,98)
        self.sig_losmin = np.nanpercentile(self.sigma,2)

        ds = gdal.Open(self.headf,gdal.GA_ReadOnly)
        band = ds.GetRasterBand(1)
        self.head = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        # self.head[np.isnan(self.los)] = float('NaN') 
        band.FlushCache()
        del ds, band

        # convert head, look to angle phi, theta in rad
        # theta: angle between LOS and horizontal
        self.theta = np.deg2rad(90.-self.look)
        # phi: horizontal angle between LOS and comp1
        self.phi = np.deg2rad(-90-self.head)

        logger.info('Average LOOK:{0:.5f}, THETA:{1:.5f} angles'.\
            format(np.nanmean(self.look),np.nanmean(np.rad2deg(self.theta))))

        logger.info('Average HEADING:{0:.5f}, PHI:{1:.5f} angles'.\
            format(np.nanmean(self.head),np.nanmean(np.rad2deg(self.phi))))

        # compute proj ENcomp3
        self.proj=[np.cos(self.phi),
               np.sin(self.phi),
               np.sin(self.theta)
              ]
        logger.info('Average horizontal LOS projection to east, north, up: {0:.5f} {1:.5f} {2:.5f}'.\
            format(np.nanmean(self.proj[0]*np.cos(self.theta)),np.nanmean(self.proj[1]*np.cos(self.theta)),np.nanmean(self.proj[2])))

        # compute proj Shortening
        self.proj=[(np.cos(rot)*self.proj[0] - np.sin(rot)*self.proj[1])*np.cos(self.theta),
               (np.sin(rot)*self.proj[0] + np.cos(rot)*self.proj[1])*np.cos(self.theta),
               self.proj[2]
               ]

        logger.info('Average LOS projection to comp1, comp2, comp3: {0:.5f} {1:.5f} {2:.5f}'.\
            format(np.nanmean(self.proj[0]),np.nanmean(self.proj[1]),np.nanmean(self.proj[2])))
        print()


################################
# MAIN
################################

def usage():
  print('invers_los2comp.py infile.py [-v] [-h]')
  print('-v Verbose mode. Show more information about the processing')
  print('-h Show this screen')

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

# init logger 
if level == 'debug':
    logging.basicConfig(level=logging.DEBUG,\
        format='%(lineno)s -- %(levelname)s -- %(message)s')
    logger = logging.getLogger('invers_los2comp.log')
    logger.info('Initialise log file {0} in DEBUG mode'.format('invers_los2comp.log'))

else:
    logging.basicConfig(level=logging.INFO,\
        format='%(lineno)s -- %(levelname)s -- %(message)s')
    logger = logging.getLogger('invers_los2comp.log')
    logger.info('Initialise log file {0} in INFO mode. Use option -v for a DEBUG mode'.format('invers_los2comp.log'))

if 1==len(argv):
  usage()
  assert False, "no input file"
  logger.critical('No input file')
  exit()

if len(argv)>1:
  try:
    fname=argv[1]
    logger.info('Read input file {0}'.format(fname))
    exec(open(path.abspath(fname)).read())
  except Exception as e: 
    logger.critical('Problem in input file')
    logger.critical(e)
    print(network.__doc__)
    exit()

# rotation angle: angle between comp1 and East
rot = np.deg2rad(rotation)

# Load data
M = len(insar)
vmax=0; vmin=0; sigmax=0; sigmin=0
for i in range(M):
    insar[i].load(rot)
    if insar[i].losmax > vmax:
        vmax = insar[i].losmax
    if insar[i].losmin < vmin:
        vmin = insar[i].losmin
    if insar[i].sig_losmax > sigmax:
        sigmax = insar[i].sig_losmax
    if insar[i].sig_losmin < sigmin:
        sigmin = insar[i].sig_losmin

    # should be able to crop within same area in the future
    nlines, ncols = insar[0].nlines, insar[0].ncols
    # check size maps identical
    if (insar[i].nlines != nlines) or (insar[i].ncols != ncols):
        logger.critical('Size input maps are not idenical. Please crop your data to the same area. Exit!')
        exit()
    # get param output files
    gt = insar[0].gt
    projref = insar[0].projref
    driver = insar[0].driver

    # Apply reference zone shift...
    if insar[i].ref_zone is not None:
      if insar[i].ref_zone[1]<=insar[i].los.shape[0] and insar[i].ref_zone[3]<=insar[i].los.shape[1]:
        for i in range(M):
            ref_z = insar[i].los[insar[i].ref_zone[0]:insar[i].ref_zone[1],insar[i].ref_zone[2]:insar[i].ref_zone[3]]
            if np.isnan(np.nanmean(ref_z)) is True:
                logger.critical('Reference zone given contains only NaN values. Re-choose [y0, y1, x0, x1]. Exit!')
                exit()
            insar[i].los = insar[i].los - np.nanmean(ref_z)
      else:
        logger.critical('Reference zone given is outside image. Re-choose [y0, y1, x0, x1]. Exit!')
        exit()
    else:
      pass

# define invert components
comp = np.array(comp) - 1 
N = len(comp)
comp_name = []
for n in range(N):
    if int(comp[n]) == 0:
        name = 'East + {} deg (clockwise)'.format(np.rad2deg(rot))
    elif int(comp[n]) == 1:
        name = 'North + {} deg (clockwise)'.format(np.rad2deg(rot))
    elif int(comp[n]) == 2:
        name = 'Up '
    else:
        logger.critical('Error defined inverted component. Exit!')
        logger.critical('[east, north, up], comp = [1,2,3]')
        logger.critical('[east, north, up], comp = [1,2,3]')
        logger.critical('[east, up], comp = [1,3]')
        exit()
    logger.info('Invert components: {}'.format(name))
    comp_name.append(name)

# crop options
if 'crop' in locals():
    if crop is not None:
        logger.info('Crop time series data between lines {}-{} and cols {}-{}'.format(int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])))
        ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])
    else:
        ibeg,iend,jbeg,jend = 0,nlines,0,ncols
else:
    ibeg,iend,jbeg,jend = 0,nlines,0,ncols

################################
# plot DATA
################################

print()
logger.info('Plot DATA ....') 
try:
    from matplotlib.colors import LinearSegmentedColormap
    cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
    cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
except:
    cmap=cm.rainbow
cmap_r = cmap.reversed()

fig=plt.figure(0, figsize=(14,12))

for i in range(M):
    d = insar[i]
    # plot LOS
    ax = fig.add_subplot(2,M,i+1)
    cax = ax.imshow(d.los,cmap=cmap_r,vmax=vmax,vmin=vmin,interpolation=None)
    ax.set_title('{}'.format(d.reduction))
    # Add the patch to the Axes
    if d.ref_zone is not None:
        rect = patches.Rectangle((d.ref_zone[2],d.ref_zone[0]),(d.ref_zone[3]-d.ref_zone[2]),(d.ref_zone[1]-d.ref_zone[0]), linewidth=1, edgecolor='grey', facecolor='none')
        ax.add_patch(rect)
    # Colorbar
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    # plot SIGMA LOS
    ax = fig.add_subplot(2,M,i+1+M)
    cax = ax.imshow(d.sigma,cmap=cmap_r,vmax=sigmax,vmin=sigmin,interpolation=None)
    ax.set_title('SIGMA {}'.format(d.reduction))
    # Add the patch to the Axes
    if d.ref_zone is not None:
        rect = patches.Rectangle((d.ref_zone[2],d.ref_zone[0]),(d.ref_zone[3]-d.ref_zone[2]),(d.ref_zone[1]-d.ref_zone[0]), linewidth=1, edgecolor='grey', facecolor='none')
        ax.add_patch(rect)
    # Colorbar
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

# fig.tight_layout()
fig.savefig('data_decomposition_{}.png'.format(output),format='PNG',dpi=300)
plt.show()

################################
# INVERSION
################################

# Initialise output matrix
disp = np.ones((nlines,ncols,3))*float('NaN')
sdisp = np.ones((nlines,ncols,3))*float('NaN')

print()
logger.info('Inversion pixel by pixel....')
logger.info('Size invers matrix: {}x{}'.format(M,N))

# lcomp1-square 2 views 
# loop over each line and cols
for i in range(ibeg,iend):
    if i % 50 == 0:
        print('Processing line: {}'.format(i)) 
    for j in range(jbeg,jend):

        # initialise matrix for each pixels
        data = np.zeros((M))
        G = np.zeros((M,N))
        rms = np.zeros((M))

        for m in range(M):
            d = insar[m]
            # build data matrix 
            data[m] = d.los[i,j]
            rms[m] = d.sigma[i,j]

            for n in range(N):
                G[m,n] = d.proj[int(comp[n])][i,j]

        rms[np.isnan(rms)] = 1.
        # remove los with NaN
        # could be only one 
        index = np.flatnonzero(~np.isnan(data))
        data = data[index]
        rms = rms[index]
        G = G[index,:]

        # Inversion
        if len(data)>1:
            try:
                pars = lst.lstsq(G,data)[0]

                if iter>1:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,pars,fprime=_fprime,iter=iter,full_output=True,iprint=0,acc=acc)[0]
                
                for n in range((N)):        
                    disp[i,j,int(comp[n])] = pars[n]

                Cd = np.diag(rms**2, k = 0)
                try:
                    sigmam = np.sqrt(np.linalg.inv(np.dot(np.dot(G.T,np.linalg.inv(Cd)),G)))
                    for n in range((N)): 
                        sdisp[i,j,int(comp[n])] = np.abs(sigmam[n,n])
                except:
                    pass
                    # print('data:',data)
                    # print('G:',G)
                    # print('rms:',rms)
                    # print()

            except:
                pars = np.zeros((N))
                #print(i,j)
                #print('data:',data)
                #print('G:',G)
                #print()


################################
# PLOT RESULTS
################################

fig=plt.figure(1, figsize=(14,12))

for n in range(N):

    data = disp[:,:,int(comp[n])]
    data[data==0] = float('NaN')

    # plot comps
    vmax =  np.nanpercentile(data,92)
    ax = fig.add_subplot(2,N,n+1)
    cax = ax.imshow(data,cmap=cmap_r,vmax=vmax,vmin=-vmax,interpolation=None)
    ax.set_title('{}'.format(comp_name[n]))


    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    # plot SIGMA LOS
    sigdata = sdisp[:,:,int(comp[n])]
    # sigdata[sigdata==0] = float('NaN')
    # sigdata[sigdata==1] = float('NaN')

    vmax=np.nanpercentile(sigdata,95)
    ax = fig.add_subplot(2,N,n+1+N)
    cax = ax.imshow(sigdata,cmap=cmap_r,vmax=vmax,vmin=0,interpolation=None)
    ax.set_title('SIGMA {}'.format(comp_name[n]))

    # Colorbar
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

# fig.tight_layout()
fig.savefig('decomposition_{}.png'.format(output), format='PNG',dpi=300)

################################
# SAVE RESULTS
################################

# Save output files
for n in range(N):
    ds = driver.Create('U{}_{}.tif'.format(comp[n],output), ncols, nlines, 1, gdal.GDT_Float32)
    band = ds.GetRasterBand(1)
    band.WriteArray(disp[:,:,int(comp[n])])
    ds.SetGeoTransform(gt)
    ds.SetProjection(projref)
    band.FlushCache()

    ds = driver.Create('sig_U{}_{}.tif'.format(comp[n],output), ncols, nlines, 1, gdal.GDT_Float32)
    band = ds.GetRasterBand(1)
    band.WriteArray(sdisp[:,:,int(comp[n])])
    ds.SetGeoTransform(gt)
    ds.SetProjection(projref)
    band.FlushCache()

plt.show()
