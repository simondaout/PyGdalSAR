#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Authors :   Hugo Watine, Louise Van Inghelandt, Simon Doaut
# Date:       03/02/2022
############################################

from sys import argv,exit,stdin,stdout
import getopt
import os, math
from os import path
import logging

import numpy as np
import scipy.ndimage
import scipy.optimize as opt
import numpy.linalg as lst
from osgeo import gdal, osr
gdal.UseExceptions()
from math import sin, cos
from numpy.lib.stride_tricks import as_strided

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

print()
print()
print('Author: Simon DAOUT')
print('Please cite:')
print('Daout, S., DAgostino, N., Pathier, E., Socquet, A., Lavé, J., Doin, M. P., ... & Benedetti, L. Along-Strike Variation of the Strain Partitioning within the Apennines as Seen from Large-Scale Multi-Temporal Insar Analysis. Available at SSRN 4429391.')
print()
print()

################################
# CLASSES
################################

class network:
    """
    Load InSAR displacements and LOS angle maps
    :name: name of raster displacement map (convention positive towards satellite)
    :reduction: reducted name 
    :wdir: relative path input files
    :incidence_file: name incidence angle file (angle between Vertical and LOS - ROIPAC convention
    :heading_file: name heading file (angle between North and LOS - ROIPAC convention)
    :scale: scale LOS displacement map
    :format: format input files: GTiff, NetCDF, ROI_PAC (default: GTiff)
    :scale_sig: scale sigma file
    :bounds: optional bounds for plot [losmin,losmax]
    :DEM: name DEM file (.dem format)
    """
    
    def __init__(self,name,wdir,reduction,incidence_file,heading_file,sigmaf=None,scale=1,scale_sig=1000,format='GTiff',bounds=None, ref_zone=None):
        self.name = name
        self.wdir = wdir
        self.reduction = reduction
        self.path = wdir + name
        self.incidence_file = wdir + incidence_file
        self.heading_file = wdir + heading_file
        if sigmaf is not None:
            self.sigmaf = wdir + sigmaf
        else:
            self.sigmaf = None
        self.scale = scale
        self.scale_sig = scale_sig
        self.bounds = bounds
        self.ref_zone = ref_zone
 
    def load(self,rot,slope):

        self.rot = rot
        self.slope = slope

        logger.info('Read track: {}'.format(self.name))
        ds = gdal.Open(self.path,gdal.GA_ReadOnly)
        band = ds.GetRasterBand(1)
        self.los = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        self.los[self.los==255] = float('NaN'); self.los *= self.scale
        band.FlushCache()
        self.nlines,self.ncols = ds.RasterYSize, ds.RasterXSize
        logger.info('Number of pixel: {}'.format(len(self.los.flatten())))
        logger.info('LOS, Nlines: {}, Ncols: {}'.format(self.nlines,self.ncols))
        del ds, band
       
        # compute losmax,losmin
        if self.bounds == None:
            self.losmax = np.nanpercentile(self.los,98)
            self.losmin = np.nanpercentile(self.los,2)
        else:
            self.losmin = self.bounds[0]
            self.losmax = self.bounds[1]
        
        ds = gdal.Open(self.incidence_file,gdal.GA_ReadOnly)
        # param output files
        self.gt = ds.GetGeoTransform()
        self.projref = ds.GetProjectionRef()
        self.driver = gdal.GetDriverByName('GTiff')

        band = ds.GetRasterBand(1)
        self.look = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        logger.info('Look, Nlines: {}, Ncols: {}'.format(ds.RasterYSize,ds.RasterXSize))
        if ds.RasterYSize != self.nlines or ds.RasterXSize != self.ncols:
            logger.critical('Look file has different size. Exit!')
            exit()
        # self.look[np.isnan(self.los)] = float('NaN'))
        del ds, band

        if self.sigmaf is not None:
            ds = gdal.Open(self.sigmaf,gdal.GA_ReadOnly)
            band = ds.GetRasterBand(1)
            logger.info('SIGMA, Nlines: {}, Ncols: {}'.format(ds.RasterYSize,ds.RasterXSize)) 
            if ds.RasterYSize != self.nlines or ds.RasterXSize != self.ncols:
                logger.critical('Sigma file has different size. Exit!')
                exit()
            self.sigma = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
            self.sigma[self.sigma==255] = float('NaN'); self.sigma = self.scale*self.sigma 
            band.FlushCache()
            del ds, band
        else:
            self.sigma = np.ones((self.nlines,self.ncols))

        # compute siglosmax,siglosmin
        self.sig_losmax = np.nanpercentile(self.sigma,98)
        self.sig_losmin = np.nanpercentile(self.sigma,2)

        ds = gdal.Open(self.heading_file,gdal.GA_ReadOnly)
        band = ds.GetRasterBand(1)
        self.head = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        logger.info('Heading, Nlines: {}, Ncols: {}'.format(ds.RasterYSize,ds.RasterXSize))
        if ds.RasterYSize != self.nlines or ds.RasterXSize != self.ncols:
            logger.critical('Look file has different size. Exit!')
            exit()
        # self.head[np.isnan(self.los)] = float('NaN') 
        band.FlushCache()
        del ds, band

        # convert head, look to angle phi, theta in rad
        # theta: angle between LOS and horizontal (positive in anti-clockwise direction)
        self.theta = np.deg2rad(90.-self.look)
        # phi: horizontal angle between LOS and East (positive in anti-clockwise direction)
        self.phi = np.deg2rad(-90-self.head)

        logger.info('Average LOOK:{0:.5f}, THETA:{1:.5f} angles'.\
            format(np.nanmean(self.look),np.nanmean(np.rad2deg(self.theta))))

        logger.info('Average HEADING:{0:.5f}, PHI:{1:.5f} angles'.\
            format(np.nanmean(self.head),np.nanmean(np.rad2deg(self.phi))))

        # compute proj ENcomp3
        self.proj=[np.cos(self.phi)*np.cos(self.theta),
               np.sin(self.phi)*np.cos(self.theta),
               np.sin(self.theta)
              ]

        logger.info('Average horizontal LOS projection to east, north, up: {0:.5f} {1:.5f} {2:.5f}'.\
            format(np.nanmean(self.proj[0]),np.nanmean(self.proj[1]),np.nanmean(self.proj[2])))

        self.proj=[
                np.cos(self.slope)*np.cos(self.theta)*(np.cos(self.rot)*np.cos(self.phi) - np.sin(self.rot)*np.sin(self.phi)) - np.sin(self.slope)*np.sin(self.theta),
                np.cos(self.theta)*(np.sin(self.rot)*np.cos(self.phi) + np.cos(self.rot)*np.sin(self.phi)),
                np.sin(self.slope)*np.cos(self.theta)*(np.cos(self.rot)*np.cos(self.phi) - np.sin(self.rot)*np.sin(self.phi)) + np.cos(self.slope)*np.sin(self.theta)
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
    try:
      sys.path.append(path.dirname(path.abspath(fname)))
      exec ("from "+path.basename(fname)+" import *")
    except:
      exec(open(path.abspath(fname)).read())
  
  except Exception as e: 
    logger.critical('Problem in input file')
    logger.critical(e)
    print(network.__doc__)
    exit()

def compute_slope_aspect(path, filter_size=2.0):
   
    #### LOAD DEM
    ds = gdal.Open(path,gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    topo = band.ReadAsArray()
    ncols, nlines = ds.RasterYSize, ds.RasterXSize

    gt = ds.GetGeoTransform()
    projref = ds.GetProjectionRef()
    drv = gdal.GetDriverByName('GTiff')
    
    # Filter to smooth data and remove noise
    filtered_topo = scipy.ndimage.gaussian_filter(topo, sigma=(filter_size, filter_size))   
 
    # Get middle latitude for geospatial calculations
    lats = gt[3] + (np.arange(nlines) * gt[5])
    lat_mean = np.mean(lats)
    print("GeoTransform:", gt)
    print("Average latitude:", lat_mean)

    # Determine resolution based on coordinate system
    res_x, res_y = gt[1], gt[5]
    srs = osr.SpatialReference(wkt=projref)

    if srs.IsGeographic():
        res_x *= 40075e3 / 360  # Convert degrees to meters
        res_y *= (40075e3 / 360) * np.cos(np.deg2rad(lat_mean))
        print("Geographic coordinates (converted to meters)")
    else:
        print("Projected coordinates (meters)")

    print(f"Resolution: dx={res_x:.2f}, dy={res_y:.2f}")

    # Compute gradient (first axis = Y, second axis = X)
    dsm_dy, dsm_dx = np.gradient(filtered_topo, res_y, res_x)

    # Calculate slope and aspect
    slope = np.sqrt(dsm_dx**2 + dsm_dy**2)
    aspect = np.arctan2(dsm_dy, -dsm_dx)  # clockwise from east
   
    # Create aspect and slope files
    dst = drv.Create('slope.tif', nlines, ncols, 1, gdal.GDT_Float32)
    bandt = dst.GetRasterBand(1)
    bandt.WriteArray(np.rad2deg(slope))
    dst.SetGeoTransform(gt)
    dst.SetProjection(projref)
    bandt.FlushCache()

    dst = drv.Create('aspect.tif', nlines, ncols, 1, gdal.GDT_Float32)
    bandt = dst.GetRasterBand(1)
    bandt.WriteArray(np.rad2deg(aspect))
    dst.SetGeoTransform(gt)
    dst.SetProjection(projref)
    bandt.FlushCache()

    # Plot DEM, Slope, Py and aspect
    fig = plt.figure(1,figsize=(11,7))
    cmap = cm.terrain
    # Plot topo
    ax = fig.add_subplot(2,2,1)
    cax = ax.imshow(filtered_topo,cmap='Greys_r',vmax=np.nanpercentile(filtered_topo,98),vmin=np.nanpercentile(filtered_topo,2))
    ax.set_title('DEM',fontsize=6)
    fig.colorbar(cax, orientation='vertical')

    ax = fig.add_subplot(2,2,2)
    cax = ax.imshow(np.rad2deg(slope),cmap='Greys_r',vmax=np.nanpercentile(np.rad2deg(slope),98),vmin=np.nanpercentile(np.rad2deg(slope),2))
    ax.set_title('Slope',fontsize=6)
    fig.colorbar(cax, orientation='vertical')

    ax = fig.add_subplot(2,2,3)
    cax = ax.imshow(dsm_dy,cmap='coolwarm',vmax=np.nanpercentile(dsm_dy,98),vmin=np.nanpercentile(dsm_dy,2))
    ax.set_title('Py',fontsize=6)
    fig.colorbar(cax, orientation='vertical')

    ax = fig.add_subplot(2,2,4)
    cax = ax.imshow(np.rad2deg(aspect),cmap='Greys_r',vmax=180,vmin=-180)
    ax.set_title('aspect',fontsize=6)
    fig.colorbar(cax, orientation='vertical')
    
    fig.savefig('dem_slope_aspect.pdf',format='PDF',dpi=300)

    if PLOT:
        plt.show()      

    return aspect, slope

#### compute rotations
# rot: tourne l'axe N vers l'axe E. 
# slope: tourne l'axe Up vers l'axe E. 
min_slope = 2

comp_name = []
rot, slope = 0, 0
# rotation angle: angle between comp1 and East
if ('DEM' in locals()) and (DEM is not None):
    logger.info('DEM is defined, compute aspect and slope for each pixel.')
    rot, slope = compute_slope_aspect(DEM)
    comp = [0, 2]
    #comp = [0]
    name1 = 'Line of Max Slope'
    comp_name.append(name1)
    logger.info('Invert components: {}'.format(name1))
    name3 = 'Normal to line of Max Slope'
    logger.info('Invert components: {}'.format(name3))
    comp_name.append(name3)
else:
    logger.info('DEM is not defined, read horizontale rotation in clockwise rotation in input file (default, rotation=0)')
    #slope = np.deg2rad(30)
    if 'rotation' not in locals():
        logger.warning('Define angle between north and new axis in the rotation variable in the input file')
        logger.warning('Set to 0 by default')
        rotation = 0
    rot = np.deg2rad(rotation)
    # define invert components
    comp = np.array(comp) - 1 
    for n in range(len(comp)):
        if int(comp[n]) == 0:
            name = 'East + {} deg (anti-clockwise)'.format(np.rad2deg(rot))
        elif int(comp[n]) == 1:
            name = 'North + {} deg (anti-clockwise)'.format(np.rad2deg(rot))
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

# Load data
N = len(comp); M = len(insar)
vmax=0; vmin=0; sigmax=0; sigmin=0

for i in range(M):
    
    insar[i].load(rot,slope)
    
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
        logger.critical('Size input maps are not idenical. Exit!')
        exit()
    # get param output files
    gt = insar[0].gt
    projref = insar[0].projref
    driver = insar[0].driver

    # Apply reference zone shift...
    if insar[i].ref_zone is not None:
      if insar[i].ref_zone[1]<=insar[i].los.shape[0] and insar[i].ref_zone[3]<=insar[i].los.shape[1]:
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

logger.info('Plot DATA ....') 
print()
try:
    from matplotlib.colors import LinearSegmentedColormap
    cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
    cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
except:
    cmap=cm.rainbow
cmap_r = cmap.reversed()

fig=plt.figure(0, figsize=(11,7))

for i in range(M):
    d = insar[i]
    # plot LOS
    ax = fig.add_subplot(2,M,i+1)
    cax = ax.imshow(d.los[ibeg:iend,jbeg:jend],cmap=cmap_r,vmax=vmax,vmin=vmin,interpolation=None)
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
    cax = ax.imshow(d.sigma[ibeg:iend,jbeg:jend],cmap=cmap_r,vmax=sigmax,vmin=sigmin,interpolation=None)
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
fig.savefig('data_decomposition_{}.pdf'.format(output),format='PDF',dpi=300)
if PLOT:
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
    if i % 200 == 0:
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
                Cd = np.diag(rms**2, k = 0)
                pars = np.dot(np.linalg.inv(np.dot(np.dot(G.T,np.linalg.inv(Cd)),G)),np.dot(np.dot(G.T,np.linalg.inv(Cd)),data))
                #pars = lst.lstsq(G,data)[0]
                #pars = np.dot(np.linalg.inv(G),data)

                if iter > 1:
                    _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
                    _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
                    pars = opt.fmin_slsqp(_func,pars,fprime=_fprime,iter=iter,full_output=True,iprint=0,acc=acc)[0]
                
                for n in range((N)):        
                    disp[i,j,int(comp[n])] = pars[n]

                Cd = np.diag(rms**2, k = 0)
                try:
                    sigmam = np.sqrt(np.linalg.inv(np.dot(np.dot(G.T,np.linalg.inv(Cd)),G)))
                except:
                    sigmam = np.ones((G.shape[1]))*float('NaN')
                
                for n in range((N)): 
                    sdisp[i,j,int(comp[n])] = np.abs(sigmam[n,n])

            except:
                pars = np.zeros((N))
                print(i,j)
                print('data:',data)
                print('G:',G)
                print()

################################
# CLEAN RESULTS
################################

mask = np.ones((nlines,ncols,3)) 
for n in range(N):
    # clean outlisers
    d = as_strided(disp[:,:,int(comp[n])])
    s = as_strided(sdisp[:,:,int(comp[n])])
    m = as_strided(mask[:,:,int(comp[n])])
    # clean based on uncertainties
    index = s>np.nanpercentile(s,98.)
    #d[index] = float('NaN') 
    m[index] = 0.
    index = s/abs(d) > 3. 
    #d[index] = float('NaN')
    m[index] = 0.
    # clean based on outiliers
    index = np.logical_or(d>np.percentile(d,99.8),d<np.percentile(d,.2))
    #d[index]= float('NaN') 
    m[index] = 0.
    
    if 'DEM' in locals():
       if DEM is not None:
         # mask based on aspect
         indice = np.nonzero(np.logical_and(np.rad2deg(rot)>60,np.rad2deg(rot)<120))
         m[indice] =  0.
         #d[indice] = float('NaN')
         indice = np.nonzero(np.logical_and(np.rad2deg(rot)<-60,np.rad2deg(rot)>-120))
         m[indice] =  0.
         #d[indice] = float('NaN')

         # clean flat regions: inversion invalid
         d[np.rad2deg(slope)<min_slope] = float('NaN')
         m[np.rad2deg(slope)<min_slope] = float('NaN') 

   
################################
# PLOT RESULTS
################################

fig=plt.figure(3, figsize=(10,12))

for n in range(N):

    data = disp[ibeg:iend,jbeg:jend,int(comp[n])]
    m = mask[ibeg:iend,jbeg:jend,int(comp[n])]
    data[data==0] = float('NaN')

    # plot comps
    vmax =  np.nanpercentile(data,98)
    vmin = np.nanpercentile(data,2) 
    ax = fig.add_subplot(3,N,n+1)
    cax = ax.imshow(data,cmap=cmap_r,vmax=vmax,vmin=vmin)
    ax.set_title('{}'.format(comp_name[n]))

    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    # plot SIGMA LOS
    sigdata = sdisp[ibeg:iend,jbeg:jend,int(comp[n])]
    # sigdata[sigdata==0] = float('NaN')
    # sigdata[sigdata==1] = float('NaN')

    vmax=np.nanpercentile(sigdata,99)
    ax = fig.add_subplot(3,N,n+1+N)
    cax = ax.imshow(sigdata,cmap=cmap_r,vmax=vmax,vmin=0)
    ax.set_title('SIGMA {}'.format(comp_name[n]))
    
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)
    
    ax = fig.add_subplot(3,N,n+1+N*2)
    cax = ax.imshow(m,cmap=cm.Greys)
    ax.set_title('MASK {}'.format(comp_name[n]))

    # Colorbar
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

fig.tight_layout()
fig.savefig('decomposition_{}.pdf'.format(output), format='PDF',dpi=300)

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

    ds = driver.Create('mask_U{}_{}.tif'.format(comp[n],output), ncols, nlines, 1, gdal.GDT_Float32)
    band = ds.GetRasterBand(1)
    band.WriteArray(mask[:,:,int(comp[n])])
    ds.SetGeoTransform(gt)
    ds.SetProjection(projref)
    band.FlushCache()

plt.show()
