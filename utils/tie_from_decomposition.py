#!/usr/bin/env python2
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Simon DAOUT (Oxford)
############################################

from sys import argv,exit,stdin,stdout
import getopt
import os, math
from os import path
import logging
import sys

import numpy as np
import scipy.optimize as opt
import scipy.linalg as lst
import gdal, pyproj
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
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

### Define GNSS class
class gpsnetwork:
    def __init__(self,network,reduction,wdir,scale=1,weight=1.):

        self.network=network
        self.reduction=reduction
        self.wdir=wdir
        self.scale=scale
        # inititialisation
        self.sigmad=1./weight
        self.N = 0
        self.x,self.y= [], []

    def load(self,rot):
        self.rot = rot
        fname=self.wdir + self.network
        if not path.isfile(fname):
            raise ValueError("invalid file name: " + fname)
        else:
            print('Load GPS time series: ', fname)
            print()
        f=open(fname,"r")
        lon, lat, ve, vn, vu, se, sn, su, cen, name = \
        np.loadtxt(f,comments='#',unpack=True,dtype='f,f,f,f,f,f,f,f,f,S4')
        
        self.lon, self.lat, self.ve, self.vn, self.vu, self.se, self.sn, self.su, self.cen,\
         self.name=np.atleast_1d(lon, lat, ve, vn, vu, se, sn, su, cen, name)

        # rotation of the axis of an angle rot in clockwise direction
        self.vep = np.cos(rot)*self.ve - np.sin(rot)*self.vn
        self.vnp = np.sin(rot)*self.ve + np.sin(rot)*self.vn

        self.sep=((self.se*np.cos(self.rot))**2 + (self.sn*np.sin(self.rot))**2)**0.5
        self.snp=((self.se*np.sin(self.rot))**2 + (self.sn*np.cos(self.rot))**2)**0.5

        self.N = len(self.name)*3
        self.d = np.column_stack([self.vep,self.vnp,self.vu])
        self.sigmad = self.sigmad*np.column_stack([self.sep,self.snp,self.su])
        self.Npoints = len(self.name)

    def info(self):
        print()
        print('GPS network:',self.network)
        print('Number of stations:', self.Npoints)
        print('Number of points:', self.N)


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
    
    def __init__(self,name,wdir,reduction,lookf,headf,sigmaf=None,scale=1,scale_sig=1,format='ROIPAC',bounds=None,outname=None):
        self.name = name
        self.wdir = wdir
        self.reduction = reduction
        self.path = wdir + name
        self.lookf = wdir + lookf
        self.headf = wdir + headf
        self.sigmaf = wdir + sigmaf
        self.scale = scale
        self.scale_sig = scale_sig
        self.bounds = bounds
        self.format=format
        if outname is None:
            self.outname = os.path.splitext(self.name)[0] + '_projected.tif'
        else:
            self.outname = outname

    def load(self,rot):
        logger.info('Read track: {}'.format(self.name))
        ds = gdal.Open(self.path,gdal.GA_ReadOnly)
        band = ds.GetRasterBand(1)
        self.los = self.scale*band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        self.ysize,self.xsize = ds.RasterYSize, ds.RasterXSize
        logger.info('Number of pixel: {}'.format(len(self.los.flatten())))
        logger.info('ysize: {}, xsize: {}'.format(self.ysize,self.xsize))
        self.left = ds.GetGeoTransform()[0]
        self.top = ds.GetGeoTransform()[3]
        self.xres = ds.GetGeoTransform()[1]
        self.yres = ds.GetGeoTransform()[5]
        self.projection = ds.GetProjection()
        self.right = self.left + self.xsize * self.xres
        self.bottom = self.top + self.ysize * self.yres
        self.pix_lin, self.pix_col = np.indices((ds.RasterYSize,ds.RasterXSize))
        self.lat,self.lon = self.top + self.yres*self.pix_lin, self.left+self.xres*self.pix_col
        # convert 0 and 255 to NaN
        self.los[self.los==0.] = np.float('NaN')
        self.los[self.los==255] = np.float('NaN')
        band.FlushCache()
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
        self.driver = gdal.GetDriverByName('GTiff')

        band = ds.GetRasterBand(1)
        self.look = np.ones((self.ysize,self.xsize))
        self.look[:ds.RasterYSize,:ds.RasterXSize] = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:self.ysize,:self.xsize]
        self.look[np.isnan(self.los)] = np.float('NaN') 
        band.FlushCache()
        del ds, band

        if self.sigmaf is not None:
            self.sigma = np.ones((self.ysize,self.xsize))
            ds = gdal.Open(self.sigmaf,gdal.GA_ReadOnly)
            band = ds.GetRasterBand(1)
            self.sigma[:ds.RasterYSize,:ds.RasterXSize] = self.scale_sig*band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:self.ysize,:self.xsize]
            self.sigma[np.isnan(self.los)] = np.float('NaN') 
            band.FlushCache()
            del ds, band
        else:
            self.sigma = np.ones((self.ysize,self.xsize))

        # compute siglosmax,siglosmin
        self.sig_losmax = np.nanpercentile(self.sigma,98)
        self.sig_losmin = np.nanpercentile(self.sigma,2)

        ds = gdal.Open(self.headf,gdal.GA_ReadOnly)
        band = ds.GetRasterBand(1)
        self.head  = np.ones((self.ysize,self.xsize))
        self.head[:ds.RasterYSize,:ds.RasterXSize] = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[:self.ysize,:self.xsize]
        self.head[np.isnan(self.los)] = np.float('NaN') 
        band.FlushCache()
        del ds, band

        if self.format=='ROIPAC':
            # convert head, look to angle phi, theta in rad
            # theta: vertical angle between LOS and vertical
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
            self.proj=[(np.cos(rot)*self.proj[0] + np.sin(rot)*self.proj[1])*np.cos(self.theta),
                   (-np.sin(rot)*self.proj[0] + np.cos(rot)*self.proj[1])*np.cos(self.theta),
                   self.proj[2]
                   ]

        logger.info('Average LOS projection to comp1, comp2, comp3: {0:.5f} {1:.5f} {2:.5f}'.\
            format(np.nanmean(self.proj[0]),np.nanmean(self.proj[1]),np.nanmean(self.proj[2])))
        print()

    def extract_pixel_value(self, lon, lat, n):

        x = int((lon-self.left)/self.xres+0.5)
        y = int((lat - self.top) / self.yres + 0.5)

        pixel_values = self.los[y - n: y + n + 1, x - n: x + n + 1]
        pixel_sigma = self.sigma[y - n: y + n + 1, x - n: x + n + 1]
        pixel_proj = [self.proj[0][y - n: y + n + 1, x - n: x + n + 1],\
        self.proj[1][y - n: y + n + 1, x - n: x + n + 1],\
        self.proj[2][y - n: y + n + 1, x - n: x + n + 1]]

        index = np.nonzero(~np.isnan(pixel_values))
        if len(index[0]) > 0:

            #m = np.nanmean(pixel_values)
            m = np.nansum(pixel_values[index] * pixel_sigma[index]) / np.sum(pixel_sigma[index])
            std = np.nanstd(pixel_values)
            proj = [np.nanmean(pixel_proj[0][index]),np.nanmean(pixel_proj[1][index]),np.nanmean(pixel_proj[2][index])]

        else:
            m = np.float('NaN')
            std = np.float('NaN')
            proj = [np.float('NaN'),np.float('NaN'),np.float('NaN')]
                            
        if m == 0:  #if only NaN nanmean is 0 
            m = np.float('NaN')
        if std == 0:
            std = np.float('NaN')
        if proj ==0:
            proj = [np.float('NaN'),np.float('NaN'),np.float('NaN')]

        return m, std, proj

    def project(self,ngps=0):
        x, y = UTM(self.lon.flatten(), self.lat.flatten())
        G = np.zeros((len(self.los.flatten()),3))
        G[:,0] = y
        G[:,1] = x
        G[:,2] = 1

        self.ramp_array = np.dot(G, self.ramp).reshape(self.ysize, self.xsize)
        self.los_corr = self.ramp_array + self.los

        # if no gps, no reference, remove a cst
        if ngps < 1 :
            cst = np.nanmean(self.los_corr)
            self.los_corr = self.los_corr - cst

    def save(self):

        if path.exists(self.outname):
            os.remove(self.outname)

        # Export merged data to tif format.
        driver = gdal.GetDriverByName("GTiff")
        outdata = driver.Create(self.outname, self.xsize, self.ysize, 1, gdal.GDT_Float32)
        outdata.SetGeoTransform([self.left, self.xres, 0, self.top, 0, self.yres])  ##sets same geotransform as input
        outdata.SetProjection(self.projection)  ##sets same projection as input
        outdata.GetRasterBand(1).WriteArray(self.los_corr)
        outdata.FlushCache()

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

    # init logger 
    if level == 'debug':
        logging.basicConfig(level=logging.DEBUG,\
        format='%(lineno)s -- %(levelname)s -- %(message)s')
        logger = logging.getLogger('tie_from_decomposition.log')
        logger.info('Initialise log file {0} in DEBUG mode'.format('tie_from_decomposition.log'))

    else:
        logging.basicConfig(level=logging.INFO,\
        format='%(lineno)s -- %(levelname)s -- %(message)s')
        logger = logging.getLogger('tie_from_decomposition.log')
        logger.info('Initialise log file {0} in INFO mode. Use option -v for a DEBUG mode'.format('tie_from_decomposition.log'))


    # rotation angle: angle between comp1 and East
    rot = np.deg2rad(-rotation)

    # load gnss
    Mgnss = len(gnss)
    for i in range(Mgnss):
        gnss[i].load(-rot)
        gnss[i].info()

    # Load insar
    M = len(insar)
    vmax=0; vmin=0; sigmax=0; sigmin=0
    for i in range(M):
        insar[i].load(rot)

        # extract bounds
        if insar[i].losmax > vmax:
            vmax = insar[i].losmax
        if insar[i].losmin < vmin:
            vmin = insar[i].losmin
        if insar[i].sig_losmax > sigmax:
            sigmax = insar[i].sig_losmax
        if insar[i].sig_losmin < sigmin:
            sigmin = insar[i].sig_losmin

    # define invert components
    comp = np.array(comp) - 1 
    Ncomp = len(comp)
    comp_name = []
    for n in range(Ncomp):
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

    ################################
    # plot DATA
    ################################

    # print()
    # logger.info('Plot DATA ....') 
    # try:
    #     from matplotlib.colors import LinearSegmentedColormap
    #     cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
    #     cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
    #     cmap_r = cmap.reversed()
    # except:
    #     cmap=cm.rainbow

    # fig=plt.figure(0, figsize=(16,9))

    # for i in range(M):
    #     d = insar[i]
    #     # plot LOS
    #     ax = fig.add_subplot(2,M,i+1)
    #     cax = ax.imshow(d.los,cmap=cmap_r,vmax=vmax,vmin=vmin,interpolation=None)
    #     ax.set_title('{}'.format(d.reduction))
    #     # Colorbar
    #     divider = make_axes_locatable(ax)
    #     c = divider.append_axes("right", size="5%", pad=0.05)
    #     plt.colorbar(cax, cax=c)
    #     plt.setp( ax.get_xticklabels(), visible=False)
    #     plt.setp( ax.get_yticklabels(), visible=False)

    #     for j in range(Mgnss):
    #         gps = gnss[j]
    #         index = np.nonzero((gps.lon>insar[i].left)&(gps.lon<insar[i].right)&(gps.lat>insar[i].bottom)&(gps.lat<insar[i].top))
    #         gpslon,gpslat = gps.lon[index], gps.lat[index]
    #         gpsve,gpsvn = gps.ve[index], gps.vn[index]

    #         x = map(int, (gpslon - insar[i].left) / insar[i].xres + 0.5)
    #         y = map(int, (gpslat - insar[i].top) / insar[i].yres + 0.5)
    #         ax.quiver(x,y,gpsve,gpsvn,scale = 100,width = 0.005,color = 'black')
        
    #     # plot SIGMA LOS
    #     ax = fig.add_subplot(2,M,i+1+M)
    #     cax = ax.imshow(d.sigma,cmap=cmap_r,vmax=sigmax,vmin=sigmin,interpolation=None)
    #     ax.set_title('SIGMA {}'.format(d.reduction))
    #     # Colorbar
    #     divider = make_axes_locatable(ax)
    #     c = divider.append_axes("right", size="5%", pad=0.05)
    #     plt.colorbar(cax, cax=c)
    #     plt.setp( ax.get_xticklabels(), visible=False)
    #     plt.setp( ax.get_yticklabels(), visible=False)

    # # fig.tight_layout()
    # fig.savefig('data_tie_decomposition{}.png'.format(output),format='PNG',dpi=300)
    # # plt.show()
    # # sys.exit()

    ################################
    # Extract pixel values in common 
    ################################
    
    n = 2 # window size for insar
    nn = 20 # window size for gnss
    d = insar[0]
    data = []
    G = []
    rms = []


    #Ngnss = sum(map((lambda x: getattr(x,'Npoints')),gnss))
    # Nmax is the number max of points to be inverted
    # Size max of the matrix is Nmax*Ncomp + Ngnss*Ncomp + Mramp
    # hardcode number max of gnss within insar0 to 20
    Mp = Nmax*Ncomp + 20*Ncomp
    Mbasis = (3*M) + Mp

    # initialisation number of points
    p = 0 # count number of insar points
    ngps = 0 # count number of gnss

    # loop over the first insar data set values
    for x,y,i,j in zip(d.lon.flatten()[::sampling],d.lat.flatten()[::sampling],\
        d.pix_lin.flatten()[::sampling],d.pix_col.flatten()[::sampling]):
        
        do = False
        for k in range(1,M):      
            
            if ((i % 100) == 0) and (j==0):
                print('Processing line: {}, data set: {}'.format(i,insar[k].reduction))
            
            m,std,proj = insar[k].extract_pixel_value(x,y,n)
            m0 = np.mean(d.los[i-n: i+n+1, j-n: j+n+1])
            rms0 = np.std(d.los[i-n: i+n+1, j-n: j+n+1])

            # if m is not NaN, we fill data and G matrix
            if ~np.isnan(m) and ~np.isnan(std) and ~np.isnan(np.sum(proj)) \
            and m0 != 0  and ~np.isnan(m0) and ~np.isnan(rms0) and p < Nmax:
                #print(m)
                do = True
                
                # value insark
                data.append(m)
                rms.append(std)

                # build G matrix for insark
                Gt = np.zeros((1,Mbasis))
                for n in range(Ncomp):  
                    Gt[0, n + p*Ncomp] = proj[int(comp[n])]
            
                # ramp param
                Gt[0, Mp + 3*k + 2] = 1               
                x, y = UTM(np.nanmean(insar[k].lon[i-n: i+n+1, j-n: j+n+1]), np.nanmean(insar[k].lat[i-n: i+n+1, j-n: j+n+1])) 
                Gt[0, Mp + 3*k ], Gt[0, Mp + 3*k + 1] = y, x
                G.append(Gt)

        # loop over gnss first and select zone within insar0
        gpsx, gpsy = [], []
        for g in range(Mgnss):
            gps = gnss[g]

            # select gnss withing the insar frame 0
            index = np.nonzero((gps.lon>d.left)&(gps.lon<d.right)&(gps.lat>d.bottom)&(gps.lat<d.top))
            gpslon,gpslat = gps.lon[index], gps.lat[index]
            gpsd = np.array([gps.vep[index], gps.vnp[index], gps.vu[index]]).T
            gpssigmad = np.array([gps.sep[index], gps.snp[index], gps.su[index]]).T

            # extract pixel values
            gpsx = map(int, (gpslon - d.left) / d.xres + 0.5)
            gpsy = map(int, (gpslat - d.top) / d.yres + 0.5)

            # if we find a gnss point, then include insar too
            for gpsxx,gpsyy in zip(gpsx,gpsy):

                index =  np.nonzero((gpsxx < i+nn+1) & (gpsxx > i-nn) & (gpsyy < j+nn+1) & (gpsyy > j-nn))
                if not np.isnan(np.mean(index)) :
                    do = True
                    
                    Gt = np.zeros((1,Mbasis))
                    for n in range(Ncomp):
                        Gt[0, n + p*Ncomp] =  1

                        # fill data and rms vector
                        data.append(gpsd[index].flatten()[int(comp[n])])
                        rms.append(gpssigmad[index].flatten()[int(comp[n])])

                        # fill G for Gnss
                        G.append(Gt)
                    
                    # count number of added gnss
                    ngps = ngps + 1
        
        if do:
            # add value insar0
            data.append(np.nanmean(d.los[i-n: i+n+1, j-n: j+n+1]))
            rms.append(np.nanstd(d.los[i-n: i+n+1, j-n: j+n+1]))

            # build G matrix for insar0
            Gt = np.zeros((1,Mbasis))
            for n in range(Ncomp):
                Gt[0, n + p*Ncomp] = np.nanmean(d.proj[int(comp[n])][i-n: i+n+1, j-n: j+n+1])
            
            Gt[0, Mp + 2] = 1
            x, y = UTM(np.nanmean(d.lon[i-n: i+n+1, j-n: j+n+1]), \
                np.nanmean(d.lat[i-n: i+n+1, j-n: j+n+1])) 
            
            Gt[0, Mp], Gt[0, Mp + 1] = y, x
            G.append(Gt)

            # one more point
            p = p + 1

            # loop over the data sets
            if p == Nmax:
                print('Number of points superior to {} ! Consider more undersampling.'.format(Nmax))


    ################################
    # Inversion
    ################################

    print()
    logger.info('Inversion ....')

    data = np.array(data).flatten()
    print('Number of inverted points: ', p)
    print('Number of inverted GNSS points: ', ngps )

    rms = np.array(rms).flatten()
    G = np.array(G).reshape(len(data),Mbasis)

    # check for NaN in data vector
    index = np.flatnonzero(~np.isnan(data))
    data = data[index]
    rms = rms[index]
    G = G[index,:]

    # check for NaN in rms vector
    index = np.flatnonzero(~np.isnan(rms))
    data = data[index]
    rms = rms[index]
    G = G[index,:]

    # check for NaN in G
    data = data[~np.isnan(G.any(axis=1))]
    G = G[~np.isnan(G.any(axis=1))]

    # Now lets remove column in G matrix with only zeros
    G = G[:,~np.all(G == 0, axis = 0)]
    #  recompute Mp
    Mp = G.shape[1] - 3*M

    pars = lst.lstsq(G,data)[0]
    # if iter > 1:
    #     _func = lambda x: np.sum(((np.dot(G,x)-data)/rms)**2)
    #     _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-data)/rms)
    #     pars = opt.fmin_slsqp(_func,pars,fprime=_fprime,iter=iter,full_output=True,iprint=0,acc=acc)[0]

    if iter > 1:
        _func = lambda x: np.sum(((np.dot(G,x)-data))**2)
        _fprime = lambda x: 2*np.dot(G.T, (np.dot(G,x)-data))
        pars = opt.fmin_slsqp(_func,pars,fprime=_fprime,iter=iter,full_output=True,iprint=0,acc=acc)[0]

    # apply ramp param
    for k in range(M):
        insar[k].ramp = pars[Mp + 3*k : Mp + 3*k + 3]
        # print(insar[i].ramp)

        # compute and save results
        insar[k].project(ngps)

    ################################
    # plot DATA
    ################################

    print()
    logger.info('Plot Corrected DATA ....') 
    try:
        from matplotlib.colors import LinearSegmentedColormap
        cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
        cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
        cmap_r = cmap.reversed()
    except:
        cmap=cm.rainbow

    fig=plt.figure(1, figsize=(16,9))

    for i in range(M):
        d = insar[i]
        # plot LOS
        ax = fig.add_subplot(2,M,i+1)
        cax = ax.imshow(d.los,cmap=cmap_r,vmax=vmax,vmin=vmin,interpolation='nearest')
        ax.set_title('{}'.format(d.reduction))
        ax.set_title('LOS {}'.format(d.reduction))

        # plot gnss
        for j in range(Mgnss):
            gps = gnss[j]
            index = np.nonzero((gps.lon>insar[i].left)&(gps.lon<insar[i].right)&(gps.lat>insar[i].bottom)&(gps.lat<insar[i].top))
            gpslon,gpslat = gps.lon[index], gps.lat[index]
            gpsve,gpsvn = gps.ve[index], gps.vn[index]

            x = map(int, (gpslon - insar[i].left) / insar[i].xres + 0.5)
            y = map(int, (gpslat - insar[i].top) / insar[i].yres + 0.5)
            ax.quiver(x,y,gpsve,gpsvn,scale = 100,width = 0.005,color = 'black')

        # Colorbar
        divider = make_axes_locatable(ax)
        c = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax, cax=c)
        plt.setp( ax.get_xticklabels(), visible=False)
        plt.setp( ax.get_yticklabels(), visible=False)

        # plot  LOS corr
        ax = fig.add_subplot(2,M,i+1+M)
        cax = ax.imshow(d.los_corr,cmap=cmap_r,vmax=vmax,vmin=vmin,interpolation='nearest')
        ax.set_title('CORR LOS {}'.format(d.reduction))
        # Colorbar
        divider = make_axes_locatable(ax)
        c = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax, cax=c)
        plt.setp( ax.get_xticklabels(), visible=False)
        plt.setp( ax.get_yticklabels(), visible=False)

    # fig.tight_layout()
    fig.savefig('data_tie_decomposition{}.png'.format(output),format='PNG',dpi=300)
    plt.show()

    ################################
    # save DATA
    ################################

    print()
    logger.info('Save Output files ....')

    # save ramp param
    for i in range(M):
        insar[i].save()




