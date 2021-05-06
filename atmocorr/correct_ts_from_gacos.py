#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT (ISTerre)
################################################################################

"""\
correct_ts_from_gacos.py
-------------
Correct Time Series data from Gacos atmospheric models. 1) convert .ztd files to .tif format, 2) re-project and re-resampel atmospheric models in geographic coordinates
3) correct data
To skip 1) and 2), use the --load=gacos_cube argument, with gacos_cube the 3-D matrix containing the time series of the gacos models 

Usage: 
    correct_ts_from_gacos.py [--cube=<file>] [--path=<path>] [--crop=<values>] [--proj=<value>] \
    [--load=<file>] [--list_images=<file>] [--imref=<value>] [--gacos2data=<value>] [--los=<file>] \
    [--plot=<yes|no>] [--rmspixel=<file>] [--threshold_rms=<value>] [--ref_zone=<values>] [--ramp=<cst|lin>] \
    [--crop_emp=<values>] [--topofile=<path>] [--fitmodel=<yes|no>] [--output_cube=<file>] [--output_gacos=<file>]\
      [--lectfile=<path>] [--meanlos=<value>]

correct_ts_from_gacos.py -h | --help

Options:
-h --help             Show this screen.
--cube=<file>         Path to the cube of TS displacements file to be corrected from atmo models [default: depl_cumul]
--path=<dir>          Path to .ztd data [default: ./GACOS/]
--list_images=<file>  Path to list images file made of 5 columns containing for each images 1) number 2) Doppler freq (not read) 3) date in YYYYMMDD format 4) numerical date 5) perpendicular baseline [default: list_images.txt]
--imref=<value>       Reference image number [default: 1]
--proj=<value>        EPSG for projection GACOS map [default: 4326]
--crop                Crop GACOS data to data extent in the output projection, eg.--crop=xmin,ymin,xmax,ymax 
--load=<file>        If a file is given, load directly GACOS cube 
--ref_zone=<lin_start,lin_end,col_start,col_end> Starting and ending lines and col numbers where phase is set to zero on the corrected data [default: None]
--ramp=<cst|lin|quad> Estimate a constant, linear or quadratic ramp in x and y in addition to the optional los/gacos relationship [default: cst].
--fitmodel=<yes|no>  If yes, then estimate the proportionlality between gacos and data in addition to the defined ramp
--crop_emp=<value>    Crop option for the empirical estimations (ie. ramp and fitmodel) (eg: 0,ncol,0,nlines)
--topo=<file>         Use DEM file to mask very low or high altitude values in the empirical estimation 
--gacos2data=<value>  Scaling value between zenithal gacos data (m) and desired output (e.g data in mm positive toward the sat.: -1000) [default: 1]
--los=<file>          LOS angle between the vertical and the Look Direction use to convert zenital delays to LOS delays [default: None]
--meanlos=<value>     Mean LOS angle between the vertical and the Look Direction. Taken into account if los argument is None [default: 0]
--output_cube=<file> Name to the output cube corrected from GACOS [default: cube_gacos_flat]
--output_gacos=<file> Name to the output GACOS flatten cube  [default: depl_cumule-gacos]
--rmspixel=<file>     Path to the RMS map that gives an error for each pixels for the empirical estimations (e.g RMSpixel, output of invers_pixel)
--threshold_rms=<value>  Thresold on rmspixel for ramp estimation [default: 2.]   
--plot=<yes|no>      Display results [default: yes]
--lectfile=<path>    Path to the lect.in file. Simple text file containing width and length and number of images of the time series cube (output of invers_pixel). By default the program will try to find an .hdr file. [default: lect.in].
"""

import gdal
from osgeo import osr
import sys,os
import numpy as np
import scipy.optimize as opt
import scipy.linalg as lst
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy.lib.stride_tricks import as_strided
import subprocess
np.warnings.filterwarnings('ignore')

try:
    from nsbas import docopt
except:
    import docopt

from contextlib import contextmanager
from functools import wraps, partial
import multiprocessing, warnings, logging
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning) 

logging.basicConfig(level=logging.INFO,\
        format='line %(lineno)s -- %(levelname)s -- %(message)s')
logger = logging.getLogger('correct_ts_from_gacos.log')

##################################################################################
###  Extras functions and context maganers
##################################################################################

# Timer for all the functions
class ContextDecorator(object):
    def __call__(self, f):
        @wraps(f)
        def decorated(*args, **kwds):
            with self:
                try:
                    return f(*args, **kwds)
                except (KeyboardInterrupt, SystemExit):
                    raise
                except:
                    Exception('{0} Failed !'.format(f))
                    raise
        return decorated

def checkinfile(file):
    if os.path.exists(file) is False:
        logger.critical("File: {0} not found, Exit !".format(file))
        print("File: {0} not found in {1}, Exit !".format(file,os.getcwd()))

# create generator for pool
@contextmanager
def poolcontext(*arg, **kargs):
    pool = multiprocessing.Pool(*arg, **kargs)
    yield pool
    pool.terminate()
    pool.join()

#####################################################################################
# FUNCTION
#####################################################################################


#####################################################################################
# END FUNCTION
#####################################################################################
    
# read arguments
arguments = docopt.docopt(__doc__)

if arguments["--path"] ==  None:
   arguments["--path"] = "./GACOS/"
if arguments["--cube"] ==  None:
    arguments["--cube"] = "depl_cumule"
if arguments["--output_cube"] ==  None:
    arguments["--output_cube"] = "depl_cumule-gacos"
if arguments["--list_images"] ==  None:
    arguments["--list_images"] = "list_images.txt"
if arguments["--lectfile"] ==  None:
    arguments["--lectfile"] = "lect.in"
if arguments["--output_gacos"] ==  None:
    arguments["--output_gacos"] = "cube_gacos_flat"
if arguments["--imref"] ==  None:
    imref = 0
else:
    imref = int(arguments["--imref"]) - 1
if arguments["--plot"] ==  None:
    arguments["--plot"] = 'yes'

if arguments["--proj"] ==  None:
    proj = False
    EPSG = 4326
else:
    proj = True
    EPSG = int(arguments["--proj"])

if arguments["--gacos2data"] ==  None:
    gacos2data = 1 
else:
    gacos2data = float(arguments["--gacos2data"])

if arguments["--ramp"] ==  None:
    arguments["--ramp"] = 'cst'
if (arguments["--ramp"] != 'cst') and (arguments["--ramp"] != 'lin'):
    arguments["--ramp"] = 'cst'

if arguments["--load"] ==  None or arguments["--load"] ==  'no':
    load = 'no'
    arguments["--load"] = 'cube_gacos'
else:
    load = 'yes'

if arguments["--topofile"] ==  None:
   radar = None
else:
   radar = arguments["--topofile"]
# read lect.in: size maps

if arguments["--crop"] ==  None:
   crop = False
else:
   crop = list(map(float,arguments["--crop"].replace(',',' ').split()))

nb,idates,dates,base=np.loadtxt(arguments["--list_images"], comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
checkinfile(arguments["--cube"])
ds = gdal.Open(arguments["--cube"])
if not ds:
  logger.info('.hdr file time series cube {0}, not found, open {1}'.format(arguments["--cube"],arguments["--lectfile"]))
  ncol, nlines = list(map(int, open(arguments["--lectfile"]).readline().split(None, 2)[0:2]))
  N = len(dates)
else:
  ncol, nlines = ds.RasterXSize, ds.RasterYSize
  N = ds.RasterCount
logger.info('nlines:{0} ncols:{1} N:{2}'.format(nlines,ncol,N))

if arguments["--los"] ==  None:
    if arguments["--meanlos"] ==  None:
        look = 0
    else:
        look = np.ones((nlines,ncol))*float(arguments["--meanlos"])
else:
    losf = arguments["--los"]
    extension = os.path.splitext(losf)[1]
    checkinfile(losf)
    if extension == ".tif":
        ds = gdal.Open(losf, gdal.GA_ReadOnly)
        band = ds.GetRasterBand(1)
        look = band.ReadAsArray()
        del ds
    else:
        fid = open(losf,'r')
        look = np.fromfile(fid,dtype=np.float32)[:nlines*ncol].reshape(nlines,ncol)
        fid.close()

if arguments["--crop_emp"] ==  None:
    empzone = [0,ncol,0,nlines]
    col_beg,col_end,line_beg,line_end = 0 , ncol, 0., nlines
else:
    empzone = list(map(float,arguments["--crop_emp"].replace(',',' ').split()))
    col_beg,col_end,line_beg,line_end = empzone[0], empzone[1], empzone[2], empzone[3]

if arguments["--ref_zone"] == None:
    refline_beg, refline_end, refcol_beg, refcol_end = None,None,None,None
else:
    ref = list(map(int,arguments["--ref_zone"].replace(',',' ').split()))
    refline_beg,refline_end, refcol_beg, refcol_end = ref[0], ref[1], ref[2], ref[3]

if arguments["--rmspixel"] ==  None:
    rms = np.ones((nlines,ncol))
else:
    rmsf = arguments["--rmspixel"]
    rms = np.fromfile(rmsf,dtype=np.float32).reshape((nlines,ncol))
    #plt.imshow(rms)
    #plt.show()

if arguments["--threshold_rms"] ==  None:
    threshold_rms = 2.
else:
    threshold_rms = float(arguments["--threshold_rms"])

if arguments["--fitmodel"] ==  None:
    arguments["--fitmodel"]  = 'no'

# Choose plot option
cmap = cm.gist_rainbow_r
cmap.set_bad('white')

if radar is not None:
    extension = os.path.splitext(radar)[1]
    if extension == ".tif":
      ds = gdal.Open(radar, gdal.GA_ReadOnly)
      band = ds.GetRasterBand(1)
      elev = band.ReadAsArray()
      del ds, band
    else:
      fid = open(radar,'r')
      elevi = np.fromfile(fid,dtype=np.float32)
      elevi = elevi[:nlines*ncol]
      elev = elevi.reshape((nlines,ncol))
      fid.close()
      del elevi
else:
    elev = np.ones((nlines,ncol))


# load cube of displacements
cubei = np.fromfile(arguments["--cube"],dtype=np.float32)
cube = as_strided(cubei[:nlines*ncol*N])
logger.info('Number of line in the cube: {0} '.format(cube.shape))
kk = np.flatnonzero(cube>9990)
cube[kk] = float('NaN')
maps = cube.reshape((nlines,ncol,N))
cst = np.copy(maps[:,:,imref])
for l in range((N)):
    maps[:,:,l] = maps[:,:,l]-cst
del cube, cubei

if load == 'no':
    gacos = np.zeros((nlines,ncol,N))
    for i in range((N)):

        logger.info('Read {0}-{1}'.format(idates[i], i))
        infile = arguments["--path"]+'{}.ztd'.format(int(idates[i]))
        rsc = arguments["--path"]+'{}.ztd.rsc'.format(int(idates[i]))
        temp1 = arguments["--path"]+'{}_temp1.tif'.format(int(idates[i]))
        outtif = arguments["--path"]+'{}_gacos.tif'.format(int(idates[i]))

        # read .rsc
        nncol = np.int(open(rsc).readlines()[0].split(None, 2)[1])
        nnlines = np.int(open(rsc).readlines()[1].split(None, 2)[1])
        xfirst = float(open(rsc).readlines()[6].split(None, 2)[1])
        yfirst = float(open(rsc).readlines()[7].split(None, 2)[1])
        xstep = float(open(rsc).readlines()[8].split(None, 2)[1])
        ystep = float(open(rsc).readlines()[9].split(None, 2)[1])
        geotransform = (xfirst, xstep, 0, yfirst, 0, ystep)
        logger.info('Read .rsc and set Geotransform: {}'.format(geotransform))
        ztd_ = np.fromfile(infile, dtype='float32').reshape(nnlines,nncol)

        # transform gacos format to samething humanily readeable
        driver = gdal.GetDriverByName('GTiff')
        ds = driver.Create(temp1, nncol, nnlines, 1, gdal.GDT_Float32)
        band = ds.GetRasterBand(1)
        band.WriteArray(ztd_)
        ds.SetGeoTransform(geotransform)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        ds.SetProjection(srs.ExportToWkt())
        band.FlushCache()
        del ds,band

        # crop, re-preoject and resample
        if proj and crop is not False:
            print("gdalwarp -overwrite -s_srs EPSG:4326 -t_srs EPSG:"+str(EPSG)+\
            " -ts "+str(ncol)+" "+str(nlines)+" -r average -te "+str(crop[0])+" "+str(crop[1])+" "+str(crop[2])+" "+str(crop[3])+" "\
            +temp1+" "+outtif+" -of GTiff")
            r = subprocess.call("gdalwarp -overwrite -s_srs EPSG:4326 -t_srs EPSG:"+str(EPSG)+\
                " -te "+str(crop[0])+" "+str(crop[1])+" "+str(crop[2])+" "+str(crop[3])+" "\
                +" -ts "+str(ncol)+" "+str(nlines)+" -r average "+temp1+" "+outtif+" -of GTiff", shell=True)
            if r != 0:
                raise Exception("gdalwarp failed")

        elif crop is not False and proj is False:
            print("gdalwarp -overwrite -s_srs EPSG:4326-ts "+str(ncol)+" "+str(nlines)+" -r average \
            -te "+str(crop[0])+" "+str(crop[1])+" "+str(crop[2])+" "+str(crop[3])+" "+temp1+" "+outtif+" -of GTiff")
            r = subprocess.call("gdalwarp -overwrite -s_srs EPSG:4326 -ts "+str(ncol)+" "+str(nlines)+" -r average \
                -te "+str(crop[0])+" "+str(crop[1])+" "+str(crop[2])+" "+str(crop[3])+" "+temp1+" "+outtif+" -of GTiff", shell=True)
            if r != 0:
                raise Exception("gdalwarp failed")

        elif proj and crop is False:
            print("gdalwarp -overwrite -s_srs EPSG:4326 -t_srs EPSG:"+str(EPSG)+" -ts "+str(ncol)+" "+str(nlines)+" -r average "+temp1+" "+outtif+" -of GTiff")
            r = subprocess.call("gdalwarp -overwrite -s_srs EPSG:4326 -t_srs EPSG:"+str(EPSG)+" -ts "+str(ncol)+" "+str(nlines)+" -r average "+temp1+" "+outtif+" -of GTiff", shell=True)
            if r != 0:
                raise Exception("gdalwarp failed")

        else:
            print("gdalwarp -overwrite -ts "+str(ncol)+" "+str(nlines)+" -r average "+temp1+" "+outtif+" -of GTiff")
            r = subprocess.call("gdalwarp -overwrite -ts "+str(ncol)+" "+str(nlines)+" -r average "+temp1+" "+outtif+" -of GTiff", shell=True)
            if r != 0:
                raise Exception("gdalwarp failed")

        # del temporary files
        cmd = 'rm -f {}'.format(temp1)
        os.system(cmd)

        # open gacos
        ds = gdal.Open(arguments["--path"]+'{}_gacos.tif'.format(int(idates[i])), gdal.GA_ReadOnly)
        band = ds.GetRasterBand(1)
        gacos[:,:,i] = band.ReadAsArray()
        del ds,band

    # Ref atmo models to the reference image
    cst = np.copy(gacos[:,:,imref])
    for l in range((N)):
        d = as_strided(gacos[:,:,l])
        gacos[:,:,l] = gacos[:,:,l] - cst

    # save new cube
    fid = open(arguments["--load"], 'wb')
    gacos.flatten().astype('float32').tofile(fid)
    fid.close()

# # load gacos cube
checkinfile(arguments["--load"])
gcubei = np.fromfile(arguments["--load"],dtype=np.float32)
gcube = as_strided(gcubei[:nlines*ncol*N])
gacos = gcube.reshape((nlines,ncol,N))*gacos2data

cst = np.copy(gacos[:,:,imref])
for l in range((N)):
    gacos[:,:,l] = (gacos[:,:,l]-cst)/np.cos(np.deg2rad(look))
del gcubei, gcube

# Plot
fig = plt.figure(0,figsize=(14,10))
fig.subplots_adjust(wspace=0.001)
vmax = np.nanpercentile(gacos,99)
vmin = np.nanpercentile(gacos,1)
for l in range((N)):
   d = as_strided(gacos[:,:,l])
   ax = fig.add_subplot(4,int(N/4)+1,l+1)
   cax = ax.imshow(d,cmap=cmap,vmax=vmax,vmin=vmin)
   ax.set_title(idates[l],fontsize=6)
   plt.setp( ax.get_xticklabels(), visible=False)
   plt.setp( ax.get_yticklabels(), visible=False)

plt.suptitle('GACOS models')
fig.colorbar(cax, orientation='vertical',aspect=10)
fig.savefig('gocos.eps', format='EPS',dpi=150)
if arguments["--plot"] == 'yes':
   plt.show()
plt.close('all')

# Apply correction
maps_flat = np.zeros((nlines,ncol,N))
model_flat = np.zeros((nlines,ncol,N))
# initialise variance
var = np.zeros((N))

for l in range((N)):

  data = as_strided(maps[:,:,l]) 
  data_flat = as_strided(maps_flat[:,:,l])
  model = as_strided(gacos[:,:,l])
  
  # print(l, imref)
  if l == imref:
    # check if data and model are set to zeros
    data_flat[np.isnan(data)] = np.float('NaN')
    data_flat[data_flat>999.]= np.float('NaN')
  
  else:

    _los_map = np.copy(data)
    _los_map[np.logical_or(data==0,data>=9990)] = float('NaN')
    losmin,losmax = np.nanpercentile(_los_map,2.),np.nanpercentile(_los_map,98.)
    gacosmin,gacosmax = np.nanpercentile(model,2),np.nanpercentile(model,98)
    del _los_map

    if radar is not None:
        maxtopo,mintopo = np.nanpercentile(elev,98), np.nanpercentile(elev,2)
    else:
        maxtopo,mintopo = 1e8, -1e8

    # index = np.nonzero(data>2)
    # data[index] = float('NaN')
    # plt.imshow(data)
    # plt.show()
    # sys.exit()

    funct = 0.
    pix_lin, pix_col = np.indices((nlines,ncol))
   
    # find proportionality between data and model
    index = np.nonzero(
        np.logical_and(~np.isnan(data),
        np.logical_and(pix_lin>line_beg,
        np.logical_and(pix_lin<line_end,
        np.logical_and(data<losmax,
        np.logical_and(data>losmin,
        np.logical_and(rms<threshold_rms,
        np.logical_and(rms>1.e-6,
        np.logical_and(model<gacosmax,
        np.logical_and(model>gacosmin,
        np.logical_and(elev<maxtopo,
        np.logical_and(elev>mintopo,
        np.logical_and(data!=0.0,
        np.logical_and(model!=0.0,
        np.logical_and(pix_col>col_beg, pix_col<col_end
        )))))))))))))))

    if refline_beg is not None:
        indexref = np.nonzero(
        np.logical_and(~np.isnan(data),
        np.logical_and(pix_lin>refline_beg,
        np.logical_and(pix_lin<refline_end,
        np.logical_and(data<losmax,
        np.logical_and(data>losmin,
        np.logical_and(rms<threshold_rms,
        np.logical_and(rms>1.e-6,
        np.logical_and(model<gacosmax,
        np.logical_and(model>gacosmin,
        np.logical_and(elev<maxtopo,
        np.logical_and(elev>mintopo,
        np.logical_and(data!=0.0,
        np.logical_and(model!=0.0,
        np.logical_and(pix_col>refcol_beg, pix_col<refcol_end
        )))))))))))))))

        if len(indexref[0]) < 2:
            print('No data in the reference zone for date: %i. Exit!'%(idates[l]))
            sys.exit()

    temp = np.array(index).T
    x_clean = temp[:,0]; y_clean = temp[:,1]
    los_clean = data[index].flatten()
    model_clean = model[index].flatten()
    rms_clean = rms[index].flatten()
    del temp
   
    if refline_beg is not None:
        # compute average phase in the ref area 
        # we want los ref area to be zero
        los_ref = data[indexref].flatten() - model[indexref].flatten()
        rms_ref = rms[indexref].flatten()
        amp_ref = 1./rms_ref
        amp_ref = amp_ref/np.nanpercentile(amp_ref,99)
        # weigth average of the phase
        cst = np.nansum(los_ref*amp_ref) / np.nansum(amp_ref)
    else:
        # compute the constante on the whole image
        # would be a better idea to define a ref_zone
        los_ref = data[index].flatten() - model[index].flatten()
        rms_ref = rms[index].flatten()
        amp_ref = 1./rms_ref
        amp_ref = amp_ref/np.nanpercentile(amp_ref,99)
        # weigth average of the phase
        cst = np.nansum(los_ref*amp_ref) / np.nansum(amp_ref)

    # digitize data in bins, compute median and std
    bins = np.arange(gacosmin,gacosmax,abs(gacosmax-gacosmin)/500.)
    inds = np.digitize(model_clean,bins)
    modelbins = []
    losbins = []
    losstd = []
    xbins, ybins = [], []
    los_clean2, model_clean2, x_clean2, y_clean2, rms_clean2 = [], [], [], [], []
    for j in range(len(bins)-1):
            uu = np.flatnonzero(inds == j)
            if len(uu)>500:
                modelbins.append(bins[j] + (bins[j+1] - bins[j])/2.)

                # do a clean within the sliding median
                indice = np.flatnonzero(np.logical_and(los_clean[uu]>np.percentile(\
                    los_clean[uu],10.),los_clean[uu]<np.percentile(los_clean[uu],90.)))

                losstd.append(np.nanstd(los_clean[uu][indice]))
                losbins.append(np.nanmedian(los_clean[uu][indice]))
                xbins.append(np.nanmedian(x_clean[uu][indice]))
                ybins.append(np.nanmedian(y_clean[uu][indice]))

                # remove outliers from data
                los_clean2.append(los_clean[uu][indice])
                x_clean2.append(x_clean[uu][indice])
                y_clean2.append(y_clean[uu][indice])
                model_clean2.append(model_clean[uu][indice])
                # new rms is clean rms time the standard deviation of the los within the bin
                rms_clean2.append(rms_clean[uu][indice]*np.nanstd(los_clean[uu][indice]))

    losbins = np.array(losbins)
    losstd = np.array(losstd)
    modelbins = np.array(modelbins)
    xbins, ybins = np.array(xbins),np.array(ybins)

    los_clean = np.concatenate(los_clean2)
    x_clean, y_clean = np.concatenate(x_clean2), np.concatenate(y_clean2)
    model_clean = np.concatenate(model_clean2)
    rms_clean = np.concatenate(rms_clean2)
    del los_clean2, x_clean2, y_clean2, model_clean2, bins, inds

    if arguments["--fitmodel"]  == "no":
      
      if (arguments["--ramp"] == 'cst'):
        
        a = cst
        logger.info('Remove cst: %f for date: %i'%(a,idates[l]))

        remove_ramp = np.ones((nlines,ncol))*cst
        remove_ramp[model==0.] = 0.
        remove_ramp[np.isnan(data)] = float('NaN')
        
        # data = gacos + cst
        remove = remove_ramp + model
       
        # Compute ramp for data = f(model) 
        functbins = a
        funct = a
        # set coef gacos to 1
        coef_model = 1

      elif (arguments["--ramp"] == 'lin'):
        # here we want to minimize data-gacos = ramp
        # invers both clean data and ref frame together
        d = los_clean - model_clean
        sigmad = rms_clean  

        G=np.zeros((len(d),3))
        G[:,0] = x_clean
        G[:,1] = y_clean
        G[:,2] = 1

        x0 = lst.lstsq(G,d)[0]
        _func = lambda x: np.sum(((np.dot(G,x)-d)/sigmad)**2)
        _fprime = lambda x: 2*np.dot(G.T/sigmad, (np.dot(G,x)-d)/sigmad)
        pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=500,full_output=True,iprint=0,acc=1.e-9)[0]
        a = pars[0]; b = pars[1]; c = pars[2]
        logger.info('Remove ramp  %f az  + %f r + %f for date: %i'%(a,b,c,idates[l]))
            
        # Compute ramp for data = f(model)
        functbins = a*xbins + b*ybins + c
        funct = a*x_clean + b*y_clean + c
        # set coef gacos to 1
        coef_model = 1

        # build total G matrix
        G=np.zeros((len(data.flatten()),3))
        for i in range(nlines):
            G[i*ncol:(i+1)*ncol,0] = i - line_beg 
            G[i*ncol:(i+1)*ncol,1] = np.arange(ncol) - col_beg
        G[:,2] = 1

        # compute ramp
        remove_ramp = np.dot(G,pars).reshape(nlines,ncol)
        remove_ramp[model==0.] = 0.
        remove_ramp[np.isnan(data)] = float('NaN')
        
        # data = gacos + (a*rg + b*az + cst) 
        remove = model + remove_ramp

    elif arguments["--fitmodel"] =='yes':
      
      if arguments["--ramp"] == 'cst': 
        # invers both digitized data and ref frame together
        # here we invert data = a*gacos + cst
        d = losbins
        sigmad = losstd   
            
        G=np.zeros((len(d),2))
        G[:,0] = 1
        G[:,1] = modelbins

        x0 = lst.lstsq(G,d)[0]
        _func = lambda x: np.sum(((np.dot(G,x)-d)/sigmad)**2)
        _fprime = lambda x: 2*np.dot(G.T/sigmad, (np.dot(G,x)-d)/sigmad)
        pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=500,full_output=True,iprint=0,acc=1.e-9)[0]
        a = pars[0]; b = pars[1]
        logger.info('Remove ramp %f + %f model for date: %i'%(a,b,idates[l]))

        #Compute ramp for data = f(model)
        funct = a
        functbins = a
        coef_model = b

        # build total G matrix
        G=np.zeros((len(data.flatten()),2))
        G[:,0] = 1
        G[:,1] = model.flatten()

        # data = a*gacos + ramp
        remove = np.dot(G,pars).reshape(nlines,ncol)
        remove[model==0.] = 0.
        remove[np.isnan(data)] = float('NaN')

        # ramp only
        remove_ramp = np.ones((nlines,ncol))*a
        remove_ramp = np.dot(G[:,:-1],pars[:-1]).reshape(nlines,ncol)
        remove_ramp[model==0.] = 0.
        remove_ramp[np.isnan(data)] = float('NaN')
   
      elif arguments["--ramp"] == "lin": 
        # invers both clean data and ref frame together
        # here we invert data = a*gacos + ramp
        d = los_clean
        sigmad = rms_clean  
            
        G=np.zeros((len(d),6))
        G[:,0] = x_clean**2
        G[:,1] = x_clean
        G[:,2] = y_clean**2
        G[:,3] = y_clean
        G[:,4] = 1
        G[:,5] = model_clean

        x0 = lst.lstsq(G,d)[0]
        _func = lambda x: np.sum(((np.dot(G,x)-d)/sigmad)**2)
        _fprime = lambda x: 2*np.dot(G.T/sigmad, (np.dot(G,x)-d)/sigmad)
        pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=500,full_output=True,iprint=0)[0]
        a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
        logger.info('Remove ramp %f az**2, %f az  + %f r**2 + %f r + %f + %f model for date: %i'%(a,b,c,d,e,f,idates[l]))

        #Compute ramp for data = f(model)
        funct = a*x_clean**2 + b*x_clean + c*y_clean**2 + d*y_clean + e
        functbins = a*xbins**2 + b*xbins + c*ybins**2 + d*ybins + e
        coef_model = f        

        # build total G matrix
        G=np.zeros((len(data.flatten()),6))
        for i in range(nlines):
                G[i*ncol:(i+1)*ncol,0] = (i - line_beg)**2
                G[i*ncol:(i+1)*ncol,1] = i - line_beg
                G[i*ncol:(i+1)*ncol,2] = np.arange((ncol - col_beg))**2
                G[i*ncol:(i+1)*ncol,3] = np.arange(ncol - col_beg)
        G[:,4] = 1
        G[:,5] = model.flatten()

        # data = a*gacos + ramp
        remove = np.dot(G,pars).reshape(nlines,ncol)
        remove[model==0.] = 0.
        remove[np.isnan(data)] = float('NaN')

        # ramp only
        remove_ramp = np.dot(G[:,:-1],pars[:-1]).reshape(nlines,ncol)
        remove_ramp[model==0.] = 0.
        remove_ramp[np.isnan(data)] = float('NaN')
   
    ###################
    ### Aply correction

    # data_flat = data - (coef_model*model + remove_ramp)    
    data_flat[:,:] = data - remove
    data_flat[np.isnan(data)] = float('NaN')
    data_flat[data_flat>999.]= float('NaN')
    # model_flat = coef_model*model + remove_ramp
    model_flat[:,:,l] = coef_model*gacos[:,:,l] + remove_ramp
 
    if refline_beg is not None:
        # Refer data again to zero in ref area 
        rms_ref = rms[indexref].flatten()
        amp_ref = 1./rms_ref
        amp_ref = amp_ref/np.nanmax(amp_ref)
        cst = np.nansum(data_flat[indexref]*amp_ref) / np.nansum(amp_ref)
        logger.info('Ref area set to zero: lines: {0}-{1}, cols: {2}-{3}'.format(refline_beg,refline_end,refcol_beg,refcol_end))
        logger.info('Average phase within ref area: {0}'.format(cst))
        data_flat = data_flat - cst
        data_flat[np.isnan(data)] = float('NaN')
        data_flat[data_flat>999.]= float('NaN')
    
    # compute variance flatten data
    # var[l] = np.sqrt(np.nanmean(data_flat**2))
    var[l] = np.nanstd(data_flat)
    logger.info('Var: {0} '.format(var[l]))
    
    fig = plt.figure(1,figsize=(12,6))
    ax = fig.add_subplot(2,2,1)
    im = ax.imshow(data,cmap=cmap,vmax=losmax,vmin=losmin)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.set_title('Data {}'.format(idates[l]),fontsize=6)

    ax = fig.add_subplot(2,2,2)
    im = ax.imshow(model,cmap=cmap)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.set_title('Model {}'.format(idates[l]),fontsize=6)
    
    ax = fig.add_subplot(2,2,3)
    im = ax.imshow(model_flat[:,:,l],vmax=losmax,vmin=losmin,cmap=cmap)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.set_title('Flatten Model {}'.format(idates[l]),fontsize=6)

    ax = fig.add_subplot(2,2,4)
    im = ax.imshow(data_flat,cmap=cmap,vmax=losmax,vmin=losmin)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.set_title('Correct Data {}'.format(idates[l]),fontsize=6)
    fig.tight_layout()
    fig.savefig('{}-data-model-maps.png'.format(idates[l]), format='PNG',dpi=150)

    fig = plt.figure(2,figsize=(7,5))
    ax = fig.add_subplot(1,1,1)
    g = np.linspace(np.nanmax(model_clean),np.nanmin(model_clean),100)
    ax.scatter(model_clean,los_clean - funct, s=0.005, alpha=0.1, rasterized=True)
    ax.plot(modelbins,losbins - functbins,'-r', lw =.5)
    ax.plot(g,coef_model*g,'-r', lw =4.,label='flatten data = {0:.4f}*model'.format(coef_model))
    ax.set_ylim([(los_clean-funct).min(),(los_clean-funct).max()])
    ax.set_xlim([model_clean.min(),model_clean.max()])
    ax.set_xlabel('GACOS ZTD')
    ax.set_ylabel('LOS delay - RAMP')
    ax.set_title('Data/Model')
    plt.legend(loc='best')
    fig.tight_layout()

    try:
        fig.savefig('{}-data-model.png'.format(idates[l]), format='PNG',dpi=150)
    except:
        pass
        
    if arguments["--plot"] == 'yes':
        plt.show()
        # sys.exit()

    plt.close('all')

# save new cube
logger.info('Create flatten time series cube: {}'.format(arguments["--output_cube"]))
fid = open(arguments["--output_cube"], 'wb')
maps_flat.flatten().astype('float32').tofile(fid)
fid.close()

# save gacos ref spatialy
logger.info('Create flatten time series model: {}'.format(arguments["--output_gacos"]))
fid = open(arguments["--output_gacos"], 'wb')
model_flat.flatten().astype('float32').tofile(fid)
fid.close()

# save rms
logger.info('Save variance for each dates in rms_gacos.txt')
np.savetxt('rms.txt', var, header='# date   |   RMS', fmt=('%.8f'))
logger.info('Plot variance in rms.eps')
fig_rms = plt.figure(4,figsize=(14,7))
plt.plot(dates, var)
fig_rms.savefig('rms.eps', format='EPS',dpi=150)

if arguments["--plot"] == 'yes':
    plt.show()

