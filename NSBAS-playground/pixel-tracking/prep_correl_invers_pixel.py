#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT 
################################################################################

"""
prep_correl_invers_pixel.py
========================

This script prepares a work directory and input files for invers_pixel.

Usage:
	prep_correl_invers_pixel.py [--correl_path=<path>] [--outputdir=<path>] [--correl_list=<path>] [--dates_list=<path>] [--do_coh=<yes/no>] [--coh_path=<value>] [--band=value] [--link=<yes/no>] 
        prep_correl_invers_pixel.py [--correl_path=<path>] [--outputdir=<path>] [--correl_list=<path>] [--dates_list=<path>] [--do_coh=<yes/no>] [--coh_path=<value>] [--sigma=<path>] [--link=<yes/no>]  
	prep_correl_invers_pixel.py [--correl_path=<path>] [--outputdir=<path>] [--correl_list=<path>] [--dates_list=<path>] [--do_coh=<yes/no>] [--coh_path=<value>] [--Bc=<values>] [--link=<yes/no>]

Options:
  --correl_path=<dir>    Absolute path to correlation directory [default: ./CORREL] 
  --outputdir=<dir>   Output directory [default: ./ts]
  --dates_list=<file>   Path to text file containing date,bp,bt,doppler_frq,date_dec [default: baseline.rsc]
  --correl_list=<file>     Text file containing list of correlation dates in two colums, $data1 $date2 [default: interf_pair.rsc]
  --do_coh<yes/no>    If yes, link coherence files used to weigth time series analysis
  --coh_path=<dir>    Absolute path to cohrence file directory  
  --sigma=<file>      Path to an uncertainty file for each interferograms (e.g rms_unwcor.txt created by invert_ramp_topo_unw.py) If not None, create a third column in list_pair file corresponding to int weight in the TS analysis  
  --Bc=<value>	      Critical temporal and perpendicular baselines for weigthing interferograms (eg. 0.5,100). Downweight small temporal baselines and large perpendicular baselines 
  --band=<value>      Band fro, the disparity map, 1: horizontal disparities, 2: vertical disparities [default: 1] 
  --link=<yes/no>     If no copy data in LN_DIR instead of doing links [default: yes]
  -h --help           Show this screen
"""

import glob, math, os, sys
import dateutil.parser
from osgeo import gdal
import numpy as np
import docopt
from parsers import gamma as gm
from osgeo import gdal
gdal.UseExceptions()
import shutil

arguments = docopt.docopt(__doc__)
if arguments["--correl_path"] == None:
    correl_path = './CORREL/'
else:
    correl_path=arguments["--correl_path"]
if arguments["--correl_list"] == None:
    correl_list = 'interf_pair.rsc'
else:
    correl_list=arguments["--correl_list"]
if arguments["--dates_list"] == None:
    baseline = 'baseline.rsc'
else:
    baseline=arguments["--dates_list"]

if arguments["--do_coh"] == None:
  do_coh = int(1)
else:
  do_coh = int(0) 

if arguments["--band"] == None:
  band = int(1)
else:
  band = int(arguments["--band"])

if arguments["--coh_path"] == None:
  coh_path = ''
else:
  coh_path = arguments["--coh_path"]

if arguments["--sigma"] == None:
    sigmaf = None
else:
    sigmaf=arguments["--sigma"]

if (arguments["--Bc"] == None):
   weight = None
else:
   bc = map(float,arguments["--Bc"].replace(',',' ').split())
   btc, bpc = bc[0], bc[1]

# A more predictable makedirs
def makedirs(name):
    if os.path.exists(name):
        return
    os.makedirs(name)

# Function to write config
def date_to_float(d):
    return (d.year + (d.month-1)/12.0 + (d.day-1)/365.0)

if arguments["--outputdir"] == None:
    tsdir = './ts/'
else:
    tsdir = os.path.join(arguments["--outputdir"]) + '/'

if arguments["--link"] == 'no':
    link = False
else:
    link = True

# Create output directories
makedirs(tsdir)

intdir = os.path.join(arguments["--correl_path"]) + '/'

# read int
date_1,date_2=np.loadtxt(correl_list,comments="#",unpack=True,dtype='i,i')
kmax=len(date_1)
print("number of interferogram: ",kmax)
# open baseline.rsc
im,bp,bt,imd=np.loadtxt(baseline,comments="#",usecols=(0,1,2,4),unpack=True,dtype='i,f,f,f')
print("image list=",baseline)
nmax=len(imd)
print("number of image: ",nmax)

# Now, write list_pair
wf = open(os.path.join(tsdir, "list_dates"), "w")
for i in range((nmax)):
    wf.write("%i %.6f %.6f %.6f\n" % (im[i], imd[i], bt[i], bp[i]))
wf.close()

# Create lndatadir variables
lndatadir = os.path.join(tsdir, "LN_DATA/")
makedirs(lndatadir)

# remove old links
unw_files = glob.glob('{0}/*_pre_inv.unw*'.format(lndatadir))
for f in unw_files:
   os.remove(f)	

# create list_pair
if (arguments["--sigma"] == None) &  (arguments["--Bc"] == None): 
    print(arguments["--sigma"])
    shutil.copy(correl_list,os.path.join(tsdir, "list_pair"))
    do_sig = int(1)
elif (arguments["--sigma"] != None) & (arguments["--Bc"] == None):
    bid,bid2,sigma = np.loadtxt(sigmaf,comments="#",unpack=True, dtype='i,i,f')
    # weight = 1./(sigma+0.001)
    # weight = np.exp(-sigma/np.percentile(sigma,80))
    weight = np.exp(-sigma)
    do_sig = int(2) # user given weigth
    # print (sigma, weight)
    if len(sigma) != kmax:
      w2 = []
      for j in range((kmax)):
        for i in range(len(sigma)):
          if (bid[i]==date_1[j]) and  (bid2[i]==date_2[j]):
              w2.append(weight[i])
      if len(w2) != kmax:
         print('Error: sigma file not the same size that the number of interferograms')
         sys.exit()
      weigth = np.array(w2)
    wf = open(os.path.join(tsdir, "list_pair"), "w")
    for i in range((kmax)):
        wf.write("%i %i %.6f\n" % (date_1[i], date_2[i], weight[i]))
    wf.close()
elif (arguments["--sigma"] == None) &  (arguments["--Bc"] != None):
     print('Weigth interferograms based on their baselines with Btc:{} and Bpc:{}'.format(btc,bpc))
     do_sig = int(2)
     weight=np.zeros((kmax))
     for i in range((kmax)):
        deltat = btc/(abs(bt[im==date_1[i]] - bt[im==date_2[i]]))
        deltap = (abs(bp[im==date_1[i]] - bt[im==date_2[i]]))/bpc
        weight[i] = (float(np.exp(-deltap)) + float(np.exp(-deltat)))/2
     wf = open(os.path.join(tsdir, "list_pair"), "w")
     for i in range((kmax)):
          wf.write("%i %i %.6f\n" % (date_1[i], date_2[i], weight[i]))
     wf.close()

cohformat = int(0)
iformat = int(1)
for kk in range((kmax)):
      date1, date2 = date_1[kk], date_2[kk]
      idate = str(date1) + '_' + str(date2)
      infile = os.path.abspath(correl_path + '/'  + idate +  '/corr-D.tif')
      ds = gdal.Open(infile, gdal.GA_ReadOnly)
      ds_band = ds.GetRasterBand(band)
      los = ds_band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
      lines,cols =  ds.RasterYSize, ds.RasterXSize
      del ds, ds_band

      if do_coh == 0:
            cohformat = int(1)
            cohfile = os.path.abspath(coh_path + '/' + prefix_coh + str(date1) + '_' + str(date2) + suffix_coh + '.cc.tif')
            ds = gdal.Open(cohfile, gdal.GA_ReadOnly)
            ds_band = ds.GetRasterBand(1)
            coh = ds_band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
            del ds, ds_band

      # save output
      outfile =  str(date1) + '-' + str(date2) + '.r4'
      print('Create  {} file, lines:{}, cols:{}'.format(outfile,lines, cols)) 
      fid = open(lndatadir +outfile,'wb')
      los.flatten().astype('float32').tofile(fid)
      fid.close()
      outrsc= lndatadir + str(date1)  + '-' + str(date2) + '.r4.rsc' 
      f = open(outrsc, "w")
      f.write("""\
      WIDTH                 %d
      FILE_LENGTH           %d
      XMIN                  0
      XMAX                  %d
      YMIN                  0
      YMAX                  %d
      """ % (cols, lines, cols-1, lines-1))
      f.close()

      if do_coh == 0:
        outcofile =  str(date1) + '-' + str(date2) + '-CC.r4'
        print('Create  {} file, lines:{}, cols:{}'.format(outcofile,lines, cols)) 
        fid = open(lndatadir +outcofile,'wb')
        coh.flatten().astype('float32').tofile(fid)
        fid.close()
        outrsc= lndatadir + str(date1)  + '-' + str(date2) + '.r4.rsc' 
        f = open(outrsc, "w")
        f.write("""\
        WIDTH                 %d
        FILE_LENGTH           %d
        XMIN                  0
        XMAX                  %d
        YMIN                  0
        YMAX                  %d
        """ % (cols, lines, cols-1, lines-1))
        f.close()

# Write the input file
name = os.path.join(tsdir, "input_inv_send")
if os.path.exists(name) is False:
    f = open(os.path.join(tsdir, "input_inv_send"), "w")
    f.write("""\
0.003  #  temporal smoothing weight, gamma liss **2 (if <0.0001, no smoothing)
0     #   mask pixels with large RMS misclosure  (y=0;n=1)
1.7    #  threshold for the mask on RMS misclosure (in same unit as input files)
1      #  range and azimuth downsampling (every n pixel)
0      #  iterations to correct unwrapping errors (y:nb_of_iterations,n:0)
3      #  iterations to weight pixels of interferograms with large residual? (y:nb_of_iterations,n:0)
0.5    #  Scaling value for weighting residuals (1/(res**2+value**2)) (in same unit as input files) (Must be approximately equal to standard deviation on measurement noise)
2      #  iterations to mask (tiny weight) pixels of interferograms with large residual? (y:nb_of_iterations,n:0)
4.     #  threshold on residual, defining clearly wrong values (in same unit as input files)
1      #  outliers elimination by the median (only if nsamp>1) ? (y=0,n=1)
list_dates
0      #  sort by date (0) ou by another variable (1) ?
list_pair
%d     #  interferogram format (RMG : 0; R4 :1) (date1-date2_pre_inv.unw or date1-date2.r4)
3100.  #  include interferograms with bperp lower than maximal baseline
%d      #  Weight input interferograms by coherence or correlation maps ? (y:0,n:1)
%d      #  coherence file format (RMG : 0; R4 :1) (date1-date2.cor or date1-date2-CC.r4)
1      #  minimal number of interferams using each image
1      #  interferograms weighting so that the weight per image is the same (y=0;n=1)
0.5    #  maximum fraction of discarded interferograms
0      #  Would you like to restrict the area of inversion ?(y=1,n=0)
1 735 1500 1585  #  Give four corners, lower, left, top, right in file pixel coord
1      #  referencing of interferograms by bands (1) or corners (2) ? (More or less obsolete)
5      #  band NW -SW(1), band SW- SE (2), band NW-NE (3), or average of three bands (4) or no referencement (5) ?
1      #  Weigthing by image quality (y:0,n:1) ? (then read quality in the list of input images)
%d     #  Weigthing by interferogram variance (y:0,n:1) or user given weight (2)?
1      #  use of covariance (y:0,n:1) ? (Obsolete)
1      #  Adjust functions to phase history ? (y:1;n:0) Require to use smoothing option (smoothing coefficient) !
0      #  compute DEM error proportional to perpendicular baseline ? (y:1;n:0)
0 2003.0     #  include a step function ? (y:1;n:0)
0      #  include a cosinus / sinus function ? (y:1;n:0)
1      #  smoothing by Laplacian, computed with a scheme at 3pts (0) or 5pts (1) ?
2      #  weigthed smoothing by the average time step (y :0 ; n : 1, int : 2) ?
1      # put the first derivative to zero (y :0 ; n : 1)?
    """ % (iformat, do_coh, cohformat, do_sig))
    f.close()


