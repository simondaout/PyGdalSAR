#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
prep_invers_pixel.py
========================

This script prepares a work directory and input files for invers_pixel.

Usage:
  prep_invers_pixel.py --int_path=<path> [--outputdir=<path>] [--int_list=<path>] [--dates_list=<path>] \
          [--format=<value>] [--prefix=<value>] [--suffix=<value>] 

Options:
  --int_path PATH     Absolute path to interferograms directory 
  --outputdir=<arg>   Output directory [default: ./ts]
  --dates_list PATH   Path to text file containing date,bp,bt,doppler_frq,date_dec [default: baseline.rsc]
  --int_list PATH     Text file containing list of interferograms dates in two colums, $data1 $date2 [default: interf_pair.rsc]
  --prefix=<value>    Prefix name $prefix$date1-$date2$suffix.unw [default: '']
  --suffix=<vaue>     Suffix name $prefix$date1-$date2$suffix.unw [default: '']
  --format VALUE       Format input files: ROI_PAC, GAMMA, GTIFF [default: ROI_PAC]
  -h --help           Show this screen

"""

from __future__ import print_function
import glob, math, os, sys
import dateutil.parser
from osgeo import gdal
import numpy as np
import docopt
import gamma as gm
import gdal
gdal.UseExceptions()
import shutil

arguments = docopt.docopt(__doc__)
int_path=arguments["--int_path"]
if arguments["--int_list"] == None:
    int_list = 'interf_pair.rsc'
else:
    int_list=arguments["--int_list"]
if arguments["--dates_list"] == None:
    baseline = 'baseline.rsc'
else:
    baseline=arguments["--dates_list"]
if arguments["--prefix"] == None:
  prefix = str(object='')
else:
  prefix=arguments["--prefix"]

if arguments["--suffix"] == None:
  suffix = str(object='')
else:
  suffix=arguments["--suffix"]

print(prefix,suffix)

if arguments["--format"] ==  None:
    sformat = 'ROI_PAC'
else:
    sformat = arguments["--format"]

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

# Create output directories
makedirs(tsdir)

intdir = os.path.join(arguments["--int_path"]) + '/'

# read int
date_1,date_2=np.loadtxt(int_list,comments="#",unpack=True,dtype='i,i')
kmax=len(date_1)
print("number of interferogram: ",kmax)
# open baseline.rsc
source2=file(baseline,'r')
im,bp,bt,imd=np.loadtxt(source2,comments="#",usecols=(0,1,2,4),unpack=True,dtype='i,f,f,f')
print("image list=",baseline)
nmax=len(imd)
print("number of image: ",nmax)

# Now, write list_pair
wf = open(os.path.join(tsdir, "list_dates"), "w")
for i in xrange((nmax)):
    wf.write("%i %.6f %.6f %.6f\n" % (im[i], imd[i], bt[i], bp[i]))
wf.close()

# Create lndatadir variables
lndatadir = os.path.join(tsdir, "LN_DATA/")
makedirs(lndatadir)

# create list_pair
shutil.copy(int_list,os.path.join(tsdir, "list_pair"))


if sformat == 'ROI_PAC':
  for kk in xrange((kmax)):
      date1, date2 = date_1[kk], date_2[kk]
      idate = str(date1) + '-' + str(date2)
      folder = int_path + 'int_'+ str(date1) + '_' + str(date2) + '/'
      rscfile='../../'+folder + prefix + str(date1) + '-' + str(date2) + suffix +  '.unw.rsc'
      outrsc= lndatadir + str(date1)  + '-' + str(date2) + '_pre_inv.unw.rsc' 
      os.symlink(rscfile,outrsc)    
      infile='../../'+folder + prefix + str(date1) + '-' + str(date2) + suffix + '.unw'
      outint=lndatadir + str(date1)  + '-' + str(date2) + '_pre_inv.unw'
      print(infile)

      os.symlink(infile,outint)

  # Last, but not least: write the input file
  f = open(os.path.join(tsdir, "input_inv_send"), "w")
  f.write("""\
0.0003 %  temporal smoothing weight, gamma liss **2 (if <0.0001, no smoothing)
0   %   mask pixels with large RMS misclosure  (y=0;n=1)
.8 %  threshold for the mask on RMS misclosure (in same unit as input files)
1  % range and azimuth downsampling (every n pixel)
3 % iterations to correct unwrapping errors (y:nb_of_iterations,n:0)
2 % iterations to weight pixels of interferograms with large residual? (y:nb_of_iterations,n:0)
2 % Scaling value for weighting residuals (1/(res**2+value**2)) (in same unit as input files) (Must be approximately equal to standard deviation on measurement noise)
2 % iterations to mask (tiny weight) pixels of interferograms with large residual? (y:nb_of_iterations,n:0)
1.2 % threshold on residual, defining clearly wrong values (in same unit as input files)
0    %   outliers elimination by the median (only if nsamp>1) ? (y=0,n=1)
list_dates
0    % sort by date (0) ou by another variable (1) ?
list_pair
0   % interferogram format (RMG : 0; R4 :1) (date1-date2_pre_inv.unw or date1-date2.r4)
3100.   %  include interferograms with bperp lower than maximal baseline
1   %Weight input interferograms by coherence or correlation maps ? (y:0,n:1)
0   %coherence file format (RMG : 0; R4 :1) (date1-date2.cor or date1-date2-CC.r4)
1   %   minimal number of interferams using each image
1     % interferograms weighting so that the weight per image is the same (y=0;n=1)
0.7 % maximum fraction of discarded interferograms
0 %  Would you like to restrict the area of inversion ?(y=1,n=0)
1 735 1500 1585  %Give four corners, lower, left, top, right in file pixel coord
1  %    referencing of interferograms by bands (1) or corners (2) ? (More or less obsolete)
5  %     band NW -SW(1), band SW- SE (2), band NW-NE (3), or average of three bands (4) or no referencement (5) ?
1   %   Weigthing by image quality (y:0,n:1) ? (then read quality in the list of input images)
1   %  Weigthing by interferogram variance (y:0,n:1) or user given weight (2)?
1    % use of covariance (y:0,n:1) ? (Obsolete)
0   % include a baseline term in inversion ? (y:1;n:0) Require to use smoothing option (smoothing coefficient) !
1   % smoothing by Laplacian, computed with a scheme at 3pts (0) or 5pts (1) ?
2   % weigthed smoothing by the average time step (y :0 ; n : 1, int : 2) ?
1 % put the first derivative to zero (y :0 ; n : 1)?
""")
  f.close()

else:
  for kk in xrange((kmax)):
      date1, date2 = date_1[kk], date_2[kk]
      idate = str(date1) + '_' + str(date2)
      if sformat == 'GTIFF':
        infile = int_path + '/' + prefix + str(date1) + '-' + str(date2) + suffix + '.tiff'
        ds = gdal.Open(infile, gdal.GA_ReadOnly)
        ds_band = ds.GetRasterBand(1)
        los = ds_band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        lines,cols =  ds.RasterYSize, ds.RasterXSize

      elif sformat == 'GAMMA':
        infile = int_path + '/' + prefix + str(date1) + '_' + str(date2) + suffix + '.unw'

        los = gm.readgamma(infile)
        los = np.float32(los)
        lines,cols = gm.readpar()

      # save output
      # outfile =  str(date1) + '-' + str(date2) + '.r4'
      # print('Create  {} file, lines:{}, cols:{}'.format(outfile,lines, cols)) 
      # fid = open(lndatadir +outfile,'wb')
      # los.flatten().astype('float32').tofile(fid)
      # fid.close()
      # outrsc= lndatadir + str(date1)  + '-' + str(date2) + '.r4.rsc' 
      # f = open(outrsc, "w")
      # f.write("""\
      # WIDTH                 %d
      # FILE_LENGTH           %d
      # XMIN                  0
      # XMAX                  %d
      # YMIN                  0
      # YMAX                  %d
      # """ % (cols, lines, cols-1, lines-1))
      # f.close()

      driver = gdal.GetDriverByName("roi_pac")
      outint=lndatadir + str(date1)  + '-' + str(date2) + '_pre_inv.unw'
      dst_ds = driver.Create(outint, cols, lines, 2, gdal.GDT_Float32)
      print('Create  {} file, lines:{}, cols:{}'.format(outint,lines, cols)) 
      dst_band1 = dst_ds.GetRasterBand(1)
      dst_band2 = dst_ds.GetRasterBand(2)
      dst_band1.WriteArray(np.ones((lines,cols)),0,0)
      dst_band2.WriteArray(los,0,0)

      outrsc= lndatadir + str(date1)  + '-' + str(date2) + '_pre_inv.unw.rsc' 
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


  # Last, but not least: write the input file
  f = open(os.path.join(tsdir, "input_inv_send"), "w")
  f.write("""\
0.0003 %  temporal smoothing weight, gamma liss **2 (if <0.0001, no smoothing)
0   %   mask pixels with large RMS misclosure  (y=0;n=1)
.8 %  threshold for the mask on RMS misclosure (in same unit as input files)
1  % range and azimuth downsampling (every n pixel)
3 % iterations to correct unwrapping errors (y:nb_of_iterations,n:0)
2 % iterations to weight pixels of interferograms with large residual? (y:nb_of_iterations,n:0)
2 % Scaling value for weighting residuals (1/(res**2+value**2)) (in same unit as input files) (Must be approximately equal to standard deviation on measurement noise)
2 % iterations to mask (tiny weight) pixels of interferograms with large residual? (y:nb_of_iterations,n:0)
1.2 % threshold on residual, defining clearly wrong values (in same unit as input files)
0    %   outliers elimination by the median (only if nsamp>1) ? (y=0,n=1)
list_dates
0    % sort by date (0) ou by another variable (1) ?
list_pair
0   % interferogram format (RMG : 0; R4 :1) (date1-date2_pre_inv.unw or date1-date2.r4)
3100.   %  include interferograms with bperp lower than maximal baseline
1   %Weight input interferograms by coherence or correlation maps ? (y:0,n:1)
0   %coherence file format (RMG : 0; R4 :1) (date1-date2.cor or date1-date2-CC.r4)
1   %   minimal number of interferams using each image
1     % interferograms weighting so that the weight per image is the same (y=0;n=1)
0.7 % maximum fraction of discarded interferograms
0 %  Would you like to restrict the area of inversion ?(y=1,n=0)
1 735 1500 1585  %Give four corners, lower, left, top, right in file pixel coord
1  %    referencing of interferograms by bands (1) or corners (2) ? (More or less obsolete)
5  %     band NW -SW(1), band SW- SE (2), band NW-NE (3), or average of three bands (4) or no referencement (5) ?
1   %   Weigthing by image quality (y:0,n:1) ? (then read quality in the list of input images)
1   %  Weigthing by interferogram variance (y:0,n:1) or user given weight (2)?
1    % use of covariance (y:0,n:1) ? (Obsolete)
0   % include a baseline term in inversion ? (y:1;n:0) Require to use smoothing option (smoothing coefficient) !
1   % smoothing by Laplacian, computed with a scheme at 3pts (0) or 5pts (1) ?
2   % weigthed smoothing by the average time step (y :0 ; n : 1, int : 2) ?
1 % put the first derivative to zero (y :0 ; n : 1)?
""")
  f.close()



