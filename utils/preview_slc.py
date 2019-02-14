#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
preview_int.py
========================
This script plot jpeg file for a list of interferograms

Usage:
  prevew_slc.py  --outputdir=<path>  [--homedir=<path>] [--dates_list=<path>] 


Options:
  --outputdir PATH    Output directory where .jpeg are saved
  --homedir PATH      Path to home directory  [default: ./]
  --dates_list PATH    Text file containing list of dates [default: baselines.rsc]
  -h --help           Show this screen
"""

import docopt
import numpy as np
import os, subprocess, glob, sys, shutil

# read arguments
arguments = docopt.docopt(__doc__)
if arguments["--homedir"] == None:
  homedir = os.path.abspath('.')
else:
  homedir=os.path.abspath(arguments["--homedir"])

outputdir=os.path.join(homedir, arguments["--outputdir"]) + '/'

if arguments["--dates_list"] == None:
  dates_list = os.path.join(homedir,'baselines.rsc')
else:
  dates_list = os.path.join(homedir,arguments["--dates_list"])

dates=np.loadtxt(dates_list,comments="#",unpack=True,usecols=(0),dtype='i')
kmax=len(dates)

# cleanif 
if os.path.exists(outputdir):
  # print '{0}/{1}*{2}{3}.jpeg'.format(outputdir, prefix, suffix, rlook)
  jpeg_files = glob.glob('{0}/*_coreg.jpeg'.format(outputdir))
  for f in jpeg_files:
    os.remove(f)
else:
  os.makedirs(outputdir)

for kk in xrange((kmax)):
    date = dates[kk]
    infile = str(date) + '/'+ str(date)+ '_coreg.slc'
    jpeg = str(date) + '/'+ str(date)+ '_coreg.jpeg'
    outjpeg = outputdir + str(date) + '/'+ str(date)+ '_coreg.jpeg'

    try:
      r = subprocess.call("nsb_preview_slc "+str(infile)+" "+str(jpeg), shell=True)
      if r != 0:
            raise Exception("nsb_preview_slc failed for date: ", infile)
      else:
            print 'Create: ', jpeg
    except:
      pass

# print '{0}/int_*/{1}*{2}{3}.jpeg'.format(int_path, prefix, suffix, rlook)
jpeg_files = glob.glob('*/*_coreg.jpeg')
print
print 'Move files into:', outputdir
for f in jpeg_files:
    shutil.move(f,outputdir)