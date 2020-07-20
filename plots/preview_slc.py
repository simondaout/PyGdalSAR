#!/usr/bin/env python2
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT 
################################################################################

"""
preview_int.py
========================
This script plot jpeg file for a list of interferograms

Usage:
  prevew_slc.py  --outputdir=<path>  [--homedir=<path>] [--dates_list=<path>] [--pixel_ratio=<value>] [--look=<value>] [--nproc=<value>] 


Options:
  --outputdir PATH    Output directory where .jpeg are saved
  --homedir PATH      Path to home directory  [default: ./]
  --dates_list PATH   Text file containing list of dates [default: baseline.rsc]
  --pixel_ratio Value Pixel ratio [default: 0.25]
  --look value        Look value [default: 4]
  --nproc Value       Number of processor [default: 4]
  -h --help           Show this screen
"""

import docopt
import numpy as np
import os, subprocess, glob, sys, shutil
import multiprocessing

# read arguments
arguments = docopt.docopt(__doc__)
if arguments["--homedir"] == None:
  homedir = os.path.abspath('.')
else:
  homedir=os.path.abspath(arguments["--homedir"])

outputdir=os.path.join(homedir, arguments["--outputdir"]) + '/'

if arguments["--dates_list"] == None:
  dates_list = os.path.join(homedir,'baseline.rsc')
else:
  dates_list = os.path.join(homedir,arguments["--dates_list"])

if arguments["--pixel_ratio"] == None:
  pixel_ratio = 0.25
else:
  pixel_ratio = float(arguments["--pixel_ratio"])

if arguments["--look"] == None:
  look = 4
else:
  look = int(arguments["--look"])

if arguments["--nproc"] == None:
  nproc = 1
else:
  nproc = int(arguments["--nproc"])

dates,bid=np.loadtxt(dates_list,comments="#",unpack=True,usecols=(0,1),dtype='i,f')
kmax=len(dates)

# cleanif 
if os.path.exists(outputdir):
  # print '{0}/{1}*{2}{3}.jpeg'.format(outputdir, prefix, suffix, rlook)
  jpeg_files = glob.glob('{0}/*_coreg*.jpeg'.format(outputdir))
  for f in jpeg_files:
    os.remove(f)
else:
  os.makedirs(outputdir)

def dolook(kk):
  date = dates[kk]
  infile = str(date) + '/'+ str(date)+ '_coreg_'+str(look)+'rlks.slc'

  if os.path.exists(infile):
    pass
  else:
    temp = str(date) + '/'+str(date)+ '_coreg.slc'
    os.system("look.pl "+str(temp)+" "+str(look*5)+" "+str(look))

pool = multiprocessing.Pool(nproc)
work = [(kk) for kk in xrange(kmax)]
pool.map(dolook, work)

def preview(kk):
    date = dates[kk]
    infile = str(date) + '/'+ str(date)+ '_coreg_'+str(look*5)+'rlks.slc'
    jpeg = str(date) + '/'+ str(date)+ '_coreg_'+str(look*5)+'rlks.jpeg'

    try:
      r = subprocess.call("nsb_preview_slc "+str(infile)+" "+str(jpeg), shell=True)
      if r != 0:
            raise Exception("nsb_preview_slc failed for date: ", infile)
      else:
            print 'Create: ', jpeg
    except:
      pass


pool = multiprocessing.Pool(nproc)
work = [(kk) for kk in xrange(kmax)]
pool.map(preview, work)

# print '{0}/int_*/{1}*{2}{3}.jpeg'.format(int_path, prefix, suffix, rlook)
jpeg_files = glob.glob('*/*_coreg*.jpeg')
print
print 'Move files into:', outputdir
for f in jpeg_files:
    shutil.move(f,outputdir)
