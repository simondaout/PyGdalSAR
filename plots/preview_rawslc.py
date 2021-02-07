#!/usr/bin/env python3
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
from contextlib import contextmanager
from os import environ

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
  look = 1
else:
  look = int(arguments["--look"])
alook=look*5

if arguments["--nproc"] == None:
  nproc = 1
else:
  nproc = int(arguments["--nproc"])

dates,bid=np.loadtxt(dates_list,comments="#",unpack=True,usecols=(0,1),dtype='i,f')
kmax=len(dates)

# cleanif 
if os.path.exists(outputdir):
  # print ('{0}/{1}*{2}{3}.jpeg'.format(outputdir, prefix, suffix, rlook))
  jpeg_files = glob.glob('{0}/*.jpeg'.format(outputdir))
  for f in jpeg_files:
    os.remove(f)
else:
  os.makedirs(outputdir)

def run(cmd):
    print(cmd)
    r = subprocess.call(cmd, shell=True, stdout=sys.stdout, stderr=subprocess.STDOUT,
        env=environ)
    if r != 0:
        print(r)
    return

@contextmanager
def poolcontext(*arg, **kargs):
    pool = multiprocessing.Pool(*arg, **kargs)
    yield pool
    pool.terminate()
    pool.join()

def dolook(kk):
  date = dates[kk]
  infile = str(date) + '/'+ str(date)+ '_'+str(look)+'rlks.slc'

  if os.path.exists(infile):
    pass
  else:
    temp = str(date) + '/'+str(date)+ '.slc'
    run("look.pl "+str(temp)+" "+str(alook)+" "+str(look)+" > log_look.txt")

if look>1:
    work = [(kk) for kk in range(kmax)]
    with poolcontext(processes=nproc) as pool:
        results = pool.map(dolook, work)

def preview(kk):
    date = dates[kk]
    if look > 1:
        infile = str(date) + '/'+ str(date)+ '_'+str(alook)+'rlks.slc'
        jpeg = str(date) + '/'+ str(date)+ '_'+str(alook)+'rlks.jpeg'
    else:
        infile = str(date) + '/'+ str(date)+ '.slc'
        jpeg = str(date) + '/'+ str(date)+ '.jpeg'

    try:
      if look>1:
          run("nsb_preview_slc "+str(infile)+" "+str(jpeg))
      else:
          run("nsb_preview_slc "+str(infile)+" "+str(jpeg)+" -l 20 -p 0.25")
    except:
      pass

work = [(kk) for kk in range(kmax)]
with poolcontext(processes=nproc) as pool:
    results = pool.map(preview, work)

# print ('{0}/int_*/{1}*{2}{3}.jpeg'.format(int_path, prefix, suffix, rlook))
jpeg_files = glob.glob('*/*.jpeg')
print()
print ('Move files into:', outputdir)
for f in jpeg_files:
    shutil.move(f,outputdir)
