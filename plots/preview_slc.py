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

dates=np.loadtxt(dates_list,comments="#",unpack=True,usecols=(0),dtype='i')
kmax=len(dates)

# cleanif 
if os.path.exists(outputdir):
  # print ('{0}/{1}*{2}{3}.jpeg'.format(outputdir, prefix, suffix, rlook))
  jpeg_files = glob.glob('{0}/*_coreg*.jpeg'.format(outputdir))
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
    return r

@contextmanager
def poolcontext(*arg, **kargs):
    pool = multiprocessing.Pool(*arg, **kargs)
    yield pool
    pool.terminate()
    pool.join()

def dolook(kk):
  date = dates[kk]
  infile = str(date) + '/'+ str(date)+ '_coreg_'+str(look)+'rlks.slc'

  if os.path.exists(infile):
    pass
  else:
    temp = str(date) + '/'+str(date)+ '_coreg.slc'
    r= run("length.pl "+str(temp))
    r= run("look.pl "+str(temp)+" "+str(alook)+" "+str(look)+" > log_look.txt")

if look>1:
    work = [(kk) for kk in range(kmax)]
    with poolcontext(processes=nproc) as pool:
        results = pool.map(dolook, work)

try:
 os.remove(os.path.join(outputdir, "dates_problems.txt"))
 os.remove(os.path.join(outputdir, "dates_success.txt"))
except:
 pass
def preview(kk):
    successf = open(os.path.join(outputdir, "dates_success.txt"), "a")
    failf =  open(os.path.join(outputdir, "dates_problems.txt"), "a")
    date = dates[kk]
    try:
      if look>1:
          infile = str(date) + '/'+ str(date)+ '_coreg_'+str(alook)+'rlks.slc'
          jpeg = str(date) + '/'+ str(date)+ '_coreg_'+str(alook)+'rlks.jpeg'
	  r= run("length.pl "+str(infile))
          r = run("nsb_preview_slc "+str(infile)+" "+str(jpeg))
      else:
          infile = str(date) + '/'+ str(date)+ '_coreg'+'.slc'
	  r= run("length.pl "+str(infile))
          jpeg = str(date) + '/'+ str(date)+ '_coreg'+'.jpeg'
          r = run("nsb_preview_slc "+str(infile)+" "+str(jpeg)+" -l 20 -p 0.25")
      if r != 0:
            raise Exception("nsb_preview_slc failed for date: ", infile)
            failf.write("%s\n" % ( str(date)))
      else:
            print ('Create: ', jpeg)
            successf.write("%s\n" % ( str(date)))
    except:
        failf.write("%s\n" % ( str(date)))
    successf.close()
    failf.close()

work = [(kk) for kk in range(kmax)]
with poolcontext(processes=nproc) as pool:
    results = pool.map(preview, work)

# print ('{0}/int_*/{1}*{2}{3}.jpeg'.format(int_path, prefix, suffix, rlook))
jpeg_files = glob.glob('*/*_coreg*.jpeg')
print()
print ('Move files into:', outputdir)
for f in jpeg_files:
    shutil.move(f,outputdir)
