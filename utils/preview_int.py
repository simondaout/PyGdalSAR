#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
preview_int.py
========================
This script plot jpeg file for a list of interferograms

Usage:
  preview_int.py  --outputdir=<path>  --radar=<path> [--homedir=<path>] [--int_list=<path>] [--int_path=<path>] \
  [--prefix=<value>] [--suffix=<value>] [--rlook=<value>] [--dlook=<value>] 

Options:
  --outputdir PATH    Output directory where .jpeg are saved
  --radar PATH        Path to the radar.hgt file
  --homedir PATH      Path to home directory  [default: ./]
  --int_list PATH     Text file containing list of interferograms dates in two colums, $data1 $date2 [default: interf_pair.rsc]
  --int_path PATH     Path to interfeorgrams [default: int]
  --prefix=<value>    Prefix name ${prefix}$date1-$date2${suffix}_${rlook}.int [default: '']
  --suffix=<vaue>     Suffix name ${prefix}$date1-$date2${suffix}_${rlook}.int  [default: '']
  --rlook VALUE       Multilook number ${prefix}$date1-$date2${suffix}_${rlook}.int [default: 2]
  --dlook VALUE       downsample look of the output jpeg [default: 1]
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
radarf=os.path.join(homedir, arguments["--radar"])

if arguments["--int_path"] == None:
    int_path= os.path.join(homedir, 'int')
else:
    int_path=os.path.join(homedir, arguments["--int_path"]) + '/'


if arguments["--int_list"] == None:
  int_list = os.path.join(homedir,'interf_pair.rsc')
else:
  int_list = os.path.join(homedir,arguments["--int_list"])
if arguments["--prefix"] == None:
    prefix = ''
else:
    prefix=arguments["--prefix"]
if arguments["--suffix"] == None:
    suffix = ''
else:
    suffix=arguments["--suffix"]
if arguments["--rlook"] == None:
    rlook = 2
else:
    rlook = int(arguments["--rlook"])
if arguments["--dlook"] == None:
    dlook = 1
else:
    dlook = int(arguments["--dlook"])
outlook = str(int(rlook)+int(dlook))
rlook = str('_' + str(rlook) + 'rlks')

date_1,date_2=np.loadtxt(int_list,comments="#",unpack=True,usecols=(0,1),dtype='i,i')
kmax=len(date_1)
# cleanif 
if os.path.exists(outputdir):
  # print '{0}/{1}*{2}{3}.jpeg'.format(outputdir, prefix, suffix, rlook)
  jpeg_files = glob.glob('{0}/{1}*{2}{3}.jpeg'.format(outputdir, prefix, suffix, rlook))
  for f in jpeg_files:
    os.remove(f)
else:
  os.makedirs(outputdir)

for kk in xrange((kmax)):
    date1, date2 = date_1[kk], date_2[kk]
    idate = str(date1) + '-' + str(date2) 
    folder =  'int_'+ str(date1) + '_' + str(date2) + '/'
    infile = int_path + folder +  prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.int'
    jpeg = int_path + folder +  prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.jpeg'
    outjpeg = outputdir + prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.jpeg'

    try:
      r = subprocess.call("nsb_preview_int -croipac -l"+str(dlook)+" "+str(radarf)\
      +" "+str(infile)+" "+str(jpeg), shell=True)
      if r != 0:
            raise Exception("nsb_preview_int failed for date: ", infile)
      else:
            print 'Create: ', jpeg
    except:
      pass

# print '{0}/int_*/{1}*{2}{3}.jpeg'.format(int_path, prefix, suffix, rlook)
jpeg_files = glob.glob('{0}/int_*/{1}*{2}{3}.jpeg'.format(int_path, prefix, suffix, rlook))
print
print 'Move files into:', outputdir
for f in jpeg_files:
    shutil.move(f,outputdir)