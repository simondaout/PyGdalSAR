#!/usr/bin/env python3
# -*- coding:utf-8 -*-

#__projet__ = "masque"
#__nom_fichier__ = "fonctions"
#__author__ = "Laure MANCEAU, Simon DAOUT"
#__date__ = "avril 2022"

"""\
gen_correl.py.py
----------------------
usage:
  gen_correl.py.py [-v] [-f] [--pairs=<path>] [--suffix=<values>] [--prefix=<value>] 
options:
  --pairs=<path>                    File that contains the pairs of images that will be correlated in two columns [Default: interf_pair.rsc]
  --suffix=<values>                 Suffix name $prefix$date1-$date2$suffix.slc [default: '_coreg']
  --prefix=<value>    Prefix name $prefix$date1-$date2$suffix.unw [default: ''] 
  -v                                Verbose mode. Show more information about the processing
  -h --help                         Show this screen.
"""

import numpy as np
from math import *
from osgeo import gdal
import os
import logging
import getopt
from sys import argv,exit,stdout
import docopt
import subprocess

def checkinfile(file):
    if os.path.exists(file) is False:
        logger.critical("File: {0} not found, Exit!".format(file))
        print("File: {0} not found in {1}, Exit!".format(file,os.os.getcwd()))

class Cd(object):
    def __init__(self,dirname):
        self.dirname = dirname
    def __enter__(self):
        self.curdir = os.getcwd()
        logger.debug('Enter {0}'.format(self.dirname))
        os.chdir(self.dirname)
    def __exit__(self, type, value, traceback):
        logger.debug('Exit {0}'.format(self.dirname))
        logger.debug('Enter {0}'.format(self.curdir))
        os.chdir(self.curdir)

def run(cmd):
    """
    Runs a shell command, and print it before running.

    Arguments:
        cmd: string to be passed to a shell

    Both stdout and stderr of the shell in which the command is run are those
    of the parent process.
    """

    logger.info(cmd)
    r = subprocess.call(cmd, shell=True, stdout=stdout, stderr=subprocess.STDOUT,
        env=os.environ)
    if r != 0:
        logger.critical(r)
    return

def correl(img1, img2):
    # asp command

    session="-t pinhole --prefilter-mode 2 --alignment-method none --prefilter-kernel-width 1.4 --nodata-value 0 "
    filtering="--filter-mode 2 --median-filter-size 3 --texture-smooth-size 13 --texture-smooth-scale 0.13  --rm-quantile-percentile 0.85 --rm-quantile-multiple 1 --rm-cleanup-passes 1 --rm-half-kernel 5 5 "
    correl="--corr-kernel 15 35 --cost-mode 2 --stereo-algorithm asp_bm --corr-tile-size 1024 --subpixel-mode 1 --subpixel-kernel 15 35 --corr-seed-mode 1 --threads-multiprocess  4 --xcorr-threshold 2 --min-xcorr-level 0 --sgm-collar-size 512"
    black_left =os.environ["PYGDALSAR"] +"/contrib/ASP/black_left.tsai"
    black_right=os.environ["PYGDALSAR"] + "/contrib/ASP/black_right.tsai"     

    run("parallel_stereo "+ session + " " +  str(img1) + " " + str(img2) + " " + black_left + " "  +  black_right + " corr --datum wgs84 --corr-timeout 300 "+ correl + " " + filtering)

def makedirs(name):
    if os.path.exists(name):
        return
    os.makedirs(name)
    logger.info('Create {} directory'.format(name))

def symlink(source, target):
    if os.path.exists(target):
        return
    os.symlink(source, target)
    logger.info('Create {} symlink'.format(source))

###########
#   MAIN  # 
###########

# Parse arguments: read arguments with dococpt
arguments = docopt.docopt(__doc__)
# init logger 
if arguments["-v"]:
    logging.basicConfig(level=logging.DEBUG,\
        format='%(asctime)s -- %(levelname)s -- %(message)s')
else:
    logging.basicConfig(level=logging.INFO,\
        format='%(asctime)s -- %(levelname)s -- %(message)s')
logger = logging.getLogger('gen_correl.log')

if arguments["--prefix"] == None:
  #prefix = str(object='')
  prefix = ''
else:
  prefix=arguments["--prefix"]
if arguments["--suffix"] == None:
  suffix = '_coreg'
else:
  suffix=arguments["--suffix"]

if arguments["--pairs"] == None:
    pairs = 'interf_pair.rsc'
else:
    pairs = arguments["--pairs"]
date_1,date_2=np.loadtxt(pairs,comments="#",unpack=True,dtype='i,i')
kmax=len(date_1)
print("number of pairs: ",kmax)

# Create CORREL DIR
correldir = os.path.join('./', "CORREL/")
makedirs(correldir)

for kk in range((kmax)):
    date1, date2 = date_1[kk], date_2[kk]
    idate = correldir + str(date1) + '_' + str(date2)
    makedirs(idate)
    folder1 = os.path.abspath(os.path.join('./', str(date1)))
    folder2 = os.path.abspath(os.path.join('./', str(date2)))
    print(folder1)
    with Cd(idate):    
        #récupération des images à corréler
        img1 = os.path.abspath(folder1 + '/' + str(date1) +  suffix +  '.tiff')
        target1 = os.path.abspath('./' + str(date1) +  suffix +  '.tiff')
        img2 = os.path.abspath(folder2 + '/' + str(date2) + suffix +  '.tiff')
        target2 = os.path.abspath('./' + str(date2) + suffix +  '.tiff')
        symlink(img1, target1)
        symlink(img2, target2)

        #réalisation de la corrélation
        outfile = 'corr-F.tif'  
        if os.path.exists(outfile):
          pass
        else: 
            correl(target1, target2)


