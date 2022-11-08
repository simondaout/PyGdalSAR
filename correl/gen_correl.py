#!/usr/bin/env python3
# -*- coding:utf-8 -*-

#__projet__ = "masque"
#__nom_fichier__ = "fonctions"
#__author__ = "Laure MANCEAU"
#__date__ = "avril 2022"

"""\
correl.py.py
----------------------
usage:
  correl.py.py [-v] [-f] [--pairs=<path>] [--suffix=<values>] [--prefix=<value>] 
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

    session="-t pinhole --prefilter-mode 0 --alignment-method none --corr-sub-seed-percent 3 "
    filtering="--rm-cleanup-passes 1 --filter-mode 2 --rm-threshold 1.5 --rm-min-matches 50 --rm-half-kernel 9 9 --edge-buffer-size -1 "
    correl="--xcorr-threshold -1 --min-xcorr-level 1 --corr-search -4 4 -4 4 "

    run("parallel_stereo "+ session +  str(img1) + " " + str(img2) + " black_left.tsai black_right.tsai corr --datum wgs84 --corr-timeout 500 --stereo-algorithm asp_mgm --corr-kernel 5 5 --force-reuse-match-files --cost-mode 3 --subpixel-mode 2 --subpixel-kernel 21 21 --threads-multiprocess 10 --processes 1 " + filtering + " " + correl)

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
        symlink('/data/work/nepal/mathillo/processing/asp/black_left.tsai', 'black_left.tsai')
        symlink('/data/work/nepal/mathillo/processing/asp/black_right.tsai', 'black_right.tsai')
        symlink('/data/work/nepal/mathillo/processing/asp/stereo.default', 'stereo.default')

        #réalisation de la corrélation
        correl(target1, target2)


