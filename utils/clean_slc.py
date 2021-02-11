#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-

"""\
clean_slc.py
-------------
Clean slc directory before running them again

Usage: check_coverage_unw.py --slc_list=<file> [--slc_path=<path>]

Options:
-h --help           Show this screen.
--slc_list=<file>   Text file containing list of slcs in one column
--slc_path=<path>   Path to input slc directories  [default ./] 
"""

from os import getcwd, path, chdir, remove
import docopt
import glob
import numpy as np

# create context manager for change dirt
class Cd(object):
    def __init__(self,dirname):
        self.dirname = dirname
    def __enter__(self):
        self.curdir = getcwd()
        print('Enter {0}'.format(self.dirname))
        chdir(self.dirname)
    def __exit__(self, type, value, traceback):
        print('Exit {0}'.format(self.dirname))
        print('Enter {0}'.format(self.curdir))
        chdir(self.curdir)

def rm(listf):
    for f in listf:
        if path.exists(f):
            remove(f)

# read arguments
arguments = docopt.docopt(__doc__)
slc_list=arguments["--slc_list"]
if arguments["--slc_path"] == None:
    slc_path='./'
else:
    slc_path=arguments["--slc_path"] + '/'

print('----------------------------------')
print('Check slc list:', slc_list)
dates=np.loadtxt(slc_list,comments="#",unpack=True,usecols=(0),dtype='S8')
Nd=len(dates)
print("number of dates to clean: {}".format(Nd))
print()

print('----------------------------------')
print('Check files: *slc*, *ampcor*,*cull*,*resamp*,fit_offsets*,*cor_topo_off*,*coreg.txt')

for j in range(Nd):
    print()
    dirname = dates[j]
    with Cd(dirname):
       #slc_files = glob.glob('*20rlks.slc*')
       slc_files = glob.glob('*slc*')
       print(slc_files)
       rm(slc_files)	
       coreg_files = glob.glob('*ampcor*')
       print(coreg_files)
       coreg_files = glob.glob('*cull*')
       print(coreg_files)
       coreg_files = glob.glob('*resamp*')
       print(coreg_files)
       coreg_files = glob.glob('fit_offsets*')
       print(coreg_files)
       coreg_files = glob.glob('*cor_topo_off*')
       print(coreg_files)
       coreg_files = glob.glob('*coreg.txt')
       print(coreg_files)

print('----------------------------------')
print()
