#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
remove_from_list.py
========================

This script prepares a work directory and input files for invers_pixel.

Usage:
    remove_from_list.py --int_list=<path> --out_list=<path>  --remove_list=<path>

Options:
  --int_list=<dir>    Absolute path to input interferograms list  
  --out_list=<dir>    Absolute path to output interferograms list
  --Bt=<value>        Absolute path to interferograms list to beremoved
  -h --help           Show this screen
"""


import docopt
import numpy as np
from datetime import datetime as datetimes
arguments = docopt.docopt(__doc__)

int_list=arguments["--int_list"]
out_list=arguments["--out_list"]
remove_list=arguments["--remove_list"]

# read int
date_1,date_2=np.loadtxt(int_list,comments="#",unpack=True,dtype='i,i')
date_11,date_22=np.loadtxt(remove_list,comments="#",unpack=True,dtype='i,i')
kmax=len(date_1)
kkmax=len(date_11)
print("Input number of interferogram: ",kmax)

# write new list_pair
k = 0 
wf = open(out_list, "w")
for i in range((kmax)):
    write=False
    for j in range((kkmax)):
        if (date_1[i] == date_11[j]) & (date_2[i] == date_22[j]):
            write=True
    if write:    
        wf.write("%i %i\n" % (date_1[i],date_2[i]))
        k = k+1
wf.close()
print("Number of output interferograms: ",k)
