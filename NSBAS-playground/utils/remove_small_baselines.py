#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
remove_small_baselines.py
========================
Remove small baselines ifgs from a pair list file

Usage:
    remove_small_baselines.py --int_list=<path> --out_list=<path>  --Bt=<values>

Options:
  --int_list=<dir>    Absolute path to inut interferograms list  
  --out_list=<dir>    Absolute path to output interferograms list
  --Bt=<value>        Critical temporal and perpendicular baselines (eg. 12) 
  -h --help           Show this screen
"""


import docopt
import numpy as np
from datetime import datetime as datetimes
arguments = docopt.docopt(__doc__)

def date2dec(dates):
    dates  = np.atleast_1d(dates)
    for date in dates:
        x = datetimes.strptime('{}'.format(date),'%Y%m%d')
        dec = float(x.strftime('%j'))/365.1
        year = float(x.strftime('%Y'))
        times = year + dec
    return times

int_list=arguments["--int_list"]
out_list=arguments["--out_list"]
bt = float(arguments["--Bt"])

# read int
date_1,date_2=np.loadtxt(int_list,comments="#",unpack=True,dtype='i,i')
kmax=len(date_1)
print("number of interferogram: ",kmax)

# write new list_pair
k = 0 
wf = open(out_list, "w")
for i in range((kmax)):
    if date2dec(date_2[i]) - date2dec(date_1[i]) >= bt/365.1:
        wf.write("%i %i\n" % (date_1[i],date_2[i]))
        k = k+1
wf.close()
print("number of output interferograms: ",k)
