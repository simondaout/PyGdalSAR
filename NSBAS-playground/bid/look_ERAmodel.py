#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from os import environ
import numpy as np
import sys
import subprocess

def run(cmd):
    """
    Runs a shell command, and print it before running.
    """
    print(cmd)
    r = subprocess.call(cmd, shell=True, stdout=sys.stdout, stderr=subprocess.STDOUT,
        env=environ)
    if r != 0:
        logger.critical(r)
    return

# read lect.in 
# read dates
fimages='liste_date_ERA5'
dates,bid=np.loadtxt(fimages, comments='#', unpack=True,dtype='i,i')
N=len(dates)
inlook = 1
outlook = 2
rlook = int(outlook/inlook)

for d in dates:
    #infile = '{}_mdel_{}rlks.unw'.format(d,inlook)
    infile = '{}_mdel.unw'.format(d,inlook)
    run("look.pl "+str(infile)+" "+str(rlook))
