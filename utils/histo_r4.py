#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Simon DAOUT (Oxford)
############################################

"""\
histo_r4.py
-------------
Compute histogram distribution of a real4 file

Usage: histo_r4.py --infile=<path> --coeffile=<path> --threshold=<value> 
histo_r4.py -h | --help

Options:
-h --help           Show this screen
--infile PATH       File for histogram computation
--coeffile PATH     Coef to compare
--threshold VALUE   Value threshold for histo
"""

import numpy as np
import docopt
import sys
import matplotlib.pyplot as plt
import scipy.stats as stat

# read arguments
arguments = docopt.docopt(__doc__)
infile = arguments["--infile"]
coeffile = arguments["--coeffile"]
threshold = float(arguments["--threshold"])

los = np.fromfile(infile,dtype=np.float32)
coef = np.fromfile(coeffile,dtype=np.float32)
if len(los) != len(coef):
    print('lenght {}: {} =! lenght {}:{}'.format(infile,len(los),coeffile,len(coef)))
    sys.exit()

indexnan = np.flatnonzero(~np.isnan(los))
los,coef = los[indexnan],coef[indexnan]
indexnan = np.flatnonzero(~np.isnan(coef))
los,coef = los[indexnan],coef[indexnan]

indicemin = np.flatnonzero(los<threshold)
indicemax = np.flatnonzero(los>threshold)
losmin,losmax = los[indicemin],los[indicemax]
coefmin,coefmax = coef[indicemin], coef[indicemax]

fig = plt.figure(4, figsize=(14,6))

ax1 = fig.add_subplot(1,2,1)
xmin,xmax = np.nanmin(los),np.nanmax(los)
histo = ax1.hist(los,range=(xmin,xmax),bins=50,histtype='step', color='black')

traces = [losmin,losmax]
labels = ['{} < {}'.format(infile,threshold),'{} > {}'.format(infile,threshold)]
colors = ['red','blue']
for trace,label,color in zip(traces,labels,colors):
    histo = ax1.hist(trace,range=(xmin,xmax),bins=50,histtype='step', color=color, label=label)
    #plt.axvline(x=trace.mean(), c=color)
plt.xlabel(infile)
plt.legend()

ax2 = fig.add_subplot(1,2,2)
xmin,xmax = np.nanmin(coef)-np.nanstd(coef),np.nanmax(coef)+np.nanstd(coef)
#histo = ax2.hist(coef,normed=True,bins=50,range=(xmin,xmax),histtype='step', color='black')

traces = [coefmin,coefmax]
labels = ['{} < {}'.format(infile,threshold),'{} > {}'.format(infile,threshold)]
for trace,label,color in zip(traces,labels,colors):
    histo = ax2.hist(trace,bins=50,range=(xmin,xmax),histtype='step', color=color, label=label)
    #plt.axvline(x=stat.moment(trace,moment=3), c=color)
    plt.axvline(x=trace.mean(), c=color)
plt.xlabel(coeffile)
plt.legend()
plt.show()
