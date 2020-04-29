#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
import sys

def quad(x,a0,a1):
    return a0+a1*x

def seas(x,a0,a1,a2):
    return a0+a1*np.cos(2*np.pi*x)+a2*np.sin(2*np.pi*x)

source=file(sys.argv[1],'r')
date,vect = np.loadtxt(source,comments="#",usecols=(0,1),unpack=True,dtype='f,f')

fig = plt.figure(10, figsize=(8,8))
fig.subplots_adjust(hspace=0.2)

# 1) temp behavior
xx = np.arange(date.min(),date.max(),0.01)

def invert_quad(x,y,ax):
    pars, cova = sp.curve_fit(quad, x, y)
    perr = 2*np.sqrt(np.diag(cova))
    reg = quad(xx,pars[0],pars[1]) 
    ax.scatter(x,y,marker='x',c='black')
    ax.plot(xx,quad(xx,pars[0],pars[1]),color='red',lw=2.,label='y = {:2.4f}x'.format(pars[1]))
    # remove quad term to the vect
    ax.set_ylim([ymin,ymax])
    out = y - quad(x,pars[0],pars[1])
    plt.legend(loc='best')
    return out

ymin,ymax = np.min(vect), np.max(vect) 

ax1 = fig.add_subplot(2,1,1)
#ax1.scatter(date,vect,marker='o')
vect1 = invert_quad(date,vect,ax1)

# 2) inversion
dd = np.mod(date,1)
xx = np.arange(dd.min(),dd.max(),0.01)

def invert_seas(x,y,ax):
    pars, cova = sp.curve_fit(seas, x, y)

# 2) inversion
dd = np.mod(date,1)
xx = np.arange(dd.min(),dd.max(),0.01)

def invert_seas(x,y,ax):
    pars, cova = sp.curve_fit(seas, x, y)
    perr = 2*np.sqrt(np.diag(cova))
    parsmax,parsmin = pars+perr, pars-perr
    reg, regmin,regmax = seas(xx,pars[0],pars[1],pars[2]), seas(xx,parsmin[0],parsmin[1],parsmin[2]),seas(xx,parsmax[0],parsmax[1],parsmax[2]) 
    ax.scatter(x,y,marker='x',color='black')
    ax.plot(xx,seas(xx,pars[0],pars[1],pars[2]),color='red',lw=2,label='y = {:0.2f}cos(wt)+{:0.2f}sin(wt)'.format(pars[1],pars[2]))
    #ax.fill_between(xx,seas(xx,parsmin[0],parsmin[1],parsmin[2]),seas(xx,parsmax[0],parsmax[1],parsmax[2]),facecolor='#CC4F1B',alpha=0.5,linewidth=0)
    ax.set_ylim([ymin,ymax])
    ax.set_xlim(0,1)
    plt.legend(loc='best')

ax2 = fig.add_subplot(2,1,2)
invert_seas(dd,vect,ax2)

fig.savefig('acp_temporal.eps', format='EPS',dpi=150)
plt.show()
