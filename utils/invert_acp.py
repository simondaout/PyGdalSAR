#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp

def quad(x,a0,a1):
    return a0+a1*x

def seas(x,a0,a1,a2):
    return a0+a1*np.cos(2*np.pi*x)+a2*np.sin(2*np.pi*x)

source=file('acp_vecteurs','r')
date,vect1,vect2,vect3 = np.loadtxt(source,comments="#",usecols=(0,1,2,3),unpack=True,dtype='f,f,f,f')

fig = plt.figure(10, figsize=(15,8))
fig.subplots_adjust(hspace=0.2)

# 1) temp behavior
xx = np.arange(date.min(),date.max(),0.01)

def invert_quad(x,y,ax):
    pars, cova = sp.curve_fit(quad, x, y)
    perr = 2*np.sqrt(np.diag(cova))
    reg = quad(xx,pars[0],pars[1]) 
    ax.scatter(x,y,marker='o')
    ax.plot(xx,quad(xx,pars[0],pars[1]),label='y = {:2.4f}x'.format(pars[1]))
    # remove quad term to the vect
    ax.set_ylim([ymin,ymax])
    out = y - quad(x,pars[0],pars[1])
    plt.legend(loc='best')
    return out

ymin,ymax = np.min([vect1,vect2,vect3]), np.max([vect1,vect2,vect3]) 

ax1 = fig.add_subplot(2,3,1)
#ax1.scatter(date,vect1,marker='o')
vect1 = invert_quad(date,vect1,ax1)

ax2 = fig.add_subplot(2,3,2)
#ax2.scatter(date,vect2,marker='o')
vect2 = invert_quad(date,vect2,ax2)

ax3 = fig.add_subplot(2,3,3)
#ax3.scatter(date,vect3,marker='o')
vect3 = invert_quad(date,vect3,ax3)

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
    ax.scatter(x,y,marker='o')
    ax.plot(xx,seas(xx,pars[0],pars[1],pars[2]),label='y = {:0.2f}cos(wt)+{:0.2f}sin(wt)'.format(pars[1],pars[2]))
    ax.fill_between(xx,seas(xx,parsmin[0],parsmin[1],parsmin[2]),seas(xx,parsmax[0],parsmax[1],parsmax[2]),facecolor='#CC4F1B',alpha=0.5,linewidth=0)
    ax.set_ylim([ymin,ymax])
    ax.set_xlim(0,1)
    plt.legend(loc='best')

ymin,ymax = np.min([vect1,vect2,vect3]), np.max([vect1,vect2,vect3]) 

ax1 = fig.add_subplot(2,3,4)
invert_seas(dd,vect1,ax1)

ax2 = fig.add_subplot(2,3,5)
invert_seas(dd,vect2,ax2)

ax3 = fig.add_subplot(2,3,6)
invert_seas(dd,vect3,ax3)

fig.savefig('acp_temp.eps', format='EPS',dpi=150)
plt.show()
