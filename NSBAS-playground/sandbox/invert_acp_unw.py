#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from numpy.lib.stride_tricks import as_strided
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
from pylab import setp
import scipy.optimize as sp
import scipy.linalg as lst
from datetime import datetime
import datetime
import sys

####################################################################################

def quad(x,a0,a1):
    return a0+a1*x

def seas(x,a0,a1,a2):
    return a0+a1*np.cos(2*np.pi*x)+a2*np.sin(2*np.pi*x)

def invert_quad(x,y,ax):
    pars, cova = sp.curve_fit(quad, x, y)
    perr = 2*np.sqrt(np.diag(cova))
    reg = quad(xx,pars[0],pars[1]) 
    ax.scatter(x,y,marker='o',alpha=0.6)
    ax.plot(xx,quad(xx,pars[0],pars[1]),'-r',lw=3,label='y = {:2.4f}x'.format(pars[1]))
    # remove quad term to the vect
    ax.set_ylim([ymin,ymax])
    out = y - quad(x,pars[0],pars[1])
    plt.legend(loc='best')
    return out

def invert_seas(x,y,ax):
    pars, cova = sp.curve_fit(seas, x, y)
    perr = 2*np.sqrt(np.diag(cova))
    parsmax,parsmin = pars+perr, pars-perr
    reg, regmin,regmax = seas(xx,pars[0],pars[1],pars[2]), seas(xx,parsmin[0],parsmin[1],parsmin[2]),seas(xx,parsmax[0],parsmax[1],parsmax[2]) 
    ax.scatter(x,y,marker='o',alpha=0.6)
    ax.plot(xx,seas(xx,pars[0],pars[1],pars[2]),'-r',lw=3, label='y = {:0.2f}cos(wt)+{:0.2f}sin(wt)'.format(pars[1],pars[2]))
    # ax.fill_between(xx,seas(xx,parsmin[0],parsmin[1],parsmin[2]),seas(xx,parsmax[0],parsmax[1],parsmax[2]),facecolor='#CC4F1B',alpha=0.5,linewidth=0)
    ax.set_ylim([ymin,ymax])
    ax.set_xlim(0,1)
    plt.legend(loc='best')

def date2dec(dates):
    dates  = np.atleast_1d(dates)
    tdatees = []
    for date in dates:
        x = datetime.datetime.strptime('{}'.format(date),'%Y%m%d')
        #x = datetime.strptime('{}'.format(date),'%Y-%m-%dT%H:%M:%S.%fZ')
        dec = float(x.strftime('%j'))/365.1
        year = float(x.strftime('%Y'))
        tdatees.append(year + dec)
    return tdatees

####################################################################################

cmap=cm.rainbow
fig = plt.figure(3,figsize=(12,8))
ncol, nlign = map(int, open('lect.in').readline().split(None, 2)[0:2])
m = np.fromfile('acp_1',dtype=np.float32)[:nlign*ncol].reshape((nlign,ncol))
vmax = np.nanpercentile(m,95)
vmin = np.nanpercentile(m,5)
ax = fig.add_subplot(1,3,1)
cax = ax.imshow(m,cmap, vmin=vmin,vmax=vmax)
setp( ax.get_xticklabels(), visible=False)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
fig.canvas.set_window_title('PCA1')
plt.colorbar(cax, cax=c)

m = np.fromfile('acp_2',dtype=np.float32)[:nlign*ncol].reshape((nlign,ncol))
vmax = np.nanpercentile(m,95)
vmin = np.nanpercentile(m,5)
ax = fig.add_subplot(1,3,2)
cax = ax.imshow(m,cmap,vmin=vmin,vmax=vmax)
setp( ax.get_xticklabels(), visible=False)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
fig.canvas.set_window_title('PCA2')
plt.colorbar(cax, cax=c)

m = np.fromfile('acp_3',dtype=np.float32)[:nlign*ncol].reshape((nlign,ncol))
vmax = np.nanpercentile(m,95)
vmin = np.nanpercentile(m,5)
ax = fig.add_subplot(1,3,3)
cax = ax.imshow(m,cmap,vmin=vmin,vmax=vmax)
setp( ax.get_xticklabels(), visible=False)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
fig.canvas.set_window_title('PCA3')
plt.colorbar(cax, cax=c)

# open acp vecteurs 
source=open('acp_vecteurs','r')
date_1,date_2,vect1,vect2,vect3 = np.loadtxt(source,comments="#",usecols=(0,1,2,3,4),unpack=True,dtype='i,i,f,f,f')
Nifg=len(date_1)
vect = np.vstack([vect1,vect2,vect3]).T

# list dates
date = []; bt = []
for date1,date2 in zip(date_1,date_2):
    if date1 not in date: date.append(date1)
    if date2 not in date: date.append(date2)
nmax=len(date)
print("number of images: {} ".format(nmax))
dated = date2dec(date)
cst = np.copy(dated[0])
# compute temporal baseline for TS inv
for i in range((nmax)):
    bt.append(dated[i]-cst)

# invert vect into TS
G_=np.zeros((Nifg,nmax))
deltat = np.zeros((Nifg))
for k in range((Nifg)):
  for n in range((nmax)):
    if (date_1[k]==date[n]): 
      G_[k,n]=-1
      t1 = bt[n]
    elif (date_2[k]==date[n]):
      G_[k,n]=1
      t2 = bt[n]
  deltat[k] = abs(t2 -t1)

w1 = np.exp(-(deltat/2.)); sig_ = 1./w1
d = np.zeros(((Nifg+1)))
sig = np.ones(((Nifg+1)))
G = np.zeros(((Nifg+1),nmax))

vect_date = np.zeros((nmax,3))
for j in range(3):
    d[:Nifg] = as_strided(vect[:,j])
    G[:Nifg,:nmax] = G_ 
    G[-1,0] = 1 # ini phi first dateage to 0
    sig[:Nifg] = sig_

    x0 = lst.lstsq(G,d)[0]
    _func = lambda x: np.sum(((np.dot(G,x)-d)/sig)**2)
    _fprime = lambda x: 2*np.dot(G.T/sig, (np.dot(G,x)-d)/sig)
    vect_date[:,j] = sp.fmin_slsqp(_func,x0,fprime=_fprime,iter=500,full_output=True,iprint=0)[0]

vect1,vect2,vect3 = vect_date[:,0],vect_date[:,1],vect_date[:,2]

####################################################################################
# PLOT

fig = plt.figure(10, figsize=(15,8))
fig.subplots_adjust(hspace=0.2)

# 1) temp behavior
dated = np.atleast_1d(dated)
xx = np.arange(dated.min(),dated.max(),0.01)
ymin,ymax = np.min([vect1,vect2,vect3]), np.max([vect1,vect2,vect3]) 

ax1 = fig.add_subplot(2,3,1)
#ax1.scatter(dated,vect1,marker='o')
vect1 = invert_quad(dated,vect1,ax1)

ax2 = fig.add_subplot(2,3,2)
#ax2.scatter(dated,vect2,marker='o')
vect2 = invert_quad(dated,vect2,ax2)

ax3 = fig.add_subplot(2,3,3)
#ax3.scatter(dated,vect3,marker='o')
vect3 = invert_quad(dated,vect3,ax3)

# 2) inversion
dd = np.mod(dated,1)
xx = np.arange(dd.min(),dd.max(),0.01)
ymin,ymax = np.min([vect1,vect2,vect3]), np.max([vect1,vect2,vect3]) 

ax1 = fig.add_subplot(2,3,4)
invert_seas(dd,vect1,ax1)

ax2 = fig.add_subplot(2,3,5)
invert_seas(dd,vect2,ax2)

ax3 = fig.add_subplot(2,3,6)
invert_seas(dd,vect3,ax3)

fig.savefig('acp_temp.eps', format='EPS',dpi=150)



plt.show()
