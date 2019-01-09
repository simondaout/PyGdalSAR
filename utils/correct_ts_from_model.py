#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT (Oxford)
################################################################################


import numpy as np
from numpy.lib.stride_tricks import as_strided
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
from pylab import *
import math
import sys
import matplotlib.dates as mdates
import datetime

def Heaviside(t):
        h=np.zeros((len(t)))
        h[t>=0]=1.0
        return h

def Box(t):
        return Heaviside(t+0.5)-Heaviside(t-0.5)

class postseismic():
    def __init__(self,date,tcar=1.):
        self.to=date
        self.tcar=tcar

    def g(self,t):
        t=(t-self.to)/self.tcar
        t[t<=0] = 0
        g = np.log10(1+t)
        return g

class coseismic():
      def __init__(self,date):
          self.to=date

      def g(self,t):
        return Heaviside(t-self.to)

# load images_retenues file
nb,idates,dates,base=np.loadtxt('images_retenues', comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
N=len(dates)
datemin, datemax = np.int(np.min(dates)), np.int(np.max(dates))+1
print 'Number images: ', N
# read lect.in
ncol, nlign = map(int, open('lect_ts.in').readline().split(None, 2)[0:2])
# lect cube
cube = np.fromfile('depl_cumule_flat',dtype=np.float32)
maps = cube[:nlign*ncol*N].reshape((nlign,ncol,N))

# compute model
co0 = 2003.29
tcar=.3
co1 = 2004.36
infile = 'post0_coeff_clean.r4'
post0 = np.fromfile(infile,dtype=np.float32).reshape((nlign,ncol))
infile = 'cos0_coeff_clean.r4'
cos0 = np.fromfile(infile,dtype=np.float32).reshape((nlign,ncol))
#infile = 'cos1_coeff_clean.r4'
#cos1 = np.fromfile(infile,dtype=np.float32).reshape((nlign,ncol))
infile = 'ref_coeff.r4'
ref = np.fromfile(infile,dtype=np.float32).reshape((nlign,ncol))

# crop
# ibeg,iend = 2180,2640
# jbeg,jend = 340,800

ibeg,iend = 0,nlign
jbeg,jend = 0,ncol

cos1=np.zeros((nlign,ncol))

model = np.zeros((nlign,ncol,N))
for i in xrange(ibeg,iend):
    for j in xrange(jbeg,jend):
        # model[i,j,:] = ref[i,j]
        model[i,j,:] = cos0[i,j]*coseismic(co0).g(dates) +post0[i,j]*postseismic(co0,tcar).g(dates) + \
        cos1[i,j]*coseismic(co1).g(dates) + ref[i,j]

# remove model
# maps_clean = np.zeros((nlign,ncol,N))
# for l in xrange((N)):
#     cst = np.nanmean(maps[:,:,l] - model[:,:,l])
#     maps_clean[:,:,l] = maps[:,:,l] - model[:,:,l] 
maps_clean = maps - model

save
fid = open('depl_cumule_nomodel', 'wb')
maps_clean.flatten().astype('float32').tofile(fid)
fid.close()

# plot TS
rad2mm = -4.4563
vmax,vmin=30,-30
maps_clean = maps_clean*rad2mm
maps = maps*rad2mm
model = model*rad2mm

ipix = 2380,2411,2250,2200
jpix = 571,640,510,400
Npix = len(ipix)
fig = plt.figure(0,figsize=(10,9))
fig.subplots_adjust(wspace=0.001)
for jj in xrange((Npix)):
    i, j = ipix[jj], jpix[jj]

    disp_clean = maps_clean[i,j,:]
    disp = maps[i,j,:]
    mdisp = model[i,j,:]

    x = [date2num(datetime.datetime.strptime('{}'.format(d),'%Y%m%d')) for d in idates]
   
    dmax = str(datemax) + '0101'
    dmin = str(datemin) + '0101'
    xmin = datetime.datetime.strptime('{}'.format(dmin),'%Y%m%d')
    xmax = datetime.datetime.strptime('{}'.format(dmax),'%Y%m%d')
    xlim=date2num(np.array([xmin,xmax]))
    xlim=date2num(np.array([xmin,xmax]))
    ax = fig.add_subplot(Npix,1,jj+1)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))

    ax.plot(x,disp,'o',label='lign:{}, column:{}'.format(i,j))
    ax.plot(x,disp_clean,'x')
    
    t = np.array([xmin + datetime.timedelta(days=d) for d in range(0, 2920)])
    tdec = np.array([float(date.strftime('%Y')) + float(date.strftime('%j'))/365.1 for date in t])
    m = ( cos0[i,j]*coseismic(co0).g(tdec) +post0[i,j]*postseismic(co0,tcar).g(tdec) + \
        cos1[i,j]*coseismic(co1).g(tdec) + ref[i,j] )*rad2mm

    # m = (np.ones(len(t))*ref[i,j])*rad2mm

    ax.plot(t,m,'-r')

    ax.set_xlim(xlim)
    fig.autofmt_xdate()
    ax.set_xlabel('Time (Year/month/day)')
    ax.set_ylabel('Displacements (mm)')

# plt.show()
# sys.exit()

# plot diplacements maps
fig = plt.figure(1,figsize=(10,9))
fig.subplots_adjust(wspace=0.001)

for l in xrange((N)):
    d = as_strided(maps_clean[ibeg:iend,jbeg:jend,l])
    ax = fig.add_subplot(4,int(N/4)+1,l+1)
    cmap = cm.jet
    cmap.set_bad('white')
    cax = ax.imshow(d,cmap=cm.rainbow,vmax=vmax,vmin=vmin)
    ax.set_title(idates[l],fontsize=6)
    setp( ax.get_xticklabels(), visible=False)
    setp( ax.get_yticklabels(), visible=False)

setp(ax.get_xticklabels(), visible=False)
setp(ax.get_yticklabels(), visible=False)
fig.tight_layout()
plt.suptitle('Corrected TS')
fig.colorbar(cax, orientation='vertical',aspect=10)
fig.savefig('maps_clean_models.eps', format='EPS',dpi=150)

# plot diplacements maps
fig = plt.figure(2,figsize=(10,9))
fig.subplots_adjust(wspace=0.001)

for l in xrange((N)):
    d = as_strided(maps[ibeg:iend,jbeg:jend,l])
    ax = fig.add_subplot(4,int(N/4)+1,l+1)
    cmap = cm.jet
    cmap.set_bad('white')
    cax = ax.imshow(d,cmap=cm.rainbow,vmax=vmax,vmin=vmin)
    ax.set_title(idates[l],fontsize=6)
    setp( ax.get_xticklabels(), visible=False)
    setp( ax.get_yticklabels(), visible=False)

setp(ax.get_xticklabels(), visible=False)
setp(ax.get_yticklabels(), visible=False)
fig.tight_layout()
plt.suptitle('TS')
fig.colorbar(cax, orientation='vertical',aspect=10)


plt.show()



