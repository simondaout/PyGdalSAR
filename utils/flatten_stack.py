#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import *

import numpy as np
from numpy.lib.stride_tricks import as_strided

import scipy
import scipy.optimize as opt
import scipy.linalg as lst



# T319
# nlign, ncol = 4078, 1420 
# ibeg, iend = 1400,2600

#T455
nlign, ncol = 4082, 1420 
ibeg, iend = 1300,2750

#T047
# nlign, ncol = 4077, 1420 
# ibeg, iend = 1700,2900

# maskfile = 'phs_sum_eq1'
# outfile = 'phs_sum_eq1_flat'  
# maskfile = 'phs_sum_eq2'
# outfile = 'phs_sum_eq2_flat'  
maskfile = 'phs_sum_eq12'
outfile = 'phs_sum_eq12_flat'  

# nlign, ncol = 4078, 1420 
# ibeg, iend = 1500,2500
# maskfile = 'test_eq12'
# outfile = 'test_eq12_flat'
 

perc=80.

# open
mask = np.zeros((nlign,ncol))
fid = open(maskfile,'r')
mask[ibeg:iend,:] = np.fromfile(fid,dtype=np.float32).reshape((nlign,ncol))[ibeg:iend,:]   
fid.close()

# flatten
los_temp = np.copy(mask)
los_temp[mask==0] = np.float('NaN')
maxlos,minlos=np.nanpercentile(los_temp,perc),np.nanpercentile(los_temp,(100-perc))
print maxlos,minlos

index = np.nonzero(
	np.logical_and(mask!=0,
	np.logical_and(~np.isnan(mask),
	np.logical_and(mask>minlos,mask<maxlos
)
)
)
)

spacial_mask = np.zeros((nlign,ncol))
spacial_mask[index] = mask[index]
temp = np.array(index).T
az = temp[:,0]; rg = temp[:,1]
los = np.copy(mask[index]).flatten()

G=np.zeros((len(los),5))
G[:,0] = rg**2
G[:,1] = rg
G[:,2] = az**2
G[:,3] = az
G[:,4] = 1
x0 = lst.lstsq(G,los)[0]
_func = lambda x: np.sum(((np.dot(G,x)-los))**2)
_fprime = lambda x: 2*np.dot(G.T, (np.dot(G,x)-los))
pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=200,full_output=True,iprint=0)[0]
a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
print 'Remove ramp %f rg**2 %f rg  + %f az**2 + %f az + %f'%(pars[0],pars[1],pars[2],pars[3],pars[4])

# build total G matrix
G=np.zeros((len(mask.flatten()),5))
for i in xrange(nlign):
    G[i*ncol:(i+1)*ncol,0] = np.arange(ncol)**2
    G[i*ncol:(i+1)*ncol,1] =  np.arange(ncol)
    G[i*ncol:(i+1)*ncol,2] = i**2  
    G[i*ncol:(i+1)*ncol,3] = i    
G[:,4] = 1

ramp = np.dot(G,pars).reshape(nlign,ncol)
mask_flat = mask - ramp
mask_flat[mask==0] = 0

# clean border 
mask_flat[:,:17] = 0

# moy = np.nanpercentile(mask_flat[ibeg:ibeg+200,:],95)
# mask_flat[ibeg-200:ibeg,:] = np.repeat(np.linspace(0,moy,200),ncol).reshape(200,ncol) 
# mask_flat[ibeg:ibeg+200,:] = mask_flat[ibeg:ibeg+200,:] - np.repeat(moy - np.linspace(0,moy,200),ncol).reshape(200,ncol)

# moy = np.nanpercentile(mask_flat[iend-200:iend,:],95)
# mask_flat[iend:iend+200,:] = np.repeat(moy-np.linspace(0,moy,200),ncol).reshape(200,ncol) 
# mask_flat[iend-200:iend,:] = mask_flat[iend-200:iend,:] - np.repeat(np.linspace(0,moy,200),ncol).reshape(200,ncol) 

mask_flat.astype('float32').tofile(outfile)

fig = plt.figure(figsize=(9,5))
ax = fig.add_subplot(1,3,1)
cax = plt.imshow(mask,cmap='plasma',vmin=-2, vmax=2)
ax.set_title('stack')
setp( ax.get_xticklabels(), visible=None)
fig.colorbar(cax, orientation='vertical',aspect=10)


ax = fig.add_subplot(1,3,2)
cax = plt.imshow(spacial_mask,cmap='plasma',vmin=-2, vmax=2)
ax.set_title('maskl')
setp( ax.get_xticklabels(), visible=None)
fig.colorbar(cax, orientation='vertical',aspect=10)


ax = fig.add_subplot(1,3,3)
cax = plt.imshow(mask_flat,cmap='plasma',vmin=-2, vmax=2)
ax.set_title('flatten stack')
setp( ax.get_xticklabels(), visible=None)
fig.colorbar(cax, orientation='vertical',aspect=10)
plt.show()