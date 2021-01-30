#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.stride_tricks import as_strided

nlign, ncol = 1821, 1828 
N = 56
im = 34

cubei = np.fromfile("depl_cumule",dtype=np.float32)
cube = as_strided(cubei[:nlign*ncol*N])
kk = np.flatnonzero(np.logical_or(cube==9990, cube==9999))
cube[kk] = float('NaN')
maps = cube.reshape((nlign,ncol,N))

eq  = np.copy(maps[:,:,im])

mask = -(np.abs(eq) - 0.8)
kk = np.nonzero(np.logical_or(mask==0,mask<0))
#mask[kk] = float('NaN')

fid = open('mask_eqs.r4', 'wb')
mask[:nlign,:ncol].flatten().astype('float32').tofile(fid)

fig, _ = plt.subplots(1)
vmax = 2
cax = fig.axes[0].imshow(mask,vmax=vmax, vmin=-vmax, cmap='seismic')
fig.axes[0].set_title('Mask')
fig.tight_layout()
fig.colorbar(cax, orientation='horizontal',aspect=10,shrink=10)
fig.savefig('create_mask.eps',format = 'EPS')
plt.show()

