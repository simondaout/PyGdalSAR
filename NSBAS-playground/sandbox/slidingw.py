#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from numpy.lib.stride_tricks import as_strided
from functools import wraps
from textwrap import dedent
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import *



#---- functions ----
def _check(a, r_c):
    """Performs the array checks necessary for stride and block.
    : a   - Array or list.
    : r_c - tuple/list/array of rows x cols.  
    :Attempts will be made to 
    :     produce a shape at least (1*c).  For a scalar, the
    :     minimum shape will be (1*r) for 1D array or (1*c) for 2D
    :     array if r<c.  Be aware
    """
    if isinstance(r_c, (int, float)):
        r_c = (1, int(r_c))
    r, c = r_c
    a = np.atleast_2d(a)
    shp = a.shape
    r, c = r_c = ( min(r, a.shape[0]), min(c, shp[1]) ) 
    a = np.ascontiguousarray(a)
    return a, shp, r, c, tuple(r_c)
    
def stride(a, r_c=(3, 3)):
    """Provide a 2D sliding/moving view of an array.  
    :  There is no edge correction for outputs.
    :
    :Requires
    :--------
    : a - array or list, usually a 2D array.  Assumes rows is >=1,
    :     it is corrected as is the number of columns.
    : r_c - tuple/list/array of rows x cols.  Attempts  to 
    :     produce a shape at least (1*c).  For a scalar, the
    :     minimum shape will be (1*r) for 1D array or 2D
    :     array if r<c.  Be aware
    """
    a, shp, r, c, r_c = _check(a, r_c)
    shape = (a.shape[0] - r + 1, a.shape[1] - c + 1) + r_c
    strides = a.strides * 2
    a_s = (as_strided(a, shape=shape, strides=strides)).squeeze()    
    return a_s


def block(a, r_c=(3, 3)):
    """See _check and/or stride for documentation.  This function
    :  moves in increments of the block size, rather than sliding
    :  by one row and column
    :
    """
    a, shp, r, c, r_c = _check(a, r_c)
    shape = (a.shape[0]/r, a.shape[1]/c) + r_c
    strides = (r*a.strides[0], c*a.strides[1]) + a.strides
    a_b = as_strided(a, shape=shape, strides=strides).squeeze()
    return a_b

n = 500
a = np.arange(n*n).reshape((n, n))
a0 = stride(a, r_c=(3, 3))  # uncomment this or below
# a0 = block(a, block=(3, 3))
print("\nArray size... {0}x{0}".format(n))
print(a0)
print(np.shape(a0))

fig = plt.figure(0,figsize=(12,8))
ax = fig.add_subplot(1,2,1)
# hax = ax.imshow(mask, cm.Greys, vmin=0, vmax=seuil)
cax = ax.imshow(a, cm.RdBu)
setp( ax.get_xticklabels(), visible=False)
cbar = fig.colorbar(cax, orientation='vertical',aspect=9)

ax = fig.add_subplot(1,2,2)
cax = ax.imshow(a0, cm.RdBu)
setp( ax.get_xticklabels(), visible=False)

plt.show()

