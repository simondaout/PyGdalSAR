####!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################

############################################
# Author        : Simon DAOUT (CRPG/ENSG)
############################################

"""invert_ramp_topo_unw.py — Empirical InSAR ramp/topo correction. Use -h for help."""

print()
print()
print('Author: Simon Daout')
print('Please cite:')
print(
    'Daout, S., Sudhaus, H., Kausch, T., Steinberg, A., & Dini, B. (2019). Interseismic and postseismic shallow creep of the North Qaidam Thrust faults detected with a multitemporal InSAR analysis. Journal of Geophysical Research: Solid Earth, 124(7), 7259-7279.')
print()
print()

# gdal — deferred: importing gdal at module level inits CoreFoundation before spawn
def _gdal():
    """Lazy gdal import — call inside worker functions only."""
    from osgeo import gdal as _g
    _g.UseExceptions()
    return _g

# system
from os import path, environ, getcwd
import os, sys

# plot — only lightweight base at module level; pyplot/mpl_toolkits deferred
import matplotlib
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap

from datetime import datetime
import logging

# numpy
import numpy as np

# scipy
import scipy.ndimage
import scipy.optimize as opt
import scipy.linalg as lst

import argparse
import shutil

from contextlib import contextmanager
from functools import wraps, partial
import multiprocessing

multiprocessing.set_start_method('spawn', force=True)

import warnings

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


##################################################################################
###  Extras functions and context maganers
##################################################################################

def makedirs(name):
    os.makedirs(name, exist_ok=True)


def _outpath(filename):
    """Return the correct output path depending on the format (set after arg parsing)."""
    if sformat == 'ROI_PAC':
        return filename
    return out_path + filename


def date2dec(dates):
    dates = np.atleast_1d(dates)
    times = []
    for date in dates:
        x = datetime.strptime(str(date), '%Y%m%d')
        dec = float(x.strftime('%j')) / 365.1
        year = float(x.strftime('%Y'))
        times.append(year + dec)
    return times


# Timer for all the functions
class ContextDecorator(object):
    def __call__(self, f):
        @wraps(f)
        def decorated(*args, **kwds):
            with self:
                try:
                    return f(*args, **kwds)
                except (KeyboardInterrupt, SystemExit):
                    raise
                except:
                    Exception('{0} Failed !'.format(f))
                    raise

        return decorated


class TimeIt(ContextDecorator):
    def __enter__(self):
        self.start = datetime.now()
        print('Starting time process: {0}'.format(self.start))

    def __exit__(self, type, value, traceback):
        print('Time process: {0}s'.format((datetime.now() - self.start).total_seconds()))


def checkinfile(file):
    if not path.exists(file):
        print("File: {0} not found in {1}, Exit!".format(file, getcwd()))
        sys.exit(1)


# create generator for pool
def _init_worker(shared):
    """Inject parent globals into each spawn worker.

    With spawn, workers start with a clean interpreter — no parent globals.
    We store everything in the module-level dict _G, then each worker function
    calls globals().update(_G) at its start to make variables visible locally.
    """
    _G.update(shared)


# Module-level store for globals shared with pool workers
# (populated by _init_worker at pool startup)
_G = {}


@contextmanager
def poolcontext(*arg, **kargs):
    pool = multiprocessing.Pool(*arg, **kargs)
    yield pool
    pool.terminate()
    pool.join()


#####################################################################################
# FUNCTIONS
#####################################################################################

def _solve_ramp(G, data, rms):
    """
    Solve the weighted least-squares ramp inversion.
    Uses lstsq as initial guess, then refines with fmin_slsqp.
    Returns parameter vector `pars`.
    """
    x0 = lst.lstsq(G, data)[0]
    try:
        _func = lambda x: np.sum(((np.dot(G, x) - data) / rms) ** 2)
        _fprime = lambda x: 2 * np.dot(G.T / rms, (np.dot(G, x) - data) / rms)
        pars = opt.fmin_slsqp(_func, x0, fprime=_fprime,
                              iter=2000, full_output=True, iprint=0, acc=1e-9)[0]
    except Exception:
        pars = x0
    return pars


def estim_ramp(los, data, topo_clean, az, rg, order, rms, nfit, ivar, los_ref, rg_ref, az_ref, topo_ref):
    """
    Empircal estmation function on flatten los vector
    """

    # initialise full vector
    sol = np.zeros((13))
    # 0:x**3 1:x**2 2:x 3:y**3 4:y**2 5:y 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2
    # initialize correction
    corr = np.zeros((mlines, mcols))

    # Offset pixel coordinates: same convention as original script.
    az_off = az - ibeg
    rg_off = rg - jbeg

    # Pre-compute vectorized coordinate arrays for the full image
    _col_idx = np.tile(np.arange(mcols) - jbeg, mlines)    # (rg - jbeg) flat
    _row_idx = np.repeat(np.arange(mlines) - ibeg, mcols)  # (az - ibeg) flat
    _elev_flat = elev_map.flatten()
    _col2 = _col_idx ** 2
    _row2 = _row_idx ** 2

    if order == 0:
        # 0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2

        if radar is None:
            pars = los_ref
            sol[8] = los_ref
            print('Remove ref frame within the ref area %f' % (los_ref))

            # build total G matrix
            G = np.zeros((len(los), 1))
            G[:, 0] = 1

        else:
            if ivar == 0 and nfit == 0:
                G = np.zeros((len(data), 2))
                G[:, 0] = 1
                G[:, 1] = topo_clean

                # ramp inversion
                pars = _solve_ramp(G, data, rms)

                sol[8] = pars[0];
                sol[9] = pars[1]
                print('Remove ref frame %f + %f z' % (pars[0], pars[1]))

                # build total G matrix
                G = np.column_stack([np.ones(len(los)), _elev_flat])  # vectorized

            elif ivar == 0 and nfit == 1:
                G = np.zeros((len(data), 3))
                G[:, 0] = 1
                G[:, 1] = topo_clean
                G[:, 2] = topo_clean ** 2

                # ramp inversion
                pars = _solve_ramp(G, data, rms)

                sol[8] = pars[0];
                sol[9] = pars[1];
                sol[10] = pars[2]
                print('Remove ref frame %f + %f z + %f z**2' % (pars[0], pars[1], pars[2]))

                # build total G matrix
                G = np.column_stack([np.ones(len(los)), _elev_flat, _elev_flat ** 2])  # vectorized

            elif ivar == 1 and nfit == 0:
                G = np.zeros((len(data), 3))
                G[:, 0] = 1
                G[:, 1] = topo_clean
                G[:, 2] = az * topo_clean

                # ramp inversion
                pars = _solve_ramp(G, data, rms)

                sol[8] = pars[0];
                sol[9] = pars[1];
                sol[11] = pars[2]
                print('Remove ref frame %f + %f z + %f az*z' % (pars[0], pars[1], pars[2]))

                # build total G matrix (vectorized)
                G = np.column_stack([np.ones(len(los)), _elev_flat, _elev_flat * _row_idx])

            elif ivar == 1 and nfit == 1:
                G = np.zeros((len(data), 5))
                G[:, 0] = 1
                G[:, 1] = topo_clean
                G[:, 2] = topo_clean ** 2
                G[:, 3] = az * topo_clean
                G[:, 4] = (az * topo_clean) ** 2

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[8] = pars[0];
                sol[9] = pars[1];
                sol[10] = pars[2];
                sol[11] = pars[3];
                sol[12] = pars[4]
                print(
                    'Remove ref frame %f + %f z + %f z**2 + %f az*z + %f (az*z)**2' % (pars[0], pars[1], pars[2],
                                                                                       pars[3], pars[4]))

                # build total G matrix (vectorized) — NOTE: was G=zeros(len,4) but needed 5 cols (bug fix)
                _ez = _elev_flat * _row_idx
                G = np.column_stack([np.ones(len(los)), _elev_flat, _elev_flat ** 2, _ez, _ez ** 2])

    elif order == 1:  # Remove a range ramp ay+b for each maps (y = col)
    # 0:x**3 1:x**2 2:x 3:y**3 4:y**2 5:y 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2

        if radar is None:
            G = np.zeros((len(data), 2))
            G[:, 0] = rg_off
            G[:, 1] = 1

            # ramp inversion
            pars = _solve_ramp(G, data, rms)
            sol[2] = pars[0];
            sol[8] = pars[1]
            print('Remove ramp %f r + %f' % (pars[0], pars[1]))

            # build total G matrix (vectorized)
            G = np.column_stack([_col_idx, np.ones(len(los))])


        else:
            if ivar == 0 and nfit == 0:
                G = np.zeros((len(data), 3))
                G[:, 0] = rg_off
                G[:, 1] = 1
                G[:, 2] = topo_clean

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[2] = pars[0];
                sol[8] = pars[1];
                sol[9] = pars[2]
                print('Remove ramp %f r + %f + %f z ' % (pars[0], pars[1], pars[2]))

                # build total G matrix (vectorized)
                G = np.column_stack([_col_idx, np.ones(len(los)), _elev_flat])

            elif ivar == 0 and nfit == 1:
                G = np.zeros((len(data), 4))
                G[:, 0] = rg_off
                G[:, 1] = 1
                G[:, 2] = topo_clean
                G[:, 3] = topo_clean ** 2

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[2] = pars[0];
                sol[8] = pars[1];
                sol[9] = pars[2];
                sol[10] = pars[3]
                print('Remove ramp %f r + %f + %f z + %f z**2' % (pars[0], pars[1], pars[2], pars[3]))

                # build total G matrix (vectorized)
                G = np.column_stack([_col_idx, np.ones(len(los)), _elev_flat, _elev_flat ** 2])

            elif ivar == 1 and nfit == 0:
                G = np.zeros((len(data), 4))
                G[:, 0] = rg  # was `y` (undefined variable) — bug fix
                G[:, 1] = 1
                G[:, 2] = topo_clean
                G[:, 3] = topo_clean * az_off

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[2] = pars[0];
                sol[8] = pars[1];
                sol[9] = pars[2];
                sol[11] = pars[3]
                print('Remove ramp %f r + %f + %f z + %f z*az' % (pars[0], pars[1], pars[2], pars[3]))

                # build total G matrix (vectorized)
                G = np.column_stack([_col_idx, np.ones(len(los)), _elev_flat, _elev_flat * _row_idx])

            elif ivar == 1 and nfit == 1:
                G = np.zeros((len(data), 6))
                G[:, 0] = rg_off
                G[:, 1] = 1
                G[:, 2] = topo_clean
                G[:, 3] = topo_clean ** 2
                G[:, 4] = topo_clean * az_off
                G[:, 5] = (topo_clean * az_off) ** 2

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[2] = pars[0];
                sol[8] = pars[1];
                sol[9] = pars[2];
                sol[10] = pars[3];
                sol[11] = pars[4];
                sol[12] = pars[5]
                print(
                    'Remove ramp %f r + %f + %f z + %f z**2 + %f z*az' % (pars[0], pars[1], pars[2], pars[3], pars[4],
                                                                          pars[5]))

                # build total G matrix (vectorized)
                _ez = _elev_flat * _row_idx
                G = np.column_stack([_col_idx, np.ones(len(los)), _elev_flat, _elev_flat ** 2, _ez, _ez ** 2])

    elif order == 2:  # Remove an azimutal ramp ax+b for each maps (x is row)
    # 0:x**3 1:x**2 2:x 3:y**3 4:y**2 5:y 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2

        if radar is None:
            G = np.zeros((len(data), 2))
            G[:, 0] = az_off
            G[:, 1] = 1

            # ramp inversion
            pars = _solve_ramp(G, data, rms)
            sol[5] = pars[0];
            sol[8] = pars[1]
            print('Remove ramp %f az + %f' % (pars[0], pars[1]))

            # build total G matrix (vectorized)
            G = np.column_stack([_row_idx, np.ones(len(los))])

        else:
            if ivar == 0 and nfit == 0:
                G = np.zeros((len(data), 3))
                G[:, 0] = az_off
                G[:, 1] = 1
                G[:, 2] = topo_clean

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[5] = pars[0];
                sol[8] = pars[1];
                sol[9] = pars[2]
                print('Remove ramp %f az + %f + %f z' % (pars[0], pars[1], pars[2]))

                # build total G matrix (vectorized)
                G = np.column_stack([_row_idx, np.ones(len(los)), _elev_flat])

            elif ivar == 0 and nfit == 1:
                G = np.zeros((len(data), 4))
                G[:, 0] = az_off
                G[:, 1] = 1
                G[:, 2] = topo_clean
                G[:, 3] = topo_clean ** 2

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[5] = pars[0];
                sol[8] = pars[1];
                sol[9] = pars[2];
                sol[10] = pars[3]
                print('Remove ramp %f az + %f + %f z + %f z**2' % (pars[0], pars[1], pars[2], pars[3]))

                # build total G matrix (vectorized)
                G = np.column_stack([_row_idx, np.ones(len(los)), _elev_flat, _elev_flat ** 2])

            elif ivar == 1 and nfit == 0:
                G = np.zeros((len(data), 4))
                G[:, 0] = az_off
                G[:, 1] = 1
                G[:, 2] = topo_clean
                G[:, 3] = topo_clean * az_off

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[5] = pars[0];
                sol[8] = pars[1];
                sol[9] = pars[2];
                sol[11] = pars[3]
                print('Remove ramp %f az + %f + %f z + %f z*az' % (pars[0], pars[1], pars[2], pars[3]))

                # build total G matrix (vectorized)
                G = np.column_stack([_row_idx, np.ones(len(los)), _elev_flat, _elev_flat * _row_idx])

            elif ivar == 1 and nfit == 1:
                G = np.zeros((len(data), 5))
                G[:, 0] = az_off
                G[:, 1] = 1
                G[:, 2] = topo_clean
                G[:, 3] = topo_clean * az_off
                G[:, 4] = (topo_clean * az_off) ** 2

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[5] = pars[0];
                sol[8] = pars[1];
                sol[9] = pars[2];
                sol[11] = pars[3];
                sol[12] = pars[4]
                print(
                    'Remove ramp %f az + %f + %f z + %f z*az + %f (z*az)**2' % (pars[0], pars[1], pars[2], pars[3],
                                                                                pars[4]))

                # build total G matrix (vectorized)
                _ez = _elev_flat * _row_idx
                G = np.column_stack([_row_idx, np.ones(len(los)), _elev_flat, _ez, _ez ** 2])

    elif order == 3:  # Remove a ramp ay+bx+c for each maps
    # 0:x**3 1:x**2 2:x 3:y**3 4:y**2 5:y 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2

        if radar is None:
            G = np.zeros((len(data), 3))
            G[:, 0] = rg_off
            G[:, 1] = az_off
            G[:, 2] = 1

            # ramp inversion
            pars = _solve_ramp(G, data, rms)
            sol[2] = pars[0];
            sol[5] = pars[1];
            sol[8] = pars[2]
            print('Remove ramp %f r  + %f az + %f' % (pars[0], pars[1], pars[2]))

            # build total G matrix (vectorized)
            G = np.column_stack([_col_idx, _row_idx, np.ones(len(los))])

        else:
            if ivar == 0 and nfit == 0:
                G = np.zeros((len(data), 4))
                G[:, 0] = rg_off
                G[:, 1] = az_off
                G[:, 2] = 1
                G[:, 3] = topo_clean

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[2] = pars[0];
                sol[5] = pars[1];
                sol[8] = pars[2];
                sol[9] = pars[3]
                print('Remove ramp %f r  + %f az + %f + %f z ' % (pars[0], pars[1], pars[2], pars[3]))

                # build total G matrix (vectorized)
                G = np.column_stack([_col_idx, _row_idx, np.ones(len(los)), _elev_flat])

            elif ivar == 0 and nfit == 1:  # was 'if' — bug fix: would run after nfit=0 block
                G = np.zeros((len(data), 5))
                G[:, 0] = rg_off
                G[:, 1] = az_off
                G[:, 2] = 1
                G[:, 3] = topo_clean
                G[:, 4] = topo_clean ** 2

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[2] = pars[0];
                sol[5] = pars[1];
                sol[8] = pars[2];
                sol[9] = pars[3];
                sol[10] = pars[4]
                print(
                    'Remove ramp %f r  + %f az + %f + %f z + %f z**2' % (pars[0], pars[1], pars[2], pars[3], pars[4]))

                # build total G matrix (vectorized)
                G = np.column_stack([_col_idx, _row_idx, np.ones(len(los)), _elev_flat, _elev_flat ** 2])

            elif ivar == 1 and nfit == 0:
                G = np.zeros((len(data), 5))
                G[:, 0] = rg_off
                G[:, 1] = az_off
                G[:, 2] = 1
                G[:, 3] = topo_clean
                G[:, 4] = topo_clean * az_off

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[2] = pars[0];
                sol[5] = pars[1];
                sol[8] = pars[2];
                sol[9] = pars[3];
                sol[11] = pars[4]
                print(
                    'Remove ramp %f r  + %f az + %f + %f z + %f z*az' % (pars[0], pars[1], pars[2], pars[3], pars[4]))

                # build total G matrix (vectorized)
                G = np.column_stack([_col_idx, _row_idx, np.ones(len(los)), _elev_flat, _elev_flat * _row_idx])

            elif ivar == 1 and nfit == 1:
                G = np.zeros((len(data), 7))
                G[:, 0] = rg_off
                G[:, 1] = az_off
                G[:, 2] = 1
                G[:, 3] = topo_clean
                G[:, 4] = topo_clean ** 2
                G[:, 5] = topo_clean * az_off
                G[:, 6] = (topo_clean * az_off) ** 2

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                # 0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2
                sol[2] = pars[0];
                sol[5] = pars[1];
                sol[8] = pars[2];
                sol[9] = pars[3];
                sol[10] = pars[4];
                sol[11] = pars[5];
                sol[12] = pars[6]
                print(
                    'Remove ramp %f r  + %f az + %f + %f z + %f z**2 + %f z*az + %f (z*az)**2' % (pars[0], pars[1],
                                                                                                  pars[2], pars[3],
                                                                                                  pars[4], pars[5],
                                                                                                  pars[6]))

                # build total G matrix (vectorized)
                _ez = _elev_flat * _row_idx
                G = np.column_stack([_col_idx, _row_idx, np.ones(len(los)), _elev_flat, _elev_flat ** 2, _ez, _ez ** 2])

    elif order == 4:
        # 0:x**3 1:x**2 2:x 3:y**3 4:y**2 5:y 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2

        if radar is None:
            G = np.zeros((len(data), 4))
            G[:, 0] = rg_off
            G[:, 1] = az_off
            G[:, 2] = rg_off * az_off
            G[:, 3] = 1

            # ramp inversion
            pars = _solve_ramp(G, data, rms)
            sol[2] = pars[0];
            sol[5] = pars[1];
            sol[7] = pars[2];
            sol[8] = pars[3]
            print('Remove ramp %f r + %f az  + %f r*az + %f' % (pars[0], pars[1], pars[2], pars[3]))

            # build total G matrix (vectorized)
            G = np.column_stack([_col_idx, _row_idx, _row_idx * _col_idx, np.ones(len(los))])

        else:
            if ivar == 0 and nfit == 0:
                G = np.zeros((len(data), 5))
                G[:, 0] = rg_off
                G[:, 1] = az_off
                G[:, 2] = rg_off * az_off
                G[:, 3] = 1
                G[:, 4] = topo_clean

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                # 0:x**3 1:x**2 2:x 3:y**3 4:y**2 5:y 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2
                sol[2] = pars[0];
                sol[5] = pars[1];
                sol[7] = pars[2];
                sol[8] = pars[3];
                sol[9] = pars[4]
                print(
                    'Remove ramp %f r, %f az  + %f r*az + %f + %f z' % (pars[0], pars[1], pars[2], pars[3], pars[4]))

                # build total G matrix (vectorized)
                G = np.column_stack([_col_idx, _row_idx, _row_idx * _col_idx, np.ones(len(los)), _elev_flat])

            elif ivar == 0 and nfit == 1:  # was 'if' — bug fix: would run after nfit=0 block
                G = np.zeros((len(data), 6))
                G[:, 0] = rg_off
                G[:, 1] = az_off
                G[:, 2] = rg_off * az_off
                G[:, 3] = 1
                G[:, 4] = topo_clean
                G[:, 5] = topo_clean ** 2

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[2] = pars[0];
                sol[5] = pars[1];
                sol[7] = pars[2];
                sol[8] = pars[3];
                sol[9] = pars[4];
                sol[10] = pars[5]
                print(
                    'Remove ramp %f r, %f az  + %f r*az + %f + %f z+ %f z**2' % (pars[0], pars[1], pars[2], pars[3],
                                                                                 pars[4], pars[5]))

                # build total G matrix (vectorized)
                G = np.column_stack(
                    [_col_idx, _row_idx, _row_idx * _col_idx, np.ones(len(los)), _elev_flat, _elev_flat ** 2])

            elif ivar == 1 and nfit == 0:
                G = np.zeros((len(data), 6))
                G[:, 0] = rg_off
                G[:, 1] = az_off
                G[:, 2] = rg_off * az_off
                G[:, 3] = 1
                G[:, 4] = topo_clean
                G[:, 5] = topo_clean * az_off

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[2] = pars[0];
                sol[5] = pars[1];
                sol[7] = pars[2];
                sol[8] = pars[3];
                sol[9] = pars[4];
                sol[11] = pars[5]
                print(
                    'Remove ramp %f r, %f az  + %f r*az + %f + %f z + %f z*az' % (pars[0], pars[1], pars[2], pars[3],
                                                                                  pars[4], pars[5]))

                # build total G matrix (vectorized) — use offset coords for consistency
                G = np.column_stack(
                    [_col_idx, _row_idx, _row_idx * _col_idx, np.ones(len(los)), _elev_flat, _elev_flat * _row_idx])

            elif ivar == 1 and nfit == 1:
                G = np.zeros((len(data), 8))
                G[:, 0] = rg_off
                G[:, 1] = az_off
                G[:, 2] = rg_off * az_off
                G[:, 3] = 1
                G[:, 4] = topo_clean
                G[:, 5] = topo_clean ** 2
                G[:, 6] = topo_clean * az
                G[:, 7] = (topo_clean * az_off) ** 2

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[2] = pars[0];
                sol[5] = pars[1];
                sol[7] = pars[2];
                sol[8] = pars[3];
                sol[9] = pars[4];
                sol[10] = pars[5];
                sol[11] = pars[6];
                sol[12] = pars[7]
                print(
                    'Remove ramp %f r, %f az  + %f r*az + %f + %f z + %f z**2 + %f z*az + %f (z*az)**2' % (pars[0],
                                                                                                           pars[1],
                                                                                                           pars[2],
                                                                                                           pars[3],
                                                                                                           pars[4],
                                                                                                           pars[5],
                                                                                                           pars[6],
                                                                                                           pars[7]))

                # build total G matrix (vectorized)
                _ez = _elev_flat * _row_idx
                G = np.column_stack(
                    [_col_idx, _row_idx, _row_idx * _col_idx, np.ones(len(los)), _elev_flat, _elev_flat ** 2, _ez,
                     _ez ** 2])

    elif order == 5:
    # 0:x**3 1:x**2 2:x 3:y**3 4:y**2 5:y 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2

        if radar is None:
            G = np.zeros((len(data), 4))
            G[:, 0] = rg_off ** 2
            G[:, 1] = rg_off
            G[:, 2] = az_off
            G[:, 3] = 1

            # ramp inversion
            pars = _solve_ramp(G, data, rms)
            sol[1] = pars[0];
            sol[2] = pars[1];
            sol[5] = pars[2];
            sol[8] = pars[3]
            print('Remove ramp %f r**2 %f r + %f az + %f' % (pars[0], pars[1], pars[2], pars[3]))

            # build total G matrix (vectorized)
            G = np.column_stack([_col2, _col_idx, _row_idx, np.ones(len(los))])

        else:
            if ivar == 0 and nfit == 0:
                G = np.zeros((len(data), 5))
                G[:, 0] = rg_off ** 2
                G[:, 1] = rg_off
                G[:, 2] = az_off
                G[:, 3] = 1
                G[:, 4] = topo_clean

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[1] = pars[0];
                sol[2] = pars[1];
                sol[5] = pars[2];
                sol[8] = pars[3];
                sol[9] = pars[4]
                print(
                    'Remove ramp %f r**2, %f r + %f az + %f + %f z' % (pars[0], pars[1], pars[2], pars[3], pars[4]))

                # build total G matrix (vectorized)
                G = np.column_stack([_col2, _col_idx, _row_idx, np.ones(len(los)), _elev_flat])

            elif ivar == 0 and nfit == 1:
                G = np.zeros((len(data), 6))
                G[:, 0] = rg_off ** 2
                G[:, 1] = rg_off
                G[:, 2] = az_off
                G[:, 3] = 1
                G[:, 4] = topo_clean
                G[:, 5] = topo_clean ** 2

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[1] = pars[0];
                sol[2] = pars[1];
                sol[5] = pars[2];
                sol[8] = pars[3];
                sol[9] = pars[4];
                sol[10] = pars[5]
                print(
                    'Remove ramp %f r**2, %f r  + %f az + %f + %f z + %f z**2' % (pars[0], pars[1], pars[2], pars[3],
                                                                                  pars[4], pars[5]))

                # build total G matrix (vectorized) — note: original had *= on row 2 (bug), should be =
                G = np.column_stack([_col2, _col_idx, _row_idx, np.ones(len(los)), _elev_flat, _elev_flat ** 2])

            elif ivar == 1 and nfit == 0:
                G = np.zeros((len(data), 6))
                G[:, 0] = rg_off ** 2  # was G[:,4] — overwritten by topo_clean (bug fix)
                G[:, 1] = rg_off
                G[:, 2] = az_off
                G[:, 3] = 1
                G[:, 4] = topo_clean
                G[:, 5] = topo_clean * az_off

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[1] = pars[0];
                sol[2] = pars[1];
                sol[5] = pars[2];
                sol[8] = pars[3];
                sol[9] = pars[4];
                sol[11] = pars[5]
                print(
                    'Remove ramp %f r**2, %f r   + %f az + %f + %f z + %f z*az' % (pars[0], pars[1], pars[2], pars[3],
                                                                                   pars[4], pars[5]))

                # build total G matrix (vectorized)
                G = np.column_stack([_col2, _col_idx, _row_idx, np.ones(len(los)), _elev_flat, _elev_flat * _row_idx])

            elif ivar == 1 and nfit == 1:
                G = np.zeros((len(data), 8))
                G[:, 0] = rg_off ** 2
                G[:, 1] = rg_off
                G[:, 2] = az_off
                G[:, 3] = 1
                G[:, 4] = topo_clean
                G[:, 5] = topo_clean ** 2
                G[:, 6] = topo_clean * az
                G[:, 7] = (topo_clean * az_off) ** 2

                # ramp inversion
                # 0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2
                pars = _solve_ramp(G, data, rms)
                sol[1] = pars[0];
                sol[2] = pars[1];
                sol[5] = pars[2];
                sol[8] = pars[3];
                sol[9] = pars[4];
                sol[10] = pars[5];
                sol[11] = pars[6];
                sol[12] = pars[7]
                print(
                    'Remove ramp %f r**2, %f r  + %f az + %f + %f z + %f z**2 + %f z*az + %f (z*az)**2' % (pars[0],
                                                                                                           pars[1],
                                                                                                           pars[2],
                                                                                                           pars[3],
                                                                                                           pars[4],
                                                                                                           pars[5],
                                                                                                           pars[6],
                                                                                                           pars[7]))

                # build total G matrix (vectorized)
                _ez = _elev_flat * _row_idx
                G = np.column_stack(
                    [_col2, _col_idx, _row_idx, np.ones(len(los)), _elev_flat, _elev_flat ** 2, _ez, _ez ** 2])

    elif order == 6:
    # 0:x**3 1:x**2 2:x 3:y**3 4:y**2 5:y 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2

        if radar is None:
            G = np.zeros((len(data), 3))
            G[:, 0] = az_off ** 2
            G[:, 1] = az_off
            G[:, 2] = 1

            # ramp inversion
            pars = _solve_ramp(G, data, rms)
            sol[4] = pars[0];
            sol[5] = pars[1];
            sol[8] = pars[2]
            print('Remove ramp %f az**2 %f az  + %f' % (pars[0], pars[1], pars[2]))

            # build total G matrix (vectorized)
            G = np.column_stack([_row2, _row_idx, np.ones(len(los))])

        else:
            if ivar == 0 and nfit == 0:
                G = np.zeros((len(data), 4))
                G[:, 0] = az_off ** 2
                G[:, 1] = az_off
                G[:, 2] = 1
                G[:, 3] = topo_clean

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[4] = pars[0];
                sol[5] = pars[1];
                sol[8] = pars[2];
                sol[9] = pars[3]
                print('Remove ramp %f az**2, %f az  + %f + %f z' % (pars[0], pars[1], pars[2], pars[3]))

                # build total G matrix (vectorized)
                G = np.column_stack([_row2, _row_idx, np.ones(len(los)), _elev_flat])

            elif ivar == 0 and nfit == 1:
                G = np.zeros((len(data), 5))
                G[:, 0] = az_off ** 2
                G[:, 1] = az_off
                G[:, 2] = 1  # was G[:,3] — bug fix: col 2 was never assigned
                G[:, 3] = topo_clean  # was G[:,4]
                G[:, 4] = topo_clean ** 2  # was G[:,5] — out of bounds on shape (N,5)

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[4] = pars[0];
                sol[5] = pars[1];
                sol[8] = pars[2];
                sol[9] = pars[3];
                sol[10] = pars[4]
                print('Remove ramp %f az**2, %f az  + %f + %f z + %f z**2' % (pars[0], pars[1], pars[2], pars[3],
                                                                                    pars[4]))

                # build total G matrix (vectorized)
                G = np.column_stack([_row2, _row_idx, np.ones(len(los)), _elev_flat, _elev_flat ** 2])

            elif ivar == 1 and nfit == 0:
                G = np.zeros((len(data), 5))
                G[:, 0] = az_off ** 2
                G[:, 1] = az_off
                G[:, 2] = 1
                G[:, 3] = topo_clean
                G[:, 4] = topo_clean * az_off

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                sol[4] = pars[0];
                sol[5] = pars[1];
                sol[8] = pars[2];
                sol[9] = pars[3];
                sol[11] = pars[4];
                print(
                    'Remove ramp %f az**2, %f az + %f + %f z + %f z*az' % (pars[0], pars[1], pars[2], pars[3], pars[4]))

                # build total G matrix (vectorized) — BUG FIX: original had G[:,4] *= i (whole-array mul each iter)
                G = np.column_stack([_row2, _row_idx, np.ones(len(los)), _elev_flat, _elev_flat * _row_idx])

            elif ivar == 1 and nfit == 1:
                G = np.zeros((len(data), 7))
                G[:, 0] = az_off ** 2
                G[:, 1] = az_off
                G[:, 2] = 1
                G[:, 3] = topo_clean
                G[:, 4] = topo_clean ** 2
                G[:, 5] = topo_clean * az_off
                G[:, 6] = (topo_clean * az_off) ** 2

                # ramp inversion
                pars = _solve_ramp(G, data, rms)
                # 0:y**3 1:y**2 2:y 3:x**3 4:x**2 5:x 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2
                sol[4] = pars[0];
                sol[5] = pars[1];
                sol[8] = pars[2];
                sol[9] = pars[3];
                sol[10] = pars[4];
                sol[11] = pars[5];
                sol[12] = pars[6]
                print(
                    'Remove ramp %f az**2, %f az + %f + %f z + %f z**2 + %f z*az + %f (z*az)**2 ' % (pars[0], pars[1],
                                                                                                     pars[2], pars[3],
                                                                                                     pars[4], pars[5],
                                                                                                     pars[6]))

                # build total G matrix (vectorized) — BUG FIX: original had G[:,5/6] *= inside loop
                _ez = _elev_flat * _row_idx
                G = np.column_stack([_row2, _row_idx, np.ones(len(los)), _elev_flat, _elev_flat ** 2, _ez, _ez ** 2])

    corr = np.dot(G, pars).reshape(mlines, mcols)
    res = los - np.dot(G, pars).flatten()
    var = np.nanstd(res)
    # plt.imshow(los.reshape(mlines,mcols))
    # plt.show()
    # plt.imshow(corr)
    # plt.show()

    return sol, corr, var, rms


def empirical_cor(kk):
    """
    Function that preapare and run empirical estimaton for each interferogram kk
    """
    # Expose variables injected by _init_worker (spawn gives a clean namespace)
    globals().update(_G)
    # deferred: avoids CoreFoundation on macOS with spawn
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    date1, date2 = date_1[kk], date_2[kk]
    idate = str(date1) + '-' + str(date2)

    if sformat == 'ROI_PAC':
        folder = 'int_' + str(date1) + '_' + str(date2) + '/'
        rscfile = int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.unw.rsc'
        infile = int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.unw'

        checkinfile(infile)
        checkinfile(rscfile)

        ds = _gdal().OpenEx(infile, allowed_drivers=["ROI_PAC"])
        # Get the band that have the data we want
        ds_band1 = ds.GetRasterBand(1)
        ds_band2 = ds.GetRasterBand(2)

        los_map = np.zeros((mlines, mcols))
        los_map[:ds.RasterYSize, :ds.RasterXSize] = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[
            :mlines, :mcols]
        # los_map[los_map==0] = float('NaN')
        lines, cols = ds.RasterYSize, ds.RasterXSize


    elif sformat == 'GTIFF':
        folder = 'int_' + str(date1) + '_' + str(date2) + '/'
        infile = int_path + folder + prefix + str(date1) + '_' + str(date2) + suffix + rlook + '.tiff'

        checkinfile(infile)

        ds = _gdal().Open(infile, _gdal().GA_ReadOnly)
        # Get the band that have the data we want
        ds_band2 = ds.GetRasterBand(1)
        lines, cols = ds.RasterYSize, ds.RasterXSize

        los_map = np.zeros((mlines, mcols))
        los_map[:ds.RasterYSize, :ds.RasterXSize] = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[
            :mlines, :mcols]
        # los_map[los_map==0] = float('NaN')

    elif sformat == 'GAMMA':
        # scfile=prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.unw.par'
        # par_file = ref
        lines, cols = gm.readpar(int_path)
        infile = int_path + prefix + str(date1) + '_' + str(date2) + suffix + rlook + '.unw'
        checkinfile(infile)
        los_map = gm.readgamma(infile, int_path)

    print('lines:{0}, cols:{1}, IFG:{2}'.format(lines, cols, idate))

    # load coherence or whatever
    spacial_mask = np.ones((mlines, mcols)) * float('NaN')

    rms_map = np.ones((mlines, mcols))
    if rmsf == 'yes':
        try:
            if sformat == 'ROI_PAC':
                rms_map[:ds.RasterYSize, :ds.RasterXSize] = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[
                    :mlines, :mcols]
                k = np.nonzero(np.logical_or(rms_map == 0.0, rms_map == 9999))
                rms_map[k] = float('NaN')
            elif sformat == 'GTIFF':
                folder = 'int_' + str(date1) + '_' + str(date2) + '/'
                rmsfile = int_path + folder + 'CNES_Coh_geo_' + str(date1) + '_' + str(date2) + rlook + '.tiff'
                checkinfile(rmsfile)
                ds = _gdal().Open(rmsfile, _gdal().GA_ReadOnly)
                ds_band1 = ds.GetRasterBand(1)
                lines, cols = ds.RasterYSize, ds.RasterXSize
                rms_map = np.zeros((mlines, mcols))
                rms_map[:ds.RasterYSize, :ds.RasterXSize] = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)[
                    :mlines, :mcols]
                k = np.nonzero(np.logical_or(rms_map == 0.0, rms_map == 9999))
                rms_map[k] = float('NaN')
            elif sformat == 'GAMMA':
                rmsfile = int_path + str(date1) + '_' + str(date2) + '.filt.cc'
                rms_map = gm.readgamma(rmsfile, int_path)
                # plt.imshow(rms_map,vmax=1,vmin=0)
                # plt.show()
        except:
            print('Coherence file cannot be read')

    # time.sleep(1.)
    # clean for estimation
    #logger.debug('Apply gaussian filter with an abitrary half-window size of 3 for estimation')
    m_filter_vals = np.copy(los_map)
    m_filter_vals[np.isnan(los_map)] = 0.
    m_lp_vals = scipy.ndimage.gaussian_filter(m_filter_vals, 3)
    # make same size array full of ones, but set to zero where there is a nan in mf
    m_filter_ones = 0 * np.copy(los_map) + 1
    m_filter_ones[np.isnan(los_map)] = 0.
    m_lp_ones = scipy.ndimage.gaussian_filter(m_filter_ones, 3)
    # find the ratio to make coefficients sum to one near nan values
    _los_map = m_lp_vals / m_lp_ones

    _los_map[np.logical_or(los_map == 0, _los_map == 0)] = float('NaN')
    _los_map[np.isnan(los_map)] = float('NaN')
    maxlos, minlos = np.nanpercentile(_los_map, perc), np.nanpercentile(_los_map, (100 - perc))

    ## CRITICAL STEP ####
    # select points for estimation only: minmax elev, los not NaN, rms<rmsthreshold ....
    index = np.nonzero(
        np.logical_and(elev_map < maxelev,
                    np.logical_and(elev_map > minelev,
                    np.logical_and(slope_map > minslope,
                    np.logical_and(los_map != 0.,
                    np.logical_and(los_map > minlos,
                    np.logical_and(los_map < maxlos,
                    np.logical_and(
                    rms_map > threshold_rms,
                    np.logical_and(
                    pix_az > ibeg,
                    np.logical_and(
                    pix_az < iend,
                    np.logical_and(
                    pix_rg > jbeg,
                    np.logical_and(
                    pix_rg < jend,
                    np.logical_and(
                    mask > threshold_mask,
                    np.logical_and(
                    ~np.isnan(
                    los_map),
                    np.logical_or(
                    pix_az < ibeg_mask,
                    pix_az > iend_mask)
                    )))))))))))))
    )

    indexref = np.nonzero(
        np.logical_and(elev_map < maxelev,
                    np.logical_and(elev_map > minelev,
                    np.logical_and(los_map != 0.,
                    np.logical_and(los_map > minlos,
                    np.logical_and(los_map < maxlos,
                    np.logical_and(
                    rms_map > threshold_rms,
                    np.logical_and(
                    pix_az > lin_start,
                    np.logical_and(
                    pix_az < lin_end,
                    np.logical_and(
                    pix_rg > col_start,
                    np.logical_and(
                    pix_rg < col_end,
                    np.logical_and(
                    mask > threshold_mask,
                    np.logical_and(
                    ~np.isnan(
                    los_map),
                    np.logical_or(
                    pix_az < ibeg_mask,
                    pix_az > iend_mask)
                    )))))))))))))

    spacial_mask[index] = np.copy(los_map[index])

    # extract range and azimuth coordinates
    temp = np.array(index).T
    az = temp[:, 0];
    rg = temp[:, 1]

    # clean maps
    los_temp = los_map.copy()
    elev_temp = elev_map.copy()
    los_clean = los_temp[index].flatten()
    los_ref = los_temp[indexref].flatten()
    rms_ref = rms_map[indexref].flatten()
    cst = np.nansum(los_ref * rms_ref) / np.nansum(rms_ref)
    print(
        'Estimation of a constant within lines {0}-{1} and cols {2}-{3}'.format(lin_start, lin_end, col_start, col_end))
    print('Average phase within ref area: {0}:'.format(cst))
    rg_ref, az_ref, topo_ref = np.nanmean(pix_rg[indexref]), np.nanmean(pix_az[indexref]), np.nanmean(
        elev_temp[indexref])

    # print('Average rg: {0}, az:{1}, topo:{2}, within ref area'.format(np.int(rg_ref), np.int(az_ref), np.int(topo_ref)))
    # sys.exit()
    elev_clean = elev_temp[index].flatten()
    rms_clean = rms_map[index].flatten()
    del los_temp, elev_temp

    # Take care to not do high polynomial estimations for short int.
    # find the begining of the image
    itemp = ibeg
    for row in range(ibeg, iend, 10):
        if np.isnan(np.nanmean(_los_map[row:row + 10, :])):
            itemp = row
        else:
            break
    del _los_map

    # 0: ref frame [default], 1: range ramp ax+b , 2: azimutal ramp ay+b,
    # 3: ax+by+c, 4: ax+by+cxy+d 5: ax**2+bx+d, 6: ay**2+by+c
    if flat > 5 and iend - itemp < .6 * (iend - ibeg):
        print('Int. too short in comparison to ref, set flat to 5')
        temp_flat = 5
    elif flat > 5 and iend - itemp < .9 * mcols:
        print('Lenght int. inferior to width, set flat to 5 and nfit to 0')
        temp_flat = 5
    else:
        temp_flat = flat

    if ivar > 0 and iend - itemp < .6 * (iend - ibeg):
        print('Int. too short in comparison to ref, set ivar and nfit to 0')
        nfit_temp = 0
        ivar_temp = 0
    else:
        nfit_temp = nfit
        ivar_temp = ivar

    # try:
    sol, corr, var, rms = estim_ramp(los_map.flatten(),
                                     los_clean[::samp], elev_clean[::samp], az[::samp], rg[::samp],
                                     temp_flat, rms_clean[::samp], nfit_temp, ivar_temp, cst, rg_ref, az_ref, topo_ref)
    # except:
    #    sol = np.zeros((13))
    #    corr = np.zeros((mlines,mcols))
    #    var = 1
    #    rg, az = az[::samp],rg[::samp]
    #    topo_clean,data,rms = elev_clean[::samp],los_clean[::samp],rms_clean[::samp]

    print('RMS: {0} '.format(var))

    if radar is not None:
        # plot phase/elevation

        fig2 = plt.figure(2, figsize=(9, 4))
        ax = fig2.add_subplot(1, 1, 1)

        z = np.linspace(np.min(elev_clean), np.max(elev_clean), 100)
        # 0:x**3 1:x**2 2:x 3:y**3 4:y**2 5:y 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2
        func = sol[0] * rg ** 3 + sol[1] * rg ** 2 + sol[2] * rg + sol[3] * az ** 3 + sol[4] * az ** 2 \
               + sol[5] * az + sol[6] * (rg * az) ** 2 + sol[7] * rg * az + sol[8] + sol[11] * az * elev_clean + \
               sol[12] * ((az * elev_clean) ** 2)

        ax.scatter(elev_clean[::5], los_clean[::5] - func[::5], s=0.05, alpha=0.1, rasterized=True)
        if nfit == 0:
            ax.plot(z, sol[8] + sol[9] * z, '-r', lw=3., label='{0:.3f}*z + {1:.3f}'.format(sol[9], sol[8]))
        else:
            ax.plot(z, sol[8] + sol[9] * z + sol[10] * z ** 2, '-r', lw=3.,
                    label='{0:.3f}*z**2 + {1:.3f}*z + {2:.3f}'.format(sol[10], sol[9], sol[8]))

        ax.set_xlabel('Elevation (m)')
        ax.set_ylabel('LOS (rad)')
        plt.legend(loc='best')
        if sformat == 'ROI_PAC':
            fig2.savefig(int_path + folder + idate + '_phase-topo.png', format='PNG')
        else:
            fig2.savefig(out_path + idate + '_phase-topo.png', format='PNG')

    _los_map = np.copy(los_map)
    _los_map[los_map == 0] = float('NaN')
    vmax = np.nanpercentile(_los_map, 98)
    vmin = np.nanpercentile(_los_map, 2)

    fig = plt.figure(3, figsize=(11, 4))

    ax = fig.add_subplot(2, 3, 1)
    cax = ax.imshow(los_map, cmap=cmap, vmax=vmax, vmin=vmin, interpolation=None)
    ax.set_title('LOS')
    plt.setp(ax.get_xticklabels(), visible=None)
    plt.setp(ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    ax = fig.add_subplot(2, 3, 2)
    cax = ax.imshow(spacial_mask, cmap=cmap, vmax=vmax, vmin=vmin, interpolation='nearest')
    ax.set_title('LOS ESTIMATION')
    plt.setp(ax.get_xticklabels(), visible=None)
    plt.setp(ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    ax = fig.add_subplot(2, 3, 3)
    cax = ax.imshow(rms_map, cmap=cmap, interpolation='nearest')
    ax.set_title('COH')
    plt.setp(ax.get_xticklabels(), visible=None)
    plt.setp(ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    # colormap correction must be same than data!!
    ax = fig.add_subplot(2, 3, 4)
    cax = ax.imshow(corr, cmap=cmap, vmax=vmax, vmin=vmin, interpolation=None)
    ax.set_title('RAMP+TOPO')
    plt.setp(ax.get_xticklabels(), visible=None)
    plt.setp(ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    # for plot we can clean
    k = np.nonzero(np.logical_or(los_map == 0., abs(los_map) > 999.))
    corr[k] = 0.

    flatlos = los_map - corr
    _los_map = np.copy(flatlos)
    _los_map[los_map == 0] = float('NaN')
    vmax = np.nanpercentile(_los_map, 98)
    vmin = np.nanpercentile(_los_map, 2)

    ax = fig.add_subplot(2, 3, 5)
    hax = ax.imshow(rms_map, cm.Greys, vmax=1, vmin=0.)
    cax = ax.imshow(flatlos, cmap=cmap, vmax=vmax, vmin=-vmax, alpha=1., interpolation=None)
    ax.set_title('CORR LOS')
    plt.setp(ax.get_xticklabels(), visible=None)
    plt.setp(ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)
    fig.tight_layout()

    if sformat == 'ROI_PAC' or sformat == 'GTIFF':
        fig.savefig(int_path + folder + idate + '_corrections.png', format='PNG')
    else:
        fig.savefig(out_path + idate + '_corrections.png', format='PNG')

    if plot == 'yes':
        plt.show()

    plt.close('all')
    del corr, flatlos

    del los_clean, rms_clean
    del elev_clean
    del az, rg
    try:
        del ds
    except:
        pass

    return iend - itemp, sol, var


def apply_cor(kk, sp, sp_inv):
    """
    Fonction that apply empirical estimatons for each interferograms kk
    """
    # Expose variables injected by _init_worker (spawn gives a clean namespace)
    globals().update(_G)
    # deferred: avoids CoreFoundation on macOS with spawn
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    # Recreate GDAL driver (SwigPyObject — not picklable, cannot be in _shared_vars)
    if sformat == 'ROI_PAC':
        driver = _gdal().GetDriverByName("roi_pac")
    elif sformat == 'GTIFF':
        driver = _gdal().GetDriverByName("GTiff")
    else:
        driver = None

    date1, date2 = date_1[kk], date_2[kk]
    folder = 'int_' + str(date1) + '_' + str(date2) + '/'
    idate = str(date1) + '-' + str(date2)

    if sformat == 'ROI_PAC':

        rscfile = int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.unw.rsc'
        infile = int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.unw'
        outfile = int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix + suffout + rlook + '.unw'
        outrsc = int_path + folder + prefix + str(date1) + '-' + str(date2) + suffix + suffout + rlook + '.unw.rsc'
        # print(infile)
        # print(outfile)

        ds = _gdal().OpenEx(infile, allowed_drivers=["ROI_PAC"])
        # Get the band that have the data we want
        ds_band1 = ds.GetRasterBand(1)
        ds_band2 = ds.GetRasterBand(2)
        los_map = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        rms_map = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        lines, cols = ds.RasterYSize, ds.RasterXSize

    elif sformat == 'GTIFF':

        infile = int_path + prefix + str(date1) + '_' + str(date2) + suffix + rlook + '.geo.unw.tif'
        outfile = out_path + prefix + str(date1) + '_' + str(date2) + suffix + suffout + rlook + '.geo.unw.tif'
        ds = _gdal().Open(infile, _gdal().GA_ReadOnly)
        ds_band2 = ds.GetRasterBand(1)
        los_map = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        rms_map = np.ones((ds.RasterYSize, ds.RasterXSize))
        lines, cols = ds.RasterYSize, ds.RasterXSize

        if rmsf == 'yes':
            rmsfile = int_path + str(date1) + '_' + str(date2) + 'geo.cc.tif'
            rmsfile = int_path + str(date1) + '_' + str(date2) + '.geo.cc.tif'
            ds = _gdal().Open(rmsfile, _gdal().GA_ReadOnly)
            ds_band1 = ds.GetRasterBand(1)
            rms_map = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
        else:
            rms_map = np.ones((lines, cols))

    elif sformat == 'GAMMA':

        infile = int_path + prefix + str(date1) + '_' + str(date2) + suffix + rlook + '.unw'
        outfile = out_path + prefix + str(date1) + '_' + str(date2) + suffix + suffout + rlook + '.unw'
        # par_file = ref
        lines, cols = gm.readpar(int_path)
        los_map = gm.readgamma(infile, int_path)

        if rmsf == 'yes':
            rmsfile = int_path + str(date1) + '_' + str(date2) + '.filt.cc'
            rms_map = gm.readgamma(rmsfile, int_path)
        else:
            rms_map = np.ones((lines, cols))

    print('mlines:{}, mcols:{}, int:{}:'.format(lines, cols, idate))

    # compute correction
    # Coordinates offset by ibeg/jbeg — must match G_full in estim_ramp
    rg = np.tile(np.arange(cols) - jbeg, (lines, 1))
    az = np.tile(np.arange(lines) - ibeg, (cols, 1)).T
    z = np.zeros((lines, cols))
    z = elev_map[:lines, :cols]

    # 0:x**3 1:x**2 2:x 3:y**3 4:y**2 5:y 6:xy**2 7:xy 8:cst 9:z 10:z**2 11:yz 12:yz**2
    # apply correction
    sol = sp_inv[kk, 3:]
    corr_inv = sol[0] * rg ** 3 + sol[1] * rg ** 2 + sol[2] * rg + sol[3] * az ** 3 + sol[4] * az ** 2 \
               + sol[5] * az + sol[6] * (rg * az) ** 2 + sol[7] * rg * az + sol[8] + sol[9] * z + sol[10] * (z ** 2) \
               + sol[11] * az * z + sol[12] * ((az * z) ** 2)

    sol = sp[kk, 3:]
    corr = sol[0] * rg ** 3 + sol[1] * rg ** 2 + sol[2] * rg + sol[3] * az ** 3 + sol[4] * az ** 2 \
           + sol[5] * az + sol[6] * (rg * az) ** 2 + sol[7] * rg * az + sol[8] + sol[9] * z + sol[10] * (z ** 2) \
           + sol[11] * az * z + sol[12] * ((az * z) ** 2)

    # apply corr
    flatlos = los_map - corr_inv

    # recompute ref frame
    zone = flatlos[lin_start:lin_end, col_start:col_end]
    amp = rms_map[lin_start:lin_end, col_start:col_end]
    minlos, maxlos = np.nanpercentile(zone[abs(zone) > 1.e-6], 5.), np.nanpercentile(zone[abs(zone) > 1.e-6], 95.)
    index = np.nonzero(
        np.logical_and(zone > minlos,
                       np.logical_and(zone < maxlos,
                                      amp > threshold_rms,
                                      )))
    # give small weights to small uncoherent vectors
    cst = np.nansum(zone[index] * amp[index]) / np.nansum(amp[index])
    if np.isnan(cst):
        pass
    else:
        flatlos = flatlos - cst
        corr_inv = corr_inv + cst
        corr = corr + cst
    print('Remove reference frame: {}'.format(cst))

    # reset to 0 areas where no data (might change after time series inversion?)
    flatlos[np.isnan(los_map)], rms_map[np.isnan(los_map)] = 0.0, 0.0
    flatlos[np.isnan(flatlos)], rms_map[np.isnan(flatlos)] = 0.0, 0.0
    rms_map[np.isnan(rms_map)], flatlos[np.isnan(rms_map)] = 0.0, 0.0
    flatlos[los_map == 0], rms_map[los_map == 0] = 0.0, 0.0

    if sformat == 'ROI_PAC':
        dst_ds = driver.Create(outfile, cols, lines, 2, _gdal().GDT_Float32)
        dst_band1 = dst_ds.GetRasterBand(1)
        dst_band2 = dst_ds.GetRasterBand(2)
        dst_band1.WriteArray(rms_map, 0, 0)
        dst_band2.WriteArray(flatlos, 0, 0)
        shutil.copy(rscfile, outrsc)
        dst_band1.FlushCache()
        dst_band2.FlushCache()
        del dst_ds, ds

    elif sformat == 'GTIFF':
        dst_ds = driver.Create(outfile, cols, lines, 1, _gdal().GDT_Float32)
        dst_band2 = dst_ds.GetRasterBand(1)
        dst_band2.WriteArray(flatlos, 0, 0)
        dst_ds.SetGeoTransform(gt)
        dst_ds.SetProjection(proj)
        dst_band2.FlushCache()
        del dst_ds, ds

    elif sformat == 'GAMMA':
        fid = open(outfile, 'wb')
        flatlos.flatten().astype('>f4').tofile(fid)
        fid.close()

    fig = plt.figure(5, figsize=(9, 6))

    _los_map = np.copy(los_map)
    _los_map[los_map == 0] = float('NaN')
    vmax, vmin = np.nanpercentile(_los_map, 99), np.nanpercentile(_los_map, 1)

    ax = fig.add_subplot(1, 4, 1)
    cax = ax.imshow(los_map, cmap=cmap, vmax=vmax, vmin=vmin, alpha=1, interpolation=None)
    ax.set_title(str(date1) + '_' + str(date2))
    plt.setp(ax.get_xticklabels(), visible=None)
    plt.setp(ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    ax = fig.add_subplot(1, 4, 2)
    cax = ax.imshow(corr, cmap=cmap, vmax=vmax, vmin=vmin, alpha=1, interpolation=None)
    ax.set_title('Best-fit Ramp')
    plt.setp(ax.get_xticklabels(), visible=None)
    plt.setp(ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    ax = fig.add_subplot(1, 4, 3)
    cax = ax.imshow(corr_inv, cmap=cmap, vmax=vmax, vmin=vmin, alpha=1, interpolation=None)
    ax.set_title('Reconstructed')
    plt.setp(ax.get_xticklabels(), visible=None)
    plt.setp(ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)

    _los_map = np.copy(flatlos)
    _los_map[_los_map == 0] = float('NaN')
    # flatlos[los_map==0] = float('NaN')
    vmax, vmin = np.nanpercentile(_los_map, 99), np.nanpercentile(_los_map, 1)

    ax = fig.add_subplot(1, 4, 4)
    cax = ax.imshow(_los_map, cmap=cmap, vmax=vmax, vmin=vmin, alpha=1, interpolation=None)
    ax.set_title('Flattened LOS')
    plt.setp(ax.get_xticklabels(), visible=None)
    plt.setp(ax.get_yticklabels(), visible=None)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c)
    fig.tight_layout()

    if sformat == 'ROI_PAC':
        fig.savefig(int_path + folder + idate + '_reconstruc_corrections.png', format='PNG')
    else:
        fig.savefig(out_path + idate + '_reconstruc_corrections.png', format='PNG')

    if plot == 'yes':
        plt.show()

    plt.close('all')
    del los_map, rms_map


#####################################################################################
# INIT LOG
#####################################################################################

if __name__ == '__main__':

    # logging.basicConfig(level=logging.INFO,\
    # logging.basicConfig(level=logging.INFO, \
    #                     format='line %(lineno)s -- %(levelname)s -- %(message)s')
    # logger = logging.getLogger('invert_ramp_topo_unw.log')

    #####################################################################################
    # READ INPUT PARAM
    #####################################################################################

    parser = argparse.ArgumentParser(
        description='Estimates atmospheric phase/elevation correlations and/or '
                    'azimuthal/range ramp polynomial coefficients on unwrapped '
                    'interferograms (ROI_PAC/GAMMA/GTIFF). Temporal inversion of '
                    'all coefficients with strong weight for small temporal baselines '
                    'and large-coverage interferograms. Reconstructs the empirical '
                    'phase correction.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument('--int_list', required=True, metavar='PATH',
                        help='Text file with list of interferogram dates in two columns: date1 date2')
    parser.add_argument('--int_path', default='./', metavar='PATH',
                        help='Relative path to input interferograms directory')
    parser.add_argument('--out_path', default=None, metavar='PATH',
                        help='Output directory path (defaults to int_path)')
    parser.add_argument('--ref_zone', default=None, metavar='l0,l1,c0,c1',
                        help='Lines and cols where phase is set to zero: lin_start,lin_end,col_start,col_end')
    parser.add_argument('--prefix', default='', metavar='VALUE',
                        help='Prefix for interferogram filename: $prefix$date1-$date2$suffix_$rlookrlks.unw')
    parser.add_argument('--suffix', default='', metavar='VALUE',
                        help='Suffix for interferogram filename')
    parser.add_argument('--rlook', default=None, metavar='VALUE',
                        help='Look factor for interferogram filename (appended as _Nrlks)')
    parser.add_argument('--flat', default=0, type=int, choices=range(7), metavar='0-6',
                        help='Spatial ramp order. 0:ref frame, 1:range ax+b, 2:azimuth ay+b, '
                             '3:ax+by+c, 4:ax+by+cxy+d, 5:ax²+bx+az+d, 6:ay²+by+c')
    parser.add_argument('--topofile', default=None, metavar='PATH',
                        help='Path to radar_look.hgt DEM file. Enables phase/elevation correlation')
    parser.add_argument('--ivar', default=0, type=int, choices=[0, 1], metavar='0/1',
                        help='Phase/elevation relationship: 0=f(elevation), 1=f(azimuth×elevation)')
    parser.add_argument('--nfit', default=0, type=int, choices=[0, 1], metavar='0/1',
                        help='Fit degree for elevation term: 0=linear, 1=quadratic')
    parser.add_argument('--ref', default=None, metavar='PATH',
                        help='Reference image to define format and dimensions (required if --topofile not given)')
    parser.add_argument('--format', default='ROI_PAC', choices=['ROI_PAC', 'GAMMA', 'GTIFF'],
                        help='Input file format')
    parser.add_argument('--tsinv', default='no', choices=['yes', 'no'],
                        help='If yes, invert corrected phase into time series')
    parser.add_argument('--estim', default='yes', choices=['yes', 'no'],
                        help='If yes, run estimation; otherwise read existing correction matrices')
    parser.add_argument('--mask', default=None, metavar='PATH',
                        help='Mask file in .r4 format — keep pixels > threshold_mask')
    parser.add_argument('--threshold_mask', default=-1, type=float, metavar='VALUE',
                        help='Threshold on mask file')
    parser.add_argument('--cohpixel', default='no', choices=['yes', 'no'],
                        help='Use amplitude/coherence image to weight and mask pixels')
    parser.add_argument('--threshold_coh', default=None, type=float, metavar='VALUE',
                        help='Threshold on cohpixel file (required when --cohpixel yes)')
    parser.add_argument('--ibeg_mask', default=None, type=int, metavar='VALUE',
                        help='Start line for the mask region')
    parser.add_argument('--iend_mask', default=None, type=int, metavar='VALUE',
                        help='End line for the mask region')
    parser.add_argument('--perc', default=99., type=float, metavar='VALUE',
                        help='Percentile threshold for LOS pixel outlier rejection')
    parser.add_argument('--perc_topo', default=99., type=float, metavar='VALUE',
                        help='Percentile threshold for elevation pixel outlier rejection')
    parser.add_argument('--min_topo', default=None, type=float, metavar='VALUE',
                        help='Minimum elevation for estimation (overrides perc_topo)')
    parser.add_argument('--max_topo', default=None, type=float, metavar='VALUE',
                        help='Maximum elevation for estimation (overrides perc_topo)')
    parser.add_argument('--perc_slope', default=98., type=float, metavar='VALUE',
                        help='Percentile threshold for slope-elevation outlier rejection')
    parser.add_argument('--samp', default=1, type=int, metavar='VALUE',
                        help='Undersampling factor for empirical estimation')
    parser.add_argument('--nproc', default=4, type=int, metavar='NB_CORES',
                        help='Number of parallel processes')
    parser.add_argument('--plot', default='no', choices=['yes', 'no'],
                        help='If yes, plot figures for each interferogram (forces nproc=1)')
    parser.add_argument('--suffix_output', default='_corrunw', metavar='VALUE',
                        help='Output filename suffix: $prefix$date1-$date2$suffix$suffix_output')
    parser.add_argument('--crop_emp', default=None, metavar='l0,l1,c0,c1',
                        help='Region of interest for spatial estimation: lin_start,lin_end,col_start,col_end')
    parser.add_argument('--cpt', default=None, metavar='PATH',
                        help='Colormap file for plots (.txt, loaded as LinearSegmentedColormap)')

    args = parser.parse_args()

    # --- Assign variables from parsed arguments ---

    int_list = args.int_list
    int_path = args.int_path if args.int_path.endswith('/') else args.int_path + '/'

    if args.out_path is None:
        out_path = int_path
    else:
        out_path = args.out_path if args.out_path.endswith('/') else args.out_path + '/'
        makedirs(out_path)

    prefix = args.prefix
    suffix = args.suffix
    rlook = ('_' + args.rlook + 'rlks') if args.rlook is not None else ''

    flat = args.flat
    ivar = args.ivar
    nfit = args.nfit

    if args.topofile is not None and os.path.exists(args.topofile):
        radar = args.topofile
    else:
        radar = None
        if args.ivar != 0 or args.nfit != 0:
            print('No valid topographic file given. Empirical phase/topo will not be performed')

    if args.ref is not None:
        ref = args.ref
    else:
        ref = None
    if ref is None and radar is None:
        print('Argument error: need --ref or --topofile')
        sys.exit(1)

    maskfile = args.mask if (args.mask is not None and os.path.exists(args.mask)) else None
    threshold_mask = args.threshold_mask

    rmsf = args.cohpixel
    if rmsf == 'yes':
        if args.threshold_coh is None:
            print('--cohpixel yes requires --threshold_coh. Exit!')
            sys.exit(1)
        threshold_rms = args.threshold_coh
    else:
        threshold_rms = -1

    tsinv = args.tsinv
    sformat = args.format
    estim = args.estim

    ibeg_mask = args.ibeg_mask if args.ibeg_mask is not None else np.inf
    iend_mask = args.iend_mask if args.iend_mask is not None else -np.inf

    perc = args.perc
    perc_topo = args.perc_topo
    perc_slope = args.perc_slope
    samp = args.samp
    nproc = args.nproc
    suffout = args.suffix_output

    if args.plot == 'yes':
        plot = 'yes'
        print('--plot yes: setting nproc to 1')
        nproc = 1
        try:
            if environ.get("TERM", "").startswith("screen"):
                matplotlib.use('Agg')
        except Exception:
            pass
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable
    else:
        plot = 'no'
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable

    if args.cpt is None:
        try:
            cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
            cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs + "roma.txt")).reversed()
        except Exception:
            cmap = cm.rainbow
    else:
        try:
            cmap = LinearSegmentedColormap.from_list(
                args.cpt.split("/")[-1].split('.')[0], np.loadtxt(args.cpt))
        except Exception:
            cmap = cm.rainbow

    #####################################################################################
    # INITIALISE
    #####################################################################################

    # print()
    # read int
    date_1, date_2 = np.loadtxt(int_list, comments="#", unpack=True, usecols=(0, 1), dtype='i,i')
    Nifg = len(date_1)
    print("number of interferogram: {}".format(Nifg))

    # list dates
    im = [];
    bt = []
    for date1, date2 in zip(date_1, date_2):
        if date1 not in im: im.append(date1)
        if date2 not in im: im.append(date2)
    nmax = len(im)
    print("number of image: {} ".format(nmax))
    imd = date2dec(im)
    cst = np.copy(imd[0])
    # compute temporal baseline for TS inv
    for i in range((nmax)):
        bt.append(imd[i] - cst)

    # load ref to define mlines, mcols and format
    if ref is not None:
        ds_extension = os.path.splitext(ref)[1]
        if sformat == 'GTIFF':
            geotiff = ref
            georef = _gdal().Open(ref)
            gt = georef.GetGeoTransform()
            proj = georef.GetProjection()
            driver = _gdal().GetDriverByName('GTiff')
            ds = _gdal().Open(ref, _gdal().GA_ReadOnly)
            mlines, mcols = ds.RasterYSize, ds.RasterXSize
        elif sformat == 'ROI_PAC':
            driver = _gdal().GetDriverByName("roi_pac")
            ds = _gdal().OpenEx(ref, allowed_drivers=["ROI_PAC"])
            mlines, mcols = ds.RasterYSize, ds.RasterXSize
        elif sformat == 'GAMMA':
            from parsers import gamma as gm

            # par_file = ref
            mlines, mcols = gm.readpar(int_path)

    # laod elevation map
    if radar is not None:
        extension = os.path.splitext(radar)[1]
        if extension == '.r4':
            fid = open(radar, 'r')
            elevi = np.fromfile(fid, dtype=np.float32)
            elev_map = elevi.reshape(mlines, mcols)
            fid.close()

        else:
            if sformat == 'GTIFF':
                geotiff = args.topofile
                georef = _gdal().Open(radar)
                gt = georef.GetGeoTransform()
                proj = georef.GetProjection()
                driver = _gdal().GetDriverByName('GTiff')
                ds = _gdal().Open(radar, _gdal().GA_ReadOnly)
                ds_band2 = ds.GetRasterBand(2)
                mlines, mcols = ds.RasterYSize, ds.RasterXSize
                elev_map = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
                del ds

            elif sformat == 'ROI_PAC':
                driver = _gdal().GetDriverByName("roi_pac")
                ds = _gdal().OpenEx(radar, allowed_drivers=["ROI_PAC"])
                ds_band2 = ds.GetRasterBand(2)
                mlines, mcols = ds.RasterYSize, ds.RasterXSize
                elev_map = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
                del ds

            elif sformat == 'GAMMA':
                from parsers import gamma as gm

                # par_file = ref
                mlines, mcols = gm.readpar()
                elev_map = gm.readgamma(radar)

            if args.max_topo != None:
                maxelev = float(args.max_topo)
            else:
                maxelev = np.nanpercentile(elev_map, perc_topo)
            if args.min_topo != None:
                minelev = float(args.min_topo)
            else:
                minelev = np.nanpercentile(elev_map, 100 - perc_topo)
            print('Max-Min topography for empirical estimation: {0:.1f}-{1:.1f}'.format(maxelev, minelev))

        # compute slope
        toposmooth = scipy.ndimage.gaussian_filter(elev_map, 3.)
        Py, Px = np.gradient(toposmooth)
        slope_map = np.sqrt(Px ** 2 + Py ** 2)
        minslope = np.nanpercentile(slope_map, 100 - perc_slope)
        print('Min relief for empirical estimation: {0}'.format(minslope))

        fig = plt.figure(0, figsize=(12, 8))

        ax = fig.add_subplot(1, 2, 1)
        cax = ax.imshow(toposmooth, cm.RdBu_r, vmin=minelev, vmax=maxelev)
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.set_title('Smoothed DEM', fontsize=6)

        ax = fig.add_subplot(1, 2, 2)
        cax = ax.imshow(slope_map, cm.RdBu_r, vmin=minslope, vmax=np.nanpercentile(slope_map, perc_slope))
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.set_title('Mask Slope bellow: {0:.2f}'.format(minslope), fontsize=8)

        if plot == 'yes':
            plt.show()

    else:
        maxelev, minelev = 1., -1
        elev_map = np.zeros((mlines, mcols))
        slope_map = np.zeros((mlines, mcols))
        minslope = -1

    # open mask file
    if maskfile is not None:
        fid = open(maskfile, 'r')
        maski = np.fromfile(fid, dtype=np.float32)[:mcols * mlines]
        mask = maski.reshape((mlines, mcols))
        k = np.nonzero(mask < threshold_mask)
        spacial_mask = np.copy(mask)
        spacial_mask[k] = float('NaN')

        if plot == 'yes':
            fig = plt.figure(0, figsize=(5, 4))
            ax = fig.add_subplot(1, 1, 1)
            cax = ax.imshow(spacial_mask, cmap=cmap)
            ax.set_title('Mask')
            plt.setp(ax.get_xticklabels(), visible=None)
            divider = make_axes_locatable(ax)
            c = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(cax, cax=c)
            plt.show()
            fid.close()
    else:
        mask = np.zeros((mlines, mcols))
        threshold_mask = -1

    if args.crop_emp is None:
        crop_emp = [0, mlines, 0, mcols]
    else:
        crop_emp = list(map(float, args.crop_emp.replace(',', ' ').split()))
        print('Crop empirical estimation between lines {}-{} and cols {}-{}'.format(
            int(crop_emp[0]), int(crop_emp[1]), int(crop_emp[2]), int(crop_emp[3])))
    ibeg, iend, jbeg, jend = int(crop_emp[0]), int(crop_emp[1]), int(crop_emp[2]), int(crop_emp[3])

    if args.ref_zone is None:
        lin_start, lin_end, col_start, col_end = 0, mlines, 0, mcols
    else:
        _rz = list(map(int, args.ref_zone.replace(',', ' ').split()))
        try:
            lin_start, lin_end, col_start, col_end = _rz[0], _rz[1], _rz[2], _rz[3]
        except ValueError:
            lin_start, lin_end = _rz[0], _rz[1]
            col_start, col_end = 0, mcols

    #####################################################################################
    # MAIN
    #####################################################################################

    # Shared header/format for coefficient files
    _coeff_hdr = ('#date1   |   dates2   |   Lenght   |   y**3   |   y**2   |   y'
                  '   |   **3   |   x**2   |   x   |   xy**2   |   xy   |   cst'
                  '   |   z   |   z**2   |   z*az   |   z**2*az')
    _coeff_fmt = ('%i', '%i', '%.10f', '%.10f', '%.10f', '%.10f', '%.10f',
                  '%.10f', '%.10f', '%.10f', '%.10f', '%.10f', '%.10f', '%.10f', '%.10f', '%.10f')

    # extract range and azimuth coordinates from ref or radar file
    pix_az, pix_rg = np.indices((mlines, mcols))

    # initialise full vector correction
    # 16 values: date1 dates2  mlines y**3 y**2 y x**3 x**2 x xy**2 xy cst z z**2 z*az az*z**2
    spint = np.zeros((Nifg, 16))
    # 16 cols: date1 date2 mlines y**3 y**2 y x**3 x**2 x xy**2 xy cst z z**2 z*az az*z**2

    # fill dates
    spint[:, 0], spint[:, 1] = date_1, date_2

    # initilise correction cube
    rmsint = np.zeros((Nifg, 3))
    rmsint[:, 0], rmsint[:, 1] = date_1, date_2

    date1_err, date2_err = [], []

    if estim == 'yes':

        print()
        #########################################
        print('#################################')
        print('Empirical estimations')
        print('#################################')
        #########################################
        print()

        output = []
        # go
        with TimeIt():
            # for kk in range(Nifg):
            work = range(Nifg)
            # Shared vars for spawn workers (excludes non-picklable SwigPyObjects like driver)
            _shared_vars = {k: globals()[k] for k in (
                'date_1', 'date_2', 'sformat', 'int_path', 'out_path', 'prefix', 'suffix',
                'rlook', 'suffout', 'mlines', 'mcols', 'elev_map', 'slope_map', 'mask',
                'pix_az', 'pix_rg', 'flat', 'ivar', 'nfit', 'samp', 'radar', 'perc',
                'maxelev', 'minelev', 'minslope', 'threshold_rms', 'threshold_mask', 'rmsf',
                'ibeg', 'iend', 'jbeg', 'jend', 'lin_start', 'lin_end', 'col_start', 'col_end',
                'ibeg_mask', 'iend_mask', 'plot', 'cmap', 'gt', 'proj',
            ) if k in globals()}
            with poolcontext(processes=nproc,
                             initializer=_init_worker,
                             initargs=(_shared_vars,)) as pool:
                results = pool.map(empirical_cor, work)
            output.append(results)
            # for work in range(Nifg):
            #     empirical_cor(work)

            for kk in range(Nifg):
                lenght, sol, rms = output[0][kk]
                # save size int to use as weight in the temporal inversion
                spint[kk, 2] = lenght
                # fill correction matrix
                spint[kk, 3:] = sol
                rmsint[kk, 2] = rms

        # print(spint)

        # save spint
        np.savetxt(_outpath('list_coeff_ramps{}.txt'.format(suffout)), spint,
                   header=_coeff_hdr, fmt=_coeff_fmt)

        # save rms
        np.savetxt(_outpath('rms{}.txt'.format(suffout)), rmsint,
                   header='# date1   |   dates2   |   RMS', fmt=('%i', '%i', '%.8f'))

    #####################################################################################

    print()
    print('read input files list_coeff_ramps.txt')

    # load spint
    date_1, date_2, length, a, b, c, d, e, f, g, h, i, j, k, l, m = np.loadtxt(
        _outpath('list_coeff_ramps{}.txt'.format(suffout)),
        comments='#', unpack=True, dtype='i,i,f,f,f,f,f,f,f,f,f,f,f,f,f,f')
    spint = np.vstack([date_1, date_2, length, a, b, c, d, e, f, g, h, i, j, k, l, m]).T
    rec_spint = np.copy(spint)

    # load rms
    rmsint = np.loadtxt(_outpath('rms{}.txt'.format(suffout)), comments='#')

    ####################################################################################

    # init inv solution
    spint_inv = np.zeros((np.shape(spint)))

    if tsinv == 'yes':

        print()
        #########################################
        print('#################################')
        print('Temporal inversion of all coefficients')
        print('#################################')
        #########################################
        print()

        G_ = np.zeros((Nifg, nmax))
        deltat = np.zeros((Nifg))
        im_arr = np.array(im)
        for k in range(Nifg):
            n1 = np.where(im_arr == date_1[k])[0][0]
            n2 = np.where(im_arr == date_2[k])[0][0]
            G_[k, n1] = -1
            G_[k, n2] = 1
            deltat[k] = abs(bt[n2] - bt[n1])

        # 1) create weight based on temporal baseline: give stronger weight to short temporal baselines
        # , where we dont expect def.
        w1 = np.exp(-(deltat / 2.))

        # 2) create a weight based on the size of the int: give a stronger weight to long interferograms
        w2 = length / mlines

        # 3) rms weight
        w3 = np.exp(-rmsint[:, 2] / np.nanpercentile(rmsint[:, 2], 80)) + 0.01

        wf = open(_outpath('list_weigths{}.txt'.format(suffout)), 'w')
        print('Weights for temporal inversion:')
        for i in range(len(w3)):
            print(int(rmsint[i, 0]), int(rmsint[i, 1]), rmsint[i, 2], w1[i], w2[i], w3[i])
            wf.write("%i %i %f %f %f %f\n" % (int(rmsint[i, 0]), int(rmsint[i, 1]), rmsint[i, 2], w1[i], w2[i], w3[i]))
        wf.close()
        print('Weights saved in: list_weights.txt')

        # compute summ of weights
        sig_ = 1. / w1 + 1. / w2 + 1. / w3

        for j in range(3, np.shape(spint_inv)[1]):

            d = np.zeros(((Nifg + 1)))
            sig = np.ones(((Nifg + 1)))
            G = np.zeros(((Nifg + 1), nmax))

            d[:Nifg] = spint[:, j]
            G[:Nifg, :nmax] = G_
            G[-1, 0] = 1  # ini phi first image to 0
            sig[:Nifg] = sig_

            try:
                x0 = lst.lstsq(G, d)[0]
                _func = lambda x: np.sum(((np.dot(G, x) - d) / sig) ** 2)
                _fprime = lambda x: 2 * np.dot(G.T / sig, (np.dot(G, x) - d) / sig)
                pars = opt.fmin_slsqp(_func, x0, fprime=_fprime, iter=20000, full_output=True, iprint=0)[0]

                # reconstruct corr for selected int
                spint_inv[:, j] = np.dot(G, pars)[:Nifg]

            except:
                pass

        spint_inv[:, :3] = spint[:, :3]
        # save spint_inv
        np.savetxt(_outpath('list_coeff_ramps{}_inv.txt'.format(suffout)), spint_inv,
                   header=_coeff_hdr,
                   fmt=('%i', '%i', '%.10f', '%.10f', '%.10f', '%.10f', '%.10f',
                        '%.10f', '%.10f', '%.10f', '%.10f', '%.10f', '%.10f', '%.12f', '%.12f', '%.12f'))
    ####################################################################################

    print()
    #########################################
    print('#################################')
    print('APPLY CORRECTION AND SAVE NEW INT.')
    print('#################################')
    #########################################
    print()

    # go
    with TimeIt():
        work = range(Nifg)
        # for j in range(Nifg):
        #     apply_cor(j, spint, spint_inv)
        # Shared vars for spawn workers (excludes non-picklable SwigPyObjects like driver)
        _shared_vars = {k: globals()[k] for k in (
            'date_1', 'date_2', 'sformat', 'int_path', 'out_path', 'prefix', 'suffix',
            'rlook', 'suffout', 'mlines', 'mcols', 'elev_map', 'slope_map', 'mask',
            'pix_az', 'pix_rg', 'flat', 'ivar', 'nfit', 'samp', 'radar', 'perc',
            'maxelev', 'minelev', 'minslope', 'threshold_rms', 'threshold_mask', 'rmsf',
            'ibeg', 'iend', 'jbeg', 'jend', 'lin_start', 'lin_end', 'col_start', 'col_end',
            'ibeg_mask', 'iend_mask', 'plot', 'cmap', 'gt', 'proj',
        ) if k in globals()}
        with poolcontext(processes=nproc,
                         initializer=_init_worker,
                         initargs=(_shared_vars,)) as pool:
            pool.map(partial(apply_cor, sp=spint, sp_inv=spint_inv), work)
