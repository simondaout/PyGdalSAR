#!/usr/bin/env python2
# -*- coding: utf-8 -*-
############################################

from __future__ import print_function
from numpy.lib.stride_tricks import as_strided
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import os, sys
import os.path as op
import glob

def check_required(required, params):
    for r in required:
        if r not in params:
            return False
    return True

def _parseParameterFile(filename):
    with open(filename, mode='r') as par:
        text = par.read().splitlines()
        raw_segs = [line.split() for line in text if ':' in line]
    
    return dict((i[0][:-1],i[1:]) for i in raw_segs)


def _getParameters(path, log=False):
    required_slc = ['nlines', 'width']
    required_int = ['range_samples','azimuth_lines']

    path = op.dirname(op.realpath(path))
    par_files = glob.glob('%s/*par' % path)

    for file in par_files:
        params = _parseParameterFile(file)
        # print (params)

        if check_required(required_int, params)\
              or check_required(required_slc, params):
                if not log:
                    self._log.info('Found parameter file %s' % file)
                return params

    raise ImportError(
                    'Parameter file does not hold required parameters')             


def readpar(par_file):
        """
        :params filename: Gamma software parameter file
        :type filename: str
        :param par_file: Corresponding parameter (:file:`*par`) file.
                         (optional)
        :type par_file: str
        :returns: Import dictionary
        :rtype: dict
        :raises: ImportError
        """
        
        params = _getParameters(par_file, log=True)
        nrows = int(params['range_samples'][0]) 
        nlines = int(params['azimuth_lines'][0])
        return nlines,nrows

def readgamma(filename,par_file):

        params = _getParameters(par_file, log=True)
        nrows = int(params['range_samples'][0]) 
        nlines = int(params['azimuth_lines'][0])
        
        return np.fromfile(filename, dtype='>f4').reshape(nlines, nrows)


