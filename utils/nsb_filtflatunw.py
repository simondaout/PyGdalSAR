#!/usr/bin/env python2
# -*- coding: utf-8 -*-
############################################

############################################
# Author        : Simon DAOUT (Oxford)
############################################

from __future__ import print_function

# gdal
import gdal, shutil
gdal.UseExceptions()
# system
from os import path, environ, system, chdir, remove
# plot
import subprocess
import matplotlib
if environ["TERM"].startswith("screen"):
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from datetime import datetime
# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided
# scipy
import scipy.optimize as opt
import scipy.linalg as lst
import logging


##################################################################################
###  INITIALISE
##################################################################################

# init logger 
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('filtcorunw')

# init figures
nfigure = 0

# init collections 
import collections
Process = collections.namedtuple('Process', 'name  do', rename=True)
Param = collections.namedtuple('Param', 'name', rename=True)
IFG = collections.namedtuple('IFG', 'date1 date2 look prefix suffix', rename=True)
Image = collections.namedtuple('Image', 'date decimal_date temporal_baseline', rename=True)

##################################################################################
###  Class and Functions
##################################################################################

def date2dec(dates):
    ''' Transform dates %Y%m%d to decimal dates'''
    dates  = np.atleast_1d(dates)
    times = []
    for date in dates:
        x = datetime.strptime('{}'.format(date),'%Y%m%d')
        dec = float(x.strftime('%j'))/365.1
        year = float(x.strftime('%Y'))
        times.append(year + dec)
    return times

class Job:
    """ Create a class of Jobs: 
    Job list is: look_int replace_amp filter flat_range flat_topo flat_model unw add_back"""
    names = 'look_int replace_amp filter flat_range flat_az flat_topo flat_model unw add_back'.split()

    def __init__(self, do_list):
        self.do_list = do_list
        try:
            self._processes = [Process(name,do) for (name,do) in zip(self.names,self.do_list)]
        except ValueError as error:
            logger.warning(error)

    # create spetial methods len() and getititem for Job class
    def __len__(self):
        return len(self._processes)
    
    def __getitem__(self ,pos):
        return self._processes[pos]

    def __call__(self, pos):
        return self._processes[post].name()

    def add_job(self,job):
        self._processes._make(job)

    def replace_job(self,pos,value):
        self._processes[pos]._replace(value)


class PileInt:
    def __init__(self,dates1, dates2, prefix, suffix, look, filterstyle ,dir):
        self.dates1, self.dates2 = dates1, dates2
        self.dir = dir
        self.filterstyle = filterstyle

        self.Nifg=len(self.dates1)
        print("number of interferogram: ",self.Nifg)

        try:
            self._ifgs = [IFG(date1,date2,look,prefix,suffix) for (date1,date2) in zip(self.dates1,self.dates2)]
        except ValueError as error:
            logger.warning(error)

    # create spetial methods len() and getititem for PileInt class
    def __len__(self):
        return len(self._ifgs)

    def __getitem__(self,kk):
        return self._ifgs[kk]

    def getlook(self,kk):
        return self._ifgs[kk].look

    def getfix(self,kk):
        return self._ifgs[kk].prefix, self._ifgs[kk].suffix

    def getname(self,kk):
        ''' Return interfergram file name '''
        return str(self._ifgs[kk].prefix) + str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + str(self._ifgs[kk].suffix) + '_' +  self._ifgs[kk].look + 'rlks.int'

    def getfilt(self,kk):
        ''' Return interfergram file name '''
        return 'filt' + str(self.filterstyle) + '_' + str(self._ifgs[kk].prefix) + str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + str(self._ifgs[kk].suffix) + '_' +  self._ifgs[kk].look + 'rlks.int'

    def getcor(self,kk):
        ''' Return cohrence file name '''
        return  str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + '_' +  self._ifgs[kk].look + 'rlks.cor'

    def getpath(self,kk):
        ''' Return path ifg dir '''
        return  self.dir + 'int_' + str(self.dates1[kk]) + '_' + str(self.dates2[kk]) 

    def updatelook(self,kk,newlook):
        self._ifgs[kk] = self._ifgs[kk]._replace(look=newlook)  

    def updatefix(self,kk,newprefix, newsuffix):
        self._ifgs[kk] = self._ifgs[kk]._replace(prefix=str(newprefix))
        self._ifgs[kk] = self._ifgs[kk]._replace(suffix=str(newsuffix))

    def info(self):
        print('List of interferograms:')
        print ([self.getname(kk) for kk in range(self.Nifg)])
        print()

class PileImages:
    def __init__(self,dates1,dates2):
        self.dates1, self.dates2 = dates1, dates2

        # define list of images 
        im = []; bt = []
        for date1,date2 in zip(self.dates1,self.dates2):
            if date1 not in im: im.append(date1)
            if date2 not in im: im.append(date2)
        self.Nimages=len(im)
        print("number of image: ",self.Nimages)
        imd = date2dec(im)
        cst = np.copy(imd[0])
        for i in xrange((self.Nimages)):
            bt.append(imd[i]-cst)
        del cst

        try:
            self._images = [Image(date,dec,baseline) for (date,dec,baseline) in zip(im,imd,bt)]
        except ValueError as error:
            logger.warning(error)
    
    # create spetial methods len() and getititem for PileImages class
    def __len__(self):
        return len(self._images)

    def __getitem__(self,kk):
        return self._images[kk]

    def info(self):
        print('List of Images:')
        print ([self._images[kk] for kk in xrange(self.Nimages)])
        print()

class FiltFlatUnw:
    """ Create a class FiltFlatUnw defining all the post-procesing functions 
    list of parameters defined in the proc file: ListInterfero, SARMasterDir, IntDir, Rlooks_int, Rlooks_unw, prefix, suffix ,
    nfit_range, hresh_amp_range, nfit_az, thresh_amp_az, filterstyle,SWwindowsize, SWamplim,
    seedx, seedy """

    def __init__(self, params):
        (self.ListInterfero, self.SARMasterDir, self.IntDir,
        self.Rlooks_int, self.Rlooks_unw, self.prefix, self.suffix, 
        self.nfit_range, self.thresh_amp_range,
        self.nfit_az, self.thresh_amp_az,
        self.filterstyle,self.SWwindowsize, self.SWamplim,
        self.seedx, self.seedy,
        ) = map(str, params)

        # initiliase number of looks
        self.look = self.Rlooks_int
        self.rlook = int(int(self.Rlooks_unw) - int(self.Rlooks_int))

        # initilise radar file
        self.dem = self.SARMasterDir + '_' + self.Rlooks_int + '.hgt'

        # define list of interferograms
        dates1,dates2=np.loadtxt(self.ListInterfero,comments="#",unpack=True,usecols=(0,1),dtype='i,i')
        self.stack = PileInt(dates1,dates2,self.prefix,self.suffix,self.look,self.filterstyle, self.IntDir)
        self.stack.info()
        self.Nifg = len(self.stack)

        # define list images
        self.images = PileImages(dates1,dates2)
        self.images.info()
        self.Nimages = len(self.images)

    def look(self):
        ''' Look Radar function '''
        chdir(self.SARMasterDir) 

        r= subprocess.call("look.pl "+str(self.dem)+" "+str(self.rlook)+" >> log_look.txt" , shell=True)
        if r != 0:
            logger.warning(r)
            logger.warning(' Can''t look file {0} in {1} look'.format(self.dem,self.rlook))

        self.dem = self.SARMasterDir + '_' + self.Rlooks_unw + '.hgt'

    def look_int(self,kk):
        ''' Look function '''

        # need to be in the corect dir for stupid perl scripts
        chdir(self.stack.getpath(kk))

        infile =  self.stack.getname(kk)
        corfile =  self.stack.getcor(kk)
        print(corfile) 
        logger.debug('Look file {0} in {1} look'.format(infile,self.rlook))

        r= subprocess.call("look.pl "+str(infile)+" "+str(self.rlook)+" >> log_look.txt" , shell=True)
        if r != 0:
            logger.warning(r)
            logger.warning(' Can''t look file {0} in {1} look'.format(infile,self.rlook))
        
        r = subprocess.call("look.pl "+str(corfile)+" "+str(self.rlook)+" >> log_look.txt", shell=True)
        if r != 0:
            logger.warning(r)
            logger.warning(' Can''t look file {0} in {1} look'.format(corfile,self.rlook))

        # update looks
        self.stack.updatelook(kk,self.Rlooks_unw)

    def replace_amp(self, kk):

        # update looks
        self.stack.updatelook(kk,self.Rlooks_unw)

        infile = self.stack.getpath(kk) + '/'+ self.stack.getname(kk)
        rscfile = infile + '.rsc'
        corfile = self.stack.getpath(kk) + '/' + self.stack.getcor(kk)
        logger.debug('Replace Amplitude by Coherence on IFG: {0}'.format(infile))

        # update names
        prefix, suffix = self.stack.getfix(kk)
        newprefix = 'coh_' + prefix
        self.stack.updatefix(kk,newprefix,suffix)
        outfile = self.stack.getpath(kk) + '/'+ self.stack.getname(kk)
        outrsc = outfile + '.rsc'

        tmp = self.stack.getpath(kk) + '/tmp'
        phs = self.stack.getpath(kk) + '/phs'
        cor = self.stack.getpath(kk) + '/cor'

        # Open
        ds_int = gdal.Open(infile, gdal.GA_ReadOnly)
        driver = ds_int.GetDriver()
        width = ds_int.RasterXSize
        print("> Width:     ", ds_int.RasterXSize)

        # check if nor done
        if path.exists(outfile) is False:
            if path.exists(corfile):

                logger.debug('Replace Amplitude by Cohrence for IFG: {}'.format(infile))
                r1 = subprocess.call("rmg2mag_phs "+str(corfile)+" "+str(tmp)+" "+str(cor)+" "+str(width), shell=True)
                r2 = subprocess.call("cpx2mag_phs "+str(infile)+" "+str(tmp)+" "+str(phs)+" "+str(width), shell=True)
                r3 = subprocess.call("mag_phs2cpx cor phs "+str(outfile)+" "+str(width), shell=True)
                if (r1 or r2 or r3) != 0:
                    logger.warning(r1, r2, r3)
                    logger.warning('Replace Amplitude by Cohrence for IFG: {}'.format(infile))
                
                remove(tmp); remove(phs); remove(cor)
                del tmp, phs, cor
                shutil.copy(rscfile,outrsc)

            else:
                logger.warning('Coherence file does not exit...')

        else:
            logger.debug('Replace Amplitude by Cohrence for IFG: {} already done'.format(infile))
        
        del ds_int 

    def filter(self, kk):
        ''' Filter function '''

        # need to be in the corect dir for stupid perl scripts
        chdir(self.stack.getpath(kk))

        infile = self.stack.getname(kk)
        inrsc = infile + '.rsc'
        inbase = path.splitext(infile)[0]
        corfile = self.stack.getcor(kk)
        corbase = path.splitext(corfile)[0]
        filtfile = self.stack.getfilt(kk)
        filtrsc = filtfile + '.rsc'
        filtbase = path.splitext(filtfile)[0]

        logger.debug('Filter filter {0} with {1} filter type'.format(infile,self.filterstyle))
        try:
            system("nsb_SWfilter.pl "+str(inbase)+" "+str(filtbase)+" "+str(corbase)+" "+str(self.SWwindowsize)+" "+str(self.SWamplim)+" "+str(self.filterstyle))
        except:
            logger.warning('Cant filter {0} with {1} filter type'.format(infile,self.filterstyle))

        if path.exists(filtrsc) == False:
            shutil.copy(inrsc,filtrsc)

    def flat_range(self,kk):
        ''' Faltten Range function '''

        # need to be in the corect dir for stupid perl scripts
        chdir(self.stack.getpath(kk))

        infile = self.stack.getname(kk)
        inrsc = infile + '.rsc'
        corfile = self.stack.getcor(kk)
        filtfile = self.stack.getfilt(kk)

        # update names
        prefix, suffix = self.stack.getfix(kk)
        newsuffix = suffix + '_flatr'
        self.stack.updatefix(kk,prefix,newsuffix)
        outfile = self.stack.getname(kk)
        outrsc = outfile + '.rsc'
        filtout = self.stack.getfilt(kk)

        if path.exists(filtfile) == False:
            logger.debug('{0} does not exist'.format(filtfile))
            # call filter function
            eval(self.filter(kk))

        if path.exists(outfile) == False:
            logger.debug('Flatten range on IFG: {0}'.format(infile))
            r = subprocess.call("flatten_range "+str(infile)+" "+str(filtfile)+" "+str(outfile)+" "+str(filtout)+" "+str(self.nfit_range)+" "+str(self.thresh_amp_range)+"  >> log_flatenrange.txt", shell=True)
            if r != 0:
                logger.warning("Flatten range failed for IFG: {0}".format(infile))
                logger.warning(r)
            else:
                shutil.copy(inrsc,outrsc) 
        else:
            logger.debug('Flatten range on IFG: {0} already done'.format(infile))

    def flat_az(self,kk):
        ''' Faltten Azimuth function '''
        infile = self.stack.getpath(kk) + '/'+ self.stack.getname(kk)
        corfile = self.stack.getpath(kk) + '/' + self.stack.getcor(kk)
        filtfile = self.stack.getpath(kk) + '/'+ self.stack.getfilt(kk)

        # update names
        prefix, suffix = self.stack.getfix(kk)
        newsuffix = suffix + '_flataz'
        self.stack.updatefix(kk,newprefix,suffix)
        outfile = self.stack.getpath(kk) + '/'+ self.stack.getname(kk)
        filtout = self.stack.getpath(kk) + '/'+ self.stack.getfilt(kk)

        if path.exists(filtfile) == False:
            logger.debug('{0} does not exist'.format(filtfile))
            # call filter function
            eval(self.filter(kk))

        if path.exists(outfile) == False:
            logger.debug('Flatten azimuth on IFG: {0}'.format(infile))
            r = subprocess.call("flatten_az "+str(infile)+" "+str(filtfile)+" "+str(outfile)+" "+str(filtout)+" "+str(fit)+" "+str(threshold), shell=True)
            if r != 0:
                logger.warning("Flatten azimuth failed for int. {0}-{1}".format(date1,date2))
            else:
                shutil.copy(rscfile,outrsc)
        else:
            logger.debug('Flatten azimuth on IFG: {0} already done'.format(infile))


    def flat_topo(kk, prefix, suffix):
        ''' Faltten topo function '''
        infile = self.stack.getpath(kk) + '/'+ self.stack.getname(kk)
        corfile = self.stack.getpath(kk) + '/' + self.stack.getcor(kk)
        filtfile = self.stack.getpath(kk) + '/'+ self.stack.getfilt(kk)

        # update names
        prefix, suffix = self.stack.getfix(kk)
        newsuffix = suffix + '_flatz'
        self.stack.updatefix(kk,newprefix,suffix)
        outfile = self.stack.getpath(kk) + '/'+ self.stack.getname(kk)
        filtout = self.stack.getpath(kk) + '/'+ self.stack.getfilt(kk)

        dem = 

        if path.exists(filtfile) == False:

        if path.exists(filtfile) == False:
            logger.debug('{0} does not exist'.format(filtfile))
            # call filter function
            eval(self.filter(kk))

        if path.exists(outfile) == False:
            logger.debug('Flatten topo on IFG: {0}'.format(infile))
            r = subprocess.call("flatten_az "+str(infile)+" "+str(filtfile)+" "+str(outfile)+" "+str(filtout)+" "+str(fit)+" "+str(threshold), shell=True)
            if r != 0:
                logger.warning("Flatten topo failed for int. {0}-{1}".format(date1,date2))
            else:
                shutil.copy(rscfile,outrsc)
        else:
            logger.debug('Flatten topo on IFG: {0} already done'.format(infile))

    def flat_model(kk, prefix, suffix):
        return

    def unwrapping(kk):
        date1, date2 = self.date_1[kk], self.date_2[kk]
        return 

    def addback(kk):
        date1, date2 = self.date_1[kk], self.date_2[kk]
        return 


##################################################################################
###  READ IMPUT PARAMETERS
##################################################################################

# # input parameters 
home='/home/cometraid14/daouts/work/tibet/qinghai/processing/Sentinel/iw1/'
IntDir=path.abspath(home)+'/'+'test/'
ListInterfero=path.abspath(home)+'/'+'interf_pair_test.rsc'
SARMasterDir=path.abspath(home)+'/'+'20160608'

Rlooks_int=int(2)
Rlooks_unw=int(4)
nfit_range = -1
thresh_amp_range = 0.3
nfit_az = 0
thresh_amp_az = 0.3
filterstyle='SWc'
SWamplim=0.05
SWwindowsize=8
seedx=268
seedy=1766
prefix = '' 
suffix = '_sd'
nproc=1

####################
# Test Process List
####################

# Job list: look_int replace_amp   filter     flat_range    flat_az    flat_topo  flat_model      unw         add_back
do_list =   [True,     True,       True,       True,       False,       False,      False,      False,      False] 
jobs = Job(do_list)

print('List of Post-Processing Jobs:')
for job in jobs:
    print(job)
print()

###########
#   MAIN 
###########

postprocess = FiltFlatUnw(
        [ListInterfero,SARMasterDir,IntDir,
        Rlooks_int, Rlooks_unw, prefix, suffix, 
        nfit_range, thresh_amp_range,
        nfit_az, thresh_amp_az,
        filterstyle,SWwindowsize, SWamplim,
        seedx,seedy]
        ) 

# loop over the processes
for p in jobs:
    # check if the process has to be done
    # print(getattr(p,'name'))
    job = getattr(p,'name')
    if p.do is True:
        print('Run {} ....'.format(p))
        [eval('postprocess.{0}({1})'.format(job,kk)) for kk in range(postprocess.Nifg)]
        print()

# [postprocess.call(job,kk) for kk in range(postprocess.Nifg)]
# work = [kk for kk in range(postprocess.Nifg)]
# pool = multiprocessing.Pool(nproc)
# pool.map(look, work)
# pool.close()      


