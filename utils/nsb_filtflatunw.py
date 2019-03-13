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
from os import path, environ, system, chdir, remove, getcwd, listdir
# plot
import subprocess
import matplotlib
if environ["TERM"].startswith("screen"):
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
from pylab import *
import matplotlib.pyplot as plt
from datetime import datetime
# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided
import logging
from multiprocessing import Pool
from contextlib import contextmanager
from functools import wraps, partial
from itertools import repeat
# import subprocess

##################################################################################
###  INITIALISE
##################################################################################

# init logger 
logging.basicConfig(level=logging.INFO,\
      format='%(asctime)s -- %(levelname)s -- %(message)s')
logger = logging.getLogger('filtflatunw_log.log')

#filename='filtflatunw_log.log'

# init collections 
import collections
Process = collections.namedtuple('Process', 'name')
IFG = collections.namedtuple('IFG', 'date1 date2 look prefix suffix width length')
Image = collections.namedtuple('Image', 'date decimal_date temporal_baseline')

##################################################################################
###  Extras functions and context maganers
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

# create context manager for change dirt
class Cd(object):
    def __init__(self,dirname):
        self.dirname = dirname
    def __enter__(self):
        self.curdir = getcwd()
        logger.debug('Enter {0}'.format(self.dirname))
        chdir(self.dirname)
    def __exit__(self, type, value, traceback):
        logger.debug('Enter {0}'.format(self.curdir))
        chdir(self.curdir)

def checkinfile(file):
    if path.exists(file) is False:
        logger.critical("File: {0} not found, Exit !".format(file))
        print("File: {0} not found, Exit !".format(file))
        # print(listdir('./'))
        sys.exit()


##################################################################################
###  RUN FONCTION
##################################################################################

# create generator for pool
# @contextmanager
# def poolcontext(*arg, **kargs):
#     pool = Pool(*arg, **kargs)
#     yield pool
#     pool.terminate()
#     pool.join()

def go(config,job,nproc):
    ''' RUN processing function '''
    
    with TimeIt():
        work = range(config.Nifg)
        
        if nproc > 1:
            pool = Pool(processes=nproc)
            results = pool.map(partial(eval(job), config), work)
            results.join()

            # with poolcontext(processes=nproc) as pool:
                # results = pool.map(partial(eval(job), config), work)
        else:
            map(eval(job), repeat(config, len(work)) , work)


##################################################################################
###  Define Job, IFG, Images and FiltFlatUnw classes  
##################################################################################

class Job():
    """ Create a class of Jobs to be run: 
    Job list is: erai look_int replace_amp filterSW filterROI flat_range flat_topo flat_model colin unwrapping add_model_back add_atmo_back add_ramp_back """

    def __init__(self, names):
        self.names = names.split()
        try:
            logging.info('Define list of processess')
            self._processes = [Process(name) for name in self.names]
        except ValueError as error:
            logger.critical(error)
            self.info()
            sys.exit()

    # create spetial methods len() and getititem for Job class
    def __len__(self):
        return len(self._processes)
    
    def __getitem__(self ,pos):
        return self._processes[pos]

    def __call__(self, pos):
        return self._processes[post].name()

    def add_job(self,job):
        self._processes = self._processes._make(job)

    def replace_job(self,pos,value):
        self._processes = self._processes[pos]._replace(value)

    def info(self):
        print('List of possible Jobs:') 
        print('erai look_int replace_amp filter flat_range flat_topo flat_model colin \
            unwrapping add_model_back add_atmo_back add_ramp_back')
        print('Choose them in the order that you want')

class PileInt:
    def __init__(self,dates1, dates2, prefix, suffix, look, filterstyle ,dir):
        self.dates1, self.dates2 = dates1, dates2
        self.dir = dir
        self.filterstyle = filterstyle

        self.Nifg=len(self.dates1)
        print("number of interferogram: ",self.Nifg)

        # init width and length to zero as long as we don't need it
        width, length = 0, 0

        try:
            logger.info('Define IFG list')
            self._ifgs = [IFG(date1,date2,look,prefix,suffix,width,length) for (date1,date2) in zip(self.dates1,self.dates2)]
        except ValueError as error:
            logger.critical(error)
            self.exit()

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
        return str(self._ifgs[kk].prefix) + str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + str(self._ifgs[kk].suffix) + '_' +  \
        self._ifgs[kk].look + 'rlks'

    def getfiltSW(self,kk):
        ''' Return interfergram file name '''
        return 'filt' + str(self.filterstyle) + '_' + str(self._ifgs[kk].prefix) + str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + \
        str(self._ifgs[kk].suffix) + '_' +  self._ifgs[kk].look + 'rlks'

    def getfiltROI(self,kk):
        ''' Return interfergram file name '''
        return 'filt_' + str(self._ifgs[kk].prefix) + str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + str(self._ifgs[kk].suffix) \
        + '_' +  self._ifgs[kk].look + 'rlks'

    def getcor(self,kk):
        ''' Return cohrence file name '''
        return  str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + '_' +  self._ifgs[kk].look + 'rlks.cor'

    def getsize(self,kk):
        ''' Return width and length IFG '''
        return  str(self._ifgs[kk].width),  str(self._ifgs[kk].length)

    def getpath(self,kk):
        ''' Return path ifg dir '''
        return  str(self.dir + 'int_' + str(self.dates1[kk]) + '_' + str(self.dates2[kk])) 

    def getstratfile(self,kk):
        ''' Return stratified file name '''
        return  str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + '_strat_' +  self._ifgs[kk].look + 'rlks'

    def updatelook(self,kk,newlook):
        self._ifgs[kk] = self._ifgs[kk]._replace(look=newlook)  

    def updatesize(self,kk,newwidth,newlength):
        self._ifgs[kk] = self._ifgs[kk]._replace(width=newwidth,length=newlength)

    def updatefix(self,kk,newprefix, newsuffix):
        self._ifgs[kk] = self._ifgs[kk]._replace(prefix=str(newprefix))
        self._ifgs[kk] = self._ifgs[kk]._replace(suffix=str(newsuffix))

    def info(self):
        print('List of interferograms:')
        print ([self._ifgs[kk] for kk in xrange(self.Nifg)])
        # print ([self.getname(kk) for kk in range(self.Nifg)])
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
            logger.info('Define list of Images')
            self._images = [Image(date,dec,baseline) for (date,dec,baseline) in zip(im,imd,bt)]
        except ValueError as error:
            logger.critical(error)
    
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
    nfit_range, hresh_amp_range, nfit_az, thresh_amp_az, filterstyle,SWwindowsize, SWamplim, filterStrength, nfit_atmo,thresh_amp_atmo, ivar, z_ref,
    seedx, seedy.threshold_unw, unw_method
    Additional parameters not in the proc file (yet?): ibeg_mask, iend_mask, jbeg_mask, jend_mask (default: 0.)
    defining the boundary of the mask zone for emprical estimations
    suffix, preffix: define name of the interferogram at the start of the processes
    model: model to be removed from wrapped interferograms (default: None)
    """

    def __init__(self, params, prefix='', suffix='_sd', ibeg_mask=0, iend_mask=0, jbeg_mask=0, jend_mask=0, model=None):
        (self.ListInterfero, self.SARMasterDir, self.IntDir,
        self.Rlooks_int, self.Rlooks_unw, 
        self.nfit_range, self.thresh_amp_range,
        self.nfit_az, self.thresh_amp_az,
        self.filterstyle,self.SWwindowsize, self.SWamplim,
        self.filterStrength,
        self.nfit_atmo,self.thresh_amp_atmo, self.ivar, self.z_ref,
        self.seedx, self.seedy,self.threshold_unw,self.unw_method,
        ) = map(str, params)

        # initialise prefix, suffix
        self.prefix, self.suffix = prefix, suffix

        # initiliase number of looks
        self.look = self.Rlooks_int
        self.rlook = int(int(self.Rlooks_unw) - int(self.Rlooks_int))

        # initialise model to be removed from wrapped int
        self.model = model
        self.strat = False

        # initilise radar file
        self.dem =  self.SARMasterDir + '/'+  'radar_' + self.Rlooks_int + 'rlks.hgt'

        # mask empirical estimations
        self.ibeg_mask, self.iend_mask, self.jbeg_mask, self.jend_mask = ibeg_mask, iend_mask, jbeg_mask, jend_mask

        # define list of interferograms
        dates1,dates2=np.loadtxt(self.ListInterfero,comments="#",unpack=True,usecols=(0,1),dtype='i,i')
        self.stack = PileInt(dates1,dates2,self.prefix,self.suffix,self.look,self.filterstyle, self.IntDir)
        self.stack.info()
        self.Nifg = len(self.stack)

        # define list images
        self.images = PileImages(dates1,dates2)
        self.images.info()
        self.Nimages = len(self.images)

    # def update(self, prefix, siffix, look):
    #     self.stack = self.stack._replace(prefix=str(prefix), suffix=str(suffix), look=str(look))

##################################################################################
###  Define Job functions 
##################################################################################

def look_file(config,file):
    ''' Look function 
    Requiered parameters:  Rlooks_int, Rlooks_unw
    '''
    
    dirname, filename = path.split(path.abspath(file)) 
    with Cd(dirname):
        logger.info("look.pl "+str(filename)+" "+str(config.rlook))
        r= subprocess.call("look.pl "+str(filename)+" "+str(config.rlook)+" >> log_look.txt" , shell=True)
        if r != 0:
            logger.critical(' Can''t look file {0} in {1} look'.format(filename,config.rlook))
            print(look_file.__doc__)

def erai(config,kk):
    return

def computesize(config,file):
    ''' Extract width anf length with gdal
    '''
    try:
        dirname, filename = path.split(path.abspath(file))
        with Cd(dirname):
            ds_int = gdal.Open(filename, gdal.GA_ReadOnly)
            driver = ds_int.GetDriver()
            return ds_int.RasterXSize, ds_int.RasterYSize
    except OSError as err:
        print("OS error: {0}".format(err))
    except ValueError as error:
        logger.critical(error)
        print(computesize.__doc__)

def replace_amp(config, kk):
    ''' Replace amplitude by coherence'''

    with Cd(config.stack.getpath(kk)):
        infile = config.stack.getname(kk)+ '.int'; checkinfile(infile)
        rscfile = infile + '.rsc'
        corfile = config.stack.getcor(kk); checkinfile(corfile)
        logger.info('Replace Amplitude by Coherence on IFG: {0}'.format(infile))

        # compute width and length
        width,length = computesize(config,infile)
        config.stack.updatesize(kk,width,length)

        # update names
        prefix, suffix = config.stack.getfix(kk)
        newprefix = 'coh_' + prefix
        config.stack.updatefix(kk,newprefix,suffix)
        outfile = config.stack.getname(kk) + '.int'
        outrsc = outfile + '.rsc' 
        shutil.copy(rscfile,outrsc)

        # check if not done
        if path.exists(outfile) is False:

            try:
                logger.info("rmg2mag_phs "+str(corfile)+" tmp cor "+str(width))
                r1 = subprocess.call("rmg2mag_phs "+str(corfile)+" tmp cor "+str(width)+"  >> log_replaceAMP.txt", shell=True)
                if r1 != 0:
                    logger.critical(r1)
                    logger.critical('Replace Amplitude by Cohrence for IFG: {} Failed!'.format(infile))
                    sys.exit()

                logger.info("cpx2mag_phs "+str(infile)+" tmp2 phs "+str(width))
                r2 = subprocess.call("cpx2mag_phs "+str(infile)+" tmp2 phs "+str(width)+"  >> log_replaceAMP.txt", shell=True)
                if (r2) != 0:
                    logger.critical(r2)
                    logger.critical('Replace Amplitude by Cohrence for IFG: {} Failed!'.format(infile))
                    sys.exit()

                logger.info("mag_phs2cpx cor phs "+str(outfile)+" "+str(width))
                r3 = subprocess.call("mag_phs2cpx cor phs "+str(outfile)+" "+str(width)+"  >> log_replaceAMP.txt", shell=True)
                if (r3) != 0:
                    logger.critical(r3)
                    logger.critical('Replace Amplitude by Cohrence for IFG: {} Failed!'.format(infile))
                    sys.exit()

                remove('tmp'); remove('tmp2'); remove('phs'); remove('cor')

            except:
                print(replace_amp.__doc__)
                sys.exit()

        else:
            logger.warning('{0} exists, assuming OK'.format(outfile))

def filterSW(config, kk):
    ''' Filter SW function form Doin et. al. 2011
    Requiered proc parameters: SWwindowsize, SWamplim, filterstyle
    '''

    with Cd(config.stack.getpath(kk)):    
        inbase = config.stack.getname(kk) 
        inrsc = inbase + '.int.rsc'
        infile = inbase + '.int'; checkinfile(infile)
        corfile = config.stack.getcor(kk); checkinfile(corfile)
        corbase = path.splitext(corfile)[0]
        filtbase = config.stack.getfiltSW(kk)
        filtrsc = filtbase + '.int.rsc'
        outfile = filtbase + '.int'

        if path.exists(outfile) == False:
            logger.info('Filter {0} with {1} filter type'.format(infile,config.filterstyle))
            r = subprocess.call("nsb_SWfilter.pl "+str(inbase)+" "+str(filtbase)+" "+str(corbase)\
                    +" "+str(config.SWwindowsize)+" "+str(config.SWamplim)+" "+str(config.filterstyle), shell=True)
            if r != 0:
                logger.critical('Filtering {0} with {1} filter type Failed!'.format(infile,config.filterstyle))
                print(filterSW.__doc__)
                sys.exit()

            if path.exists(filtrsc) == False:
                shutil.copy(inrsc,filtrsc)

        else:
            logger.warning('{0} exists, assuming OK'.format(outfile))

def filterROI(config, kk):
    ''' ROI-PAC Filter function
    Requiered proc file parameter: filterStrength
    '''

    with Cd(config.stack.getpath(kk)):

        infile = config.stack.getname(kk) + '.int'; checkinfile(infile)
        inrsc = infile + '.rsc'
        # get width and compute if not already done
        width,length =  config.stack.getsize(kk)
        if (int(width) == 0) or (int(length) == 0):
            width,length = computesize(config,infile)
            config.stack.updatesize(kk,width,length)

        filtfile = config.stack.getfiltROI(kk) + '.int'
        filtrsc = filtfile + '.rsc'
        if path.exists(filtrsc) == False:
            shutil.copy(inrsc,filtrsc)

        if path.exists(filtfile) == False:
            logger.info("myadapt_filt "+str(infile)+" "+str(filtfile)+" "+str(width)+" 0.25"+" "+str(config.filterStrength))
            r = subprocess.call("myadapt_filt "+str(infile)+" "+str(filtfile)+" "\
                    +str(width)+" 0.25"+" "+str(config.filterStrength)+"  >> log_filtROI.txt", shell=True)
            if r != 0:
                logger.critical('Failed filtering {0} with ROI-PAC adaptative filter Failed!'.format(infile))
                print(filterROI.__doc__)
                sys.exit()
        else:
            logger.warning('{0} exists, assuming OK'.format(filtfile))
        
def flat_range(config,kk):
    ''' Function flatten range  on wrapped phase  (See Doin et al., 2015)
    Requiered proc file parameters: nfit_range, thresh_amp_range
    Estimation done on filterSW file
    '''

    with Cd(config.stack.getpath(kk)):
        infile = config.stack.getname(kk) + '.int'; checkinfile(infile)
        inrsc = infile + '.rsc'
        corfile = config.stack.getcor(kk); checkinfile(corfile)
        filtfile = config.stack.getfiltSW(kk) + '.int'; checkinfile(filtfile)

        # update names
        prefix, suffix = config.stack.getfix(kk)
        newsuffix = suffix + '_flatr'
        config.stack.updatefix(kk,prefix,newsuffix)
        outfile = config.stack.getname(kk) + '.int' 
        outrsc = outfile + '.rsc'
        filtout = config.stack.getfiltSW(kk)
        shutil.copy(inrsc,outrsc)

        if path.exists(filtfile) == False:
            logger.warning('{0} does not exist'.format(filtfile))
            # call filter function
            config.filterSW(kk)

        if path.exists(outfile) == False:
            logger.info("flatten_range "+str(infile)+" "+str(filtfile)+" "+str(outfile)+" "+str(filtout)+\
                " "+str(config.nfit_range)+" "+str(config.thresh_amp_range))
            r = subprocess.call("flatten_range "+str(infile)+" "+str(filtfile)+" "+str(outfile)+" "+str(filtout)+\
                " "+str(config.nfit_range)+" "+str(config.thresh_amp_range)+"  >> log_flatenrange.txt", shell=True)
            if r != 0:
                logger.critical("Flatten range failed for IFG: {0} Failed!".format(infile))
                print(flat_range.__doc__) 
                sys.exit()
        else:
            logger.warning('{0} exists, assuming OK'.format(outfile))

def flat_az(config,kk):
    ''' Function flatten azimuth  on wrapped phase  (See Doin et al., 2015)
        Requiered proc file parameters: nfit_az, thresh_amp_az
        Estimation done on filterSW file
    '''

    with Cd(config.stack.getpath(kk)):
        ''' Faltten Azimuth function '''
        infile = config.stack.getname(kk) + '.int'; checkinfile(infile)
        inrsc = infile + '.rsc'
        corfile = config.stack.getcor(kk); checkinfile(corfile)
        filtfile = config.stack.getfiltSW(kk) + '.int'

        # update names
        prefix, suffix = config.stack.getfix(kk)
        newsuffix = suffix + '_flataz'
        config.stack.updatefix(kk,prefix,newsuffix)
        outfile = config.stack.getname(kk) + '.int' 
        filtout = config.stack.getfiltSW(kk) + '.int'
        outrsc = outfile + '.rsc'
        shutil.copy(inrsc,outrsc)

        if path.exists(filtfile) == False:
            logger.warning('{0} does not exist'.format(filtfile))
            # call filter function
            eval(config.filterSW(kk))

        if path.exists(outfile) == False:
            logger.info("flatten_az "+str(infile)+" "+str(filtfile)+" "+str(outfile)+" "+str(filtout)+\
                " "+str(nfit_az)+" "+str(thresh_amp_az))
            r = subprocess.call("flatten_az "+str(infile)+" "+str(filtfile)+" "+str(outfile)+" "+str(filtout)+\
                " "+str(nfit_az)+" "+str(thresh_amp_az), shell=True)
            if r != 0:
                logger.critical("Flatten azimuth failed for int. {0}-{1} Failed!".format(date1,date2))
                print(flat_az.__doc__)
                sys.exit()
        else:
            logger.warning('{0} exists, assuming OK'.format(outfile))

def flat_topo(config, kk):
    ''' Function flatten atmosphere on wrapped phase  (See Doin et al., 2015)
    Requiered proc file parameters: nfit_atmo, ivar, z_ref, thresh_amp_atmo
    Estimation done on filterSW file
    Plot phase/topo in *_phase-topo.png file
    '''

    with Cd(config.stack.getpath(kk)):
        infile = config.stack.getname(kk) + '.int';  checkinfile(infile)
        corfile = config.stack.getcor(kk)
        filtfile = config.stack.getfiltSW(kk) + '.int';  checkinfile(filtfile)
        stratfile = config.stack.getstratfile(kk) + '.unw'

        # filt must be done before changing name
        if path.exists(filtfile) == False:
            logger.info('{0} does not exist'.format(filtfile))
            # call filter function
            config.filterSW(kk)

        # update names
        prefix, suffix = config.stack.getfix(kk)
        newsuffix = suffix + '_flatz'
        config.stack.updatefix(kk,prefix,newsuffix)
        outfile = config.stack.getname(kk) + '.int'
        filtout = config.stack.getfiltSW(kk) + '.int'

        # check width IFG
        width,length = config.stack.getsize(kk)
        if (int(width) == 0) or (int(length) == 0):
            width,length = computesize(config,infile)
            config.stack.updatesize(kk,width,length)

    # look dem if necessary
    w,l = computesize(config,config.dem)

    if int(w) != int(width):
        logger.warning('IFG:{0} and DEM file are not the same size: {0}'.format(infile))
        look_file(config,config.dem)
        # update DEM
        config.dem = config.SARMasterDir + '/'+  'radar_' + config.Rlooks_unw + 'rlks.hgt'
        
    with Cd(config.stack.getpath(kk)):

        if path.exists(outfile) == False:
            logger.info("flatten_topo "+str(infile)+" "+str(filtfile)+" "+str(config.dem)+" "+str(outfile)+" "+str(filtout)\
                +" "+str(config.nfit_atmo)+" "+str(config.ivar)+" "+str(config.z_ref)+" "+str(config.thresh_amp_atmo)+" "+\
                str(stratfile))
            r = subprocess.call("flatten_topo "+str(infile)+" "+str(filtfile)+" "+str(config.dem)+" "+str(outfile)+" "+str(filtout)\
                +" "+str(config.nfit_atmo)+" "+str(config.ivar)+" "+str(config.z_ref)+" "+str(config.thresh_amp_atmo)+" "+\
                str(stratfile)+" >> log_flattopo.txt", shell=True)
            if r != 0:
                logger.critical("Flatten topo failed for int. {0} Failed!".format(infile))
                print(flat_topo.__doc__)
                sys.exit()
        else:
            print('{0} exists, assuming OK'.format(outfile))
        
        inrsc = infile + '.rsc'
        outrsc = outfile + '.rsc'
        # print(inrsc,outrsc)
        filtrsc = filtout + '.rsc'
        shutil.copy(inrsc,outrsc)
        shutil.copy(inrsc,filtrsc)

        # select points
        i, j, z, phi, coh, deltaz = np.loadtxt('ncycle_topo',comments='#', usecols=(0,1,2,3,5,10), unpack=True,dtype='f,f,f,f,f,f')
        z = z - float(config.z_ref)
        phi = phi*0.00020944

        topfile = path.splitext(infile)[0] + '.top'
        b1, b2, b3, b4, b5 =  np.loadtxt(topfile,usecols=(0,1,2,3,4), unpack=True, dtype='f,f,f,f,f')

        if ((config.jend_mask > config.jbeg_mask) or (config.iend_mask > config.ibeg_mask)) and config.ivar<2 :
            sys.exit(0)
            b1, b2, b3, b4, b5 = 0, 0, 0, 0, 0

            index = np.nonzero(
            np.logical_and(coh>config.thresh_amp_atmo,
            np.logical_and(deltaz>75.,
            np.logical_and(np.logical_or(i<config.ibeg_mask,pix_az>config.iend_mask),
            np.logical_or(j<config.jbeg_mask,j>config.jend_mask),
            ))))

            phi_select = phi[index]
            z_select = z[index]

            if config.nfit_atmo == -1:
                b1 = np.nanmedian(phi_select)
                fit = z_select*b1
            elif config.nfit_atmo == 0:
                b1 = np.nanmean(phi_select)
                fit = z_select*b1
            elif config.nfit_atmo == 1:
                from sklearn.linear_model import LinearRegression
                model = LinearRegression()
                model.fit(z_select, phi_select)
                fit = model.predict(z_select)
            else:
                from sklearn.preprocessing import PolynomialFeatures
                polynomial_features= PolynomialFeatures(degree=config.nfit_atmo)
                x_poly = polynomial_features.fit_transform(z_select)
                model = LinearRegression()
                model.fit(x_poly, phi_select)
                fit = model.predict(x_poly)
            
            # save median phase/topo
            strattxt = path.splitext(infile)[0] + '_strat.top'
            np.savetxt(strattxt, np.array([b1, b2, b3, b4]), header='# z   |   z**2   |   z**3   |   z**4  ' ,  fmt=('%.8f','%.8f','%.8f','%.8f'))
            
            # clean and prepare for write strat
            infileunw = path.splitext(infile)[0] + '.unw'
            remove(infileunw)

            logger.info("cpx2rmg.pl "+str(infile)+" "+str(infileunw))
            r = subprocess.call("cpx2rmg.pl "+str(infile)+" "+str(infileunw), shell=True)
            if r != 0:
                logger.critical("Failed to convert {0} to .unw".format(infile))
                sys.exit()

            inrsc = infile + '.rsc'
            outrsc = infileunw + '.rsc'
            shutil.copy(inrsc,outrsc)

            # remove strat file created by flatten_topo and write
            remove(stratfile)
            logger.info("write_strat_unw "+str(strattxt)+" "+str(infileunw)+" "+str(config.dem)+" "+str(stratfile)\
                    +" "+str(nfit_atmo)+" "+str(ivar)+" "+str(config.z_ref))
            r = subprocess.call("write_strat_unw "+str(strattxt)+" "+str(infileunw)+" "+str(config.dem)+" "+str(stratfile)\
                    +" "+str(nfit_atmo)+" "+str(ivar)+" "+str(config.z_ref)+" >> log_flattopo.txt", shell=True)
            if r != 0:
                logger.critical("Failed creating stratified file: {0}".format(stratfile))
                sys.exit()

            outrsc = stratfile + '.rsc'
            shutil.copy(inrsc,outrsc)

            # remove model created by flatten_topo and run
            remove(outfile)
            remove(filtout)
            logger.info("removeModel.pl "+str(infile)+" "+str(stratfile)+" "+str(outfile))
            r = subprocess.call("removeModel.pl "+str(infile)+" "+str(stratfile)+" "+str(outfile), shell=True)
            if r != 0:
                logger.critical("Failed removing stratified file: {0} from IFG {1}: ".format(stratfile,infile))
                sys.exit()

            corfile = config.stack.getcor(kk)
            corbase = path.splitext(corfile)[0]
            logger.info("nsb_SWfilter.pl "+str(path.splitext(outfile)[0])+" "+str(path.splitext(filtout)[0])+" "+str(corbase)\
                +" "+str(config.SWwindowsize)+" "+str(config.SWamplim)+" "+str(config.filterstyle))
            r = subprocess.call("nsb_SWfilter.pl "+str(path.splitext(outfile)[0])+" "+str(path.splitext(filtout)[0])+" "+str(corbase)\
                +" "+str(config.SWwindowsize)+" "+str(config.SWamplim)+" "+str(config.filterstyle), shell=True)
            if r != 0:
                logger.critical("Failed filtering IFG {0}: ".format(outfile))
                sys.exit()

        else:
            z_select = z; phi_select = phi
            # I dont understand ivar=0
            # ivar=2 needs to be implemented
            fit = b1*z_select + (b2/2.)*z_select**2 + (b3/3.)*z_select**2 + (b4/4.)*z_select**2

        # What about the cst??
        cst = np.nanmean(fit-phi_select*z_select)

        # plot phase/topo
        fig = plt.figure(0)
        ax = fig.add_subplot(1,1,1)
        # lets not plot dphi but phi
        ax.plot(z,phi*z,'.',alpha=.6)
        ax.plot(z_select,phi_select*z_select,'.',label='selected points',alpha=.6)
        ax.plot(z_select,fit-cst,'-r',lw=3,label='Fit: {0:.3f}z + {1:.3f}z**2 + {2:.3f}z**3 + {3:.3f}z**4'.format(b1,b2,b3,b4))
        ax.set_xlabel('Elevation (m)')
        ax.set_ylabel('Phase (rad)')
        plt.legend(loc='best')
        plotfile = path.splitext(infile)[0] + '_phase-topo.png'
        fig.savefig(plotfile, format='PNG')
        # plt.show()
        plt.close()
        del fig, ax

        # update strat
        config.strat = True

def flat_model(config,kk):
    return

def colin(config,kk):
    ''' Compute and replace amplitude by colinearity (See Pinel-Puyssegur et al., 2012)'''

    with Cd(config.stack.getpath(kk)):

        infile = config.stack.getname(kk) + '.int';  checkinfile(infile)
        inrsc = infile + '.rsc'

        # update names
        prefix, suffix = config.stack.getfix(kk)
        newprefix = 'col_'
        config.stack.updatefix(kk,newprefix,suffix)
        outfile = config.stack.getname(kk) + '.int'
        outrsc = outfile + '.rsc'

        # Retrieve length and width
        width,length =  config.stack.getsize(kk)
        if (int(width) == 0) or (int(length) == 0):
            width,length = computesize(config,infile)
            config.stack.updatesize(kk,width,length)

        if path.exists(outfile) == False:
            logger.info('Replace Amplitude by colinearity on IFG: {0}'.format(infile))
            shutil.copy(infile,'temp')
            logger.info("colin "+str(infile)+" temp "+str(outfile)+" "+str(width)+" "+str(length)+\
                " 3 0.0001 2")
            r = subprocess.call("colin "+str(infile)+" temp "+str(outfile)+" "+str(width)+" "+str(length)+\
                " 3 0.0001 2  >> log_flatenrange.txt", shell=True)
            if r != 0:
                logger.critical('Failed replacing Amplitude by colinearity on IFG: {0}'.format(infile))
                print(colin.__doc__)
                sys.exit()
            # clean
            remove('temp')

        else:
            logger.warning('{0} exists, assuming OK'.format(outfile))

def look_int(config,kk):
    ''' Look function for IFG, coherence, strat, and radar files
    Requiered parameters:  Rlooks_int, Rlooks_unw
    '''

    # need to be in the corect dir for stupid perl scripts
    with Cd(config.stack.getpath(kk)):

        infile =  config.stack.getname(kk) + '.int';  checkinfile(infile)
        corfile =  config.stack.getcor(kk);  checkinfile(corfile) 

        # look radar file if not done
        dem = config.SARMasterDir + '/'+  'radar_' + config.Rlooks_unw + 'rlks.hgt'
        if path.exists(config.dem) is False:
            look_file(config,config.dem)
            config.dem = dem

        logger.info('Look file {0} in {1} look'.format(infile,config.rlook))
        chdir(config.stack.getpath(kk))

        # update looks
        config.stack.updatelook(kk,config.Rlooks_unw)
        outfile =  config.stack.getname(kk) + '.int'
        outcor = config.stack.getcor(kk) 

        # look strat file
        stratfile = config.stack.getstratfile(kk) + '.unw'
        if (path.exists(stratfile) == False) and (config.strat == True):
            # config.strat check if flatten_topo has been done
            look_file(config,stratfile)
            config.strat = stratfile

        chdir(config.stack.getpath(kk))

        if path.exists(outfile) is False:
            logger.info("look.pl "+str(infile)+" "+str(config.rlook))
            r= subprocess.call("look.pl "+str(infile)+" "+str(config.rlook)+" >> log_look.txt" , shell=True)
            if r != 0:
                logger.critical(' Can''t look file {0} in {1} look'.format(infile,config.rlook))
                print(look_int.__doc__)
        else:
            logger.warning('{0} exists, assuming OK'.format(outfile))
        
        if path.exists(outcor) is False:
            logger.info("look.pl "+str(corfile)+" "+str(config.rlook))
            r = subprocess.call("look.pl "+str(corfile)+" "+str(config.rlook)+" >> log_look.txt", shell=True)
            if r != 0:
                logger.critical(' Can''t look file {0} in {1} look'.format(corfile,config.rlook))
                print(look_int.__doc__)
        else:
            logger.warning('{0} exists, assuming OK'.format(outfile))        
                
        # update size
        width,length = computesize(config,outfile)
        config.stack.updatesize(kk,width,length)
        # print(outfile)
        # print(config.stack.getlook(kk))

def unwrapping(config,kk):
    ''' Unwrap function from strating seedx, seedy
    if unw_method: mpd, MP.DOIN algorthim (Grandin et al., 2012). Requiered: threshold_unw, filterSW, filterROI
    if unw_method: roi, ROIPAC algorthim. Requiered threshold_unw, filterROI
    '''

    config.stack.updatelook(kk,config.Rlooks_unw)

    with Cd(config.stack.getpath(kk)):

        infile = config.stack.getname(kk)+ '.int';  checkinfile(infile)
        inrsc = infile + '.rsc'
        filtSWfile = config.stack.getfiltSW(kk)+ '.int'
        filtROIfile = config.stack.getfiltROI(kk)+ '.int'

        unwfile = config.stack.getname(kk)+ '.unw'
        unwrsc = unwfile + '.rsc'
        unwfiltSW = config.stack.getfiltSW(kk)+ '.unw'
        unwSWrsc = unwfiltSW + '.rsc'
        unwfiltROI = config.stack.getfiltROI(kk)+ '.unw'
        unwROIrsc = unwfiltROI + '.rsc'
        shutil.copy(inrsc,unwrsc)
        shutil.copy(inrsc,unwSWrsc)
        shutil.copy(inrsc,unwROIrsc)

        # Filter with colinearity
        if path.exists(filtROIfile) == False:
            filterROI(config, kk)
            checkinfile(filtROIfile)

        if path.exists(unwfiltROI) == False:

            logger.info('Unwraped IFG:{0} with strating point col:{1} line:{2} and filterd coherence threshold {3}'.\
                format(unwfile,config.seedx,config.seedy,config.threshold_unw))
            if config.unw_method == 'mpd':
                logger.info("Unwraped IFG:{0} with MP.DOIN algorthim (Grandin et al., 2012) ".format(unwfile))
                if path.exists(unwfiltSW) == False:
                    filterSW(config, kk)
                    checkinfile(unwfiltROI)

                # my_deroul_interf has ana additional input parameter for threshold on amplitude infile (normally colinearity)
                # unwrapped firt filtSWfile and then add high frequency of filtROIfile
                logger.info("my_deroul_interf_filt "+str(filtSWfile)+" cut "+str(infile)+" "+str(unwfiltROI)\
                    +" "+str(config.seedx)+" "+str(config.seedy)+" "+str(0.04)+" "+str(config.threshold_unw)+" 0")
                r = subprocess.call("my_deroul_interf_filt "+str(filtSWfile)+" cut "+str(infile)+" "+str(unwfiltROI)\
                    +" "+str(config.seedx)+" "+str(config.seedy)+" "+str(0.04)+" "+str(config.threshold_unw)+" 0  >> log_unw.txt", shell=True)
                if r != 0:
                    print(_unwrapping.__doc__)
                    logger.critical("Failed unwrapping with MP.DOIN algorthim (Grandin et al., 2012)".format(unwfile))
                    sys.exit()
                # remove('cut')

            if config.unw_method == 'roi':

                logger.info("Unwraped IFG:{0} with ROIPAC algorithm ".format(unwfile))
                mask = path.splitext(filtSWfile)[0] + '_msk'

                logger.info("make_mask.pl "+str(path.splitext(filtROIfile)[0])+" "+str(mask)+" "+str(0.02))
                r = subprocess.call("make_mask.pl "+str(path.splitext(filtROIfile)[0])+" "+str(mask)+" "+str(0.02)+"  >> log_unw.txt", shell=True)
                if r != 0:
                    print(_unwrapping.__doc__)
                    logger.critical("Failed unwrapping IFG {0} with ROIPAC algorithm ".format(unwfile))
                    sys.exit()

                logger.info("new_cut.pl "+str(path.splitext(filtROIfile)[0]))
                r = subprocess.call("new_cut.pl "+str(path.splitext(filtROIfile)[0])+"  >> log_unw.txt", shell=True)
                if r != 0:
                    print(_unwrapping.__doc__)
                    logger.critical("Failed unwrapping IFG {0} with ROIPAC algorithm ".format(unwfile))
                    sys.exit()

                logger.info("unwrap.pl "+str(path.splitext(filtROIfile)[0])+" "+str(mask)+" "+str(path.splitext(filtROIfile)[0])\
                    +" "+str(config.threshold_unw)+" "+str(config.seedx)+" "+str(config.seedy))
                r = subprocess.call("unwrap.pl "+str(path.splitext(filtROIfile)[0])+" "+str(mask)+" "+str(path.splitext(filtROIfile)[0])\
                    +" "+str(config.threshold_unw)+" "+str(config.seedx)+" "+str(config.seedy)+"  >> log_unw.txt",shell=True)
                if r != 0:
                    print(_unwrapping.__doc__)
                    logger.critical("Failed unwrapping IFG {0} with ROIPAC algorithm ".format(unwfile))
                    sys.exit()
        else:
            logger.warning('{0} exists, assuming OK'.format(unwfiltROI))

def add_atmo_back(config,kk):
    ''' Add back stratified model computed by flatten_topo'''

    with Cd(config.stack.getpath(kk)):

        # the final product is always filtROI
        unwfile = config.stack.getfiltROI(kk) + '.unw'; checkinfile(unwfile)
        stratfile = config.stack.getstratfile(kk) + '.unw'; checkinfile(stratfile)

        # update names
        prefix, suffix = config.stack.getfix(kk)
        newsuffix = suffix.replace("_flatz", "")
        config.stack.updatefix(kk,prefix,newsuffix)
        outfile = config.stack.getname(kk) + '.unw'

        if path.exists(outfile) == False:
            logger.info("add_rmg.py --infile="+str(unwfile)+" --outfile="+str(outfile)+" --add="+str(stratfile))
            r = subprocess.call("add_rmg.py --infile="+str(unwfile)+" --outfile="+str(outfile)+" --add="+str(stratfile)+\
                 " >> log_flatenrange.txt", shell=True)
            if r != 0:
                logger.critical('Failed adding back {0} on IFG: {1}'.format(stratfile,unwfile))
                logger.critical(r) 
                sys.exit()
        else:
            logger.warning('{0} exists, assuming OK'.format(outfile))

def add_ramp_back(config,kk):
    return

def add_model_back(config,kk):
    return

##################################################################################
###  READ IMPUT PARAMETERS
##################################################################################

nproc=1

# # input parameters (not in the proc file)
home='/home/cometraid14/daouts/work/tibet/qinghai/processing/Sentinel/iw1/'
seedx=336 ## iw1
seedy=1840

# home='/home/cometraid14/daouts/work/tibet/qinghai/processing/Sentinel/iw2/'
# seedx=300 ## iw2
# seedy=2384

prefix = 'col_' 
suffix = '_sd_flatz'
iend_mask=0 # mask for empirical estimations
jend_mask=0
jbeg_mask=0
jend_mask=0
model=None # model to be removed from wrapped int

# proc file parameters
IntDir=path.abspath(home)+'/'+'int/'
ListInterfero=path.abspath(home)+'/'+'interf_pair.rsc'
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
filterStrength=2.
threshold_unw=0.35
unw_method='mpd'

nfit_topo=-1
thresh_amp_topo=0.2
ivar=1
z_ref=8000.


#### TEST DIR
prefix = '' 
suffix = '_sd'
home='/home/cometraid14/daouts/work/tibet/qinghai/processing/Sentinel/iw1/'
IntDir=path.abspath(home)+'/'+'test/'
ListInterfero=path.abspath(home)+'/'+'interf_pair_test.rsc'
nproc=2

####################
# Test Process List
####################

""" Job list is: erai look_int replace_amp filterSW filterROI flat_range flat_topo flat_model colin unwrapping add_model_back add_atmo_back add_ramp_back """
print(Job.__doc__)
# do_list =  'unwrapping add_model_back'  
do_list =  'replace_amp filterSW flat_topo colin look_int' 
jobs = Job(do_list)

print('List of Post-Processing Jobs:')
for p in jobs:
    print(p)
print()

###########
#   MAIN 
###########

print(FiltFlatUnw.__doc__)
print()

postprocess = FiltFlatUnw(
        [ListInterfero,SARMasterDir,IntDir,
        Rlooks_int, Rlooks_unw, 
        nfit_range, thresh_amp_range,
        nfit_az, thresh_amp_az,
        filterstyle,SWwindowsize, SWamplim,
        filterStrength,
        nfit_topo,thresh_amp_topo,ivar,z_ref,
        seedx,seedy,threshold_unw,unw_method], 
        prefix=prefix, suffix=suffix,
        ) 

# RUN
for p in jobs:
    
    print()
    job = getattr(p,'name')
    # print ifg names at the begining of each process
    postprocess.stack.info()

    print('----------------------------------')
    print('Run {} ....'.format(job))
    
    # run process
    go(postprocess, job, nproc)

    print('----------------------------------')
    print()

print("That's all folks")

