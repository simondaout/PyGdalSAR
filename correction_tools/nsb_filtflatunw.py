#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################

############################################
# Author        : Simon DAOUT (Oxford)
############################################

"""\
nsb_filtflatunw.py
----------------------

usage:
  nsb_filtflatunw.py [-v] [-f] [--nproc=<nb_cores>] [--prefix=<value>] [--suffix=<value>] [--jobs=<job1/job2/...>] [--list_int=<path>] [--look=<value>] \
  [--model=<path>] [--cutfile=<path>] [--ibeg_mask=<value>] [--iend_mask=<value>] [--jbeg_mask=<value>] [--jend_mask=<value>]  <proc_file> 
  nsb_filtflatunw.py -h | --help

options:
  --nproc=<nb_cores>    Use <nb_cores> local cores to create delay maps [Default: 4]
  --prefix=<value>      Prefix of the IFG at the starting of the processes $prefix$date1-$date2$suffix_$rlookrlks.int [default: '']
  --suffix=<value>      Suffix of the IFG at the starting of the processes $prefix$date1-$date2$suffix_$rlookrlks.int [default: '_sd']
  --jobs<job1/job2/...> List of Jobs to be done (eg. check_look ecmwf look_int replace_amp filterSW filterROI flatr flat_atmo flat_model colin unwrapping add_model_back add_atmo_back add_era_back add_flata_back add_flatr_back refer) 
Job list is: check_look ecmwf look_int replace_amp filterSW filterROI flatr flat_atmo flat_model colin unwrapping add_model_back add_era_back add_atmo_back add_flata_back add_flatr_back
  --list_int=<path>     Overwrite liste ifg in proc file            
  --look=<value>        starting look number, default is Rlooks_int
  --model=<path>        Model to be removed from wrapped IFG [default: None]
  --cutfile=<path>      Cut file for unwrappin [default: None]
  --ibeg_mask,iend_mask Starting and Ending columns defining mask for empirical estimations [default: 0,0]
  --jbeg_mask,jend_mask Starting and Ending lines defining mask for empirical estimations [default: 0,0]
  -v                    Verbose mode. Show more information about the processing
  -f                    Force mode. Overwrite output files
  -h --help             Show this screen.
"""

from __future__ import print_function
import shutil, sys
from os import path, environ, system, chdir, remove, getcwd, listdir, symlink
import matplotlib
if environ["TERM"].startswith("screen"):
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
from pylab import *
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
from numpy.lib.stride_tricks import as_strided
import logging
import multiprocessing
from contextlib import contextmanager
from functools import wraps, partial
# from nsbas import docopt, gdal, procparser, subprocess
import subprocess, gdal, procparser
gdal.UseExceptions()
import filecmp
from operator import methodcaller
import collections
try:
    from nsbas import docopt
except:
    import docopt

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
        logger.debug('Exit {0}'.format(self.dirname))
        logger.debug('Enter {0}'.format(self.curdir))
        chdir(self.curdir)

def checkoutfile(config,file):
    do = True
    try:
        w,l = computesize(config,file)
        if w > 0 :
            do = False
            logger.warning('{0} exists, assuming OK'.format(file))
    except:
        pass
    return do

def force_link(src,dest):
    try:
        symlink(src,dest)
    except:
        rm(dest)
        symlink(src,dest)

# check if rsc file exist and good size
def copyrsc(src,dest):
    if path.exists(dest):
        if filecmp.cmp(src,dest):
            pass
        else:
            shutil.copy(src,dest)
    else:
            shutil.copy(src,dest)

def rm(file):
    if path.exists(file):
        remove(file)

def checkinfile(file):
    if path.exists(file) is False:
        logger.critical("File: {0} not found, Exit!".format(file))
        print("File: {0} not found in {1}, Exit!".format(file,getcwd()))
        # print(listdir('./'))

##################################################################################
###  RUN FONCTION
##################################################################################

# create generator for pool
@contextmanager
def poolcontext(*arg, **kargs):
    pool = multiprocessing.Pool(*arg, **kargs)
    yield pool
    pool.terminate()
    pool.join()

def go(config,job,nproc):
    ''' RUN processing function '''

    with TimeIt():

        work = range(config.Nifg)
        pool = multiprocessing.Pool(processes=nproc)

        with poolcontext(processes=nproc) as pool:
            results = pool.map(partial(eval(job), config), work)
        
    return results

def run(cmd):
    """
    Runs a shell command, and print it before running.

    Arguments:
        cmd: string to be passed to a shell

    Both stdout and stderr of the shell in which the command is run are those
    of the parent process.
    """

    logger.info(cmd)
    r = subprocess.call(cmd, shell=True, stdout=sys.stdout, stderr=subprocess.STDOUT,
        env=environ)
    if r != 0:
        logger.critical(r)
    return

##################################################################################
###  Define Job, IFG, Images and FiltFlatUnw classes  
##################################################################################

class Job():
    """ Create a class of Jobs to be run: 
    Job list is: check_look ecmwf look_int replace_amp filterSW filterROI flatr flat_atmo flat_model colin unwrapping add_model_back add_atmo_back add_era_back add_flata_back add_flatr_back """

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
        print('ecmwf look_int replace_amp filter flatr flat_atmo flat_model colin \
            unwrapping add_model_back add_atmo_back add_ramp_back add_era_back')
        print('Choose them in the order that you want')

class PileInt:
    def __init__(self,dates1, dates2, prefix, suffix, Rlooks_int, Rlooks_unw, look, filterstyle ,dir, EraDir):
        self.dates1, self.dates2 = dates1, dates2
        self.dir = dir
        self.EraDir = EraDir
        self.filterstyle = filterstyle
        self.prefix = prefix
        self.suffix = suffix
        self.look = look
        self.Rlooks_int = Rlooks_int
        self.Rlooks_unw = Rlooks_unw

        self.Nifg=len(np.atleast_1d(self.dates1))
        print("number of interferogram: ",self.Nifg)

        # init width and length to zero as long as we don't need it
        width, length = 0, 0

        try:
            logger.info('Define IFG list')
            self._ifgs = [IFG(date1,date2,look,prefix,suffix,width,length,1) for (date1,date2) in zip(self.dates1,self.dates2)]
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

    def get_default_name(self,kk):
        ''' Return interfergram file name'''
        _f = str(self._ifgs[kk].prefix) + str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + str(self._ifgs[kk].suffix) 
        temp_rlook = 1 
        print(_f+'.int')
        print()
        if path.exists(_f+'.int') == False:
            _f = str(self._ifgs[kk].prefix) + str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + str(self._ifgs[kk].suffix) + '_2rlks'
            temp_rlook = 2
        return _f, temp_rlook

    def get_default_cor(self,kk):
        ''' Return cohrence file name'''
        _f = str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + '.cor'
        temp_rlook = 1 
        if path.exists(_f) == False:
            _f = str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + '_2rlks.cor'
        return  _f

    def getname(self,kk):
        ''' Return interfergram file name '''
        if int(self._ifgs[kk].look) > 1:
          _f = str(self._ifgs[kk].prefix) + str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + str(self._ifgs[kk].suffix) + '_' + self._ifgs[kk].look + 'rlks'
        else:
          _f = str(self._ifgs[kk].prefix) + str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + str(self._ifgs[kk].suffix) 
        return _f

    def getfiltSW(self,kk):
        ''' Return interfergram file name '''
        if int(self._ifgs[kk].look) > 1:
          _f = 'filt' + str(self.filterstyle) + '_' + str(self._ifgs[kk].prefix) + str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + str(self._ifgs[kk].suffix) + '_' +  self._ifgs[kk].look + 'rlks'
        else:
          _f = 'filt' + str(self.filterstyle) + '_' + str(self._ifgs[kk].prefix) + str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + str(self._ifgs[kk].suffix) 
        return _f

    def getfiltROI(self,kk):
        ''' Return interfergram file name '''
        if int(self._ifgs[kk].look) > 1:
          _f = 'filt_' + str(self._ifgs[kk].prefix) + str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + str(self._ifgs[kk].suffix) + '_' +  self._ifgs[kk].look + 'rlks'
        else:
          _f = 'filt_' + str(self._ifgs[kk].prefix) + str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + str(self._ifgs[kk].suffix) 
        return _f

    def getcor(self,kk):
        ''' Return cohrence file name '''
        if int(self._ifgs[kk].look) > 1:
          _f = str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + '_' +  self._ifgs[kk].look + 'rlks.cor'
        else:
          _f = str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + '.cor'
        return _f

    def getsize(self,kk):
        ''' Return width and length IFG '''
        return  str(self._ifgs[kk].width),  str(self._ifgs[kk].length)

    def getpath(self,kk):
        ''' Return path ifg dir '''
        return  str(self.dir + '/'  + 'int_' + str(self.dates1[kk]) + '_' + str(self.dates2[kk])) 

    def getstratfile(self,kk):
        ''' Return stratified file name '''
        if int(self._ifgs[kk].look) > 1:
          _f = str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + '_strat_' +  self._ifgs[kk].look + 'rlks'
        else:
          _f = str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + '_strat' 
        return _f

    def getmodelfile(self,kk):
        ''' Return model file name 
        Assume model computed on the Rlooks_unw IFG...
        '''
        if int(self.Rlooks_unw) > 1:
          _f = str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) +  '_' + self.Rlooks_unw + 'rlks' + '.acp'
        else:
          _f = str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + '.acp'
        return _f

    def geterafiles(self,kk):
        ''' Return ERA model file names '''
        if int(self.Rlooks_int) > 1:
          mdel1 = str(self.EraDir) + '/'+ str(self._ifgs[kk].date1) + '_mdel' + '_' + self.Rlooks_int + 'rlks' + '.unw'
          mdel2 = str(self.EraDir) + '/'+ str(self._ifgs[kk].date2) + '_mdel' + '_' + self.Rlooks_int + 'rlks' + '.unw'
          mdel12 = str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + '_mdel' + '_' + self.Rlooks_int + 'rlks' + '.unw'
        else:
          mdel1 = str(self.EraDir) + '/'+ str(self._ifgs[kk].date1) + '_mdel' + '.unw'
          mdel2 = str(self.EraDir) + '/'+ str(self._ifgs[kk].date2) + '_mdel' + '.unw'
          mdel12 = str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + '_mdel' + '.unw'
        return  mdel1, mdel2, mdel12

    def getflatrfile(self,kk):
        ''' Return stratified file name 
        Assume estimation computed on the Rlooks_int IFG...
        '''
        if int(self.Rlooks_int) > 1:
          _f=str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + '_' +  self.Rlooks_int + 'rlks' + '.flatr'
        else:
          _f=str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + '.flatr'
        return _f

    def getflatafile(self,kk):
        ''' Return stratified file name 
        Supposed estimation computed on the Rlooks_int IFG...
        '''
        if int(self.Rlooks_int) > 1:
          _f=str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + '_' +  self.Rlooks_int + 'rlks' + '.flata'
        else:
          _f=str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + '.flata'
        return _f

    def updatelook(self,kk,newlook):
        self._ifgs[kk] = self._ifgs[kk]._replace(look=newlook)  
        self.look = newlook

    def updatesize(self,kk,newwidth,newlength):
        self._ifgs[kk] = self._ifgs[kk]._replace(width=newwidth,length=newlength)

    def updatesuccess(self,kk):
        self._ifgs[kk] = self._ifgs[kk]._replace(success=0)

    def updatefix(self,kk,newprefix, newsuffix):
        self._ifgs[kk] = self._ifgs[kk]._replace(prefix=str(newprefix), suffix=str(newsuffix))
        self.prefix, self.suffix = newprefix, newsuffix

    def info(self):
        print('List of interferograms to be process:')
        # print ([self._ifgs[kk] for kk in range(self.Nifg)])
        print ([self.getname(kk) for kk in range(self.Nifg)])
        print()

# class PileImages:
#     def __init__(self,dates1,dates2):
#         self.dates1, self.dates2 = dates1, dates2

#         # define list of images 
#         im = []; bt = []
#         for date1,date2 in zip(self.dates1,self.dates2):
#             if date1 not in im: im.append(date1)
#             if date2 not in im: im.append(date2)
#         self.Nimages=len(im)
#         print("number of image: ",self.Nimages)
#         print(im)
#         imd = date2dec(im)
#         cst = np.copy(imd[0])
#         for i in range((self.Nimages)):
#             bt.append(imd[i]-cst)
#         del cst

#         try:
#             logger.info('Define list of Images')
#             self._images = [Image(date,dec,baseline) for (date,dec,baseline) in zip(im,imd,bt)]
#         except ValueError as error:
#             logger.critical(error)
    
#     # create spetial methods len() and getititem for PileImages class
#     def __len__(self):
#         return len(self._images)

#     def __getitem__(self,kk):
#         return self._images[kk]

#     def info(self):
#         print('List of Images:')
#         print ([self._images[kk] for kk in range(self.Nimages)])
#         print()

class FiltFlatUnw:
    """ Create a class FiltFlatUnw defining all the post-procesing functions 
    list of parameters defined in the proc file: ListInterfero, SARMasterDir, IntDir, Rlooks_int, Rlooks_unw, prefix, suffix ,
    nfit_range, hresh_amp_range, nfit_az, thresh_amp_az, filterstyle, SWwindowsize, SWamplim, FilterStrength, nfit_atmo,thresh_amp_atmo, ivar, z_ref, min_z, detla_z,
    seedx, seedy.threshold_unw, unw_method, ref_top, ref_left, ref_width, ref_length
    Additional parameters not in the proc file (yet?): ibeg_mask, iend_mask, jbeg_mask, jend_mask (default: 0.)
    defining the boundary of the mask zone for emprical estimations
    suffix, preffix: define name of the interferogram at the start of the processes
    model: model to be removed from wrapped interferograms (default: None)
    cutfile: cut file fro unwrapping, decrease coh threshold for pixel with cut>0
    """

    def __init__(self, params, prefix='', suffix='_sd', look=2, ibeg_mask=0, iend_mask=0, jbeg_mask=0, jend_mask=0, model=None,force=False,cutfile=None):
        (self.ListInterfero, self.SARMasterDir, self.IntDir, self.EraDir,
        self.Rlooks_int, self.Rlooks_unw, 
        self.nfit_range, self.thresh_amp_range,
        self.nfit_az, self.thresh_amp_az,
        self.filterstyle,self.SWwindowsize, self.SWamplim,
        self.FilterStrength,self.Filt_method,
        self.nfit_atmo,self.thresh_amp_atmo,self.ivar,self.z_ref,self.min_z,self.delta_z,
        self.seedx, self.seedy,self.threshold_unw,self.threshold_unfilt,self.unw_method,self.ref_top,        self.ref_left,self.ref_width,self.ref_length
        ) = map(str, params)

        # initialise prefix, suffix
        self.prefix, self.suffix = prefix, suffix

        # initiliase number of looks
        self.look = look
        self.rlook = int(int(self.Rlooks_unw) / int(self.Rlooks_int))

        # initialise model to be removed from wrapped int
        self.model = model
        self.strat = False
        self.cutfile = cutfile

        # initilise radar file
        if int(self.Rlooks_int) > 1:
          self.dem =  self.SARMasterDir + '/'+  'radar_' + self.Rlooks_int + 'rlks.hgt'
        else:
          self.dem =  self.SARMasterDir + '/'+  'radar.hgt'

        # mask empirical estimations
        self.ibeg_mask, self.iend_mask, self.jbeg_mask, self.jend_mask = ibeg_mask, iend_mask, jbeg_mask, jend_mask

        # define list of interferograms
        dates1, dates2=np.loadtxt(self.ListInterfero,comments="#",unpack=True,usecols=(0,1),dtype='i,i') 
        self.stack = PileInt(dates1,dates2,self.prefix,self.suffix,self.Rlooks_int, self.Rlooks_unw, self.look,self.filterstyle, self.IntDir, self.EraDir)
        self.stack.info()
        self.Nifg = len(self.stack)

        # # define list images
        # self.images = PileImages(dates1,dates2)
        # self.images.info()
        # self.Nimages = len(self.images)

    def getconfig(self,kk):
         return self.stack[kk].prefix, self.stack[kk].suffix, self.stack[kk].look, self.stack[kk].success

##################################################################################
###  Define Job functions 
##################################################################################

def look_file(config,file):
    ''' Look function 
    Requiered parameters:  Rlooks_int, Rlooks_unw
    '''
    
    dirname, filename = path.split(path.abspath(file)) 
    with Cd(dirname):
        try:
            run("look.pl "+str(filename)+" "+str(config.rlook)+" > log_look.txt")
        except Exception as e:
            logger.critical(e)
            config.stack.updatesuccess(kk)
            print(look_file.__doc__)

def check_look(config,kk):
    ''' This function aims at multilooking all your IFGs and radar file to Rlooks_int. 
    To be used at a first stage of the processing if you want to corrrect your wrapped Ifgs to a higher look that the look they have been generated.
    Attention: Default hardcoding factor of Rlooks_int/2. Igs are assumed to be generated at 2looks.. need to add a new argument to change that.
    Requiered parameters:  Rlooks_int
    '''

    with Cd(config.stack.getpath(kk)):

        infile, temp_rlook = config.stack.get_default_name(kk) 
        infile= infile + '.int'; checkinfile(infile)
        corfile = config.stack.get_default_cor(kk); checkinfile(corfile)

        outfile =  config.stack.getname(kk) + '.int'
        outcor =  config.stack.getcor(kk)
        if force:
            rm(outfile); rm(outcor)

        look = int(int(config.Rlooks_int)/temp_rlook)

        do = checkoutfile(config,outfile)
        if do:
            try:
                run("look.pl "+str(infile)+" "+str(look)+" > log_look.txt")
            except Exception as e:
                logger.critical(e)
                logger.critical(' Can''t look file {0} in {1} look'.format(infile,look))
                print(check_look.__doc__)
                config.stack.updatesuccess(kk)
        
        do = checkoutfile(config,outcor)
        if do:
            try:
                run("look.pl "+str(corfile)+" "+str(look)+" >> log_look.txt")
            except Exception as e:
                logger.critical(e)
                logger.critical(' Can''t look file {0} in {1} look'.format(corfile,look))
                print(check_look.__doc__)
                config.stack.updatesuccess(kk)
                
        # update size
        width,length = computesize(config,outfile)
        config.stack.updatesize(kk,width,length)

        return config.getconfig(kk)

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

def ecmwf(config,kk):
    '''ERA Atmospheric corrections applied before filtering and unwrapping
    Requiered proc parameter: EraDir '''
    
    with Cd(config.stack.getpath(kk)):
        infile = config.stack.getname(kk)+ '.int'; checkinfile(infile)
        rscfile = infile + '.rsc'
        erafiles = config.stack.geterafiles(kk)
        checkinfile(erafiles[0]); checkinfile(erafiles[1]); 

        # update names
        prefix, suffix = config.stack.getfix(kk)
        newsuffix = suffix + '_era'
        config.stack.updatefix(kk,prefix,newsuffix)
        outfile = config.stack.getname(kk) + '.int'
        outrsc = outfile + '.rsc' 
        copyrsc(rscfile,outrsc)

        if force:
            rm(outfile), rm(erafiles[2])
        # check if not done
        do = checkoutfile(config,outfile)
        if do:
            try:
                run("add_rmg.pl "+str(erafiles[1])+" "+str(erafiles[0])+" "+str(erafiles[2])+" -1 0 > log_erai.txt" )
                checkinfile(erafiles[2])
                run("removeModel.pl "+str(infile)+" "+str(erafiles[2])+" "+str(outfile)+" >> log_erai.txt")

            except Exception as e:
                logger.critical(e)
                config.stack.updatesuccess(kk)
                print(ecmwf.__doc__)

    return config.getconfig(kk)

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
        copyrsc(rscfile,outrsc)

        if force:
            rm(outfile)
        # check if not done
        do = checkoutfile(config,outfile)
        if do:
            try:
                run("rmg2mag_phs "+str(corfile)+" tmp cor "+str(width)+" > log_replaceAMP.txt")
                run("cpx2mag_phs "+str(infile)+" tmp2 phs "+str(width)+" >> log_replaceAMP.txt")
                run("mag_phs2cpx cor phs "+str(outfile)+" "+str(width)+" >> log_replaceAMP.txt")
                rm('tmp'); rm('tmp2'); rm('phs'); rm('cor')

            except Exception as e:
                logger.critical(e)
                config.stack.updatesuccess(kk)
                print(replace_amp.__doc__)

    return config.getconfig(kk)

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

        if force:
            rm(outfile)
        do = checkoutfile(config,outfile)
        if do:
         try:
            run("nsb_SWfilter.pl "+str(inbase)+" "+str(filtbase)+" "+str(corbase)\
                    +" "+str(config.SWwindowsize)+" "+str(config.SWamplim)+" "+str(config.filterstyle)+"> log_filtSW.txt")
            if path.exists(filtrsc) == False:
                copyrsc(inrsc,filtrsc)
         except Exception as e:
            logger.critical(e)
            logger.critical('Filtering {0} with {1} filter type Failed!'.format(infile,config.filterstyle))
            config.stack.updatesuccess(kk)
            print(filterSW.__doc__)
            

    return config.getconfig(kk)

def filterROI(config, kk):
    ''' ROI-PAC Filter function
    Requiered proc file parameter: FilterStrength, Filt_method
    '''

    with Cd(config.stack.getpath(kk)):

        inbase = config.stack.getname(kk)
        infile = config.stack.getname(kk) + '.int'; checkinfile(infile)
        inrsc = infile + '.rsc'
        # get width and compute if not already done
        width,length =  config.stack.getsize(kk)
        if (int(width) == 0) or (int(length) == 0):
            width,length = computesize(config,infile)
            config.stack.updatesize(kk,width,length)

        filtbase = config.stack.getfiltROI(kk)
        filtfile = config.stack.getfiltROI(kk) + '.int'
        filtrsc = filtfile + '.rsc'
        if path.exists(filtrsc) == False:
            copyrsc(inrsc,filtrsc)

        if force:
            rm(filtfile)
        do = checkoutfile(config,filtfile)
        if do:
          try:
            #run("adapt_filt "+str(infile)+" "+str(filtfile)+" "+str(width)+" 0.25"+" "+str(config.FilterStrength)+"> log_filtROI.txt")
            run("filter.pl "+str(inbase)+" "+str(config.FilterStrength)+" "+str(config.Filt_method)+" "+str(filtbase)+"> log_filtROI.txt")
          except Exception as e:
            logger.critical(e)
            logger.critical('Failed filtering {0} with ROI-PAC filter Failed!'.format(infile))
            config.stack.updatesuccess(kk)
            print(filterROI.__doc__)
            
    return config.getconfig(kk)
        
def flatr(config,kk):
    ''' Function flatten range  on wrapped phase  (See Doin et al., 2015)
    Requiered proc file parameters: nfit_range, thresh_amp_range
    Estimation done on filterSW file
    '''

    with Cd(config.stack.getpath(kk)):
        infile = config.stack.getname(kk) + '.int'; checkinfile(infile)
        inrsc = infile + '.rsc'
        corfile = config.stack.getcor(kk); checkinfile(corfile)
        filtfile = config.stack.getfiltSW(kk) + '.int'; checkinfile(filtfile)

        if path.exists(filtfile) == False:
            logger.warning('{0} does not exist'.format(filtfile))
            # call filter function
            filterSW(config,kk)

        # parameter file
        param = config.stack.getname(kk) + '.flatr'
        newparam = config.stack.getflatrfile(kk)

        # update names
        prefix, suffix = config.stack.getfix(kk)
        newsuffix = suffix + '_flatr'
        config.stack.updatefix(kk,prefix,newsuffix)
        outfile = config.stack.getname(kk) + '.int' 
        outrsc = outfile + '.rsc'
        filtout = config.stack.getfiltSW(kk)
        copyrsc(inrsc,outrsc)

        if force:
            rm(outfile)
        do = checkoutfile(config,outfile)
        if do:
            try:
                run("flatten_range "+str(infile)+" "+str(filtfile)+" "+str(outfile)+" "+str(filtout)+\
                " "+str(config.nfit_range)+" "+str(config.thresh_amp_range)+" > log_flatr.txt")
            except Exception as e:
                logger.critical(e)
                logger.critical("Flatten range failed for IFG: {0} Failed!".format(infile))
                print(flatr.__doc__) 
                config.stack.updatesuccess(kk)
        force_link(param,newparam)

    return config.getconfig(kk)

def flata(config,kk):
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

        # parameter file
        param = config.stack.getname(kk) + '.flata'
        newparam = config.stack.getflatafile(kk)

        if path.exists(filtfile) == False:
            logger.warning('{0} does not exist'.format(filtfile))
            # call filter function
            filterSW(config,kk)

        # update names
        prefix, suffix = config.stack.getfix(kk)
        newsuffix = suffix + '_flataz'
        config.stack.updatefix(kk,prefix,newsuffix)
        outfile = config.stack.getname(kk) + '.int' 
        filtout = config.stack.getfiltSW(kk) + '.int'
        outrsc = outfile + '.rsc'
        print(inrsc,outrsc)
        copyrsc(inrsc,outrsc)

        if force:
            try:
                rm(outfile)
            except:
                pass
        print(outfile)
        do = checkoutfile(config,outfile)
        if do:
            try:
                run("flatten_az "+str(infile)+" "+str(filtfile)+" "+str(outfile)+" "+str(filtout)+\
                " "+str(config.nfit_az)+" "+str(config.thresh_amp_az)+" > log_flata.txt")
            except Exception as e:
                logger.critical(e)
                logger.critical("Flatten azimuth failed for int. {0}-{1} Failed!".format(date1,date2))
                print(flata.__doc__)
                config.stack.updatesuccess(kk)

        force_link(param,newparam)

    return config.getconfig(kk)

def flat_atmo(config, kk):
    ''' Function flatten atmosphere on wrapped phase  (See Doin et al., 2015)
    Requiered proc file parameters: nfit_atmo, ivar, z_ref, thresh_amp_atmo, min_z, delta_z
    Estimation done on filterSW file
    Plot phase/topo in *_phase-topo.png file
    '''

    with Cd(config.stack.getpath(kk)):
        infile = config.stack.getname(kk) + '.int';  checkinfile(infile)
        corfile = config.stack.getcor(kk)
        filtfile = config.stack.getfiltSW(kk) + '.int'
        stratfile = config.stack.getstratfile(kk) + '.unw' 
        topfile = path.splitext(infile)[0] + '.top'

        # filt must be done before changing name
        if path.exists(filtfile) == False:
            logger.info('{0} does not exist'.format(filtfile))
            # call filter function
            filterSW(config,kk)
            checkinfile(filtfile)

        width,length = computesize(config,filtfile)
        if (int(width) == 0) or (int(length) == 0):
            logger.info('{0} has zero size'.format(filtfile))
            filterSW(config,kk)
            checkinfile(filtfile)

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
        if int(config.Rlooks_unw) > 1:
          config.dem = config.SARMasterDir + '/'+  'radar_' + config.Rlooks_unw + 'rlks.hgt'
        else:
          config.dem = config.SARMasterDir + '/'+  'radar.hgt'

    if (np.int(config.ivar) == 2) and int(config.nfit_atmo) < 1:
        logger.warning("ivar = 2, nfit must be > 0. Set nfit to 1")
        config.nfit_atmo = 1

    with Cd(config.stack.getpath(kk)):

        if force:
            rm(outfile); rm(topfile)
        do = checkoutfile(config,outfile)
        if do:
            try:
                 run("flatten_topo "+str(infile)+" "+str(filtfile)+" "+str(config.dem)+" "+str(outfile)+" "+str(filtout)\
                +" "+str(config.nfit_atmo)+" "+str(config.ivar)+" "+str(config.z_ref)+" "+str(config.thresh_amp_atmo)+" "+\
                str(stratfile)+" > log_flatatmo.txt")
            except Exception as e:
                logger.critical(e)
                logger.critical("Flatten topo failed for int. {0} Failed!".format(infile))
                print(flat_atmo.__doc__)
                config.stack.updatesuccess(kk)
        
        inrsc = infile + '.rsc'
        outrsc = outfile + '.rsc'
        filtrsc = filtout + '.rsc'
        copyrsc(inrsc,outrsc)
        copyrsc(inrsc,filtrsc)
        stratrsc =  stratfile + '.rsc'
        copyrsc(inrsc,stratrsc)

        # select points
        i, j, z, dphi, coh, az, deltaz = np.loadtxt('ncycle_topo',comments='#', usecols=(0,1,2,3,5,9,10), unpack=True,dtype='f,f,f,f,f,f,f')
        dphi = dphi*0.00020944

        # open parameters for plot 
        b1, b2, b3, b4, b5 =  np.loadtxt(topfile,usecols=(0,1,2,3,4), unpack=True, dtype='f,f,f,f,f')

        # I dont understand ivar=0
        # ivar=2 needs to be implemented with mask
        if int(config.ivar)>1 :
            logger.warning("ivar=2, no masking implemented !!!")

        # print(int(config.jend_mask),int(config.jbeg_mask),int(config.iend_mask),int(config.ibeg_mask),np.float(config.min_z),config.delta_z)       
        if ((int(config.jend_mask) > int(config.jbeg_mask)) or (int(config.iend_mask) > int(config.ibeg_mask)) or np.float(config.min_z) > 0. or  np.float(config.delta_z) != 75.)  and int(config.ivar)<2 :

            import scipy.optimize as opt
            import scipy.linalg as lst

            b1, b2, b3, b4, b5 = 0, 0, 0, 0, 0
            
            #print(config.thresh_amp_atmo, config.min_z, config.delta_z)
            # print(config.ibeg_mask,config.iend_mask)
            # print(config.jbeg_mask,config.jend_mask)
            # print(config.thresh_amp_atmo,config.min_z,config.delta_z)

            index = np.nonzero(
            np.logical_and(coh>np.float(config.thresh_amp_atmo),
            np.logical_and(z>np.float(config.min_z),
            np.logical_and(deltaz>np.float(config.delta_z),
            np.logical_and(np.logical_or(i<int(config.ibeg_mask),i>int(config.iend_mask)),
            np.logical_or(j<int(config.jbeg_mask),j>int(config.jend_mask)),
            )))))
            # )))

            dphi_select = dphi[index]; z_select = z[index]; az_select = az[index]
            rms = 1./coh[index]

            # new top file
            strattxt = path.splitext(infile)[0] + '_strat.top'
            rm(strattxt)
            if config.nfit_atmo == str(-1):
                b1 = np.nanmedian(dphi_select)

                # save median phase/topo
                wf = open(strattxt, "w")
                wf.write("%16.8E" % (b1))
                wf.close()

            elif config.nfit_atmo == str(0):
                b1 = np.nanmean(dphi_select)

                # save mean phase/topo
                wf = open(strattxt, "w")
                wf.write("%16.8E" % (b1))
                wf.close()

            elif config.nfit_atmo == str(1):
                G=np.zeros((len(z_select),2))
                G[:,0] = z_select 
                G[:,1] = 1
                x0 = lst.lstsq(G,dphi_select)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-dphi_select)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-dphi_select)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0)[0]
                b2 = pars[0]; b1 = pars[1]
                # print(pars)

                # save mean phase/topo
                wf = open(strattxt, "w")
                wf.write("%15.8E %15.8E" % (b1, b2))
                wf.close()

            elif config.nfit_atmo == str(2):
                G=np.zeros((len(z_select),3))
                G[:,0] = z_select**2 
                G[:,1] = z_select
                G[:,2] = 1
                x0 = lst.lstsq(G,dphi_select)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-dphi_select)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-dphi_select)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0)[0]
                b3 = pars[0]; b2 = pars[1]; b1 = pars[2]
                # print(pars)

                # save mean phase/topo
                wf = open(strattxt, "w")
                wf.write("%15.8E %15.8E %15.8E" % (b1, b2, b3))
                wf.close()

            elif config.nfit_atmo == str(3):
                G=np.zeros((len(z_select),4))
                G[:,0] = z_select**3 
                G[:,1] = z_select**2 
                G[:,2] = z_select
                G[:,3] = 1
                x0 = lst.lstsq(G,dphi_select)[0]
                _func = lambda x: np.sum(((np.dot(G,x)-dphi_select)/rms)**2)
                _fprime = lambda x: 2*np.dot(G.T/rms, (np.dot(G,x)-dphi_select)/rms)
                pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0)[0]
                b4 = pars[0]; b3 = pars[1]; b2 = pars[2]; b1 = pars[3]

                # save mean phase/topo
                wf = open(strattxt, "w")
                wf.write("%15.8E %15.8E %15.8E %15.8E" % (b1, b2, b3, b4))
                wf.close()

            # clean and prepare for write strat
            infileunw = path.splitext(infile)[0] + '.unw'
            rm(infileunw)

            try: 
                run("cpx2rmg.pl "+str(infile)+" "+str(infileunw)+" >> log_flatatmo.txt")
            except Exception as e:
                logger.critical(e)
                logger.critical("Failed to convert {0} to .unw".format(infile))
                print(flat_atmo.__doc__)
                config.stack.updatesuccess(kk)

            inrsc = infile + '.rsc'
            outrsc = infileunw + '.rsc'
            copyrsc(inrsc,outrsc)

            # remove strat file created by flatten_topo and write the new one
            rm(stratfile)
            try:
                run("write_strat_unw "+str(strattxt)+" "+str(infileunw)+" "+str(config.dem)+" "+str(stratfile)\
                    +" "+str(config.nfit_atmo)+" "+str(config.ivar)+" "+str(config.z_ref)+" >> log_flatatmo.txt")
            except:
                logger.critical("Failed creating stratified file: {0}".format(stratfile))
                print(flat_atmo.__doc__)
                config.stack.updatesuccess(kk)
                sys.exit()

            outrsc = stratfile + '.rsc'
            copyrsc(inrsc,outrsc)

            # remove corrected ifg created by flatten_topo and create new one
            rm(outfile)
            rm(filtout)
            try:
                run("removeModel.pl "+str(infile)+" "+str(stratfile)+" "+str(outfile)+" >> log_flatatmo.txt")
            except Exception as e:
                logger.critical(e)
                logger.critical("Failed removing stratified file: {0} from IFG {1}: ".format(stratfile,infile))
                print(flat_atmo.__doc__)
                config.stack.updatesuccess(kk)

            corfile = config.stack.getcor(kk)
            corbase = path.splitext(corfile)[0]
            try:
                run("nsb_SWfilter.pl "+str(path.splitext(outfile)[0])+" "+str(path.splitext(filtout)[0])+" "+str(corbase)\
                +" "+str(config.SWwindowsize)+" "+str(config.SWamplim)+" "+str(config.filterstyle)+" >> log_flatatmo.txt")
            except Exception as e:
                logger.critical(e)
                logger.critical("Failed filtering IFG {0}: ".format(outfile))
                config.stack.updatesuccess(kk)
                print(flat_atmo.__doc__)

            # rm(infileunw)
                
        else:
            z_select = np.copy(z); dphi_select = np.copy(dphi); az_select = np.copy(az)

        # lets not plot dphi but phi
        phi = dphi*z
        
        index = z_select.argsort()
        z_select = z_select[index]; dphi_select = dphi_select[index]; az_select = az_select[index]
        phi_select = dphi_select*z_select

        if int(config.ivar)<2:
            fit = b1*(z_select - np.float(config.z_ref)) + (b2/2.)*((z_select)**2 -np.float(config.z_ref)**2) + \
            (b3/3.)*((z_select)**3-np.float(config.z_ref)**3) + (b4/4.)*((z_select)**4 - np.float(config.z_ref)**4)
        elif (config.nfit_atmo) == 3 :
            fit = (b1+b2*az_select)/2.*(z_select-np.float(config.z_ref))**2 + (b3+b4*az_select)/3.*(z_select-np.float(config.z_ref))**3
        else:
            fit = b1*(z_select-np.float(config.z_ref)) +(b2+b3*az_select)/2.*(z_select-np.float(config.z_ref))**2 + (b4+b5*az_select)/3.*(z_select-np.float(config.z_ref))**3

        # compute crappy constant for the plot
        cst = np.nanmedian(fit-phi_select)

        # plot phase/topo
        fig = plt.figure(0, figsize=(11,4))
        ax1 = fig.add_subplot(1,2,1)
        ax1.plot(z,dphi*z,'.',label='all points',alpha=.2)
        ax1.set_xlabel('Elevation (m)')
        ax1.set_ylabel('Phase (rad)')
        ymin,ymax = ax1.get_ylim() 
        # print(ymin,ymax)

        ax = fig.add_subplot(1,2,2)
        ax.plot(z,dphi*z,'.',label='all points',alpha=.2)
        ax.plot(z_select,phi_select,'.',label='selected points',alpha=.4)
        if int(config.ivar)<2:
            ax.plot(z_select,fit-cst,'-r',lw=4,label='Fit: {0:.2E}z + {1:.2E}z**2 + {2:.2E}z**3 + {3:.2E}z**4'.format(b1,b2,b3,b4))
        else:
            ax.plot(z_select,fit-cst,'-r',lw=4,label='Fit: {0:.2E}*z + ({1:.2E}+{2:.2E}az)/2.*z**2 + ({3:.2E}+{4:.2E}*az)/3.*z**3'.format(b1,b2,b3,b4,b5))

        ax.set_xlabel('Elevation (m)')
        ax.set_ylabel('Phase (rad)')
        ax.set_ylim(ymin,ymax)
        plt.legend(loc='best',fontsize = 'x-small')
        plotfile = path.splitext(infile)[0] + '_phase-topo.png'
        fig.savefig(plotfile, format='PNG')
       

        # plt.show()
        plt.close()
        del fig, ax

        return config.getconfig(kk)

def flat_model(config,kk):
    ''' Function flatten model on wrapped phase  (See Daout et al., 2017)
        Requiered model file
        Requiered proc file parameters: thresh_amp_atmo
        Estimation done on filterSW file
    '''

    with Cd(config.stack.getpath(kk)):

        # lets asume we do this in look unw
        config.stack.updatelook(kk,config.Rlooks_unw)

        infile = config.stack.getname(kk) + '.int'; checkinfile(infile)
        inrsc = infile + '.rsc'

        filtfile = config.stack.getfiltSW(kk) + '.int'
        # filt must be done before changing name
        if path.exists(filtfile) == False:
            logger.info('{0} does not exist'.format(filtfile))
            # call filter function
            filterSW(config,kk)
            checkinfile(filtfile)

        width,length = computesize(config,filtfile)
        if (int(width) == 0) or (int(length) == 0):
            logger.info('{0} has zero size'.format(filtfile))
            filterSW(config,kk)
            checkinfile(filtfile)

        # param file
        param = config.stack.getname(kk) + '.acp'
        newparam = config.stack.getmodelfile(kk)

        if path.exists(filtfile) == False:
            logger.warning('{0} does not exist'.format(filtfile))
            # call filter function
            filterSW(config, kk)

        # update names
        prefix, suffix = config.stack.getfix(kk)
        newsuffix = suffix + '_nomodel'
        config.stack.updatefix(kk,prefix,newsuffix)
        outfile = config.stack.getname(kk) + '.int' 
        filtout = config.stack.getfiltSW(kk) + '.int'
        outrsc = outfile + '.rsc'
        filtoutrsc = filtout + '.rsc'
        copyrsc(inrsc,outrsc)
        copyrsc(inrsc,filtoutrsc)

        if force:
            rm(outfile); rm(param); rm(newparam)
        if config.model != None:
            do = checkoutfile(config,outfile)
            if do:
                try:
                    run("flatten_stack "+str(infile)+" "+str(filtfile)+" "+str(config.model)+" "+str(outfile)+" "+str(filtout)\
                    +" "+str(config.thresh_amp_atmo)+" > log_flatmodel.txt")
                except Exception as e:
                    logger.critical(e)
                    logger.critical("Flatten model failed for int. {0} Failed!".format(infile))
                    print(flat_model.__doc__)
                    config.stack.updatesuccess(kk)

        else:
            logger.critical('Model file is not defined. Exit!')
            sys.exit()

        # move param file into a file name independent of prefix and suffix
        force_link(param,newparam)

    return config.getconfig(kk)

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
        copyrsc(inrsc,outrsc)

        # Retrieve length and width
        width,length =  config.stack.getsize(kk)
        if (int(width) == 0) or (int(length) == 0):
            width,length = computesize(config,infile)
            config.stack.updatesize(kk,width,length)


        if force:
            rm(outfile)
        do = checkoutfile(config,outfile)
        if do:
            try:
                copyrsc(infile,'temp')
                #run("colin "+str(infile)+" temp "+str(outfile)+" "+str(width)+" "+str(length)+\
                #" 3 0.0001 2 > log_colin.txt")
                run("colin "+str(infile)+" temp "+str(outfile)+" "+str(width)+" "+str(length)+\
                " 3 0.1 0 > log_colin.txt")
                rm('temp')
            except Exception as e:
                logger.critical(e)
                logger.critical('Failed replacing Amplitude by colinearity on IFG: {0}'.format(infile))
                print(colin.__doc__)
                config.stack.updatesuccess(kk)

    return config.getconfig(kk)

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

        # look strat file
        stratfile = config.stack.getstratfile(kk) + '.unw'

        # update looks
        config.stack.updatelook(kk,config.Rlooks_unw)
        outfile =  config.stack.getname(kk) + '.int'
        outcor = config.stack.getcor(kk) 

        chdir(config.stack.getpath(kk))

        if force:
            rm(outfile)
        do = checkoutfile(config,outfile)
        if do:
            try:
                run("look.pl "+str(infile)+" "+str(config.rlook)+" > log_look.txt")
            except Exception as e:
                logger.critical(e)
                logger.critical(' Can''t look file {0} in {1} look'.format(infile,config.rlook))
                print(look_int.__doc__)
                config.stack.updatesuccess(kk)
        
        do = checkoutfile(config,outcor)
        if do:
            try:
                run("look.pl "+str(corfile)+" "+str(config.rlook)+" >> log_look.txt")
            except Exception as e:
                logger.critical(e)
                logger.critical(' Can''t look file {0} in {1} look'.format(corfile,config.rlook))
                print(look_int.__doc__)
                config.stack.updatesuccess(kk)
                
        # update size
        width,length = computesize(config,outfile)
        config.stack.updatesize(kk,width,length)

        return config.getconfig(kk)

def unwrapping(config,kk):
    ''' Unwrap function from strating seedx, seedy
    if unw_method: mpd, MP.DOIN algorthim (Grandin et al., 2012). Requiered: threshold_unw, filterSW, filterROI, threshold_unfilt
    Unwrapped firt filtSWfile and then add high frequency of filtROIfile
    if unw_method: roi, ROIPAC algorthim. Requiered threshold_unw, filterROI
    '''

    config.stack.updatelook(kk,config.Rlooks_unw)

    with Cd(config.stack.getpath(kk)):

        infile = config.stack.getname(kk)+ '.int';  checkinfile(infile)
        inrsc = infile + '.rsc'
        corfile = config.stack.getcor(kk); checkinfile(corfile)
        filtSWfile = config.stack.getfiltSW(kk)+ '.int'
        filtROIfile = config.stack.getfiltROI(kk)+ '.int'

        unwfile = config.stack.getname(kk)+ '.unw'
        unwrsc = unwfile + '.rsc'
        unwfiltSW = config.stack.getfiltSW(kk)+ '.unw'
        unwSWrsc = unwfiltSW + '.rsc'
        unwfiltROI = config.stack.getfiltROI(kk)+ '.unw'
        unwROIrsc = unwfiltROI + '.rsc'

        bridgefile = 'bridge.in'

        if force: 
            rm(filtROIfile); rm(filtSWfile); rm(unwfiltROI); rm(unwSWrsc)

        # Filter with colinearity
        if path.exists(filtROIfile) == False:
            filterROI(config, kk)
            checkinfile(filtROIfile)
            logger.info("length.pl "+str(filtROIfile))
            r = subprocess.call("length.pl "+str(filtROIfile), shell=True)

        do = checkoutfile(config,unwfiltROI)
        if do:
          try:
            logger.info('Unwraped IFG:{0} with strating point col:{1} line:{2} and filterd coherence threshold {3}'.\
                format(unwfile,config.seedx,config.seedy,config.threshold_unw))
            
            if config.unw_method == 'mpd':
                logger.info("Unwraped IFG:{0} with MP.DOIN algorthim (Grandin et al., 2012) ".format(unwfile))
                if path.exists(filtSWfile) == False:
                    filterSW(config, kk)
                    checkinfile(filtSWfile)
                    logger.info("length.pl "+str(filtSWfile))
                    run("length.pl "+str(filtSWfile))

                if path.exists(bridgefile) == False:
                    wf = open(bridgefile,"w")
                    wf.write("1  1  1  1  0 ")
                    wf.close()

                if config.cutfile==None:
                    opt=0
                else:
                    opt=1

                # my_deroul_interf has an additional input parameter for threshold on amplitude infile (normally colinearity)
                run("my_deroul_interf_filt "+str(filtSWfile)+" "+str(config.cutfile)+" "+str(infile)+" "+str(filtROIfile)\
                    +" "+str(config.seedx)+" "+str(config.seedy)+" "+str(config.threshold_unfilt)+" "+str(config.threshold_unw)+" "+str(opt)+" > log_unw.txt")

            if config.unw_method == 'roi':

                logger.info("Unwraped IFG:{0} with ROIPAC algorithm ".format(unwfile))
                mask = path.splitext(filtROIfile)[0] + '_msk'
                cut = path.splitext(filtROIfile)[0] + '_cut.flg'

                if force:
                    rm(mask); rm(cut)

                run("make_mask.pl "+str(path.splitext(filtROIfile)[0])+" "+str(mask)+" "+str(config.threshold_unw)+" > log_unw_roi.txt")
                run("new_cut.pl "+str(path.splitext(filtROIfile)[0])+" >> log_unw_roi.txt")
                run("unwrap.pl "+str(path.splitext(filtROIfile)[0])+" "+str(mask)+" "+str(path.splitext(filtROIfile)[0])\
                    +" "+str(config.threshold_unw)+" "+str(config.seedx)+" "+str(config.seedy)+" >> log_unw_roi.txt")

            if config.unw_method == 'icu':

                logger.info("Unwraped IFG:{0} with ICU algorithm ".format(unwfile))
                mask = path.splitext(filtROIfile)[0] + '_msk'
                if force:
                    rm(mask)

                run("make_mask.pl "+str(path.splitext(filtROIfile)[0])+" "+str(mask)+" "+str(0.02)+" > log_unw_icu.txt")

                run("icu.pl "+str(path.splitext(infile)[0])+" "+str(path.splitext(filtROIfile)[0])\
                    +" "+str(path.splitext(corfile)[0])+" "+str(mask)+" "+str(path.splitext(filtROIfile)[0])+" "+str(config.Filt_method)+" "+str(config.FilterStrength)+" "+"Phase_Sigma"+" 5 "+str(config.threshold_unw)+" >> log_unw_icu.txt")

            if config.unw_method == 'snaphu':

                logger.info("Unwraped IFG:{0} with SNAPHU algorithm ".format(unwfile))
                mask = path.splitext(filtROIfile)[0] + '_msk'
                if force:
                    rm(mask)

                run("snaphu_mcf.pl smooth"+str(path.splitext(infile)[0])+" "+str(path.splitext(infile)[0])+" "+str(mask)+" "+str(path.splitext(corfile)[0])+" "+str(config.threshold_unw)+ " > log_unw_snaphu.txt")
                snaphu_file = "lp_" + path.splitext(filtSWfile)[0] 
                run("lowpass.pl "+str(path.splitext(filtROIfile)[0])+" "+str(path.splitext(corfile)[0])+" "+str(snaphu_file)+" 3 Phase_Sigma > log_unw_icu.txt")
                snaphu_file = "lp_" + path.splitext(filtSWfile)[0] + "_phsig"
                run("snaphu_mcf.pl smooth "+str(path.splitext(filtROIfile)[0])+" "+str(path.splitext(unwfile)[0])+" "+str(mask)+" "+str(snaphu_file)+" "+str(config.threshold_unw)+" > log_unw_snaphu.txt")

          except Exception as e:
            logger.critical(e)
            print(unwrapping.__doc__)
            logger.critical("Failed unwrapping IFG {0} ".format(unwfile))
            config.stack.updatesuccess(kk)
        
        copyrsc(inrsc,unwrsc)
        copyrsc(inrsc,unwSWrsc)
        copyrsc(inrsc,unwROIrsc)

    return config.getconfig(kk)

def add_era_back(config,kk):
    '''Add back ERA Atmospheric corrections'''
    
    with Cd(config.stack.getpath(kk)):
        # look era file
        erafile = str(config.stack.geterafiles(kk)[2])
        look_file(config,str(erafile))
        
        # update look unw in case not done already
        config.stack.updatelook(kk,config.Rlooks_unw)

        # the final product is always filtROI
        unwfile = config.stack.getfiltROI(kk) + '.unw'; checkinfile(unwfile)
        if "_era" in unwfile:
            unwrsc = unwfile + '.rsc'

            # update names
            prefix, suffix = config.stack.getfix(kk)
            newsuffix = suffix.replace("_era", "")
            config.stack.updatefix(kk,prefix,newsuffix)
            outfile = config.stack.getfiltROI(kk) + '.unw'
            outrsc = outfile + '.rsc'
            copyrsc(unwrsc,outrsc)

            if force:
                rm(outfile)
            do = checkoutfile(config,outfile)
            if do:
                try:
                    run("length.pl "+str(unwfile))
                    run("add_rmg.py --infile="+str(unwfile)+" --outfile="+str(outfile)+" --add="+str(erafile)+" >> log_flatatmo.txt")
                except Exception as e:
                    logger.critical(e)
                    logger.critical('Failed adding back {0} on IFG: {1}'.format(erafile,unwfile))
                    config.stack.updatesuccess(kk)
                    print(add_atmo_back.__doc__)

        else:
            logger.critical('ecmwf() does not seem to have been done... Exit!')
            config.stack.updatesuccess(kk)
            print(add_era_back.__doc__)
 
    return config.getconfig(kk)

def add_atmo_back(config,kk):
    ''' Add back stratified model computed by flatten_topo'''

    with Cd(config.stack.getpath(kk)):

        # look strat file
        if int(config.Rlooks_int) > 1:
          stratfile = str(config.stack[kk].date1) + '-' + str(config.stack[kk].date2) + '_strat_' + config.Rlooks_int + 'rlks.unw'
        else:
          stratfile = str(config.stack[kk].date1) + '-' + str(config.stack[kk].date2) + '_strat.unw'
        look_file(config,stratfile)
        
        # update look unw in case not done already
        config.stack.updatelook(kk,config.Rlooks_unw)

        # the final product is always filtROI
        unwfile = config.stack.getfiltROI(kk) + '.unw'; checkinfile(unwfile)
        if "_flatz" in unwfile:
            stratfile = config.stack.getstratfile(kk) + '.unw'; checkinfile(stratfile)
            unwrsc = unwfile + '.rsc'

            # update names
            prefix, suffix = config.stack.getfix(kk)
            newsuffix = suffix.replace("_flatz", "")
            config.stack.updatefix(kk,prefix,newsuffix)
            outfile = config.stack.getfiltROI(kk) + '.unw'
            outrsc = outfile + '.rsc'
            copyrsc(unwrsc,outrsc)

            if force:
                rm(outfile)
            do = checkoutfile(config,outfile)
            if do:
                try:
                    run("length.pl "+str(unwfile))
                    run("add_rmg.py --infile="+str(unwfile)+" --outfile="+str(outfile)+" --add="+str(stratfile)+" >> log_flatatmo.txt")
                except Exception as e:
                    logger.critical(e)
                    logger.critical('Failed adding back {0} on IFG: {1}'.format(stratfile,unwfile))
                    config.stack.updatesuccess(kk)
                    print(add_atmo_back.__doc__)

        else:
            logger.critical('flat_atmo() does not seem to have been done... Exit!')
            config.stack.updatesuccess(kk)
            print(add_atmo_back.__doc__)
 

    return config.getconfig(kk)

def add_flatr_back(config,kk):
    ''' Add range ramp estimated on wrapped IFG back on unwrapped IFG
    Requiered .flatr parameter file containing polynomial fit
    !!! Supposed flattening estimated on Rlooks_int file. 
    '''

    with Cd(config.stack.getpath(kk)):
   
        # update look unw in case not done already
        config.stack.updatelook(kk,config.Rlooks_unw)

        # the final product is always filtROI
        unwfile = config.stack.getfiltROI(kk) + '.unw'; checkinfile(unwfile)
        if "_flatr" in unwfile:
            param = config.stack.getflatrfile(kk)
            unwrsc = unwfile + '.rsc'

            # update names
            prefix, suffix = config.stack.getfix(kk)
            newsuffix = suffix.replace("_flatr", "")
            config.stack.updatefix(kk,prefix,newsuffix)
            outfile = config.stack.getfiltROI(kk) + '.unw'
            outrsc = outfile + '.rsc'
            copyrsc(unwrsc,outrsc)

            if path.exists(param):
                if force:
                    rm(outfile)
                do = checkoutfile(config,outfile)
                if do:
                    try:
                        run("length.pl "+str(unwfile))
                        run("correct_rgaz_unw.py --infile="+str(unwfile)+" --outfile="+str(outfile)+" --param="+str(param)+" --rlook_factor="+str(config.rlook)+" >> log_flatr.txt")
                    except Exception as e:
                        logger.critical(e)
                        logger.critical('Failed adding back {0} on IFG: {1}'.format(param,unwfile))
                        config.stack.updatesuccess(kk)
                        print(add_flatr_back.__doc__)

            else:
                logger.critical('Param file {0} does not exist. Exit!'.format(param))
                config.stack.updatesuccess(kk)
                print(add_flatr_back.__doc__)

        else:
            logger.critical('flatr() does not seem to have been done... Exit!')
            config.stack.updatesuccess(kk)
            print(add_flatr_back.__doc__)

    return config.getconfig(kk)

def add_flata_back(config,kk):
    ''' Add azimutal ramp estimated on wrapped IFG back on unwrapped IFG
    Requiered .flata parameter file containing polynomial fit
    !!! Supposed flattening estimated on Rlooks_int file. 
    '''

    with Cd(config.stack.getpath(kk)):
   
        # update look unw in case not done already
        config.stack.updatelook(kk,config.Rlooks_unw)

        # the final product is always filtROI
        unwfile = config.stack.getfiltROI(kk) + '.unw'; checkinfile(unwfile)
        if "_flata" in unwfile:

            param = config.stack.getflatafile(kk)
            unwrsc = unwfile + '.rsc'

            # assume flatr done on Rlooks_int...but it is always the case?
            # look_factor = int(config.Rlooks_unw) / int(config.Rlooks_int)

            # update names
            prefix, suffix = config.stack.getfix(kk)
            newsuffix = suffix.replace("_flata", "")
            config.stack.updatefix(kk,prefix,newsuffix)
            outfile = config.stack.getfiltROI(kk) + '.unw'
            outrsc = outfile + '.rsc'
            copyrsc(unwrsc,outrsc)

            if force:
                rm(outfile)
            do = checkoutfile(config,outfile)
            if do:
                if path.exists(outfile) == False:
                    try:
                        run("length.pl "+str(unwfile))
                        run("correct_rgaz_unw.py --infile="+str(unwfile)+" --outfile="+str(outfile)+" --param="+str(param)+" --rlook_factor="+str(config.rlook)+" >> log_flata.txt")
                    except Exception as e:
                        logger.critical(e)
                        logger.critical('Failed adding back {0} on IFG: {1}'.format(param,unwfile))
                        config.stack.updatesuccess(kk)
                        print(add_flata_back.__doc__)

            else:
                logger.critical('Param file {0} does not exist. Exit!'.format(param))
                config.stack.updatesuccess(kk)
                print(add_flata_back.__doc__)  
        else:
            logger.critical('flata() does not seem to have been done... Exit!')
            config.stack.updatesuccess(kk)
            print(add_flata_back.__doc__) 
                    

    return config.getconfig(kk)

def refer(config,kk):
    "  Refer unwrapped interferograms to a rectangular area  "
    " Define ref_left, ref_top, ref_width, ref_length in the proc file !"

    with Cd(config.stack.getpath(kk)):
        config.stack.updatelook(kk,config.Rlooks_unw)
        unwfile = config.stack.getfiltROI(kk) + '.unw'; checkinfile(unwfile)
        unwrsc = unwfile + '.rsc'
        
        # update names
        prefix, suffix = config.stack.getfix(kk)
        newsuffix = suffix + '_refer'
        config.stack.updatefix(kk,prefix,newsuffix)
        outfile = config.stack.getfiltROI(kk) + '.unw'
        outrsc = outfile + '.rsc'
        copyrsc(unwrsc,outrsc)
        
        if force:
            rm(outfile)
        
        try:
          run("length.pl "+str(unwfile))
          run("refer_interf "+str(unwfile)+" "+str(outfile)+" "+str(config.ref_left)+" "+str(config.ref_top)+" "+str(config.ref_width)+" "+str(config.ref_length)+" >> log_refer.txt")
          run("length.pl "+str(outfile))
        except Exception as e:
          logger.critical(e)
          logger.critical("Refer failed for int. {0} Failed!".format(unwfile))
          print(refer.__doc__)
          config.stack.updatesuccess(kk)

    return config.getconfig(kk)

def add_model_back(config,kk):
    ''' Function adding model on unwrapped IFG previously removed on wrapped IFG  (See Daout et al., 2017)
        Requiered model file and parameter file .stack
    '''

    with Cd(config.stack.getpath(kk)):

        # update look unw in case not done already
        config.stack.updatelook(kk,config.Rlooks_unw)

        # the final product is always filtROI
        unwfile = config.stack.getfiltROI(kk) + '.unw'; checkinfile(unwfile)
        if "_nomodel" in unwfile:
            unwrsc = unwfile + '.rsc'
            param = config.stack.getmodelfile(kk)

            # update names
            prefix, suffix = config.stack.getfix(kk)
            newsuffix = suffix.replace("_nomodel", "")
            config.stack.updatefix(kk,prefix,newsuffix)
            outfile = config.stack.getfiltROI(kk) + '.unw'
            outrsc = outfile + '.rsc'
            copyrsc(unwrsc,outrsc)

            if force:
                rm(outfile)
            if config.model != None:
                if path.exists(param) is True:
                    do = checkoutfile(config,outfile)
                    if do:
                        try:
                            run("length.pl "+str(unwfile))
                            run("unflatten_stack "+str(unwfile)+" "+str(outfile)+" "+str(config.model)+" "+str(param)+" >> log_flatmodel.txt")
                        except Exception as e:
                            logger.critical(e)
                            logger.critical("Unflatten model failed for int. {0} Failed!".format(unwfile))
                            print(add_model_back.__doc__)
                            config.stack.updatesuccess(kk)
                else:
                        logger.critical("Unflatten model for int. {0} Failed!".format(unwfile))
                        logger.critical("Param file {0}, does not exist!".format(param))
                        print(add_model_back.__doc__)        
                        config.stack.updatesuccess(kk)
            else:
                logger.critical("Model file not found. Exit!".format(unwfile))
                print(add_model_back.__doc__)  
                config.stack.updatesuccess(kk)
                print(add_model_back.__doc__)

        else:
            logger.critical('flat_model() does not seem to have been done... Exit!')
            config.stack.updatesuccess(kk)
            print(add_model_back.__doc__)         

    return config.getconfig(kk)

##################################################################################
###  READ IMPUT PARAMETERS
##################################################################################

# Parse arguments
arguments = docopt.docopt(__doc__)
home = getcwd() + '/'


if arguments["--nproc"] == None:
    nproc = 4
elif arguments["--nproc"] == 'all':
    nproc = multiprocessing.cpu_count() - 4
    print('Using all {} cpu available with all option'.format(multiprocessing.cpu_count()))
    time.sleep(1)
else:
    nproc = int(arguments["--nproc"])

if arguments["--prefix"] == None:
    prefix = ''
else:
    prefix=arguments["--prefix"]
if arguments["--suffix"] == None:
    suffix = '_sd'
else:
    suffix=arguments["--suffix"]
if arguments["--jobs"] ==  None:
    print('--jobs list is empty. Nothing to be done. Exit!')
    print(Job.__doc__)
    sys.exit()
else:
    split = methodcaller('replace','/',' ')
    do_list = split(arguments["--jobs"])

if arguments["--model"] == None:
    model = None
else:
    model = path.abspath(home)+'/'+arguments["--model"]

if arguments["--cutfile"] == None:
    cutfile = None
else:
    cutfile = path.abspath(home)+'/'+arguments["--cutfile"]

if arguments["--ibeg_mask"] == None:
    ibeg_mask = 0
else:
    ibeg_mask = int(arguments["--ibeg_mask"])

if arguments["--iend_mask"] == None:
    iend_mask = 0
else:
    iend_mask = int(arguments["--iend_mask"])

if arguments["--jbeg_mask"] == None:
    jbeg_mask = 0
else:
    jbeg_mask = int(arguments["--jbeg_mask"])

if arguments["--jend_mask"] == None:
    jend_mask = 0
else:
    jend_mask = int(arguments["--jend_mask"])

if arguments["-f"]:
    force = True
else:
    force = False

##################################################################################
###  READ PROC FILE
##################################################################################

# Read proc file
# default procfile: I dont know if that is a good idea as we maybe want error messages
# and process to stop when procfile parameters not well defined for run process
# but like this it run fine all the time
proc_defaults = {
    "IntDir": "int",
    "EraDir": "ERA",
    "ListInterfero": "interf_pair.rsc",
    "Rlooks_int": "2",
    "Rlooks_unw": "4",
    "nfit_range": "-1", # median
    "thresh_amp_range": "0.3", # threshold on coherence
    "nfit_az": "-1", # median
    "thresh_amp_az": "0.3", # threshold on coherence
    "nfit_topo": "-1", # median 
    "thresh_amp_topo": "0.2",
    "ivar": "1", # fct of topography only
    "z_ref": "8000.", # reference
    "min_z": "0.", # min elevation
    "delta_z": "75", 
    "filterstyle": "SWc",
    "SWamplim": "0.05",
    "SWwindowsize": "8",
    "FilterStrength": "0.75", # strenght filter roi
    "Filt_method": "adapt_filt", # could be: adapt_filt
    "unw_method": "roi",
    "threshold_unw": "0.35", # threshold on filtered colinearity 
    "threshold_unfilt": "0.02", # threshold on colinearity 
    "seedx": "50", # starting col for unw
    "seedy": "50", # starting line for unw
    "ref_left": "0", # left corner for refer
    "ref_top": "0", # top corner for  refer
    "ref_width": "200", # box width for refer
    "ref_length": "200", # box length for refer
    }

proc = procparser.ProcParser(proc_defaults)
if arguments["-v"]:
    print("reading proc file "+arguments["<proc_file>"])
if path.exists(arguments["<proc_file>"]):
    proc.read(arguments["<proc_file>"])
else:
    print('{0} does not exist. Exit!'.format(arguments["<proc_file>"]))
    sys.exit()

# is that usefull? not sure...
IntDir = path.abspath(home)+'/'+proc.get("IntDir")
if arguments["--list_int"] == None:
    ListInterfero = path.abspath(home)+'/'+proc["ListInterfero"]
else:
    ListInterfero = path.abspath(home)+'/'+arguments["--list_int"]
SARMasterDir = path.abspath(home)+'/'+proc["SarMasterDir"]
EraDir = path.abspath(home)+'/'+proc["EraDir"]

print('Proc File parameters:')
print('ListInterfero: {0}\n SARMasterDir: {1}\n IntDir: {2}\n  EraDir: {3}\n\
    Rlooks_int: {4}, Rlooks_unw: {5}\n\
    nfit_range: {6}, thresh_amp_range: {7}\n\
    nfit_az: {8}, thresh_amp_az: {9}\n\
    filterstyle : {10}, SWwindowsize: {11}, SWamplim: {12}\n\
    FilterStrength : {13}, Filt_method : {14}\n\
    nfit_topo: {15}, thresh_amp_topo: {16}, ivar: {17}, z_ref: {18}, min_z: {19}, delta_z: {20}\n\
    seedx: {21}, seedy: {22}, threshold_unw: {23}, threshold_unfilt: {24}, unw_method: {25}, ref_top: {26}, ref_left: {27}, ref_width: {28}, ref_length: {29}'.format(ListInterfero,SARMasterDir,IntDir,EraDir,\
    proc["Rlooks_int"], proc["Rlooks_unw"],\
    proc["nfit_range"], proc["thresh_amp_range"],\
    proc["nfit_az"], proc["thresh_amp_az"],\
    proc["filterstyle"], proc["SWwindowsize"], proc["SWamplim"],\
    proc["FilterStrength"],proc["Filt_method"],\
    proc["nfit_topo"], proc["thresh_amp_topo"], proc["ivar"], proc["z_ref"], proc["min_z"], proc["delta_z"],\
    proc["seedx"], proc["seedy"], proc["threshold_unw"],proc["threshold_unfilt"], proc["unw_method"],\
    proc["ref_top"], proc["ref_left"],proc["ref_width"],proc["ref_length"]
    ))
print()

if arguments["-v"]:
    # give me the time to read terminal
    time.sleep(3)


###########
#   MAIN 
###########

# init collections 
Process = collections.namedtuple('Process', 'name')
IFG = collections.namedtuple('IFG', 'date1 date2 look prefix suffix width length success')
Image = collections.namedtuple('Image', 'date decimal_date temporal_baseline')

# init logger 
if arguments["-v"]:
    logging.basicConfig(level=logging.DEBUG,\
        format='%(asctime)s -- %(levelname)s -- %(message)s')
else:
    logging.basicConfig(level=logging.INFO,\
        format='%(asctime)s -- %(levelname)s -- %(message)s')
logger = logging.getLogger('filtflatunw_log.log')

# initialise job list
jobs = Job(do_list)
print('List of Post-Processing Jobs:')
for p in jobs:
    print(p)
    # test if job in available list
    job = getattr(p,'name')
    if callable(eval(job)):
        pass
    else:
        print('Job: {0} does not exist'.format(job))
        print(Job.__doc__)
        sys.exit()
print()

if arguments["-v"]:
    print(FiltFlatUnw.__doc__)
    print()

if arguments["-v"]:
    # give me the time to read terminal
    time.sleep(3)

if arguments["--look"] == None:
    # init look to Rlooks_int
    look = str(np.copy(proc["Rlooks_int"]))
else:
    look = str(int(arguments["--look"]))

# RUN
for p in jobs:
    postprocess = FiltFlatUnw(
        [ListInterfero,SARMasterDir,IntDir,EraDir,
        proc["Rlooks_int"], proc["Rlooks_unw"], 
        proc["nfit_range"], proc["thresh_amp_range"],
        proc["nfit_az"], proc["thresh_amp_az"],
        proc["filterstyle"], proc["SWwindowsize"], proc["SWamplim"],
        proc["FilterStrength"], proc["Filt_method"], 
        proc["nfit_topo"], proc["thresh_amp_topo"], proc["ivar"], proc["z_ref"], proc["min_z"], proc["delta_z"],
        proc["seedx"], proc["seedy"], proc["threshold_unw"], proc["threshold_unfilt"], proc["unw_method"],proc["ref_top"], proc["ref_left"],proc["ref_width"],proc["ref_length"]], 
        prefix=prefix, suffix=suffix, look=look, model=model, cutfile=cutfile, force=force, 
        ibeg_mask=ibeg_mask, iend_mask=iend_mask, jbeg_mask=jbeg_mask, jend_mask=jend_mask,
        ) 

    print()
    job = getattr(p,'name')
    # print ifgs list at the begining of each process
    postprocess.stack.info()

    print('----------------------------------')
    print('Run {} ....'.format(job))
    
    # run process
    output = []
    output.append(go(postprocess, job, nproc))

    # update name
    prefix, suffix, look = output[0][0][:3]

    # load list of dates
    dates1, dates2 = np.loadtxt(ListInterfero,comments="#",unpack=True,usecols=(0,1),dtype='i,i') 
    dates = np.vstack([dates1,dates2]).T; Nifg = len(dates1)

    # update list ifg
    success = []; [success.append(output[0][i][3]) for i in range(Nifg)]; success = np.array(success)
    index = np.flatnonzero(success==1); newdates =  dates[index,:] 

    # save new list
    ListInterfero = path.join(path.abspath(home) + '/' + "interf_pair_success.txt")
    logger.info("Save successfull list of interferograms in {}".format(ListInterfero))
    wf = open(ListInterfero, 'w')
    for i in range(len(index)):
        wf.write("%i  %i\n" % (newdates[i][0], newdates[i][1]))
    wf.close()

    print('----------------------------------')
    print()

print("That's all folks")

