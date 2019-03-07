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
from os import path, environ, system, chdir, remove, getcwd
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
import logging
import multiprocessing

##################################################################################
###  INITIALISE
##################################################################################

# init logger 
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger('filtcorunw')

# init collections 
import collections
Process = collections.namedtuple('Process', 'name')
IFG = collections.namedtuple('IFG', 'date1 date2 look prefix suffix width length')
Image = collections.namedtuple('Image', 'date decimal_date temporal_baseline')

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

class Job():
    """ Create a class of Jobs to be run: 
    Job list is: erai look_int replace_amp filterSW filterROI flat_range flat_topo flat_model colin unwrapping add_model_back add_atmo_back add_ramp_back """

    def __init__(self, names):
        self.names = names.split()
        try:
            self._processes = [Process(name) for name in self.names]
        except ValueError as error:
            logger.warning(error)
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
            self._ifgs = [IFG(date1,date2,look,prefix,suffix,width,length) for (date1,date2) in zip(self.dates1,self.dates2)]
        except ValueError as error:
            logger.warning(error)
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
        return  self.dir + 'int_' + str(self.dates1[kk]) + '_' + str(self.dates2[kk]) 

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
    nfit_range, hresh_amp_range, nfit_az, thresh_amp_az, filterstyle,SWwindowsize, SWamplim, filterStrength, nfit_atmo,thresh_amp_atmo, ivar, z_ref,
    seedx, seedy.threshold_unw, unw_method
    Additional parameters not in the proc file (yet?): ibeg_mask, iend_mask, jbeg_mask, jend_mask (default: 0.)
    defining the boundary of the mask zone for emprical estimations
    suffix, preffix: define name of the interferogram at the start of the processes
    model: model to be removed from wrapped interferograms (default: None)
    """

    def __init__(self, params, prefix='', siffix='_sd', ibeg_mask=0, iend_mask=0, jbeg_mask=0, jend_mask=0, model=None):
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
        self.strat = None

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

    def look_file(self,file):
        ''' Look function 
        Requiered parameters:  Rlooks_int, Rlooks_unw
        '''
        
        dirname, filename = path.split(path.abspath(file))
        chdir(dirname) 

        print("look.pl "+str(filename)+" "+str(self.rlook))
        r= subprocess.call("look.pl "+str(filename)+" "+str(self.rlook)+" >> log_look.txt" , shell=True)
        if r != 0:
            logger.warning(' Can''t look file {0} in {1} look'.format(filename,self.rlook))
            print(self.look_file.__doc__)

    def erai(self,kk):
        return

    def computesize(self,file):
        ''' Open file with gdal and retrieve width anf length 
        '''
        
        try:
            dirname, filename = path.split(path.abspath(file))
            chdir(dirname) 

            ds_int = gdal.Open(filename, gdal.GA_ReadOnly)
            driver = ds_int.GetDriver()
            return ds_int.RasterXSize, ds_int.RasterYSize
        except ValueError as error:
            logger.warning(error)
            print(self.computesize.__doc__)

    def replace_amp(self, kk):
        ''' Replace amplitude by coherence'''

        chdir(self.stack.getpath(kk))

        infile = self.stack.getname(kk)+ '.int'
        rscfile = infile + '.rsc'
        corfile = self.stack.getcor(kk) 
        logger.info('Replace Amplitude by Coherence on IFG: {0}'.format(infile))

        # compute width and length
        width,length = self.computesize(infile)
        self.stack.updatesize(kk,width,length)

        # update names
        prefix, suffix = self.stack.getfix(kk)
        newprefix = 'coh_' + prefix
        self.stack.updatefix(kk,newprefix,suffix)
        outfile = self.stack.getname(kk) + '.int'
        outrsc = outfile + '.rsc' 
        shutil.copy(rscfile,outrsc)

        # check if not done
        if path.exists(outfile) is False:
            if path.exists(corfile):

                logger.info('Replace Amplitude by Cohrence for IFG: {}'.format(infile))
                print("rmg2mag_phs "+str(corfile)+" tmp cor "+str(width))
                r1 = subprocess.call("rmg2mag_phs "+str(corfile)+" tmp cor "+str(width)+"  >> log_replaceAMP.txt", shell=True)
                if r1 != 0:
                    logger.warning(r1)
                    logger.warning('Replace Amplitude by Cohrence for IFG: {}'.format(infile))
                    sys.exit()

                print("cpx2mag_phs "+str(infile)+" tmp2 phs "+str(width))
                r2 = subprocess.call("cpx2mag_phs "+str(infile)+" tmp2 phs "+str(width)+"  >> log_replaceAMP.txt", shell=True)
                if (r2) != 0:
                    logger.warning(r2)
                    logger.warning('Replace Amplitude by Cohrence for IFG: {}'.format(infile))
                    sys.exit()

                print("mag_phs2cpx cor phs "+str(outfile)+" "+str(width))
                r3 = subprocess.call("mag_phs2cpx cor phs "+str(outfile)+" "+str(width)+"  >> log_replaceAMP.txt", shell=True)
                if (r3) != 0:
                    logger.warning(r3)
                    logger.warning('Replace Amplitude by Cohrence for IFG: {}'.format(infile))
                    sys.exit()

                remove('tmp'); remove('tmp2'); remove('phs'); remove('cor')

            else:
                logger.warning('Coherence file does not exit...')
                print(self.replace_amp.__doc__)
                sys.exit()

        else:
            logger.debug('Replace Amplitude by Cohrence for IFG: {} already done'.format(infile))
            print('{0} exists, assuming OK'.format(outfile))

    def filterSW(self, kk):
        ''' Filter SW function form Doin et. al. 2011
        Requiered proc parameters: SWwindowsize, SWamplim, filterstyle
        '''

        chdir(self.stack.getpath(kk))

        inbase = self.stack.getname(kk) 
        inrsc = inbase + '.int.rsc'
        infile = inbase + '.int'
        corfile = self.stack.getcor(kk) 
        corbase = path.splitext(corfile)[0]
        filtbase = self.stack.getfiltSW(kk)
        filtrsc = filtbase + '.int.rsc'

        logger.info('Filter {0} with {1} filter type'.format(infile,self.filterstyle))
        r = subprocess.call("nsb_SWfilter.pl "+str(inbase)+" "+str(filtbase)+" "+str(corbase)\
                +" "+str(self.SWwindowsize)+" "+str(self.SWamplim)+" "+str(self.filterstyle), shell=True)
        if r != 0:
            logger.warning('Failed filtering {0} with {1} filter type'.format(infile,self.filterstyle))
            print(self.filterSW.__doc__)
            sys.exit()

        if path.exists(filtrsc) == False:
            shutil.copy(inrsc,filtrsc)

    def filterROI(self, kk):
        ''' ROI-PAC Filter function
        Requiered proc file parameter: filterStrength
        '''

        chdir(self.stack.getpath(kk))

        infile = self.stack.getname(kk) + '.int'
        inrsc = infile + '.rsc'
        
        # get width and compute if not already done
        width,length =  self.stack.getsize(kk)
        if (int(width) == 0) or (int(length) == 0):
            width,length = self.computesize(infile)
            self.stack.updatesize(kk,width,length)

        filtfile = self.stack.getfiltROI(kk) + '.int'
        filtrsc = filtfile + '.rsc'
        if path.exists(filtrsc) == False:
            shutil.copy(inrsc,filtrsc)

        logger.info('Filter {0} with ROI-PAC adaptative filter'.format(infile))
        print("myadapt_filt "+str(infile)+" "+str(filtfile)+" "+str(width)+" 0.25"+" "+str(self.filterStrength))
        r = subprocess.call("myadapt_filt "+str(infile)+" "+str(filtfile)+" "\
                +str(width)+" 0.25"+" "+str(self.filterStrength)+"  >> log_filtROI.txt", shell=True)
        if r != 0:
            logger.warning('Failed filtering {0} with ROI-PAC adaptative filter'.format(infile))
            print(self.filterROI.__doc__)
            sys.exit()
        
    def flat_range(self,kk):
        ''' Function flatten range  on wrapped phase  (See Doin et al., 2015)
        Requiered proc file parameters: nfit_range, thresh_amp_range
        Estimation done on filterSW file
        '''

        chdir(self.stack.getpath(kk))

        infile = self.stack.getname(kk) + '.int'
        inrsc = infile + '.rsc'
        corfile = self.stack.getcor(kk)
        filtfile = self.stack.getfiltSW(kk) + '.int'

        # update names
        prefix, suffix = self.stack.getfix(kk)
        newsuffix = suffix + '_flatr'
        self.stack.updatefix(kk,prefix,newsuffix)
        outfile = self.stack.getname(kk) + '.int' 
        outrsc = outfile + '.rsc'
        filtout = self.stack.getfiltSW(kk)
        shutil.copy(inrsc,outrsc)

        if path.exists(filtfile) == False:
            logger.warning('{0} does not exist'.format(filtfile))
            # call filter function
            self.filterSW(kk)

        if path.exists(outfile) == False:
            logger.info('Flatten range on IFG: {0}'.format(infile))
            print("flatten_range "+str(infile)+" "+str(filtfile)+" "+str(outfile)+" "+str(filtout)+\
                " "+str(self.nfit_range)+" "+str(self.thresh_amp_range))
            r = subprocess.call("flatten_range "+str(infile)+" "+str(filtfile)+" "+str(outfile)+" "+str(filtout)+\
                " "+str(self.nfit_range)+" "+str(self.thresh_amp_range)+"  >> log_flatenrange.txt", shell=True)
            if r != 0:
                logger.warning("Flatten range failed for IFG: {0}".format(infile))
                print(self.flat_range.__doc__) 
                sys.exit()
        else:
            logger.debug('Flatten range on IFG: {0} already done'.format(infile))
            print('{0} exists, assuming OK'.format(outfile))

    def flat_az(self,kk):
        ''' Function flatten azimuth  on wrapped phase  (See Doin et al., 2015)
            Requiered proc file parameters: nfit_az, thresh_amp_az
            Estimation done on filterSW file
        '''

        # need to be in the corect dir for stupid perl scripts
        chdir(self.stack.getpath(kk))

        ''' Faltten Azimuth function '''
        infile = self.stack.getname(kk) + '.int'
        inrsc = infile + '.rsc'
        corfile = self.stack.getcor(kk) 
        filtfile = self.stack.getfiltSW(kk) + '.int'

        # update names
        prefix, suffix = self.stack.getfix(kk)
        newsuffix = suffix + '_flataz'
        self.stack.updatefix(kk,prefix,newsuffix)
        outfile = self.stack.getname(kk) + '.int' 
        filtout = self.stack.getfiltSW(kk) + '.int'
        outrsc = outfile + '.rsc'
        shutil.copy(inrsc,outrsc)

        if path.exists(filtfile) == False:
            logger.info('{0} does not exist'.format(filtfile))
            # call filter function
            eval(self.filterSW(kk))

        if path.exists(outfile) == False:
            logger.info('Flatten azimuth on IFG: {0}'.format(infile))
            print("flatten_az "+str(infile)+" "+str(filtfile)+" "+str(outfile)+" "+str(filtout)+\
                " "+str(nfit_az)+" "+str(thresh_amp_az))

            r = subprocess.call("flatten_az "+str(infile)+" "+str(filtfile)+" "+str(outfile)+" "+str(filtout)+\
                " "+str(nfit_az)+" "+str(thresh_amp_az), shell=True)
            if r != 0:
                logger.warning("Flatten azimuth failed for int. {0}-{1}".format(date1,date2))
                print(self.flat_az.__doc__)
                sys.exit()
        else:
            logger.debug('Flatten azimuth on IFG: {0} already done'.format(infile))
            print('{0} exists, assuming OK'.format(outfile))

    def flat_topo(self, kk):
        ''' Function flatten atmosphere on wrapped phase  (See Doin et al., 2015)
        Requiered proc file parameters: nfit_atmo, ivar, z_ref, thresh_amp_atmo
        Estimation done on filterSW file
        Plot phase/topo in *_phase-topo.png file
        '''

        chdir(self.stack.getpath(kk))
        infile = self.stack.getname(kk) + '.int'
        corfile = self.stack.getcor(kk)
        filtfile = self.stack.getfiltSW(kk) + '.int'
        stratfile = self.stack.getstratfile(kk) + '.unw'

        # filt must be done before changing name
        if path.exists(filtfile) == False:
            logger.info('{0} does not exist'.format(filtfile))
            # call filter function
            self.filterSW(kk)

        # update names
        prefix, suffix = self.stack.getfix(kk)
        newsuffix = suffix + '_flatz'
        self.stack.updatefix(kk,prefix,newsuffix)
        outfile = self.stack.getname(kk) + '.int'
        filtout = self.stack.getfiltSW(kk) + '.int'

        # check width IFG
        width,length = self.stack.getsize(kk)
        if (int(width) == 0) or (int(length) == 0):
            width,length = self.computesize(infile)
            self.stack.updatesize(kk,width,length)

        # look dem if necessary
        w,l = self.computesize(self.dem)

        if int(w) != int(width):
            logger.warning('IFG:{0} and DEM file are not the same size: {0}'.format(infile))
            self.look_file(self.dem)
            # update DEM
            self.dem = self.SARMasterDir + '/'+  'radar_' + self.Rlooks_unw + 'rlks.hgt'
            
        # ca c'est pourri...
        chdir(self.stack.getpath(kk))

        if path.exists(outfile) == False:
            logger.info('Flatten topo on IFG: {0}'.format(infile))
            print("flatten_topo "+str(infile)+" "+str(filtfile)+" "+str(self.dem)+" "+str(outfile)+" "+str(filtout)\
                +" "+str(self.nfit_atmo)+" "+str(self.ivar)+" "+str(self.z_ref)+" "+str(self.thresh_amp_atmo)+" "+\
                str(stratfile))
            r = subprocess.call("flatten_topo "+str(infile)+" "+str(filtfile)+" "+str(self.dem)+" "+str(outfile)+" "+str(filtout)\
                +" "+str(self.nfit_atmo)+" "+str(self.ivar)+" "+str(self.z_ref)+" "+str(self.thresh_amp_atmo)+" "+\
                str(stratfile)+" >> log_flattopo.txt", shell=True)
            if r != 0:
                logger.warning("Flatten topo failed for int. {0}".format(infile))
                print(self.flat_topo.__doc__)
                sys.exit()
        else:
            print('{0} exists, assuming OK'.format(outfile))
        
        inrsc = infile + '.rsc'
        outrsc = outfile + '.rsc'
        print(inrsc,outrsc)
        filtrsc = filtout + '.rsc'
        shutil.copy(inrsc,outrsc)
        shutil.copy(inrsc,filtrsc)

        # select points
        i, j, z, phi, coh, deltaz = np.loadtxt('ncycle_topo',comments='#', usecols=(0,1,2,3,5,10), unpack=True,dtype='f,f,f,f,f,f')
        z = z - float(self.z_ref)
        phi = phi*0.00020944

        topfile = path.splitext(infile)[0] + '.top'
        b1, b2, b3, b4, b5 =  np.loadtxt(topfile,usecols=(0,1,2,3,4), unpack=True, dtype='f,f,f,f,f')

        if ((self.jend_mask > self.jbeg_mask) or (self.iend_mask > self.ibeg_mask)) and self.ivar<2 :
            sys.exit(0)
            b1, b2, b3, b4, b5 = 0, 0, 0, 0, 0

            index = np.nonzero(
            np.logical_and(coh>self.thresh_amp_atmo,
            np.logical_and(deltaz>75.,
            np.logical_and(np.logical_or(i<self.ibeg_mask,pix_az>self.iend_mask),
            np.logical_or(j<self.jbeg_mask,j>self.jend_mask),
            ))))

            phi_select = phi[index]
            z_select = z[index]

            if self.nfit_atmo == -1:
                b1 = np.nanmedian(phi_select)
                fit = z_select*b1
            elif self.nfit_atmo == 0:
                b1 = np.nanmean(phi_select)
                fit = z_select*b1
            elif self.nfit_atmo == 1:
                from sklearn.linear_model import LinearRegression
                model = LinearRegression()
                model.fit(z_select, phi_select)
                fit = model.predict(z_select)
            else:
                from sklearn.preprocessing import PolynomialFeatures
                polynomial_features= PolynomialFeatures(degree=self.nfit_atmo)
                x_poly = polynomial_features.fit_transform(z_select)
                model = LinearRegression()
                model.fit(x_poly, phi_select)
                fit = model.predict(x_poly)
            
            # save median phase/topo
            strattxt = path.splitext(infile)[0] + '_strat.top'
            np.savetxt(strattxt, median, fmt=('%.8f'))

            ax.plot(z_select,phi_selct*z_select,'.r',label='selected points')
            av = np.median(phi_select*z_select)
            
            # clean and prepare for write strat
            infileunw = path.splitext(infile)[0] + '.unw'
            remove(infileunw)

            logger.debub("Convert {0} to .unw".format(infile))
            print("cpx2rmg.pl "+str(infile)+" "+str(infileunw))
            r = subprocess.call("cpx2rmg.pl "+str(infile)+" "+str(infileunw), shell=True)
            if r != 0:
                logger.warning("Failed to convert {0} to .unw".format(infile))
                sys.exit()

            inrsc = infile + '.rsc'
            outrsc = infileunw + '.rsc'
            shutil.copy(inrsc,outrsc)

            # remove strat file created by flatten_topo and write
            remove(stratfile)
            logger.debub("Create stratified file {0}: ".format(strat))
            print("write_strat_unw "+str(strattxt)+" "+str(infileunw)+" "+str(self.dem)+" "+str(stratfile)\
                    +" "+str(nfit_atmo)+" "+str(ivar)+" "+str(self.z_ref))
            r = subprocess.call("write_strat_unw "+str(strattxt)+" "+str(infileunw)+" "+str(self.dem)+" "+str(stratfile)\
                    +" "+str(nfit_atmo)+" "+str(ivar)+" "+str(self.z_ref)+" >> log_flattopo.txt", shell=True)
            if r != 0:
                logger.warning("Failed creating stratified file: {0}".format(stratfile))
                sys.exit()

            outrsc = stratfile + '.rsc'
            shutil.copy(inrsc,outrsc)

            # remove model created by flatten_topo and run
            remove(outfile)
            remove(filtout)
            logger.debub("Remove stratified file {0} from IFG {1}: ".format(stratfile,infile))
            print("removeModel.pl "+str(infile)+" "+str(stratfile)+" "+str(outfile))
            r = subprocess.call("removeModel.pl "+str(infile)+" "+str(stratfile)+" "+str(outfile), shell=True)
            if r != 0:
                logger.warning("Failed removing stratified file: {0} from IFG {1}: ".format(stratfile,infile))
                sys.exit()

            corfile = self.stack.getcor(kk)
            corbase = path.splitext(corfile)[0]
            logger.debub("Filter IFG {0}: ".format(outfile))
            print("nsb_SWfilter.pl "+str(path.splitext(outfile)[0])+" "+str(path.splitext(filtout)[0])+" "+str(corbase)\
                +" "+str(self.SWwindowsize)+" "+str(self.SWamplim)+" "+str(self.filterstyle))
            r = subprocess.call("nsb_SWfilter.pl "+str(path.splitext(outfile)[0])+" "+str(path.splitext(filtout)[0])+" "+str(corbase)\
                +" "+str(self.SWwindowsize)+" "+str(self.SWamplim)+" "+str(self.filterstyle), shell=True)
            if r != 0:
                logger.warning("Failed filtering IFG {0}: ".format(outfile))
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
        self.strat = True

    def look_int(self,kk):
        ''' Look function for IFG, coherence, strat, and radar files
        Requiered parameters:  Rlooks_int, Rlooks_unw
        '''

        # need to be in the corect dir for stupid perl scripts
        chdir(self.stack.getpath(kk))

        infile =  self.stack.getname(kk) + '.int'
        corfile =  self.stack.getcor(kk) 

        # look strat file
        stratfile = self.stack.getstratfile(kk) + '.unw'
        if (path.exists(stratfile) == False) and (self.strat == True):
            # self.strat check if flatten_topo has been done
            self.look_file(stratfile)
            self.strat = stratfile
 
        chdir(self.stack.getpath(kk))

        # look radar file if not done
        dem = self.SARMasterDir + '/'+  'radar_' + self.Rlooks_unw + 'rlks.hgt'
        if path.exists(self.dem) is False:
            self.look_file(self.dem)
            self.dem = dem

        logger.info('Look file {0} in {1} look'.format(infile,self.rlook))
        chdir(self.stack.getpath(kk))

        # update looks
        self.stack.updatelook(kk,self.Rlooks_unw)
        outfile =  self.stack.getname(kk) + '.int'
        outcor = self.stack.getcor(kk) 

        if path.exists(outfile) is False:
            print("look.pl "+str(infile)+" "+str(self.rlook))
            r= subprocess.call("look.pl "+str(infile)+" "+str(self.rlook)+" >> log_look.txt" , shell=True)
            if r != 0:
                logger.warning(' Can''t look file {0} in {1} look'.format(infile,self.rlook))
                print(self.look_int.__doc__)
        else:
            print('{0} exists, assuming OK'.format(outfile))
        
        if path.exists(outcor) is False:
            print("look.pl "+str(corfile)+" "+str(self.rlook))
            r = subprocess.call("look.pl "+str(corfile)+" "+str(self.rlook)+" >> log_look.txt", shell=True)
            if r != 0:
                logger.warning(' Can''t look file {0} in {1} look'.format(corfile,self.rlook))
                print(self.look_int.__doc__)
                
            r = subprocess.call("look.pl "+str(self.dem)+" "+str(self.rlook)+" >> log_look.txt", shell=True)
            if r != 0:
                logger.warning(' Can''t look file {0} in {1} look'.format(corfile,self.rlook))
                print(self.look_int.__doc__)
        else:
            print('{0} exists, assuming OK'.format(outfile))
                
        # update size
        width,length = self.computesize(outcor)
        self.stack.updatesize(kk,width,length)

    def flat_model(self,kk):
        return

    def colin(self,kk):
        ''' Compute and replace amplitude by colinearity (See Pinel-Puyssegur et al., 2012)'''

        chdir(self.stack.getpath(kk))

        infile = self.stack.getname(kk) + '.int'
        inrsc = infile + '.rsc'
        filtfile = self.stack.getfiltSW(kk) + '.int'

        # update names
        prefix, suffix = self.stack.getfix(kk)
        newprefix = 'col_'
        self.stack.updatefix(kk,newprefix,suffix)
        outfile = self.stack.getname(kk) + '.int'
        outrsc = outfile + '.rsc'
        filtout = self.stack.getfiltSW(kk) + '.int'
        filtrsc = filtout + '.rsc'
        filtoutroi = self.stack.getfiltROI(kk)+ '.int'
        filtroirsc = filtoutroi + '.rsc'

        shutil.copy(inrsc,outrsc)
        shutil.copy(inrsc,filtrsc)
        shutil.copy(inrsc,filtroirsc)

        # Retrieve length and width
        width,length =  self.stack.getsize(kk)
        if (int(width) == 0) or (int(length) == 0):
            width,length = self.computesize(infile)
            self.stack.updatesize(kk,width,length)

        if path.exists(outfile) == False:
            logger.info('Replace Amplitude by colinearity on IFG: {0}'.format(infile))
            shutil.copy(infile,'temp')
            print("colin "+str(infile)+" temp "+str(outfile)+" "+str(width)+" "+str(length)+\
                " 3 0.0001 2")
            r = subprocess.call("colin "+str(infile)+" temp "+str(outfile)+" "+str(width)+" "+str(length)+\
                " 3 0.0001 2  >> log_flatenrange.txt", shell=True)
            if r != 0:
                logger.warning('Failed replacing Amplitude by colinearity on IFG: {0}'.format(infile))
                print(self.colin.__doc__)
                sys.exit()
            # clean
            remove('temp')

        else:
            logger.debug('Colinearity on IFG {0} already computed'.format(infile))
            print('{0} exists, assuming OK'.format(outfile))

        # # Filter with colinearity for unwrapping
        # if path.exists(filtout) == False:
        #     self.filterSW(kk)
        # if path.exists(filtoutroi) == False:
        #     self.filterROI(kk)

    def unwrapping(self,kk):
        ''' Unwrap function from strating seedx, seedy
        if unw_method: mpd, MP.DOIN algorthim (Grandin et al., 2012). Requiered: threshold_unw, filterSW, filterROI
        if unw_method: roi, ROIPAC algorthim. Requiered threshold_unw, filterROI
        '''

        chdir(self.stack.getpath(kk))

        infile = self.stack.getname(kk)+ '.int'
        inrsc = infile + '.rsc'
        filtSWfile = self.stack.getfiltSW(kk)+ '.int'
        filtROIfile = self.stack.getfiltROI(kk)+ '.int'

        unwfile = self.stack.getname(kk)+ '.unw'
        unwrsc = unwfile + '.rsc'
        unwfiltSW = self.stack.getfiltSW(kk)+ '.unw'
        unwSWrsc = unwfiltSW + '.rsc'
        unwfiltROI = self.stack.getfiltROI(kk)+ '.unw'
        unwROIrsc = unwfiltROI + '.rsc'
        shutil.copy(inrsc,unwrsc)
        shutil.copy(inrsc,unwSWrsc)
        shutil.copy(inrsc,unwROIrsc)

        # Filter with colinearity
        if path.exists(filtROIfile) == False:
            self.filterROI(kk)

        if path.exists(unwfiltROI) == False:

            print('Unwraped IFG:{0} with strating point col:{1} line:{2} and filterd coherence threshold {3}'.\
                format(unwfile,self.seedx,self.seedy,self.threshold_unw))
            if self.unw_method == 'mpd':
                logger.info("Unwraped IFG:{0} with MP.DOIN algorthim (Grandin et al., 2012) ".format(unwfile))
                if path.exists(unwfiltSW) == False:
                    self.filterSW(kk)

                # my_deroul_interf has ana additional input parameter for threshold on amplitude infile (normally colinearity)
                # unwrapped firt filtSWfile and then add high frequency of filtROIfile
                print("my_deroul_interf_filt "+str(filtSWfile)+" cut "+str(infile)+" "+str(unwfiltROI)\
                    +" "+str(self.seedx)+" "+str(self.seedy)+" "+str(0.04)+" "+str(self.threshold_unw)+" 0")
                r = subprocess.call("my_deroul_interf_filt "+str(filtSWfile)+" cut "+str(infile)+" "+str(unwfiltROI)\
                    +" "+str(self.seedx)+" "+str(self.seedy)+" "+str(0.04)+" "+str(self.threshold_unw)+" 0  >> log_unw.txt", shell=True)
                if r != 0:
                    print(self._unwrapping.__doc__)
                    logger.warning("Failed unwrapping with MP.DOIN algorthim (Grandin et al., 2012)".format(unwfile))
                    sys.exit()
                # remove('cut')

            if self.unw_method == 'roi':

                logger.info("Unwraped IFG:{0} with ROIPAC algorithm ".format(unwfile))
                mask = path.splitext(filtSWfile)[0] + '_msk'

                print("make_mask.pl "+str(path.splitext(filtROIfile)[0])+" "+str(mask)+" "+str(0.02))
                r = subprocess.call("make_mask.pl "+str(path.splitext(filtROIfile)[0])+" "+str(mask)+" "+str(0.02)+"  >> log_unw.txt", shell=True)
                if r != 0:
                    print(self._unwrapping.__doc__)
                    logger.warning("Failed unwrapping IFG {0} with ROIPAC algorithm ".format(unwfile))
                    sys.exit()

                print("new_cut.pl "+str(path.splitext(filtROIfile)[0]))
                r = subprocess.call("new_cut.pl "+str(path.splitext(filtROIfile)[0])+"  >> log_unw.txt", shell=True)
                if r != 0:
                    print(self._unwrapping.__doc__)
                    logger.warning("Failed unwrapping IFG {0} with ROIPAC algorithm ".format(unwfile))
                    sys.exit()

                print("unwrap.pl "+str(path.splitext(filtROIfile)[0])+" "+str(mask)+" "+str(path.splitext(filtROIfile)[0])\
                    +" "+str(self.threshold_unw)+" "+str(self.seedx)+" "+str(self.seedy))
                r = subprocess.call("unwrap.pl "+str(path.splitext(filtROIfile)[0])+" "+str(mask)+" "+str(path.splitext(filtROIfile)[0])\
                    +" "+str(self.threshold_unw)+" "+str(self.seedx)+" "+str(self.seedy)+"  >> log_unw.txt",shell=True)
                if r != 0:
                    print(self._unwrapping.__doc__)
                    logger.warning("Failed unwrapping IFG {0} with ROIPAC algorithm ".format(unwfile))
                    sys.exit()
        else:
            logger.debug("Unwraped IFG:{0} already done  ".format(unwfiltROI))
            print('{0} exists, assuming OK'.format(unwfiltROI))

    def add_atmo_back(self,kk):
        ''' Add back stratified model computed by flatten_topo'''

        chdir(self.stack.getpath(kk))

        # the final product is always filtROI
        unwfile = self.stack.getfiltROI(kk) + '.unw'
        stratfile = self.stack.getstratfile(kk) + '.unw'

        # update names
        prefix, suffix = self.stack.getfix(kk)
        newsuffix = suffix.replace("_flatz", "")
        self.stack.updatefix(kk,prefix,newsuffix)
        outfile = self.stack.getname(kk) + '.unw'

        if path.exists(outfile) == False:
            logger.info('Adding back {0} on IFG: {1}'.format(stratfile,unwfile))
            print("add_rmg.py --infile="+str(unwfile)+" --outfile="+str(outfile)+" --add="+str(stratfile))
            r = subprocess.call("add_rmg.py --infile="+str(unwfile)+" --outfile="+str(outfile)+" --add="+str(stratfile)+\
                 " >> log_flatenrange.txt", shell=True)
            if r != 0:
                logger.warning('Failed adding back {0} on IFG: {1}'.format(stratfile,unwfile))
                logger.warning(r) 
                sys.exit()
        else:
            logger.debug("Adding back {0} on IFG already done: {1}".format(stratfile,unwfile))
            print('{0} exists, assuming OK'.format(outfile))

    def add_ramp_back(self,kk):
        return

    def add_model_back(self,kk):
        return

##################################################################################
###  READ IMPUT PARAMETERS
##################################################################################

# # input parameters (not in the proc file)
home='/home/cometraid14/daouts/work/tibet/qinghai/processing/Sentinel/iw1/'
prefix = '' 
suffix = '_sd'
iend_mask=0 # mask for empirical estimations
jend_mask=0
jbeg_mask=0
jend_mask=0
model=None # model to be removed from wrapped int
nproc=2

# proc file parameters
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
filterStrength=2.
seedx=268
seedy=1766
threshold_unw=0.35
unw_method='roi'

nfit_topo=-1
thresh_amp_topo=0.2
ivar=1
z_ref=8000.


####################
# Test Process List
####################

""" Job list is: erai look_int replace_amp filterSW filterROI flat_range flat_topo flat_model colin unwrapping add_model_back add_atmo_back add_ramp_back """
print(Job.__doc__)
do_list =  'replace_amp filterSW flat_topo colin look_int filterSW unwrapping add_atmo_back'  
jobs = Job(do_list)

print('List of Post-Processing Jobs:')
for job in jobs:
    print(job)
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
        seedx,seedy,threshold_unw,unw_method]
        ) 

# loop over the processes
for p in jobs:
    # check if the process has to be done
    # print(getattr(p,'name'))
    job = getattr(p,'name')
    print('----------------------------------')
    print('Run {} ....'.format(p))
    # [eval('postprocess.{0}({1})'.format(job,kk)) for kk in range(postprocess.Nifg)]
    work = [kk for kk in range(postprocess.Nifg)]
    pool = multiprocessing.Pool(nproc)
    pool.map(eval('postprocess.{0}({1})'.format(job,kk)), work)
    pool.close() 
    print('----------------------------------')
    print()
print("That's all folks")
 


