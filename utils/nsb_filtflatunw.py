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
        self._processes = self._processes._make(job)

    def replace_job(self,pos,value):
        self._processes = self._processes[pos]._replace(value)

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

    def getstratfile(self,kk):
        ''' Return stratified file name '''
        return  str(self._ifgs[kk].date1) + '-' + str(self._ifgs[kk].date2) + '_strat_' +  self._ifgs[kk].look + 'rlks.unw'

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

    def __init__(self, params, ibeg_mask=0, iend_mask=0, jbeg_mask=0, jend_mask=0):
        (self.ListInterfero, self.SARMasterDir, self.IntDir,
        self.Rlooks_int, self.Rlooks_unw, self.prefix, self.suffix, 
        self.nfit_range, self.thresh_amp_range,
        self.nfit_az, self.thresh_amp_az,
        self.filterstyle,self.SWwindowsize, self.SWamplim,
        self.nfit_atmo,self.thresh_amp_atmo, self.ivar, self.z_ref,
        self.seedx, self.seedy,
        ) = map(str, params)

        # initiliase number of looks
        self.look = self.Rlooks_int
        self.rlook = int(int(self.Rlooks_unw) - int(self.Rlooks_int))

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
        ''' Look Radar function '''
        
        dirname, filename = path.split(path.abspath(file))
        chdir(dirname) 

        r= subprocess.call("look.pl "+str(filename)+" "+str(self.rlook)+" >> log_look.txt" , shell=True)
        if r != 0:
            logger.warning(r)
            logger.warning(' Can''t look file {0} in {1} look'.format(filename,self.rlook))

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

        logger.debug('Filter {0} with {1} filter type'.format(infile,self.filterstyle))
        r = subprocess.call("nsb_SWfilter.pl "+str(inbase)+" "+str(filtbase)+" "+str(corbase)\
                +" "+str(self.SWwindowsize)+" "+str(self.SWamplim)+" "+str(self.filterstyle), shell=True)
        if r != 0:
            logger.warning('Failed filtering {0} with {1} filter type'.format(infile,self.filterstyle))

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
            r = subprocess.call("flatten_range "+str(infile)+" "+str(filtfile)+" "+str(outfile)+" "+str(filtout)+\
                " "+str(self.nfit_range)+" "+str(self.thresh_amp_range)+"  >> log_flatenrange.txt", shell=True)
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
        inrsc = infile + '.rsc'
        corfile = self.stack.getpath(kk) + '/' + self.stack.getcor(kk)
        filtfile = self.stack.getpath(kk) + '/'+ self.stack.getfilt(kk)

        # update names
        prefix, suffix = self.stack.getfix(kk)
        newsuffix = suffix + '_flataz'
        self.stack.updatefix(kk,prefix,newsuffix)
        outfile = self.stack.getpath(kk) + '/'+ self.stack.getname(kk)
        filtout = self.stack.getpath(kk) + '/'+ self.stack.getfilt(kk)
        outrsc = outfile + '.rsc'

        if path.exists(filtfile) == False:
            logger.debug('{0} does not exist'.format(filtfile))
            # call filter function
            eval(self.filter(kk))

        if path.exists(outfile) == False:
            logger.debug('Flatten azimuth on IFG: {0}'.format(infile))
            r = subprocess.call("flatten_az "+str(infile)+" "+str(filtfile)+" "+str(outfile)+" "+str(filtout)+\
                " "+str(nfit_az)+" "+str(thresh_amp_az), shell=True)
            if r != 0:
                logger.warning("Flatten azimuth failed for int. {0}-{1}".format(date1,date2))
            else:
                shutil.copy(inrsc,outrsc)
        else:
            logger.debug('Flatten azimuth on IFG: {0} already done'.format(infile))

    def flat_topo(self, kk):
        ''' Faltten topo function '''
        chdir(self.stack.getpath(kk))

        infile = self.stack.getname(kk)
        corfile = self.stack.getcor(kk)
        filtfile = self.stack.getfilt(kk)
        stratfile = self.stack.getstratfile(kk)

        # update names
        prefix, suffix = self.stack.getfix(kk)
        newsuffix = suffix + '_flatz'
        self.stack.updatefix(kk,prefix,newsuffix)
        outfile = self.stack.getpath(kk) + '/'+ self.stack.getname(kk)
        filtout = self.stack.getpath(kk) + '/'+ self.stack.getfilt(kk)

        # look dem
        rscin = self.dem + '.rsc'
        self.look_file(self.dem)
        self.dem =  self.SARMasterDir + '/'+  'radar_' + self.Rlooks_unw + 'rlks.hgt'
        rscout = self.dem + '.rsc'
        shutil.copy(rscin,rscout)
        del rscin,rscout

        if path.exists(filtfile) == False:
            logger.debug('{0} does not exist'.format(filtfile))
            # call filter function
            eval(self.filter(kk))

        if path.exists(outfile) == False:
            logger.debug('Flatten topo on IFG: {0}'.format(infile))
            r = subprocess.call("flatten_topo "+str(infile)+" "+str(filtfile)+" "+str(self.dem)+" "+str(outfile)+" "+str(filtout)\
                +" "+str(self.nfit_atmo)+" "+str(self.ivar)+" "+str(self.z_ref)+" "+str(self.thresh_amp_atmo)+" "+str(stratfile)+" >> log_flattopo.txt", shell=True)
            if r != 0:
                logger.warning("Flatten topo failed for int. {0}".format(infile))
            else:
                inrsc = infile + '.rsc'
                outrsc = outfile + '.rsc'
                filtrsc = filtout + '.rsc'
                shutil.copy(inrsc,outrsc)
                shutil.copy(inrsc,filtrsc)

            # select points
            i, j, z, phi, coh, deltaz = np.loadtxt('ncycle_topo',comments='#', usecols=(0,1,2,3,5,10), unpack=True,dtype='f,f,f,f,f,f')
            z = z - self.z_ref
            phi = phi*0.00020944

            topfile = path.splitext(infile)[0] + '.top'
            b1, b2, b3, b4, b5 =  np.loadtxt(topfile,usecols=(0,1,2,3,4), unpack=True, dtype='f,f,f,f,f')

            if ((self.jend_mask > self.jbeg_mask) or (self.iend_mask > self.ibeg_mask)) and self.ivar<2 :
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
                r = subprocess.call("cpx2rmg.pl "+str(infile)+" "+str(infileunw), shell=True)
                if r != 0:
                    logger.warning("Failed to convert {0} to .unw".format(infile))
                inrsc = infile + '.rsc'
                outrsc = infileunw + '.rsc'
                shutil.copy(inrsc,outrsc)

                # remove strat file created by flatten_topo and write
                remove(stratfile)
                logger.debub("Create stratified file {0}: ".format(strat))
                r = subprocess.call("write_strat_unw "+str(strattxt)+" "+str(infileunw)+" "+str(self.dem)+" "+str(stratfile)\
                        +" "+str(nfit_atmo)+" "+str(ivar)+" "+str(self.z_ref)+" >> log_flattopo.txt", shell=True)
                if r != 0:
                    logger.warning("Failed creating stratified file: {0}".format(stratfile))
                outrsc = stratfile + '.rsc'
                shutil.copy(inrsc,outrsc)

                # remove model created by flatten_topo and run
                remove(outfile)
                remove(filtout)
                logger.debub("Remove stratified file {0} from IFG {1}: ".format(stratfile,infile))
                r = subprocess.call("removeModel.pl "+str(infile)+" "+str(stratfile)+" "+str(outfile), shell=True)
                if r != 0:
                    logger.warning("Failed removing stratified file: {0} from IFG {1}: ".format(stratfile,infile))

                corfile = self.stack.getcor(kk)
                corbase = path.splitext(corfile)[0]
                logger.debub("Filter IFG {0}: ".format(outfile))
                r = subprocess.call("nsb_SWfilter.pl "+str(path.splitext(outfile)[0])+" "+str(path.splitext(filtout)[0])+" "+str(corbase)\
                    +" "+str(self.SWwindowsize)+" "+str(self.SWamplim)+" "+str(self.filterstyle))
                if r != 0:
                    logger.warning("Failed filtering IFG {0}: ".format(outfile))
            else:
                z_select = z; phi_select = phi
                # I dont understand ivar=0
                if self.ivar == 1:
                    fit = b1*z_select + (b2/2.)*z_select**2 + (b3/3.)*z_select**2 + (b4/4.)*z_select**2

            # plot phase/topo
            fig = plt.figure(nfigure)
            nfigure =+ 1
            ax = fig.add_subplot(1,1,1)
            # lets not plot dphi but phi
            ax.plot(z,phi*z,'.',alpha=.6)
            ax.plot(z_select,fit,'-r',lw=4,label='Fit: {0:.3f}z + {1:.3f}z2 + {2:.3f}z3 + {3:.3f}z4'.format(b1,b2,b3,b4))
            ax.set_xlabel('Elevation (m)')
            ax.set_ylabel('Phase (rad)')
            plt.legend(loc='best')
            plotfile = path.splitext(infile)[0] + '_phase-topo.png'
            fig.savefig(plotfile, format='PNG')
            plt.show()

        else:
            logger.debug('Flatten topo on IFG: {0} already done'.format(infile))

    def flat_model(kk):
        return

    def colin(kk):
        return 

    def unwrapping(kk):
        return 

    def addback(kk):
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
nfit_topo=-1
thresh_amp_topo=0.2
ivar=1
z_ref=8000.

####################
# Test Process List
####################

# Job list: look_int replace_amp   filter     flat_range    flat_az    flat_topo  flat_model      unw         add_back
do_list =   [True,     True,       True,       True,       False,       True,      False,      False,      False] 
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
        nfit_topo,thresh_amp_topo,ivar,z_ref,
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


