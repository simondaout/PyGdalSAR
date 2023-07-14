#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Simon DAOUT (CRPG)
############################################

"""\
invers_disp_gps.py
-------------
Temporal inversions of the gps time series displacements

Usage: invers_disp_gps.py --network=<path> --reduction=<path> [--dim=<value>] [--wdir=<path>] [--extension=<value>] [--coseismic=<value>][--postseismic=<value>] [--seasonal=<yes/no>] [--cond=<value>] [--ineq=<value>] [--proj=<value>] [--scale=<value>]  [--bianual=<yes/no>]
invers_disp_gps.py -h | --help

Options:
-h --help           Show this screen 
--network PATH      Text file of name lat, lon for each stations
--reduction PATH    Directory name containing the TS (each TS in a text file untitled by the GPS name: East, North, Down, sigma_East, sigma_North, sigma Down format)      
--dim VALUE         dimension of the gps time series [default:2]
--wdir PATH         Path to the network file and reduction directory [default: ./]
--extension Value   Extention of the the GPS TS files [default: .dat]
--coseismic PATH    Add heaviside functions to the inversion, indicate coseismic time (e.g 2004.,2006.)
--postseismic PATH  Add logarithmic transients to each coseismic step, indicate characteristic time of the log function, must be a serie of values of the same lenght than coseismic (e.g 1.,1.). To not associate postseismic function to a give coseismic step, put None (e.g None,1.) 
--seasonal PATH     If yes, add seasonal terms in the inversion
--bianual=<yes/no>       If yes, add bianual terms in the decomposition [default: no]
--cond VALUE        Condition value for optimization: Singular value smaller than cond*largest_singular_value are considered zero [default: 1.0e-10]
--ineq VALUE        If yes, add inequqlity constrained in the inversion: use lself.east square result to iterate the inversion. Force postseismic to be the same sign than coseismic [default: no].  
--proj VALUE        LOS projection [default for Envisat satellite: 0.38690,-0.09645,0.91706]
--scale VALUE       Scale factor for the three components [default: 1,1,1]
"""

print()
print()
print('Author: Simon Daout')
print()
print()

import numpy as np
from numpy.lib.stride_tricks import as_strided

# scipy
import scipy
import scipy.optimize as opt
import scipy.linalg as lst

import matplotlib.pyplot as plt
import matplotlib.dates as dates
import datetime, math, time

from os import path
import sys

# docopt (command line parser)
import docopt

# read arguments
arguments = docopt.docopt(__doc__)
network = arguments["--network"]
reduction = arguments["--reduction"]
if arguments["--dim"] ==  None:
    dim = 2
else:
    dim = float(arguments["--dim"])
if arguments["--wdir"] ==  None:
    wdir = './'
else:
    wdir = arguments["--wdir"]
if arguments["--extension"] ==  None:
    ext = '.dat'
else:
    ext = arguments["--extension"]    
if arguments["--ineq"] ==  None:
    ineq = 'no'
else:
    ineq = arguments["--ineq"]
if arguments["--cond"] ==  None:
    rcond = 1.0e-10
else:
    rcond = float(arguments["--cond"])
if arguments["--seasonal"] ==  None:
    seasonal = 'no'
else:
    seasonal = arguments["--seasonal"]
if arguments["--bianual"] ==  None:
    arguments["--bianual"] = 'no'
if arguments["--coseismic"] ==  None:
    cos = []
else:
    cos = list(map(float,arguments["--coseismic"].replace(',',' ').split()))
if arguments["--postseismic"] ==  None:
    pos = []
else:
    pos = list(map(float,arguments["--postseismic"].replace('None','-1').replace(',',' ').split()))
if arguments["--proj"] ==  None:
    proj = [0.318,-0.0827,0.9396]
else:
    proj = list(map(float,arguments["--proj"].replace(',',' ').split()))
if arguments["--scale"] ==  None:
    scale = [1.,1.,1.]
else:
    scale = list(map(float,arguments["--scale"].replace(',',' ').split()))

### Define GPS class
class point:
    def __init__(self,x,y):
        self.x=x
        self.y=y

class gpspoint(point):
    def __init__(self,x,y,name):
        point.__init__(self,x,y)
        self.name=name
        self.t=[]
        self.d=[] #self.east,self.north,self.down
        self.sigmad=[] #self.east,self.north, self.down

        # lenght of the time series
        self.Nt = 0

    def info(self):
        print('gps station: {}, self.east(km): {}, self.north(km): {}'.format(self.name,self.x,self.y))

class gpstimeseries:
    def __init__(self,network,reduction,dim,wdir,extension,proj,scale,weight=1.):

        self.network=network
        self.reduction=reduction
        self.dim=dim
        self.wdir=wdir
        self.extension=extension
        self.proj=proj
        self.scale=scale
    
        # inititialisation
        self.sigmad=1./weight
        self.points=[]
        self.Npoints=0
        self.N=0
        self.x,self.y=[],[]
        self.tmin,self.tmax=10000000,0

    def load(self):

        fname=self.wdir + self.network
        if not path.isfile(fname):
            raise ValueError("invalid file name: " + fname)
        else:
            print('Load GPS time series: ', fname)
            print() 
        
        # name, self.east(km), self.north(km)
        name,x,y=np.loadtxt(fname,comments='#',unpack=True,dtype='S4,f,f')

        self.name,self.x,self.y=np.atleast_1d(name,x,y)

        self.Npoints=len(self.name)
        self.N=0
        self.d = []
        
        print('Load time series... ')
        for i in range(self.Npoints):
            station=self.wdir+self.reduction+'/'+self.name[i].decode('utf-8')+self.extension
            print(station)
            if not path.isfile(station):
                raise ValueError("invalid file name: " + station)
                pass
            else:
                self.points.append(gpspoint(self.x[i],self.y[i],self.name[i]))

            if 3==self.dim:
                dated,east,north,up,esigma,nsigma,upsigma=\
                np.loadtxt(station,comments='#',usecols=(0,1,2,3,4,5,6), unpack=True,\
                    dtype='f,f,f,f,f,f,f')

                self.points[i].dated,self.points[i].east,self.points[i].north,self.points[i].up,\
                self.points[i].esigma,self.points[i].nsigma,self.points[i].upsigma=\
                np.atleast_1d(dated,east,north,\
                    up,esigma,nsigma,upsigma)

                self.points[i].d=[self.points[i].east*1000,self.points[i].north*1000,self.points[i].up*1000]
                self.points[i].sigmad=[self.points[i].esigma*1000,self.points[i].nsigma*1000,self.points[i].upsigma*1000]
                self.points[i].comp=['east','north','up']
                self.points[i].t=np.atleast_1d(self.points[i].dated)

            if 2==self.dim:
                dated,east,north,esigma,nsigma=\
                np.loadtxt(station,comments='#',usecols=(0,1,2,3,4), unpack=True,\
                    dtype='f,f,f,f,f')

                self.points[i].dated,self.points[i].east,self.points[i].north,\
                self.points[i].esigma,self.points[i].nsigma,\
                np.atleast_1d(dated,east,north,\
                    esigma,nsigma)

                self.points[i].d=[self.points[i].east*1000,self.points[i].north*1000,self.points[i].up*1000]
                self.points[i].sigmad=[self.points[i].esigma*1000,self.points[i].nsigma*1000,self.points[i].upsigma*1000]
                self.points[i].comp=['east','north']
                self.points[i].t=np.atleast_1d(self.points[i].dated)


            self.points[i].tmin,self.points[i].tmax=min(self.points[i].dated),max(self.points[i].dated)
            if self.points[i].tmin < self.tmin:
                self.tmin = self.points[i].tmin
            if self.points[i].tmax > self.tmax:
                self.tmax = self.points[i].tmax
            
            self.points[i].Nt=len(self.points[i].t)
            self.N += int(len(self.points[i].t)*self.dim)
            self.d.append(self.points[i].d)

        self.d = np.array(self.d).flatten()
        self.sigmad = self.sigmad*np.ones(self.N)

    def info(self):
        print()
        print('GPS time series from network:',self.network)
        print('Number of stations:', self.Npoints)
        print('Lenght data vector:', self.N)


manifold = gpstimeseries(network=network,reduction=reduction,dim=dim,wdir=wdir,extension=ext,proj=proj,scale=scale)
manifold.load()
manifold.info()

# SVD inversion with cut-off eigenvalues
def invSVD(A,b,cond):
    try:
        U,eignv,V = lst.svd(A, full_matrices=False)
        s = np.diag(eignv)
        print(s)
        index = np.nonzero(s<cond)
        inv = lst.inv(s)
        inv[index] = 0.
        fsoln = np.dot( V.T, np.dot( inv , np.dot(U.T, b) ))
    except:
        fsoln = lst.lstsq(A,b)[0]
        #fsoln = lst.lstsq(A,b,rcond=cond)[0]
    
    return fsoln

## inversion procedure 
def consInvert(A,b,sigmad,ineq='yes',cond=1.0e-3, iter=2000,acc=1e-12, eguality=False):
    '''Solves the constrained inversion problem.
    Minimize:
    
    ||Ax-b||^2
    Subject to:
    mmin < m < mmax
    '''

    if A.shape[0] != len(b):
        raise ValueError('Incompatible dimensions for A and b')

    if ineq == 'no':
        print('ineq=no: SVD decomposition neglecting small eigenvectors inferior to {} (cond)'.format(cond))
        fsoln = invSVD(A,b,cond)
        print('SVD solution:', fsoln)

    else:
        print('ineq=yes: Iterative least-square decomposition. Prior obtained with SVD.')
        if len(indexpo>0):
          # invert first without post-seismic
          Ain = np.delete(A,indexpo,1)
          try:
              U,eignv,V = lst.svd(Ain, full_matrices=False)
              s = np.diag(eignv) 
              print('Eigenvalues:', eignv)
              index = np.nonzero(s<cond)
              inv = lst.inv(s)
              inv[index] = 0.
              mtemp = np.dot( V.T, np.dot( inv , np.dot(U.T, b) ))
          except:
              mtemp = lst.lstsq(Ain,b,rcond=cond)[0]
          print('SVD solution:', mtemp)

          # rebuild full vector
          for z in range(len(indexpo)):
            mtemp = np.insert(mtemp,indexpo[z],0)
          minit = np.copy(mtemp)
          # # initialize bounds
          mmin,mmax = -np.ones(len(minit))*np.inf, np.ones(len(minit))*np.inf 

          # We here define bounds for postseismic to be the same sign than coseismic
          # and coseismic inferior or egual to the coseimic initial 
          print('ineq=yes: Impose postseismic to be the same sign than coseismic')
          for i in range(len(indexco)):
            if (pos[i] > 0.) and (minit[int(indexco[i])]>0.):
                mmin[int(indexpofull[i])], mmax[int(indexpofull[i])] = 0, np.inf 
                mmin[int(indexco[i])], mmax[int(indexco[i])] = 0, minit[int(indexco[i])] 
            if (pos[i] > 0.) and (minit[int(indexco[i])]<0.):
                mmin[int(indexpofull[i])], mmax[int(indexpofull[i])] = -np.inf , 0
                mmin[int(indexco[i])], mmax[int(indexco[i])] = minit[int(indexco[i])], 0
          bounds=list(zip(mmin,mmax))

        else:
          minit=invSVD(A,b,cond)
          print('SVD solution:', minit)
          bounds=None
        
        def eq_cond(x, *args):
           return math.atan2(x[indexseast+1],x[[indexseast]]) - math.atan2(x[indexseas+1],x[[indexseas]])
       
        ####Objective function and derivative
        _func = lambda x: np.sum(((np.dot(A,x)-b)/sigmad)**2)
        _fprime = lambda x: 2*np.dot(A.T/sigmad, (np.dot(A,x)-b)/sigmad)
        if eguality:
            res = opt.fmin_slsqp(_func,minit,bounds=bounds,fprime=_fprime,eqcons=[eq_cond], \
                iter=iter,full_output=True,iprint=0,acc=acc)  
        else:
            res = opt.fmin_slsqp(_func,minit,bounds=bounds,fprime=_fprime, \
                iter=iter,full_output=True,iprint=0,acc=acc)  
        fsoln = res[0]
        print('Optimization:', fsoln)

    # tarantola:
    # Cm = (Gt.Cov.G)-1 --> si sigma=1 problems
    # sigma m **2 =  misfit**2 * diag([G.TG]-1)
    try:
       varx = np.linalg.inv(np.dot(A.T,A))
       # res2 = np.sum(pow((b-np.dot(A,fsoln))/sigmad,2))
       res2 = np.sum(pow((b-np.dot(A,fsoln)),2))
       scale = 1./(A.shape[0]-A.shape[1])
       # scale = 1./A.shape[0]
       sigmam = np.sqrt(scale*res2*np.diag(varx))
    except:
       sigmam = np.ones((A.shape[1]))*float('NaN')
    print('model errors:', sigmam)

    return fsoln,sigmam


### Define basis functions for plot
class pattern:
    def __init__(self,name,reduction,date):
        self.name=name
        self.reduction=reduction
        self.date=date
    
    def info(self):
        print(self.name, self.date)

def Heaviside(t):
        h=np.zeros((len(t)))
        h[t>=0]=1.0
        return h

class coseismic(pattern):
      def __init__(self,name,reduction,date):
          pattern.__init__(self,name,reduction,date)
          self.to=date

      def g(self,t):
        return (Heaviside(t-self.to))

class postseismic(pattern):
      def __init__(self,name,reduction,date,tcar=1):
          pattern.__init__(self,name,reduction,date)
          self.to=date
          self.tcar=tcar

      def g(self,t):
        t=(t-self.to)/self.tcar
        t[t<=0] = 0
        g = np.log10(1+t)
        return g

class reference(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
    def g(self,t):
        return np.ones((t.size))

class interseismic(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
        self.to=date

    def g(self,t):
        func=(t-self.to)
        return func

class sinvar(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
        self.to=date

    def g(self,t):
        func=np.zeros(t.size)
        for i in range(t.size):
            func[i]=math.sin(2*math.pi*(t[i]-self.to))
        return func

class cosvar(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
        self.to=date

    def g(self,t):
        func=np.zeros(t.size)
        for i in range(t.size):
            func[i]=math.cos(2*math.pi*(t[i]-self.to))
        return func

class sin5var(pattern):
     def __init__(self,name,reduction,date):
         pattern.__init__(self,name,reduction,date)
         self.to=date

     def g(self,t):
         func=np.zeros(t.size)
         for i in range(t.size):
             func[i]=math.sin(math.pi*(t[i]-self.to))
         return func

class cos5var(pattern):
     def __init__(self,name,reduction,date):
         pattern.__init__(self,name,reduction,date)
         self.to=date

     def g(self,t):
         func=np.zeros(t.size)
         for i in range(t.size):
             func[i]=math.cos(math.pi*(t[i]-self.to))
         return func

#datemin, datemax = int(np.min(manifold.tmin)), int(np.max(manifold.tmax))+1 
datemin, datemax= 2015, 2022

basis=[
      reference(name='reference',date=datemin,reduction='ref'),
      interseismic(name='interseismic',reduction='lin',date=datemin),
      ]

index = 2
if seasonal=='yes':
   basis.append(cosvar(name='seas. var (cos)',reduction='coswt',date=datemin))
   basis.append(sinvar(name='seas. var (sin)',reduction='sinwt',date=datemin))
   index = index +2

if arguments["--bianual"]=='yes':
     indexbi = index
     basis.append(cos5var(name='bi-anual var (cos)',reduction='cos5wt',date=datemin))
     basis.append(sin5var(name='bi-anual var (sin)',reduction='sin5wt',date=datemin))
     index = index + 2

indexco = np.zeros(len(cos))
for i in range(len(cos)):
   basis.append(coseismic(name='coseismic {}'.format(i),reduction='cos{}'.format(i),date=cos[i])),
   indexco[i] = index
   index = index + 1

indexpo = np.zeros(len(pos))
for i in range(len(pos)):
   if pos[i] > 0. :
      basis.append(postseismic(name='postseismic {}'.format(i),reduction='post{}'.format(i),date=cos[i],tcar=pos[i])),
      indexpo[i] = index
      index = index + 1

indexpo = indexpo.astype(int)
indexco = indexco.astype(int)

print()
M=len(basis)
print('Number of basis functions:', M)
for i in range((M)):
    basis[i].info()

# plot diplacements maps
nfigure = 10

colors = ['royalblue','darkgreen','darkorchid','firebrick']

for jj in range((manifold.Npoints)):
    fig = plt.figure(nfigure+1,figsize=(12,8))
    name, i, j = manifold.name[jj], manifold.x[jj], manifold.y[jj]
    print('station:', name,i,j)
    fig.suptitle('{0} {1:.3f} {2:.3f}'.format(name.decode('utf-8'),i,j))
    pt = manifold.points[jj]
    
    pt.md, pt.md_lin, pt.d_lin, md, md_lin = [], [], [], [], []
    # iter over all components
    i = 0
    for disp, sigma, name in zip(pt.d, pt.sigmad, pt.comp):
        print(name)

        # inversion model
        mdisp=np.ones((pt.Nt))*float('NaN')

        G=np.zeros((pt.Nt,M))
        for l in range((M)):
            G[:,l]=basis[l].g(pt.t)

        # Inisilize m
        m = np.zeros((M))

        t = time.time()
        # inversion
        m,sigmam = consInvert(G,disp,sigma,cond=rcond,ineq=ineq)
        
        print() 
        print('computation time:', time.time() - t)
        print()

        # forward model in original order
        mdisp = np.dot(G,m)
        pt.md.append(mdisp)
        
        tdec = np.arange(datemin, datemax, 0.01)
        G=np.zeros((len(tdec),M))
        for l in range((M)):
            G[:,l]=basis[l].g(tdec)
        model = np.dot(G[:,:M],m[:M])
        md.append(model)

        # plot date
        ax = fig.add_subplot(len(pt.comp),1,i+1)
        ax.plot(pt.t, disp, 'o', color=colors[i], ms=2)
        ax.errorbar(pt.t,disp,yerr=sigma,ecolor=colors[i],fmt='none', alpha=0.1)
        ax.plot(tdec,model,'-r')
        ax.set_ylabel('{} (mm)'.format(name))
        # ax.legend(loc='best',fontsize='small')
        ax.set_ylim([np.nanpercentile(disp,.2),np.nanpercentile(disp,99.8)])
        ax.set_xlim([datemin,datemax])
        ax.grid(True)
        i = i+1

    ax.set_xlabel('Time (Year/month/day)')
    fig.savefig('Disp_{}.eps'.format(pt.name), format='EPS')
    plt.show()

# plt.legend(loc='best')
sys.exit()
