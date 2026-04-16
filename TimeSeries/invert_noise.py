#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT (ISTerre)
################################################################################

"""\
invert_noise.py
-------------
Time series inversion: solve noise/spectrum for each acquisition dates using phase/noise/spectrum of all interferograms
Solve problem : sigma(AB)^2 = sigma(A)^2 + sigma(B)^2

Usage: invert_noise.py [--basefile=<path>] --input=<path> --output=<path> [--cons=<yes/no>]
invert_phi.py  -h | --help

Options:
-h --help           Show this screen.
--basefile PATH     list images file [default: baseline.rsc]
--input PATH        Infile containing the phase for each pair of int.
--output PATH       Outfile
--cons VALUE        Add postive constrain to the inversion
"""

import scipy.linalg as lst
import scipy.optimize as opt
import numpy as np
import math

from nsbas import docopt

def consInvert(A,b,sigmad=1,ineq=[None,None], cond=1.0e-10, iter=250,acc=1e-06):
    '''Solves the constrained inversion problem.

    Minimize:

    ||Ax-b||^2

    Subject to:
    Ex >= f
    '''

    Ain = A
    bin = b

    if Ain.shape[0] != len(bin):
        raise ValueError('Incompatible dimensions for A and b')

    Ein = ineq[0]
    fin = ineq[1]

    if Ein is not None:
        if Ein.shape[0] != len(fin):
            raise ValueError('Incompatible shape for E and f')
        if Ein.shape[1] != Ain.shape[1]:
            raise ValueError('Incompatible shape for A and E')

    ####Objective function and derivative
    _func = lambda x: np.sum(((np.dot(Ain,x)-bin)/sigmad)**2)
    #_fprime = lambda x: 2*np.dot(Ain.T, (np.dot(Ain,x)-bin))
    _fprime = lambda x: 2*np.dot(Ain.T, ((np.dot(Ain,x)-bin)/sigmad**2))

    ######Inequality constraints and derivative
    if Ein is not None:
        _f_ieqcons = lambda x: np.dot(Ein,x)-fin
        _fprime_ieqcons = lambda x: Ein

    ######Actual solution of the problem
    temp = lst.lstsq(Ain,bin,cond=cond)   ####Initial guess.
    x0 = temp[0]

    if Ein is None:
        res = temp
    else:
        res = opt.fmin_slsqp(_func,x0,f_ieqcons=_f_ieqcons,fprime=_fprime, fprime_ieqcons=_fprime_ieqcons,iter=iter,full_output=True,acc=acc)
        if res[3] != 0:
            print('Exit mode %d: %s \n'%(res[3],res[4]))

    fsoln = res[0]
    return fsoln

# read arguments
arguments = docopt.docopt(__doc__)
liste_int = arguments["--input"]
outfile = arguments["--output"]
if arguments["--cons"] == None:
    cons = 'no'
else:
    cons = arguments["--cons"]
if arguments["--basefile"] ==  None:
    basefile = "baseline.rsc"
else:
    basefile=arguments["--basefile"]

#Data loading
print("int list=",liste_int)
source1=file(liste_int,'r')
date1,date2,spint=np.loadtxt(source1,comments="#",unpack=True,dtype='i,i,f')
kmax=len(date1)
print("number of interferogram: ",kmax)

print("image list=",basefile)
source2=file(basefile,'r')
im,imd=np.loadtxt(source2,comments="#",usecols=(0,4),unpack=True,dtype='i,f')
nmax=len(im)
print("number of image: ",nmax)

#build G
G=np.zeros((kmax+1,nmax))
for k in range((kmax)):
    for n in range((nmax)):
        if (date1[k]==im[n]):
            G[k,n]=1
        elif (date2[k]==im[n]):
            G[k,n]=1
# ini noise one image to 0
G[-1,0] = 1

#build d
d=np.zeros((kmax+1))
d[:kmax]=spint

print()
# Constrain
if cons=='yes':
    print("Add positive constrain to the inversion")
    f = np.zeros(nmax)
    E = np.diag(np.ones(nmax))

    print("Inversion....")
    sp = consInvert(G,d,ineq=[E,f])

else:
    print("Inversion....")
    sp = consInvert(G,d)

#temp = lst.lstsq(G,d,cond=1.0e-10)
#sp = temp[0]

## check resolution of the inversion
#Res = np.diag(np.dot(np.dot(G,np.linalg.pinv(np.dot(G.T,G))),G.T))
#print 'date , phase, res: '
#for n in xrange(len(imd)):
#    print imd[n], sp[n], Res[n]
#print

#for n in xrange(len(imd)):
#    print imd[n], sp[n]
#print

print('Saving in the output file', outfile)
# save in output file
np.savetxt(outfile, np.vstack([imd,sp]).T, fmt='%.6f')
