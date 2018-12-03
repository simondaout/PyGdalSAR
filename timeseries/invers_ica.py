#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        :    Louise Maubant (Grenoble)
############################################

"""\
invers_ica.py
InSAR Time Series ICA decomposition 

Usage: invers_ica.py [--cube=<path>] [--lectfile=<path>] [--list_images=<path>] [--imref=<value>]  [--n_comp=<values>]  [--crop=<values>] [--resize=<values>] [--type=<space/time>] [--demfile=<path>] [--events=<values>] [--save_resize=<yes/no>] [--save_matrix=<yes/no>] [--smooth=<yes/no>] [--plot=<yes/no>] [--outdir=<path>]

invers_ica.py -h | --help

Options:
-h --help               Show this screen
--cube PATH             Path to displacement file [default: depl_cumul]
--lectfile PATH         Path to the lect.in file (output of invers_pixel) [default: lect.in]
--list_images PATH      Path to list images file made of 5 columns containing for each images 1) number 2) Doppler freq (not read) 3) date in YYYYMMDD format 4) numerical date 5) perpendicular baseline [default: images_retenues]
--imref VALUE           Reference image number [default: 1]
--n_comp VALUE          Number of eigenvector 
--crop VALUE            Define a region of interest for the temporal decomposition [default: 0,nlign,0,ncol]
--type space/time       Space or time decomposition [default:space]
--demfile               Optinal DEM error coefficient to be removed before ICA decomposition [default: no] 
--events VALUES         List of event dates 
--resize VALUE          Integer resize factor [default: 0]
--save_resize YES/NO    Save the resized cube [default:no]
--save_matrix YES/NO    Save the A and S matrix for each components [default:no] 
--smooth YES/NO         Smooth the vector [default:no]
--plot YES/NO           Display plots [default: yes]
--outdir PATH           Path to output dir [default: ICA]

"""


from __future__ import division
from sklearn.decomposition import FastICA, PCA
from numpy.lib.stride_tricks import as_strided

import matplotlib
matplotlib.use('TkAgg')

from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib.dates as mdates
from matplotlib.dates import date2num
from mpl_toolkits.axes_grid1 import make_axes_locatable

from datetime import datetime
import sys, os

import scipy.signal

import argparse
# docopt (command line parser)
import docopt

#====================
def date2dec(date):
    """
    date could be a datetime list
    """
    imax=len(date)
    dates_dec=np.zeros((imax))
    i = 0
    for i in range(len(date)):
        y = float(datetime.strftime(date[i],'%Y'))
        dec = float(datetime.strftime(date[i],'%j'))/365.1
        dates_dec[i] = y + dec
        i = i+1
    return dates_dec
#=====================

def resize_2d_nonan(array,factor):
    """
    Resize a 2D array by different factor on two axis sipping NaN values.
    If a new pixel contains only NaN, it will be set to NaN

    Parameters
    ----------

    array : 2D np array

    factor : int or tuple. If int x and y factor wil be the same

    Returns
    -------
    array : 2D np array scaled by factor

    @author: damo_ma
    """
    xsize, ysize = array.shape

    if isinstance(factor,int):
        factor_x = factor
        factor_y = factor
    elif isinstance(factor,tuple):
        factor_x , factor_y = factor[0], factor[1]
    else:
        raise NameError('Factor must be a tuple (x,y) or an integer')

    if not (xsize %factor_x == 0 or ysize % factor_y == 0) :
        raise NameError('Factors must be intger multiple of array shape')

    new_xsize, new_ysize = xsize/factor_x, ysize/factor_y

    new_array = np.empty([new_xsize, new_ysize])
    new_array[:] = np.nan # this saves us an assignment in the loop below

    # submatrix indexes : is the average box on the original matrix
    subrow, subcol  = np.indices((factor_x, factor_y))

     # new matrix indexs
    row, col  = np.indices((new_xsize, new_ysize))


    for i, j, ind in zip(row.reshape(-1), col.reshape(-1),range(row.size)) :
        # define the small sub_matrix as view of input matrix subset
        sub_matrix = array[subrow+i*factor_x,subcol+j*factor_y]
        # modified from any(a) and all(a) to a.any() and a.all()
        # see https://stackoverflow.com/a/10063039/1435167
        if not (np.isnan(sub_matrix)).all(): # if we haven't all NaN
            if (np.isnan(sub_matrix)).any(): # if we haven no NaN at all
                msub_matrix = np.ma.masked_array(sub_matrix,np.isnan(sub_matrix))
                (new_array.reshape(-1))[ind] = np.nanmean(msub_matrix)
            else: # if we haven some NaN
                (new_array.reshape(-1))[ind] = np.nanmean(sub_matrix)
        # the case assign NaN if we have all NaN is missing due
        # to the standard values of new_array
    return new_array
#========================
def ica_spatial(cube, ncol, nrow, ncomp, nbr_dates):
    print('Decomposition spatiale, {} components'.format(ncomp))
    N = nbr_dates
    d = cube.reshape(((nrow)*(ncol)), N)
    #c.reshape((N3570,N))
    n_comp = ncomp
    nrow = nrow
    ncol = ncol

    ica = FastICA(n_components=n_comp, max_iter=10000, tol=0.0001)

    ##### Remove les nans values
    X = d[~np.any(np.isnan(d), axis=1)]
    Snew_ = ica.fit_transform(X)
    m = ica.mixing_
    #mT_, STnew_ = ica.fit_transform(X.T)

    ### Decomposition spatiale, recrer le cube avec les bonnes dimensions
    S = np.zeros((d.shape[0], n_comp))
    S[~np.any(np.isnan(d), axis=1)] = Snew_
    S = S.reshape((nrow,ncol,n_comp))
    S[S==0] = 'nan'

    return S, m

def ica_temporel(cube, ncol, nrow, ncomp, nbr_dates):
    print('Decomposition temporelle, {} components'.format(ncomp))
    N = nbr_dates
    d = cube.reshape(((nrow)*(ncol)), N)
    n_comp = ncomp
    nrow = nrow
    ncol = ncol
    ica = FastICA(n_components=n_comp, max_iter=10000, tol=0.0001)

    ##### Remove les nans values
    X = d[~np.any(np.isnan(d), axis=1)]
    m = ica.fit_transform(X.T)
    Snew_ = ica.mixing_
    print('shape of m is {}, shape of S is {}'.format(m.shape, Snew_.shape))
    ### Decomposition spatiale, recrer le cube avec les bonnes dimensions
    S = np.zeros((d.shape[0], n_comp))
    S[~np.any(np.isnan(d), axis=1)] = Snew_
    S = S.reshape((nrow,ncol,n_comp))
    S[S==0] = 'nan'

    return S, m

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

#=========================


################ Load arguments

if __name__ == "__main__":
  
  arguments = docopt.docopt(__doc__)
  if arguments["--list_images"] ==  None:
    listim = "images_retenues"
  else:
    listim = arguments["--list_images"]
  if arguments["--lectfile"] ==  None:
    infile = "lect.in"
  else:
    infile = arguments["--lectfile"] 
  
  ncol, nlign = map(int, open(infile).readline().split(None, 2)[0:2]) 
  nb,idates,dt,base=np.loadtxt(listim, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
  N=len(dt)
  print 'Number images: ', N
  
  if arguments["--resize"] ==  None:
    resize = 0
  else:
    resize = np.int(arguments["--resize"])   
  if arguments["--n_comp"] ==  None:
    n_comp = 3
  else:
    n_comp = np.int(arguments["--n_comp"])   
  if arguments["--type"] ==  None:
      type_decomp = 'space'
  else:
      type_decomp = arguments["--type"]
  if arguments["--imref"] ==  None:
    imref = 0
  elif arguments["--imref"] < 1:
    print '--imref must be between 1 and Nimages'
  else:
    imref = int(arguments["--imref"]) - 1
  if arguments["--crop"] ==  None:
    crop = [0,nlign,0,ncol]
  else:
    crop = map(float,arguments["--crop"].replace(',',' ').split())
    docrop = 'yes'
  ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])
  print 'Compute ICA between ibeg:{} - iend:{} and jbeg:{} - jend:{} '.format(ibeg,iend,jbeg,jend)
  
  if arguments["--plot"] ==  None:
    plot = 'yes'
  else:
    plot = arguments["--plot"]
  if arguments["--cube"] ==  None:
    cubef = "depl_cumule"
  else:
    cubef = arguments["--cube"]
   
  if arguments["--save_resize"] ==  None:
    save_resize = 'no'
  else:
    save_resize = arguments["--save_resize"]
  if arguments["--save_matrix"] ==  None:
    save_matrix = 'no'
  else:
    save_matrix = arguments["--save_matrix"]
  if arguments["--smooth"] ==  None:
    smooth = 'no'
  else:
    smooth = arguments["--smooth"]

  if arguments["--events"] ==  None:
      events = []
  else:
      events = arguments["--events"].replace(',',' ').split()

  if arguments["--demfile"] ==  None:
    demf = 'no'
    dem = np.zeros((nlign,ncol))
  else:
    demf = arguments["--demfile"]
    extension = os.path.splitext(demf)[1]
    if extension == ".tif":
        ds = gdal.Open(demf, gdal.GA_ReadOnly)
        dem = ds.GetRasterBand(1).ReadAsArray()
    else:
        dem = np.fromfile(demf,dtype=np.float32).reshape((nlign,ncol))
  
  if arguments["--outdir"] ==  None:
    outdir = './ICA/'
  else:
    outdir = './' + arguments["--outdir"] + '/'

#====================================
#========= Code =====================

if not os.path.exists(outdir):
    os.makedirs(outdir)

cubei = np.fromfile(cubef,dtype=np.float32)
cube = as_strided(cubei[:nlign*ncol*N])
print 'Number of line in the cube: ', cube.shape 
kk = np.flatnonzero(cube>9990)
cube[kk] = float('NaN')
maps = cube.reshape((nlign,ncol,N))
print 'Reshape cube: ', maps.shape
# ref
cst = np.copy(maps[:,:,imref])
for l in xrange((N)):
    maps[:,:,l] = maps[:,:,l] - cst - dem*(base[l] - base[imref])
 
if resize !=0:
  maps_liss = np.zeros((np.int(nlign/resize), np.int(ncol/resize), N))
  for i in range(N):
      maps_liss[:,:,i] = resize_2d_nonan(maps[:,:,i], (resize,resize))
  ncol = np.int(ncol/resize)
  nlign = np.int(nlign/resize)
  maps = maps_liss.reshape((nlign,ncol,N))
  if save_resize == 'yes':
      fid = open(outdir+'{}_{}'.format(cube, resize),'wb')
      maps_liss.flatten().astype('float32').tofile(fid)
      fid.close()

      fid = open(outdir+'lect_ica_{}.in'.format(resize),'w')
      np.savetxt(fid, (ncol,nlign),fmt='%6i',newline='\t')
      fid.close()

if docrop == 'yes':
    ncol, nlign = iend-ibeg, jend-jbeg
    crop_maps = np.zeros((nlign,ncol,N))
    for j in xrange((N)):
        crop_maps[:,:,j] = maps[jbeg:jend,ibeg:iend,j]
    maps = np.copy(crop_maps)
    print 'Crop cube: ', maps.shape
    
# remove NaN
maps[np.isnan(maps)] = 0.

if type_decomp =='space':
   S, m = ica_spatial(maps.flatten(), ncol, nlign, n_comp, N)

if type_decomp =='time':
   print('decomposition temporelle')
   S, m = ica_temporel(maps.flatten(), ncol, nlign, n_comp, N)

fig=plt.figure(1,figsize=(14,12))
for i in range(n_comp):
    ax = fig.add_subplot(2,int(n_comp/2)+1,1+i)
    vmax = np.nanpercentile(S[:,:,i], 98)
    vmin = -vmax
    c = ax.imshow(S[:,:,i], vmax=vmax, vmin=vmin, cmap=cm.jet)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(c, cax=cax)
fig.suptitle('Spatial Eigenvectors')
fig.savefig(outdir+'decomp_{}_{}.pdf'.format(n_comp,type_decomp), format='PDF')

dates_sar = []
for i in range((N)):
    d = (datetime.strptime(str(int(idates[i])),'%Y%m%d'))
    dates_sar.append(d)

x_d = [date2num(datetime.strptime('{}'.format(d),'%Y%m%d')) for d in idates]
x_eq = [date2num(datetime.strptime('{}'.format(d),'%Y%m%d')) for d in events]

datemin, datemax = np.int(np.nanmin(dt)), np.int(np.nanmax(dt))+1
dmax = str(datemax) + '0101'
dmin = str(datemin) + '0101'
xmin = datetime.strptime('{}'.format(dmin),'%Y%m%d')
xmax = datetime.strptime('{}'.format(dmax),'%Y%m%d')
xlim=date2num(np.array([xmin,xmax]))

fig=plt.figure(2,figsize=(14,18))
fig.autofmt_xdate()
for i in range(n_comp):
    ax = fig.add_subplot(n_comp,1,i+1)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
    ax.plot(x_d, m[:,i])
    ax.scatter(x_d, m[:,i])
    for j in xrange(len(x_eq)):
        ax.axvline(x=x_eq[j],linestyle='--', color = 'r', linewidth = 1)
    ax.set_xlim(xlim)
fig.autofmt_xdate()
fig.suptitle('Temporal eigenvectors'.format(type_decomp))
fig.savefig(outdir+'decomp_{}_{}_mixing.pdf'.format(type_decomp,n_comp), format='PDF')

fid = open(outdir+'matrix_ica_{}_{}_{}'.format(cube, n_comp, type_decomp), 'wb')
S.flatten().astype('float32').tofile(fid)
fid.close()

S2 = S.reshape((nlign*ncol,n_comp))
S2 = S2[~np.any(np.isnan(S2), axis=1)]

fig=plt.figure(3,figsize=(8,4))
for i in range(n_comp):
    ax = fig.add_subplot(1,n_comp,1+i)
    xmin = np.nanpercentile(S[:,i], 1)
    xmax = np.nanpercentile(S[:,i], 99)
    plt.xlim(xmin= xmin, xmax = xmax)
    ax.hist(S2[:,i], 200,density=True)

fig.suptitle('PDFs eigenvectors'.format(type_decomp))
fig.savefig(outdir+'PDF_{}_{}.pdf'.format(n_comp, type_decomp), format='PDF')

if smooth == 'yes':
        m_smooth =np.zeros((N, n_comp))
        m_smooth = np.array(m_smooth)
        for i in range(n_comp):
            y_smooth = scipy.signal.savgol_filter(m[:,i],21,3)
            m_smooth[:,i] = y_smooth
#        np.savetxt(outdir+'vector_smooth_ica_{}_{}.txt'.format(type_decomp,n_comp),  np.column_stack((idates,m_smooth)))
        fig=plt.figure(22,figsize=(12,18))
        fig.autofmt_xdate()
        for i in range(n_comp):
            ax = fig.add_subplot(n_comp,1,i+1)
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
            ax.plot(x_d, m_smooth[:,i])
            ax.scatter(x_d, m_smooth[:,i])
            for j in xrange(len(x_eq)):
                ax.axvline(x=x_eq[j],linestyle='--', color = 'r', linewidth = 1)
            ax.set_xlim(xlim)
        fig.autofmt_xdate()
        fig.suptitle('Decomposition {}'.format(type_decomp))
        fig.savefig(outdir+'decomp_{}_{}_smooth_mixing.pdf'.format(type_decomp,n_comp), format='PDF')
else:
    m_smooth = m
np.savetxt(outdir+'vector_smooth_ica_{}_{}.txt'.format(type_decomp,n_comp),  np.column_stack((idates,m)))

if save_matrix == 'yes':
    fid = open(outdir+'maps_ica_{}_{}.r4'.format(type_decomp,n_comp), 'wb')
    S.flatten().astype('float32').tofile(fid)
    fid.close()
    np.savetxt(outdir+'vector_ica_{}_{}.txt'.format(type_decomp,n_comp),  np.column_stack((idates,m)))

###### Compute RMS 
depl_cumule_ica = np.dot(S, m.T).flatten()
models = as_strided(depl_cumule_ica.reshape((nlign, ncol, N)))
res = as_strided(maps) - as_strided(models)
res_ts = maps.reshape((nlign*ncol,N)) - depl_cumule_ica.reshape((nlign*ncol,N))

rms_pixel = np.sqrt(np.nansum(res_ts**2, axis=1)/N).reshape((nlign,ncol))
rmsd = np.sqrt(np.nansum((res_ts)**2)/((nlign)*(ncol)))
print('RMS : {}'.format(rmsd))
perc_tot = (1 - np.nansum((res)**2)/np.nansum((maps)**2))*100

S_reshape = S.reshape((nlign*ncol,n_comp))
perc = np.array(range(n_comp))
for i in range((n_comp)):
    comp = np.array([m[:,i]])
    S_def = np.array([S_reshape[:,i]]).T
    cube_def = np.dot(S_def, comp)
    cube_def = cube_def.reshape((nlign, ncol, len(m[:,0])))
    perc[i] = (1 - np.nansum((maps - cube_def)**2)/np.nansum(maps**2))*100

sum = np.sum(perc)
perc_2 = np.array(range(n_comp))
for i in range(len(perc)):
    perc_2[i] = (perc[i]/sum)*100
    print('Percentile {}th component: {}'.format(i, perc_2[i]))

fig_perc = plt.figure(500)
plt.bar((range(len(perc))), perc_2)
plt.xlabel('n comp', fontsize=18)
plt.ylabel("Percentile", fontsize=16)
fig_perc.suptitle('Tot. percentile for {} components: {:.3f}'.format(n_comp, perc_tot))
fig_perc.savefig(outdir+'percentil_{}_{}.pdf'.format(type_decomp, n_comp), format='PDF')


######### Calcule l'amplitude de chaque eigenvector et le multiplie a la matrice S
comp_ampli = np.zeros((S_reshape.shape))
comp_norm = np.zeros((m.shape))
for i in range((n_comp)):
    amplitude = ((np.nanmax(m[:,i]))-(np.nanmin(m[:,i])))
    comp_ampli[:,i] = S_reshape[:,i]*amplitude
    comp_norm[:,i] = m[:,i]/amplitude
comp_ampli = comp_ampli.reshape((nlign, ncol, n_comp))

fig_ampli=plt.figure(31,figsize=(14,12))
for i in range(n_comp):
    ax_norm_2 = fig_ampli.add_subplot(2,int(n_comp/2)+1,1+i)
    vmax = np.nanpercentile(comp_ampli[:,:,i], 98)
    vmin = -vmax
    c = ax_norm_2.imshow(comp_ampli[:,:,i], vmax=vmax, vmin=vmin, cmap=cm.jet)
    plt.setp(ax_norm_2.get_xticklabels(), visible=False)
    plt.setp(ax_norm_2.get_yticklabels(), visible=False)
    ax_norm_2.set_title('IC{}'.format(i+1),fontsize=20)
    divider = make_axes_locatable(ax_norm_2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(c, cax=cax)
fig_ampli.suptitle('Normalised spatial eigenvectors')
fig_ampli.savefig(outdir+'decomp_norm_carte_{}_{}.pdf'.format(n_comp,type_decomp), format='PDF')

# plt.show()

fig = plt.figure(figsize=(14,8))
ax_n = fig.add_subplot(111)
fig.autofmt_xdate()
for i in range(n_comp):
    #ax_n = fig.add_subplot(n_comp,1,1)
    ax_n.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
    ax_n.plot(x_d, comp_norm[:,i], label='IC{}'.format(i+1))
    for j in xrange(len(x_eq)):
        ax.axvline(x=x_eq[j],linestyle='--', color = 'r', linewidth = 1) 
    ax_n.legend()
fig.autofmt_xdate()
fig.suptitle('Normalised temporal eigenvectors')
fig.savefig(outdir+'decomp_norm_{}_{}_mixing.pdf'.format(type_decomp,n_comp), format='PDF')

for i in range(n_comp):
    fid = open(outdir+'IC{}_{}_{}.r4'.format(i+1, type_decomp, n_comp), 'wb')
    comp_ampli[:,:,i].flatten().astype('float32').tofile(fid)
    fid.close()

fig_rms = plt.figure(112)
vmax = np.nanpercentile(rms_pixel, 98)
vmin = np.nanpercentile(rms_pixel, 2)
c = plt.imshow(rms_pixel[:,:], vmax=vmax, vmin = vmin, cmap=cm.jet)
cbar = fig_rms.colorbar(c, orientation='vertical',aspect=10)
fig_rms.suptitle('Tot. RMS: {}'.format(rmsd))
fig_rms.savefig(outdir+'RMS_pixel_ica_{}_{}.pdf'.format(type_decomp, n_comp), format='PDF')

# save foward model
fid = open(outdir+'{}_ica_{}'.format(cube, n_comp), 'wb')
depl_cumule_ica.flatten().astype('float32').tofile(fid)
fid.close()

figd = plt.figure(4,figsize=(14,10))
vmax = np.nanpercentile(depl_cumule_ica, 98)
for l in range((N)):
    axd = figd.add_subplot(4,int(N/4)+1,l+1)
    axd.set_title(idates[l],fontsize=6)
    caxd = axd.imshow(maps[:,:,l],cmap=cm.jet,vmax=vmax,vmin=-vmax)
    plt.setp(axd.get_xticklabels(), visible=False)
    plt.setp(axd.get_yticklabels(), visible=False)
cbar =figd.colorbar(caxd, orientation='vertical',aspect=10)
figd.tight_layout()
figd.suptitle('Time series displacement maps')
figd.savefig(outdir+'{}_ica_{}_{}.pdf'.format(cube, type_decomp, n_comp), format='PDF')

figd = plt.figure(5,figsize=(14,10))
for l in range((N)):
    axd = figd.add_subplot(4,int(N/4)+1,l+1)
    axd.set_title(idates[l],fontsize=6)
    caxd = axd.imshow(models[:,:,l],cmap=cm.jet,vmax=vmax,vmin=-vmax)
    plt.setp(axd.get_xticklabels(), visible=False)
    plt.setp(axd.get_yticklabels(), visible=False)
cbar =figd.colorbar(caxd, orientation='vertical',aspect=10)
figd.tight_layout()
figd.suptitle('Time series displacement models')
figd.savefig(outdir+'{}_ica_{}_{}.pdf'.format(cube, type_decomp, n_comp), format='PDF')


figd = plt.figure(6,figsize=(14,10))
vmax = np.nanpercentile(res[:,:,:], 98)
for l in range((N)):
    axd = figd.add_subplot(4,int(N/4)+1,l+1)
    axd.set_title(idates[l],fontsize=6)
    caxd = axd.imshow(res[:,:,l],cmap=cm.jet,vmax=vmax,vmin=-vmax)
    plt.setp(axd.get_xticklabels(), visible=False)
    plt.setp(axd.get_yticklabels(), visible=False)
cbar =figd.colorbar(caxd, orientation='vertical',aspect=10)
figd.tight_layout()
figd.suptitle('RMS time series')
figd.savefig(outdir+'{}_ica_{}_{}_rms.pdf'.format(cube, type_decomp, n_comp), format='PDF')

fid = open(outdir+'{}_rms_ica_{}'.format(cube, n_comp), 'wb')
res.flatten().astype('float32').tofile(fid)
fid.close()

if plot == 'yes':
    plt.show()

# S_comp = S.reshape((nlign*ncol, n_comp))
# for i in range(n_comp):
#     X = np.dot(np.array([S_comp[:,i]]).T, np.array([m[:,i]]))
#     X = X.reshape((nlign, ncol, N))
#     figd = plt.figure(220+i,figsize=(14,10))
#     vmax = np.nanpercentile(X, 98)
#     for l in range((N)):
#         axd = figd.add_subplot(4,int(N/4)+1,l+1)
#         axd.set_title(date[l],fontsize=6)
#         caxd = axd.imshow(X[:,:,l],cmap=cm.jet,vmax=vmax,vmin=-vmax)
#         plt.setp(axd.get_xticklabels(), visible=False)
#         plt.setp(axd.get_yticklabels(), visible=False)
#     cbar = figd.colorbar(caxd, orientation='vertical',aspect=10)
#     figd.suptitle('Decomposition spatiale cube')
#     figd.savefig(outdir+'{}_ica_{}_{}.pdf'.format(cube, type_decomp, i), format='PDF')
#     plt.close()
