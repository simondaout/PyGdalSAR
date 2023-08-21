#!/usr/bin/env python3
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
invers_pca.py
InSAR Time Series PCA decomposition

Usage: invers_pca.py [--cube=<path>] [--lectfile=<path>] [--list_images=<path>] [--imref=<value>]  [--n_comp=<values>]  [--crop=<values>] [--type=<space/time>] [--demfile=<path>] [--events=<values>] [--save_matrix=<yes/no>] [--plot=<yes/no>] [--outdir=<path>]

invers_ica.py -h | --help

Options:
-h --help               Show this screen
--cube PATH             Path to displacement file [default: depl_cumul]
--lectfile PATH         Path to the lect.in file (output of invers_pixel) [default: lect.in]
--list_images PATH      Path to list images file made of 5 columns containing for each images 1) number 2) Doppler freq (not read) 3) date in YYYYMMDD format 4) numerical date 5) perpendicular baseline [default: images_retenues]
--imref VALUE           Reference image number [default: 1]
--n_comp VALUE          Number of eigenvector [default: 10]
--crop VALUE            Define a region of interest for the temporal decomposition [default: 0,nlign,0,ncol]
--type space/time       Space or time decomposition [default:space]
--demfile               Optional DEM error coefficient to be removed before ICA decomposition [default: no]
--events VALUES         List of event dates
--save_matrix YES/NO    Save the A and S matrix for each components [default:no]
--plot YES/NO           Display plots [default: yes]
--outdir PATH           Path to output dir [default: ICA]

"""

print()
print()
print('Please cite:')
print('Maubant, L., Pathier, E., Daout, S., Radiguet, M., Doin, M. P., Kazachkina, E., ... & Walpersdorf, A. (2020). Independent component analysis and parametric approach for source separation in InSAR time series at regional scale: application to the 2017–2018 Slow Slip Event in Guerrero (Mexico). Journal of Geophysical Research: Solid Earth, 125(3), e2019JB018187.')
print()
print()

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
from nsbas import docopt

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


def pca_spatial(cube, ncol, nrow, n_comp, nbr_dates):
    print('Decomposition spatiale, {} components'.format(n_comp))
    N = nbr_dates
    d = cube.reshape(((nrow)*(ncol)), N)
    #c.reshape((N3570,N))
    nrow = nrow
    ncol = ncol
    pca = PCA(n_components=n_comp)

    ##### Remove les nans values
    X = d[~np.any(np.isnan(d), axis=1)]
    Snew_ = pca.fit_transform(X)
    m = pca.components_



    ### Decomposition spatiale, recrer le cube avec les bonnes dimensions
    S = np.zeros((d.shape[0], n_comp))
    S[~np.any(np.isnan(d), axis=1)] = Snew_
    S = S.reshape((nrow,ncol,n_comp))
    S[S==0] = 'nan'

    var = pca.explained_variance_ratio_*100
    return S, m, var



def pca_temporel(cube, ncol, nrow, n_comp, nbr_dates):
    print('Decomposition temporelle, {} components'.format(n_comp))
    N = nbr_dates
    d = cube.reshape(((nrow)*(ncol)), N)
    #c.reshape((N3570,N))
    nrow = nrow
    ncol = ncol
    pca = PCA(n_components=n_comp)

    ##### Remove les nans values
    X = d[~np.any(np.isnan(d), axis=1)]
    m = pca.fit_transform(X.T)
    Snew_ = pca.components_

    ### Decomposition spatiale, recrer le cube avec les bonnes dimensions
    S = np.zeros((d.shape[0], n_comp))
    S[~np.any(np.isnan(d), axis=1)] = Snew_.T
    S = S.reshape((nrow,ncol,n_comp))
    S[S==0] = 'nan'
    var = pca.explained_variance_ratio_*100

    return S, m, var
#==================================
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
  print('Number images: ', N)

  if arguments["--n_comp"] ==  None:
    n_comp = 10
  else:
    n_comp = np.int(arguments["--n_comp"])
  if arguments["--type"] ==  None:
      type_decomp = 'space'
  else:
      type_decomp = arguments["--type"]
  if arguments["--imref"] ==  None:
    imref = 0
  elif arguments["--imref"] < 1:
    print('--imref must be between 1 and Nimages')
  else:
    imref = int(arguments["--imref"]) - 1
  if arguments["--crop"] ==  None:
    crop = [0,nlign,0,ncol]
    docrop = 'no'
  else:
    crop = map(float,arguments["--crop"].replace(',',' ').split())
    docrop = 'yes'
  ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])
  print('Compute PCA between ibeg:{} - iend:{} and jbeg:{} - jend:{} '.format(ibeg,iend,jbeg,jend))

  if arguments["--plot"] ==  None:
    plot = 'yes'
  else:
    plot = arguments["--plot"]
  if arguments["--cube"] ==  None:
    cubef = "depl_cumule"
    print("Carefull: Data should be detrend")
  else:
    cubef = arguments["--cube"]

  if arguments["--save_matrix"] ==  None:
    save_matrix = 'no'
  else:
    save_matrix = arguments["--save_matrix"]

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
    outdir = './PCA/'
  else:
    outdir = './' + arguments["--outdir"] + '/'


if not os.path.exists(outdir):
    os.makedirs(outdir)

cubei = np.fromfile(cubef,dtype=np.float32)
cube = as_strided(cubei[:nlign*ncol*N])
print('Number of line in the cube: ', cube.shape)
kk = np.flatnonzero(cube>9990)
cube[kk] = float('NaN')
maps = cube.reshape((nlign,ncol,N))
print('Reshape cube: ', maps.shape)
# ref
cst = np.copy(maps[:,:,imref])
for l in range((N)):
    maps[:,:,l] = maps[:,:,l] - cst - dem*(base[l] - base[imref])


if docrop == 'yes':
    ncol, nlign = iend-ibeg, jend-jbeg
    crop_maps = np.zeros((nlign,ncol,N))
    for j in range((N)):
        crop_maps[:,:,j] = maps[jbeg:jend,ibeg:iend,j]
    maps = np.copy(crop_maps)
    print('Crop cube: ', maps.shape)


if type_decomp =='space':
   S, m, var = pca_spatial(maps.flatten(), ncol, nlign, n_comp, N)

if type_decomp =='time':
   print('decomposition temporelle')
   S, m, var = pca_temporel(maps.flatten(), ncol, nlign, n_comp, N)

N_Var = len(var)+1
fig_var=plt.figure(101,figsize=(14,12))
cumulative = np.cumsum(var)
plt.bar(range(1,N_Var), var)
plt.plot(range(1,N_Var), cumulative, 'red')
plt.scatter(range(1,N_Var), cumulative)
plt.yticks(np.arange(0, 100, 10))
plt.grid()
fig_var.suptitle('Variance')
fig_var.savefig(outdir+'variance_{}_{}.pdf'.format(n_comp,type_decomp), format='PDF')
np.savetxt(outdir+'variance_{}_{}.txt'.format(n_comp,type_decomp), var)


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

print(m.T.shape)

if type_decomp == 'space':
    m = m.T

fig=plt.figure(2,figsize=(14,18))
fig.autofmt_xdate()
for i in range(n_comp):
    ax = fig.add_subplot(n_comp,1,i+1)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
    ax.plot(x_d, m[:,i])
    ax.scatter(x_d, m[:,i])
    for j in range(len(x_eq)):
        ax.axvline(x=x_eq[j],linestyle='--', color = 'r', linewidth = 1)
    ax.set_xlim(xlim)
fig.autofmt_xdate()
fig.suptitle('Temporal eigenvectors'.format(type_decomp))
fig.savefig(outdir+'decomp_{}_{}_mixing.pdf'.format(type_decomp,n_comp), format='PDF')

fid = open(outdir+'matrix_pca_{}_{}_{}'.format(cubef, n_comp, type_decomp), 'wb')
S.flatten().astype('float32').tofile(fid)
fid.close()

S2 = S.reshape((nlign*ncol,n_comp))
S2 = S2[~np.any(np.isnan(S2), axis=1)]

fig=plt.figure(3,figsize=(12,5))
for i in range(n_comp):
    ax = fig.add_subplot(1,n_comp,1+i)
    xmin = np.nanpercentile(S[:,i], 1)
    xmax = np.nanpercentile(S[:,i], 99)
    plt.xlim(xmin= xmin, xmax = xmax)
    ax.hist(S2[:,i], 200,density=True)

fig.suptitle('PDFs eigenvectors'.format(type_decomp))
fig.savefig(outdir+'PDF_{}_{}.pdf'.format(n_comp, type_decomp), format='PDF')


np.savetxt(outdir+'vector_pca_{}_{}.txt'.format(type_decomp,n_comp),  np.column_stack((idates,m)))

if save_matrix == 'yes':
    fid = open(outdir+'maps_pca_{}_{}.r4'.format(type_decomp,n_comp), 'wb')
    S.flatten().astype('float32').tofile(fid)
    fid.close()
    np.savetxt(outdir+'vector_pca_{}_{}.txt'.format(type_decomp,n_comp),  np.column_stack((idates,m)))


depl_cumule_pca = np.dot(S, m.T).flatten()
models = as_strided(depl_cumule_pca.reshape((nlign, ncol, N)))
res = as_strided(maps) - as_strided(models)
res_ts = maps.reshape((nlign*ncol,N)) - depl_cumule_pca.reshape((nlign*ncol,N))

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
    for j in range(len(x_eq)):
        ax.axvline(x=x_eq[j],linestyle='--', color = 'r', linewidth = 1)
    ax_n.legend()
fig.autofmt_xdate()
fig.suptitle('Normalised temporal eigenvectors')
fig.savefig(outdir+'decomp_norm_{}_{}_mixing.pdf'.format(type_decomp,n_comp), format='PDF')

# for i in range(n_comp):
#     fid = open(outdir+'PC{}_{}_{}.r4'.format(i+1, type_decomp, n_comp), 'wb')
#     comp_ampli[:,:,i].flatten().astype('float32').tofile(fid)
#     fid.close()

fig_rms = plt.figure(112)
vmax = np.nanpercentile(rms_pixel, 98)
vmin = np.nanpercentile(rms_pixel, 2)
c = plt.imshow(rms_pixel[:,:], vmax=vmax, vmin = vmin, cmap=cm.jet)
cbar = fig_rms.colorbar(c, orientation='vertical',aspect=10)
fig_rms.suptitle('Tot. RMS: {}'.format(rmsd))
fig_rms.savefig(outdir+'RMS_pixel_pca_{}_{}.pdf'.format(type_decomp, n_comp), format='PDF')

# save foward model
fid = open(outdir+'{}_pca_{}'.format(cubef, n_comp), 'wb')
depl_cumule_pca.flatten().astype('float32').tofile(fid)
fid.close()

figd = plt.figure(4,figsize=(14,10))
vmax = np.nanpercentile(depl_cumule_pca, 98)
for l in range((N)):
    axd = figd.add_subplot(4,int(N/4)+1,l+1)
    axd.set_title(idates[l],fontsize=6)
    caxd = axd.imshow(maps[:,:,l],cmap=cm.jet,vmax=vmax,vmin=-vmax)
    plt.setp(axd.get_xticklabels(), visible=False)
    plt.setp(axd.get_yticklabels(), visible=False)
cbar =figd.colorbar(caxd, orientation='vertical',aspect=10)
#figd.tight_layout()
figd.suptitle('Time series displacement maps')
figd.savefig(outdir+'{}_pca_{}_{}.pdf'.format(cubef, type_decomp, n_comp), format='PDF')

figd = plt.figure(5,figsize=(14,10))
for l in range((N)):
    axd = figd.add_subplot(4,int(N/4)+1,l+1)
    axd.set_title(idates[l],fontsize=6)
    caxd = axd.imshow(models[:,:,l],cmap=cm.jet,vmax=vmax,vmin=-vmax)
    plt.setp(axd.get_xticklabels(), visible=False)
    plt.setp(axd.get_yticklabels(), visible=False)
cbar =figd.colorbar(caxd, orientation='vertical',aspect=10)
#figd.tight_layout()
figd.suptitle('Time series displacement models')
figd.savefig(outdir+'{}_pca_{}_{}.pdf'.format(cubef, type_decomp, n_comp), format='PDF')


figd = plt.figure(6,figsize=(14,10))
vmax = np.nanpercentile(res[:,:,:], 98)
for l in range((N)):
    axd = figd.add_subplot(4,int(N/4)+1,l+1)
    axd.set_title(idates[l],fontsize=6)
    caxd = axd.imshow(res[:,:,l],cmap=cm.jet,vmax=vmax,vmin=-vmax)
    plt.setp(axd.get_xticklabels(), visible=False)
    plt.setp(axd.get_yticklabels(), visible=False)
cbar =figd.colorbar(caxd, orientation='vertical',aspect=10)
#figd.tight_layout()
figd.suptitle('RMS time series')
figd.savefig(outdir+'{}_pca_{}_{}_rms.pdf'.format(cubef, type_decomp, n_comp), format='PDF')

fid = open(outdir+'{}_rms_pca_{}'.format(cubef, n_comp), 'wb')
res.flatten().astype('float32').tofile(fid)
fid.close()

if plot == 'yes':
    plt.show()
