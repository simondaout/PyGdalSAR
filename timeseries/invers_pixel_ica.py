#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 09:57:15 2018

author: maubantl

invers_pixel_ica.py
-------------------
inversion en composantes independantes d'une serie temporelle insar, dans un premier temps supprime toutes les valeurs ayant un nan

Usage: invers_pixel_ica.py -c <path> -n_compo <value> -type <value> [-i <path> -p yes/no -lect <path> -resize <value>]

Options:
-h --help		 Show this screen.
-c --cube		 Cube a decomposer
-n_compo		 Nombre de composantes
-type			 type de decomposition: Spatiale ou temporelle
-i --images		 default: images_retenues
-p --plot		 show plot, default: no
-lect			 lect_file default: lect.in
-resize			 if value != 0 yes, default: no, la valeur doit être un multiple du nombre ce colonnes et de lignes
"""
from __future__ import division
from sklearn.decomposition import FastICA, PCA

import matplotlib
matplotlib.use('TkAgg')

from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib.dates as mdates
from matplotlib.dates import date2num
from datetime import datetime
import sys, os

import scipy.signal

import argparse

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
def ica_spatial(cube, ncol, nrow, ncompo, nbr_dates):
    print('Decomposition spatiale, {} composantes'.format(ncompo))
    N = nbr_dates
    d = cube.reshape(((nrow)*(ncol)), N)
    #c.reshape((N3570,N))
    n_compo = ncompo
    nrow = nrow
    ncol = ncol


    ica = FastICA(n_components=n_compo, max_iter=10000, tol=0.0001)

    ##### Remove les nans values
    X = d[~np.any(np.isnan(d), axis=1)]
    Snew_ = ica.fit_transform(X)
    m = ica.mixing_
    #mT_, STnew_ = ica.fit_transform(X.T)

    ### Decomposition spatiale, recrer le cube avec les bonnes dimensions
    S = np.zeros((d.shape[0], n_compo))
    S[~np.any(np.isnan(d), axis=1)] = Snew_
    S = S.reshape((nrow,ncol,n_compo))
    S[S==0] = 'nan'

    return S, m

def ica_temporel(cube, ncol, nrow, ncompo, nbr_dates):
    print('Decomposition temporelle, {} composantes'.format(ncompo))
    N = nbr_dates
    d = cube.reshape(((nrow)*(ncol)), N)
    n_compo = ncompo
    nrow = nrow
    ncol = ncol
    ica = FastICA(n_components=n_compo, max_iter=10000, tol=0.0001)

    ##### Remove les nans values
    X = d[~np.any(np.isnan(d), axis=1)]
    m = ica.fit_transform(X.T)
    Snew_ = ica.mixing_
    print('shape of m is {}, shape of S is {}'.format(m.shape, Snew_.shape))
    ### Decomposition spatiale, recrer le cube avec les bonnes dimensions
    S = np.zeros((d.shape[0], n_compo))
    S[~np.any(np.isnan(d), axis=1)] = Snew_
    S = S.reshape((nrow,ncol,n_compo))
    S[S==0] = 'nan'

    return S, m

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

#=========================


################ Load arguments

if __name__ == "__main__":
 # if len(sys.argv) == 0:
 #   sys.argv.append('-h')
  parser = argparse.ArgumentParser(description = "ICA decomposition")
  parser.add_argument("-c" , type=str, help="Cube, OBLIGATOIRE")
  parser.add_argument("-n_compo" , type=int, help="nombre de composantes, OBLIGATOIRE")
  parser.add_argument("-type_decom" , type=str, help="spatiale ou temporelle, OBLIGATOIRE")
  parser.add_argument("-i" , type=str, help="file of images retenues",default='images_retenues')
  parser.add_argument("-plot", type=str, help="plot yes or no",default='no')
  parser.add_argument("-lect", type=str, help="lect_file", default='lect.in')
  parser.add_argument("-resize" ,  type=int, help="resize the cube, value", default=0)
  parser.add_argument("-save_resize" ,  type=str, help="save the resized cube", default='no')
  parser.add_argument("-save_matrix" ,  type=str, help="save the A and S matrix for each components", default='no')
  parser.add_argument("-smooth" ,  type=str, help="smooth the vector", default='no')

  #if "-H" in sys.argv:
   #  print(__doc__)
   #  sys.exit(0)

  args = parser.parse_args()


#====================================
#========= Code =====================

cube = (args.c)
n_compo = (args.n_compo)
type_decom = (args.type_decom)
images = (args.i)
lect = (args.lect)
resize = (args.resize)
plot = (args.plot)
save_resize = (args.save_resize)
save_matrix = (args.save_matrix)
smooth = (args.smooth)
outdir = './ICA/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

maps=np.fromfile(cube, dtype='float32')
date, idates, dt =np.loadtxt(images, comments='#', usecols=(1,3,4),unpack=True,dtype='i,f,f')
N = len(date)
ncol, nlign = map(int, open(lect).readline().split(None, 2)[0:2])
print(ncol, nlign, N)
if resize !=0:
  maps_reshape = maps.reshape((nlign,ncol,N))
  maps_liss = np.zeros((nlign/resize, ncol/resize, N))
  for i in range(N):
      maps_liss[:,:,i] = resize_2d_nonan(maps_reshape[:,:,i], (resize,resize))
  ncol = ncol/resize
  nlign = nlign/resize
  maps = maps_liss.reshape((nlign*ncol,N))
  if save_resize == 'yes':
      fid = open(outdir+'{}_{}'.format(cube, resize),'wb')
      maps_liss.flatten().astype('float32').tofile(fid)
      fid.close()

      fid = open(outdir+'lect_ica_{}.in'.format(resize),'w')
      np.savetxt(fid, (ncol,nlign),fmt='%6i',newline='\t')
      fid.close()
maps = maps.reshape((nlign,ncol,N))
if type_decom =='spatiale':
   S, m = ica_spatial(maps, ncol, nlign, n_compo, N)

if type_decom =='temporelle':
   print('decomposition temporelle')
   S, m = ica_temporel(maps, ncol, nlign, n_compo, N)

fig=plt.figure(1,figsize=(16,12))
for i in range(n_compo):
    ax = fig.add_subplot(2,int(n_compo/2)+1,1+i)
    vmax = np.nanpercentile(S[:,:,i], 98)
    vmin = -vmax
    c = ax.imshow(S[:,:,i], vmax=vmax, vmin=vmin, cmap=cm.jet)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    cbar = fig.colorbar(c, orientation='vertical',shrink=0.5)
fig.suptitle('Decomposition {}, composantes'.format(type_decom))
fig.savefig(outdir+'decompo_{}_{}.pdf'.format(n_compo,type_decom), format='PDF')



dates_sar = []
for i in range(len(date)):
    d = (datetime.strptime(str(int(date[i])),'%Y%m%d'))
    dates_sar.append(d)

x_d = [date2num(datetime.strptime('{}'.format(d),'%Y%m%d')) for d in date]
dates_EQ = ['20170917', '20170921', '20180205', '20170501']

x_eq = [date2num(datetime.strptime('{}'.format(d),'%Y%m%d')) for d in dates_EQ]

fig=plt.figure(2,figsize=(12,18))
fig.autofmt_xdate()
for i in range(n_compo):
    ax = fig.add_subplot(n_compo,1,i+1)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
    ax.plot(x_d, m[:,i])
    ax.scatter(x_d, m[:,i])
      #ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
    ax.axvline(x=x_eq[0],linestyle='--', color = 'r', linewidth = 1)
    ax.axvspan(x_eq[3], x_eq[0], alpha=0.1, color='red')
    ax.axvline(x=x_eq[1],linestyle='--', color = 'g', linewidth = 1)
    ax.axvline(x=x_eq[2],linestyle='--', color = 'orange', linewidth = 1)
fig.autofmt_xdate()
fig.suptitle('Decomposition {}'.format(type_decom))
fig.savefig(outdir+'decompo_{}_{}_mixing.pdf'.format(type_decom,n_compo), format='PDF')


print S.shape
fid = open(outdir+'matrix_ica_{}_{}_{}'.format(cube, n_compo, type_decom), 'wb')
S.flatten().astype('float32').tofile(fid)
fid.close()

S2 = S.reshape((nlign*ncol,n_compo))
S2 = S2[~np.any(np.isnan(S2), axis=1)]

fig=plt.figure(3,figsize=(16,2))
for i in range(n_compo):
    ax = fig.add_subplot(1,n_compo,1+i)
    xmin = np.nanpercentile(S[:,i], 1)
    xmax = np.nanpercentile(S[:,i], 99)
    plt.xlim(xmin= xmin, xmax = xmax)
    ax.hist(S2[:,i], 200,density=True)

fig.suptitle('Decomposition {}, PDF'.format(type_decom))
fig.savefig(outdir+'PDF_{}_{}.pdf'.format(n_compo, type_decom), format='PDF')

if smooth == 'yes':
        m_smooth =np.zeros((N, n_compo))
        m_smooth = np.array(m_smooth)
        for i in range(n_compo):
            y_smooth = scipy.signal.savgol_filter(m[:,i],21,3)
            m_smooth[:,i] = y_smooth
#        np.savetxt(outdir+'vector_smooth_ica_{}_{}.txt'.format(type_decom,n_compo),  np.column_stack((idates,m_smooth)))
        fig=plt.figure(22,figsize=(12,18))
        fig.autofmt_xdate()
        for i in range(n_compo):
            ax = fig.add_subplot(n_compo,1,i+1)
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
            ax.plot(x_d, m_smooth[:,i])
            ax.scatter(x_d, m_smooth[:,i])
            ax.axvline(x=x_eq[0],linestyle='--', color = 'r', linewidth = 1)
            ax.axvspan(x_eq[3], x_eq[0], alpha=0.1, color='red')
            ax.axvline(x=x_eq[1],linestyle='--', color = 'g', linewidth = 1)
            ax.axvline(x=x_eq[2],linestyle='--', color = 'orange', linewidth = 1)
        fig.autofmt_xdate()
        fig.suptitle('Decomposition {}'.format(type_decom))
        fig.savefig(outdir+'decompo_{}_{}_smooth_mixing.pdf'.format(type_decom,n_compo), format='PDF')
else:
    m_smooth = m
np.savetxt(outdir+'vector_smooth_ica_{}_{}.txt'.format(type_decom,n_compo),  np.column_stack((idates,m)))


if save_matrix == 'yes':
    fid = open(outdir+'maps_ica_{}_{}.r4'.format(type_decom,n_compo), 'wb')
    S.flatten().astype('float32').tofile(fid)
    fid.close()

    np.savetxt(outdir+'vector_ica_{}_{}.txt'.format(type_decom,n_compo),  np.column_stack((idates,m)))


###### Calcul du RMS et du pourcentage d'explication du signal
depl_cumule_ica = np.dot(S, m.T).reshape((nlign, ncol, N))
data = maps.reshape((nlign, ncol, N))
cube_rms = data - depl_cumule_ica
rms = maps.reshape((nlign*ncol,N)) - depl_cumule_ica.reshape((nlign*ncol,N))

rms_pixel = np.sqrt(np.nansum(rms**2, axis=1)/N)
rms_pixel = rms_pixel.reshape((nlign,ncol))
rmsd = np.sqrt(np.nansum((rms)**2)/((nlign)*(ncol)))
print('RMS : {}'.format(rmsd))


pourcentage_tot = (1 - np.nansum((cube_rms)**2)/np.nansum((data)**2))*100
print(np.nansum(pourcentage_tot))

S_reshape = S.reshape((nlign*ncol,n_compo))
perc = np.array(range(n_compo))
for i in range((n_compo)):
    compo = np.array([m[:,i]])
    S_def = np.array([S_reshape[:,i]]).T
    cube_def = np.dot(S_def, compo)
    cube_def = cube_def.reshape((nlign, ncol, len(m[:,0])))
    pourcentage = (1 - np.nansum((data - cube_def)**2)/np.nansum((data)**2))*100
    perc[i] = pourcentage





sum = np.sum(perc)
perc_2 = np.array(range(n_compo))

for i in range(len(perc)):
    perc_2[i] = (perc[i]/sum)*100

    print('pour la composante {} le pourcentage est de {}'.format(i, perc_2[i]))

fig_perc = plt.figure(500)
plt.bar((range(len(perc))), perc_2)
plt.xlabel('n compo', fontsize=18)
plt.ylabel("percentage", fontsize=16)
fig_perc.suptitle('percentage total for {} componante = {}'.format(n_compo, pourcentage_tot))
fig_perc.savefig(outdir+'percentil_{}_{}.pdf'.format(type_decom, n_compo), format='PDF')

#plt.close('all')
#plt.show()


######### Calcule l'amplitude de chaque eigenvector et le multiplie a la matrice S
compo_ampli = np.zeros((S_reshape.shape))
compo_norm = np.zeros((m.shape))
for i in range((n_compo)):
    amplitude = ((np.nanmax(m[:,i]))-(np.nanmin(m[:,i])))
    compo_ampli[:,i] = S_reshape[:,i]*amplitude
    compo_norm[:,i] = m[:,i]/amplitude
compo_ampli = compo_ampli.reshape((nlign, ncol, n_compo))

fig_ampli=plt.figure(31)
for i in range(n_compo):
    ax_norm_2 = fig_ampli.add_subplot(2,int(n_compo/2)+1,1+i)
    vmax = np.nanpercentile(compo_ampli[:,:,i], 98)
    vmin = -vmax
    c = ax_norm_2.imshow(compo_ampli[:,:,i], vmax=vmax, vmin=vmin, cmap=cm.jet)
    plt.setp(ax_norm_2.get_xticklabels(), visible=False)
    plt.setp(ax_norm_2.get_yticklabels(), visible=False)
    ax_norm_2.set_title('IC{}'.format(i+1),fontsize=20)
    cbar = fig_ampli.colorbar(c, orientation='vertical',shrink=0.5)
fig_ampli.suptitle('Decomposition {}, composantes'.format(type_decom))
fig_ampli.savefig(outdir+'decompo_norm_carte_{}_{}.pdf'.format(n_compo,type_decom), format='PDF')


# plt.show()

#fig=plt.figure(28,figsize=(10,4))
fig = plt.figure(figsize=(4,3))
ax_n = fig.add_subplot(111)
fig.autofmt_xdate()
for i in range(n_compo):
    #ax_n = fig.add_subplot(n_compo,1,1)
    ax_n.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
    ax_n.plot(x_d, compo_norm[:,i], label='IC{}'.format(i+1))
    #ax_n.scatter(x_d, compo_norm[:,i])
      #ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
    ax_n.axvline(x=x_eq[0],linestyle='--', color = 'r', linewidth = 1)
    ax_n.axvspan(x_eq[3], x_eq[0], alpha=0.01, color='red')
    ax_n.axvline(x=x_eq[1],linestyle='--', color = 'g', linewidth = 1)
    ax_n.axvline(x=x_eq[2],linestyle='--', color = 'orange', linewidth = 1)
    ax_n.legend()
fig.autofmt_xdate()
fig.suptitle('Decomposition {}'.format(type_decom))
fig.savefig(outdir+'decompo_norm_{}_{}_mixing.pdf'.format(type_decom,n_compo), format='PDF')

for i in range(n_compo):
    fid = open(outdir+'IC{}_{}_{}.r4'.format(i+1, type_decom, n_compo), 'wb')
    compo_ampli[:,:,i].flatten().astype('float32').tofile(fid)
    fid.close()

fig_rms = plt.figure(112)
vmax = np.nanpercentile(rms_pixel, 98)
vmin = np.nanpercentile(rms_pixel, 2)
#print(vmax, vmin)
c = plt.imshow(rms_pixel[:,:], vmax=vmax, vmin = vmin, cmap=cm.jet)
cbar = fig_rms.colorbar(c, orientation='vertical',aspect=10)
fig_rms.suptitle('RMS total = {}'.format(rmsd))
fig_rms.savefig(outdir+'RMS_pixel_ica_{}_{}.pdf'.format(type_decom, n_compo), format='PDF')





depl_cumule_ica = np.dot(S.reshape((nlign*ncol, n_compo)), m.T)
depl_cumule_ica = depl_cumule_ica.reshape((nlign, ncol, N))
fid = open(outdir+'{}_ica_{}'.format(cube, n_compo), 'wb')
depl_cumule_ica.flatten().astype('float32').tofile(fid)
fid.close()


figd = plt.figure(4,figsize=(14,10))
vmax = np.nanpercentile(depl_cumule_ica[:,:,:], 98)
for l in range((N)):
    axd = figd.add_subplot(4,int(N/4)+1,l+1)
    axd.set_title(date[l],fontsize=6)
    caxd = axd.imshow(depl_cumule_ica[:,:,l],cmap=cm.jet,vmax=vmax,vmin=-vmax)
    plt.setp(axd.get_xticklabels(), visible=False)
    plt.setp(axd.get_yticklabels(), visible=False)
cbar =figd.colorbar(caxd, orientation='vertical',aspect=10)
figd.suptitle('Decomposition spatiale cube')
figd.savefig(outdir+'{}_ica_{}_{}.pdf'.format(cube, type_decom, n_compo), format='PDF')


figd = plt.figure(5,figsize=(14,10))
vmax = np.nanpercentile(cube_rms[:,:,:], 98)
for l in range((N)):
    axd = figd.add_subplot(4,int(N/4)+1,l+1)
    axd.set_title(date[l],fontsize=6)
    caxd = axd.imshow(cube_rms[:,:,l],cmap=cm.jet,vmax=vmax,vmin=-vmax)
    plt.setp(axd.get_xticklabels(), visible=False)
    plt.setp(axd.get_yticklabels(), visible=False)
cbar =figd.colorbar(caxd, orientation='vertical',aspect=10)
figd.suptitle('Decomposition spatiale cube RMS')
figd.savefig(outdir+'{}_ica_{}_{}_rms.pdf'.format(cube, type_decom, n_compo), format='PDF')

fid = open(outdir+'{}_rms_ica_{}'.format(cube, n_compo), 'wb')
cube_rms.flatten().astype('float32').tofile(fid)
fid.close()

if plot == 'yes':
    plt.show()

# S_compo = S.reshape((nlign*ncol, n_compo))
# for i in range(n_compo):
#     X = np.dot(np.array([S_compo[:,i]]).T, np.array([m[:,i]]))
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
#     figd.savefig(outdir+'{}_ica_{}_{}.pdf'.format(cube, type_decom, i), format='PDF')
#     plt.close()
