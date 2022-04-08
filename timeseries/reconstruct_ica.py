#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : LM
################################################################################
"""
Reconstruct the deformation signal of ICA decomposition.

invers_pixel_ica.py should run before
"""


from sklearn.decomposition import FastICA, PCA

import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib.dates as mdates
from matplotlib.dates import date2num
from datetime import datetime
import sys, os

import scipy.signal

import argparse

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description = "ICA decomposition")
  parser.add_argument("-i" , type=str, help="file of images retenues",default='images_retenues')
  parser.add_argument("-S" , type=str, help="Matrix S")
  parser.add_argument("-A" , type=str, help="Vector A")
  parser.add_argument("-plot", type=str, help="plot yes or no",default='no')
  parser.add_argument("-lect", type=str, help="lect_file", default='lect_ts.in')
  parser.add_argument("-smooth" ,  type=str, help="smooth the vector", default='yes')
  parser.add_argument("-dir" ,  type=str, help="Directory of ICA", default='./ICA/')
  parser.add_argument("-compo_choose" , nargs = '*', type=int, help="composantes contenant la deformation")

  if "-H" in sys.argv:
     print(__doc__)
     sys.exit(0)

  # if len(sys.argv) == 0:
  #     print("no argument given!")
  #     parser.print_help()
  #     sys.exit(1)

  args = parser.parse_args()


images = (args.i)
lect = (args.lect)
plot = (args.plot)
smooth = (args.smooth)
outdir = './ICA/'
compo_choose = (args.compo_choose)
print(compo_choose)
S =(args.S)
A =(args.A)


A = np.loadtxt(A)
dates = (A[:,0])
print(dates)

n_compo = A.shape[1] -1
print(n_compo)


S = np.fromfile(S, dtype='float32')
ncol, nlign = map(int, open(lect).readline().split(None, 2)[0:2])
date, idates, dt =np.loadtxt(images, comments='#', usecols=(1,3,4),unpack=True,dtype='i,f,f')
N = len(A)
dates_EQ = ['20170917', '20170921', '20180205', '20170501']
x_eq = [date2num(datetime.strptime('{}'.format(d),'%Y%m%d')) for d in dates_EQ]
x_d = [date2num(datetime.strptime('{}'.format(d),'%Y%m%d')) for d in date]


try:
    S = S.reshape((nlign*ncol,n_compo))
except ValueError:
    print('S and A should have the same compositions dimensions !!!')

if smooth =='yes':
    A_bis = np.zeros((A.shape))
    for i in range(n_compo):
        A_bis[:,i+1] = scipy.signal.savgol_filter(A[:,i+1],21,3)
    #print(A_bis-A)
    fig=plt.figure(22,figsize=(12,18))
    fig.autofmt_xdate()
    for i in range(n_compo):
        ax = fig.add_subplot(n_compo,1,i+1)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
        ax.plot(x_d, A_bis[:,i+1])
        ax.scatter(x_d, A[:,i+1],c='r')
        ax.axvline(x=x_eq[0],linestyle='--', color = 'r', linewidth = 1)
        ax.axvspan(x_eq[3], x_eq[0], alpha=0.1, color='red')
        ax.axvline(x=x_eq[1],linestyle='--', color = 'g', linewidth = 1)
        ax.axvline(x=x_eq[2],linestyle='--', color = 'orange', linewidth = 1)

    # fig.autofmt_xdate()
    # fig.suptitle('composantes smooth')
    # fig=plt.figure(23,figsize=(12,18))
    # fig.autofmt_xdate()
    # for i in range(n_compo):
    #     ax = fig.add_subplot(n_compo,1,i+1)
    #     ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
    #     ax.plot(x_d, (A_bis[:,i+1] - A[:,i+1]))
    #     ax.scatter(x_d, (A_bis[:,i+1] - A[:,i+1]))
    #     ax.axvline(x=x_eq[0],linestyle='--', color = 'r', linewidth = 1)
    #     ax.axvspan(x_eq[3], x_eq[0], alpha=0.1, color='red')
    #     ax.axvline(x=x_eq[1],linestyle='--', color = 'g', linewidth = 1)
    #     ax.axvline(x=x_eq[2],linestyle='--', color = 'orange', linewidth = 1)
    #
    # fig.autofmt_xdate()
    # fig.suptitle('composantes smooth - composante raw')
    #A = A_bis

        #fig.savefig(outdir+'decompo_{}_{}_smooth_mixing.pdf'.format(type_decom,n_compo), format='PDF')
print(nlign, ncol, len(A[:,0]))
if len(compo_choose) == 5:
    #compo = np.vstack((A[:,2], A[:,7])).T
    compo_def = np.vstack((A[:,compo_choose[0]+1], A[:,compo_choose[1]+1],A[:,compo_choose[2]+1],A[:,compo_choose[3]+1],A[:,compo_choose[4]+1])).T
    S_def = np.vstack((S[:,compo_choose[0]],S[:,compo_choose[1]],S[:,compo_choose[2]],S[:,compo_choose[3]],S[:,compo_choose[4]])).T
    print(S_def.shape)

    cube_def = np.dot(S_def, compo_def.T)
    print(np.nanmax(cube_def))
    cube_def = cube_def.reshape((nlign, ncol, len(A[:,0])))

    figd = plt.figure(4,figsize=(14,10))
    vmax = np.nanpercentile(cube_def[:,:,:], 98)
    for l in range((N)):
        axd = figd.add_subplot(4,int(N/4)+1,l+1)
        axd.set_title(date[l],fontsize=6)
        caxd = axd.imshow(cube_def[:,:,l],cmap=cm.jet,vmax=vmax,vmin=-vmax)
        plt.setp(axd.get_xticklabels(), visible=False)
        plt.setp(axd.get_yticklabels(), visible=False)
    cbar =figd.colorbar(caxd, orientation='vertical',aspect=10)
    plt.show()
if len(compo_choose) == 4:
    #compo = np.vstack((A[:,2], A[:,7])).T
    compo_def = np.vstack((A[:,compo_choose[0]+1], A[:,compo_choose[1]+1],A[:,compo_choose[2]+1],A[:,compo_choose[3]+1])).T
    S_def = np.vstack((S[:,compo_choose[0]],S[:,compo_choose[1]],S[:,compo_choose[2]],S[:,compo_choose[3 ]])).T
    print(S_def.shape)

    cube_def = np.dot(S_def, compo_def.T)
    print(np.nanmax(cube_def))
    cube_def = cube_def.reshape((nlign, ncol, len(A[:,0])))

    figd = plt.figure(4,figsize=(14,10))
    vmax = np.nanpercentile(cube_def[:,:,:], 98)
    for l in range((N)):
        axd = figd.add_subplot(4,int(N/4)+1,l+1)
        axd.set_title(date[l],fontsize=6)
        caxd = axd.imshow(cube_def[:,:,l],cmap=cm.jet,vmax=vmax,vmin=-vmax)
        plt.setp(axd.get_xticklabels(), visible=False)
        plt.setp(axd.get_yticklabels(), visible=False)
    cbar =figd.colorbar(caxd, orientation='vertical',aspect=10)
    plt.show()

if len(compo_choose) == 3:
    #compo = np.vstack((A[:,2], A[:,7])).T
    compo_def = np.vstack((A[:,compo_choose[0]+1], A[:,compo_choose[1]+1],A[:,compo_choose[2]+1])).T
    S_def = np.vstack((S[:,compo_choose[0]],S[:,compo_choose[1]],S[:,compo_choose[2]])).T
    print(S_def.shape)

    cube_def = np.dot(S_def, compo_def.T)
    print(np.nanmax(cube_def))
    cube_def = cube_def.reshape((nlign, ncol, len(A[:,0])))

    figd = plt.figure(4,figsize=(14,10))
    vmax = np.nanpercentile(cube_def[:,:,:], 98)
    for l in range((N)):
        axd = figd.add_subplot(4,int(N/4)+1,l+1)
        axd.set_title(date[l],fontsize=6)
        caxd = axd.imshow(cube_def[:,:,l],cmap=cm.jet,vmax=vmax,vmin=-vmax)
        plt.setp(axd.get_xticklabels(), visible=False)
        plt.setp(axd.get_yticklabels(), visible=False)
    cbar =figd.colorbar(caxd, orientation='vertical',aspect=10)
    plt.show()


if len(compo_choose) == 2:
    #compo = np.vstack((A[:,2], A[:,7])).T
    compo_def = np.vstack((A[:,compo_choose[0]+1], A[:,compo_choose[1]+1])).T
    S_def = np.vstack((S[:,compo_choose[0]],S[:,compo_choose[1]])).T
    print(S_def.shape)

    cube_def = np.dot(S_def, compo_def.T)
    print(np.nanmax(cube_def))
    cube_def = cube_def.reshape((nlign, ncol, len(A[:,0])))

    figd = plt.figure(4,figsize=(14,10))
    vmax = np.nanpercentile(cube_def[:,:,:], 98)
    for l in range((N)):
        axd = figd.add_subplot(4,int(N/4)+1,l+1)
        axd.set_title(date[l],fontsize=6)
        caxd = axd.imshow(cube_def[:,:,l],cmap=cm.jet,vmax=vmax,vmin=-vmax)
        plt.setp(axd.get_xticklabels(), visible=False)
        plt.setp(axd.get_yticklabels(), visible=False)
    cbar =figd.colorbar(caxd, orientation='vertical',aspect=10)
    plt.show()

    fid = open('depl_cumule_def_{}'.format(n_compo), 'wb')
    cube_def.flatten().astype('float32').tofile(fid)
    fid.close()

if len(compo_choose) == 1:
    #compo = np.vstack((A[:,2], A[:,7])).T
    compo_def = np.array([A[:,compo_choose[0]+1]])
    S_def = np.array([S[:,compo_choose[0]]]).T
    print(S_def.shape)

    cube_def = np.dot(S_def, compo_def)
    cube_def = cube_def.reshape((nlign, ncol, len(A[:,0])))

    figd = plt.figure(4,figsize=(14,10))
    vmax = np.nanpercentile(cube_def[:,:,:], 98)
    for l in range((N)):
        axd = figd.add_subplot(4,int(90/4)+1,l+1)
        axd.set_title(date[l],fontsize=6)
        caxd = axd.imshow(cube_def[:,:,l],cmap=cm.jet,vmax=vmax,vmin=-vmax)
        plt.setp(axd.get_xticklabels(), visible=False)
        plt.setp(axd.get_yticklabels(), visible=False)
    cbar =figd.colorbar(caxd, orientation='vertical',aspect=10)
    #plt.show()

fid = open('depl_cumule_def_{}'.format(n_compo), 'wb')
cube_def.flatten().astype('float32').tofile(fid)
fid.close()

# vmax = np.nanmax(compo_def)
# vmin = np.nanmin(compo_def)
# ampl = vmax - vmin
# print(ampl)
# defo = S_def.reshape((nlign, ncol))*-ampl*-0.44133
# figd = plt.figure(5,figsize=(14,10))
# axd = figd.add_subplot(1,1,1)
# caxd = axd.imshow(defo[:,:], vmax=(np.nanpercentile(defo, 98)), vmin=-(np.nanpercentile(defo, 98)), cmap=cm.jet)
# cbar =figd.colorbar(caxd, orientation='vertical',aspect=10)
# plt.show()
