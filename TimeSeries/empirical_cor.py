#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Simon DAOUT (CRPG, Nancy)
############################################

import numpy as np
from numpy.lib.stride_tricks import as_strided
import scipy as sp
import scipy.optimize as opt
import scipy.linalg as lst
from os import path, environ, getcwd
import matplotlib
import matplotlib.cm as cm
import matplotlib.dates as mdates

def linear_inv(G, data, sigma):
      'Iterative linear inversion'

      x0 = lst.lstsq(G,data)[0]
      _func = lambda x: np.sum(((np.dot(G,x)-data)/sigma)**2)
      _fprime = lambda x: 2*np.dot(G.T/sigma, (np.dot(G,x)-data)/sigma)
      pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0,acc=1.e-9)[0]

      return pars

def estim_ramp(los,los_clean,topo_clean,az,rg,order,sigma,nfit,ivar,l,ax_dphi,new_lines,new_cols):
      'Ramp/Topo estimation and correction. Estimation is performed on sliding median'

      # initialize topo
      topo = np.zeros((new_lines,new_cols))
      ramp = np.zeros((new_lines,new_cols))      

      # y: range, x: azimuth
      if arguments["--topofile"] is None:
          topobins = topo_clean
          rgbins, azbins = rg, az
          data = np.copy(los_clean)

      else:
          # lets try to digitize to improve the fit
          # digitize data in bins, compute median and std
          bins = np.arange(mintopo,maxtopo,abs(maxtopo-mintopo)/500.)
          inds = np.digitize(topo_clean,bins)
          topobins = []
          losbins = []
          losstd = []
          azbins, rgbins = [], []
          los_clean2, topo_clean2, az_clean2, rg_clean2, sigma_clean2 = [], [], [], [], []
          for j in range(len(bins)-1):
                  uu = np.flatnonzero(inds == j)
                  if len(uu)>200:
                      topobins.append(bins[j] + (bins[j+1] - bins[j])/2.)

                      # do a small clean within the bin
                      indice = np.flatnonzero(np.logical_and(los_clean[uu]>np.percentile(\
                          los_clean[uu],100-float(arguments["--perc_los"])),los_clean[uu]<np.percentile(los_clean[uu],float(arguments["--perc_los"]))))

                      losstd.append(np.nanstd(los_clean[uu][indice]))
                      losbins.append(np.nanmedian(los_clean[uu][indice]))
                      azbins.append(np.nanmedian(az[uu][indice]))
                      rgbins.append(np.nanmedian(rg[uu][indice]))

                      # remove outliers from data
                      los_clean2.append(los_clean[uu][indice])
                      az_clean2.append(az[uu][indice])
                      rg_clean2.append(rg[uu][indice])
                      topo_clean2.append(topo_clean[uu][indice])
                      # new sigma is clean sigma time the standard deviation of the los within the bin
                      sigma_clean2.append(sigma[uu][indice]*np.nanstd(los_clean[uu][indice]))

          topobins = np.array(topobins)
          rgbins, azbins = np.array(rgbins),np.array(azbins)

          # data and model cleaned a second time by sliding median
          los_clean = np.concatenate(los_clean2)
          az, rg = np.concatenate(az_clean2), np.concatenate(rg_clean2)
          topo_clean = np.concatenate(topo_clean2)
          sigma = np.concatenate(sigma_clean2)
          del los_clean2, topo_clean2, az_clean2, rg_clean2, sigma_clean2 
          
          if order == 0 and ivar == 0:
              data = np.array(losbins)
              sigma = np.array(losstd)

          else:
              data = np.array(los_clean)

          if len(data) < 10:
            logger.critical('Too small area for empirical phase/topo relationship. Re-defined crop values Exit!')
            sys.exit()
    
      if order==0:

        if arguments["--topofile"] is None:

            a = 0.
            ramp = np.zeros((new_lines,new_cols))
            res = los 

        else:

            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),2))
                G[:,0] = 1
                G[:,1] = topobins

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]
                logger.info('Remove ref frame %f + %f z for date: %i'%(a,b,idates[l]))

                # plot phase/elev
                funct = a
                funcbins = a
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,b*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),2))
                G[:,0] = 1
                G[:,1] = elevi

                res = los - np.dot(G,pars)
                topo = np.dot(G,pars).reshape(new_lines,new_cols)


            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),3))
                G[:,0] = 1
                G[:,1] = topobins
                G[:,2] = topobins**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c=pars[2]
                print ('Remove ref frame %f + %f z + %f z**2 for date: %i'%(a,b,c,idates[l]))

                # plot phase/elev
                funct = a
                funcbins = a
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,b*x+c*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),3))
                G[:,0] = 1
                G[:,1] = elevi
                G[:,2] = elevi**2

                res = los - np.dot(G,pars)
                topo = np.dot(G,pars).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),3))
                G[:,0] = 1
                G[:,1] = topo_clean
                G[:,2] = az*topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]
                print ('Remove ref frame %f + %f z + %f az*z for date: %i'%(a,b,c,idates[l]))

                # plot phase/elev
                funct = a + c*topo_clean*az
                funcbins = a + c*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,b*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),3))
                G[:,0] = 1
                G[:,1] = elevi
                G[:,2] = elevi
                for i in range(nlines):
                    G[i*new_cols:(i+1)*new_cols,2] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)

                topo = np.dot(G,pars).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),4))
                G[:,0] = 1
                G[:,1] = az*topo_clean
                G[:,2] = topo_clean
                G[:,3] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
                print ('Remove ref frame %f + %f az*z + %f z + %f z**2 for date: %i'%(a,b,c,d,idates[l]))

                # plot phase/elev
                funct = a + b*topo_clean*az
                funcbins = a + b*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,c*x+d*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                G[:,0] = 1
                G[:,1] = elevi
                G[:,2] = elevi
                G[:,3] = elevi**2
                for i in range(nlines):
                    G[i*new_cols:(i+1)*new_cols,1] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                rms = np.sqrt(np.nanmean(res**2))
                logger.info('RMS dates %i: %f'%(idates[l], rms))
                topo = np.dot(G,pars).reshape(new_lines,new_cols)

      elif order==1: # Remove a range ramp ay+b for each maps (y = col)

        if arguments["--topofile"] is None:
            G=np.zeros((len(data),2))
            G[:,0] = rg
            G[:,1] = 1

            # ramp inversion
            x0 = lst.lstsq(G,data)[0]
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]
            print ('Remove ramp %f r + %f for date: %i'%(a,b,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),2))
            for i in range(nlines):
                G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
            G[:,1] = 1

            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),3))
                G[:,0] = rg
                G[:,1] = 1
                G[:,2] = topo_clean

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]
                print ('Remove ramp %f r + %f + %f z for date: %i'%(a,b,c,idates[l]))

                # plot phase/elev
                funct = a*rg + b
                funcbins = a*rgbins + b
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,c*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),3))
                for i in range(nlines):
                    G[i*ncol:(i+1)*ncol,0] = np.arange((ncol)) - jbeg_emp
                G[:,1] = 1
                G[:,2] = elevi

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),4))
                G[:,0] = rg
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d=pars[3]
                print ('Remove ramp %f r + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,idates[l]))

                # plot phase/elev
                funct = a*rg+ b
                funcbins = a*rgbins+ b
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,c*x+d*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                for i in range(nlines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                G[:,1] = 1
                G[:,2] = elevi
                G[:,3] = elevi**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),4))
                G[:,0] = rg
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
                print ('Remove ramp %f r + %f + %f z + %f z*az for date: %i'%(a,b,c,d,idates[l]))

                # plot phase/elev
                funct = a*rg+ b + d*topo_clean*az
                funcbins = a*rgbins+ b + d*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,c*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                G[:,1] = 1
                G[:,2] = elevi
                G[:,3] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,3] *= (i - ibeg_emp)


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),5))
                G[:,0] = rg
                G[:,1] = 1
                G[:,2] = topo_clean*az
                G[:,3] = topo_clean
                G[:,4] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
                print ('Remove ramp %f r + %f +  %f z*az + %f z + %f z**2 for date: %i'%(a,b,c,d,e,idates[l]))

                # plot phase/elev
                funct = a*rg+ b + c*topo_clean*az
                funcbins = a*rgbins+ b + c*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,d*x+e*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                G[:,1] = 1
                G[:,2] = elevi
                G[:,3] = elevi
                G[:,4] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)


      elif order==2: # Remove a azimutal ramp ax+b for each maps (x is lign)
        if arguments["--topofile"] is None:
            G=np.zeros((len(data),2))
            G[:,0] = az
            G[:,1] = 1

            # ramp inversion
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]
            print ('Remove ramp %f az + %f for date: %i'%(a,b,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),2))
            for i in range(new_lines):
                G[i*new_cols:(i+1)*new_cols,0] =(i - ibeg_emp)
            G[:,1] = 1


            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),3))
                G[:,0] = az
                G[:,1] = 1
                G[:,2] = topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]
                print ('Remove ramp %f az + %f + %f z for date: %i'%(a,b,c,idates[l]))

                # plot phase/elev
                funct = a*az + b
                funcbins = a*azbins + b
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,c*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),3))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] =(i - ibeg_emp)
                G[:,1] = 1
                G[:,2] = elevi


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),4))
                G[:,0] = az
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
                print ('Remove ramp %f az + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,idates[l]))

                # plot phase/elev
                funct = a*az + b
                funcbins = a*azbins + b
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,c*x + d*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] =(i - ibeg_emp)
                G[:,1] = 1
                G[:,2] = elevi
                G[:,3] = elevi**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),4))
                G[:,0] = az
                G[:,1] = 1
                G[:,2] = topo_clean
                G[:,3] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
                print ('Remove ramp %f az + %f + %f z + %f z*az for date: %i'%(a,b,c,d,idates[l]))

                # plot phase/elev
                funct = a*az + b + d*topo_clean*az
                funcbins = a*azbins + b + d*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,c*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                G[:,1] = 1
                G[:,2] = elevi
                G[:,3] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,3] *= (i - ibeg_emp)


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),5))
                G[:,0] = az
                G[:,1] = 1
                G[:,2] = topo_clean*az
                G[:,3] = topo_clean
                G[:,4] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
                print ('Remove ramp %f az + %f + %f z*az + %f z + %f z**2 for date: %i'%(a,b,c,d,e,idates[l]))

                # plot phase/elev
                funct = a*az + b + c*topo_clean*az
                funcbins = a*azbins + b + c*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,d*x+e*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                G[:,1] = 1
                G[:,2] = elevi
                G[:,3] = elevi
                G[:,4] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)

      elif order==3: # Remove a ramp ay+bx+c for each maps
        if arguments["--topofile"] is None:
            G=np.zeros((len(data),3))
            G[:,0] = rg
            G[:,1] = az
            G[:,2] = 1

            # ramp inversion
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]; c = pars[2]
            print ('Remove ramp %f r  + %f az + %f for date: %i'%(a,b,c,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),3))
            for i in range(new_lines):
                G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                G[i*new_cols:(i+1)*new_cols,1] = (i - ibeg_emp)
            G[:,2] = 1

            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),4))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = 1
                G[:,3] = topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
                print ('Remove ramp %f r  + %f az + %f + %f z for date: %i'%(a,b,c,d,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c
                funcbins = a*rgbins+ b*azbins + c
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,d*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),4))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] =(i - ibeg_emp)
                G[:,2] = 1
                G[:,3] = elevi


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),5))
                G[:-1,0] = rg
                G[:-1,1] = az
                G[:,2] = 1
                G[:-1,3] = topo_clean
                G[:-1,4] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
                print ('Remove ramp %f r  + %f az + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,e,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c
                funcbins = a*rgbins+ b*azbins + c
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,d*x+e*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] =(i - ibeg_emp)
                G[:,2] = 1
                G[:,3] = elevi
                G[:,4] = elevi**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),5))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = 1
                G[:,3] = topo_clean
                G[:,4] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e=pars[4]
                print ('Remove ramp %f r  + %f az + %f + %f z +  %f z*az for date: %i'%(a,b,c,d,e,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c + e*topo_clean*az
                funcbins = a*rgbins+ b*azbins + c + e*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,d*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                G[:,2] = 1
                G[:,3] = elevi
                G[:,4] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] =(i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,4] *= (i - ibeg_emp)


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),6))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = 1
                G[:,3] = topo_clean*az
                G[:,4] = topo_clean
                G[:,5] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e=pars[4]; f=pars[5]
                print ('Remove ramp %f r  + %f az + %f +  %f z*az + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c + d*topo_clean*az
                funcbins = a*rgbins+ b*azbins + c + d*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.1, alpha=0.01, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x+f*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,2] = 1
                G[:,3] = elevi
                G[:,4] = elevi
                G[:,5] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] =(i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,3] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)

      elif order==4:
        if arguments["--topofile"] is None:
            G=np.zeros((len(data),4))
            G[:,0] = rg
            G[:,1] = az
            G[:,2] = rg*az
            G[:,3] = 1

            # ramp inversion
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
            print ('Remove ramp %f r %f az  + %f r*az + %f for date: %i'%(a,b,c,d,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),4))
            for i in range(new_lines):
                G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                G[i*new_cols:(i+1)*new_cols,2] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
            G[:,3] = 1

            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),5))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = rg*az
                G[:,3] = 1
                G[:,4] = topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]

                print ('Remove ramp %f r, %f az  + %f r*az + %f + %f z for date: %i'%(a,b,c,d,e,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c*az*rg+ d
                funcbins = a*rgbins+ b*azbins + c*azbins*rgbins+ d
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
                G[:,3] = 1
                G[:,4] = elevi


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),6))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = rg*az
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]

                print ('Remove ramp %f r, %f az  + %f r*az + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c*az*rg+ d
                funcbins = a*rgbins+ b*azbins + c*azbins*rgbins+ d
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x+f*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),6))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = rg*az
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]

                print ('Remove ramp %f r, %f az  + %f r*az + %f + %f z + %f az*z for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c*az*rg+ d + f*topo_clean*az
                funcbins = a*rgbins+ b*azbins + c*azbins*rgbins+ d + f*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
                    G[i*new_cols:(i+1)*new_cols,5] *=  (i - ibeg_emp)


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),6))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = rg*az
                G[:,3] = 1
                G[:,4] = topo_clean*az
                G[:,5] = topo_clean
                G[:,6] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]

                print ('Remove ramp %f r, %f az  + %f r*az + %f + + %f az*z +  %f z + %f z**2  for date: %i'%(a,b,c,d,e,f,g,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c*az*rg+ d + e*topo_clean*az
                funcbins = a*rgbins+ b*azbins + c*azbins*rgbins+ d + e*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,f*x+g*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi
                G[:,6] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
                    G[i*new_cols:(i+1)*new_cols,4] *=  (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)

      elif order==5:

        if arguments["--topofile"] is None:
            G=np.zeros((len(data),4))
            G[:,0] = rg**2
            G[:,1] = rg
            G[:,2] = az
            G[:,3] = 1

            # ramp inversion
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
            print ('Remove ramp %f r**2 + %f r  + %f az + %f for date: %i'%(a,b,c,d,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),4))
            for i in range(new_lines):
                G[i*new_cols:(i+1)*new_cols,0] = (np.arange((new_cols)) - jbeg_emp)**2
                G[i*new_cols:(i+1)*new_cols,1] = np.arange((new_cols)) - jbeg_emp
                G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
            G[:,3] = 1


            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:

            if (ivar==0 and nfit==0):

                G=np.zeros((len(data),5))
                G[:,0] = rg**2
                G[:,1] = rg
                G[:,2] = az
                G[:,3] = 1
                G[:,4] = topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
                print ('Remove ramp %f r**2 + %f r  + %f az + %f + %f z for date: %i'%(a,b,c,d,e,idates[l]))

                # plot phase/elev
                funct = a*rg**2 + b*rg+ c*az + d
                funcbins = a*rgbins**2 + b*rgbins+ c*azbins + d
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = np.arange((new_cols)) -  jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                G[:,3] = 1
                G[:,4] = elevi


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),6))
                G[:,0] = rg**2
                G[:,1] = rg
                G[:,2] = az
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
                print ('Remove ramp %f r**2 + %f r  + %f az + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funct = a*rg**2 + b*rg+ c*az + d
                funcbins = a*rgbins**2 + b*rgbins+ c*azbins + d
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x+f*x**2,'-r', lw =4.)


                # build total G matrix
                G=np.zeros((len(los),6))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = np.arange((new_cols)) -  jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi**2

                res = los - np.dot(G,pars)

                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==0):

                G=np.zeros((len(data),6))
                G[:,0] = rg**2
                G[:,1] = rg
                G[:,2] = az
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
                print ('Remove ramp %f r**2 + %f r  + %f az + %f + %f z + %f z*az for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funct = a*rg**2 + b*rg+ c*az + d + f*topo_clean*az
                funcbins = a*rgbins**2 + b*rgbins+ c*azbins + d + f*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = np.arange((new_cols)) -  jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,5] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)


            elif (ivar==1 and nfit==1):

                G=np.zeros((len(data),7))
                G[:,0] = rg**2
                G[:,1] = rg
                G[:,2] = az
                G[:,3] = 1
                G[:,4] = topo_clean*az
                G[:,5] = topo_clean
                G[:,6] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]
                print ('Remove ramp %f r**2 + %f r  + %f az + %f + + %f z*az + %f z +%f z**2 for date: %i'%(a,b,c,d,e,f,g,idates[l]))

                # plot phase/elev
                funct = a*rg**2 + b*rg+ c*az + d + e*topo_clean*az
                funcbins = a*rgbins**2 + b*rgbins+ c*azbins + d + e*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,f*x+g*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi
                G[:,6] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = np.arange((new_cols)) -  jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,4] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)

            else:
                pass

      elif order==6:
        if arguments["--topofile"] is None:
            G=np.zeros((len(data),4))
            G[:,0] = az**2
            G[:,1] = az
            G[:,2] = rg
            G[:,3] = 1

            # ramp inversion
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]
            print ('Remove ramp %f az**2 + %f az  + %f r + %f for date: %i'%(a,b,c,d,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),4))
            for i in range(new_lines):
                G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                G[i*new_cols:(i+1)*new_cols,1] = (i - ibeg_emp)
                G[i*new_cols:(i+1)*new_cols,2] = np.arange((new_cols)) - jbeg_emp
            G[:,3] = 1


            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:
            if (ivar==0 and nfit==0) :
                G=np.zeros((len(data),5))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = rg
                G[:,3] = 1
                G[:,4] = topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
                print ('Remove ramp %f az**2 + %f az  + %f r + %f + %f z for date: %i'%(a,b,c,d,e,idates[l]))

                # plot phase/elev
                funct = a*az**2 + b*az + c*rg+ d
                funcbins = a*azbins**2 + b*azbins + c*rgbins+ d
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),5))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = i -  ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = np.arange((new_cols)) - jbeg_emp
                G[:,3] = 1
                G[:,4] = elevi


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            elif (ivar==0 and nfit==1):
                G=np.zeros((len(data),6))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = rg
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
                print ('Remove ramp %f az**2 + %f az  + %f r + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funct = a*az**2 + b*az + c*rg+ d
                funcbins = a*azbins**2 + b*azbins + c*rgbins+ d
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x+f*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = i -  ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = np.arange((new_cols)) - jbeg_emp
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),6))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = rg
                G[:,3] = 1
                G[:,4] = topo_clean
                G[:,5] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
                print ('Remove ramp %f az**2 + %f az  + %f r + %f + %f z + %f z*az for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funct = a*az**2 + b*az + c*rg+ d + f*topo_clean*az
                funcbins = a*azbins**2 + b*azbins + c*rgbins+ d + f*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = i -  ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,5] *= (i - ibeg_emp)


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),8))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = rg
                G[:,3] = 1
                G[:,4] = topo_clean*az
                G[:,5] = topo_clean
                G[:,6] = topo_clean**2
                G[:,7] = (topo_clean*az)**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g=pars[6]; h = pars[7]
                print ('Remove ramp %f az**2 + %f az  + %f r + %f + %f z*az + %f z + %f z**2 + %f (z*az)**2 for date: %i'%(a,b,c,d,e,f,g,h,idates[l]))

                # plot phase/elev
                funct = a*az**2 + b*az + c*rg+ d + e*topo_clean*az + h*(topo_clean*az)**2
                funcbins = a*azbins**2 + b*azbins + c*rgbins+ d + e*topobins*azbins + h*(topobins*azbins)**2
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,f*x+g*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),8))
                G[:,3] = 1
                G[:,4] = elevi
                G[:,5] = elevi
                G[:,6] = elevi**2
                G[:,7] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = i -  ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,4] *= (i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,7] *= (i - ibeg_emp)**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)


      elif order==7:
        if arguments["--topofile"] is None:
            G=np.zeros((len(data),5))
            G[:,0] = az**2
            G[:,1] = az
            G[:,2] = rg**2
            G[:,3] = rg
            G[:,4] = 1

            # ramp inversion
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
            print ('Remove ramp %f az**2 + %f az  + %f r**2 + %f r + %f for date: %i'%(a,b,c,d,e,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),5))
            for i in range(new_lines):
                G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                G[i*new_cols:(i+1)*new_cols,2] = (np.arange((new_cols)) - jbeg_emp)**2
                G[i*new_cols:(i+1)*new_cols,3] = np.arange((new_cols)) - jbeg_emp
            G[:,4] = 1


            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:
            if (ivar==0 and nfit ==0):
                G=np.zeros((len(data),6))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = rg**2
                G[:,3] = rg
                G[:,4] = 1
                G[:,5] = topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
                print ('Remove ramp %f az**2 + %f az  + %f r**2 + %f r + %f + %f z for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funct = a*az**2 + b*az + c*rg**2 + d*az+ e
                funcbins = a*azbins**2 + b*azbins + c*rgbins**2 + d*rgbins+ e
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,f*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,3] = np.arange((new_cols)) - jbeg_emp
                G[:,4] = 1
                G[:,5] = elevi


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            if (ivar==0 and nfit ==1):
                G=np.zeros((len(data),7))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = rg**2
                G[:,3] = rg
                G[:,4] = 1
                G[:,5] = topo_clean
                G[:,6] = topo_clean**2

                # ramp inversion
                x0 = lst.lstsq(G,data)[0]
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]
                print ('Remove ramp %f az**2 + %f az  + %f r**2 + %f r + %f + %f z + %f z**2  for date: %i'%(a,b,c,d,e,f,g,idates[l]))

                # plot phase/elev
                funct = a*az**2 + b*az + c*rg**2 + d*rg+ e
                funcbins = a*azbins**2 + b*azbins + c*rgbins**2 + d*rgbins+ e
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,f*x+g*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,3] = np.arange((new_cols)) - jbeg_emp
                G[:,4] = 1
                G[:,5] = elevi
                G[:,6] = elevi**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit ==0):
                G=np.zeros((len(data),7))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = rg**2
                G[:,3] = rg
                G[:,4] = 1
                G[:,5] = topo_clean
                G[:,6] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g=pars[6]
                print ('Remove ramp %f az**2 + %f az  + %f r**2 + %f r + %f + %f z + %f az*z for date: %i'%(a,b,c,d,e,f,g,idates[l]))

                # plot phase/elev
                funct = a*az**2 + b*az + c*rg**2 + d*rg+ e + g*topo_clean*az
                funcbins = a*azbins**2 + b*azbins + c*rgbins**2 + d*rgbins+ e + g*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,f*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                G[:,4] = 1
                G[:,5] = elevi
                G[:,6] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,3] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,6] *= (i - ibeg_emp)


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),8))
                G[:,0] = az**2
                G[:,1] = az
                G[:,2] = rg**2
                G[:,3] = rg
                G[:,4] = 1
                G[:,5] = topo_clean*az
                G[:,6] = topo_clean
                G[:,7] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g=pars[6]; h=pars[7]
                print ('Remove ramp %f az**2 + %f az  + %f r**2 + %f r + %f +  %f az*z + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,g,h,idates[l]))

                # plot phase/elev
                funct = a*az**2 + b*az + c*rg**2 + d*rg+ e + f*topo_clean*az
                funcbins = a*azbins**2 + b*azbins + c*rgbins**2 + d*rgbins+ e + f*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,g*x+h*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),8))
                G[:,4] = 1
                G[:,5] = elevi
                G[:,6] = elevi
                G[:,7] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,3] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,5] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)

      elif order==8:
        if arguments["--topofile"] is None:
            G=np.zeros((len(data),6))
            G[:,0] = az**3
            G[:,1] = az**2
            G[:,2] = az
            G[:,3] = rg**2
            G[:,4] = rg
            G[:,5] = 1

            # ramp inversion
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]
            print ('Remove ramp %f az**3 + %f az**2  + %f az + %f r**2 + %f r + %f for date: %i'%(a,b,c,d,e,f,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),6))
            for i in range(new_lines):
                G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**3
                G[i*new_cols:(i+1)*new_cols,1] = (i - ibeg_emp)**2
                G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                G[i*new_cols:(i+1)*new_cols,3] = (np.arange((new_cols)) - jbeg_emp)**2
                G[i*new_cols:(i+1)*new_cols,4] = (np.arange((new_cols)) - jbeg_emp)
            G[:,5] = 1


            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),7))
                G[:,0] = az**3
                G[:,1] = az**2
                G[:,2] = az
                G[:,3] = rg**2
                G[:,4] = rg
                G[:,5] = 1
                G[:,6] = topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]
                print ('Remove ramp %f az**3 + %f az**2  + %f az + %f r**2 + %f r + %f + %f z for date: %i'%(a,b,c,d,e,f,g,idates[l]))

                # plot phase/elev
                funct = a*az**3 + b*az**2 + c*az + d*rg**2 + e*rg+ f
                funcbins = a*azbins**3 + b*azbins**2 + c*azbins + d*rgbins**2 + e*rgbins+ f
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,g*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**3
                    G[i*new_cols:(i+1)*new_cols,1] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,3] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,4] = np.arange((new_cols)) - jbeg_emp
                G[:,5] = 1
                G[:,6] = elevi


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            if (ivar==0 and nfit==1):
                G=np.zeros((len(data),8))
                G[:,0] = az**3
                G[:,1] = az**2
                G[:,2] = az
                G[:,3] = rg**2
                G[:,4] = rg
                G[:,5] = 1
                G[:,6] = topo_clean
                G[:,7] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]; h = pars[7]
                print ('Remove ramp %f az**3 + %f az**2  + %f az + %f r**2 + %f r + %f + %f z + %f z**2 for date: %i'%(a,b,c,d,e,f,g,h,idates[l]))

                # plot phase/elev
                funct = a*az**3 + b*az**2 + c*az + d*rg**2 + e*rg+ f
                funcbins = a*azbins**3 + b*azbins**2 + c*azbins + d*rgbins**2 + e*rgbins+ f
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,g*x+h*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),8))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**3
                    G[i*new_cols:(i+1)*new_cols,1] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,3] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,4] = np.arange((new_cols)) - jbeg_emp
                G[:,5] = 1
                G[:,6] = elevi
                G[:,7] = elevi**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)


            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),8))
                G[:,0] = az**3
                G[:,1] = az**2
                G[:,2] = az
                G[:,3] = rg**2
                G[:,4] = rg
                G[:,5] = 1
                G[:,6] = topo_clean
                G[:,7] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]; h=pars[7]
                print ('Remove ramp %f az**3 + %f az**2  + %f az + %f r**2 + %f r + %f + %f z + %f z*az for date: %i'%(a,b,c,d,e,f,g,h,idates[l]))

                # plot phase/elev
                funct = a*az**3 + b*az**2 + c*az + d*rg**2 + e*rg + f + h*topo_clean*az
                funcbins = a*azbins**3 + b*azbins**2 + c*azbins + d*rgbins**2 + e*rgbins+ f + h*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,g*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),8))
                G[:,5] = 1
                G[:,6] = elevi
                G[:,7] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**3
                    G[i*new_cols:(i+1)*new_cols,1] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,3] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,4] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,7] *= (i - ibeg_emp)


                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),10))
                G[:,0] = az**3
                G[:,1] = az**2
                G[:,2] = az
                G[:,3] = rg**2
                G[:,4] = rg
                G[:,5] = 1
                G[:,6] = topo_clean*az
                G[:,7] = topo_clean
                G[:,8] = topo_clean**2
                G[:,9] = (topo_clean*az)**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]; h=pars[7]; i=pars[8]; k=pars[9]
                print ('Remove ramp %f az**3 + %f az**2  + %f az + %f r**2 + %f r + %f z*az + %f + %f z + %f z**2 + %f (z*az)**2 for date: %i'%(a,b,c,d,e,f,g,h,i,k,idates[l]))

                # plot phase/elev
                funct = a*az**3 + b*az**2 + c*az + d*rg**2 + e*rg + f + g*topo_clean*az + k*(topo_clean*az)**2
                funcbins = a*azbins**3 + b*azbins**2 + c*azbins + d*rgbins**2 + e*rgbins+ f + g*topobins*azbins + k*(topobins*azbins)**2
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,h*x+i*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),10))
                G[:,5] = 1
                G[:,6] = elevi
                G[:,7] = elevi
                G[:,8] = elevi**2
                G[:,9] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = (i - ibeg_emp)**3
                    G[i*new_cols:(i+1)*new_cols,1] = (i - ibeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,2] =(i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,3] = (np.arange((new_cols)) - jbeg_emp)**2
                    G[i*new_cols:(i+1)*new_cols,4] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,6] *= (i - ibeg_emp)
                    G[i*new_cols:(i+1)*new_cols,9] *= (i - ibeg_emp)**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)

      elif order==9:
        if arguments["--topofile"] is None:
            G=np.zeros((len(data),5))
            G[:,0] = rg
            G[:,1] = az
            G[:,2] = (rg*az)**2
            G[:,3] = rg*az
            G[:,4] = 1

            # ramp inversion
            pars = linear_inv(G, data, sigma)
            a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]
            print ('Remove ramp %f r + %f az  + %f r*az**2 + %f r*az + %f for date: %i'%(a,b,c,d,e,idates[l]))

            # build total G matrix
            G=np.zeros((len(los),5))
            for i in range(new_lines):
                G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                G[i*new_cols:(i+1)*new_cols,2] = ((i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp))**2
                G[i*new_cols:(i+1)*new_cols,3] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
            G[:,4] = 1

            res = los - np.dot(G,pars)
            ramp = np.dot(G,pars).reshape(new_lines,new_cols)

        else:
            if (ivar==0 and nfit==0):
                G=np.zeros((len(data),6))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = (rg*az)**2
                G[:,3] = rg*az
                G[:,4] = 1
                G[:,5] = topo_clean

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]

                print ('Remove ramp %f r + %f az  + %f (r*az)**2 + %f r*az + %f + %f z for date: %i'%(a,b,c,d,e,f,idates[l]))

                # plot phase/elev
                funcbins = a*rgbins+ b*azbins + c*(azbins*rgbins)**2 + d*azbins*rgbins+ e
                funct = a*rg+ b*az + c*(az*rg)**2 + d*az*rg+ e
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,e*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),6))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = ((i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp))**2
                    G[i*new_cols:(i+1)*new_cols,3] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
                G[:,4] = 1
                G[:,5] = elevi

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-1)],pars[:nparam-1]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-1):],pars[(nparam-1):]).reshape(new_lines,new_cols)

            if (ivar==0 and nfit==1):
                G=np.zeros((len(data),7))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = (rg*az)**2
                G[:,3] = rg*az
                G[:,4] = 1
                G[:,5] = topo_clean
                G[:,6] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5]; g = pars[6]

                print ('Remove ramp %f r + %f az  + %f (r*az)**2 + %f r*az + %f + %f z + %f z**2  for date: %i'%(a,b,c,d,e,f,g,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c*(rg*az)**2 + d*rg*az+ e
                funcbins = a*rgbins+ b*azbins + c*(azbins*rgbins)**2 + d*azbins*rgbins+ e
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,f*x+g*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = ((i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp))**2
                    G[i*new_cols:(i+1)*new_cols,3] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
                G[:,4] = 1
                G[:,5] = elevi
                G[:,6] = elevi**2

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==0):
                G=np.zeros((len(data),7))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = (rg*az)**2
                G[:,3] = rg*az
                G[:,4] = 1
                G[:,5] = topo_clean
                G[:,6] = topo_clean*az

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5] ; g = pars[6]

                print ('Remove ramp %f r + %f az  + %f (r*az)**2 + %f r*az + %f + %f z + %f az*z for date: %i'%(a,b,c,d,e,f,g,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c*(rg*az)**2 + d*rg*az+ e + g*topo_clean*az
                funcbins = a*rgbins+ b*azbins + c*(azbins*rgbins)**2 + d*azbins*rgbins+ e + g*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,f*x,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),7))
                G[:,4] = 1
                G[:,5] = elevi
                G[:,6] = elevi
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = ((i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp))**2
                    G[i*new_cols:(i+1)*new_cols,3] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
                    G[i*new_cols:(i+1)*new_cols,6] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-2)],pars[:nparam-2]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-2):],pars[(nparam-2):]).reshape(new_lines,new_cols)

            elif (ivar==1 and nfit==1):
                G=np.zeros((len(data),8))
                G[:,0] = rg
                G[:,1] = az
                G[:,2] = (rg*az)**2
                G[:,3] = rg*az
                G[:,4] = 1
                G[:,5] = topo_clean*az
                G[:,6] = topo_clean
                G[:,7] = topo_clean**2

                # ramp inversion
                pars = linear_inv(G, data, sigma)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5] ; g = pars[6]; h=pars[7]

                print ('Remove ramp %f r + %f az  + %f (r*az)**2 + %f r*az + %f + %f az*z + %f z + %f z**2  for date: %i'%(a,b,c,d,e,f,g,h,idates[l]))

                # plot phase/elev
                funct = a*rg+ b*az + c*(rg*az)**2 + d*rg*az+ e + f*topo_clean*az
                funcbins = a*rgbins+ b*azbins + c*(azbins*rgbins)**2 + d*azbins*rgbins+ e + f*topobins*azbins
                x = np.linspace(mintopo, maxtopo, 100)
                ax_dphi.scatter(topo_clean,los_clean-funct, s=0.01, alpha=0.3, rasterized=True)
                ax_dphi.plot(topobins,losbins - funcbins,'-r', lw =1., label='sliding median')
                ax_dphi.plot(x,g*x+h*x**2,'-r', lw =4.)

                # build total G matrix
                G=np.zeros((len(los),8))
                G[:,4] = 1
                G[:,5] = elevi
                G[:,6] = elevi
                G[:,7] = elevi**2
                for i in range(new_lines):
                    G[i*new_cols:(i+1)*new_cols,0] = np.arange((new_cols)) - jbeg_emp
                    G[i*new_cols:(i+1)*new_cols,1] = i - ibeg_emp
                    G[i*new_cols:(i+1)*new_cols,2] = ((i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp))**2
                    G[i*new_cols:(i+1)*new_cols,3] = (i-ibeg_emp) * (np.arange((new_cols))-jbeg_emp)
                    G[i*new_cols:(i+1)*new_cols,5] *= (i - ibeg_emp)

                res = los - np.dot(G,pars)
                nparam = G.shape[1]
                ramp = np.dot(G[:,:(nparam-3)],pars[:nparam-3]).reshape(new_lines,new_cols)
                topo = np.dot(G[:,(nparam-3):],pars[(nparam-3):]).reshape(new_lines,new_cols)

      rms = np.sqrt(np.nanmean(res**2))
      logger.info('RMS dates %i: %f'%(idates[l], rms))
      flata = los.reshape(new_lines,new_cols) - ramp - topo
      
      try:
         del G; del los
      except:
         pass
      return ramp, flata, topo, rms
 
