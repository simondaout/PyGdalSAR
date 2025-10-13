#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Model coefficient estimation and removal for wrapped interferograms

Usage:
    flatten_stack.py <phase> <phase_filt> <model> [--outfile=<outfile>] [--nreg=<float>] [--thresh_amp=<float>] [--thresh_cohreg=<float>] [--thresh_model=<float>] [--thresh_std_model=<float>] [--thresh_min_pixel=<int>] [--cyclmax=<float>] [--plot=<yes/no>] [--plot_reg=<yes/no>]
    flatten_stack.py -h | --help
    flatten_stack.py estimate <phase_filt> <model> [--outfile=<outfile>] [--nreg=<nreg>] [--thresh_amp=<float>] [--thresh_cohreg=<float>] [--thresh_model=<float>] [--thresh_std_model=<float>] [--thresh_min_pixel=<int>] [--cyclmax=<float>] [--plot_reg=<yes/no>]
    flatten_stack.py add <unwrapped> <model> --coeff=<coeff> [--outfile=<outfile>]
    flatten_stack.py remove <phase> <model> --coeff=<coeff> [--outfile=<outfile>]

Options:
    -h --help           Display command details
    phase               None filtered interferogram (.int)       
    phase_filt          Filtered interferogram (.int)
    model               Model for the propotionality estimation (.r4)
    outfile             Outfile name
    nreg                Number of region in range, put 1 for 1 global region (default: 28)
    thresh_amp          Minimal value of amplitude for pixel selection (default: 0.1)
    thresh_cohreg       Minimal value of coherence for region selection (default: 0.2)
    thresh_model        Minimal value of model for pixel selection (default: 0.4)
    thresh_std_model    Minimal standard deviation within a window for estimation (default: 0.3)
    thresh_min_pixel    Minimal percent of pixel per windows for estimation (default: 10)
    cyclmax             Maximum number of phase cycle (default: 3)
    plot                plot intermediate results (default: no)
    plot_reg            plot each phase model relation per region (default: no)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from osgeo import gdal
from numba import jit
import docopt
import sys
import os
from scipy.optimize import minimize
gdal.UseExceptions()

def open_gdal(file, band=1, supp_ndv=None, complex=False):
    """
    Use GDAL to open band as real value
    """
    print('-----')
    print(file)
    if not os.path.isfile(file):
        raise FileNotFoundError('File does not exists: {}'.format(file))
    ds = gdal.Open(file)

    print('dims', ds.RasterXSize, ds.RasterYSize, ds.RasterCount)
    band = ds.GetRasterBand(band)
    ndv = band.GetNoDataValue()
    data = band.ReadAsArray()
    Xsize = ds.RasterXSize
    Ysize = ds.RasterYSize
    if complex:
        amp = np.absolute(data)
        phi = np.angle(data)
        data = [amp, phi]
        ndv = 0.0
        if ndv is not None and ndv != np.nan:
            amp[amp==ndv] = np.nan
            phi[phi==ndv] = np.nan
        if supp_ndv is not None and supp_ndv != np.nan:
            amp[amp==supp_ndv] = np.nan
            phi[phi==supp_ndv] = np.nan
    else:
        if ndv is not None and ndv != np.nan:
            data[data==ndv] = np.nan
        if supp_ndv is not None and supp_ndv != np.nan:
            data[data==supp_ndv] = np.nan
    return data, Xsize, Ysize

def save_gdal(outfile, data, template=None, band=1, alt_ndv=None, roi_pac=True):
    """
    Create a new raster with same dims as the template
    """
    if outfile == 'nomodel':
        pos_rlks = template.rfind('rlks')
        pos_underscore = template.rfind('_', 0, pos_rlks)
        outfile = template[:pos_underscore+1] + outfile + '_' + template[pos_underscore+1:]
    if template is not None:
        ds_template = gdal.Open(template)
        if ds_template is None:
            raise RuntimeError(f"Imposible to open the template {template}")
        driver_name = "ROI_PAC" if roi_pac else "GTiff"
        driver = gdal.GetDriverByName(driver_name)
        ds = driver.CreateCopy(outfile, ds_template)
    else:
        ds = gdal.OpenEx(outfile, gdal.GA_Update)
    band = ds.GetRasterBand(band)
    ndv = band.GetNoDataValue()
    if alt_ndv is not None:
        ndv = alt_ndv
        band.SetNoDataValue(ndv)
    if ndv is not None:
        data[np.isnan(data)] = ndv
    band.WriteArray(data)
    ds.FlushCache()
    print('save : ', outfile)


def former_flatten_stack(phase, phase_filt, model, **kwargs):
    """
    Remove a pattern from a wrapped signal, as in flatten_stack.f from NSBAS
    """

    # nregs = [1, 5, 10, 20, 40, 80, 160, 200]
    # thresh_model = [-0.6, -0.3, 0, 0.3, 0.4, 0.6]
    # #nregs = [1, 200]
    # #thresh_model = [0, 0.6]
    # results_mediane = {t: [] for t in thresh_model}
    # for nreg in nregs:
    #     kwargs["nreg"] = nreg
    #     for thresh in thresh_model:
    #         kwargs["thresh_model"] = thresh
    #         median = estimate_coeff(phase_filt, model, **kwargs)
    #         results_mediane[thresh].append(median)
    # colors = plt.cm.viridis(np.linspace(0, 1, len(thresh_model)))
    # for color, thresh in zip(colors, thresh_model):
    #     plt.plot(nregs, results_mediane[thresh], '-o', color=color,
    #             label=f'médiane (thresh={thresh})', linewidth=2)
        
    # plt.xlabel('nreg', fontsize=12)
    # plt.ylabel('Coefficient', fontsize=12)
    # plt.title('Évolution des coefficients selon nreg et thresh_model', fontsize=13)
    # plt.grid(alpha=0.3)
    # plt.legend(frameon=False)
    # plt.tight_layout()
    # plt.show()

    coeff = estimate_coeff(phase_filt, model, **kwargs)
    nomodel_phifilt = remove_model(phase_filt[1], model, coeff)
    nomodel_phi = remove_model(phase[1], model, coeff)

    nomodel_phifilt[np.isnan(nomodel_phifilt)] = 0.0
    nomodel_phi[np.isnan(nomodel_phi)] = 0.0
    phase_filt[0][np.isnan(phase_filt[0])] = 0.0
    phase[0][np.isnan(phase[0])] = 0.0

    complex_nomodel_phifilt = phase_filt[0] * np.exp(1j*nomodel_phifilt)
    complex_nomodel_phi = phase[0] * np.exp(1j*nomodel_phi)

    if kwargs['plot'] == 'yes':
        fig, axes = plt.subplots(1, 3, figsize=(12, 5))

        im0 = axes[0].imshow(phase_filt[1], cmap='twilight')
        axes[0].set_title('phase')
        plt.colorbar(im0, ax=axes[0])

        im1 = axes[1].imshow(nomodel_phifilt, cmap='twilight')
        axes[1].set_title('corrected_phase')
        plt.colorbar(im1, ax=axes[1])
        
        im1 = axes[2].imshow(model, cmap='RdBu_r',vmin=np.nanpercentile(model, 2), vmax=np.nanpercentile(model, 98))
        axes[2].set_title('Model')
        plt.colorbar(im1, ax=axes[2])
        plt.show()
    
    print('')
    save_gdal(kwargs['outfile'], complex_nomodel_phifilt, template=arguments["<phase_filt>"])
    save_gdal(kwargs['outfile'], complex_nomodel_phi, template=arguments["<phase>"])

    dates = next(part for part in kwargs['outfile'].split('_') if '-' in part)
    with open(f"{dates}.stack", 'w') as f:
        f.write(f"{coeff:.10f}\n")

def estimate_coeff(phase_filt, model, **kwargs) -> float:
    """
    Estimate the proportionnality coefficient between the model and the wrapped phase values
    minimizing the complex product (maximizing the coherence)
    """
    amp_filt = np.copy(phase_filt[0])
    phi_filt = np.copy(phase_filt[1])
    model_copy = np.copy(model)
    
    index = np.isnan(amp_filt)| np.isnan(phi_filt) | np.isnan(model) | (amp_filt<kwargs['thresh_amp']) | (model<kwargs['thresh_model'])

    amp_filt[index] = np.nan
    phi_filt[index] = np.nan
    model_copy[index] = np.nan

    ## Creation of sub-area
    nreg = kwargs["nreg"]
    print('---')
    print('Number of region :', nreg)
    if nreg != 1:
        amp_filt = reg_pad(amp_filt, nreg)
        phi_filt = reg_pad(phi_filt, nreg)
        model_split = reg_pad(model_copy, nreg)
    else:
        amp_filt = np.array([[amp_filt]])
        phi_filt = np.array([[phi_filt]])
        model_split = np.array([[model_copy]])
        kwargs['thresh_cohreg'] = 0.           

    # Initialisation
    cyclmax = kwargs['cyclmax']
    list_coeff = np.full(shape=(model_split.shape[0], model_split.shape[1]), fill_value=np.nan)
    list_coh = np.full(shape=(model_split.shape[0], model_split.shape[1]), fill_value=np.nan)

    print('\nStarting maximisation')
    for i in range(model_split.shape[0]):     
        for j in range(model_split.shape[1]):
            # for each region
            phi_filt_region = phi_filt[i,j]
            model_region = model_split[i,j]

            if np.count_nonzero(~np.isnan(phi_filt_region)) < kwargs['thresh_min_pixel']:
                # print("skip not enough px", i, j)
                continue
            elif np.nanstd(model_region) < kwargs['thresh_std_model']:
                # print("skip not enough std", i, j)
                continue

            res = minimize(misfit, x0=[0.0], args=(phi_filt_region, model_region), 
                          bounds=[(-cyclmax, cyclmax) ] )
            coeff = res.x[0]

            coh = np.abs(np.nanmean(np.exp(1j * phi_filt_region) * np.exp(-1j * (coeff * model_region ))))
            
            if coh > kwargs['thresh_cohreg']:
                list_coeff[i,j] = coeff
                list_coh[i,j] = coh 

            if kwargs['plot_reg'] == 'yes':
                plt.figure(figsize=(5, 3))
                plt.scatter(np.angle(np.exp(1j*model_region)), np.angle(np.exp(1j*phi_filt_region)), s=0.1, c='k', alpha=0.05, label="Data", rasterized=True)
                
                x_fit = np.linspace(np.nanmin(np.angle(np.exp(1j*model_region))), np.nanmax(np.angle(np.exp(1j*model_region))), 100)
                y_fit = np.angle(np.exp(1j*coeff*x_fit))

                x = np.angle(np.exp(1j*model_region))
                y = np.angle(np.exp(1j*phi_filt_region))

                y_fit2 = np.angle(np.exp(1j*coeff*model_region))
                cst =  np.nanmean(y - y_fit2)

                plt.plot(x_fit, y_fit + cst, 'r-', lw=2, label=f"Fit: coeff={list_coeff[i,j]:.3f}, coh={list_coh[i,j]:.3f}")
                plt.xlabel("Model")
                plt.ylabel("Phase filtered")
                plt.legend()
                plt.title(f"Window ({i}, {j})")
                plt.grid(True, ls='--', alpha=0.4)
                plt.tight_layout()
                plt.show()
    
    med = np.nanmedian(list_coeff)
    print(f'--> coeff median: {med:.3f}')

    if kwargs['plot'] == 'yes':
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 7))
        im1 = ax1.imshow(list_coeff, cmap='RdBu_r',
                    vmax=cyclmax,
                    vmin=-cyclmax)
        ax1.set_title('Coefficients')
        plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
        
        im2 = ax2.imshow(list_coh, cmap='viridis')
        ax2.set_title('Coherence')
        plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)

        im3 = ax3.imshow(model_copy, cmap='RdBu_r', vmin=np.nanpercentile(model, 2), vmax=np.nanpercentile(model, 98))
        ax3.set_title('Model')
        plt.colorbar(im3, ax=ax3, fraction=0.046, pad=0.04)
        plt.tight_layout()
        plt.show()

    return med

def misfit(params, phi_filt_region, model_region):
    c = params  
    coh = np.nanmean(np.exp(1j * phi_filt_region) *
            np.exp(-1j * (c * model_region )))
    return 1./np.abs(coh)
    
def reg_pad(array, block_nb_x):
    """
    Regionalize an array from a given block number in the x dimension.
    The squared block size is determined from the x dimension, by padding the array
    to ensure nx % block_size = 0.
    The number of blocks in the y dimension is determined by conserving the block size 
    constant between x and y axis.
    """
    ny, nx = array.shape
    # taille de bloc à nombre de blocs fixé
    dxy = nx // block_nb_x + 1
    pad_x = dxy * block_nb_x - nx

    # Nombre de blocs à taille de bloc fixée
    nby = ny // dxy + 1
    pad_y = dxy * nby - ny
    
    # Pad to ensure that it can perfectly be divided into blocks
    pad_array = np.pad(array, ((0, pad_y), (0, pad_x)), constant_values=np.nan)
    # Cut the different blocks into different sub-arrays
    y_blocks = np.array(np.array_split(pad_array, nby, axis=0))
    blocks = np.array([np.array_split(b, block_nb_x, axis=1) for b in y_blocks])
    
    return blocks

#@jit(nopython=True)
def remove_model(phase, model, coeff):
    """
    Remove a model from wrapped values with a given proportionnality coefficient
    """
    corrected_phase = np.angle(np.exp(1j*phase) * np.exp(-1j*(coeff*model)))
    return corrected_phase

@jit
def weighted_median(values, weights):
    mask = ~np.isnan(values)
    values = values[mask] ; weights = weights[mask]

    sorter = np.argsort(values)
    values = values[sorter]
    weights = weights[sorter]

    cumulative_weight = np.cumsum(weights)

    return values[cumulative_weight >= cumulative_weight[-1] / 2.0][0]

def add_model(phase, model, coeff):
    """
    Add a model to values given a proportionnality coefficient
    """
    return (phase + coeff*model)

def arg2value(value, conversion=None, default=None):
    """Convert a string argument if exists otherwise use default value"""
    if value is None:
        return default
    elif conversion is not None:
        return conversion(value)
    return value


if __name__ == '__main__':
    arguments = docopt.docopt(__doc__)
    print(arguments)
    
    # Additionnal parameters
    outfile = arg2value(arguments["--outfile"], str, 'nomodel')
    nreg = arg2value(arguments["--nreg"], int, 28)
    thresh_amp = arg2value(arguments["--thresh_amp"], float, 0.1)
    thresh_cohreg = arg2value(arguments["--thresh_cohreg"], float, 0.2)
    thresh_model = arg2value(arguments["--thresh_model"], float, 0.4)
    thresh_std_model = arg2value(arguments["--thresh_std_model"], float, 0.3)
    thresh_min_pixel = arg2value(arguments["--thresh_min_pixel"], float, 10) # nb px (TODO en %)
    cyclmax = arg2value(arguments["--cyclmax"], float, 3.)
    plot = arg2value(arguments["--plot"], str, 'no')
    plot_reg = arg2value(arguments["--plot_reg"], str, 'no')

    if arguments["estimate"]:
        # Find the optimal coefficient of proportionality between the phase and the model
        phase, pXsize, pYsize = open_gdal(arguments["<phase>"], band=1, complex=True)
        phase_filt, pfXsize, pfYsize = open_gdal(arguments["<phase_filt>"], band=1, complex=True)
        model, mXsize, mYsize = open_gdal(arguments["<model>"])
        if not (pXsize == pfXsize == mXsize and pYsize == pfYsize == mYsize):
            sys.exit("Error: input rasters (phase, phase_filt, model) do not have the same dimensions.")
        coeff = estimate_coeff(phase_filt, model, nreg=nreg, thresh_amp=thresh_amp, thresh_cohreg=thresh_cohreg,thresh_model=thresh_model, 
                       thresh_std_model=thresh_std_model, thresh_min_pixel=thresh_min_pixel, cyclmax=cyclmax, plot=plot, plot_reg=plot_reg)
        print(coeff)
    elif arguments["add"]:
        # Add the model to an unwrapped file
        unwrapped, unwXsize, unwYsize = open_gdal(arguments["<unwrapped>"], band=2)
        model, mXsize, mYsize = open_gdal(arguments["<model>"])
        if not (unwXsize == mXsize and unwYsize == mYsize):
            sys.exit("Error: input rasters (unw, model) do not have the same dimensions.")
        coeff = float(arguments["--coeff"])
        outfile = arguments["--outfile"]
        unw_plus_model = add_model(unwrapped, model, coeff)
        save_gdal(outfile, unw_plus_model, template=arguments["<unwrapped>"], band=2)
    elif arguments["remove"]:
        # Remove the model from wrapped signal
        phase, pXsize, pYsize = open_gdal(arguments["<phase>"], band=1, complex=True)
        model, mXsize, mYsize = open_gdal(arguments["<model>"])
        if not (pXsize == mXsize and pYsize == mYsize):
            sys.exit("Error: input rasters (phase, model) do not have the same dimensions.")
        
        coeff = float(arguments["<coeff>"])
        outfile = arguments["--outfile"]
        phase_minus_model = remove_model(phase.astype(np.float64), modelastype(np.float64), float(coeff))
        save_gdal(outfile, phase, phase_minus_model, band=2)
    else:
        # Remove a pattern from a wrapped signal, as in flatten_stack.f from NSBAS
        phase, pXsize, pYsize = open_gdal(arguments["<phase>"], complex=True)
        phase_filt, pfXsize, pfYsize = open_gdal(arguments["<phase_filt>"], complex=True)
        model, mXsize, mYsize = open_gdal(arguments["<model>"])

        if not (pXsize == pfXsize == mXsize and pYsize == pfYsize == mYsize):
            sys.exit("Error: input rasters (phase, phase_filt, model) do not have the same dimensions.")
        phase_minus_model = former_flatten_stack(phase, phase_filt, model, 
                                outfile=outfile, nreg=nreg, thresh_amp=thresh_amp, thresh_cohreg=thresh_cohreg, thresh_model=thresh_model, 
                                thresh_std_model=thresh_std_model, thresh_min_pixel=thresh_min_pixel, cyclmax=cyclmax, plot=plot, plot_reg=plot_reg)

