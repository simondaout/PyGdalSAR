#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################################
# Authors        : Hugo Watine, Leo Letellier (CRPG)
# inspired from NSBAS code flatten_stack 
############################################

"""
Model coefficient estimation and removal for wrapped interferograms

Usage:
    flatten_stack.py estimate <phase_filt> <model> [--outfile=<outfile>] [--nreg=<nreg>] [--thresh_amp=<float>] [--thresh_cohreg=<float>] [--thresh_model=<float>] [--mask_center=<yes/no>] [--thresh_std_model=<float>] [--thresh_min_pixel=<int>] [--cyclmax=<float>] [--plot_reg=<yes/no>] [--plot=<yes/no>] [--weightmedian=<yes/no>]
    flatten_stack.py add <unwrapped> <model> --coeff=<coeff> [--outfile=<outfile>]
    flatten_stack.py remove <phase> <model> --coeff=<coeff> [--outfile=<outfile>]
    flatten_stack.py 2model <phase> <phase_filt> <model1> <model2> [--outfile=<outfile>] [--nreg=<float>] [--thresh_amp=<float>] [--thresh_cohreg=<float>] [--thresh_model=<float>] [--mask_center=<yes/no>] [--thresh_std_model=<float>] [--thresh_model2=<float>] [--mask_center2=<yes/no>] [--thresh_std_model2=<float>] [--thresh_min_pixel=<int>] [--cyclmax=<float>] [--plot=<yes/no>] [--plot_reg=<yes/no>]
    flatten_stack.py <phase> <phase_filt> <model> [--outfile=<outfile>] [--nreg=<float>] [--thresh_amp=<float>] [--thresh_cohreg=<float>] [--thresh_model=<float>] [--mask_center=<yes/no>] [--thresh_std_model=<float>] [--thresh_min_pixel=<int>] [--cyclmax=<float>] [--plot=<yes/no>] [--plot_reg=<yes/no>] [--weightmedian=<yes/no>]
    flatten_stack.py -h | --help

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
    mask_center         Mask asb(model) < thresh_model insted of model < thresh_model (default: no)
    thresh_model2       Minimal value for model2 pixel selection (default: 0.4)
    mask_center2         Mask asb(model2) < thresh_model insted of model2 < thresh_model (default: no)
    thresh_std_model    Minimal standard deviation within a window for estimation (default: 0.3)
    thresh_std_model2   Minimal standard deviation for model2 in a window (default: 0.3)
    thresh_min_pixel    Minimal percent of pixel per windows for estimation (default: 10)
    cyclmax             Maximum number of phase cycle (default: 3)
    plot                plot intermediate results (default: no)
    plot_reg            plot each phase model relation per region (default: no)
    weightmedian        if yes do a weighted median based on the increase of coherence before and after removing the model (default: no)
"""

print()
print('Author: Hugo Watine (CRPG)')
print('inspired from NSBAS code flatten_stack')
print()

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

def plot_reg_model1(model1_region,phi_filt_region, list_coeff, list_coh, i, j ):
    plt.figure(figsize=(5, 3))
    plt.scatter(np.angle(np.exp(1j*model1_region)), np.angle(np.exp(1j*phi_filt_region)), s=0.1, c='k', alpha=0.05, label="Data", rasterized=True)
    
    x_fit = np.linspace(np.nanmin(np.angle(np.exp(1j*model1_region))), np.nanmax(np.angle(np.exp(1j*model1_region))), 100)
    y_fit = np.angle(np.exp(1j*coeff*x_fit))

    x = np.angle(np.exp(1j*model1_region))
    y = np.angle(np.exp(1j*phi_filt_region))

    y_fit2 = np.angle(np.exp(1j*coeff*model1_region))
    cst =  np.nanmean(y - y_fit2)

    plt.plot(x_fit, y_fit + cst, 'r-', lw=2, label=f"Fit: coeff={list_coeff[i,j]:.3f}, coh={list_coh[i,j]:.3f}")
    plt.xlabel("Model")
    plt.ylabel("Phase filtered")
    plt.legend()
    plt.title(f"Window ({i}, {j})")
    plt.grid(True, ls='--', alpha=0.4)
    plt.tight_layout()
    plt.show()

def plot_reg_model2(alpha, model1_region, beta, model2_region, phi_filt_region):
    try:
        plt.figure(figsize=(5, 3))
        x = np.angle(np.exp(1j * (alpha * model1_region + beta * model2_region)))
        y = np.angle(np.exp(1j * phi_filt_region))
        plt.scatter(x, y, s=0.5, c='k', alpha=0.05, rasterized=True)
        plt.xlabel("Combined model angle")
        plt.ylabel("Phase filtered")
        plt.title(f"Window ({i}, {j})")
        plt.grid(True, ls='--', alpha=0.4)
        plt.tight_layout()
        plt.show()
    except Exception:
        pass

def plot_final(list_coeff, list_coh, model1_copy):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 7))

    im1 = ax1.imshow(list_coeff, cmap='RdBu_r',
                vmax=cyclmax,
                vmin=-cyclmax)
    ax1.set_title('Coefficients')
    plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
    
    im2 = ax2.imshow(list_coh, cmap='viridis')
    ax2.set_title('Coherence')
    plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)

    im3 = ax3.imshow(model1_copy, cmap='RdBu_r',
                        vmin=np.nanpercentile(model1_copy, 2),
                        vmax=np.nanpercentile(model1_copy, 98))
    ax3.set_title('Model')
    plt.colorbar(im3, ax=ax3, fraction=0.046, pad=0.04)

    plt.tight_layout()
    plt.show()


def open_gdal(file, band=1, supp_ndv=None, complex=False):
    """
    Use GDAL to open band as real value or complex interferogram.
    Returns (data, Xsize, Ysize)
    If complex=True: data = [amplitude, phase]
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
    Create a new raster with same dims as the template or overwrite if exists.
    If data contains complex values (complex dtype), will write magnitude/phase splitted is not implemented.
    We assume template is a complex ROI_PAC-like stack when writing complex arrays.
    """

    if "." in outfile:
        final_outfile = outfile
    else:
        if template is None:
            raise RuntimeError("Template required when outfile has no extension")
        pos_rlks = template.rfind('rlks')
        pos_underscore = template.rfind('_', 0, pos_rlks)
        final_outfile = template[:pos_underscore+1] + outfile + '_' + template[pos_underscore+1:]


    # Create dataset
    ds_template = gdal.Open(template)
    driver = gdal.GetDriverByName("ROI_PAC" if roi_pac else "GTiff")
    ds = driver.CreateCopy(final_outfile, ds_template)

    band_obj = ds.GetRasterBand(band)
    ndv = alt_ndv if alt_ndv is not None else band_obj.GetNoDataValue()
    if ndv is not None:
        band_obj.SetNoDataValue(ndv)
        data = np.where(np.isnan(data), ndv, data)

    band_obj.WriteArray(data)
    ds.FlushCache()

    print("save:", final_outfile)


def former_flatten_stack(phase, phase_filt, model1, model2=None, **kwargs):
    """
    Remove a pattern from a wrapped signal, as in flatten_stack.f from NSBAS
    """
    coeffs = estimate_coeff(phase_filt, model1, model2=model2, **kwargs)

    if model2 is None:
        coeff = float(coeffs)
        nomodel_phifilt = remove_model_single(phase_filt[1], model1, coeff)
        nomodel_phi = remove_model_single(phase[1], model1, coeff)
    else:
        alpha, beta = coeffs
        nomodel_phifilt = remove_model_double(phase_filt[1], model1, model2, alpha, beta)
        nomodel_phi = remove_model_double(phase[1], model1, model2, alpha, beta)

    nomodel_phifilt[np.isnan(nomodel_phifilt)] = 0.0
    nomodel_phi[np.isnan(nomodel_phi)] = 0.0
    phase_filt[0][np.isnan(phase_filt[0])] = 0.0
    phase[0][np.isnan(phase[0])] = 0.0

    complex_nomodel_phifilt = phase_filt[0] * np.exp(1j*nomodel_phifilt)
    complex_nomodel_phi = phase[0] * np.exp(1j*nomodel_phi)

    if kwargs.get('plot', 'no') == 'yes':
        try:
            fig, axes = plt.subplots(1, 3, figsize=(12, 5))
            im0 = axes[0].imshow(phase_filt[1], cmap='twilight')
            axes[0].set_title('phase_filt (orig)')
            plt.colorbar(im0, ax=axes[0])
            im1 = axes[1].imshow(nomodel_phifilt, cmap='twilight')
            axes[1].set_title('phase_filt (corrected)')
            plt.colorbar(im1, ax=axes[1])
            if model2 is None:
                im2 = axes[2].imshow(model1, cmap='RdBu_r', vmin=np.nanpercentile(model1, 2), vmax=np.nanpercentile(model1, 98))
            else:
                im2 = axes[2].imshow(model1 + 0.0, cmap='RdBu_r', vmin=np.nanpercentile(model1, 2), vmax=np.nanpercentile(model1, 98))
            axes[2].set_title('Model')
            plt.colorbar(im2, ax=axes[2])
            plt.tight_layout()
            plt.show()
        except Exception:
            pass
    
    print('')
    save_gdal(kwargs['outfile'], complex_nomodel_phifilt, template=arguments["<phase_filt>"])
    save_gdal(kwargs['outfile'], complex_nomodel_phi, template=arguments["<phase>"])

    pos_rlks = arguments["<phase_filt>"].rfind('rlks')
    pos_underscore = arguments["<phase_filt>"].rfind('_', 0, pos_rlks)
    outfile = arguments["<phase_filt>"][:pos_underscore+1] + kwargs['outfile'] + '_' + arguments["<phase_filt>"][pos_underscore+1:]

    dates = next(part for part in outfile.split('_') if '-' in part)
    if kwargs['outfile'] != 'nomodel':
        dates = dates + '_' + kwargs['outfile']
    with open(f"{dates}.stack", 'w') as f:
        f.write(f"{coeff:.10f}\n")


def estimate_coeff(phase_filt, model, model2=None, **kwargs) -> float:
    """
    Estimate coefficient(s) by maximizing coherence (minimizing misfit).
    If model2 is None -> returns single float coeff
    If model2 provided -> returns tuple (alpha_med, beta_med)
    """
    amp_filt = np.copy(phase_filt[0])
    phi_filt = np.copy(phase_filt[1])
    model1_copy = np.copy(model)
    model2_copy = None
    if model2 is not None:
        model2_copy = np.copy(model2)
    
    if model2 is None:  
        if kwargs["mask_center"] == 'no':
            index = (np.isnan(amp_filt) | np.isnan(phi_filt)| np.isnan(model1_copy) | (amp_filt < kwargs['thresh_amp']) | (model1_copy < kwargs['thresh_model']))
        else:
            index = (np.isnan(amp_filt) | np.isnan(phi_filt)| np.isnan(model1_copy) | (amp_filt < kwargs['thresh_amp']) | (np.abs(model1_copy) < kwargs['thresh_model']))
    else:
        if kwargs["mask_center"] == 'no' and kwargs["mask_center2"] == 'no':
            index = (np.isnan(amp_filt) | np.isnan(phi_filt)| np.isnan(model1_copy) | np.isnan(model2_copy) | (amp_filt < kwargs['thresh_amp']) | ((model1_copy < kwargs['thresh_model']) & (model2_copy < kwargs['thresh_model2'])))
        else:
            cond1 = np.abs(model1_copy) < kwargs['thresh_model'] if kwargs["mask_center"] == 'yes' else model1_copy < kwargs['thresh_model']
            cond2 = np.abs(model2_copy) < kwargs['thresh_model2'] if kwargs["mask_center2"] == 'yes' else model2_copy < kwargs['thresh_model2']
            index = (np.isnan(amp_filt) | np.isnan(phi_filt)| np.isnan(model1_copy) | np.isnan(model2_copy) | (amp_filt < kwargs['thresh_amp']) | (cond1 & cond2))

    # if kwargs["mask_center"] == 'no':
    #     if model2 is not None:
    #         index = np.isnan(amp_filt) | np.isnan(phi_filt) | np.isnan(model1_copy) | np.isnan(model2_copy) | (amp_filt < kwargs['thresh_amp']) | (model1_copy < kwargs['thresh_model']) & (model2_copy < kwargs['thresh_model2'])
    #     else:
    #         index = np.isnan(amp_filt) | np.isnan(phi_filt) | np.isnan(model1_copy) | (amp_filt < kwargs['thresh_amp']) | (model1_copy < kwargs['thresh_model'])
    # else :
    #     index = np.isnan(amp_filt) | np.isnan(phi_filt) | np.isnan(model1_copy) | (amp_filt < kwargs['thresh_amp']) | (np.abs(model1_copy) < kwargs['thresh_model'])

    amp_filt[index] = np.nan
    phi_filt[index] = np.nan
    model1_copy[index] = np.nan
    if model2 is not None:
        model2_copy[index] = np.nan

    ## Creation of sub-area
    nreg = kwargs["nreg"]
    print('---')
    print('Number of region :', nreg)
    if nreg != 1:
        amp_filt = reg_pad(amp_filt, nreg)
        phi_filt = reg_pad(phi_filt, nreg)
        model1_split = reg_pad(model1_copy, nreg)
        if model2 is not None:
            model2_blocks = reg_pad(model2_copy, nreg)
    else:
        amp_filt = np.array([[amp_filt]])
        phi_filt = np.array([[phi_filt]])
        model1_split = np.array([[model1_copy]])
        if model2 is not None:
            model2_split = reg_pad(model2_copy, nreg)
        kwargs['thresh_cohreg'] = 0.0           

    # Initialisation
    cyclmax = kwargs['cyclmax']
    ny = model1_split.shape[0]
    nx = model1_split.shape[1]

    if model2 is None:
        list_coeff = np.full((ny, nx), np.nan)
        list_coh = np.full((ny, nx), np.nan)
    else:
        list_alpha = np.full((ny, nx), np.nan)
        list_beta = np.full((ny, nx), np.nan)
        list_coh = np.full((ny, nx), np.nan)

    if kwargs['weightmedian'] == 'yes':
        list_coh_initial = np.full((ny, nx), np.nan)

    print('\nStarting maximisation')
    for i in range(ny):     
        for j in range(nx):
            # for each region
            phi_filt_region = phi_filt[i,j]
            model1_region = model1_split[i,j]

            if kwargs['weightmedian'] == 'yes':
                list_coh_initial[i,j] = np.abs(np.nanmean(np.exp(1j * phi_filt_region)))

            #kwargs['thresh_std_model1'] = 0.0
            if np.count_nonzero(~np.isnan(phi_filt_region)) < kwargs['thresh_min_pixel']:
                # print("skip not enough px", i, j)
                continue
            elif np.nanstd(model1_region) < kwargs['thresh_std_model']:
                # print("skip not enough std", i, j)
                continue
            
            #### THINK about what to do if the std of model 2 is low ###
            

            if model2 is None:
                res = minimize(misfit_single, x0=[0.0], args=(phi_filt_region, model1_region), 
                            bounds=[(-cyclmax, cyclmax) ] )
                coeff = res.x[0]
                coh = np.abs(np.nanmean(np.exp(1j * phi_filt_region) * np.exp(-1j * (coeff * model1_region ))))
                
                if coh > kwargs['thresh_cohreg']:
                    list_coeff[i,j] = coeff
                    list_coh[i,j] = coh 

                if kwargs['plot_reg'] == 'yes':
                    plot_reg_model1(model1_region,phi_filt_region, list_coeff, list_coh, i, j)
                   
            else :
                model2_region = model2_blocks[i, j]
                #kwargs['thresh_std_model2'] = 0.0
                if np.nanstd(model2_region) < kwargs['thresh_std_model2']:
                    continue

                res = minimize(misfit_double, x0=[0.0, 0.0], args=(phi_filt_region, model1_region, model2_region),
                               bounds=[(-cyclmax, cyclmax), (-cyclmax, cyclmax)])
                alpha = float(res.x[0]); beta = float(res.x[1])
                coh = np.abs(np.nanmean(np.exp(1j * phi_filt_region) * np.exp(-1j * (alpha * model1_region + beta * model2_region))))
                if coh > kwargs['thresh_cohreg']:
                    list_alpha[i, j] = alpha
                    list_beta[i, j] = beta
                    list_coh[i, j] = coh

                if kwargs['plot_reg'] == 'yes':
                    plot_reg_model2(alpha, model1_region, beta, model2_region, phi_filt_region)
        
    if model2 is None:
        if kwargs['weightmedian'] == 'no':
            med = np.nanmedian(list_coeff)
        else:
            mask = ~np.isnan(list_coh_initial) & ~np.isnan(list_coh)
            weights = np.zeros_like(list_coeff)
            valid_mask = mask & (list_coh_initial != 0)
            weights[valid_mask] = (list_coh[valid_mask] - list_coh_initial[valid_mask]) / list_coh_initial[valid_mask]
            zero_mask = mask & (list_coh_initial == 0)
            weights[zero_mask] = list_coh[zero_mask] - list_coh_initial[zero_mask]
            weights = np.abs(weights)
            if np.sum(weights) > 0:
                weights = weights / np.sum(weights)
            else:
                weights = np.ones_like(weights) / len(weights)

            med = weighted_median(list_coeff, weights)

        if kwargs['plot'] == 'yes':
            plot_final(list_coeff, list_coh, model1_copy)

        print(f'--> coeff median: {med:.3f}')
        return med
    
    else:
        valid_alpha = list_alpha[~np.isnan(list_alpha)]
        valid_beta  = list_beta[~np.isnan(list_beta)]

        print("alpha (valid):", valid_alpha)
        print("beta  (valid):",  valid_beta)

        alpha_med = np.nanmedian(valid_alpha)
        beta_med  = np.nanmedian(valid_beta)
        print(f'--> coeff median: alpha={alpha_med:.3f}, beta={beta_med:.3f}')

        # === PLOTS ==========================================================
        if kwargs['plot'] == 'yes':

            fig, axes = plt.subplots(1, 5, figsize=(15, 6))

            # α map
            im0 = axes[0].imshow(list_alpha, cmap='RdBu_r',
                                 vmin=-cyclmax, vmax=cyclmax)
            axes[0].set_title("Alpha coefficient")
            plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

            # β map
            im1 = axes[1].imshow(list_beta, cmap='RdBu_r',
                                 vmin=-cyclmax, vmax=cyclmax)
            axes[1].set_title("Beta coefficient")
            plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

            # Coherence
            im2 = axes[2].imshow(list_coh, cmap='viridis')
            axes[2].set_title("Coherence")
            plt.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

            im3 = axes[3].imshow(model1_copy, cmap='RdBu_r',
                             vmin=np.nanpercentile(model1_copy, 2),
                             vmax=np.nanpercentile(model1_copy, 98))
            axes[3].set_title('Model 1')
            plt.colorbar(im3, ax=axes[3], fraction=0.046, pad=0.04)

            im4 = axes[4].imshow(model2_copy, cmap='RdBu_r',
                             vmin=np.nanpercentile(model2_copy, 2),
                             vmax=np.nanpercentile(model2_copy, 98))
            axes[4].set_title('Model 2')
            plt.colorbar(im4, ax=axes[4], fraction=0.046, pad=0.04)

            plt.tight_layout()
            plt.show()
        sys.exit()

        return alpha_med, beta_med

#    return med

def misfit_single(params, phi_filt_region, model_region):
    """
    params: [c]
    """
    c = params  
    coh = np.nanmean(np.exp(1j * phi_filt_region) *
            np.exp(-1j * (c * model_region )))
    return 1./np.abs(coh)

def misfit_double(params, phi_filt_region, model1_region, model2_region):
    """
    params: [alpha, beta]
    """
    alpha = params[0]
    beta = params[1]
    total = alpha * model1_region + beta * model2_region
    coh = np.nanmean(np.exp(1j * phi_filt_region) * np.exp(-1j * total))
    return 1.0 / np.abs(coh)
    
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
def remove_model_single(phase, model, coeff):
    """
    Remove a single model from wrapped phase values using coefficient coeff.
    phase : wrapped phase (radians)
    model : model raster
    """
    corrected_phase = np.angle(np.exp(1j*phase) * np.exp(-1j*(coeff*model)))
    return corrected_phase

def remove_model_double(phase, model1, model2, alpha, beta):
    """
    Remove combined model alpha*model1 + beta*model2 from wrapped phase values.
    """
    total = alpha * model1 + beta * model2
    corrected_phase = np.angle(np.exp(1j * phase) * np.exp(-1j * total))
    return corrected_phase


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
    thresh_model2 = arg2value(arguments.get("--thresh_model2"), float, 0.4)
    thresh_std_model = arg2value(arguments["--thresh_std_model"], float, 0.3)
    thresh_std_model2 = arg2value(arguments.get("--thresh_std_model2"), float, 0.3)
    thresh_min_pixel = arg2value(arguments["--thresh_min_pixel"], float, 10) # nb px (TODO en %)
    cyclmax = arg2value(arguments["--cyclmax"], float, 3.)
    plot = arg2value(arguments["--plot"], str, 'no')
    plot_reg = arg2value(arguments["--plot_reg"], str, 'no')
    mask_center = arg2value(arguments["--mask_center"], str, 'no')
    mask_center2 = arg2value(arguments["--mask_center2"], str, 'no')
    weightmedian = arg2value(arguments["--weightmedian"], str, 'no')

    if arguments["estimate"]:
        # Find the optimal coefficient of proportionality between the phase and the model
        #phase, pXsize, pYsize = open_gdal(arguments["<phase>"], band=1, complex=True)
        
        phase_filt, pfXsize, pfYsize = open_gdal(arguments["<phase_filt>"], band=1, complex=True)
        model, mXsize, mYsize = open_gdal(arguments["<model>"])
        if not (pfXsize == mXsize and pfYsize == mYsize):
            sys.exit("Error: input rasters (phase, phase_filt, model) do not have the same dimensions.")
        coeff = estimate_coeff(phase_filt, model, nreg=nreg, thresh_amp=thresh_amp, thresh_cohreg=thresh_cohreg,thresh_model=thresh_model, 
                       thresh_std_model=thresh_std_model, thresh_min_pixel=thresh_min_pixel, cyclmax=cyclmax, plot=plot, plot_reg=plot_reg, mask_center=mask_center)
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

        coeff = float(arguments["--coeff"])
        outfile = arguments["--outfile"]
        phase_minus_model = remove_model_single(phase[1].astype(np.float64), model.astype(np.float64), float(coeff))

        complex_nomodel_phi = phase[0] * np.exp(1j*phase_minus_model)
        save_gdal(outfile, complex_nomodel_phi, template=arguments["<phase>"])

    elif arguments.get('2model'):
        # 2model mode: <phase> <phase_filt> <model1> <model2>

        phase, pX, pY = open_gdal(arguments["<phase>"], complex=True)
        phase_filt, pfX, pfY = open_gdal(arguments["<phase_filt>"], complex=True)
        model1, m1X, m1Y = open_gdal(arguments["<model1>"])
        model2, m2X, m2Y = open_gdal(arguments["<model2>"])

        if not (pX == pfX == m1X == m2X and pY == pfY == m1Y == m2Y):
            sys.exit("Error: input rasters do not have the same dimensions.")

        # run flatten for two models (silent: no printing of coeff)
        former_flatten_stack(phase, phase_filt, model1, model2=model2, outfile=outfile, nreg=nreg, thresh_amp=thresh_amp, thresh_cohreg=thresh_cohreg, thresh_model=thresh_model, 
                                thresh_std_model=thresh_std_model, thresh_min_pixel=thresh_min_pixel, cyclmax=cyclmax, plot=plot, plot_reg=plot_reg, thresh_model2=thresh_model2, thresh_std_model2=thresh_std_model2, mask_center=mask_center, mask_center2=mask_center2, weightmedian=weightmedian)

    else:
        # Remove a pattern from a wrapped signal, as in flatten_stack.f from NSBAS
        phase, pXsize, pYsize = open_gdal(arguments["<phase>"], complex=True)
        phase_filt, pfXsize, pfYsize = open_gdal(arguments["<phase_filt>"], complex=True)
        model, mXsize, mYsize = open_gdal(arguments["<model>"])

        if not (pXsize == pfXsize == mXsize and pYsize == pfYsize == mYsize):
            sys.exit("Error: input rasters (phase, phase_filt, model) do not have the same dimensions.")
        phase_minus_model = former_flatten_stack(phase, phase_filt, model,
                                outfile=outfile, nreg=nreg, thresh_amp=thresh_amp, thresh_cohreg=thresh_cohreg, thresh_model=thresh_model, 
                                thresh_std_model=thresh_std_model, thresh_min_pixel=thresh_min_pixel, cyclmax=cyclmax, plot=plot, plot_reg=plot_reg, thresh_model2=None,  mask_center=mask_center, mask_center2=mask_center2, weightmedian=weightmedian)

