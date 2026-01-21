#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Derived from plot_raster.py
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Mathieu Volat 
#                 Simon Daout (CRPG-ENSG)
# Revision      : Leo Letellier
############################################

"""
rplot.py
-------------
Display and Cut image file (.unw/.int/.r4/.tiff)

Usage: rplot.py <infile> [--cpt=<values>] [--crop=<values>] \
[--dim=<dim> | --gdal | --lectfile=<lectfile> | --lectcube=<lectcube> | --parfile=<parfile> | --amfile=<amfile>] \
[--rad2mm=<rad2mm>] [--title=<title>] [--wrap=<wrap>] [--vmin=<vmin>] [--vmax=<vmax>] [--band=<band>] \
[--cols=<cols>] [--lines=<lines>] [--zoom=<zoom>] [--histo] [--save] [--ndv=<ndv>] [--stats] \
[--bg=<bg>] [--alpha=<alpha>] [--vario] [--samples=<samples>] [--dlag=<dlag>] [--nlag=<nlag>] [--model=<model>] \
[--res=<res>]


Options:
-h --help               Show this screen.
<infile>                Raster to be displayed
--dim=<dim>             Indicate the raster's dimension (x,y or x,y,band)
--gdal                  Force the openning with GDAL
--lectfile=<lectfile>   Path of the lect.in file for r4 format
--lectcube=<lectcube>   Path to lect.in file containing band metadata
--parfile=<parfile>     Path of the .par file of GAMMA
--amfile=<amfile>       Path of the AMSTer InsarParameter file
--crop=<crop>           Crop option ("xmin,xmax,ymin,ymax")
--cpt=<cpt>             Indicate colorscale for phase
--wrap=<wrap>           Wrapped phase between value for unwrapped files 
--rad2mm=<rad2mm>       Convert data [default: 1]
--title=<title>         Title plot 
--band=<band>           Select band number [default: 1] 
--vmax=<vmax>           Max colorscale [default: 98th percentile]
--vmin=<vmin>           Min colorscale [default: 2th percentile]
--cols=<cols>           Add marker on pixel column numbers (eg. 200,400,450)
--lines=<lines>         Add marker on pixel lines numbers  (eg. 1200,1200,3000)
--ndv=<ndv>             Use an additionnal no data value
--zoom=<zoom>           Additionnaly display a zoom of the raster ("xmin,xmax,ymin,ymax")
--histo                 Additionnaly display the raster histogram
--stats                 Display the raster and zoom statistics
--save                  Save the display to pdf
--bg                    Path to a background raster of same dimension as infile [must be GDAL raster]
--alpha                 Alpha value to apply to the infile to show the background behind [default: 0.8]
"""

print()
print()
print('Author: Simon Daout')
print()
print('revised version September 2025 (Leo Letellier)')
print()

try:
    from nsbas import docopt
except:
    import docopt

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from osgeo import gdal
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

gdal.UseExceptions()

EXT = {
    'ROIPAC': [
        '.unw',
        '.hgt',
    ],
    'REAL4': [
        '.r4',
    ],
    'GDAL': [
        '.tif',
        '.tiff',
        '.bil',
        '.int',
        '.slc',
        '.flat',
    ],
    'AMSTER': [
    ],
    'GAMMA': [
        '.diff',
    ],
}


def arg2value(value, conversion=None, default=None):
    """Convert a string argument if exists otherwise use default value"""
    if value is None:
        return default
    elif conversion is not None:
        return conversion(value)
    return value


def resolve_format(infile):
    """Resolve format based on extension and automatic parameter file detection"""
    ext = os.path.splitext(infile)[1]
    infile = os.path.abspath(infile)

    file_format = None
    for key, items in EXT.items():
        if ext in items:
            file_format = key
            break

    maybe_real4_param = os.path.join(os.path.dirname(infile), "lect.in")
    maybe_amster_param = os.path.join(os.path.dirname(os.path.dirname(infile)), "TextFiles", "InSARParameters.txt")
    maybe_hdr = os.path.splitext(infile)[0] + '.hdr'
    maybe_rsc = infile + '.rsc'
    has_real4_param = os.path.isfile(maybe_real4_param)
    has_amster_param = os.path.isfile(maybe_amster_param)
    has_hdr = os.path.isfile(maybe_hdr)
    has_rsc = os.path.isfile(maybe_rsc)
    
    if has_hdr:
        return 'GDAL', None
    elif file_format in [None, 'REAL4'] and has_rsc:
        return 'REAL4', maybe_rsc
    elif file_format == 'REAL4' or (file_format is None and has_real4_param):
        return 'REAL4', maybe_real4_param
    elif file_format == 'AMSTER' or (file_format is None and has_amster_param):
        return 'AMSTER', maybe_amster_param
    elif file_format == 'ROIPAC':
        return 'ROIPAC', None
    elif file_format == 'GAMMA':
        raise ValueError('To use GAMMA file please provide the par header file')
    
    try:
        gdal.Open(infile)
    except:
        raise ValueError('Unsupported file')
    return 'GDAL', None


def open_band_gdal(file, band, crop):
    """ Open as GDAL raster
    crop: xmin, xmax, ymin, ymax
    return: data, driver, x, y, b, dtype
    """
    ds = gdal.Open(file)
    band = ds.GetRasterBand(band)
    ndv = band.GetNoDataValue()
    if crop is None:
        crop = [0, band.XSize, 0, band.YSize]
    x_dim = crop[1] - crop[0]
    y_dim = crop[3] - crop[2]
    array = band.ReadAsArray(crop[0], crop[2], x_dim, y_dim)
    if ndv is not None and ndv != np.nan:
        try:
            array[array == ndv] = np.nan
        except:
            pass
    return [array], ds.GetDriver().ShortName, ds.RasterXSize, ds.RasterYSize, ds.RasterCount, gdal.GetDataTypeName(band.DataType)


def open_band_real4(file, band, params, crop, cube):
    """Open as REAL4 raster"""
    fid = open(file, 'r')
    if cube:
        if type(param_file) is list:
            x_dim, y_dim, band_nb = param_file
        else:
            x_dim, y_dim, band_nb = list(map(int, open(params).readline().split(None, 3)[0:3]))
    else:
        band_nb = 1
        if type(param_file) is not str:
            x_dim, y_dim = param_file[:2]
        elif os.path.splitext(params)[1] == '.rsc':
            lines = open(params).read().strip().split('\n')
            x_dim, y_dim = None, None
            for l in lines:
                if 'WIDTH' in l:
                    x_dim = int(''.join(filter(str.isdigit, l)))
                elif 'FILE_LENGTH' in l:
                    y_dim = int(''.join(filter(str.isdigit, l)))
                if x_dim is not None and y_dim is not None:
                    break
        else:
            x_dim, y_dim = list(map(int, open(params).readline().split(None, 2)[0:2]))
    
    phase = np.fromfile(fid, dtype=np.float32)[:y_dim * x_dim * band_nb].reshape((y_dim, x_dim, band_nb))
    phase = phase[:, :, band - 1]

    if crop is not None:
        phase = phase[crop[2]:crop[3], crop[0]:crop[1]]

    data = [phase]
    data_type = np.float32
    driver = 'REAL4'
    return data, driver, x_dim, y_dim, band_nb, data_type


def open_band_roipac(file, crop):
    """Open as custom ROIPAC raster (amplitude / phase)"""
    ds = gdal.OpenEx(file, allowed_drivers=["ROI_PAC"])
    if crop is None:
        crop = [0, ds.RasterXSize, 0, ds.RasterYSize]
    x_dim = crop[1] - crop[0]
    y_dim = crop[3] - crop[2]
    driver = ds.GetDriver().ShortName
    band_nb = ds.RasterCount

    phase_band = ds.GetRasterBand(2)
    phase_ndv = phase_band.GetNoDataValue()
    phase_data = phase_band.ReadAsArray(crop[0], crop[2], x_dim, y_dim)
    amp_band = ds.GetRasterBand(1)
    amp_ndv = amp_band.GetNoDataValue()
    amp_data = amp_band.ReadAsArray(crop[0], crop[2], x_dim, y_dim)
    
    if phase_ndv is not None and phase_ndv != np.nan:
        phase_data[phase_data == phase_ndv] = np.nan
    if amp_ndv is not None and amp_ndv != np.nan:
        amp_data[amp_data == amp_ndv] = np.nan
    data = [phase_data, amp_data]
    data_type = gdal.GetDataTypeName(phase_band.DataType)
    
    return data, driver, ds.RasterXSize, ds.RasterYSize, band_nb, data_type


def open_band_gamma(file, params, crop):
    """Open as GAMMA raster"""
    try:
        from parsers import gamma as gm
    except:
        ModuleNotFoundError("GAMMA parser not found in python installation. Need gamma from module parsers.")

    if params is not None:
        y_dim, x_dim = gm.readpar(par=params)
        phase = gm.readgamma_int(file, par=params)
    else:
        y_dim, x_dim = gm.readpar()
        phase = gm.readgamma_int(file)

    if crop is not None:
        phase = phase[crop[2]:crop[3], crop[0]:crop[1]]
    
    data = [phase]
    driver = 'GAMMA'
    band_nb = 1
    data_type = 'unknown'

    return data, driver, x_dim, y_dim, band_nb, data_type


def open_band_amster(file, params, crop):
    """Open as AMSTer raster"""
    with open(params, 'r') as pfile:
        lines = [''.join(l.strip().split('\t\t')[0]) for l in pfile.readlines()]
        jump_index = lines.index('/* -5- Interferometric products computation */')
        img_dim = lines[jump_index + 2: jump_index + 4]
        y_dim, x_dim = (int(img_dim[1].strip()), int(img_dim[0].strip()))
        band_nb = 1
    
    array = np.fromfile(file, dtype=np.float32)
    data = [array[:x_dim * y_dim].reshape((y_dim, x_dim))]
    driver = 'AMSTer'
    data_type = np.float32

    if crop is not None:
        data[0] = data[0][crop[2]:crop[3], crop[0]:crop[1]]

    return data, driver, x_dim, y_dim, band_nb, data_type


def correct_values_phase(phase, ext, rad2mm, wrap, supp_ndv):
    """Apply corrections to phase values"""
    if supp_ndv is not None:
        phase[phase == supp_ndv] = np.nan

    if rad2mm is not None:
        # scale the values
        phase = phase * rad2mm

    if ext in ['.slc']:
        phase = np.absolute(phase)
    if ext in ['.int', '.flat']:
        phase = np.angle(phase)

    if wrap is not None:
        # simulate wrapped values
        phase = np.mod(phase + wrap, 2 * wrap) - wrap
    
    return phase


def correct_values_amp(amp, ext):
    """Cpply corrections to amplitude values"""
    if ext in ['.int', '.flat', '.diff']:
        amp = np.absolute(amp)
    return amp


def resolve_plot(data, arguments, crop, do_save, bg, alpha):
    """Manage all displays to be plotted"""
    vmin = arg2value(arguments["--vmin"], float)
    vmax = arg2value(arguments["--vmax"], float)
    if (vmax is None) ^ (vmin is None):
        vmin = -vmax if vmax is not None else vmin
        vmax = -vmin if vmin is not None else vmax
    elif (vmax is None and vmin is None):
        vmin = np.nanpercentile(data[0], 2)
        vmax = np.nanpercentile(data[0], 98)
    
    cpt = arguments["--cpt"]
    if cpt is None:
        try:
            from matplotlib.colors import LinearSegmentedColormap
            cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
            cpt = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
            cpt = cpt.reversed()
        except:
            cpt=cm.rainbow

    cols = arguments["--cols"]
    lines = arguments["--lines"]
    cross = None
    if cols is not None and lines is not None:
        cross = [[int(k) for k in cols.split(",")], [int(k) for k in lines.split(",")]]
        if crop is not None:
            cross[0] = [k - crop[0] for k in cross[0]]
            cross[1] = [k - crop[2] for k in cross[1]]

    title = arg2value(arguments["--title"], default=arguments["<infile>"])

    zoom = arguments["--zoom"]
    if zoom is not None:
        zoom = [int(z) for z in zoom.split(',')]
    
    origin = None
    if crop is not None:
        origin = (crop[0], crop[2])

    # Plot the main dislay (phase)
    plot_raster(data[0], cpt, vmin, vmax, cross, title, zoom, origin, bg, alpha)

    if do_save:
        print("Saving figure...")
        plt.savefig(infile + '.pdf', format='PDF', dpi=180)
    
    if len(data) > 1:
        # Plot the secondary display (amplitude)
        vmin = np.nanpercentile(data[1], 2)
        vmax = np.nanpercentile(data[1], 98)
        plot_raster(data[1], 'Greys_r', vmin, vmax, cross, title + " [Amplitude]", zoom, origin, bg, alpha)
        if do_save:
            print("Saving amplitude...")
            plt.savefig(infile + '_amplitude.pdf', format='PDF', dpi=180)

    if arguments["--histo"]:
        # Plot all histograms
        plot_histo(data, title, crop, zoom)
        if do_save:
            print("Saving histo...")
            plt.savefig(infile + '_histo.pdf', format='PDF', dpi=180)

    if zoom is not None:
        plot_zoom(data, crop, zoom, cpt, vmin, vmax, title, bg, alpha)
        if do_save:
            print("Saving zoom...")
            plt.savefig(infile + '_zoom.pdf', format='PDF', dpi=180)
    
    if arguments["--stats"]:
        display_stats(data, zoom, crop)


def plot_raster(raster, cpt, vmin, vmax, cross, title, zoom, origin, bg, alpha):
    """Construct the raster display"""
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1,1,1)
    extent = None
    if origin is not None:
        extent = (origin[0], origin[0] + raster.shape[1], origin[1] + raster.shape[0], origin[1])
    if bg is not None:
        bax = ax.imshow(bg, 'Greys_r', interpolation='nearest', extent=extent)
    hax = ax.imshow(raster, cpt, interpolation='nearest', vmin=vmin, vmax=vmax, extent=extent, alpha=alpha)
    ax.set_title(title)
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(hax, cax=c)

    if cross is not None:
        for i in range(len(cross[0])):
            ax.scatter(cross[0][i], cross[1][i], marker='x', color='black', s=150.0)
    
    if zoom is not None:
        ax.plot([zoom[0], zoom[0], zoom[1], zoom[1], zoom[0]], 
                 [zoom[2], zoom[3], zoom[3], zoom[2], zoom[2]],
                 "-", color='black', linewidth=1)
    
    plt.tight_layout()


def plot_histo(data, title, crop, zoom):
    """Construct the histogram display"""
    fig = plt.figure(figsize=(5, 5))

    histo_data = [data[0]]
    histo_label = ['Main']
    if len(data) > 1:
        histo_data.append(data[1])
        histo_label.append('Secondary')
    if crop is None or len(crop) != 4:
        crop = [0, data[0].shape[0], 0, data[0].shape[1]]
    if zoom is not None:
        histo_data.append(data[0][zoom[2] - crop[2]:zoom[3] - crop[2], zoom[0] - crop[0]:zoom[1] - crop[0]])
        histo_label.append('Zoom')
    
    for d, l in zip(histo_data, histo_label):
        lower = np.nanpercentile(d, 1)
        upper = np.nanpercentile(d, 99)
        hist_values, bin_edges = np.histogram(d[~np.isnan(d)].flatten(), bins=50, range=(lower, upper))
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        plt.plot(bin_centers, hist_values / np.sum(hist_values), label=l)
    plt.title(title + " [Histogram]")
    plt.legend()

    # TODO add a box plot view ?
    # min min2% 


def plot_zoom(data, crop, zoom, cpt, vmin, vmax, title, bg, alpha):
    """Construct the zoom display"""
    if crop is None or len(crop) != 4:
        crop = [0, data[0].shape[0], 0, data[0].shape[1]]
    zdata = data[0][zoom[2] - crop[2]:zoom[3] - crop[2], zoom[0] - crop[0]:zoom[1] - crop[0]]
    if bg is not None:
        zbg = bg[zoom[2] - crop[2]:zoom[3] - crop[2], zoom[0] - crop[0]:zoom[1] - crop[0]]
    else:
        zbg = None
    plot_raster(zdata, cpt, vmin, vmax, None, title + " [ZOOM]", None, (zoom[0], zoom[2]), zbg, alpha)


def plot_vario(data, dims, model='spherical', samples=10000, maxlag=400, n_lags=20):
    print(">  Variogram:")
    rows, cols = np.indices((data.shape[0], data.shape[1]))
    data, rows, cols = data.flatten(), rows.flatten() * dims[1], cols.flatten() * dims[0]
    valid_data = ~np.isnan(data)
    data = data[valid_data]
    rows = rows[valid_data]
    cols = cols[valid_data]
    v = Variogram(np.column_stack((cols , rows)), data, model=model, use_nugget=True, samples=samples, maxlag=maxlag, n_lags=n_lags)
    v.plot()
    stats = v.parameters
    print("Range:\t{:.4}".format(stats[0]))
    print("Sill:\t{:.4}".format(stats[1]))
    print("Nugget:\t{:.4}".format(stats[2]))


def display_raster_format(infile, driver, x, y, b, dtype):
    """Display information about the raster reading"""
    print(">  File:", infile)
    print(">  Driver:", driver)
    print(">  Size:", x, y, b)
    print(">  DataType:", dtype)


def display_stats(data, zoom, crop):
    from scipy.stats import describe

    STATS = ['MIN', 'MAX', 'MEAN', 'VAR', 'MED', 'SKW', 'KRT', 'VAL']
    
    print(">  Stats:")
    # print("\tMin\tMax\tMean\tVariance\tMedian\tSkewness\tKurtosis\tValid", end='\n\n')
    # print("Main:")
    desc = describe(data[0], axis=None, nan_policy="omit")
    med = np.nanmedian(data[0])
    nb_values = data[0].shape[0] * data[0].shape[1]
    nb_nans = np.count_nonzero(np.isnan(data[0]))
    region = ['MAIN']
    stats = [['{:.4}'.format(desc[1][0]), '{:.4}'.format(desc[1][1]), '{:.4}'.format(desc[2]), '{:.4}'.format(desc[3]), 
              '{:.4}'.format(med), '{:.4}'.format(desc[4]), '{:.4}'.format(desc[5]), '{:.3}%'.format((1 - nb_nans/nb_values) * 100)]]
    # stats = [[desc[1][0], desc[1][1], desc[2], desc[3], med, desc[4], desc[5], (1 - nb_nans/nb_values) * 100]]
    # print(f"\t{desc[1][0]}\t{desc[1][1]}\t{desc[2]}\t{desc[3]}\t{med}\t{desc[4]}\t{desc[5]}\t{1 - nb_nans/nb_values}", end='\n\n')

    if len(data) > 1:
        # print("Secondary: ")
        desc = describe(data[1], axis=None, nan_policy="omit")
        med = np.nanmedian(data[1])
        nb_values = data[1].shape[0] * data[1].shape[1]
        nb_nans = np.count_nonzero(np.isnan(data[1]))
        region.append('SECOND')
        stats.append(['{:.4}'.format(desc[1][0]), '{:.4}'.format(desc[1][1]), '{:.4}'.format(desc[2]), '{:.4}'.format(desc[3]), 
              '{:.4}'.format(med), '{:.4}'.format(desc[4]), '{:.4}'.format(desc[5]), '{:.3}%'.format((1 - nb_nans/nb_values) * 100)])
        # print(f"\t{desc[1][0]}\t{desc[1][1]}\t{desc[2]}\t{desc[3]}\t{med}\t{desc[4]}\t{desc[5]}\t{1 - nb_nans/nb_values}", end='\n\n')
    
    if zoom is not None:
        if crop is None or len(crop) != 4:
            crop = [0, data[0].shape[0], 0, data[0].shape[1]]
        # print("Zoom: ")
        zdata = data[0][zoom[2] - crop[2]:zoom[3] - crop[2], zoom[0] - crop[0]:zoom[1] - crop[0]]
        desc = describe(zdata, axis=None, nan_policy="omit")
        med = np.nanmedian(zdata)
        nb_values = zdata.shape[0] * zdata.shape[1]
        nb_nans = np.count_nonzero(np.isnan(zdata))
        region.append('ZOOM')
        stats.append(['{:.4}'.format(desc[1][0]), '{:.4}'.format(desc[1][1]), '{:.4}'.format(desc[2]), '{:.4}'.format(desc[3]), 
              '{:.4}'.format(med), '{:.4}'.format(desc[4]), '{:.4}'.format(desc[5]), '{:.3}%'.format((1 - nb_nans/nb_values) * 100)])
        # print(f"\t{desc[1][0]}\t{desc[1][1]}\t{desc[2]}\t{desc[3]}\t{med}\t{desc[4]}\t{desc[5]}\t{1 - nb_nans/nb_values}", end='\n\n')
    
    print('\t' + '\t\t'.join(region))
    for i, s in enumerate(STATS):
        print(s + '\t', end='')
        for r in range(len(region)):
            print(stats[r][i], end='\t\t')
        print()
    

if __name__ == "__main__":
    arguments = docopt.docopt(__doc__)
    infile = arguments["<infile>"]
    if not os.path.isfile(infile):
        raise FileNotFoundError("No such raster file: {}".format(infile))
    ext = os.path.splitext(infile)[1]
    if arguments["--vario"]:
        try:
            from skgstat import Variogram
        except:
            raise ImportError("Need module scikit-gstats for variogram features")

    crop = arguments["--crop"]
    if crop is not None:
        crop = [int(k) for k in crop.split(",")]
    
    roicube = False

    if arguments["--dim"] is not None:
        file_format = 'REAL4'
        param_file = [int(d) for d in arguments["--dim"].split(',', 3)]
        if len(param_file) == 2:
            param_file.append(1)
    elif arguments["--gdal"]:
        file_format = 'GDAL'
        param_file = None
    elif arguments["--lectfile"] is not None:
        file_format = 'REAL4'
        param_file = arguments["--lectfile"]
    elif arguments["--lectcube"] is not None:
        file_format = 'REAL4'
        param_file = arguments["--lectcube"]
        roicube = True
    elif arguments["--parfile"] is not None:
        file_format = 'GAMMA'
        param_file = arguments["--parfile"]
    elif arguments["--amfile"] is not None:
        file_format = 'AMSTER'
        param_file = arguments["--amfile"]
    else:
        file_format = None
        param_file = None

    if type(param_file) is not list:
        if param_file is not None and not os.path.isfile(param_file):
            raise FileNotFoundError("No such parameter file: {}".format(param_file))

    if file_format is None:
        file_format, param_file = resolve_format(infile)

    band = arg2value(arguments["--band"], int, 1)
    supp_ndv = arg2value(arguments["--ndv"], float)
    
    if file_format == 'REAL4':
        data, driver, x, y, b, dtype = open_band_real4(infile, band, param_file, crop, cube=roicube)
    elif file_format == 'GDAL' or (file_format == 'ROIPAC' and arguments["--band"] is not None):
        data, driver, x, y, b, dtype = open_band_gdal(infile, band, crop)
    elif file_format == 'ROIPAC':
        data, driver, x, y, b, dtype = open_band_roipac(infile, crop)
    elif file_format == 'AMSTER':
        data, driver, x, y, b, dtype = open_band_amster(infile, param_file, crop)
    elif file_format == 'GAMMA':
        data, driver, x, y, b, dtype = open_band_gamma(infile, param_file, crop)

    display_raster_format(infile, driver, x, y, b, dtype)

    bg = arguments["--bg"]
    alpha = arg2value(arguments["--alpha"], float, 0.8)
    if bg is not None:
        try:
            bg = open_band_gdal(bg, 1, crop)[0][0]
        except:
            raise ValueError('Background file is not valid: {}'.format(bg))
    else:
        alpha = 1

    rad2mm = arg2value(arguments["--rad2mm"], float)
    wrap = arg2value(arguments["--wrap"], float)

    data[0] = correct_values_phase(data[0], ext, rad2mm, wrap, supp_ndv)
    if len(data) > 1:
        data[1] = correct_values_amp(data[1], ext)

    do_save = arguments["--save"]
    
    resolve_plot(data, arguments, crop, do_save, bg, alpha)

    if arguments["--vario"]:
        samples = arg2value(arguments["--samples"], conversion=int, default=10000)
        dlag = arg2value(arguments["--dlag"], conversion=int, default=500)
        nlag = arg2value(arguments["--nlag"], conversion=int, default=20)
        model = arg2value(arguments["--model"], default="spherical")
        res = arg2value(arguments["--res"], conversion=lambda x: x.split(',')[:2], default=[1, 1])
        plot_vario(data[0], [float(res[0]), float(res[1])], model, samples ,dlag, nlag)

    plt.show()
    