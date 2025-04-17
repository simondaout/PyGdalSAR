#!/usr/bin/env python
# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""\
filter_cube.py
-------------
Clean a time series cube with a HP or LP filter band by band

Usage: filter_cube.py --infile=<path> [--outfile=<path>] [--plot=<yes/no>] [--filter=<HP/LP>] [--fwindsize=<value>] [--mask=<path>] [--threshold=<value>] \

Options:
-h --help           Show this screen.
--infile PATH       Time series cube to filter
--outfile PATH      Output file
--filter=<HP/LP>    Apply a high pass (HP) or a low pass (LP) gaussian filter to the image [default: HP]
--fwindsize=<value> Filter window size [default: 30]
--plot=<yes/no>     Plot intermediate result [default:no]
--mask PATH         File used as mask before filtering (e.g. RMSpixel)
--threshold VALUE   Threshold value on mask file (Keep pixel with mask > threshold) [default: 2] 
"""

print()
print()
print('Author: Hugo WATINE /  Simon DAOUT')
print()
print()

import time 
from scipy import ndimage
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from osgeo import gdal
import docopt
import os

# read arguments
arguments = docopt.docopt(__doc__)
if arguments["--fwindsize"] == None:
    arguments["--fwindsize"] = 30
if arguments["--filter"] == None:
    arguments["--filter"] = 'HP'
if arguments["--plot"] == None:
    arguments["--plot"] = 'no'
if arguments["--mask"] ==  None:
    maskf = 'no'
else:
    maskf = arguments["--mask"]
if arguments["--threshold"] ==  None:
    arguments["--threshold"] = 2.
else:
    arguments["--threshold"] = float(arguments["--threshold"])
if arguments["--outfile"] == None:
    basename = os.path.splitext(arguments["--infile"])[0]
    arguments["--outfile"] = basename + '_' + arguments["--filter"] + arguments["--fwindsize"]    

# mask
if maskf != 'no':
  ds_extension = os.path.splitext(maskf)[1]
  if (ds_extension == ".tif" or ds_extension ==".tiff" or ds_extension ==".grd"):
    from osgeo import gdal
    ds = gdal.Open(maskf, gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    # Attributes
    print("> Driver:   ", ds.GetDriver().ShortName)
    print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
    print("> Datatype: ", gdal.GetDataTypeName(band.DataType))
    mask = band.ReadAsArray(0, 0,
                               ds.RasterXSize, ds.RasterYSize,
                               ds.RasterXSize, ds.RasterYSize)
  else:
    ncols, nlines = list(map(int, open('lect.in').readline().split(None, 2)[0:2]))
    fid2 = open(maskf, 'r')
    mask = np.fromfile(fid2,dtype=np.float32)[:nlines*ncols].reshape((nlines,ncols))
    mask =  mask
    mask[np.isnan(mask)] = 0
else:
    mask = np.zeros((nlines,ncols))*float('nan')

def clean_raster(band, arguments, nodata):
    if arguments["--filter"] == 'HP':
        m_filter = np.copy(band)
        sum_coef = 0*np.copy(band) +1 # The goal is to take into acount nodata value in the filter

        index = (band == nodata) | np.isnan(band)
        
        m_filter[index] = 0.
        sum_coef[index] = 0.

        band = band - ndimage.gaussian_filter(m_filter, int(arguments["--fwindsize"]))/ndimage.gaussian_filter(sum_coef, int(arguments["--fwindsize"]))

        band[index] = float('nan')
    
    elif arguments["--filter"] == 'LP':
        m_filter = np.copy(band)
        index = np.isnan(band)
        m_filter[index] = 0.
        band = ndimage.gaussian_filter(m_filter, int(arguments["--fwindsize"]))
        band[index] = float('nan')
   
    return band

# Ouvrir le fichier raster cube
dataset = gdal.Open(arguments["--infile"], gdal.GA_ReadOnly)
if dataset is None:
    print(f"Erreur lors de l'ouverture du fichier {arguments['--infile']}")

print("> Driver:   ", dataset.GetDriver().ShortName)
print("> Size:     ", dataset.RasterXSize,'x',dataset.RasterYSize,'x',dataset.RasterCount)

# Extraire les dimensions du raster
ncols = dataset.RasterXSize
nrows = dataset.RasterYSize
nbands = dataset.RasterCount
gt = dataset.GetGeoTransform()
proj = dataset.GetProjectionRef()
driver = gdal.GetDriverByName('GTiff')

# Créer un nouveau cube vide pour stocker les bandes traitées
processed_cube = np.zeros((nbands, nrows, ncols))

# Traiter chaque bande individuellement
for band_index in range(1, nbands + 1):
    print(f"Traitement de la bande {band_index}/{nbands}...")
    start_time= time.time()
    # Lire la bande actuelle
    
    band = dataset.GetRasterBand(band_index)
    nodata_value = band.GetNoDataValue()
    band = band.ReadAsArray(0, 0, ncols, nrows, ncols, nrows)
    # replace 9999 by NaN
    kk = np.nonzero(np.logical_or(band==9990, band==9999))
    band[kk] = float('NaN')
    # clean based on mask
    if maskf != 'no':
        band[mask>arguments["--threshold"]] = float('NaN')        

    # Appliquer un traitement à la bande (par exemple, nettoyer ou filtrer)
    processed_band = clean_raster(band, arguments, nodata_value)
    
    if arguments["--plot"] == 'yes' and band_index % 5 == 0:
        fig = plt.figure(0, figsize=(12,6))    
        ax = fig.add_subplot(1,2,1)

        try:
            import matplotlib.cm as cm
            from matplotlib.colors import LinearSegmentedColormap
            cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
            cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
            cmap = cmap.reversed()
        except:
            cmap=cm.rainbow
        
        # Afficher la bande non filtrée
        vmax, vmin = np.nanpercentile(band,98), np.nanpercentile(band,2)
        cax = ax.imshow(band, cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
        ax.set_title('Band')
        
        # Ajouter une colorbar
        divider = make_axes_locatable(ax)
        c = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax, cax=c)

        ax = fig.add_subplot(1,2,2)
        # Afficher la bande filtrée
        #vmax, vmin = np.nanpercentile(processed_band,98), np.nanpercentile(band,2)
        cax = ax.imshow(processed_band, cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
        ax.set_title('Filter Band')
        
        # Masquer les labels des ticks x
        plt.setp(ax.get_xticklabels(), visible=False)

        # Ajouter une colorbar
        divider = make_axes_locatable(ax)
        c = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax, cax=c)

        # Afficher la figure
        plt.show()
 
    # Ajouter la bande traitée au cube final
    processed_cube[band_index - 1, :, :] = processed_band
    
    elapsed_time = time.time() - start_time
    print(f"Temps de traitement pour la bande : {elapsed_time:.2f} secondes")

# Réassembler le cube et l'enregistrer
driver = gdal.GetDriverByName('GTiff')  # Choisir le format de sortie (ex: GeoTIFF)
out_dataset = driver.Create(arguments["--outfile"], ncols, nrows, nbands, gdal.GDT_Float32)
out_dataset.SetGeoTransform(gt)
out_dataset.SetProjection(proj)

# Sauvegarder chaque bande traitée dans le fichier de sortie
for band_index in range(1, nbands + 1):
    out_band = out_dataset.GetRasterBand(band_index)
    out_band.WriteArray(processed_cube[band_index - 1, :, :])
    band_metadata = dataset.GetRasterBand(band_index).GetMetadata()
    if band_metadata:
        out_band.SetMetadata(band_metadata)
# Fermer les fichiers
out_dataset.FlushCache()
dataset = None
out_dataset = None
print(f"Cube traité et enregistré dans {arguments['--outfile']}")
