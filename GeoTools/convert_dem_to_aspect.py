#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy.ndimage
from osgeo import gdal, osr
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import sys
import getopt

def compute_slope_aspect(dem_path, filter_size=2.0, plot=True):
    """ 
    Compute slope and aspect from a DEM file.

    Args:
        dem_path (str): Path to the DEM file.
        filter_size (float): Gaussian filter size (default=2.0).
        plot (bool): Whether to display plots.

    Returns:
        tuple: (-aspect, slope)
    """

    if not os.path.isfile(dem_path):
        print(f"Error: File '{dem_path}' not found.")
        sys.exit(1)

    # Load DEM
    ds = gdal.Open(dem_path, gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    topo = band.ReadAsArray().astype(float)
    ncols, nlines = ds.RasterXSize, ds.RasterYSize
    gt = ds.GetGeoTransform()
    projref = ds.GetProjection()
    drv = gdal.GetDriverByName('GTiff')

    # Filter to smooth data and remove noise
    filtered_topo = scipy.ndimage.gaussian_filter(topo, sigma=(filter_size, filter_size))

    # Get middle latitude for geospatial calculations
    lats = gt[3] + (np.arange(nlines) * gt[5])
    lat_mean = np.mean(lats)
    print("GeoTransform:", gt)
    print("Average latitude:", lat_mean)

    # Determine resolution based on coordinate system
    res_x, res_y = gt[1], gt[5]
    srs = osr.SpatialReference(wkt=projref)

    if srs.IsGeographic():
        res_x *= 40075e3 / 360  # Convert degrees to meters
        res_y *= (40075e3 / 360) * np.cos(np.deg2rad(lat_mean))
        print("Geographic coordinates (converted to meters)")
    else:
        print("Projected coordinates (meters)")

    print(f"Resolution: dx={res_x:.2f}, dy={res_y:.2f}")

    # Compute gradient (first axis = Y, second axis = X)
    dsm_dy, dsm_dx = np.gradient(filtered_topo, res_y, res_x)

    # Calculate slope and aspect
    slope = np.sqrt(dsm_dx**2 + dsm_dy**2)
    aspect = np.rad2deg(np.arctan2(dsm_dy, -dsm_dx))  # clockwise from east

    crop_slice = (slice(450, 850), slice(500, 900))
    aspect_crop = aspect[crop_slice] 
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    cax= ax.imshow(aspect_crop, cmap='Greys_r', origin='upper', alpha=0.5, vmax=180, vmin=-180)
    ax.set_xticks([])
    ax.set_yticks([])
    divider = make_axes_locatable(ax)
    c = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=c) 
    #plt.show()
    #sys.exit()

    # Save outputs to GeoTIFF
    def save_raster(filename, array):
        dst = drv.Create(filename, ncols, nlines, 1, gdal.GDT_Float32)
        dst.SetGeoTransform(gt)
        dst.SetProjection(projref)
        dst.GetRasterBand(1).WriteArray(array)
        dst.FlushCache()
        print(f"Saved: {filename}")

    save_raster('slope.tif', np.rad2deg(slope))
    save_raster('aspect.tif', aspect)

    # Plot results
    if plot:
        fig, axes = plt.subplots(2, 2, figsize=(11, 7))
        titles = ['Slope', 'Gradient X', 'Gradient Y', 'Aspect']
        datasets = [slope[450:850, 500:900], dsm_dx[450:850, 500:900], dsm_dy[450:850, 500:900], aspect[450:850, 500:900]]
        cmaps = ['Greys_r', 'coolwarm', 'coolwarm', 'Greys_r']

        for ax, title, data, cmap in zip(axes.flatten(), titles, datasets, cmaps):
            im = ax.imshow(data, cmap=cmap, origin='upper', vmax=np.percentile(data,90), vmin=np.percentile(data,10))
            ax.set_title(title)
            divider = make_axes_locatable(ax)
            cbar = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(im, cax=cbar)

        fig.savefig('dem_slope_aspect.png', dpi=300)
        plt.show()

    return -aspect, slope

def usage():
    print("Usage: convert_dem_to_aspect.py demfile [-v] [-h]")
    print("-v Verbose mode. Show more information")
    print("-h Show this help message")

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hv", ["help", "verbose"])
    except getopt.GetoptError:
        print("Error in arguments. Use -h for help.")
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-v", "--verbose"):
            print("Verbose mode enabled")

    if len(args) < 1:
        print("Error: No DEM file provided.")
        usage()
        sys.exit(1)

    demfile = args[0]
    compute_slope_aspect(demfile, plot=True)

