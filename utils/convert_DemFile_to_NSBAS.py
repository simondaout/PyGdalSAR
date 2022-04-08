#!/usr/bin/env python3

"""
Build the DEM for use in NSBAS:

* remove the geoid
* convert it to the roi_pac format

You should have access to the internet for this script to work.
Heavily rely on gdal
"""

import subprocess
import re
import os
import sys
import logging
from pathlib import Path

import nsbas.nsbio as nsbio
import dem.fetch_raw as demfetch
if "GDAL_SYS" in os.environ and os.environ["GDAL_SYS"] == 'True':
    from osgeo import gdal
else:
    import nsbas.gdal as gdal

logger = logging.getLogger(__name__)

SRTM_DEM_CACHE = 'SRTM_DEM_CACHE'

# ############## Definition of functions ###################


def clean_dem_dir(directory):
    """
    clean directory, removing temporal files
    """
    for fic in os.listdir(directory):
        abs_file = os.path.join(directory, fic)
        logger.debug("checking %s", abs_file)
        if re.search(r"(\.tiff|\.zip|\.hgt)$", abs_file):
            logger.debug("removing %s", abs_file)
            os.unlink(os.path.join(directory, abs_file))


def clean_and_die(err_num, err_str, directory):
    """
    clean directory and exit

    :param err_num: the error number to exit with
    :type err_num: int
    :param err_str: the message to print before exiting
    :type err_str: str
    :param directory: the directory to clean
    :type directory: str
    """
    clean_dem_dir(directory)
    if err_num == 0:
        logger.info("quitting %s", err_str)
    else:
        logger.critical("aborting: %s", err_str)
    sys.exit(err_num)


def sample_srtm(currdir, in_file, out_file):
    """
    The downloaded SRTM tiles should be referenced to the WGS84 ellipsoid.
    If not, this step allows to assign elevation values where information
    is lacking.
    :param in_file: merged tiff srtm
    :param out_file: resampled tiff srtm
    :type in_file: str
    :type out_file:str
    """

    gdalwarp = "gdalwarp"
    try:
        logger.info('Resample merged tiles.')
        proc = subprocess.Popen([gdalwarp, "-t_srs",
                                 "+proj=longlat +datum=WGS84 +no_defs",
                                 "-srcnodata", "-32768",
                                 "-dstnodata", "-32768",
                                 os.path.join(currdir, in_file),
                                 out_file],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate()
        exitcode = proc.returncode

    except OSError:
        raise OSError("{} command not found".format(gdalwarp))
    if exitcode != 0:
        raise RuntimeError(("command FAILED\nEXIT CODE={exitcode}\n"
                            f"stdout={out}\nstderr={err}"))
    logger.debug(f"result of {gdalwarp} is: {out_file}")


def get_geographic_coordinates(currdir, in_file):
    """
    Get the geographic coordinates of the created tiff file.
    :param in_file: resampled tiff srtm
    :type: str
    """
    src = gdal.Open(os.path.join(currdir, in_file))
    ulx, xres, _, uly, _, yres = src.GetGeoTransform()
    lrx = ulx + (src.RasterXSize * xres)
    lry = uly + (src.RasterYSize * yres)
    coordinates = [ulx, lry, lrx, uly]
    del src
    return coordinates


def get_latLong_coordinates(in_file):
    """
    Get the lattitude and longitude coordinates of the radar data.
    :param in_file: path to tiff radar data
    :type: str
    """
    # en test ...
    logger.debug("get_lat_long_coordinates: checking %s", in_file)
    src = gdal.Open(in_file)
    ground_control_points = src.GetGCPs()
    p0 = ground_control_points[0]
    plast = ground_control_points[-1]
    del src
    return [p0.GCPX, p0.GCPY, plast.GCPX, plast.GCPY]


def crop_egm96(geo_coordinates, currdir,egm_file, out_file, xres, yres):
    """
    Extract the geoid undulations for the study area from a global map.
    :param geo_coordinates: geocoordinates of the resampled tiff file
    :param egm_file: geoid file
    :param out_file: cropped geoid file
    :type geo_coordinates: array of str
    :type egm_file: str
    :type out_file:str
    """
    str_geo_coordinates = [str(coord) for coord in geo_coordinates]
    if os.path.exists(out_file):
        logger.info("crop egm96: file %s exist, remove it to be gdalwrap compliant", out_file)
        os.unlink(out_file)
    try:
        logger.info('Extract the geoid undulations for the study area from a global map.')
        proc = subprocess.Popen(["gdalwarp", "-t_srs",
                                 "+proj=longlat +datum=WGS84 +no_defs", "-te",
                                 str_geo_coordinates[0], str_geo_coordinates[1],
                                 str_geo_coordinates[2], str_geo_coordinates[3],
                                 "-tr", str(xres), str(yres),
                                 "-r", "cubic",
                                 os.path.join(currdir, egm_file),
                                 out_file],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate()
        exitcode = proc.returncode
    except OSError:
        raise OSError("{} command not found".format("gdalwarpd"))
    if exitcode != 0:
        raise RuntimeError((f"command failed\nExit code={exitcode}\n"
                            f"stdout={out}\nstderr={err}"))
    logger.debug(f"result of crop_egm96 is: {out_file}")


def add_egm96_to_dem(egm_file, currdir, in_file, out_file):
    """
    Add the EGM96 to the tiff DEM file.
    :param egm_file: cropped geoid file
    :param in_file: resampled tiff srtm
    :param out_file: cropped tiff srtm
    :type egm_file: str
    :type in_file: str
    :type out_file:str
    """

    gdalcalculate = "gdal_calc.py"

    try:
        logger.info('Add geoid to the resampled file.')
        proc = subprocess.Popen([gdalcalculate,
                                 "-A",
                                 os.path.join(currdir, in_file),
                                 "-B",
                                 os.path.join(currdir, egm_file),
                                 "--outfile",
                                 out_file,
                                 "--calc",
                                 "A+B"],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate()
        exitcode = proc.returncode
    except OSError:
        raise OSError(f"{gdalcalculate}: command not found")
    if exitcode != 0:
        raise RuntimeError((f"command failed\nExit code={exitcode}\n"
                            f"stdout={out}\nstderr={err}"))
    logger.debug(f"result of {gdalcalculate} is: {out_file}")


def convert_to_dem(currdir, in_file, out_file):
    """
    Convert the .tiff file to a .dem file readable by NSBAS/ROI_PAC.
    :param in_file: cropped tiff srtm
    :param out_file: dem file readable by NSBAS/ROI_PAC
    :type in_file:str
    :type out_file: str
    """

    gdaltranslate = "gdal_translate"

    try:
        logger.info('Convert the tiff file to a file in NSBAS/ROI_PAC format')
        proc = subprocess.Popen(["gdal_translate", "-of",
                                 "roi_pac", "-ot", "int16",
                                 os.path.join(currdir, in_file),
                                 out_file],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate()
        exitcode = proc.returncode
    except OSError:
        raise OSError("{} command not found".format(gdaltranslate))
    if exitcode != 0:
        raise RuntimeError((f"command failed\nExit code={exitcode}\n"
                            f"stdout={out}\nstderr={err}"))

    logger.debug("result of {} is: {}".format(gdaltranslate, out_file))


def file_check(file_name, file_min_size):
    """
    :param file_name: the file name to check
    :type file_name: str
    :param file_min_size: the minimal file size (default 0,
           ie file_size should be >=0)
    :type file_min_size: int
    :return: error message, or None if ok
    :rtype: str or None
    """
    if not os.path.exists(file_name):
        return 1, f"file {file_name} does not exist"
    if os.stat(file_name).st_size < file_min_size:
        return 2, (f"file {file_name} exist but has a "
                   f"weird size ({os.stat(file_name).st_size})")
    return None


def path_check(path_check):
    """
    Check if path starts with "/"
    :return: string error code
    :rtype: an error string if not starting with /, None else
    """
    if not path_check.startswith('/'):
        return "Error: {} is not an absolute path!".format(path_check)
    return None


if __name__ == "__main__":

    import argparse
    # =====================
    # Parse arguments
    parser = argparse.ArgumentParser(description=("Download the srtm 30m tiles"
                                                  ", merge and convert them "
                                                  "to .dem format"))
     
    parser.add_argument("-t", type=str, default=".",
                        help=("defines the target dir, ie the working dir. "
                              "Target dir is created if it does not exist"))

    parser.add_argument("--dem", type=str,
                          help=("DEM file in EPSG:4326"
                                "set nodata to -32768"))    
    parser.add_argument("-v", type=int, default=3,
                        help=("set logging level: 0 critical, 1 error, "
                              "2 warning, 3 info, 4 debug, default=info"))
    
    args = parser.parse_args()
    merged_tiles_file_name = args.dem 
    dirname=args.t

    file_prefix = os.path.splitext(merged_tiles_file_name)[0]
    percent = 0.0
    step_percent = 100.0/4.0
    keep= False
    
    logging_translate = [logging.CRITICAL, logging.ERROR, logging.WARNING,
                         logging.INFO, logging.DEBUG]
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                        level=logging_translate[args.v])
    logger = logging.getLogger(__name__)
    logger.info("General start")

    # =====================
    # Sample srtm
    # =====================

    logger.info('start resampling srtm')
    resample_file_name = file_prefix + "_resampled.tiff"
    try:
        sample_srtm(dirname, merged_tiles_file_name, resample_file_name)
    except Exception as e:
        logger.critical("cannot sample srtm: {}".format(str(e)))
        sys.exit(2)
    err = file_check(resample_file_name, 0)
    if err is not None and not keep:
        logger.critical("file checking failed, aborting")
        clean_and_die(2, err[1], args.t)

    percent += step_percent
    logger.info("percent %d", percent)
    logger.info('end resampling srtm')

    # =====================
    # Get geocoordinates
    # =====================
    logger.info('start geting geocoordinates')

    geo_coordinates = get_geographic_coordinates(dirname, resample_file_name)
    logger.info('geting geocoordinates done')
    percent += step_percent
    logger.info('end geting geocoordinates')
    logger.info("percent %d", percent)
    
    ds = gdal.Open(resample_file_name, gdal.GA_ReadOnly)
    xres = ds.GetGeoTransform()[1]
    yres =  ds.GetGeoTransform()[5]

    # =====================
    # Crop geoid
    # =====================
    logger.info("start cropping geoid")
    egm_file = dirname + "/egm96-5-flipped.tiff"
    if not os.path.exists(egm_file):
        # backoff to env variable
        flipped_dem = None
        if "NSBAS" in os.environ and os.path.exists(os.environ["NSBAS"] + '/eap/EGM96/egm96-5-flipped.tiff'):
            flipped_dem = os.environ["NSBAS"] + '/eap/EGM96/egm96-5-flipped.tiff'
        elif "EGM96FLIPPED" in os.environ and os.path.exists(os.environ["EGM96FLIPPED"]):
            flipped_dem = os.environ["EGM96FLIPPED"]
        else:
            logger.critical(("cannot find egm96-5-flipped.tiff, "
                             "neither in %s nor in as env variables", dirname))
            clean_and_die(3, err[1], dirname)
        os.symlink(flipped_dem, egm_file)
    cropped_egm = os.path.join(dirname, "egm96-5-flipped-cropped.tiff")
    try:
        crop_egm96(geo_coordinates, dirname, egm_file, cropped_egm, xres, yres)
    except Exception as excpt:
        logger.critical("cannot crop egm96: %s", str(excpt))
        clean_and_die(3, err[1], dirname)
    err = file_check(cropped_egm, 0)
    if err is not None and not keep:
        clean_and_die(3, err[1], dirname)
    percent += step_percent
    logger.info("percent %d", percent)
    logger.info("end cropping geoid")

    # =====================
    # Add geoid
    # =====================
    logger.info("start adding geoid")

    cropped_file_name = file_prefix + "_cropped.tiff"
    try:
        add_egm96_to_dem(cropped_egm, dirname, resample_file_name,
                         cropped_file_name)
    except Exception as excpt:
        logger.critical("cannot add egm96: %s", str(excpt))
        clean_and_die(4, err[1], dirname)
    err = file_check(cropped_file_name, 0)
    if err is not None and not keep:
        clean_and_die(4, err[1], dirname)
    percent += step_percent
    logger.info("percent %d", percent)
    logger.info("end adding geoid")

    # =====================
    # Convert to dem
    # =====================
    logger.info("start convert to dem")
    dem_file_name = file_prefix+".dem"
    try:
        convert_to_dem(dirname, cropped_file_name, dem_file_name)
    except Exception as excpt:
        logger.critical("cannot add egm86: %s", str(excpt))
        sys.exit(5)
    rscDemFile = open((dem_file_name+'.rsc'), "r")
    found = 0
    for line in rscDemFile:
        if line.startswith("WIDTH"):
            width = int(line.split()[1])
            found += 1
            if found == 3:
                break
        if line.startswith("FILE_LENGTH"):
            length = int(line.split()[1])
            found += 2
            if found == 3:
                break
    err = file_check(dem_file_name, width*length*2)
    if err is not None and not keep:
        clean_and_die(5, err[1], dirname)
    percent += step_percent
    logger.info("percent %d", percent)
    logger.info("end convert to dem")

    # =======================
    # cleaning dem dir
    # =======================
    if not keep:
        logger.warning("cleaning working dir %s", dirname)
        clean_dem_dir(dirname)
    percent = 100.0
    logger.info("percent %d", percent)
