#!/usr/bin/env python3

"""
Build the DEM for use in NSBAS:

* download the tiles
* merge them
* do resampling
* cropping to the bounding box
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

gdal.UseExceptions()

logger = logging.getLogger(__name__)

SRTM_DEM_CACHE = 'SRTM_DEM_CACHE'
COP_DEM_CACHE = 'COP_DEM_CACHE'

# ############## Definition of functions ###################


def clean_dem_dir(directory):
    """
    clean directory, removing temporal files
    """
    for fic in os.listdir(directory):
        abs_file = os.path.join(directory, fic)
        logger.debug("checking %s", abs_file)
        #if re.search(r"(\.tiff|\.zip|\.hgt|\.dt2|\.tif)$", abs_file):
        if re.search(r"(\.hgt|\.dt2)$", abs_file):
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


def merge_tiles(extension, currdir, out_file):
    """
    Merge the SRTM tiles with .hgt extension to a tiff file
    :param extension: the extension of the downloaded SRTM tiles
    :param currdir: the working directory containing "extension" files
    :param out_file: merged tiff file
    :type extension: str
    :type currdir:str
    :type out_file: str
    """
    logger.debug('merging tiles')

    gdalmerge = "gdal_merge.py"
    try:
        logger.info('Merge tiles.')
        cmd_array = [gdalmerge, "-o", out_file, "-of", "GTIFF"]
        file_array = [f"{currdir}/{f}" for f in os.listdir(currdir)
                      if f.endswith(extension)]
        proc = subprocess.Popen(cmd_array + file_array,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate()
        exitcode = proc.returncode

    except OSError:
        err_msg = "{} command not found".format(gdalmerge)
        logger.error(err_msg)
        raise OSError(err_msg)
    if exitcode != 0:
        err_msg = (f"command failed \n EXIT CODE={exitcode}\n"
                   f"stdout={out}\nstderr={err}")
        logger.error(err_msg)
        raise RuntimeError(err_msg)
    logger.debug("result of %s is %s", gdalmerge, out_file)
    logger.debug('merging tiles done')


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

    logger.debug(f"sample_srtm opening {in_file}")
    gdalwarp = "gdalwarp"
    try:
        logger.info('Resample merged tiles.')
        proc = subprocess.Popen([gdalwarp, "-t_srs",
                                 "+proj=longlat +datum=WGS84 +no_defs",
                                 "-srcnodata", "-32768",
                                 "-dstnodata", "-32768",
                                 "-tr", "0.000277778", "-0.000277778", "-r", "cubic",
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
    logger.debug(f"get_lat_long_coordinates: checking {in_file}")
    src = gdal.Open(in_file)
    ground_control_points = src.GetGCPs()
    p0 = ground_control_points[0]
    plast = ground_control_points[-1]
    del src
    return [p0.GCPX, p0.GCPY, plast.GCPX, plast.GCPY]


def crop_egm(geo_coordinates, currdir, egm_file, out_file):
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
        logger.info("crop EGM: file %s exist, remove it to be gdalwrap compliant", out_file)
        os.unlink(out_file)
    try:
        logger.info(('Extract the geoid undulations for the study area from a global map.'
                     f' Filename = {os.path.join(currdir, egm_file)}'))
        proc = subprocess.Popen(["gdalwarp",  "-t_srs",
                                 "+proj=longlat +datum=WGS84 +no_defs", "-te",
                                 str_geo_coordinates[0], str_geo_coordinates[1],
                                 str_geo_coordinates[2], str_geo_coordinates[3],
                                 "-tr", "0.000277778", "-0.000277778", "-r", "cubic",
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


def add_egm_to_dem(egm_file, currdir, in_file, out_file):
    """
    Add the EGM to the tiff DEM file.
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
    epilog_str = f'''Dem cache is first checked in env variable
                    {SRTM_DEM_CACHE} and {COP_DEM_CACHE}. If not found, the -c option
                    is used. If none is available the option -t
                    is used.
                 '''
    # =====================
    # Parse arguments
    parser = argparse.ArgumentParser(description=("Download the SRTM or Copernicus 30m tiles"
                                                  ", merge and convert them "
                                                  "to .dem format"),
                                     epilog=epilog_str)
    parser.add_argument("-v", type=int, default=3,
                        help=("set logging level: 0 critical, 1 error, "
                              "2 warning, 3 info, 4 debug"
                              "default=%(default)s"))
    parser.add_argument("-k", action="store_true",
                        help=("if set, keep tiles already downloaded, "
                              "do not re-download if present in directory. "
                              "default=%(default)s"))
    parser.add_argument("-c", type=str, default=None,
                        help=("cache dir: look target files in cache dir."
                              "if not present in cache dir, download tiles in "
                              "cache dir. If option is not provided, download"
                              "locally"))
    parser.add_argument("-t", type=str, default=".",
                        help=("defines the target dir, ie the working dir. "
                              "Target dir is created if it does not exist."
                              "default=%(default)s"))
    parser.add_argument("-s", type=str, default="SRTMGL1.003",
                        choices=["SRTMGL1.003", "NASADEM_HGT.001", "COP_DEM"],
                        help=("defines version of DEM to be used, "
                              "default=%(default)s"))
    coord_in = parser.add_mutually_exclusive_group(required=True)
    coord_in.add_argument("--bbox", type=str,
                          help=("bounding box of the data, in the form "
                                "minlong, maxlong, minlat, maxlat"))
    coord_in.add_argument("--safe", type=str,
                          help=("directory that contains the images from"
                                "which we extract the bounding box"))
    parser.add_argument("-p", type=str, default="vv",
                          help=("polarization of the images "
                                "can be vv, vh, hv, or hh"))

    # search for the -V flag in order to avoid pb with mandatory arguments

    args = parser.parse_args()
    logging_translate = [logging.CRITICAL, logging.ERROR, logging.WARNING,
                         logging.INFO, logging.DEBUG]
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                        level=logging_translate[args.v])
    logger = logging.getLogger(f"{os.path.basename(__file__).replace('.py','')}:main")
    logger.info("General start")

    # =====================
    # normalize paths
    if not args.t.startswith('/'):
        args.t = os.path.abspath(args.t)
        logger.debug("target dir: %s", args.t)

    if args.c is not None:
        cache_dir = os.path.abspath(args.c)
    else:
        if SRTM_DEM_CACHE in os.environ:
            logging.debug(f"using cache dir from env var: {SRTM_DEM_CACHE}={os.environ[SRTM_DEM_CACHE]}")
            cache_dir = os.path.abspath(os.environ[SRTM_DEM_CACHE])
        else:
            cache_dir = args.t
    logging.debug(f"using cache dir {cache_dir}")
    if not os.path.exists(cache_dir):
        logging.debug(f"create dir {cache_dir}")
        os.makedirs(cache_dir)
    percent = 0.0
    step_percent = 100.0/9.0

    logger.info("start extract coordinates")

    # =====================
    # extracting polarization
    polarization = args.p.lower()
    if polarization not in ["vv", "hh", "vh", "hv"]:
        logger.warning(f"polarization should be either vv, hh, vh, or hv, not {polarization}!"
                       f" we will set polarization to vv!")
        polarization = "vv"
    logger.info(f"polarization set to {polarization}")


    # ======================
    # extracting coordinates
    if args.bbox is not None:
        logger.info("extracting bbox from command line option")
        (minLong, maxLong, minLat, maxLat) = (float(x) for x in args.bbox.split(","))
    if args.safe is not None:
        if not args.safe.startswith('/'):
            args.safe = os.path.abspath(args.safe)
        logger.info("extracting bbox from safes: safe dir: %s", args.safe)
        maxLong = -180
        maxLat = -90
        minLat = 1000000
        minLong = 1000000
        for abs_safe in Path(args.safe).glob(f"*.SAFE/measurement/*{polarization}*.tiff"):
            latlong_coord = get_latLong_coordinates(str(abs_safe))
            logger.info(f"{str(abs_safe)}: {latlong_coord}")
            minLong = min(latlong_coord[0], latlong_coord[2], minLong)
            minLat = min(latlong_coord[1], latlong_coord[3], minLat)
            maxLong = max(latlong_coord[0], latlong_coord[2], maxLong)
            maxLat = max(latlong_coord[1], latlong_coord[3], maxLat)
        logger.info("bounding box (minlong, maxlong, minlat, maxlat) = (%f, %f, %f,%f)",
                    minLong, maxLong, minLat, maxLat)
        if minLat == 1000000 or minLong == 1000000:
            logger.warning("Bounding box ill-defined, please check SAFE directory!")
    percent += step_percent
    logger.info("percent %d", percent)
    logger.info("end extract coordinates")

    # =====================
    # Fetch tiles
    DEM_SOURCE = args.s
    logger.info(f"start building tile list and fetching tiles using {DEM_SOURCE}")

    # Search for .netrc file (Earthdata USGS credientials) in standard locations
    netrcfilepath = None
    kwargs_fetcher = {}
    if 'HOME' in os.environ: # First, home directory
        netrcfilepathsearch = os.path.join(os.environ['HOME'] , '.netrc')
        if os.path.exists(netrcfilepathsearch):
            netrcfilepath = netrcfilepathsearch
        else:
            if 'PATH' in os.environ: # Second, in the PATH
                for mypath in os.environ['PATH'].split(os.pathsep):
                    netrcfilepathsearch = os.path.join(mypath , '.netrc')
                    if os.path.exists(netrcfilepathsearch):
                        netrcfilepath = netrcfilepathsearch
                        break # Just use the first file that is found

    if netrcfilepath is None:
        logger.info('Could not locate file containing Earthdata USGS credientials in $HOME or $PATH. Using default credentials...')
    else:
        logger.info("Found .netrc file (Earthdata USGS credientials) in %s." % (netrcfilepath))
        # Attempt to read login / password from that file
        with open(netrcfilepath) as fp:
            login = None
            pwd = None
            line = fp.readline()
            while len(line.split())==2: # Only read lines that include exactly two words
                if re.match(line.split()[0], 'login'):
                    login = line.split()[-1]
                if re.match(line.split()[0], 'password'):
                    pwd = line.split()[-1]
                line = fp.readline()
        if login is not None and pwd is not None:
            logger.info("Earthdata USGS credientials read with success from %s." % (netrcfilepath))
            kwargs_fetcher = {'login': login, 'pwd': pwd}
        else:
            logger.info("Earthdata USGS credientials could not be understood in %s. Using default credentials..." % (netrcfilepath))

    fetcher = demfetch.DemFetcher(keep=args.k, srtm_source=DEM_SOURCE, **kwargs_fetcher)

    if not os.path.exists(args.t):
        os.makedirs(args.t)
    logger.info('Fetching tiles minlat=%s minlong=%s maxlat=%s maxlong=%s',
                minLat, minLong, maxLat, maxLong)
    try:
        tile_list = fetcher.build_tile_list(minLat, minLong, maxLat, maxLong)
        logger.info("tile list={}".format(tile_list))
        #fetcher.download_unzip(tile_list, cache_dir)
    except Exception as ex:
        msg = "ERROR while fetching tiles ABORTING: %s" % str(ex)
        logger.critical(msg)
        if not args.k:
            clean_and_die(1, msg, args.t)
    if cache_dir is not None:
        if DEM_SOURCE == 'SRTMGL1.003':
            tile_names = (x.replace("SRTMGL1.hgt.zip", "hgt") for x in tile_list)
        else:
            tile_names = (x.replace("NASADEM_HGT_", "").replace("zip", "hgt") for x in tile_list)
        nsbio.make_dir_and_links(cache_dir, args.t, tile_names)

    percent += step_percent
    logger.info("end building tile list and fetching tiles")
    logger.info("percent %d", percent)

    logger.info("start merging tiles")
    # =====================
    # Merge tiles
    str_int_coordinates = [str(int(x)) for x in [minLong, maxLong, minLat, maxLat]]
    if DEM_SOURCE == 'SRTMGL1.003':
        file_prefix = os.path.join(args.t, "_".join(["srtm30"] + str_int_coordinates))
        extension = ".hgt"
    elif DEM_SOURCE == 'NASADEM_HGT.001':
        file_prefix = os.path.join(args.t, "_".join(["nasadem30"] + str_int_coordinates))
        extension = ".hgt"
    else:
        file_prefix = os.path.join(args.t, "_".join(["cop_dem30"] + str_int_coordinates))
        extension = ".dt2"
    merged_tiles_file_name = file_prefix + "_merged.tiff"
    try:
        merge_tiles(extension, args.t, merged_tiles_file_name)
    except Exception as e:
        logger.critical(f"cannot merge tiles: {e}")
        sys.exit(1)
    err = file_check(merged_tiles_file_name, 0)
    if err is not None and not args.k:
        logger.critical("file checking failed, aborting")
        clean_and_die(1, err[1], args.t)
    percent += step_percent
    logger.info("percent %d", percent)
    logger.info("end merging tiles")
    # =====================
    # Sample DEM
    # =====================

    logger.info('start resampling DEM')
    resample_file_name = file_prefix + "_resampled.tiff"
    try:
        sample_srtm(args.t, merged_tiles_file_name, resample_file_name)
    except Exception as e:
        logger.critical(f"cannot sample DEM: {e}")
        sys.exit(2)
    err = file_check(resample_file_name, 0)
    if err is not None and not args.k:
        logger.critical("file checking failed, aborting")
        clean_and_die(2, err[1], args.t)

    percent += step_percent
    logger.info("percent %d", percent)
    logger.info('end resampling DEM')

    # =====================
    # Get geocoordinates
    # =====================
    logger.info('start geting geocoordinates')

    geo_coordinates = get_geographic_coordinates(args.t, resample_file_name)
    logger.info('geting geocoordinates done')
    percent += step_percent
    logger.info('end geting geocoordinates')
    logger.info("percent %d", percent)

    # =====================
    # Crop geoid
    # =====================
    logger.info("start cropping geoid")
    if DEM_SOURCE == "COP_DEM":
        egm_model = 'EGM2008'
        egm_tif = "us_nga_egm2008_25.tif"
        # egm_model = 'EGM96'
        # egm_tif = "egm96-5-flipped.tiff"
    else:
        egm_model = 'EGM96'
        egm_tif = "egm96-5-flipped.tiff"
    egm_file = args.t + "/" + egm_tif
    if not os.path.exists(egm_file):
        # backoff to env variable
        flipped_dem = None
        if "NSBAS" in os.environ and os.path.exists(os.environ["NSBAS"] + '/eap/' + egm_model + '/' + egm_tif):
            flipped_dem = os.environ["NSBAS"] + '/eap/' + egm_model + '/' + egm_tif
        elif "EGM96FLIPPED" in os.environ and os.path.exists(os.environ["EGM96FLIPPED"]):
            flipped_dem = os.environ["EGM96FLIPPED"]
        else:
            logger.critical(("cannot find EGM geoid TIFF file, "
                             "neither in %s nor in as env variables", args.t))
            clean_and_die(3, err[1], args.t)
        os.symlink(flipped_dem, egm_file)
    if egm_tif == "us_nga_egm2008_25.tif":
        cropped_egm = os.path.join(args.t, "us_nga_egm2008_25-cropped.tif")
    else:
        cropped_egm = os.path.join(args.t, "egm96-5-flipped-cropped.tiff")
    try:
        crop_egm(geo_coordinates, args.t, egm_file, cropped_egm)
    except Exception as excpt:
        logger.critical("cannot crop EGM: %s", str(excpt))
        clean_and_die(3, err[1], args.t)
    err = file_check(cropped_egm, 0)
    if err is not None and not args.k:
        clean_and_die(3, err[1], args.t)
    percent += step_percent
    logger.info("percent %d", percent)
    logger.info("end cropping geoid")

    # =======================
    # cleaning dem dir
    # =======================
    if not args.k:
        logger.warning("cleaning working dir %s", args.t)
    clean_dem_dir(args.t)
    percent = 100.0
    logger.info("percent %d", percent)
