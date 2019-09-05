#!/usr/bin/env python2.7

from __future__ import print_function

import subprocess
import re
import os
import sys
import logging
import dem
import nsbas.nsbio as nsbio

import nsbas.utils as nsb_utils
import dem.fetch_raw as demfetch

# managing version
revision = "$Id: nsb_getDemFile.py 1993 2017-06-25 20:22:56Z rgrandin $"
__version__ = revision.split()[2]


# ################## Logging  ############################
#logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
#                    level=logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# ############## Definition of functions ###################


def clean_dem_dir(directory):
    """
    clean directory, removing temporal files
    """
    for fic in os.listdir(directory):
        abs_file = os.path.join(directory, fic)
        logger.debug("checking %s", abs_file)
        if re.search("(\.tiff|\.zip|\.hgt)$", abs_file):
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
        proc = subprocess.Popen([gdalmerge, "-o",
                                out_file,
                                "-of", "GTIFF"] +
                                [currdir + '/' + f for f in os.listdir(currdir)
                                 if f.endswith(extension)],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate()
        exitcode = proc.returncode

    except OSError:
        err_msg = "{} command not found".format(gdalmerge)
        logger.error(err_msg)
        raise OSError(err_msg)
    if exitcode is not 0:
        err_msg = "command failed \n EXIT CODE={}\nstdout={}\nstderr={}".format(exitcode, out, err)
        logger.error(err_msg)
        raise RuntimeError(err_msg)
    logger.debug("result of {} is: {}".format(gdalmerge, out_file))
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
    if exitcode is not 0:
        raise RuntimeError("command FAILED\nEXIT CODE={}\nstdout={}\nstderr={}"
                           .format(exitcode, out, err))
    logger.debug("result of {} is: {}".format(gdalwarp, out_file))


def get_geographic_coordinates(currdir, in_file):
    """
    Get the geographic coordinates of the created tiff file.
    :param in_file: resampled tiff srtm
    :type: str
    """
    coordinates = []
    gdalinfo = "gdalinfo"

    try:
        logger.info('Get the geographic coordinates of the resampled file.')
        proc = subprocess.Popen([gdalinfo,
                                 os.path.join(currdir, in_file)],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate()
        exitcode = proc.returncode
    except OSError:
        raise OSError("{} command not found".format(gdalinfo))
    if exitcode is not 0:
        raise RuntimeError("command failed\nExit code={}\nstdout={}\nstderr={}"
                           .format(exitcode, out, err))

    match_min = re.search(r"Lower Left\s*\(\s*(-?[0-9]*\.[0-9]*),\s*(-?[0-9]+\.?[0-9]*)", out, re.M)
    if match_min:
        coordinates += match_min.groups()
    else:
        raise RuntimeError("re: match fail for lower left corner")
    match_max = re.search(r"Upper Right\s*\(\s*(-?[0-9]*\.[0-9]*),\s*(-?[0-9]*\.[0-9]*)", out, re.M)
    if match_max:
        coordinates += match_max.groups()
    else:
        raise RuntimeError("re: match fail for upper right corner")
    logger.info(coordinates)
    logger.debug("result of {} is: {}".format(gdalinfo, coordinates))
    return coordinates


def get_latLong_coordinates(in_file):
    """
    Get the lattitude and longitude coordinates of the radar data.
    :param in_file: path to tiff radar data
    :type: str
    """
    coordinates = []
    gdalinfo = "gdalinfo"
    logger.debug("get_lat_long_coordinates: checking %s", in_file)
    try:
        proc = subprocess.Popen([gdalinfo, in_file],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate()
        exitcode = proc.returncode
    except OSError:
        raise OSError("{} command not found".format(gdalinfo))
    if exitcode is not 0:
        raise RuntimeError("command failed\nExit code={}\nstdout={}\nstderr={}"
                           .format(exitcode, out, err))
    match_min = re.search(r"\(0,0\)\s*->\s*\(\s*(-?[0-9]*\.[0-9]*), *(-?[0-9]+\.[0-9]*)?", out)
    if match_min:
        coordinates += match_min.groups()
    else:
        raise RuntimeError("re: match fail for long min")
    match_max = re.search(r"\(\s*(-?[0-9]*\.[0-9]*),\s*(-?[0-9]+\.[0-9]*),\s*-?[0-9]*\.?[0-9]* *\) *\nMetadata ?", out, re.M)
    if match_max:
        coordinates += match_max.groups()
    else:
        raise RuntimeError("re: match fail for long max")
    logger.info("get_latLong_coordinates: coordinates {}".format(coordinates))
    return [float(x) for x in coordinates]


def crop_egm96(geo_coordinates, currdir, egm_file, out_file):
    """
    Extract the geoid undulations for the study area from a global map.
    :param geo_coordinates: geocoordinates of the resampled tiff file
    :param egm_file: geoid file
    :param out_file: cropped geoid file
    :type geo_coordinates: array of str
    :type egm_file: str
    :type out_file:str
    """
    if os.path.exists(out_file):
        logger.info("crop egm96: file %s exist, remove it to be gdalwrap compliant", out_file)
        os.unlink(out_file)
    try:
        logger.info('Extract the geoid undulations for the study area from a global map.')
        proc = subprocess.Popen(["gdalwarp", "-t_srs",
                                 "+proj=longlat +datum=WGS84 +no_defs", "-te",
                                 geo_coordinates[0], geo_coordinates[1],
                                 geo_coordinates[2], geo_coordinates[3],
                                 "-tr", "0.000277778", "-0.000277778",
                                 "-r", "cubic", os.path.join(currdir,egm_file), out_file],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate()
        exitcode = proc.returncode
    except OSError:
        raise OSError("{} command not found".format("gdalwarpd"))
    if exitcode is not 0:
        raise RuntimeError("command failed\nExit code={}\nstdout={}\nstderr={}"
                           .format(exitcode, out, err))

    logger.debug("result of {} is: {}".format("cropGeoid", out_file))


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
                    raise OSError("{} command not found".format(gdalcalculate))
    if exitcode is not 0:
        raise RuntimeError("command failed\nExit code={}\nstdout={}\nstderr={}"
                           .format(exitcode, out, err))
    logger.debug("result of {} is: {}".format(gdalcalculate, out_file))


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
    if exitcode is not 0:
        raise RuntimeError("command failed\nExit code={}\nstdout={}\nstderr={}"
                           .format(exitcode, out, err))

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
        return 1, "file {} does not exist".format(file_name)
    if os.stat(file_name).st_size < file_min_size:
        return 2, "file {} exist but has a weird size ({})"\
            .format(file_name, os.stat(file_name).st_size)
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
    parser.add_argument("-v", type=int, default=3,
                        help=("set logging level: 0 critical, 1 error, "
                              "2 warning, 3 info, 4 debug, default=info"))
    parser.add_argument("-V", action="store_true",
                        help="print version (ie svn revision) and exit")
    parser.add_argument("-k", action="store_true",
                        help=("if set, keep tiles already downloaded, "
                              "do not re-download if present in directory"))
    parser.add_argument("-c", type=str, default=None,
                        help=("cache dir: look target files in cache dir."
                              "if not present in cache dir, download tiles in "
                              "cache dir. If option is not provided, download"
                              "locally"))
    parser.add_argument("-t", type=str, default=".",
                        help=("defines the target dir, ie the working dir. "
                              "Target dir is created if it does not exist"))

    coord_in = parser.add_mutually_exclusive_group(required=True)
    coord_in.add_argument("--bbox", type=str,
                          help=("bounding box of the data, in the form "
                                "longLowerLeft,latLowerLeft,longUpperRight,"
                                "latUpperRight"))
    coord_in.add_argument("--safe", type=str,
                          help=("directory that contains the images from"
                                "which we extract the bounding box"))

    # search for the -V flag in order to avoid pb with mandatory arguments
    if "-V" in sys.argv:
        print('Version (ie SVN revision): {}'.format(__version__))
        sys.exit(0)

    args = parser.parse_args()
    logging_translate = [logging.CRITICAL, logging.ERROR, logging.WARNING,
                         logging.INFO, logging.DEBUG]
    logger = logging.getLogger(__name__)
    logger.setLevel(logging_translate[args.v])
    # reload imported module with the correct log level
    logger = nsb_utils.change_log_level(args.v,
                                        [nsb_utils, nsbio, dem, demfetch])
    # =====================
    # normalize paths
    if not args.t.startswith('/'):
        args.t = os.path.abspath(args.t)
        logger.debug("target dir: %s", args.t)

    if args.c is not None:
        if not args.c.startswith('/'):
            args.c = os.path.abspath(args.c)
            logger.debug("cache dir: %s", args.c)
    percent = 0.0
    step_percent =100.0/9.0

    # ======================
    # extracting coordinates
    if args.bbox is not None:
        (minLong, minLat, maxLong, maxLat) = (float(x) for x in args.bbox.split(","))
    if args.safe is not None:
        if not args.safe.startswith('/'):
            args.safe = os.path.abspath(args.safe)
        logger.debug("safe dir: %s", args.safe)
        maxLong = -180
        maxLat = -90
        minLat = 1000000
        minLong = 1000000

        for root, dirs, files in os.walk(args.safe,followlinks=True):
            for name in files:
                if ".SAFE" in root \
                        and "measurement" in root \
                        and name.endswith(".tiff"):
                    abs_tiff = os.path.join(root, name)
                    logger.info("processing %s", abs_tiff)
                    latlong_coord = None
                    try:
                        latlong_coord = get_latLong_coordinates(abs_tiff)
                    except Exception as excpt:
                        msg = "Fail to extract lat/long coordinates from file:"
                        msg = msg + str(excpt)
                        logger.critical(msg)
                        clean_and_die(1, msg, args.t)
                    minLong = min(latlong_coord[0], latlong_coord[2], minLong)
                    minLat = min(latlong_coord[1], latlong_coord[3], minLat)
                    maxLong = max(latlong_coord[2], latlong_coord[0], maxLong)
                    maxLat = max(latlong_coord[3], latlong_coord[1], maxLat)
        logger.info("bounding box (minlong, maxlong, minlat, maxlat) = (%f, %f, %f,%f)",
                    minLong, maxLong, minLat, maxLat)
    percent += step_percent
    logger.info("percent %d",  percent)

    # =====================
    # Fetch tiles
    fetcher = demfetch.DemFetcher(url=("http://e4ftl01.cr.usgs.gov/"
                                       "MEASURES/SRTMGL1.003/"
                                       "2000.02.11/"),
                                  keep=True,
                                  log_level=logger.getEffectiveLevel())
    target_dir = args.c if args.c is not None else args.t
    if not os.path.exists(args.t):
        os.makedirs(args.t)
    logger.info('Fetching tiles minlat=%s minlong=%s maxlat=%s maxlong=%s',
                minLat, minLong, maxLat, maxLong)
    try:
        tile_list = fetcher.build_tile_list(minLat, minLong, maxLat, maxLong)
        logger.info("tile list={}".format(tile_list))
        fetcher.download_unzip(tile_list, target_dir)
    except Exception as ex:
        msg = "ERROR while fetching tiles ABORTING: %s", str(ex)
        logger.critical(msg)
        if not args.k:
            clean_and_die(1, msg, args.t)
    if args.c is not None:
        tile_names = (x.replace("SRTMGL1.hgt.zip", "hgt") for x in tile_list)
        nsbio.make_dir_and_links(args.c, args.t, tile_names)

    percent += step_percent
    logger.info("percent %d",  percent)

    # =====================
    # Merge tiles
    extension = ".hgt"
    str_int_coordinates = [str(int(x)) for x in [minLong, maxLong, minLat, maxLat]]
    file_prefix = os.path.join(args.t, "_".join(["srtm30"] + str_int_coordinates))
    merged_tiles_file_name = file_prefix + "_merged.tiff"
    try:
        merge_tiles(extension, args.t,  merged_tiles_file_name)
    except Exception as e:
        logger.critical("cannot merge tiles: {}".format(str(e)))
        sys.exit(1)
    err = file_check(merged_tiles_file_name, 0)
    if err is not None and not args.k:
        logger.critical("file checking failed, aborting")
        clean_and_die(1, err[1], args.t)
    percent += step_percent
    logger.info("percent %d",  percent)

    # =====================
    # Sample srtm
    # =====================

    logger.info('resampling srtm')
    resample_file_name = file_prefix + "_resampled.tiff"
    try:
        sample_srtm(args.t, merged_tiles_file_name, resample_file_name)
    except Exception as e:
        logger.critical("cannot sample srtm: {}".format(str(e)))
        sys.exit(2)
    err = file_check(resample_file_name, 0)
    if err is not None and not args.k:
        logger.critical("file checking failed, aborting")
        clean_and_die(2, err[1], args.t)

    percent += step_percent
    logger.info("percent %d",  percent)

    # =====================
    # Get geocoordinates
    # =====================
    logger.info('geting geocoordinates')

    geo_coordinates = get_geographic_coordinates(args.t, resample_file_name)
    logger.info('geting geocoordinates done')
    percent += step_percent
    logger.info("percent %d",  percent)


    # =====================
    # Crop geoid
    # =====================

    #egm_file = args.t + "/egm96-5-flipped.tiff"
    egm_file = "/home/cometraid1/geoide/egm96-5-flipped.tiff"
    if not os.path.exists(egm_file):
        # backoff to env variable
        flipped_dem = None
        if "NSBAS" in os.environ and os.path.exists(os.environ["NSBAS"] + '/eap/EGM96/egm96-5-flipped.tiff'):
            flipped_dem = os.environ["NSBAS"] + '/eap/EGM96/egm96-5-flipped.tiff'
        elif "EGM96FLIPPED" in os.environ and os.path.exists(os.environ["EGM96FLIPPED"]):
            flipped_dem = os.environ["EGM96FLIPPED"]
        else:
            logger.critical("cannot find egm96-5-flipped.tiff, neither in %s nor in as env variables",
                    args.t)
            clean_and_die(3, err[1], args.t)
        os.symlink(flipped_dem, egm_file)
    cropped_egm = os.path.join(args.t, "egm96-5-flipped-cropped.tiff")
    try:
        crop_egm96(geo_coordinates, args.t, egm_file, cropped_egm)
    except Exception as e:
        logger.critical("cannot crop egm96: {}".format(str(e)))
        clean_and_die(3, err[1], args.t)
    err = file_check(cropped_egm, 0)
    if err is not None and not args.k:
        clean_and_die(3, err[1], args.t)
    percent += step_percent
    logger.info("percent %d",  percent)


    # =====================
    # Add geoid
    # =====================

    cropped_file_name = file_prefix + "_cropped.tiff"
    try:
        add_egm96_to_dem(cropped_egm, args.t, resample_file_name,
                         cropped_file_name)
    except Exception as e:
        logger.critical("cannot add egm96: {}".format(str(e)))
        clean_and_die(4, err[1], args.t)
    err = file_check(cropped_file_name, 0)
    if err is not None and not args.k:
        clean_and_die(4, err[1], args.t)
    percent += step_percent
    logger.info("percent %d",  percent)


    # =====================
    # Convert to dem
    # =====================

    dem_file_name = file_prefix+".dem"
    try:
        convert_to_dem(args.t, cropped_file_name, dem_file_name)
    except Exception as e:
        logger.critical("cannot add egm86: {}".format(str(e)))
        sys.exit(5)
    rscDemFile = open((dem_file_name+'.rsc'), "r")
    found = 0
    for l in rscDemFile:
        if l.startswith("WIDTH"):
            width = int(l.split()[1])
            found += 1
            if found == 3:
                break
        if l.startswith("FILE_LENGTH"):
            length = int(l.split()[1])
            found += 2
            if found == 3:
                break
    err = file_check(dem_file_name, width*length*2)
    if err is not None and not args.k:
        clean_and_die(5, err[1], args.t)
    percent += step_percent
    logger.info("percent %d",  percent)


    # =======================
    # cleaning dem dir
    # =======================
    if not args.k:
        logger.warning("cleaning working dir %s", args.t)
        clean_dem_dir(args.t)
    percent = 100.0
    logger.info("percent %d",  percent)



