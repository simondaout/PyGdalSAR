#!/usr/bin/env python2
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Authors        : Mathieu Volat
#                  Simon DAOUT (Oxford)
############################################

"""\
fetch_erai.py
-----------------
Retrieve ERAI data relative to a area of interest.

Usage:
  fetch_erai.py [-v] <proc_file>
  fetch_erai.py -h | --help

Options:
  -v         Verbose
  -h --help  Show this screen.

"""

from __future__ import print_function
import datetime, os, itertools, string

import docopt
import ecmwfapi

from nsbas.printfunctions import *
from nsbas import procparser, rscparser, utils

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)

# Parse arguments
arguments = docopt.docopt(__doc__)

# Read proc file
proc_defaults = { "Rlooks_sim": 1,
                  "Rlooks_int": 1,
                  "pres_levels": "1000/975/950/925/900/875/850/825/800/775/750/700/650/600/550/500/450/400/350/300/250/225/200/175/150/125/100/70/50/30/20/10/7/5/3/2/1",
                  # temperature 2m
                  "param_ecmwf": "167.128",
                  # precipitation
                  #"param_ecmwf": "228.128",
                  "grid_step_ecmwf": "0.75/0.75",
                  "time_steps_ecmwf": "0/6/12/18",
                  "grid_limits_ecmwf": "NULL",
                  "dates_list_ecmwf": "NULL",
                  "ncfile_ecmwf": "NULL" }
proc = procparser.ProcParser(proc_defaults)
proc.read(arguments["<proc_file>"])

# Create EraDir if needed
if not os.path.exists(proc.get("EraDir")):
    if arguments["-v"]:
        print("create "+proc.get("EraDir"))
    os.makedirs(proc.get("EraDir"))

# Find closest time available
if arguments["-v"]:
    print("time_of_day_ecmwf = "+proc.get("time_of_day_ecmwf")+", day_shift = "+str(day_shift))

if arguments["-v"]:
    print("dates_list_ecmwf = "+proc.get("dates_list_ecmwf"))

# Fetch the data
if arguments["-v"]:
    print("fetching data...")

server = ecmwfapi.ECMWFDataServer()
# First, data.nc
server.retrieve({
    "dataset"  : "interim",
    "stream"   : "oper",
    "step"     : "0",
    "levtype"  : "sfc",
    "date"     : proc.get("dates_list_ecmwf"),
    "time"     : proc.get("time_steps_ecmwf"),
    "type"     : "an",
    "param"    : proc.get("param_ecmwf"),
    "area"     : proc.get("grid_limits_ecmwf"),
    "grid"     : proc.get("grid_step_ecmwf"),
    "target"   : os.path.join(proc.get("EraDir"), "data.nc"),
    "format"   : "netcdf"})

if arguments["-v"]:
    print("work done")
exit(0)
