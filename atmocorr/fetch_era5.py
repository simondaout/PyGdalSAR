#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Simon DAOUT (Oxford)
############################################

import ecmwfapi

server = ecmwfapi.ECMWFDataServer()
server.retrieve({
    "class": "ea",
    "dataset": "era5",
    "date": "2007-01-01/to/2010-12-31", # Time period
    "expver": "1",
    "levtype": "pl",
    "levelist": "850/700/600/500",
    "param": "130",           # Parameters. Here we use 2m Temperature (2t) and Surface Pressure (sp). See the ECMWF parameter database, at http://apps.ecmwf.int/codes/grib/param-db
    "stream": "oper",
    "type": "an",
    "time": "0/12",
    "step": "0",
    "area": "27.25/89.00/28.25/90.00",    # Subset or clip to an area, here to Europe. Specify as North/West/South/East in Geographic lat/long degrees. Southern latitudes and Western longitudes must be given as negative numbers.
    "grid": "0.25/0.25",        # Regrid from the default grid to a regular lat/lon with specified resolution. The first number is east-west resolution (longitude) and the second is north-south (latitude).
    "format": "netcdf",         # Convert the output file from the default GRIB format to NetCDF format. Requires "grid" to be set to a regular lat/lon grid.
    "target": "data.nc",    # The output file name. Set this to whatever you like.
})
