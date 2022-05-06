#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
destripe.py
===============

Usage:
   destripe.py [-d=<direction>] [-a=<angle>] [-m=<mask_file>] <img1_file> <im2_file> <output_file>

   Usage:
   destripe.py [-d=<direction>] [-a=<angle>] [-m=<mask_file>] <input_file> <output_file>

   Options:
   -h --help       Show this screen.
   -a=<angle>      Angle: can be orbit, manual or numeric value [default: orbit].
   -d=<direction>  Direction of the destriping: 'v'(ertical) or 'h'(orizontal)
   [default: v].
   -m=<mask_file>  Mask file (raster or .mat file).

"""
