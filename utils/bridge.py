#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################################
# Author        : Simon DAOUT (Oxford)
############################################

"""
bridge.py
========================
Create bridge.in file given output form mdx wirtten in a tex file

Usage:
  bridge.py  <file>  

Options:
  file         input text file from mdx
"""

import docopt
import numpy as np
from collections import namedtuple
import re

arguments = docopt.docopt(__doc__)
infile = arguments["<file>"]


LINE_EXPR = r"(?P<option>.*?)\s*(?P<vi>:)\s*(?P<value>.*)$"
_re = re.compile(LINE_EXPR.format(delim=":"), re.VERBOSE)

lines = []
with open(infile) as in_file:
	for line in in_file:
		lines.append(line.strip())
print(lines)

