#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################################
# Author        : Simon DAOUT (Oxford)
############################################

"""
bridge.py
========================
Create bridge.in file given output form mdx wirtten in a text file
To create the text file: right clic on the two points to bridge on mdx windows and copy past terminal in a text file

Usage:
  bridge.py  <file>  

Options:
  file         input text file from mdx
"""

import docopt
import re
import sys

arguments = docopt.docopt(__doc__)
infile = arguments["<file>"]


# there is probably a better way to do that re search... 
#_re = re.compile(r"(.*) COL:(.*?) ROW:(.*) ", re.IGNORECASE)
_re = re.compile(r"(.*) COL: (.*) ROW: {1,6}([0-9]*) ([0-9]*)", re.IGNORECASE)

file = open('bridge.in', "w")
lines = []
count = 0
with open(infile) as in_file:
    for i, line in enumerate(in_file):
        if line.strip().startswith("#"):
            continue
        line = line.strip()
        if (not line): 
            continue

        m = _re.match(line)
        if m:
           # print(m.group(0))
           # print(m.group(2),m.group(3))
            col, row =  m.group(2).strip(), m.group(3).strip()
            print(count, col, row)
            if count % 2:
                file.write("%s %s 0\n" % (str(col), str(row)))
            else:
                file.write("%s %s " % (str(col), str(row)))
            count  += 1

