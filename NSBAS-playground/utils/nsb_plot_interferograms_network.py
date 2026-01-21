#!/usr/bin/env python3
##########################################################################
#
#   This file is part of NSBAS.
#
#   NSBAS is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   NSBAS is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with NSBAS.  If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################

"""\
nsb_plot_interferograms_network

Usage:
  nsb_plot_interferograms_network [-o <out_path>] <pair_file> <baseline_file> [-w <y/n>]
  nsb_plot_interferograms_network -h | --help

Options:
  -o OUT_PATH  Write figure to file.
  -h --help    Show this screen.
  -w <y/n>     show the weight used for TS analysis in the third column of pair_file (default: no)

"""

import string, os
from datetime import datetime

import numpy as np
import matplotlib as mpl
#if "DISPLAY" not in os.environ:
#    mpl.use('Agg')
#else:
#    mpl.use('pdf')

from nsbas import docopt
arguments = docopt.docopt(__doc__)
if arguments["-o"]:
   mpl.use('pdf')

if arguments["-w"] == None:
    w = 'n'
else :
    w = 'y'

import matplotlib.pyplot as plt
from matplotlib import dates
from matplotlib.dates import date2num, num2date
import matplotlib.cm as cm
import matplotlib.colors as colors

import pydot
#except ModuleNotFoundError as e:
#from nsbas import pydot

fig = plt.figure(figsize=(13,6))
ax = fig.add_subplot(111)

# Load graph from <pair_file> and <baseline_file>
graph = pydot.Dot("interferogram_network", graph_type="digraph")
weight = []
for line in open(arguments["<baseline_file>"], "r"):
    lines = line.split()
    date, pbaseline, rest = lines[0], lines[1], lines[2]
    graph.add_node(pydot.Node(date.strip(), label=date.strip(), bperp=float(pbaseline)))
for line in open(arguments["<pair_file>"], "r"):
    if w == 'n':
        date1, date2 = line.split()[0], line.split()[1]
        graph.add_edge(pydot.Edge(date1.strip(), date2.strip()))
        weight.append(1)
    else :
        date1, date2, wei = line.split()[0], line.split()[1], float(line.split()[2])
        graph.add_edge(pydot.Edge(date1.strip(), date2.strip()))
        weight.append(wei)

if w == 'y':
    norm = colors.Normalize(vmin=min(weight), vmax=max(weight))
    cmap = cm.PuBu

# Draw nodes
x, y = [], []
for n in graph.get_nodes():
    x.append(date2num(datetime.strptime(n.get_label(), "%Y%m%d")))
    y.append(float(n.get_attributes()["bperp"]))
ax.plot(x, y, "o", color='dodgerblue', mec='black', markersize=4, picker=5)
# Draw arrows
for i, edge in enumerate(graph.get_edges()):
    master = graph.get_node(edge.get_source())[0]
    slave = graph.get_node(edge.get_destination())[0]
    x = date2num(datetime.strptime(master.get_label(), "%Y%m%d"))
    y = float(master.get_attributes()["bperp"])
    dx = date2num(datetime.strptime(slave.get_label(), "%Y%m%d")) - x
    dy = float(slave.get_attributes()["bperp"]) - y

    if w == 'y':
        wei = weight[i]
        color = cmap(norm(wei))
    else:
        color = 'black'
    
    ax.arrow(x, y, dx, dy, linewidth=.5, color=color, alpha=.5)

if w == 'y':
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Weight')

# Register click usage to display date of nearest point
def onpick(event):
    global ax, an
    pt = event.artist
    ind = event.ind
    x = np.take(pt.get_xdata(), ind)[0]
    y = np.take(pt.get_ydata(), ind)[0]
    an.xy = (x, y)
    if hasattr(an, "xyann"): # matplotlib >= 1.4
        an.xyann = (x + abs(ax.get_xlim()[1]-ax.get_xlim()[0])/30,
                    y + abs(ax.get_ylim()[1]-ax.get_ylim()[0])/15)
    else:
        an.xytext = (x + abs(ax.get_xlim()[1]-ax.get_xlim()[0])/30,
                     y + abs(ax.get_ylim()[1]-ax.get_ylim()[0])/15)
    an.set_text(num2date(x).strftime("%Y/%m/%d"))
    an.set_visible(True)
    fig.canvas.draw()
#an = ax.annotate("",
#                 xy=(0, 0), xycoords="data",
#                 xytext=(0, 0), textcoords="data",
#                 arrowprops=dict(facecolor="black", width=1, frac=0.3),
#                 bbox=dict(boxstyle="round", fc="w"))
#an.set_visible(False)
fig.canvas.mpl_connect("pick_event", onpick)
# Show
#fig.canvas.set_window_title("Interferogram network")
fig.suptitle("Interferogram network")
fig.autofmt_xdate()
ax.set_xlabel("Acquisition Date")
ax.set_ylabel("Perpendicular Baseline (m)")
ax.xaxis.set_major_formatter(dates.DateFormatter("%Y/%m/%d"))
#plt.margins(0.1, 0.1)
if arguments["-o"]:
    plt.savefig(arguments["-o"], bbox_inches="tight")
else:
    plt.show()
