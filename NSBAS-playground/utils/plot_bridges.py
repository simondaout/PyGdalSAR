#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

import os, getopt, posixpath as path, sys, warnings

import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
import gdal

mpl.rc("axes", titlesize="medium")

bridges_path = "bridge.in"
bridges = []

def update_ax(ds):
    global ax, already_updating

    b1 = ds.GetRasterBand(1)
    b2 = ds.GetRasterBand(2) if ds.RasterCount > 1 else None

    if b1.DataType == gdal.GDT_CFloat32:
        data = b1.ReadAsArray()
        amp, pha = np.absolute(data), np.angle(data)
        pha_min, pha_max = -np.pi, np.pi
    else:
        amp, pha = b1.ReadAsArray(), b2.ReadAsArray()
        pha_min, pha_max, pha_mean, pha_stddev = b2.GetStatistics(0, 1)
        pha_min, pha_max = pha_mean - pha_stddev*2.0, pha_mean + pha_stddev*2.0

    pha[pha == 0] = np.nan

    ax.imshow(amp, cmap="gray", interpolation="bicubic")
    im = ax.imshow(pha, cm.gist_rainbow, interpolation="none", vmin=pha_min, vmax=pha_max, alpha=0.5)

opts, argv = getopt.getopt(sys.argv[1:], "f:h")
for o, a in opts:
    if o == "-f":
        bridges_path = a
    elif o == "-h":
        print("%s [-f bridge_file] dataset1 [dataset2]" % (sys.argv[0]))
        exit(0)
    else:
        pass

# Open the interferogram
ds1 = gdal.Open(argv[0], gdal.GA_ReadOnly)
if ds1 == None:
    exit(1)
if len(argv) > 1:
    ds2 = gdal.Open(argv[1], gdal.GA_ReadOnly)
    if ds2 == None:
        exit(1)
else:
    ds2 = None

# Initialize the figure
fig, ax = plt.subplots(1, 1)
plt.subplots_adjust(left=0.3)

# Initialize ax image
update_ax(ds1)
# Initialize ax bridges
if path.exists(bridges_path):
    for line in open(bridges_path, "r"):
        cols = [int(x) for x in line.strip().split()]
        if cols[4] == 0:
            curstyle = dict(arrowstyle="-")
        elif cols[4] == 1:
            curstyle = dict(arrowstyle="->")
        elif cols[4] == -1:
            curstyle = dict(arrowstyle="<-")
        newbridge = ax.annotate("",
                                xy=(cols[2], cols[3]), xycoords="data", 
                                xytext=(0, 0), textcoords="data",
                                arrowprops=curstyle)
        newbridge.xyann = newbridge.xytext = (cols[0], cols[1])
        newbridge.bridgetype = int(cols[4])
        bridges.append(newbridge)

# Box to select the image
already_updating = False
def on_dsswitch_clicked(label):
    global already_updating
    if already_updating:
        return
    already_updating = True
    update_ax(ds1 if label == "input 1" else ds2)
    fig.canvas.draw()
    already_updating = False
if ds2 != None:
    iax = plt.axes([0.05, 0.7, 0.15, 0.15])
    iax.set_title("Dataset")
    dsswitch = widgets.RadioButtons(iax, ("input 1", "input 2"))
    dsswitch.on_clicked(on_dsswitch_clicked)

# Box to select the bridge type
rbridgetype_value_selected = "0"
def on_rbridgetype_clicked(label):
    global rbridgetype_value_selected
    rbridgetype_value_selected = label
rax = plt.axes([0.05, 0.5, 0.15, 0.15])
rax.set_title("Bridge type")
rbridgetype = widgets.RadioButtons(rax, ("0", "1", "-1"))
rbridgetype.on_clicked(on_rbridgetype_clicked)

# Box to save bridges
def dosave(event):
    bridgesf = open(bridges_path, "w")
    for b in bridges:
        bridgesf.write("%d %d %d %d %d\n" % (b.xyann[0], b.xyann[1], b.xy[0], b.xy[1], b.bridgetype))
    bridgesf.close()
bax = plt.axes([0.05, 0.1, 0.15, 0.05])
bsave = widgets.Button(bax, "Save to file")
bsave.on_clicked(dosave)

# Global event handler functions
shift_is_held = False
def on_key_press(event):
    global shift_is_held
    if event.key == "shift":
        shift_is_held = True
def on_key_release(event):
    global shift_is_held
    if event.key == "shift":
        shift_is_held = False
def on_button_press(event):
    if not shift_is_held or event.inaxes != ax:
        return
    if rbridgetype_value_selected == "0":
        curstyle = dict(arrowstyle="-")
    elif rbridgetype_value_selected == "1":
        curstyle = dict(arrowstyle="->")
    elif rbridgetype_value_selected == "-1":
        curstyle = dict(arrowstyle="<-")
    newbridge = ax.annotate("",
                            xy=(0, 0), xycoords="data", 
                            xytext=(0, 0), textcoords="data",
                            arrowprops=curstyle)
    if hasattr(newbridge, "xyann"): # matplotlib >= 1.4
        newbridge.xyann = (event.xdata, event.ydata)
    else:
        newbridge.xytext = (event.xdata, event.ydata)
    newbridge.bridgetype = int(rbridgetype_value_selected)
    bridges.append(newbridge)
def on_button_release(event):
    if not shift_is_held or event.inaxes != ax:
        return
    bridges[-1].xy = (event.xdata, event.ydata)
    fig.canvas.draw()
# Register the event handler functions
fig.canvas.mpl_connect("key_press_event", on_key_press)
fig.canvas.mpl_connect("key_release_event", on_key_release)
fig.canvas.mpl_connect("button_press_event", on_button_press)
fig.canvas.mpl_connect("button_release_event", on_button_release)

def pick_handler(event):
    mouseevent = event.mouseevent
    artist = event.artist
    print(artist)

plt.show()
