#!/bin/csh

# Check input arguments
if ($#argv < 2 || $#argv > 3) then
    echo " "
    echo "Usage: grd2kml.csh grd_file_stem cptfile [-R<west>/<east>/<south>/<north>] "
    echo " "
    echo "Example: grd2kml.csh phase phase.cpt "
    echo " "
    exit 1
endif

# Extract grid resolution (cell size)
set DX = `gmt grdinfo $1.grd -C | awk '{print $8}'`
if ($status != 0) then
    echo "Error: Unable to extract grid information from $1.grd"
    exit 1
endif

# Define output resolution
set DPI = 2000
echo "Using DPI: $DPI"

# forcer GMT à recalculer les limites et résolutions de la grille.
gmt grdedit $1.grd -A

# Set GMT parameters for visualization
gmt gmtset COLOR_MODEL=hsv
gmt gmtset PS_MEDIA=letter

# Uncomment to calculate gradients (currently not used)
# gmt grdgradient $1.grd -Ggrad.grd -Nt0.7 -A60 -V

# Generate the image with optional region specification
if ($#argv == 3) then
    gmt grdimage $1.grd -C$2 $3 -Jx1id -P -Y2 -X2 -Q -V > $1.ps
else
    gmt grdimage $1.grd -C$2 -Jx1id -P -Y2 -X2 -Q -V > $1.ps
endif

# Convert the PostScript file to PNG and embed into KML
#gmt psconvert $1.ps -W+k+t"$1" -TG -P -S -V -F$1.png
#-A+m0.5c → Adds margins to prevent clipping (adjustable if needed).
#-E$DPI → Sets the resolution to the defined DPI.
#-S → Strips unnecessary metadata for a cleaner output.
#-F$1.png → Specifies the output file name.
gmt psconvert $1.ps -W+k+t"$1" -TG -P -S -A+m0.5c -E$DPI -F$1.png

# Clean up temporary files
rm -f $1.ps grad.grd  gmt.history gmt.conf
rm psconvert_*

echo "Conversion complete: $1.kml and $1.png generated."

