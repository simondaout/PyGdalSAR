#!/bin/sh

# Initial checkings
if [ "$1" = ""] || ["$2" = ""] || ["$3" = ""] || [ "$4" = ""] || ["$5" = ""] || ["$6" = ""] ; then
	echo "usage : geotiff.sh home_dir liste_interfero.rsc suffix prefix rlook geomaptrans "
	echo "home_dir: full path to T??? dir "
	echo "liste_interfero.rsc: path from home dir to rsc file"
    echo "suffix prefix: ${prefix}$date1-$date2_${suffix}"
	echo "geomaptrans: full path from home dir to geomap trans file"
	echo ; exit         
fi 

home=$1
list_int=$home/$2
suffix=$3
prefix=$4
rlook=$5
geo=$home/$6

cat < $list_int | while true
do
   read ligne
   if [ "$ligne" = "" ]; then break; fi
   set -- $ligne ; image1=$1 ; image2=$2 ;
   cd int_$1_$2/
   
   cp filt_$1-$2_cor_flatrange_flatcutbox_4rlks.unw.rsc ${prefix}$1-$2_${suffix}_${rlook}rlks.unw.rsc
   rm -f geo_$1-$2_${rlook}rlks.unw*
   geocode.pl $geo ${prefix}$1-$2_${suffix}_${rlook}rlks.unw geo_$1-$2_${rlook}rlks.unw
   gdal_translate -ot Float32 -b 2 -co COMPRESS=DEFLATE -co COMPRESS=PREDICTOR2 geo_$1-$2_${rlook}rlks.unw geo_$1-$2_${rlook}rlks.tiff 
   # gdal_translate -of JPEG -b 2 geo_$1-$2_${rlook}rlks.unw geo_$1-$2_${rlook}rlks.jpg

   cd ..
done
