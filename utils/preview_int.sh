#!/bin/sh

# Initial checkings
if [ "$1" = ""] || ["$2" = ""] || ["$5" = ""] || [ "$6" = ""] || ["$7" = ""] \
  || ["$8" = ""] ; then
	echo "usage : geotiff.sh home_dir liste_interfero.rsc suffix prefix rlook dlook radar dir"
	echo "home_dir: full path to T??? dir "
	echo "liste_interfero.rsc: path from home dir to rsc file"
  echo "suffix prefix rlook: ${prefix}$date1-$date2_${suffix}_${rlook}rlks"
  echo "dlook: downsample look of the output jpeg"
	echo "radar: full path from home dir to radar.hgt file"
	echo "dir: path to output dir form home"
  echo ; exit         
fi 

home=$1
list_int=$home/$2
suffix=$3
prefix=$4
rlook=$5
dlook=$6
radar=$home/$7
outdir=$home/$8
outlook=$(echo "scale=3; $rlook + $dlook" | bc )

cat < $list_int | while true
do
   read ligne
   if [ "$ligne" = "" ]; then break; fi
   set -- $ligne ; image1=$1 ; image2=$2 ;
   cd int_$1_$2/
   echo $1 $2
   
   # length.pl ${prefix}$1-$2${suffix}_${rlook}rlks.int
   nsb_preview_int -croipac -l$dlook $radar \
   ${prefix}$1-$2${suffix}_${rlook}rlks.int ${prefix}$1-$2${suffix}_${outlook}rlks.jpeg
   # nsb_preview_int -croipac -l$dlook $radar $1-$2_${rlook}rlks.int $1-$2_${outlook}rlks.jpeg

   cd ..
done

mkdir $outdir
mv int_*/*jpeg $outdir/.
