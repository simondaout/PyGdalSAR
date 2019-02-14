#!/bin/sh
# daout simon

# Initial checkings
if [ "$1" = ""] || ["$2" = ""] ||  [ "$4" = ""] || ["$5" = ""] \
  || ["$6" = ""] || ["$7" = ""] || ["$8" = ""] ; then
	echo "usage : geotiff.sh home_dir liste_interfero.rsc suffix prefix rlook radar wrap dir"
	echo "home_dir: full path to T??? dir "
	echo "liste_interfero.rsc: path from home dir to rsc file"
  echo "suffix prefix: ${prefix}$date1-$date2_${suffix}"
	echo "radar: full path from home dir to radar.hgt file"
  echo "wrap: wrapped cycle of the phase"
	echo "dir: path to output dir form home"
  echo ; exit         
fi 

home=$1
list_int=$home/$2
suffix=$3
prefix=$4
rlook=$5
radar=$home/$6
min_cb=-$7
max_cb=$7
outdir=$home/$8

clean
rm -f $outdir/${prefix}????????-????????${suffix}_${rlook}rlks.jpeg

cat < $list_int | while true
do
   read ligne
   if [ "$ligne" = "" ]; then break; fi
   set -- $ligne ; image1=$1 ; image2=$2 ;
   cd ${home}/int/int_$1_$2
   echo $1 $2
	
   #length.pl ${prefix}$1-$2${suffix}_${rlook}rlks.unw    
   nsb_preview_unw -croipac -m${min_cb} -M${max_cb} $radar \
   ${prefix}$1-$2${suffix}_${rlook}rlks.unw  ${prefix}$1-$2${suffix}_${rlook}rlks.jpeg

   cd ..
done


mkdir $outdir
mv ${home}/int/int_*/${prefix}????????-????????${suffix}_${rlook}rlks.jpeg $outdir/.
