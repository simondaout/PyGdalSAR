#!/bin/bash

#home=/data/ATF/WATF/Envisat/T391
#rsc=$home/interf_pair_200_4_700_1_r_ts3.rsc
#jbeg=2000


# Initial checkings
if [ "$1" = ""] || ["$2" = ""] || ["$3" = ""] ; then
	echo "usage : cut_unw.sh path_to_home liste_interfero.rsc jbeg "
	echo "absolute path to dates"
	echo "liste_interfero.rsc: rsc file"
	echo "jbeg: define the cutting zone"
	echo ; exit         
fi 

home=$1
rsc=$home/$2
jbeg=$3
jend=$4
rep_post=$home/SBAS

cat < $rsc | while true
do
   read ligne
   echo $ligne
   if [ "$ligne" = "" ]; then break; fi
   set -- $ligne ; image1=$1 ; image2=$2 ;

cd  $rep_post/int_$1_$2

#rm -f "filtSW_cor_$1-$2_col_era_flatr_cut_8rlks.unw"
#cut_unw.py --infile="filtSW_cor_$1-$2_col_era_flatr_8rlks.unw" --outfile="filtSW_cor_$1-$2_col_era_flatr_cut_8rlks.unw" --plot=no $jbeg $jend
#cp filtSW_cor_$1-$2_col_era_flatr_8rlks.unw.rsc filtSW_cor_$1-$2_col_era_flatr_cut_8rlks.unw.rsc

rm -f "filtSW_cor_$1-$2_col_era1pt_flatr_cut_8rlks.unw"
cut_unw.py --infile="filtSW_cor_$1-$2_col_era1pt_flatr_8rlks.unw" --outfile="filtSW_cor_$1-$2_col_era1pt_flatr_cut_8rlks.unw" --plot=no $jbeg $jend
cp filtSW_cor_$1-$2_col_era1pt_flatr_8rlks.unw.rsc filtSW_cor_$1-$2_col_era1pt_flatr_cut_8rlks.unw.rsc

done
