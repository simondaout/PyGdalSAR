#!/bin/bash

# Initial checkings
if [ "$1" = ""] || ["$2" = ""] || [ "$3" = ""] || ["$4" = ""] || ["$5" = ""] || ["$8" = ""] ; then
   echo "usage : make_unwcol.sh home_dir liste_interfero.rsc master olook ilook suffix prefix width seedx seedy threshold"
   echo "home_dir: full path to T??? dir "
   echo "liste_interfero.rsc: path fom home dir to rsc file"
   echo "master: master date"
   echo "olook: look unwrapped interferograms"
   echo "ilook: look input wrapped interferograms"
   echo "int name: ${prefix}$1-$2${suffix}_${ilook}rlks.int"
   echo "width: number col interferogram"
   echo "seedx, seedy: col and line starting point for unwrapping"
   echo "threshold: threshold cohrence filter int. for unwrapping"
   echo ; exit         
fi 

# # input parameters 
# home='/home/cometraid14/daouts/work/tibet/qinghai/processing/Sentinel/iw3'
# master='20160608'
# suffix='_sd_flatrange_flatz'
# prefix='coh_'
# olook=4
# ilook=2
# width=1526
# list_int=$1
# seedx=268
# seedy=1766
# thresholdfiltSW=0.45 # firt round

home=$1
list_int=$2
master=$3
olook=$4
ilook=$5
suffix=$6
prefix=$7
width=$9
seedx=$10
seedy=$11
thresholdfiltSW=$12

################################################################################
# harcoding parameters

filterstyle=SWc
SWamplim=0.1
SWwindowsize=8
thresholdcol=0.04

################################################################################

look=$(echo "scale=3; $olook - $ilook" | bc )
dem="$home/$master/radar_${olook}rlks.hgt"

if [ ! -e $home/$master/radar_${olook}rlks.hgt ]; then
   look.pl  $home/$master/radar_${ilook}rlks.hgt ${look}
fi

cat < $list_int | while true
do
   read ligne
   if [ "$ligne" = "" ]; then break; fi
   set -- $ligne ; image1=$1 ; image2=$2 ;
   cd $home/int/int_$1_$2
   echo $1 $2

   if [ ! -e ${prefix}${image1}-${image2}${suffix}_${olook}rlks.int ]; then
   look.pl  ${prefix}${image1}-${image2}${suffix}_${ilook}rlks.int ${look}
   fi
   
   if [ ! -e $${image1}-${image2}_${olook}rlks.cor ]; then
   look.pl ${image1}-${image2}_${ilook}rlks.cor ${look}
   fi

   cor=${image1}-${image2}_${olook}rlks

   if [ -e bridge.asc ]; then
      bridge.sh bridge.asc
   fi

   if [ ! -e bridge.in ]; then
      echo "1 1 1 1 0" > bridge.in
   fi

   echo '......................'
   echo 'Filtering'
   echo '......................'

   # clean filters and colin
   rm -f col_${image1}-${image2}${suffix}_${olook}rlks.int
   rm -f filtSW_${prefix}${image1}-${image2}${suffix}_${olook}rlks.int
   rm -f filt_col_${image1}-${image2}${suffix}_${olook}rlks.int
   length.pl ${prefix}${image1}-${image2}${suffix}_${olook}rlks.int

   # compute colinearity
   cp ${prefix}${image1}-${image2}${suffix}_${olook}rlks.int col_${image1}-${image2}${suffix}_${olook}rlks.int
   cp ${prefix}${image1}-${image2}${suffix}_${olook}rlks.int.rsc \
   col_${image1}-${image2}${suffix}_${olook}rlks.int.rsc
   cp ${prefix}${image1}-${image2}${suffix}_${olook}rlks.int temp
   length=`/bin/grep FILE_LENGTH ${prefix}${image1}-${image2}${suffix}_${olook}rlks.int.rsc  | awk '{print $2}'`
   
   colin ${prefix}${image1}-${image2}${suffix}_${olook}rlks.int temp \
   col_${image1}-${image2}${suffix}_${olook}rlks.int $width $length 3 0.0001 2

   # adapt filt with col
   myadapt_filt col_${image1}-${image2}${suffix}_${olook}rlks.int \
    filt_col_${image1}-${image2}${suffix}_${olook}rlks.int  \
    $width 0.25 2 
   cp col_${image1}-${image2}${suffix}_${olook}rlks.int.rsc \
   filt_col_${image1}-${image2}${suffix}_${olook}rlks.int.rsc

   # SW filt with col
   rm -f filtSW_col_${image1}-${image2}${suffix}_${olook}rlks.int
   # threshold on col: 0.1, filter windows: 12
   nsb_SWfilter.pl col_${image1}-${image2}${suffix}_${olook}rlks filtSW_col_${image1}-${image2}${suffix}_${olook}rlks \
   $cor $SWwindowsize $SWamplim $filterstyle
   cp ${prefix}${image1}-${image2}${suffix}_${olook}rlks.int.rsc filtSW_col_${image1}-${image2}${suffix}_${olook}rlks.int.rsc

   length.pl filt_col_${image1}-${image2}${suffix}_${olook}rlks.int
   length.pl filtSW_col_${image1}-${image2}${suffix}_${olook}rlks.int

   echo '......................'
   echo 'Unwrapping'
   echo '......................'

   ## unwrapping

   my_deroul_interf_filt filtSW_col_${image1}-${image2}${suffix}_${olook}rlks.int \
   cut \
   col_${image1}-${image2}${suffix}_${olook}rlks.int \
   filt_col_${image1}-${image2}${suffix}_${olook}rlks.int \
   $seedx $seedy $thresholdcol $thresholdfiltSW 0

   # deroul_interf_filt filtSW_col_${image1}-${image2}${suffix}_${olook}rlks.int \
   # cut \
   # col_${image1}-${image2}${suffix}_${olook}rlks.int \
   # filt_col_${image1}-${image2}${suffix}_${olook}rlks.int \
   # $seedx $seedy $thresholdfiltSW 0

   # rm -f *msk *flg
   # echo " " > cut
   # make_mask.pl filt_${prefix}${image1}-${image2}${suffix}_${olook}rlks filt_${prefix}${image1}-${image2}${suffix}_${olook}rlks_msk 0.02
   # new_cut.pl filt_${prefix}${image1}-${image2}${suffix}_${olook}rlks
   # unwrap.pl filt_${prefix}${image1}-${image2}${suffix}_${olook}rlks filt_${prefix}${image1}-${image2}${suffix}_${olook}rlks_msk \
   # filt_${prefix}${image1}-${image2}${suffix}_${olook}rlks 0.2 $seedx $seedy

   cp filt_col_${image1}-${image2}${suffix}_${olook}rlks.int.rsc \
   filt_col_${image1}-${image2}${suffix}_${olook}rlks.unw.rsc
   length.pl filt_col_${image1}-${image2}${suffix}_${olook}rlks.unw

   rm -f temp filt_col_${image1}-${image2}${suffix}_${olook}rlks_bridge.int

   echo '......................'
   echo 'Add empirical stratified model'
   echo '......................'

   add_rmg.py --infile="filt_col_${image1}-${image2}${suffix}_${olook}rlks.unw" --outfile="filt_col_${image1}-${image2}_sd_flatrange_4rlks.unw" --add="${image1}-${image2}_strat_${olook}rlks.unw"
   # add_rmg.pl filt_col_${image1}-${image2}${suffix}_${olook}rlks.unw ${image1}-${image2}_strat_${olook}rlks.unw filt_col_${image1}-${image2}_sd_flatrange_4rlks.unw  1
   cd ..

done





