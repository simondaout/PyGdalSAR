#!/bin/bash

# Initial checkings
#if [ "$1" = ""] || ["$2" = ""] || [ "$3" = ""] || ["$4" = ""] || \
#   ["$5" = ""] [ "$6" = ""] || ["$7" = ""] [ "$8" = ""] || ["$9" = ""]; then
#	echo "usage : refer.sh home_dir liste_interfero.rsc suffix prefix rlook \
#   left top lenthr width "
#	echo "home_dir: full path to T??? dir "
#	echo "liste_interfero.rsc: path fom home dir to rsc file"
#   echo "suffix prefix int: ${prefix}$date1-$date2_${suffix}"
#	echo "left top lenthr width: define the referencing band"
#	echo ; exit         
#fi 

home=$1
list_int=$home/$2
suffix=$3
prefix=$4
rlook=$5
left=$6
top=$7
lengthr=$8
widthr=$9

cat < $list_int | while true
do
   read ligne
   if [ "$ligne" = "" ]; then break; fi
   set -- $ligne ; image1=$1 ; image2=$2 ;
   cd $home/INTERFERO/int_$1_$2/

   length.pl ${prefix}$1-$2_${suffix}_${rlook}rlks.unw
   refer_interf "${prefix}$1-$2_${suffix}_${rlook}rlks.unw" "${prefix}$1-$2_${suffix}_refer_${rlook}rlks.unw" $left $top $widthr $lengthr
   cp ${prefix}$1-$2_${suffix}_${rlook}rlks.unw.rsc ${prefix}$1-$2_${suffix}_refer_${rlook}rlks.unw.rsc
   length.pl ${prefix}$1-$2_${suffix}_refer_${rlook}rlks.unw

done
