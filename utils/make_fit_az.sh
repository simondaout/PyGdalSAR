#!/bin/sh

# Initial checkings
if [ "$1" = ""] || ["$2" = ""] || [ "$4" = ""] || ["$5" = ""] || \
   ["$6" = ""] [ "$7" = ""] || ["$8" = ""] [ "$9" = ""] || ["$10" = ""] || \
   ["$11" = ""] || ["$12" = ""]; then
	echo "usage : make_fit_az.sh home_dir liste_interfero.rsc suffix prefix rlook \
   jstart jend left top lenthr width r"
	echo "home_dir: full path to T??? dir "
	echo "liste_interfero.rsc: path fom home dir to rsc file"
   echo "suffix prefix int: ${prefix}$date1-$date2_${suffix}"
	echo "jstart jend: define a zone at zero for the flattening"
	echo "left top lenthr width: define the referencing band"
   echo "r: limits for quad. and cubic fit in az are (1+r)*width and (1+2*r)*width"
	echo ; exit         
fi 

home=$1
list_int=$home/$2
suffix=$3
prefix=$4
rlook=$5
jstart=$6
jend=$7
left=$8
top=$9
lengthr=$10
widthr=$11
r=$12

right=$(echo "scale=3; $left + $width" | bc )

#cat < $list_int | while true
#do
#   read ligne
#   if [ "$ligne" = "" ]; then break; fi
#   set -- $ligne ; image1=$1 ; image2=$2 ;
#   cd $home/int/int_$1_$2/
#   
#   echo int_$1_$2
#   length.pl ${prefix}$1-$2_${suffix}_${rlook}rlks.unw
#done

baseline="$home/baseline.rsc"

#the limits for quad. and cubic fit in az are (1+r)*nx and (1+2*r)*nx
#r:10 strong to avoid quadratic ramps 
invert_rampe_az $list_int $baseline ${prefix} _${suffix}_${rlook}rlks $jstart $jend $left $right $r no  
invers_fit_az << END
$baseline
liste_fit_az
1
0
END

cat < $list_int | while true
do
   read ligne
   if [ "$ligne" = "" ]; then break; fi
   set -- $ligne ; image1=$1 ; image2=$2 ;
   cd $home/int/int_$1_$2/
   
   cp ${prefix}$1-$2_${suffix}_${rlook}rlks.unw.rsc  ${prefix}$1-$2_${suffix}_flata_${rlook}rlks.unw.rsc
   correct_fit_az ${prefix}$1-$2_${suffix}_${rlook}rlks.unw ${prefix}$1-$2_${suffix}_flata_${rlook}rlks.unw
   length.pl ${prefix}$1-$2_${suffix}_flata_${rlook}rlks.unw

   refer_interf "${prefix}$1-$2_${suffix}_flata_${rlook}rlks.unw" "${prefix}$1-$2_${suffix}_flata_refer_${rlook}rlks.unw"  $left $top $widthr $lengthr
   cp ${prefix}$1-$2_${suffix}_${rlook}rlks.unw.rsc ${prefix}$1-$2_${suffix}_flata_refer_${rlook}rlks.unw.rsc
   length.pl ${prefix}$1-$2_${suffix}_flata_refer_${rlook}rlks.unw

   cd ..
done
