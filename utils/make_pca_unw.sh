#!/bin/sh

if [ "$1" = "" ] || [ "$2" = "" ] ; then
  echo "utilisation : make_acp_unw.sh home list_int prefix suffix jstrat jend istart iend"
  echo; exit
fi

home=$1
list_int=$home/$2
prefix=$3
suffix=$4
jstart=$5 
jend=$6
istart=$7
iend=$8

##
baseline="$home/baseline_top_bot.rsc"

## compute ACP
ACP_unw $list_int $prefix $suffix $jstart $jend $istart $iend

# temp invers
invers_fit_frange  << END
$baseline
acp_vecteurs
0
1
1
1
END
mv inverted_cst inverted_acp1
mv inverted_lin inverted_acp2
mv inverted_quad inverted_acp3
mv inverted_cub inverted_acp4
rm inverted_cst_interf
rm inverted_lin_interf
rm inverted_quad_interf
rm inverted_cub_interf
