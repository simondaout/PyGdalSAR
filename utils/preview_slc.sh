#!/bin/sh

# Initial checkings
if [ "$1" = ""] || ["$2" = ""] || ["$3" = ""] ; then
	echo "usage : previeuw_slc.sh home_dir baselines.rsc dir"
	echo "home_dir: full path to T??? dir "
	echo "baslines.rsc: path from list images in first columnn"
  echo "dlook: downsample look of the output jpeg"
	echo "dir: path to output dir form home"
  echo ; exit         
fi 

home=$1
list_im=$home/$2
outdir=$home/$3

cat < $list_im | while true
do
   read ligne
   if [ "$ligne" = "" ]; then break; fi
   set -- $ligne ; image=$1
   echo $1

   # look.pl $1/$1_coreg.slc $dlook
   nsb_preview_slc $1/$1_coreg.slc $1/$1_coreg.jpeg
done

mkdir $outdir
mv 20*/*jpeg $outdir/.

