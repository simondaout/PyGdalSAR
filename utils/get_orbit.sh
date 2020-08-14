#!/bin/bash

PROC_FILE="nsbas.proc"
num_images=`grep SarDir $PROC_FILE | wc | awk '{print $1}'`
max_try=1
count_images=1

while [ $count_images -le $num_images ]
  do
  NameDate=`grep SarDir $PROC_FILE | awk 'NR=='$count_images'' | awk 'BEGIN{FS="="} {print $2}'`
  number_try=0
  while [ $number_try -le $max_try ]
  do
	  echo $NameDate/${NameDate}.slc
	  touch $NameDate/${NameDate}.slc
          
	  /home/cometsoft/Ubuntu-18.04/nsbas//scripts/nsb_fetch_s1_orb.py -s $NameDate/${NameDate}.slc -t  RESORB -v 4
	  /home/cometsoft/Ubuntu-18.04/nsbas//scripts/nsb_fetch_s1_orb.py -s $NameDate/${NameDate}.slc -t  POEORB -v 4
	  number_try=$(( $number_try+1 ))
	  echo $number_try
  done
  count_images=$(( $count_images+1 ))
done


