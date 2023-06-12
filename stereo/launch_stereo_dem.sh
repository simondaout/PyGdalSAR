#!/usr/bin/bash

module load asp
. ./asp_parameters.txt

# set input DEM
DEM=$DEM_FILE
DEM_UTM_FILE=$(basename $DEM .tif)"_utm.tif"

##################
# TILED IMAGES   #
##################
# With some Airbus Pleiades data, each of the left and right images may arrive broken up into .TIF or .JP2 tiles, with names ending in R1C1.tif, R2C1.tif, etc.
# These need to be mosaicked before being used
#gdalbuildvrt $DATA_DIR"/"$NAME1"/IMG_PHR1A_P_001/vrt.tif" $DATA_DIR"/"$NAME1"/IMG_PHR1A_P_001/*R*C*.TIF"
#gdalbuildvrt $DATA_DIR"/"$NAME2"/IMG_PHR1A_P_001/vrt.tif" $DATA_DIR"/"$NAME2"/IMG_PHR1A_P_001/*R*C*.TIF"

#actual self-contained image can be produced with:
#gdal_translate -co TILED=YES -co BLOCKXSIZE=256 -co BLOCKYSIZE=256 -co BIGTIFF=IF_SAFER $DATA_DIR"/"$NAME1"/IMG_PHR1A_P_001/vrt.tif" $IMG1
#gdal_translate -co TILED=YES -co BLOCKXSIZE=256 -co BLOCKYSIZE=256 -co BIGTIFF=IF_SAFER $DATA_DIR"/"$NAME2"/IMG_PHR1A_P_001/vrt.tif" $IMG2

mkdir $OUTPUT_DIR
# copy parameter file to new working dir
cp asp_parameters.txt $OUTPUT_DIR
# change to working dir; TIMESTAMP bc this is the name of the current working dir
cd $OUTPUT_DIR

##################
# BUNDLE ADJUST #
##################

bundle_adjust $IMG1 $IMG2 $Lrpc $Rrpc -t $SESSION_TYPE --camera-weight 0 --tri-weight 0.1  --datum wgs84 -o ba/run 

# ###############
# # Map Project #
# ###############

mapproject -t $SESSION_TYPE --t_srs EPSG:$UTM --tr $RESMP $DEM $IMG1 $Lrpc $IMG1_MP --bundle-adjust-prefix ba/run --nodata-value 0
mapproject -t $SESSION_TYPE --t_srs EPSG:$UTM --tr $RESMP $DEM $IMG2 $Rrpc $IMG2_MP --bundle-adjust-prefix ba/run --nodata-value 0

###########
# REF DEM #
###########

gdalwarp -tr $GDAL_OUT_RES -t_srs EPSG:$UTM $DEM $DEM_UTM_FILE -r $RESAMP_M -overwrite

##############
# RUN STEREO #
##############

# if filtering is selected in asp_parameters.txt, the following parameters will be included for parallel_stereo: rm-quantile-percentile, rm-quantile-multiple, rm-cleanup-passes, filter-mode, rm-half-kernel, rm-min-matches, rm-threshold

if [ $FILTERING = true ]
then
parallel_stereo -t $SESSION_TYPE --alignment-method $A_M $IMG1_MP $IMG2_MP $Lrpc $Rrpc demPleiades/dem $DEM_UTM_FILE  --bundle-adjust-prefix ba/run --nodata-value $NO_DATA_S --prefilter-mode $PREF_MODE --prefilter-kernel-width $PREF_KER_M --corr-kernel $CORR_KERNEL --cost-mode $COST_MODE --stereo-algorithm $ST_ALG --corr-tile-size $CORR_T_S --subpixel-mode $SUBP_MODE --subpixel-kernel $SUBP_KERNEL --corr-seed-mode $CORR_S_MODE --processes $THREADS --xcorr-threshold $XCORR_TH --min-xcorr-level $MIN_XCORR_LVL --sgm-collar-size $SGM_C_SIZE --rm-quantile-percentile $RM_QUANT_PC --rm-quantile-multiple $RM_QUANT_MULT --rm-cleanup-passes $RM_CLEAN_PASS --filter-mode $FILTER_MODE --rm-half-kernel $RM_HALF_KERN --rm-min-matches $RM_MIN_MAT --rm-threshold $RM_TH

else
parallel_stereo -t $SESSION_TYPE --alignment-method $A_M $IMG1_MP $IMG2_MP $Lrpc $Rrpc demPleiades/dem $DEM_UTM_FILE --bundle-adjust-prefix ba/run --nodata-value $NO_DATA_S --corr-kernel $CORR_KERNEL --cost-mode $COST_MODE --stereo-algorithm $ST_ALG --corr-tile-size $CORR_T_S --subpixel-mode $SUBP_MODE --subpixel-kernel $SUBP_KERNEL --corr-seed-mode $CORR_S_MODE --processes $THREADS --xcorr-threshold $XCORR_TH --min-xcorr-level $MIN_XCORR_LVL --sgm-collar-size $SGM_C_SIZE --prefilter-mode $PREF_MODE --prefilter-kernel-width $PREF_KER_M
fi

#############
# POINT2DEM #
#############

# the dem-PC file is stored in demPleiades; the final dem-DEM.tif will be created within demPleiades

cd demPleiades

point2dem --t_srs EPSG:$UTM --tr $RES dem-PC.tif --median-filter-params $MED_F_PAR --dem-hole-fill-len $DEM_HOLE_F_L --erode-length $ERODE_L --nodata-value $NO_DATA_DEM --tif-compress $TIF_COMPR --max-valid-triangulation-error $MAX_V_TRIANG_ERR --remove-outliers-param $RM_OUTL_PARA --threads $THREADS

exit

