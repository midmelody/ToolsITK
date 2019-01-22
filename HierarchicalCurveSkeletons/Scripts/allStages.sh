#! /bin/bash
#
# This script goes through all the stages of the skeletonization process
# using all the different modules separatelly
# It uses WORK_DIR as a temporary directory for files that nedd to be exchanged
# between the modules as well as for the output 
#
# Calling this script: 
#	./allStages \
#            <path> <volumefilename> <extension> <sizeX> <sizeY> <sizeZ> \
#            <chargesDistance> <fieldStrength> <percHighDiverg> <outputDir>
# Example:
#	./allStages /home/volumes cow.100x46x68 vol 100 46 68 \
#            0 7 10 /home/outputdir
#
#

if [ $# -lt 10 ] 
then
	echo "Usage: $0 <path> <volumefilename> <extension> <sizeX> <sizeY> <sizeZ> <chargesDistance> <fieldStrength> <percHighDiverg> <outputDir>"
	exit 1
fi


#
# paths are set separatelly in this small script
#
source setPaths.sh

##
## Run all the tools
##

WORK_DIR=${10}

## make solid volume
echo Executing: $MakeSolidVol_DIR/driver $1/$2.$3 $4 $5 $6 $WORK_DIR/$2.vol
$MakeSolidVol_DIR/driver $1/$2.$3 $4 $5 $6 $WORK_DIR/$2.vol


## expand the volume with the necessary number of layers.
echo Executing: $ExpandVol_DIR/driver $WORK_DIR/$2.vol $4 $5 $6 $7 $WORK_DIR/$2.vol
$ExpandVol_DIR/driver $WORK_DIR/$2.vol $4 $5 $6 $7 $WORK_DIR/$2.vol

## make solid volume again - to make sure there are no holes in the new volume
echo Executing: $MakeSolidVol_DIR/driver $1/$2.$3 $4 $5 $6 $WORK_DIR/$2.vol
$MakeSolidVol_DIR/driver $WORK_DIR/$2.vol $4 $5 $6 $WORK_DIR/$2.vol

## potential field
echo Executing: $PotField_DIR/driver $WORK_DIR/$2.vol $4 $5 $6 $8 $WORK_DIR/$2.pf
$PotField_DIR/driver $WORK_DIR/$2.vol $4 $5 $6 $8 $WORK_DIR/$2.pf

## find critical points
echo Executing: $CritPts_DIR/driver $WORK_DIR/$2.vol $4 $5 $6 $WORK_DIR/$2.pf $WORK_DIR/$2.critpts
$CritPts_DIR/driver $WORK_DIR/$2.vol $4 $5 $6 $WORK_DIR/$2.pf $WORK_DIR/$2.critpts


## find high divergence points
echo Executing: $HighDiverg_DIR/driver $WORK_DIR/$2.vol $4 $5 $6 $WORK_DIR/$2.pf $9 $WORK_DIR/$2-$9.hd
$HighDiverg_DIR/driver $WORK_DIR/$2.vol $4 $5 $6 $WORK_DIR/$2.pf $9 $WORK_DIR/$2-$9.hd

## find streamlines
echo Executing: $StreamLn_DIR/driver $WORK_DIR/$2.pf $4 $5 $6 $WORK_DIR/$2.critpts $WORK_DIR/$2-$9.hd $WORK_DIR/$2-$9.skel
$StreamLn_DIR/driver $WORK_DIR/$2.pf $4 $5 $6 $WORK_DIR/$2.critpts $WORK_DIR/$2-$9.hd $WORK_DIR/$2-$9.skel

