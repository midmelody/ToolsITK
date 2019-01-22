#! /bin/bash
#
# This script calls the skeletonization program
# that includes all the different modules in one executable
#
# Calling this script: 
#	./pfSkel.sh <volumefilename> <sizeX> <sizeY> <sizeZ> <distCharges> <fieldStrength> <percHighDiverg> <outputFilename>
# Example:
#	./pfSkel.sh /home/user/cow.100x46x68.vol 100 46 68 2 7 50 /home/user/cow.100x46x68.skel
#
#

if [ $# -lt 8 ] 
then
	echo "Usage: $0 <volumeFilename> <sizeX> <sizeY> <sizeZ> <distCharges> <fieldStrength> <percHighDiverg> <outputFilename>"
	exit 1
fi

#
# paths are set separatelly in this small script
#
source setPaths.sh

##
## Run the program
##
$pfSkel_DIR/pfSkel $1 $2 $3 $4 $5 $6 $7 $8

