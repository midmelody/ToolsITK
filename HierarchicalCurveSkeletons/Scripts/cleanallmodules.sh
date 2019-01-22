#! /bin/bash
#
#
# This script updates all the executable modules
#
#

source setPaths.sh

for module in $Modules
do
  echo "Entering directory $module"
  pushd $module > res.txt
  make clean || exit 1
  echo "Leaving directory $module"
  popd > res.txt
  rm $module/res.txt
  rm res.txt
done

## pfSkel
pushd $pfSkel_DIR
make clean
popd
