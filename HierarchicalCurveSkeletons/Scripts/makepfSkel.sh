#! /bin/bash
#
#
# This script updates the pfSkel executable
#
#

source setPaths.sh

## pfSkel
pushd $pfSkel_DIR
make clean
make pfSkel
popd

