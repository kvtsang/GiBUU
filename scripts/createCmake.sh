#!/bin/bash
#*******************************************************************************
#****e* /createCMake.sh
# NAME
# createCmake.sh
# PURPOSE
# This script sets up a parallel directory, where the compilation is
# steered my cmake.
#
# NOTES
# It uses the Makefile targets 'renew' and 'subbdir' to build up the sources
# in the directory ./objects. The it copies these files to a parallel directory
# tree and installs CMakeLists.txt files to enable the configuration of
# the build process via cmake.
#
# run it as
#   ./scripts/createCmake.sh
#
# Open topics:
# * not all possibilities of the original Makefile a carried over.
# * setting up version.f90 with accurate svn info is not supported yet.
#*******************************************************************************

##### 1) Setup the source files in the 'objects' directory:
make renew
make subdirs

##### 2) create the parallel source directory and copy the files:
# (we are actually copying the source files the soft links are pointing to)

sDir=../release2025_cmake/trunk/src
install -d ${sDir}

for i in objects/*90 objects/*f objects/*inc;
do
    # echo $i
    f=`readlink -f $i`
    cp -a $f ${sDir}
done

cp -a version.txt ${sDir}/..

##### 3) set up the CMakeLists.txt files
cp scripts/forCmake/CMakeLists.txt.trunk  ${sDir}/../CMakeLists.txt
cp scripts/forCmake/CMakeLists.txt.src  ${sDir}/CMakeLists.txt

##### 4) set up the libraries

install -d ${sDir}/../libs

rsync -avP ../libraries/RootTuple ${sDir}/../libs
# the provided findROOT can not cope with new version numbers anymore:
rm ${sDir}/../libs/RootTuple/RootTuple-master/cmake/Modules/FindROOT.cmake
# In future we want to update the file also directly in the original tree:
cp -a scripts/forCmake/CMakeLists.txt.RootTuple_trunk  ${sDir}/../libs/RootTuple/RootTuple-master/CMakeLists.txt
cp -a scripts/forCmake/CMakeLists.txt.RootTuple_src  ${sDir}/../libs/RootTuple/RootTuple-master/src/CMakeLists.txt

rsync -avP ../libraries/HEPMC3/HEPMC3event ${sDir}/../libs
cp -a scripts/forCmake/CMakeLists.txt.HepMC3event_trunk ${sDir}/../libs/HEPMC3event/CMakeLists.txt
cp -a scripts/forCmake/CMakeLists.txt.HepMC3event_src ${sDir}/../libs/HEPMC3event/src/CMakeLists.txt


##### 99) do the compile

install -d ${sDir}/../build

# now you should 'cd' into the directory '${sDir}/../build' and run 'cmake ..'
