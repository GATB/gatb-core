#!/bin/bash
#--------------------------------------------------------------#
#         Continuous integration script for Jenkins            #
#--------------------------------------------------------------#
#
# Default mode :
# This script will exit with error (exit code 1) if any of its steps fails.
# To change this behaviour, choose DO_NOT_STOP_AT_ERROR in Jenkins (see below).
#--------------------------------------------------------------#
set +xv

echo "
-----------------------------------------
 Miscellaneous information 
-----------------------------------------
date      : `date`
hostname  : `hostname`
pwd       : `pwd`

-----------------------------------------
 Jenkins build parameters (user defined)
-----------------------------------------
VERSION_TO_TEST      : ${VERSION_TO_TEST}
INRIA_FORGE_LOGIN    : ${INRIA_FORGE_LOGIN}
DO_NOT_STOP_AT_ERROR : ${DO_NOT_STOP_AT_ERROR}

-----------------------------------------
 Jenkins build parameters (built in)
-----------------------------------------
BUILD_NUMBER         : ${BUILD_NUMBER}


"

error_code () { [ "$DO_NOT_STOP_AT_ERROR" = "true" ] && { return 0 ; } }


[ "$DO_NOT_STOP_AT_ERROR" != "true" ] && { set -e ; } || { echo "(!) DEBUG mode, the script will NOT stop..." ; echo; }
set -xv

# quick look at resources
#---------------------------------------------------------------
free -h
#---------------------------------------------------------------
lstopo
#---------------------------------------------------------------
df -kh
#---------------------------------------------------------------


################################################################
#                       DIRECTORIES                            #
################################################################

JENKINS_TASK=test-bin-debian7-64bits-gcc-4.7
BUILD_DIR=/scratchdir/$JENKINS_TASK/gatb-core/build

rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR

cd $BUILD_DIR

################################################################
#                       UNPACKING                              #
################################################################
# Upload bin bundle to the forge; source bundle is made by OSX Jenkins task
if [ $? -eq 0 ] && [ "$INRIA_FORGE_LOGIN" != none ] && [ "$DO_NOT_STOP_AT_ERROR" != true ]; then
   echo "Getting a binary archive... "
   scp ${INRIA_FORGE_LOGIN}@scm.gforge.inria.fr:/home/groups/gatb-core/htdocs/versions/bin/gatb-core-${VERSION_TO_TEST}-bin-Linux.tar.gz .
fi

################################################################
#                       COMPILATION                            #
################################################################
gcc --version
g++ --version

[ `gcc -dumpversion` = 4.7 ] && { echo "GCC 4.7"; } || { echo "GCC version is not 4.7, we exit"; exit 1; }

gunzip gatb-core-${VERSION_TO_TEST}-bin-Linux.tar.gz
tar -xf gatb-core-${VERSION_TO_TEST}-bin-Linux.tar
cd gatb-core-${VERSION_TO_TEST}-bin-Linux

code_snippets=($(find ./examples -name "*1.cpp"))
for code_snippet in ${code_snippets[*]}
do
	prg_name=`echo $code_snippet | cut -d'/' -f4 | cut -d'.' -f1`
    g++ $code_snippet -Iinclude -Llib -lgatbcore -lhdf5 -ldl -lz -lpthread  -std=c++0x -O3 -o $prg_name
done

# do some cleanup to save disk space
cd ..
rm -rf gatb-core-${VERSION_TO_TEST}-bin-Linux*

