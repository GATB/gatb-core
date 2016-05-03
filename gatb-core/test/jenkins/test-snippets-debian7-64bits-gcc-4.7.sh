#!/bin/bash
#--------------------------------------------------------------#
#         Continuous integration script for Jenkins            #
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
BRANCH_TO_BUILD      : ${BRANCH_TO_BUILD}
RELEASE_TO_BUILD     : ${RELEASE_TO_BUILD}
INRIA_FORGE_LOGIN    : ${INRIA_FORGE_LOGIN}
TEST_VARIABLE        : ${TEST_VARIABLE}
DO_NOT_STOP_AT_ERROR : ${DO_NOT_STOP_AT_ERROR}

-----------------------------------------
 Jenkins build parameters (built in)
-----------------------------------------
BUILD_NUMBER         : ${BUILD_NUMBER}


"


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
#                       COMPILATION                            #
################################################################

gcc --version
g++ --version

[ `gcc -dumpversion` = 4.7 ] && { echo "GCC 4.7"; } || { echo "GCC version is not 4.7, we exit"; exit 1; }

JENKINS_TASK=test-snippets-debian7-64bits-gcc-4.7
GIT_DIR=/builds/workspace/$JENKINS_TASK/gatb-core
BUILD_DIR=/scratchdir/$JENKINS_TASK/gatb-core/build

rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR

cd $BUILD_DIR

#---------------------------------------------------------------
cmake -Wno-dev $GIT_DIR

#---------------------------------------------------------------
make -j 2 || error_code

################################################################
#                   COMPILE SNIPPETS                           #
################################################################

cd $BUILD_DIR/examples

make

################################################################
#                      RUN SNIPPETS                            #
################################################################

cd $BUILD_DIR/bin/

./kmer1



################################################################
#                         END                                  #
################################################################
