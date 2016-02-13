#!/bin/bash
#--------------------------------------------------------------#
#         Continuous integration script for Jenkins            #
#--------------------------------------------------------------#
#
# Default mode :
# This script will exit with error (exit code 1) if any of its steps fails.·
# To change this behaviour, choose DO_NOT_STOP_AT_ERROR in Jenkins (see below).
#--------------------------------------------------------------#
set +xv
echo "
date      : `date`
hostname  : `hostname`
pwd       : `pwd`

--------------------------
 Jenkins build parameters
 --------------------------
 BRANCH_TO_BUILD      : ${BRANCH_TO_BUILD}
 RELEASE_TO_BUILD     : ${RELEASE_TO_BUILD}
 INRIA_FORGE_LOGIN    : ${INRIA_FORGE_LOGIN}
 TEST_VARIABLE        : ${TEST_VARIABLE}
 DO_NOT_STOP_AT_ERROR : ${DO_NOT_STOP_AT_ERROR}
 "

 [ "$DO_NOT_STOP_AT_ERROR" != "true" ] && { set -e ; } || { echo "DEBUG mode, the script will NOT stop..." ; }
 set -xv

################################################################
#                       COMPILATION                            #
################################################################

gcc --version
g++ --version

[ `gcc -dumpversion` = 4.2.1 ] && { echo "GCC 4.2.1"; } || { echo "GCC version is not 4.2.1, we exit"; exit 1; }

sw_vers -productVersion
#system_profiler SPSoftwareDataType

cd gatb-core

JENKINS_TASK=test-suite-macos-10.9.5-gcc-4.2.1
GIT_DIR=/builds/workspace/$JENKINS_TASK/gatb-core
#BUILD_DIR=/scratchdir/$JENKINS_TASK/gatb-core/build
BUILD_DIR=$GIT_DIR/build         #N.B. /scratchdir not yet mounted on the osx slave (ciosx)

rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR

cd $BUILD_DIR
 

 cmake -Wno-dev $GIT_DIR


 make 


################################################################
#                       UNIT TESTS                             #
################################################################
export CPPUNIT_VERBOSE=1

# Specify single unit tests
#$BUILD_DIR/bin/gatb-core-cppunit TestBag
#$BUILD_DIR/bin/gatb-core-cppunit TestMap

# Launch the full test suite
$BUILD_DIR/bin/gatb-core-cppunit

################################################################
#    CHECK FUNCTIONS (with precomputed reference results)      #
################################################################

# Not ready

[ "$DO_NOT_STOP_AT_ERROR" = "true" ] && { exit 0 ; }
