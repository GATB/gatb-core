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

error_code () { [ "$DO_NOT_STOP_AT_ERROR" = "true" ] && { return 0 ; } }

[ "$DO_NOT_STOP_AT_ERROR" != "true" ] && { set -e ; } || { echo "DEBUG mode, the script will NOT stop..." ; }
set -xv

# quick look at resources
#-----------------------------------------------
sw_vers -productVersion
#-----------------------------------------------
system_profiler SPSoftwareDataType
#-----------------------------------------------
lstopo
#-----------------------------------------------
top -l 1|head -15
#-----------------------------------------------


################################################################
#                       COMPILATION                            #
################################################################

clang --version
clang++ --version
cmake --version

export CC=/usr/bin/clang
export CXX=/usr/bin/clang++

cd gatb-core

JENKINS_TASK=${JOB_NAME}
GIT_DIR=/builds/workspace/$JENKINS_TASK/gatb-core
#BUILD_DIR=/scratchdir/$JENKINS_TASK/gatb-core/build
BUILD_DIR=$GIT_DIR/build         #N.B. /scratchdir not yet mounted on the osx slave (ciosx)

rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR

cd $BUILD_DIR

cmake -Wno-dev -DJENKINS_TAG=${BRANCH_TO_BUILD} $GIT_DIR

make -j 2 || error_code

################################################################
#                       UNIT TESTS                             #
################################################################
export CPPUNIT_VERBOSE=1

# Copy database for unit tests
cp -r $GIT_DIR/test/db $BUILD_DIR/test/

# Specify single unit tests
#$BUILD_DIR/bin/gatb-core-cppunit TestBag
#$BUILD_DIR/bin/gatb-core-cppunit TestMap

# Launch the full test suite
cd $BUILD_DIR/bin
ls ../test/db/           # default directory for test db
./gatb-core-cppunit || error_code

################################################################
#    CHECK FUNCTIONS (with precomputed reference results)      #
################################################################

# Not ready

################################################################
#                   VALGRIND CHECK                             #
################################################################

# not ready
