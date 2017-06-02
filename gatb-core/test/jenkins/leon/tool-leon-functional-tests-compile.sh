#!/bin/bash
#--------------------------------------------------------------#
#         Continuous integration script for Jenkins            #
#--------------------------------------------------------------#
#
# !!! This script has only been tested using Koriscale CI slave.
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
BRANCH_TO_BUILD      : ${BRANCH_TO_BUILD}
RELEASE_TO_BUILD     : ${RELEASE_TO_BUILD}
INRIA_FORGE_LOGIN    : ${INRIA_FORGE_LOGIN}
TEST_VARIABLE        : ${TEST_VARIABLE}
DO_NOT_STOP_AT_ERROR : ${DO_NOT_STOP_AT_ERROR}

-----------------------------------------
 Jenkins build parameters (built in)
-----------------------------------------
BUILD_NUMBER         : ${BUILD_NUMBER}
JOB_NAME             : ${JOB_NAME}

"

#---------------------------------------------------------------
# quick look at resources
free -h

################################################################
#                       COMPILATION                            #
################################################################

if $DO_NOT_COMPILE; then
  echo "SKIP COMPILE PHASE"
  exit 0
fi

# Make sure, we use the appropriate cmake on Koriscale
export PATH=/home/ci-gatb/cmake-3.7.2-Linux-x86_64/bin:$PATH

# dump compiler information
gcc --version
g++ --version
gcc -dumpversion
cmake --version

JENKINS_TASK=${JOB_NAME}
GIT_DIR=/home/ci-gatb/workspace/$JENKINS_TASK/gatb-core
BUILD_DIR=/home/ci-gatb/scratchdir/$JENKINS_TASK/gatb-core/build

mkdir -p $BUILD_DIR
cd $BUILD_DIR

#---------------------------------------------------------------
# compile in default mode: Release (otherwise Leon compressor will
# be very, very slow using Debug mode)
cmake -Wno-dev  \
      -DCPPUNIT_INCLUDE_DIR=/usr/include/cppunit/ \
      -DCPPUNIT_LIBRARY=/usr/lib64/libcppunit.so \
      $GIT_DIR

# we compile all GATB-Core: library, tools (including leon) and snippets
# (snippets bank26 to bank28 are used to test leon)
make -j8
make -j8 examples

