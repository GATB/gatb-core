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
#                       COMPILATION                            #
################################################################

gcc --version
g++ --version

[ `gcc -dumpversion` = 4.7 ] && { echo "GCC 4.7"; } || { echo "GCC version is not 4.7, we exit"; exit 1; }

JENKINS_TASK=test-suite-debian7-64bits-gcc-4.7
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
#                       PACKAGING                              #
################################################################
# Upload bin bundle to the forge; source bundle is made by OSX Jenkins task
if [ $? -eq 0 ] && [ "$INRIA_FORGE_LOGIN" != none ] && [ "$DO_NOT_STOP_AT_ERROR" != true ]; then
   echo "Creating a binary archive... "
   echo "N.B. this is NOT an official binary release"
   make package
   scp gatb-core-${BRANCH_TO_BUILD}-bin-Linux.tar.gz ${INRIA_FORGE_LOGIN}@scm.gforge.inria.fr:/home/groups/gatb-core/htdocs/ci-inria
fi

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

#---------------------------------------------------------------
./gatb-core-cppunit || error_code
# end of unit tests

################################################################
#    CHECK FUNCTIONS (with precomputed reference results)      #
################################################################

# Note: if "dgbh5 -check" fails, exit code will be 1 (0 otherwise), and the Jenkins build will be reported as FAILED

$BUILD_DIR/bin/dbgh5 -verbose 0 -in $HOME/reference/fastq/aphid_662451seq.fa               -check $HOME/reference/check/aphid_662451.props

$BUILD_DIR/bin/dbgh5 -verbose 0 -in $HOME/reference/fastq/aphid_662451seq.album/album.txt  -check $HOME/reference/check/aphid_662451.props

$BUILD_DIR/bin/dbgh5 -verbose 0 -in $HOME/reference/fastq/SRR959239_clean.fastq.gz         -check $HOME/reference/check/SRR959239_clean.props


################################################################
#                         END                                  #
################################################################