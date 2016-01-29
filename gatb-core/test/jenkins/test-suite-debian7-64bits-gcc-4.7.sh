#!/bin/bash
#--------------------------------------------------------------#
#         Continuous integration script for Jenkins            #
#--------------------------------------------------------------#
#
# Default mode : 
# This script will exit with error (exit code 1) if any of its steps fails. 
# To change this behaviour, launch the script with the DEBUG argument.
#--------------------------------------------------------------#
 
[ "$1" != "DEBUG" ] && { set -e ; } || { echo "DEBUG mode, the script will NOT stop..." ; }
set -xv

date
hostname
pwd

# Jenkins Build parameters
echo "GIT_TAG : $GIT_TAG"

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

# Note: if "dgbh5 -check" fails, exit code will be 1 (0 otherwise), and the Jenkins build will be reported as FAILED

$BUILD_DIR/bin/dbgh5 -verbose 0 -in $HOME/reference/fastq/aphid_662451seq.fa               -check $HOME/reference/check/aphid_662451.props

$BUILD_DIR/bin/dbgh5 -verbose 0 -in $HOME/reference/fastq/aphid_662451seq.album/album.txt  -check $HOME/reference/check/aphid_662451.props

$BUILD_DIR/bin/dbgh5 -verbose 0 -in $HOME/reference/fastq/SRR959239_clean.fastq.gz         -check $HOME/reference/check/SRR959239_clean.props

