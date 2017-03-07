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
cmake --version

[ `gcc -dumpversion` = 4.7 ] && { echo "GCC 4.7"; } || { echo "GCC version is not 4.7, we exit"; exit 1; }

JENKINS_TASK=test-valgrind-debian7-64bits-gcc-4.7
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
#                   VALGRIND CHECKS                            #
################################################################


#>>>>>>>>>>>> 1. Diagnostics for TestDebruijn >>>>>>>>>>>>>>>>>>

export CPPUNIT_VERBOSE=1

# Copy database for unit tests
cp -r $GIT_DIR/test/db $BUILD_DIR/test/

XMLFILE=$BUILD_DIR/bin/valgrind_TestDebruijn_${BUILD_NUMBER}.xml
cd $BUILD_DIR/bin

#------------ TestDebruijn without valgrind --------------------
./gatb-core-cppunit TestDebruijn || error_code


#------------ TestDebruijn with valgrind -----------------------

valgrind --xml=yes --xml-file=$XMLFILE \
	./gatb-core-cppunit TestDebruijn || error_code

xmllint --xpath "string(//error)" $XMLFILE

xmllint --xpath "string(//errorcounts)" $XMLFILE

#>>>>>>>>>>>> 2. Diagnostics for dbgh5 >>>>>>>>>>>>>>>>>>>>>>>>>

# XMLFILE=$BUILD_DIR/bin/valgrind_dbgh5_${BUILD_NUMBER}.xml
# valgrind --xml=yes --xml-file=$XMLFILE \
# $BUILD_DIR/bin/dbgh5 -verbose 0 -in $HOME/reference/fastq/aphid_662451seq.fa -check $HOME/reference/check/aphid_662451.props || error_code

# xmllint --xpath "string(//error)" $XMLFILE

# xmllint --xpath "string(//errorcounts)" $XMLFILE


################################################################
#                         END                                  #
################################################################
