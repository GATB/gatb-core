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
DO_NOT_STOP_AT_ERROR : ${DO_NOT_STOP_AT_ERROR}

-----------------------------------------
 Jenkins build parameters (built in)
-----------------------------------------
BUILD_NUMBER         : ${BUILD_NUMBER}
JOB_NAME             : ${JOB_NAME}


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
#df -kh
#---------------------------------------------------------------


################################################################
#                       COMPILATION                            #
################################################################

JENKINS_TASK=${JOB_NAME}

MACHINE="`hostname`"
case $MACHINE in
koriscale*)
  echo $MACHINE
  GIT_DIR=/home/ci-gatb/workspace/$JENKINS_TASK/gatb-core
  BUILD_DIR=/home/ci-gatb/scratchdir/$JENKINS_TASK/gatb-core/build
  ;;
gatb-core-ubuntu16-docker)
  echo $MACHINE
  GIT_DIR=/builds/workspace/$JENKINS_TASK/gatb-core
  BUILD_DIR=/scratchdir/$JENKINS_TASK/gatb-core/build
  ;;
*)
  echo Erreur
  exit 1
  ;;
esac

rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR

cd $BUILD_DIR

#---------------------------------------------------------------
cmake -Wno-dev -DJENKINS_TAG=${BRANCH_TO_BUILD} -DJENKINS_GFORGE_USER=${INRIA_FORGE_LOGIN} $GIT_DIR

#---------------------------------------------------------------
make -j 2 doc || error_code
make deploy-doc || error_code

################################################################
#                         END                                  #
################################################################
