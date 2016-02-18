#!/bin/bash
set +xv
echo "
date      : `date`
hostname  : `hostname`
pwd       : `pwd`

--------------------------
Jenkins build parameters
--------------------------
BRANCH_TO_BUILD   : ${BRANCH_TO_BUILD}
RELEASE_TO_BUILD  : ${RELEASE_TO_BUILD}
INRIA_FORGE_LOGIN : ${INRIA_FORGE_LOGIN}
TEST_VARIABLE     : ${TEST_VARIABLE}
"

set -xv

################################################################
#       CHECKOUT $BRANCH_TO_BUILD FROM SCRATCH                 #
################################################################

gcc --version
g++ --version

[ `gcc -dumpversion` = 4.7 ] && { echo "GCC 4.7"; } || { echo "GCC version is not 4.7, we exit"; exit 1; }

JENKINS_TASK=finalize-release-debian
ROOT_DIR=/builds/workspace/$JENKINS_TASK/

rm -rf $ROOT_DIR
mkdir -p $ROOT_DIR

cd $ROOT_DIR
git clone git+ssh://${INRIA_FORGE_LOGIN}@scm.gforge.inria.fr/gitroot/gatb-core/gatb-core.git 
cd gatb-core
git checkout $BRANCH_TO_BUILD
cd gatb-core   # this level should be removed soon

################################################################
#                       MAKE DELIVERY                          #
################################################################

MAJOR="`echo $RELEASE_TO_BUILD |cut -d. -f1`"
MINOR="`echo $RELEASE_TO_BUILD |cut -d. -f2`"
PATCH="`echo $RELEASE_TO_BUILD |cut -d. -f3`"

/bin/bash -xv ./scripts/make_official_release.sh -s -M $MAJOR -m $MINOR -p $PATCH -c delivery
