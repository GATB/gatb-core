#!/bin/bash

#***************************************************************************
# GATB-Core management script: compile source code using a Docker container.
#
# When running inside its Docker container, we expect to find
# two working paths (inside the container):
#  /tmp/gatb-core-code: will contain the git clone of GATB-Core
#  /tmp/gatb-core-build: will contain the compiled code 
#
# The script also uses the environment variable: GIT_BRANCH, than can
# be set to an appropriate GATB-Core branch. If not set, we work on
# the master branch.
#
#
# Author: Patrick Durand, Inria
# Created: February 2017
#****************************************************************************

#set -xv

GCORE_SOURCE=/tmp/gatb-core-code
GCORE_BUILD=/tmp/gatb-core-build

[ ! -d ${GCORE_SOURCE} ] && { echo "${GCORE_SOURCE} does not exist. Abort."; exit 1; }
[ ! -d ${GCORE_BUILD} ] && { echo "${GCORE_BUILD} does not exist. Abort."; exit 1; }

if [ -z ${GIT_PROVIDER} ]; then
  GIT_PROVIDER="hub"
fi

# git management not done for Jenkins/CI (Inria only)
if [ ! ${GIT_PROVIDER} == "ci" ]; then
  # Figure out whether or not we have to get GATB-Core source code
  cd ${GCORE_SOURCE}
  if [ ! -d "gatb-core" ]; then
      git clone https://github.com/GATB/gatb-core.git 
  fi

  # Update GATB-Core then checkout appropriate branch
  cd gatb-core
  git pull --all

  # set default branch to master if not specified otherwise using
  # GIT_BRANCH environment variable
  if [ -z ${GIT_BRANCH} ]; then
    GIT_BRANCH=master
  fi

  git checkout ${GIT_BRANCH}
fi

# Prepare a fresh build directory
cd ${GCORE_BUILD}
if [ -d "build" ]; then
    rm -rf build
fi
mkdir build 
cd build 

# Compile source code
cmake ${GCORE_SOURCE}/gatb-core/gatb-core 
make 
