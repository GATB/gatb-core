#!/bin/bash

#**************************************************************************
# GATB-Core management script: test compiled code using a Docker container.
#
# When running inside its Docker container, we expect to find
# a working path (inside the container):
#  /tmp/gatb-core-build: will contain the compiled code 
#
# Author: Patrick Durand, Inria
# Created: February 2017
#**************************************************************************

#set -xv

GCORE_BUILD=/tmp/gatb-core-build/build
GCORE_SOURCE=/tmp/gatb-core-code

[ ! -d ${GCORE_BUILD} ] && { echo "${GCORE_BUILD} does not exist. Abort."; exit 1; }

# Run unit tests
cd ${GCORE_BUILD}

cp -r ${GCORE_SOURCE}/gatb-core/gatb-core/test/db ./test/ 
export CPPUNIT_VERBOSE=1 
cd bin 
./gatb-core-cppunit
