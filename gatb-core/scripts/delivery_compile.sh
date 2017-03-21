#!/bin/bash

#*****************************************************************************************
# GATB-CORE management script. Only for use from Delivery.cmake script.
#
# Compile the binary version of gatb-core from current source. 
# 
# Author: Patrick Durand, Inria
# Created: January 2016
#*****************************************************************************************

# ARGUMENTS:
# $1:   silent mode (true|flase)

# SHA 1 MANAGEMENT
export set GIT_SHA1=`git rev-parse HEAD`

# temporarely change the config_sha1.hpp file
export set CONFIG_FILE_IN="../src/gatb/system/api/config_sha1.hpp"
echo "#define STR_GIT_SHA1 " \"$GIT_SHA1\" > sha1.tmp
\mv sha1.tmp $CONFIG_FILE_IN
cat $CONFIG_FILE_IN

# clean, compile and package library
make clean
make package  

# get back the to official config_sha1.hpp
git checkout $CONFIG_FILE_IN

#
if [ "$1" == "false" ]; then
  echo "***"
  echo "*** SUCCESS: source code compiled and package created."
  echo "***"
  echo "/!\\  /!\\  /!\\  /!\\  /!\\  /!\\  /!\\  /!\\  /!\\  /!\\"
  echo " "
  echo "Now, you are about to:"
  echo "  1/ create an official tag on the Inria Forge"
  echo "  2/ create an official release on Github"
  echo "  3/ upload library binary on Github"
  echo " "
  echo "If ok, press enter to continue... (CTRL+C to abort)"
  echo " "
  read 
fi

