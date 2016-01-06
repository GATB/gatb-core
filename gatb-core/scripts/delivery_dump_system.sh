#!/bin/bash

#*****************************************************************************************
# GATB-CORE management script. Only for use from Delivery.cmake script.
#
# Dump some system properties. 
# 
# Author: Patrick Durand, Inria
# Created: January 2016
#*****************************************************************************************

# ARGUMENTS:
# $1:   file where to dump properties
# $2:   cmake version
# $3:   cmake_system_name
# $4:   cmake_system
# $5:   cmake_system_processor
# $6:   compiler_id
# $7:   compiler_version
# $8:   cxx_flags
# $9:   lib_flags

echo "The library for $3 has been compiled as follows:" > $1
echo "CMake $2" >> $1
echo "OS: $4 running on $5 processor" >> $1
echo "Compiler: $6 $7" >> $1
echo "Compiler flags: $8" >> $1
echo "Library flags: $9" >> $1
