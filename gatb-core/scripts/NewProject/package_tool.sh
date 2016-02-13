#!/bin/bash

#*****************************************************************************************
# GATB-Tool management script.
#
# This script is used to prepare a release of a tool on your local computer.
#
# Usage: run script with -h as argument.
#
# Author: Patrick Durand, Inria
# Created: January 2016
#*****************************************************************************************

MAJOR_V=""
MINOR_V=""
PATCH_V=""

# ========================================================================================
# Section: utility function declarations
# --------
# FUNCTION: display help message
function help(){
	printf "\n[$0]: a script to prepare a release of your tool.\n\n"
  printf "usage: [-h] -M <version_major> -m <version_minor> -p <version_patch> \n\n"
  printf "  -M <tag>             : major release number.\n"
  printf "  -m <tag>             : minor release number.\n"
  printf "  -p <tag>             : patch number.\n"
  printf "  -h                   : this message.\n\n"
  printf "This tool automatically makes a fresh release of your tool: it compiles source\n"
  printf "codes on your computer, then package the binary within a tar.gz file. It also\n"
  printf "packages source codes within a separate tar.gz file.\n"
  exit 0
}

# --------
# FUNCTION: check if mandatory arguments are provided
#  args: none
#  return: nothing
#  exit application if value is empty
function checkMandatoryArg(){
  local params=( "-M" "-m" "-p")
  local vars=( "$MAJOR_V" "$MINOR_V" "$PATCH_V" )

  for ((i=0;i<${#params[@]};++i)); do
    if [ -z "${vars[i]}" ]; then
      printf "\n/!\ Missing mandatory argument: ${params[i]}\n\n" >&2
      exit 1
    fi
  done
}


# ========================================================================================
# Section : main

# Prepare arguments for processing
while getopts hM:m:p: opt
do
    case "$opt" in
      M)  MAJOR_V="$OPTARG";;
      m)  MINOR_V="$OPTARG";;
      p)  PATCH_V="$OPTARG";;
      h)  help;;
      \?)	help;;
    esac
done
shift `expr $OPTIND - 1`

#check provided values. Exit script if something is wrong.
checkMandatoryArg

#now, start the whole job
script_dir=$( cd -P -- "$(dirname -- "$(command -v -- "$0")")" && pwd -P )

cd $script_dir/..

# prepare a clean build
rm -rf build
mkdir build
cd build/
# configure the build with a release number
cmake -DMAJOR=$MAJOR_V -DMINOR=$MINOR_V -DPATCH=$PATCH_V ..
# compile code and package binary
make package -j8
# package source 
make package_source -j8

