#!/bin/bash

#*****************************************************************************************
# GATB-CORE management script.
#
# This script is used to prepare an official release of GATB-CORE library.
#
# Usage: run script with -h as argument.
#
# Author: Patrick Durand, Inria
# Created: December 2015
#*****************************************************************************************

CRED_FILE="$HOME/.gatbcore_manager"
SILENT="false"
GH_OWNER="GATB"
GH_REPO="gatb-core"
GH_LOGIN=""
GH_TOKEN=""
IF_LOGIN="$INRIA_FORGE_LOGIN"          # known inside Jenkins tasks, empty otherwise
MAJOR_V=""
MINOR_V=""
PATCH_V=""
COMMAND=""

# ========================================================================================
# Section: utility function declarations
# --------
# FUNCTION: display help message
function help(){
	printf "\n[$0]: a tool to prepare an official GATB-CORE release.\n\n"
  printf "usage: [-h] [-s] [-l <credential_file>] -M <version_major> -m <version_minor> -p <version_patch> -c <command> \n\n"
  printf "  -l <credential_file> : optional credential file (see below).\n"
  printf "  -M <tag>             : major release number.\n"
  printf "  -m <tag>             : minor release number.\n"
  printf "  -p <tag>             : patch number.\n"
  printf "  -h                   : this message.\n"
  printf "  -s                   : silent mode (not set means interactive mode).\n"
  printf "  -c <command>         : command to execute.\n"
  printf "   'command' is one of delivery, upload or package.\n"
  printf "      delivery: makes a full release processing as follows:\n"
	printf "         1. compile the library on the hosted system using curent source code \n"
	printf "         2. package library within a tar.gz file \n"
	printf "         3. tag Inria Forge using release number \n"
	printf "         4. create a corresponding Github release \n"
	printf "         5. upload the library binary on Github \n"
  printf "      upload: runs only steps 1, 2 and 5.\n"
  printf "      package: runs only steps 1 and 2.\n\n"
  printf "About the credential file. It is intended to contain login:token pair of values\n"
  printf "required to connect to remote github repository. Only authorized users are able\n"
  printf "to use this script to make a release on github.\n"
  printf "File format is quite simple: it contains two lines:\n"
  printf "  - first line is:  login=<value>\n"
  printf "  - second line is: token=<value>\n"
  printf "where values have to be replaced accordingly.\n"
  printf "When '-l' argument is not provided, this script looks for a file called \$HOME/.gatbcore_manager\n\n"
  exit 0
}

# --------
# FUNCTION: print out a simple message on stderr 
function errorMsg(){
  printf "$* \n" >&2
}

# --------
# FUNCTION: check if mandatory arguments are provided
#  args: none
#  return: nothing
#  exit application if value is empty
function checkMandatoryArg(){
  local params=( "-M" "-m" "-p" "-c" )
  local vars=( "$MAJOR_V" "$MINOR_V" "$PATCH_V" "$COMMAND" )

  for ((i=0;i<${#params[@]};++i)); do
    if [ -z "${vars[i]}" ]; then
      errorMsg "/!\ Missing mandatory argument: ${params[i]}"
      exit 1
    fi
  done
}

# --------
# FUNCTION: read the credential file to retrieve github login:token pair of values
#  args: none
#  return: nothing
#  exit application if credential file is missing or if login/token cannot be retrieved
function readCrendentialFile(){
  local value=""
  if [ ! -e $CRED_FILE ];then
    errorMsg "/!\ Crendential file not found: $CRED_FILE\n"
    exit 2
  fi
  GH_LOGIN=$( cat $CRED_FILE | grep "login" | cut -d = -s -f 2)
  if [ -z $GH_LOGIN ]; then
    errorMsg "/!\ Missing 'login' field. Check your credential file: $CRED_FILE"
    exit 3
  fi
  GH_TOKEN=$( cat $CRED_FILE | grep "token" | cut -d = -s -f 2)
  if [ -z $GH_TOKEN ]; then
    errorMsg "/!\ Missing 'token' field. Check your credential file: $CRED_FILE"
    exit 4
  fi
}

# --------
# FUNCTION: check if command is valid
#  arg1: command name (a string)
#  return: nothing
#   exit application if value is unknown
function checkCommand(){
  local cmds=( "delivery" "upload" "package")
  
  if [[ "${cmds[*]}" =~ "$1" ]]; then
    return 0
  fi
  errorMsg "/!\ Unknown command: '$1'. Valid command is one of: ${cmds[*]}"
  exit 5
}

# ========================================================================================
# Section : main

# Prepare arguments for processing
while getopts hlsM:m:p:c: opt
do
    case "$opt" in
      s)  SILENT="true";;
      l)  CRED_FILE="$OPTARG";;
      M)  MAJOR_V="$OPTARG";;
      m)  MINOR_V="$OPTARG";;
      p)  PATCH_V="$OPTARG";;
      c)  COMMAND="$OPTARG";;
      h)  help;;
      \?)	help;;
    esac
done
shift `expr $OPTIND - 1`

# Ask users if we can go ahead
if [ "$SILENT" == "false" ]; then
  echo "***"
  echo "*** Making official GATB-CORE library"
  echo "***"
  echo " "
  echo ">>>   DID YOU CHECK THAT EVERYTHING IS COMMITED TO GIT ?   <<<"
  echo " "
  echo " "
  echo "If ok, press enter to continue... (CTRL+C to abort)"
  echo " "
  read 
fi

#check provided values. Exit script if something is wrong.
readCrendentialFile
checkMandatoryArg
checkCommand $COMMAND

#now, start the whole job
script_dir=$( cd -P -- "$(dirname -- "$(command -v -- "$0")")" && pwd -P )

cd $script_dir/..
rm -rf build
mkdir build
cd build

cmake -DGH_LOGIN=$GH_LOGIN -DGH_TOKEN=$GH_TOKEN -DGH_OWNER=$GH_OWNER \
  -DGH_REPO=$GH_REPO -DMAJOR=$MAJOR_V -DMINOR=$MINOR_V -DPATCH=$PATCH_V \
  -DSILENT_MODE=$SILENT -DCPACK_USER_NAME=$IF_LOGIN ..

make $COMMAND
