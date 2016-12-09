#!/bin/bash

#*****************************************************************************************
# GATB-CORE: new tool script.
#
# This script is used to create a new software relying on GATB-CORE library.
#
# Usage: run script with -h as argument.
#
# Author: Patrick Durand, Inria
# Created: February 2016
#*****************************************************************************************

# This script directory (absolute path)
_scripts_dir=$( cd -P -- "$(dirname -- "$(command -v -- "$0")")" && pwd -P )

# Some global variables
INCLUDE_GATB_COPY="false"
INSTALL_DIRECTORY=""
TOOL_NAME=""
GATB_CORE_PATH="${_scripts_dir}/../.."
NB_TOOLS=1
PROJECT_DIR=""
GATB_CMAKE_PATH=""
ERR=0

# Some error messages and corresponding exit codes
ERR_MSG_2=("/!\ Error: tool installation directory not provided"      "2")
ERR_MSG_3=("/!\ Error: tool name not provided"                        "3")
ERR_MSG_4=("/!\ Error: Directory for GATB-CORE submodule not found:"  "4")
ERR_MSG_5=("/!\ Error: Directory for new project already exists:"     "5")
ERR_MSG_6=("/!\ Error: unable to prepare source template"             "6")
ERR_MSG_7=("/!\ Error: unable to copy scripts"                        "7")
ERR_MSG_8=("/!\ Error: unable to prepare CMake files"                 "8")
ERR_MSG_9=("/!\ Error: unable to prepare GATB-Core"                   "9")
ERR_MSG_10=("/!\ Error: unable to copy other files"                  "10")
ERR_MSG_11=("/!\ Error: unable to create project directory:"         "11")

# ========================================================================================
# Section: utility function declarations
# --------
# FUNCTION: display help message
function help(){
	printf "\n[$0]: a tool to create a new software relying on GATB-Core library.\n\n"
  printf "usage: [-h] [-t <number_of_binaries>] -d <install_directory> -n <tool_name>\n\n"
  printf "  -d <install_directory>   : the directory where to place the new software; use absolute path.\n"
  printf "  -n <tool_name>           : the tool name; use alpha-num characters only.\n"
  printf "  -t <number_of_binaries>  : number of binaries making the tool (default is 1).\n"
  printf "  -h                       : this message.\n\n"
  printf "This script will automatically create the path made of <install_directory>/<tool_name>. For these\n"
  printf "two values, we recommand you do not use space characters.\n\n"
  printf "From time to time, you may need to organize your tool using more than one binaries. In such a case\n"
  printf "you may be interested in organizing your source codes using a dedicated source folder, one for each\n"
  printf "binary. In such a case, use -t <number_of_binaries>.\n"
  exit 0
}

# --------
# FUNCTION: creates directory to host tool source codes
#  args: step number
#  return: nothing
function createToolSourceDirectory(){
  ERR=0
  printf "$1- Preparing template source code...\n"
  if [[ $NB_TOOLS == 1 ]]; then
    mkdir -p $PROJECT_DIR/tools/${TOOL_NAME}/src  || ERR=1
    cat $_scripts_dir/main.cpp | sed s/XXX/${TOOL_NAME}/g  > $PROJECT_DIR/tools/${TOOL_NAME}/src/main.cpp  || ERR=1
    cat $_scripts_dir/XXX.hpp  | sed s/XXX/${TOOL_NAME}/g  > $PROJECT_DIR/tools/${TOOL_NAME}/src/${TOOL_NAME}.hpp  || ERR=1
    cat $_scripts_dir/XXX.cpp  | sed s/XXX/${TOOL_NAME}/g  > $PROJECT_DIR/tools/${TOOL_NAME}/src/${TOOL_NAME}.cpp  || ERR=1
  else
    for ((i=0;i<${NB_TOOLS};++i)); do
      # We create the tools directory and two sub directory for two different tools
      mkdir -p $PROJECT_DIR/tools/${TOOL_NAME}_${i}/src  || ERR=1
      cat $_scripts_dir/main.cpp | sed s/XXX/${TOOL_NAME}_${i}/g  > $PROJECT_DIR/tools/${TOOL_NAME}_${i}/src/main.cpp  || ERR=1
      cat $_scripts_dir/XXX.hpp  | sed s/XXX/${TOOL_NAME}_${i}/g  > $PROJECT_DIR/tools/${TOOL_NAME}_${i}/src/${TOOL_NAME}_${i}.hpp  || ERR=1
      cat $_scripts_dir/XXX.cpp  | sed s/XXX/${TOOL_NAME}_${i}/g  > $PROJECT_DIR/tools/${TOOL_NAME}_${i}/src/${TOOL_NAME}_${i}.cpp  || ERR=1
    done
  fi
  printf "   done\n"
}

# --------
# FUNCTION: creates directory to host tool scripts
#  args: step number
#  return: nothing
function createToolScriptDirectory(){
  ERR=0
  printf "$1- Copying utility scripts...\n"
  mkdir -p $PROJECT_DIR/scripts  || ERR=1
  # only copy required scripts (do not use *.sh: we cannot include NewProject.sh into new tool!)
  cp  $_scripts_dir/package_tool.sh  $PROJECT_DIR/scripts/  || ERR=1
  printf "   done\n"
}

# --------
# FUNCTION: copy and update CMake files
#  args: step number
#  return: nothing
function createCMakeFiles(){
  ERR=0
  printf "$1- Preparing CMake files...\n"
  # We copy the top level CMakeLists.txt  
  cat $_scripts_dir/CMakeLists.txt | sed s/XXX/${TOOL_NAME}/g | sed "s|GGGGG|${GATB_CMAKE_PATH}|g" > $PROJECT_DIR/CMakeLists.txt  || ERR=1

  # We copy the 'tools' CMakeLists.txt
  cat $_scripts_dir/CMakeLists_tools.txt | sed s/XXX/${TOOL_NAME}/g  > $PROJECT_DIR/tools/CMakeLists.txt  || ERR=1
  printf "   done\n"
}

# --------
# FUNCTION: creates directory to host a hard copy of gatb-core 
#  args: step number
#  return: nothing
function createGatbCoreDirectory(){
    ERR=0
    if [ "$INCLUDE_GATB_COPY" == "true" ] ; then
      printf "$1- Including hard copy of GATB-Core library...\n"
      mkdir -p $PROJECT_DIR/thirdparty/gatb-core  || ERR=1
      cp -r $GATB_CORE_PATH/cmake             $PROJECT_DIR/thirdparty/gatb-core/  || ERR=1
      cp -r $GATB_CORE_PATH/doc               $PROJECT_DIR/thirdparty/gatb-core/  || ERR=1
      cp -r $GATB_CORE_PATH/examples          $PROJECT_DIR/thirdparty/gatb-core/  || ERR=1
      cp -r $GATB_CORE_PATH/src               $PROJECT_DIR/thirdparty/gatb-core/  || ERR=1
      cp -r $GATB_CORE_PATH/thirdparty        $PROJECT_DIR/thirdparty/gatb-core/  || ERR=1
      cp -r $GATB_CORE_PATH/tools             $PROJECT_DIR/thirdparty/gatb-core/  || ERR=1
      cp    $GATB_CORE_PATH/CMakeLists.txt    $PROJECT_DIR/thirdparty/gatb-core/  || ERR=1
      cp    $GATB_CORE_PATH/LICENCE           $PROJECT_DIR/thirdparty/gatb-core/  || ERR=1
      cp    $GATB_CORE_PATH/../README.md      $PROJECT_DIR/thirdparty/gatb-core/  || ERR=1
      cp    $GATB_CORE_PATH/RELEASES.md       $PROJECT_DIR/thirdparty/gatb-core/  || ERR=1
      cp    $GATB_CORE_PATH/THIRDPARTIES.md   $PROJECT_DIR/thirdparty/gatb-core/  || ERR=1
    else
      printf "$1- Preparing GATB-Core link...\n"
    fi
    printf "   done\n"
}

# --------
# FUNCTION: copy and/or update other files
#  args: step number
#  return: nothing
function createOtherFiles(){
  ERR=0
  local MSG1="Place here the documentation of your tool, then remove this file."
  local MSG2="Place here any other third-party librairies this tool relies on."
  local MSG3="Place here test suites: scripts and/or small/medium sized data files."
  
  printf "$1- Copying other files...\n"
  # We copy the default README
  cat $_scripts_dir/README.md | sed s/XXX/${TOOL_NAME}/g > $PROJECT_DIR/README.md  || ERR=1
  cp $GATB_CORE_PATH/LICENCE                               $PROJECT_DIR/           || ERR=1
  mkdir $PROJECT_DIR/doc                                                           || ERR=1
  echo $MSG1 | tee $PROJECT_DIR/doc/README                                         || ERR=1
  if [[ ! -e $PROJECT_DIR/thirdparty ]]; then
     mkdir $PROJECT_DIR/thirdparty                                                 || ERR=1
  fi
  echo $MSG2 | tee $PROJECT_DIR/thirdparty/README                                  || ERR=1
  cp $GATB_CORE_PATH/THIRDPARTIES.md                       $PROJECT_DIR/thirdparty || ERR=1
  mkdir $PROJECT_DIR/tests                                                         || ERR=1
  echo $MSG3 | tee $PROJECT_DIR/tests/README                                       || ERR=1
  printf "   done\n"
}

# --------
# FUNCTION: display final msg to provide user with some help to start with the project
#  args: none
#  return: nothing
function displayEndingMsg(){
  printf "\n*** SUCCESS: project '$TOOL_NAME' has been created in '$PROJECT_DIR'\n\n"
  printf "    To edit source code: have a look at $PROJECT_DIR/tools/$TOOL_NAME/src\n\n"
  printf "    To compile source code:\n"
  printf "       cd $PROJECT_DIR \n"
  printf "       mkdir build ; cd build ; cmake .. \n"
  printf "           (check that cmake did not report any errors)\n"
  printf "       make -j8 \n"
  printf "           (adapt -j8 to your system; here, '8' means '8-core processor')\n\n"
  printf "    Binary of your tool will be located in $PROJECT_DIR/tools\n\n"
  printf "    More information: read $PROJECT_DIR/README.md file\n\n"
}

# ========================================================================================
# Section : main

# Prepare arguments for processing
while getopts ht:d:n: opt
do
    case "$opt" in
      h)  help;;
      d)  INSTALL_DIRECTORY="$OPTARG";;
      n)  TOOL_NAME="$OPTARG";;
      t)  NB_TOOLS="$OPTARG";;
      \?)	help;;
    esac
done
shift `expr $OPTIND - 1`

# Define the new project directory
PROJECT_DIR="$INSTALL_DIRECTORY/$TOOL_NAME"

printf "\n>>> Start making new tool '$TOOL_NAME' within '$PROJECT_DIR'\n\n"

# Check mandatory arguments
[ -z "$INSTALL_DIRECTORY" ]   && { echo ${ERR_MSG_2[0]}; exit ${ERR_MSG_2[1]}; }
[ -z "$TOOL_NAME" ]           && { echo ${ERR_MSG_3[0]}; exit ${ERR_MSG_3[1]}; }

# Check that we found the GATB-CORE submodule directory
[ ! -d "$GATB_CORE_PATH" ]    && { echo "${ERR_MSG_4[0]} $GATB_CORE_PATH"; exit ${ERR_MSG_4[1]}; }

# Check that the new project does not already exist
[ -d "$PROJECT_DIR" ]         && { echo "${ERR_MSG_5[0]} $PROJECT_DIR";    exit ${ERR_MSG_5[1]}; }

# Create the project directory
(mkdir -p "$PROJECT_DIR")     || { echo "${ERR_MSG_11[0]} $PROJECT_DIR";   exit ${ERR_MSG_11[1]}; }

# Check how this script is called to figure out how to handle GATB-Core dependency
ABS_PATH=`cd "$INSTALL_DIRECTORY"; pwd` 
if [[ "$ABS_PATH" =~ "gatb-tools/tools" ]]; then
    # if the new tool is located in gatb-tools: "link" tool to the gatb-core included 
    # within gatb-tools
    INCLUDE_GATB_COPY="false"
    GATB_CMAKE_PATH="/../../thirdparty/gatb-core/gatb-core"
else
    # we the new tool is outside gatb-tools: make a hard copy of gatb-core
    INCLUDE_GATB_COPY="true"
    GATB_CMAKE_PATH="/thirdparty/gatb-core"
fi

# Prepare the new tool
createToolSourceDirectory 1
[ $ERR == 1 ] && { echo ${ERR_MSG_6[0]}; exit ${ERR_MSG_6[1]}; }

createToolScriptDirectory 2
[ $ERR == 1 ] && { echo ${ERR_MSG_7[0]}; exit ${ERR_MSG_7[1]}; }

createCMakeFiles 3
[ $ERR == 1 ] && { echo ${ERR_MSG_8[0]}; exit ${ERR_MSG_8[1]}; }

createGatbCoreDirectory 4
[ $ERR == 1 ] && { echo ${ERR_MSG_9[0]}; exit ${ERR_MSG_9[1]}; }

createOtherFiles 5
[ $ERR == 1 ] && { echo ${ERR_MSG_10[0]}; exit ${ERR_MSG_10[1]}; }

# success!
displayEndingMsg
