#!/usr/bin/env bash

#*****************************************************************************************
# Git tag management script.
#
# This script can be used to tag a gatb-core release. When a tag is added, it is 
# automatically pushed to the remote server. 
#
# Author: Patrick Durand, Inria
# Created: December 2015
#*****************************************************************************************

# --------
# FUNCTION: print a spacer
function printspacer(){
	sname=$0 ; length=${#sname} ; ((length+=2))
	for ((i=0;i<length;i++)); do 
	  printf " " 
	done
}
# --------
# FUNCTION: display help message
function help(){
	printf "\n$0: a tool to handle git release tag creation using canonical \n"
	printspacer
	printf "numbering including major, minor and patch numbers.\n\n"
  printf "usage: $0 [-h] [-t <message>] -M <major> -m <minor> -p <patch>\n"
  printf "\n"
  printf "  -t <message> : new tag message (optional). Default message: 'new release'.\n"
  printf "  -M <tag>     : major release number.\n"
  printf "  -m <tag>     : minor release number.\n"
  printf "  -p <tag>     : patch number.\n"
  printf "  -h           : this message.\n"
  exit 0
}

# Some variables
release_msg=""
curTag=""
MAJOR=""
MINOR=
PATCH=
TAG=""
OUT=""

# Prepare arguments for processing
while getopts hM:m:p:t: opt
do
    case "$opt" in
      M)  MAJOR="$OPTARG";;
      m)  MINOR="$OPTARG";;
      p)  PATCH="$OPTARG";;
      t)  CREATE_MESSAGE="$OPTARG";;
      h)  help;;
      \?)	help;;
    esac
done
shift `expr $OPTIND - 1`

#Do we have all mandatory arguments ?
mandatory_params=( "-M" "-m" "-p" )
mandatory_values=( "$MAJOR" "$MINOR" "$PATCH" )

for ((i=0;i<${#mandatory_params[@]};++i)); do
  #printf "  %s %s\n" "${mandatory_params[i]}" "${mandatory_values[i]}"
  if [ -z "${mandatory_values[i]}" ]; then
    printf "/!\ Missing mandatory argument: ${mandatory_params[i]}\n" >&2
    exit 1
  fi
done

# Prepare tag
TAG=$(echo "v${MAJOR}.${MINOR}.${PATCH}")

curTag=`git tag -l $TAG`

# Check whether provided tag already exists on repository
if [ ! -z "$curTag" ]; then
  printf "/!\ git tag '$TAG' already exist.\n" >&2
  exit 1
fi

# Prepare message
if [ ! -z "$CREATE_MESSAGE" ]; then
  release_msg=CREATE_MESSAGE
else
  release_msg="new release"
fi


# Create tag and check if it's ok
printf "Tagging git repository with '$TAG':'$release_msg'\n"
git tag -m "$release_msg" $TAG

OUT=$?
if [ $OUT -eq 0 ];then
  printf "git tag '$TAG' created.\n"
else
  printf "/!\ unable to create git tag '$TAG'.\n" >&2
  exit 1
fi

# Push tab to the server and check if it's ok
git push origin $TAG
OUT=$?
if [ $OUT -eq 0 ];then
  printf "git tag '$TAG' pushed to remote repository.\n"
  exit 0
else
  printf "/!\ git tag '$TAG' not pushed to remote repository.\n" >&2
  exit 1
fi
