#!/usr/bin/env bash

#*****************************************************************************************
# Git tag management script.
#
# This script can be used to tag a gatb-core release. When a tag is added, it is 
# automatically pushed to the remote server. 
# 
# Usage:
#   use option -h
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
  printf "usage: $0 [-h] [-D] [-t <message>] -M <major> -m <minor> -p <patch>\n"
  printf "\n"
  printf "  -t <message> : new tag message (optional). Default message: 'new release'.\n"
  printf "  -M <tag>     : major release number.\n"
  printf "  -m <tag>     : minor release number.\n"
  printf "  -p <tag>     : patch number.\n"
  printf "  -D           : delete existing tag.\n"
  printf "  -h           : this message.\n"
  printf "\n"
  printf "Using values from '-M A', '-m B' and '-p C' arguments, the script creates a tag\n"
  printf "named 'vA.B.C' (without quotes). The tag is automatically pushed to the remote\n"
  printf "server.\n"
  printf "\n"
  printf "Return value is one of:\n"
  printf "  0: ok\n"
  printf "  1: missing mandatory argument\n"
  printf "  2: tag version already exists\n"
  printf "  3: failed to create tag version\n"
  printf "  4: failed to push tag to remote repository\n"
  printf "  5: failed to delete tag\n"
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
CMD=""

# Prepare arguments for processing
while getopts hDM:m:p:t: opt
do
    case "$opt" in
      M)  MAJOR="$OPTARG";;
      m)  MINOR="$OPTARG";;
      p)  PATCH="$OPTARG";;
      t)  CREATE_MESSAGE="$OPTARG";;
      D)  CMD="delete";;
      h)  help;;
      \?)	help;;
    esac
done
shift `expr $OPTIND - 1`

#Do we have all mandatory arguments ?
mandatory_params=( "-M" "-m" "-p" )
mandatory_values=( "$MAJOR" "$MINOR" "$PATCH" )

for ((i=0;i<${#mandatory_params[@]};++i)); do
  if [ -z "${mandatory_values[i]}" ]; then
    printf "/!\ Missing mandatory argument: ${mandatory_params[i]}\n" >&2
    printf "    use option -h to get help.\n" >&2
    exit 1
  fi
done

# Prepare tag
TAG=$(echo "v${MAJOR}.${MINOR}.${PATCH}")

# Check whether provided tag already exists on repository
curTag=`git tag -l $TAG`
if [ ! -z "$curTag" ]; then
  #if tag exists and we do not want to delete it: error
  if [ ! "$CMD" == "delete" ]; then
    printf "/!\ git tag '$TAG' already exist.\n" >&2
    exit 2
  fi
  #delete tag from local repository...
  git tag --delete $curTag
  #... then from remote repository
  git push --delete origin $curTag
  exit 0
fi

#do we have to delete tag?
if [ "$CMD" == "delete" ]; then
  printf "/!\ git tag '$TAG' does not exist: nothing to delete.\n" >&2
  exit 5
fi


# Prepare message
if [ ! -z "$CREATE_MESSAGE" ]; then
  release_msg=$CREATE_MESSAGE
else
  release_msg="new release"
fi


# Create tag on local repository and check if it's ok
printf "Tagging git repository with \n"
printf "   tag: $TAG\n"
printf "   msg: $release_msg\n"
git tag -m "$release_msg" $TAG

OUT=$?
if [ $OUT -eq 0 ];then
  printf "git tag '$TAG' created.\n"
else
  printf "/!\ unable to create git tag '$TAG'.\n" >&2
  exit 3
fi

# Push tag to the remote server and check if it's ok
git push origin $TAG
OUT=$?
if [ $OUT -eq 0 ];then
  printf "git tag '$TAG' pushed to remote repository.\n"
  exit 0
else
  printf "/!\ git tag '$TAG' not pushed to remote repository.\n" >&2
  exit 4
fi

