#!/usr/bin/env bash

#*****************************************************************************************
# Basic test suite for Github management script.
#
# Author: Patrick Durand, Inria
# Created: December 2015
#*****************************************************************************************

# ========================================================================================
# Include the API
script_dir=$( cd -P -- "$(dirname -- "$(command -v -- "$0")")" && pwd -P )
. $script_dir/github_release_api.sh

# ========================================================================================
# Section: utility function declarations
# --------
# FUNCTION: display help message
function help(){
	printf "\n$0: a tool to test github_release_api.\n\n"
  printf "usage: $0 [-h] -l <login> -t <token> -o <owner> -r <repository> -d <git_tab> -c <command> file ...\n\n"
	exit 1
}

# ========================================================================================
# Section: main

# Prepare arguments for processing
while getopts shl:t:o:r:c:d:m: opt
do
    case "$opt" in
      l)  LOGIN="$OPTARG";;
      t)  TOKEN="$OPTARG";;
      o)  OWNER="$OPTARG";;
      r)  REPOSITORY="$OPTARG";;
      c)  COMMAND="$OPTARG";;
      d)  TAG="$OPTARG";;
      m)  CREATE_MESSAGE="$OPTARG";;
      h)  help;;
      \?)	help;;
    esac
done
shift `expr $OPTIND - 1`

# remaining arguments, if any, are supposed to be the [file ...] part of the command-line
FILES=$@

echo "> Test createRelease()"
createRelease $TAG "my new description for $TAG"

echo "> Test getGithubReleaseId()"
release_id=$(getGithubReleaseId $TAG)
echo "$release_id"

echo "> Test getGithubReleaseDescription()"
echo $(getGithubReleaseDescription $release_id)

echo "> Test getAssetsDescription()"
assetList=$(getAssetsDescription $release_id)
echo $assetList

echo "> Test listAssetSummary()"
listAssetSummary $release_id

echo "> Test uploadAsset()"
for fname in $FILES
  do
    uploadAsset $release_id "$fname"
  done

echo "> Test deleteAsset()"
for fname in $FILES
  do
    deleteAsset $release_id "$fname"
  done

echo "> Test deleteRelease"
deleteRelease $release_id

