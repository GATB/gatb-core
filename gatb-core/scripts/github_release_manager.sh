#!/bin/bash

#*****************************************************************************************
# Github management script.
#
# Script relying on Github API v3, as described in:
# https://developer.github.com/v3/repos/releases/
#
# This script can be used either as a command-line tool (use -h to get help) or as an API.
# For the latter use, see github_release_api.sh script.
#
# Author: Patrick Durand, Inria
# Created: December 2015
#*****************************************************************************************

# ========================================================================================
# Section: include API
script_dir=$( cd -P -- "$(dirname -- "$(command -v -- "$0")")" && pwd -P )
. $script_dir/github_release_api.sh

# ========================================================================================
# Section: utility function declarations
# --------
# FUNCTION: display help message
function help(){
	printf "\n$0: a tool to handle github release management using github/api v3.\n\n"
  printf "usage: $0 [-h] [-s] -l <login> -t <token> -o <owner> -r <repository> [-d <git_tab>] [-m <message>] -c <command> [file ...]\n\n"
  printf "Credential parameters used to access remote github repository:\n"
  printf "  -l <login> -t <token> -o <owner> -r <repository> \n"
  printf "\n"
  printf "Release identification:\n"
  printf "  -d <git_tab>\n"
  printf "\n"
  printf "Release managment commands are provided using:\n"
  printf "  -c <command> command to execute.\n"
  printf "   'command' is one of create, list, upload, delete, info.\n"
  printf "   create: create a new release.\n"
  printf "    flist: list files available for an existing release.\n"
  printf "   upload: upload file(s) to an existing release.\n"
  printf "   delete: permanently delete remote file(s) from an existing release.\n"
  printf "     info: print out some information about an existing release.\n"
  printf "    erase: permanently delete an existing release. Use with extreme caution!\n"
  printf "    rlist: list existing releases for a repository.\n"
  printf "\n"
  printf "   All commands but 'rlist' require the -d <git_tag>.\n"
  printf "\n"
  printf "   Commands 'upload' and 'delete' expect files as remaining command line arguments:\n"
  printf "       -c upload file1.tgz file2.tgz\n"
  printf "\n"
  printf "Notice:\n"
  printf "   /!\ this script does not handle file name/path containing space characters.\n"
  printf "\n"
  printf "Other arguments:\n"
  printf "  -s turn script to silent mode\n"
  printf "  -m release message: only used when creating a new release\n"
  printf "  -h display this message\n"
	exit 1
}

# ========================================================================================
# Section : main

# Prepare arguments for processing
while getopts shl:t:o:r:c:d:m: opt
do
    case "$opt" in
      s)  SILENT="true";;
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

#check we have all required arguments
checkMandatoryArg

#check the validity of the command name
checkCommand "$COMMAND"

#we need at least one file for these commands: upload, delete
checkFileList "$COMMAND" "$FILES"

#execute command
case "$COMMAND" in
  info) 
    checkTag
    release_id=$(getGithubReleaseId $TAG)
    infoMsg "Git tag $TAG refers to github release ID: $release_id"
    echo $(getGithubReleaseDescription $release_id)
    ;;
  create)
    checkTag
    createRelease $TAG "$CREATE_MESSAGE"
    ;;
  upload) 
    checkTag
    release_id=$(getGithubReleaseId $TAG)
    infoMsg "Git tag $TAG refers to github release ID: $release_id"
    for fname in $FILES
      do
        uploadAsset $release_id "$fname"
      done
    ;;
  delete) 
    checkTag
    release_id=$(getGithubReleaseId $TAG)
    infoMsg "Git tag $TAG refers to github release ID: $release_id"
    for fname in $FILES
      do
        deleteAsset $release_id "$fname"
      done
    ;;
  rlist)
    listReleaseSummary
    ;;
  flist) 
    checkTag
    release_id=$(getGithubReleaseId $TAG)
    infoMsg "Git tag $TAG refers to github release ID: $release_id"
    listAssetSummary $release_id
    ;;
  erase) 
    checkTag
    release_id=$(getGithubReleaseId $TAG)
    infoMsg "Git tag $TAG refers to github release ID: $release_id"
    deleteRelease $release_id
    ;;
esac

clean

exit 0
