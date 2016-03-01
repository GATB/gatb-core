##  /usr/bin/env bash

#*****************************************************************************************
# Github management script.
#
# Script relying on Github API v3, as described in:
# https://developer.github.com/v3/repos/releases/
#
# This script contains an API and is intended to be used as an include within other 
# scripts. Methods of interest are:
#     getGithubReleaseId()
#     getGithubReleaseDescription()
#     uploadAsset()
#     deleteAsset()
#     listReleaseSummary()
#     listAssetSummary()
#     deleteRelease()
# For more information, see methods description below, as well as their use within the
# main section of the script.
#
# Author: Patrick Durand, Inria
# Created: December 2015
#*****************************************************************************************

# ========================================================================================
# Section: prepare access to JSON script
s_dir=$( cd -P -- "$(dirname -- "$(command -v -- "$0")")" && pwd -P )
JSON_SH=$s_dir/json-v2.sh

# ========================================================================================
# Section: variable declarations

#for internal use, to store github api answers
github_answer="_ga_.json"
#some global variables to store credentials and github repos
LOGIN=""
TOKEN=""
OWNER=""
REPOSITORY=""
#command to execute
COMMAND=""
#git tag name
TAG=""
#list of files from command-line ([file ...] part if any)
FILES=""
#default message used to create a new release
CREATE_MESSAGE="new release"
#by default this script talk to you
SILENT="off"

# ========================================================================================
# Section: utility function declarations

# --------
# FUNCTION: print out an simple message on stderr (only if SILENT mode is off)
function errorMsg(){
  if [ "$SILENT" == "off" ]; then
    printf "$* \n" >&2
  fi
}

# --------
# FUNCTION: print out an simple message on stdout (only if SILENT mode is off)
function infoMsg(){
  if [ "$SILENT" == "off" ]; then
    printf "$* \n"
  fi
}

# --------
# FUNCTION: print out an error message on stderr and exit from application
function throw () {
  errorMsg "$* \n"
  exit 77
}

# --------
# FUNCTION: check whether or not github returns an error
#   arg1: file containing json answer
#   return: error msg to display
function isError(){
  local github_error_msg=$(cat "$1" | $JSON_SH -b | grep -F -e "[message]" | cut -s -f 2 | sed -e "s/\"//g" | tr [:upper:] [:lower:])

  if [ "$github_error_msg" == "not found" ]; then 
    echo "true"
  else
    echo "false"
  fi
}

# --------
# FUNCTION: clean any resources created
function clean(){
  if [ -e $github_answer ]; then
    rm -f $github_answer
  fi
}

# --------
# FUNCTION: return a key associated to a particular key from a json file
#   arg1: file containing json answer from github
#   arg2: key. Can be a simple string or structured keys should have the form A.B.C when
#         a json file has a structured contents. 
#   return: a string containing the value or an empty string if not found
function getDataField(){
  echo $(cat "$1" | $JSON_SH -b | grep -F -e "[$2]" | cut -s -f 2 | sed -e "s/\"//g")
}

# ========================================================================================
# Section: github_API_v3 based function declarations

# --------
# FUNCTION: return a github release id given a git release tag
#   arg1: git release tag
#   return: a string containing release id 
#   throw: exit application if release does not exist
function getGithubReleaseId(){
  # Connect github to check whether or not release already exists
  curl --user ${LOGIN}:${TOKEN} \
     --request GET \
     --output "$github_answer" \
     --silent \
     --data @- \
     https://api.github.com/repos/${OWNER}/${REPOSITORY}/releases/tags/${1} <<END
END

  #do we have an error ?
  if [ $(isError "$github_answer") == "true" ]; then
	  throw "/!\ Cannot find release github ID for git tag ${1}"
  fi

  echo $(getDataField "$github_answer" "id")

}

# --------
# FUNCTION: return a github release description
#   arg1: github release ID
#   return: a json containing release description 
#   throw: exit application if release does not exist
function getGithubReleaseDescription(){
  # Connect github to check whether or not release already exists
  curl --user ${LOGIN}:${TOKEN} \
     --request GET \
     --output "$github_answer" \
     --silent \
     --data @- \
     https://api.github.com/repos/${OWNER}/${REPOSITORY}/releases/${1} <<END
END

  #do we have an error ?
  if [ $(isError "$github_answer") == "true" ]; then
	  throw "/!\ Cannot find github ID ${1}"
  fi

  echo $(cat "$github_answer")
}


# --------
# FUNCTION: return a list of assets that are contained in a specific release
#   arg1: github release id
#   return: a json structured contents
function getAssetsDescription(){
  curl --user ${LOGIN}:${TOKEN} \
     --request GET \
     --output "$github_answer" \
     --silent \
     --data @- \
     https://api.github.com/repos/${OWNER}/${REPOSITORY}/releases/${1}/assets <<END
END
  echo $(cat "$github_answer")
}

# --------
# FUNCTION: return a list of release that are contained in a specific repository
#   args: none
#   return: a json structured contents
function getReleasesDescription(){
  curl --user ${LOGIN}:${TOKEN} \
     --request GET \
     --output "$github_answer" \
     --silent \
     --data @- \
     https://api.github.com/repos/${OWNER}/${REPOSITORY}/releases <<END
END
  echo $(cat "$github_answer")
}

# --------
# FUNCTION: display the release available on github for a particular repository
#   args: none
#   return: nothing
function listReleaseSummary(){
  local dRelName="" 
  local dRelDate="" 
  local dRelId="" 

  infoMsg "Releases available on github for ${OWNER}/${REPOSITORY}:"
  local rels=$($(echo getReleasesDescription $1) | $JSON_SH -b | cut -s -d "." -f 1 | sed -e "s/\[//g" | uniq)
  if [ ! -z "$rels" ]; then
    for key in $rels; 
      do
        dRelName=$(getDataField "$github_answer" "${key}.tag_name")
        dRelId=$(getDataField "$github_answer" "${key}.id")
        dRelDate=$(getDataField "$github_answer" "${key}.created_at")
        infoMsg " asset $key: $dRelName ($dRelId), $dRelDate";
      done
  else
    infoMsg "  none."
  fi
}

# --------
# FUNCTION: display the assets (id and name) contained in a release
#   arg1: github release id
#   return: nothing
function listAssetSummary(){
  local dFileName="" 
  local dFileSize="" 
  local dFileDate="" 
  local dFileDownload=""

  infoMsg "File(s) for release $1:"
  local files=$($(echo getAssetsDescription $1) | $JSON_SH -b | cut -s -d "." -f 1 | sed -e "s/\[//g" | uniq)
  if [ ! -z "$files" ]; then
    for key in $files; 
      do
        dFileName=$(getDataField "$github_answer" "${key}.name")
        dFileSize=$(getDataField "$github_answer" "${key}.size")
        dFileDate=$(getDataField "$github_answer" "${key}.created_at")
        dFileDownload=$(getDataField "$github_answer" "${key}.download_count")
        
        infoMsg " asset $key: $dFileName ($dFileSize bytes), $dFileDate; $dFileDownload downloads.";
      done
  else
    infoMsg "  none."
  fi
}

# --------
# FUNCTION: delete an asset by its file name. Utility method, use deleteAsset() instead.
#   arg1: github release id
#   arg2: file name to delete
#   return: 0 if success (file found and deleted), 1 otherwise
function deleteAssetbyName(){
  local dFileName="" 
  local dFileId="" 
  for key in $($(echo getAssetsDescription $1) | $JSON_SH -b | cut -s -d "." -f 1 | sed -e "s/\[//g" | uniq); 
  do
    dFileName=$(getDataField "$github_answer" "${key}.name")
    if [ "$dFileName" == "$2" ]; then
      dFileId=$(getDataField "$github_answer" "${key}.id")
      curl --user ${LOGIN}:${TOKEN} \
           --request DELETE \
           --output "$github_answer" \
           --silent \
           --data @- \
           https://api.github.com/repos/${OWNER}/${REPOSITORY}/releases/assets/$dFileId <<END
END
      return 0
    fi
  done
  return 1
}

# --------
# FUNCTION: delete a file on remote github repository
#   arg1: github release id
#   arg2: file name to delete
#   return: nothing
function deleteAsset(){
  infoMsg "Deleting $2"
  if $(deleteAssetbyName "$1" "$2"); 
    then 
      infoMsg "  deleted."; 
    else 
      errorMsg "  not deleted (maybe not found)"; 
  fi
}

# --------
# FUNCTION: upload a file to github
#   arg1: github release id
#   arg2: file to upload
#   return: nothing
#   throw: exit application if file cannot be found of if upload failed
function uploadAsset(){
  local release_desc=""
  local tmpFile=""
  local dField=""
  local dField2=""
  local upURL=""
  local upFile=""
  local upBaseName=""
  
  # Connect github to get upload_url from release description
  release_desc=$(getGithubReleaseDescription $1)
  tmpFile="zmtp.json"
  echo $release_desc > $tmpFile
  dField=$(getDataField "$tmpFile" "upload_url")
  rm -f $tmpFile
  # Connect github to upload a file (binary mode only)
  upURL=$(echo $dField | cut -d "{" -f 1)
  upURL="${upURL}?name="
  upFile=$2
  infoMsg "Uploading files on github for release $1"
  infoMsg "Upload URL: ${upURL}"
  infoMsg "  uploading file: ${upFile}"
  upBaseName=$(basename $upFile)
  curl --user ${LOGIN}:${TOKEN} \
     --request POST \
     --silent \
     --output "$github_answer" \
     --header "Content-Type: application/octet-stream" \
     --data-binary @${upFile} \
     ${upURL}${upBaseName} <<END
END

  dField=$(getDataField "$github_answer" "state")
  if [ ! "$dField" == "uploaded" ]; then
	  dField=$(getDataField "$github_answer" "message")
	  dField2=$(getDataField "$github_answer" "errors.0.code")
	  throw "/!\ Cannot upload file. ${dField}. ${dField2}."
  fi
  infoMsg "    uploaded."
}

# --------
# FUNCTION: upload a file to github
#   arg1: tag name
#   arg2: release description (optional)
#   return: nothing
#   throw: exit if tag already exist as a release on github side, or if release 
#          creation failed.
function createRelease(){
  local release_desc="new release"
  local dField=""
  # Connect github to create release
  infoMsg "Connecting github to create release: $1"

  if [ $# -eq 2 ]; then
    release_desc=$2 
  fi
  curl --user ${LOGIN}:${TOKEN} \
     --request POST \
     --output "$github_answer" \
     --silent \
     --data @- \
     https://api.github.com/repos/${OWNER}/${REPOSITORY}/releases <<END
{
 "tag_name": "$1",
 "target_commitish": "master",
 "name": "$1",
 "body": "$release_desc",
 "draft": false,
 "prerelease": false
}
END

  #do we find that release ?
  dField=$(getDataField "$github_answer" "errors.0.code")
  if [ "$dField" == "already_exists" ]; then
	  throw "/!\ Release ${1} already exists, cannot create it again."
  fi

  #if created, we should have a github ID, show it!
  dField=$(getDataField "$github_answer" "id")
  if [ ! -z "$dField" ]; then
	  infoMsg "  created with github ID: $dField"
  else
    dField=$(getDataField "$github_answer" "message")
	  throw "  Failed. $dField"
  fi
}

# --------
# FUNCTION: delete a github release
# !!!! se with caution: there is no undo !!!!
# arg1: github release id
# return: nothing
function deleteRelease(){
  infoMsg "Deleting github release: $1"
  curl --user ${LOGIN}:${TOKEN} \
       --request DELETE \
       --output "$github_answer" \
       --silent \
       --data @- \
       https://api.github.com/repos/${OWNER}/${REPOSITORY}/releases/${1} <<END
END
  infoMsg "  done."
}

# --------
# FUNCTION: check if mandatory arguments are provided
#  args: none
#  return: nothing
#  throw: exit application if value is empty
function checkMandatoryArg(){
  local params=( "-l" "-t" "-o" "-r" "-c" )
  local vars=( "$LOGIN" "$TOKEN" "$OWNER" "$REPOSITORY" "$COMMAND" )

  for ((i=0;i<${#params[@]};++i)); do
    #printf "  %s %s\n" "${params[i]}" "${vars[i]}"
    if [ -z "${vars[i]}" ]; then
      throw "/!\ Missing mandatory argument: ${params[i]}"
    fi
  done
}

# --------
# FUNCTION: check if command is valid
#  arg1: command name (a string)
#  return: nothing
#   throw: exit application if value is unknown
function checkCommand(){
  local cmds=( "create" "flist" "rlist" "upload" "delete" "info" "erase")
  
  if [[ "${cmds[*]}" =~ "$1" ]]; then
    return 0
  fi
  throw "/!\ Unknown command: '$1'. Valid command is one of: ${cmds[*]}"
}

# --------
# FUNCTION: check if we have a git tag for commands that require it
function checkTag(){
  if [ -z "$TAG" ]; then
    throw "/!\ Missing mandatory argument: -d"
  fi
}

# --------
# FUNCTION: check if we have some files for command delete/upload
#  arg1: command name (a string)
#  arg2: a string containing some file names
#  return: nothing
#   throw: exit application if we do not provide some files for commands delete/upload
function checkFileList(){
  local cmds=( "upload" "delete")
  
  if [[ "${cmds[*]}" =~ "$1" ]]; then
    if [ -z "$2" ]; then
      throw "/!\ no file(s) provided for command: $1"
    fi
    for fname in $2
    do
      if [ ! -f "$fname" ]; then
        throw "/!\ not a file: $fname"
      fi
    done
  fi
}

# ========================================================================================
# Section : tips and tricks...

# these two lines from solution 12 on http://unix.stackexchange.com/questions/48533/exit-shell-script-from-a-subshell
# to enable exit to work within the throw() method.
set -E
trap '[ "$?" -ne 77 ] || exit 77' ERR

