#!/bin/bash

#*****************************************************************************************
# GATB-CORE management script. Only for use from Delivery.cmake script.
#
# Check that we are preparing a gatb-core build targeting an Inria Forge repository. 
# 
# Author: Patrick Durand, Inria
# Created: January 2016
#*****************************************************************************************

l_script_dir=$( cd -P -- "$(dirname -- "$(command -v -- "$0")")" && pwd -P )
git_config="$l_script_dir/../../.git/config"

if [ ! -e $git_config ]; then
  echo "Git configuration file not found: $git_config"
  exit 1
fi

cat $git_config | grep --quiet "scm.gforge.inria.fr"

if [ $? -eq 1 ];then
  echo "This Git repository is not targeting Inria Forge."
  echo "TODO_migration_gitlab: we do not stop however"
  #exit 2
fi


