#!/bin/bash

echo ""
echo "--------------------------------------------------------------------------------"
echo "          CREATE PROJECT '"$1"' in directory '"$2"'"
echo "--------------------------------------------------------------------------------"
echo ""


# We get the directory of the script (absolute path)
_scripts_dir=$( cd -P -- "$(dirname -- "$(command -v -- "$0")")" && pwd -P )

# We define the new project directory
_project_dir="$2/$1"

# We set the GATB-CORE directory
_gatb_core_dir="$_scripts_dir/../.."

# We check that we found the GATB-CORE submodule directory
[ ! -d "$_gatb_core_dir" ] && { echo "Error: Directory for GATB-CORE submodule not found"; exit 2; }

# We check that the new project doesn't already exist
[ -d "$_project_dir" ] && { echo "Error: Directory '$_project_dir' for new project already exists"; exit 2; }

# We create the project directory
mkdir $_project_dir

# We create the project source directory
mkdir $_project_dir/src

# We copy gatb-core as thirdparty
mkdir $_project_dir/thirdparty
mkdir $_project_dir/thirdparty/gatb-core
cp -r $_gatb_core_dir/cmake             $_project_dir/thirdparty/gatb-core/
cp -r $_gatb_core_dir/src               $_project_dir/thirdparty/gatb-core/
cp -r $_gatb_core_dir/tools             $_project_dir/thirdparty/gatb-core/
cp -r $_gatb_core_dir/thirdparty        $_project_dir/thirdparty/gatb-core/
cp -r $_gatb_core_dir/examples          $_project_dir/thirdparty/gatb-core/
cp -r $_gatb_core_dir/CMakeLists.txt    $_project_dir/thirdparty/gatb-core/

# We init the CMakeLists.txt
touch $_project_dir/CMakeLists.txt

# Note: we do it this way because it is cumbersome to do it with sed because of special characters in _gatb_core_dir

echo "project($1)"                                                                          >> $_project_dir/CMakeLists.txt
echo ""                                                                                     >> $_project_dir/CMakeLists.txt
echo "cmake_minimum_required(VERSION 2.6)"                                                  >> $_project_dir/CMakeLists.txt
echo ""                                                                                     >> $_project_dir/CMakeLists.txt
echo "################################################################################"     >> $_project_dir/CMakeLists.txt
echo "# GATB-CORE"                                                                          >> $_project_dir/CMakeLists.txt
echo "################################################################################"     >> $_project_dir/CMakeLists.txt
echo ""                                                                                     >> $_project_dir/CMakeLists.txt
echo "# we depend on gatb-core; here, we define where to find all the required material"    >> $_project_dir/CMakeLists.txt
echo "add_subdirectory (thirdparty/gatb-core \"\${CMAKE_CURRENT_BINARY_DIR}/ext/gatb-core\") "   >> $_project_dir/CMakeLists.txt

# We copy the remaining of the file
cat $_scripts_dir/CMakeLists.txt | sed 's/__PROJECT_NAME__/'$1'/g' >> $_project_dir/CMakeLists.txt

# We copy the default main.cpp file
cp $_scripts_dir/main.cpp $_project_dir/src/main.cpp

# We copy the default README
cp $_scripts_dir/README.md $_project_dir/README.md

# We go to the new project
echo "PROJECT CREATED IN DIRECTORY '$_project_dir'"
echo ""
