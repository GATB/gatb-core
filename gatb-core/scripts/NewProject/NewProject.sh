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

# We create the tools directory and two sub directory for two different tools
mkdir $_project_dir/tools
mkdir $_project_dir/tools/$1_1
mkdir $_project_dir/tools/$1_2

# We create the source directory for two tools
mkdir $_project_dir/tools/$1_1/src
mkdir $_project_dir/tools/$1_2/src

# We create the script directory and copy the delivery script
mkdir $_project_dir/scripts
cp  $_gatb_core_dir/scripts/delivery.sh  $_project_dir/scripts/

# If no $3 argument is provided, we copy GATB-CORE as thirdparty.
# => in GATB-TOOLS context, we can provide a $3 argument and so
#    we don't have a copy of GATB-CORE
if test -z "$3"
then
    # We copy gatb-core as thirdparty
    mkdir $_project_dir/thirdparty
    mkdir $_project_dir/thirdparty/gatb-core
    cp -r $_gatb_core_dir/cmake             $_project_dir/thirdparty/gatb-core/
    cp -r $_gatb_core_dir/src               $_project_dir/thirdparty/gatb-core/
    cp -r $_gatb_core_dir/tools             $_project_dir/thirdparty/gatb-core/
    cp -r $_gatb_core_dir/thirdparty        $_project_dir/thirdparty/gatb-core/
    cp -r $_gatb_core_dir/examples          $_project_dir/thirdparty/gatb-core/
    cp -r $_gatb_core_dir/doc               $_project_dir/thirdparty/gatb-core/
    cp -r $_gatb_core_dir/CMakeLists.txt    $_project_dir/thirdparty/gatb-core/
fi

# We copy the top level CMakeLists.txt
cat $_scripts_dir/CMakeLists.txt | sed s/XXX/$1/g  > $_project_dir/CMakeLists.txt

# We copy the 'tools' CMakeLists.txt
cat $_scripts_dir/CMakeLists_tools.txt | sed s/XXX/$1/g  > $_project_dir/tools/CMakeLists.txt

# We copy the default main.cpp file
cat $_scripts_dir/main.cpp | sed s/XXX/$1_1/g  > $_project_dir/tools/$1_1/src/main.cpp
cat $_scripts_dir/XXX.hpp  | sed s/XXX/$1_1/g  > $_project_dir/tools/$1_1/src/$1_1.hpp
cat $_scripts_dir/XXX.cpp  | sed s/XXX/$1_1/g  > $_project_dir/tools/$1_1/src/$1_1.cpp

cat $_scripts_dir/main.cpp | sed s/XXX/$1_2/g  > $_project_dir/tools/$1_2/src/main.cpp
cat $_scripts_dir/XXX.hpp  | sed s/XXX/$1_2/g  > $_project_dir/tools/$1_2/src/$1_2.hpp
cat $_scripts_dir/XXX.cpp  | sed s/XXX/$1_2/g  > $_project_dir/tools/$1_2/src/$1_2.cpp


# We copy the default README
cp $_scripts_dir/README.md $_project_dir/README.md

# We go to the new project
echo "PROJECT CREATED IN DIRECTORY '$_project_dir'"
echo ""
