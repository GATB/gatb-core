#!/bin/bash

# ARGUMENTS:
# $1:   KIND ('BIN' or 'SRC')
# $2:   PROJECT_NAME
# $3:   PACKAGE_VERSION
# $4:   UPLOAD_VERSIONS
# $5:   VERSIONS_FILENAME
# $6    INFO_XXX            (XXX=BIN or SRC)
# $7:   URI_XXX             (XXX=BIN or SRC)
# $8:   UPLOAD_URI_XXX      (XXX=BIN or SRC)

echo "-----------------------------------------------------------"
echo "DELIVERY $1 FOR $2, VERSION $3"
echo "-----------------------------------------------------------"

echo ""
echo ""
echo ""
echo ""
echo "WARNING ! WARNING ! WARNING ! WARNING ! "
echo ""
echo "DID YOU CHECK THAT YOUR SOURCE IS COMPATIBLE WITH OLD VERSIONS OF GCC ?"
echo ""
echo "REMBEMBER THAT A LOT OF SERVERS STILL USE OLD LINUX DISTRIBUTIONS..."
echo ""
echo ""
echo ""
echo "press enter to continue..."
read 
    
# We get the versions.txt file from the server
scp -q $4 $5

# We check that the delivery doesn't already exist.
if [[ ! "$*" =~ "--override" ]] ; then if grep -q "$6" $5 ; then echo"" ; echo '===> THIS VERSION ALREADY EXISTS' ; echo"" ; exit 0 ; fi ; fi

# SHA 1 MANAGEMENT
export set GIT_SHA1=`git rev-parse HEAD`

# We temporarely change the config_sha1.hpp file
export set CONFIG_FILE_IN="../src/gatb/system/api/config_sha1.hpp"
echo "#define STR_GIT_SHA1 " \"$GIT_SHA1\" > sha1.tmp
\mv sha1.tmp $CONFIG_FILE_IN
cat $CONFIG_FILE_IN

# A clean won't hurt
make clean

# We build the package
if [[ "$1" = "SRC" ]]
    then make -j8 package_source 
    else make -j8 package  
fi

# We get back the official config_sha1.hpp
git checkout $CONFIG_FILE_IN

# We set the file rights
chmod a+r $7

# We copy the archive to the server
scp -q $7 $8

# We add information to the versions file
echo $6 >> $5  

# We set the versions.txt file to the server
scp -q $5  $4 
