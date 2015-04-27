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
    
# We get the versions.txt file from the server
scp -q $4 $5

# We check that the delivery doesn't already exist.
if grep -q "$6" $5 ; then echo"" ; echo '===> THIS VERSION ALREADY EXISTS' ; echo"" ; exit 0 ; fi

# SHA 1 MANAGEMENT
export set GIT_SHA1=`git rev-parse HEAD`
export set CONFIG_FILE_IN="../src/gatb/system/api/config.hpp.in"
echo ""                    >> $CONFIG_FILE_IN
echo "#undef STR_GIT_SHA1"                 >> $CONFIG_FILE_IN
echo "#define STR_GIT_SHA1 " \"$GIT_SHA1\" >> $CONFIG_FILE_IN

# We build the package
if [[ "$1" = "SRC" ]]
    then make -j8 package_source 
    else make -j8 package  
fi

git checkout $CONFIG_FILE_IN

# We set the file rights
chmod a+r $7

# We copy the archive to the server
scp -q $7 $8

# We add information to the versions file
echo $6 >> $5  

# We set the versions.txt file to the server
scp -q $5  $4 
