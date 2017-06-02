#!/bin/bash

# This script tests integrity of gzipped files.
# It can be used to test SRA files retrieved from NCBI
# usind download.sh script.
#
# Patrick G. Durand, May 2017 
#

i=0
while read line
do
    array[ $i ]="$line"
    (( i++ ))
done < <(ls -1 *.fastq.gz)

echo "Nb. SRA files to test: ${#array[*]}"

echo "SRA are:"
for item in ${array[*]}
do
  echo "> Checking ${item} ..."
  FILESIZE=$(stat -c%s "${item}")
  FILESIZE=$(( FILESIZE/(1024*1024) ))
  echo "  size= $FILESIZE Mb."
  gunzip -t ${item}
  if [ $? -ne 0 ]; then
    echo "   ${item}: invalid"
  fi
done

