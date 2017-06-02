#!/bin/bash

# This script compares original SRA FastQ files with Leon compressed ones.
#
# Patrick G. Durand, May 2017 
#

# !!! This script only works with bash --version 4+ on Linux systems.
#     (it uses hastable, stat -c, amon others)
# !!! This script has only been tested using Koriscale CI slave.

# ==================================================================
# Where are the refaerence SRA files?
DATA_DIR=/mnt/windows/ci-gatb/data/leon
# Where are the binaries: leon and bank snippets
JENKINS_TASK=${JOB_NAME}
BIN_DIR=/home/ci-gatb/scratchdir/$JENKINS_TASK/gatb-core/build/bin
# How many cores to use? 
CORES=8
# Do we have to use leon in verbose mode?
VERBOSE=0
# default k-mer size
KMER_SIZE=31
# Bank snippets used to test Leon generated data files
BANK_SNIPPET_1=bank25
BANK_SNIPPET_2=bank26
BANK_SNIPPET_3=bank28


# ==================================================================
# == Control values
# 5 data fields: nb. letters, nb. sequences, max seq size, 
#                min seq size,  nb. sequences < k-mer size
ERR385912_valids=(139175685 2728935 51 51 0)
SRR2916693_valids=(204077299 405343 1201 67 0)
SRR959239_valids=(526537536 5372832 98 98 0)
ERR386131_valids=(544138147 3601856 371 8 14165)
ERR739513_valids=(450281548 122240 246140 5 1235)
SRR1870605_valids=(543194861 1119218 502 70 0)
SRR857303_valids=(505750634 2581532 368 6 32095)
SRR3211986_valids=(917247967 163477 62746 2 1762)
SRR2994368_valids=(2258185851 5054526 502 70 0)
SRR034509_valids=(2091430836 10353618 202 202 0)
SRR445718_valids=(3294366500 32943665 100 100 0)
SRR3190692_valids=(1300585440 9314994 602 70 0)
SRR065390_valids=(2466741904 67617092 100 100 0)
SRR1519083_valids=(1734577366 59698462 101 101 0)
ERR174310_valids=(3276346670 207579467 202 202 0)

declare -A control_map
control_map[ERR385912]=${ERR385912_valids[*]}
control_map[SRR2916693]=${SRR2916693_valids[*]}
control_map[SRR959239]=${SRR959239_valids[*]}
control_map[ERR386131]=${ERR386131_valids[*]}
control_map[ERR739513]=${ERR739513_valids[*]}
control_map[SRR1870605]=${SRR1870605_valids[*]}
control_map[SRR857303]=${SRR857303_valids[*]}
control_map[SRR3211986]=${SRR3211986_valids[*]}
control_map[SRR2994368]=${SRR2994368_valids[*]}
control_map[SRR034509]=${SRR034509_valids[*]}
control_map[SRR445718]=${SRR445718_valids[*]}
control_map[SRR3190692]=${SRR3190692_valids[*]}
control_map[SRR065390]=${SRR065390_valids[*]}
control_map[SRR1519083]=${SRR1519083_valids[*]}
control_map[ERR174310]=${ERR174310_valids[*]}

# ==================================================================
# == Usefull variables
# store size of SRA reference files (.fastq.gz)
declare -A ref_file_size_map
# use to check whether or not some tests failed
#  (0: ok ; !=0: not ok)
TEST_RESULT=0
# set a dedicated time format
TIMEFORMAT='    Time - real:%3lR | user:%3lU | sys:%3lS'

# ==================================================================
# == Usefull methods
checkDataFile(){ 
  # takes two argument: 
  #  $1 : the file base name without its extension
  #  $2 : the file name to test
  item=$1
  srafile=$2
  echo "    $BIN_DIR/$BANK_SNIPPET_3 -in ${srafile} -kmer-size $KMER_SIZE"
  RESULT=`time $BIN_DIR/$BANK_SNIPPET_3 -in ${srafile} -kmer-size $KMER_SIZE`
  echo "    answer : $RESULT"
  echo "    control: ${control_map[${item}]}"
  if [ -z "${control_map[${item}]}" ]; then
      echo "      ERROR: no control value"
      (( TEST_RESULT++ ))
  else
    i=0
    RESULT_ARRAY=($RESULT)
    for VALUE in ${control_map[${item}]}
    do
      if (( VALUE != RESULT_ARRAY[i] )); then
        echo "      ERROR on value $i: $VALUE != ${RESULT_ARRAY[$i]}"
        (( TEST_RESULT++ ))
      fi
      (( i++ ))
    done
  fi
}

# ==================================================================
# == Start processing

echo "
-----------------------------------------
 Data information 
-----------------------------------------
slave : $NODE_NAME
dir   : $DATA_DIR"

cd $DATA_DIR

# dump the full list of '*.fastq.gz' files available 
# in working directory
echo "> SRA FastQ files available in ${DATA_DIR}:"
i=0
j=0
while read line
do
    FILESIZE=$(stat -c%s "${line}")
    FILESIZE=$(( FILESIZE/(1024*1024) ))  
    echo "  ${line}: $FILESIZE Mb"
    fname=`echo $line | cut -d'.' -f 1`
    ref_file_size_map[$fname]=$FILESIZE
    if (( FILESIZE <= MAX_FILE_SIZE )); then
      array[ $i ]="${fname}"
      (( i++ ))
    fi
    allfiles[ $j ]="$fname"        
    (( j++ ))
done < <(ls -Sr -1 *.fastq.gz)

# Filter file by size retaining only those < MAX_FILE_SIZE
echo "> Total SRA files: ${#allfiles[*]}"
echo "> SRA files to handle now: ${#array[*]} (size <= $MAX_FILE_SIZE Mb)"
for item in ${array[*]}
do
  srafile="${item}.fastq.gz"
  echo "  ${srafile}: ${ref_file_size_map[$item]} Mb"
done

echo "-----------------------"
echo "Running content test..."
echo "-----------------------"
date

# Check file content with GATB Bank API
echo "> SRA file content (using $BANK_SNIPPET_3 snippet):"
for item in ${array[*]}
do
  srafile="${item}.fastq.gz"
  echo "  ${srafile}: ${ref_file_size_map[$item]} Mb"
  checkDataFile $item $srafile
done

if (( TEST_RESULT != 0 )); then
  echo "FAILURE: check $TEST_RESULT error(s), above"
  exit 1
fi

echo "----------------------"
echo "Running compression..."
echo "----------------------"
date
for item in ${array[*]}
do
  srafile="${item}.fastq.gz"
  echo "> compress ${srafile}: ${ref_file_size_map[$item]} Mb"
  echo "  leon -file ${srafile} -c -lossless -nb-cores $CORES -verbose $VERBOSE -kmer-size $KMER_SIZE"
  time $BIN_DIR/leon -file ${srafile} -c -lossless -nb-cores $CORES -verbose $VERBOSE -kmer-size $KMER_SIZE
  leonfile="${item}.fastq.leon"
  FILESIZE=$(stat -c%s "${leonfile}")
  FILESIZE=$(( FILESIZE/(1024*1024) ))  
  echo "  ${leonfile}: $FILESIZE Mb"
done

echo "-----------------------------------"
echo "Running compression content test..."
echo "-----------------------------------"
date
TEST_RESULT=0
echo "> LEON file content (using $BANK_SNIPPET_3 snippet):"
for item in ${array[*]}
do
  leonfile="${item}.fastq.leon"
  echo "  ${leonfile}:"
  checkDataFile $item $leonfile
done

if (( TEST_RESULT != 0 )); then
  echo "FAILURE: check $TEST_RESULT error(s), above"
  exit 1
fi

echo "------------------------"
echo "Running decompression..."
echo "------------------------"
date
for item in ${array[*]}
do
  leonfile="${item}.fastq.leon"
  echo "> decompress ${leonfile}"
  echo "  leon -file ${leonfile} -d -nb-cores $CORES -verbose $VERBOSE"
  time $BIN_DIR/leon -file ${leonfile} -d -nb-cores $CORES -verbose $VERBOSE
done

echo "------------------------------------"
echo "Running decompressed content test..."
echo "------------------------------------"
date
TEST_RESULT=0
echo "> FastQ/Leon file content (using $BANK_SNIPPET_3 snippet):"
for item in ${array[*]}
do
  leonfile="${item}.fastq.d"
  echo "  ${leonfile}:"
  checkDataFile $item $leonfile
done

if (( TEST_RESULT != 0 )); then
  echo "FAILURE: check $TEST_RESULT error(s), above"
  exit 1
fi

echo "----------------------------------------"
echo "Running comparison '.gz' vs. '.leon'..."
echo "----------------------------------------"
date
for item in ${array[*]}
do
  srafile="${item}.fastq.gz"
  leonfile="${item}.fastq.leon"
  echo "> compare ${srafile} vs ${leonfile} ..."
  echo "  $BANK_SNIPPET_1 ${srafile} ${leonfile}"
  time $BIN_DIR/$BANK_SNIPPET_1 ${srafile} ${leonfile}
done

echo "------------------------------------"
echo "Running comparison '.gz' vs. '.d'..."
echo "------------------------------------"
date
for item in ${array[*]}
do
  srafile="${item}.fastq.gz"
  leonfile="${item}.fastq.d"
  echo "> compare ${srafile} vs ${leonfile} ..."
  echo "  $BANK_SNIPPET_1 ${srafile} ${leonfile}"
  time $BIN_DIR/$BANK_SNIPPET_1 ${srafile} ${leonfile}
done

echo "-----------"
echo "Cleaning..."
echo "-----------"
for item in ${array[*]}
do
  rm -f ${item}.fastq.d
  rm -f ${item}.fastq.h5
  rm -f ${item}.fastq.leon
done
rm -rf trashme*

