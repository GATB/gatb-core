#!/bin/bash

# This script relies on the use of NCBI SRA Tools to fetch FastQ files from 
# SRA entry IDs.
#
# To install SRA Tools on Genocluster, simply use conda as follows:
#   . /local/env/envconda.sh
#   conda config --add channels bioconda
#   conda create -p ~/sra sra-tools=2.8.1
#
# The above procedure has to de bone one times.
# Then, to activate the conda 'sra' tools, simply use:
#   source activate ~/sra
#
# Patrick G. Durand, May 2017 
#

# === CONDA environment: direct use on Genocluster ===
#. /local/env/envconda.sh
#source activate ~/sra
# === CONDA environment: use from GoDocker/Genocluster ===
. /softs/local/env/envconda.sh
source activate $GODOCKER_HOME/sra

# === DATA directory ===
cd /omaha-beach/pdurand/sra-for-leon

SRA_TOOL_PATH=`which fastq-dump`
if [ -z "$SRA_TOOL_PATH" ] ; then
  echo "NCBI SRA Tools not found. Please activate Conda sratools..."
  echo " use: source activate ~/sra"
  exit 1
fi

# This list taken from http://csse.szu.edu.cn/staff/zhuzx/LWFQZip2/SupplementaryFile.pdf
# Table S10.
# array_ok: SRA files on which Leon is ok
# array_nok: SRA files on which Leon is not ok 
#   (according to article authors: lose fidelity after decompression)
array_ok=(SRR2916693 SRR2994368 SRR3190692 ERR385912 SRR034509 ERR174310)
# ok files should include: ERR194147, but failed to download from NCBI after several trials!
array_nok=(SRR3211986 ERR739513 ERR386131)

# data set from Leon publication
array_leon=(SRR065390 SRR959239 SRR857303 SRR1870605 SRR445718 SRR1519083)

array=("${array_ok[@]}" "${array_nok[@]}" "${array_leon[@]}")

echo "Nb. SRA files to download: ${#array[*]}"

echo "SRA are:"
for item in ${array[*]}
do
  srafile="${item}.fastq.gz"
  if [ -e $srafile ]; then
    echo "   ${srafile}: exists"
  else
    echo "   ${srafile}: will be downloaded"
  fi  
done

for item in ${array[*]}
do
  srafile="${item}.fastq.gz"
  if [ ! -e $srafile ]; then
    echo "> downloading ${item} to ${srafile} ..."
    echo "  time fastq-dump --gzip $item"
    time fastq-dump --gzip $item
  fi
done


