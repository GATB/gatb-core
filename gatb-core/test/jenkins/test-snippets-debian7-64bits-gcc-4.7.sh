#!/bin/bash
#--------------------------------------------------------------#
#         Continuous integration script for Jenkins            #
#--------------------------------------------------------------#

set +xv

echo "
-----------------------------------------
 Miscellaneous information·
-----------------------------------------
date      : `date`
hostname  : `hostname`
pwd       : `pwd`

-----------------------------------------
 Jenkins build parameters (user defined)
-----------------------------------------
BRANCH_TO_BUILD      : ${BRANCH_TO_BUILD}
RELEASE_TO_BUILD     : ${RELEASE_TO_BUILD}

-----------------------------------------
 Jenkins build parameters (built in)
-----------------------------------------
BUILD_NUMBER         : ${BUILD_NUMBER}

"

set -xv

# quick look at resources
#---------------------------------------------------------------
free -h
#---------------------------------------------------------------
lstopo
#---------------------------------------------------------------
df -kh
#---------------------------------------------------------------


################################################################
#                       COMPILATION                            #
################################################################

gcc --version
g++ --version
cmake --version

[ `gcc -dumpversion` = 4.7 ] && { echo "GCC 4.7"; } || { echo "GCC version is not 4.7, we exit"; exit 1; }

JENKINS_TASK=test-snippets-debian7-64bits-gcc-4.7
GIT_DIR=/builds/workspace/$JENKINS_TASK/gatb-core
BUILD_DIR=/scratchdir/$JENKINS_TASK/gatb-core/build

#>>>>>>>>>>>>>>>>>>>>>
#if false; then 
#>>>>>>>>>>>>>>>>>>>>>

rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR

cd $BUILD_DIR

#---------------------------------------------------------------
cmake -DGATB_CORE_INCLUDE_EXAMPLES=True -Wno-dev $GIT_DIR

#---------------------------------------------------------------
make -j2 examples

#>>>>>>>>>>>>>>>>>>>>>
#fi
#>>>>>>>>>>>>>>>>>>>>>

################################################################
#                      RUN SNIPPETS                            #
################################################################

cd $BUILD_DIR/bin/

# list all compiled snippets
ls -l

# prepare environment
DB_DIR=${GIT_DIR}/test/db
TMP_DIR=${BUILD_DIR}/tmp
[ ! -d ${TMP_DIR} ] && { mkdir ${TMP_DIR}; }

# Note:
#   we redirect stderr and stdout to dev/null to avoir poluting log
#   with lots of data. In case of an error on a snippet:
#   1. un-comment above 'if-then' to avoid compiling code
#   2. remove stream redirection on the snippet of interest

###################### iterator snippets ###########################
./iterators1 > /dev/null 2>&1 # no args expected
./iterators2 > /dev/null 2>&1 # no args expected
./iterators3 > /dev/null 2>&1 # no args expected
./iterators4 > /dev/null 2>&1 # no args expected
./iterators5 > /dev/null 2>&1 # no args expected
./iterators6 > /dev/null 2>&1 # no args expected
./iterators7 > /dev/null 2>&1 # no args expected
./iterators8 > /dev/null 2>&1 # no args expected
./iterators9 > /dev/null 2>&1 # no args expected

###################### bank snippets ###########################
./bank1 ${DB_DIR}/reads1.fa > /dev/null 2>&1
./bank2 ${DB_DIR}/reads1.fa > /dev/null 2>&1
./bank3 ${DB_DIR}/reads1.fa > /dev/null 2>&1
./bank4 ${DB_DIR}/reads1.fa > /dev/null 2>&1
./bank5 ${DB_DIR}/reads1.fa > /dev/null 2>&1
./bank6 ${DB_DIR}/reads1.fa > /dev/null 2>&1
./bank7 ${DB_DIR}/reads1.fa 500 > /dev/null 2>&1
./bank8 ${DB_DIR}/reads1.fa ${TMP_DIR}/reads1.fa.bin > /dev/null 2>&1
./bank9 ${DB_DIR}/reads1.fa ${TMP_DIR}/reads1_filtered.fa 5 > /dev/null 2>&1
./bank10 -in ${DB_DIR}/reads1.fa -out-dir ${TMP_DIR} -max-size 5000 > /dev/null 2>&1
./bank11 -kmer-size 4 -out ${TMP_DIR}/kmer.fa > /dev/null 2>&1
./bank12 -in ${DB_DIR}/reads1.fa -seq-ids ${DB_DIR}/reads1_list_ids.txt
./bank13 -in ${DB_DIR}/reads1.fa -filter-ratio 0.8 > /dev/null 2>&1
./bank14 -in ${DB_DIR}/reads1.fa > /dev/null 2>&1
./bank15 ${DB_DIR}/reads1.fa > /dev/null 2>&1
./bank16 -in ${DB_DIR}/reads1.fa > /dev/null 2>&1
./bank17 -in ${DB_DIR}/album.txt  > /dev/null 2>&1
./bank18 -one ${DB_DIR}/reads1.fa -two ${DB_DIR}/reads2.fa > /dev/null 2>&1
./bank19 -in ${DB_DIR}/reads1.fa > /dev/null 2>&1
./bank20 -in ${DB_DIR}/reads1.fa > /dev/null 2>&1

###################### kmer snippets ###########################
./kmer1 > /dev/null 2>&1 # no args expected
./kmer2 > /dev/null 2>&1 # no args expected
./kmer3 > /dev/null 2>&1 # no args expected
./kmer4 -in ${DB_DIR}/reads1.fa -kmer-size 11 > /dev/null 2>&1
./kmer5 -in ${DB_DIR}/reads1.fa -kmer-size 11 -minimizer-size 11 > /dev/null 2>&1
# kmer6 will report an EXCEPTION: this is the expected behavior
./kmer6 > /dev/null 2>&1 # no args expected
./kmer7 -in ${DB_DIR}/reads1.fa -kmer-size 11 > /dev/null 2>&1
./kmer8 -in ${DB_DIR}/reads1.fa -kmer-size 11 -minimizer-size 11 > /dev/null 2>&1
./kmer9 -in ${DB_DIR}/reads1.fa > /dev/null 2>&1
./kmer10 -in ${DB_DIR}/reads1.fa > /dev/null 2>&1
./kmer11 -in ${DB_DIR}/reads1.fa > /dev/null 2>&1
./kmer12 -in "${DB_DIR}/reads1.fa,${DB_DIR}/reads2.fa" > /dev/null 2>&1
./kmer13 -in "${DB_DIR}/reads1.fa,${DB_DIR}/reads2.fa" > /dev/null 2>&1
./kmer14 -in ${DB_DIR}/reads1.fa > /dev/null 2>&1
./kmer15 -in ${DB_DIR}/reads1.fa > /dev/null 2>&1
./kmer16 -in ${DB_DIR}/reads1.fa > /dev/null 2>&1
./kmer17 -in ${DB_DIR}/reads1.fa -kmer-size 11 -minimizer-size 11 > /dev/null 2>&1

###################### debruijn snippets #######################
./debruijn1 -in ${DB_DIR}/reads1.fa > /dev/null 2>&1
./debruijn2 ${DB_DIR}/reads1.fa > /dev/null 2>&1
./debruijn3 ${DB_DIR}/reads1.fa > /dev/null 2>&1
./debruijn4 > /dev/null 2>&1 # no args expected
./debruijn5 ${DB_DIR}/celegans_reads.h5 > /dev/null 2>&1
./debruijn6 ${DB_DIR}/celegans_reads.h5 > /dev/null 2>&1
./debruijn7 ${DB_DIR}/celegans_reads.h5 > /dev/null 2>&1
./debruijn8 ${DB_DIR}/celegans_reads.h5 > /dev/null 2>&1
./debruijn9 > /dev/null 2>&1  # no args expected
./debruijn10 > /dev/null 2>&1 # no args expected
./debruijn11 > /dev/null 2>&1 # no args expected
./debruijn12 > /dev/null 2>&1 # no args expected
./debruijn13 > /dev/null 2>&1 # no args expected
./debruijn14 > /dev/null 2>&1 # no args expected
./debruijn15 > /dev/null 2>&1 # no args expected
./debruijn16 > /dev/null 2>&1 # no args expected
./debruijn17 > /dev/null 2>&1 # no args expected
./debruijn18 -graph ${DB_DIR}/celegans_reads.h5 > /dev/null 2>&1
./debruijn19 -graph ${DB_DIR}/celegans_reads.h5 > /dev/null 2>&1
./debruijn20 -graph ${DB_DIR}/celegans_reads.h5 > /dev/null 2>&1
./debruijn21 -graph ${DB_DIR}/celegans_reads.h5 > /dev/null 2>&1
./debruijn22 -graph ${DB_DIR}/celegans_reads.h5 -out ${TMP_DIR}/data.h5 > /dev/null 2>&1
./debruijn23 > /dev/null 2>&1 # no args expected
./debruijn24 -graph ${DB_DIR}/celegans_reads.h5 > /dev/null 2>&1
./debruijn25 -in ${DB_DIR}/reads1.fa > /dev/null 2>&1
./debruijn26 > /dev/null 2>&1 # no args expected
./debruijn27 ${DB_DIR}/celegans_reads.h5 "TTTGCCCATTTCCTGCCATTTGTC" > /dev/null 2>&1 

./traversal1 -graph ${DB_DIR}/celegans_reads.h5 > /dev/null 2>&1
./traversal2 -traversal unitig > /dev/null 2>&1
./traversal2 -traversal contig > /dev/null 2>&1

###################### storage snippets ########################
./storage1 > /dev/null 2>&1 # no args expected
./storage2 > /dev/null 2>&1 # no args expected
./storage3 > /dev/null 2>&1 # no args expected
./storage4 > /dev/null 2>&1 # no args expected
./storage5 > /dev/null 2>&1 # no args expected
./storage6 -graph ${DB_DIR}/celegans_reads.h5 > /dev/null 2>&1
./storage7 > /dev/null 2>&1 # no args expected
./storage8 > /dev/null 2>&1 # no args expected
./storage9 ${DB_DIR}/celegans_reads.h5 > /dev/null 2>&1

###################### tools snippets ##########################
./ToyTool > /dev/null 2>&1 # no args expected
./multithreading1 > /dev/null 2>&1 # no args expected
./multithreading2 > /dev/null 2>&1 # no args expected
./multithreading3 > /dev/null 2>&1 # no args expected
./multithreading4 > /dev/null 2>&1 # no args expected
./multithreading5 > /dev/null 2>&1 # no args expected
./multithreading6 ${DB_DIR}/reads1.fa > /dev/null 2>&1
./observer1 > /dev/null 2>&1 # no args expected
./optionsparser1 > /dev/null 2>&1 # no args expected
./optionsparser2 > /dev/null 2>&1 # no args expected

###################### protos snippets #########################
./MicroSNP -in ${DB_DIR}/microsnp.fa -kmer-size 7 -abundance-min 1 > /dev/null 2>&1

