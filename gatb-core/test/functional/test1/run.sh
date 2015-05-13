#! /bin/bash

################################################################################
# ARG 1 : login on gforge INRIA
# ARG 2 : email address to send a report 
# ARG 3 : dbgh5 binary; if not set, we do a git clone and build the binary 
################################################################################

set -e

if [ -z "$1" ]; 
then
    export set TEST_USER=$USER
else
    export set TEST_USER=$1
fi  

if [ -z "$2" ]; 
then
    export set TEST_MAIL=$TEST_USER@irisa.fr
else
    export set TEST_MAIL=$2
fi  

if [ -z "$3" ]; 
then
    export set TEST_BINARY=bin/dbgh5
    export set TEST_GIT_CLONE=1
else
    export set TEST_BINARY=$3
    export set TEST_GIT_CLONE=0
fi  

# we get the check directory
TEST_CHECK_TMP=`readlink -f $0`
export set TEST_CHECK=`dirname $TEST_CHECK_TMP`/check
echo $TEST_CHECK

################################################################################
# we dump the configuration
################################################################################

echo "----------------------------------------------------------------------"
echo "                  NON REGRESSION TEST FOR dbgh5                       "
echo "----------------------------------------------------------------------"

echo "TEST_USER      : " $TEST_USER
echo "TEST_MAIL      : " $TEST_MAIL
echo "TEST_BINARY    : " $TEST_BINARY
echo "TEST_GIT_CLONE : " $TEST_GIT_CLONE
echo "TEST_CHECK     : " $TEST_CHECK
echo ""


################################################################################
# we may have to clone the git directory if no dbgh5 is provided
################################################################################

if [ $TEST_GIT_CLONE -eq 1 ]; then
    ################################################################################
    # we clone the git repository
    ################################################################################
    git clone git+ssh://$TEST_USER@scm.gforge.inria.fr//gitroot/gatb-core/gatb-core.git

    ################################################################################
    # we build the project
    ################################################################################
    cd gatb-core/gatb-core

    # we build the project
    mkdir build; cd build; cmake ..; make 
fi


################################################################################
# we download some banks
################################################################################
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR039/ERR039477/ERR039477.fastq.gz"
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR387/SRR387476/SRR387476.fastq.gz"
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/000/SRR1785130/SRR1785130.fastq.gz"


################################################################################
# we define a function that runs dbgh5 with different sets of parameters, and with a check file
################################################################################
launch () {
    echo "--------------------------------------------------------- $1  $3  ---------------------------------------------------------"
    
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2
    
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -bloom basic     -debloom original   -debloom-impl basic
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -bloom cache     -debloom original   -debloom-impl basic
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -bloom neighbor  -debloom original   -debloom-impl basic
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -bloom basic     -debloom cascading  -debloom-impl basic
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -bloom cache     -debloom cascading  -debloom-impl basic
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -bloom neighbor  -debloom cascading  -debloom-impl basic
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -bloom neighbor  -debloom original   -debloom-impl minimizer
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -bloom neighbor  -debloom cascading  -debloom-impl minimizer
    
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -nb-cores 1
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -nb-cores 2
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -nb-cores 4
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -nb-cores 8

    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -max-memory 1000
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -max-memory 1000 -nb-cores 1
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -max-memory 1000 -nb-cores 2
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -max-memory 1000 -nb-cores 4
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -max-memory 1000 -nb-cores 8

    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -max-memory 1000 -max-disk 2000 
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -max-memory 1000 -max-disk 2000 -nb-cores 1
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -max-memory 1000 -max-disk 2000 -nb-cores 2
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -max-memory 1000 -max-disk 2000 -nb-cores 4
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -max-memory 1000 -max-disk 2000 -nb-cores 8

    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -minimizer-size 6
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -minimizer-size 7
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -minimizer-size 8
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -minimizer-size 9
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -minimizer-size 10

    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -minimizer-type 0 
    $TEST_BINARY  -kmer-size $3  -in $1 -check $2 -minimizer-type 1 
}



################################################################################
# we run the tests
################################################################################

launch ./ERR039477.fastq.gz   $TEST_CHECK/k31/ERR039477.props     31
launch ./ERR039477.fastq.gz   $TEST_CHECK/k63/ERR039477.props     63
launch ./ERR039477.fastq.gz   $TEST_CHECK/k95/ERR039477.props     95

launch ./SRR387476.fastq.gz   $TEST_CHECK/k31/SRR387476.props     31
launch ./SRR387476.fastq.gz   $TEST_CHECK/k63/SRR387476.props     63
launch ./SRR387476.fastq.gz   $TEST_CHECK/k95/SRR387476.props     95

launch ./SRR1785130.fastq.gz  $TEST_CHECK/k31/SRR1785130.props    31
launch ./SRR1785130.fastq.gz  $TEST_CHECK/k63/SRR1785130.props    63
launch ./SRR1785130.fastq.gz  $TEST_CHECK/k95/SRR1785130.props    95
launch ./SRR1785130.fastq.gz  $TEST_CHECK/k127/SRR1785130.props  127

launch ./ERR039477.fastq.gz,./SRR387476.fastq.gz   $TEST_CHECK/k31/ERR039477_SRR387476.props     31
launch ./ERR039477.fastq.gz,./SRR387476.fastq.gz   $TEST_CHECK/k63/ERR039477_SRR387476.props     63
launch ./ERR039477.fastq.gz,./SRR387476.fastq.gz   $TEST_CHECK/k95/ERR039477_SRR387476.props     95

################################################################################
# clean up
################################################################################

if [ $TEST_GIT_CLONE -eq 1 ]; then
    cd ../../..
    \rm -rf gatb-core
fi

if [ $? -eq 0 ]; then

   echo "TEST OK"

   echo $0 | mail -s "[gatb] TEST OK" $TEST_MAIL

else
   echo "TEST KO"
   echo $0 | mail -s "[gatb] TEST KO" $TEST_MAIL
fi