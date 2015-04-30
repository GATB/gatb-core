#! /bin/bash

set -e

if [ -z "$1" ]; 
then
    export set TEST_USER=$USER
else
    export set TEST_USER=$1
fi  

# we clone the git repository
git clone git+ssh://$TEST_USER@scm.gforge.inria.fr//gitroot/gatb-core/gatb-core.git

# we go into the directory
cd gatb-core/gatb-core

# we build the project
mkdir build; cd build; cmake ..; make 

# we download some banks
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR039/ERR039477/ERR039477.fastq.gz"
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR387/SRR387476/SRR387476.fastq.gz"
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/000/SRR1785130/SRR1785130.fastq.gz"

# we define a function that runs dbgh5 with different sets of parameters, and with a check file
launch () {
    echo "--------------------------------------------------------- $1  $3  ---------------------------------------------------------"
    bin/dbgh5  -kmer-size $3  -in $1 -check $2
    bin/dbgh5  -kmer-size $3  -in $1 -check $2 -bloom basic     -debloom original   -debloom-impl basic
    bin/dbgh5  -kmer-size $3  -in $1 -check $2 -bloom cache     -debloom original   -debloom-impl basic
    bin/dbgh5  -kmer-size $3  -in $1 -check $2 -bloom neighbor  -debloom original   -debloom-impl basic
    bin/dbgh5  -kmer-size $3  -in $1 -check $2 -bloom basic     -debloom cascading  -debloom-impl basic
    bin/dbgh5  -kmer-size $3  -in $1 -check $2 -bloom cache     -debloom cascading  -debloom-impl basic
    bin/dbgh5  -kmer-size $3  -in $1 -check $2 -bloom neighbor  -debloom cascading  -debloom-impl basic
    bin/dbgh5  -kmer-size $3  -in $1 -check $2 -bloom neighbor  -debloom original   -debloom-impl minimizer
    bin/dbgh5  -kmer-size $3  -in $1 -check $2 -bloom neighbor  -debloom cascading  -debloom-impl minimizer
}

# we get the check directory  
export set TEST_CHECK="../test/functional/test1/check"


################################################################################
# we run the tests
################################################################################

launch ./ERR039477.fastq.gz $TEST_CHECK/check/k31/ERR039477.props  31

launch ./ERR039477.fastq.gz $TEST_CHECK/check/k63/ERR039477.props  63
launch ./ERR039477.fastq.gz $TEST_CHECK/check/k95/ERR039477.props  95

launch ./SRR387476.fastq.gz $TEST_CHECK/check/k31/SRR387476.props  31
launch ./SRR387476.fastq.gz $TEST_CHECK/check/k63/SRR387476.props  63
launch ./SRR387476.fastq.gz $TEST_CHECK/check/k95/SRR387476.props  95

launch ./SRR387476.album5/album.txt  $TEST_CHECK/check/k31/SRR387476.props  31
launch ./SRR387476.album5/album.txt  $TEST_CHECK/check/k63/SRR387476.props  63
launch ./SRR387476.album5/album.txt  $TEST_CHECK/check/k95/SRR387476.props  95

launch ./SRR387476.album246/album.txt  $TEST_CHECK/check/k31/SRR387476.props  31
launch ./SRR387476.album246/album.txt  $TEST_CHECK/check/k63/SRR387476.props  63
launch ./SRR387476.album246/album.txt  $TEST_CHECK/check/k95/SRR387476.props  95
launch ./SRR387476.album491/album.txt  $TEST_CHECK/check/k31/SRR387476.props  31
launch ./SRR387476.album491/album.txt  $TEST_CHECK/check/k63/SRR387476.props  63
launch ./SRR387476.album491/album.txt  $TEST_CHECK/check/k95/SRR387476.props  95

launch ./SRR387476.album1227/album.txt  $TEST_CHECK/check/k31/SRR387476.props  31
launch ./SRR387476.album1227/album.txt  $TEST_CHECK/check/k63/SRR387476.props  63
launch ./SRR387476.album1227/album.txt  $TEST_CHECK/check/k95/SRR387476.props  95

launch ./SRR387476.album4904/album.txt  $TEST_CHECK/check/k31/SRR387476.props  31
launch ./SRR387476.album4904/album.txt  $TEST_CHECK/check/k63/SRR387476.props  63
launch ./SRR387476.album4904/album.txt  $TEST_CHECK/check/k95/SRR387476.props  95

launch ./SRR1785130.fastq.gz  $TEST_CHECK/check/k31/SRR1785130.props   31
launch ./SRR1785130.fastq.gz  $TEST_CHECK/check/k63/SRR1785130.props   63
launch ./SRR1785130.fastq.gz  $TEST_CHECK/check/k95/SRR1785130.props   95
launch ./SRR1785130.fastq.gz  $TEST_CHECK/check/k127/SRR1785130.props  127
