This folder contains data set used by GATB-Core CPPUnit test codes.

## Leon files

* Leon*.fastq: DO NOT modify

* leon1.fastq.leon-ref
  and leon2.fastq.leon-ref created as follows

    leon -c -file ./gatb-core/gatb-core/test/db/leon1.fastq -lossless -verbose 0 -kmer-size 31 -abundance 1

* NIST7035* files from: 

    ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz 

* giab.hg002* files from:

    ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/CORNELL_Oxford_Nanopore/giab.hg002.2D.fastq

* Note: loading files from ftp server can be none as follows:

    curl --user anonymous:YOUR-EMAIL ftp://ftp-trace.../.../NIST7035.fastq.gz -o NIST7035.fastq.gz
