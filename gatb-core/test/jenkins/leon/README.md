| **Functional tests** |
|----------------------|
| [![Build Status](https://ci.inria.fr/gatb-core/view/Leon/job/tool-leon-functional-tests/badge/icon)](https://ci.inria.fr/gatb-core/view/Leon/job/tool-leon-functional-tests/) |

# Introduction

This directory contains reference material used to test Leon compressor
by means of Jenkins tasks.

It is intended to be used by Genscale team.

# Scripts

Tests are actually run by 'tool-leon-functional-tests' Jenkins task
from GATB project on INRIA CI platform; please refer to:
[https://ci.inria.fr/gatb-core/view/Leon/job/tool-leon-functional-tests/](https://ci.inria.fr/gatb-core/view/Leon/job/tool-leon-functional-tests/)

this task runs in the row the following scripts:

* first: ```tool-leon-functional-tests.sh``` to compile GATB-Core binaries
* then: ```tool-leon-functional-tests-test.sh``` to run Leon on 15 SRA files
    totalizing more than 50 Gb of gzipped reads data.

# Retrieving data files

NCBI SRA files have to be retrieved from Genocluster computing nodes 
using the following script: ```download.sh```. That script can be used from 
a terminal session or, even better, using GoDocker.
 
Then, files integrity can be checked using ```test_integrity.sh``` script.

# Data files

Files used to test Leon are as follows:

```Leon original publication, G. Benoit et al., 2015```:
[https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0709-7](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0709-7)
cf. Table ST1, supplementary file:

* SRR065390: C. elegans WGS Illumina
* SRR959239: E. coli WGS Illumina
* SRR1519083: Metagenome Illumina
* SRR1870605: E. coli WGS Illumina Miseq
* SRR445718: Human RNA-seq Illumina
* SRR857303: E. coli WGS Ion Torrent

```Evaluation of Leon by Y. Zhang et al, 2017```:
[http://csse.szu.edu.cn/staff/zhuzx/LWFQZip2/SupplementaryFile.pdf](http://csse.szu.edu.cn/staff/zhuzx/LWFQZip2/SupplementaryFile.pdf)
Cf. table S10:

* SRR2916693: 454 GS
* SRR2994368: Illumina Miseq
* SRR3211986: Pacbio RS
* ERR739513: MiniION
* SRR3190692: Illumina Miseq
* ERR385912: Illumina Hiseq 2000
* ERR386131: Ion Torrent PGM
* SRR034509: Illumina Analyzer II
* ERR174310: Illumina Hiseq 2000



