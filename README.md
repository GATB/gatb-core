# GATB - The Genome Analysis Toolbox with de-Bruijn graph

&nbsp;&nbsp;&nbsp;[![](https://img.shields.io/badge/release-1.3.0-orange.svg?style=plastic)](https://github.com/GATB/gatb-core/releases)&nbsp;-&nbsp;[![](https://img.shields.io/badge/build--Linux-passing-green.svg?style=plastic)]()&nbsp;&nbsp;[![](https://img.shields.io/badge/build--OSX-passing-green.svg?style=plastic)]()

--------------------------------------------------------------------------------

&nbsp;&nbsp;&nbsp;[![License](http://img.shields.io/:license-Affero--GPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)&nbsp;&nbsp;[![](https://tokei.rs/b1/github/GATB/gatb-core?category=code)](https://github.com/GATB/gatb-core)&nbsp;&nbsp;[![](https://img.shields.io/badge/platform-c++/11-yellow.svg)](https://isocpp.org/wiki/faq/cpp11)&nbsp;&nbsp;[![](https://img.shields.io/badge/run_on-Linux--Mac_OSX-yellowgreen.svg)]()

--------------------------------------------------------------------------------
**Continuous integration on master branch - build status:**

|**Linux**| **gcc 4.7** | **gcc 4.8** | **gcc 4.9** | **clang 3.6** | **clang 3.9**| **Valgrind** |
|---------|-------------|-------------|-------------|---------------|--------------|--------------|
|*Debian 8*|    | [![Build Status](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-docker-gatb-core-compile-gcc48/badge/icon)](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-docker-gatb-core-compile-gcc48/) | [![Build Status](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-docker-gatb-core-compile-gcc49/badge/icon)](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-docker-gatb-core-compile-gcc49/) | [![Build Status](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-docker-gatb-core-compile-clang36/badge/icon)](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-docker-gatb-core-compile-clang36/) | [![Build Status](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-docker-gatb-core-compile-clang39/badge/icon)](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-docker-gatb-core-compile-clang39/) |  - | 
|*Debian 7*| [![Build Status](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-suite-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-suite-debian7-64bits-gcc-4.7/) | - | - | - | - | [![Build Status](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-valgrind-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-valgrind-debian7-64bits-gcc-4.7/) | 
|*Fedora 20*| - | [![Build Status](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-suite-fedora20-gcc-4.8/badge/icon)](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-suite-fedora20-gcc-4.8/) | - | - | - | - | 

| **Mac OSX** | **clang-600** | **gcc 4.2.1** |
|    :--:     |---------------|---------------|
| *10.9* | [![Build Status](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-suite-macos-10.9.5-clang-6.0/badge/icon)](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-suite-macos-10.9.5-clang-6.0/) | [![Build Status](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-suite-macos-10.9.5-gcc-4.2.1/badge/icon)](https://ci.inria.fr/gatb-core/view/GATB-CORE/job/test-suite-macos-10.9.5-gcc-4.2.1/) | 


--------------------------------------------------------------------------------
# What is GATB?

GATB is made of two master projects: 

* The **GATB-CORE project** provides a set of highly efficient algorithms to analyse NGS data sets. These 
methods enable the analysis of data sets of any size on multi-core desktop computers, including very 
huge amount of reads data coming from any kind of organisms such as bacteria, plants, animals and 
even complex samples (e.g. metagenomes). Read more about GATB at <a href="https://gatb.inria.fr/">https://gatb.inria.fr/</a>.
 
* The **GATB-TOOLS project** contains a set of ready-to-use softwares relying on GATB-CORE algorithms. You can have a look at available tools at <a href="https://gatb.inria.fr/software/">https://gatb.inria.fr/software/</a>.

--------------------------------------------------------------------------------
# What is GATB-CORE ?

GATB-CORE is a high-performance and low memory footprint C++ library.

GATB-Core natively provides the following operations:

* Reads handling: 	
 * FASTA/FASTQ parsing
 * Parallel iteration of sequences

* K-mer: 	

 * K-mer counting
 * Minimizer computation of k-mers, partitioning of datasets by minimizers
 * Bloom data structure of k-mers
 * Hash table of k-mers
 * Minimal perfect hash function of k-mers
 * Arbitrarily large k-mers representations

* de Bruijn graph: 	

 * graph construction
 * graph traversal operations (contigs, unitigs)
 * graph simplifications for assembly (tip removal, bulge removal)

By itself GATB-CORE is not an NGS data analysis tool. However, it can be 
used to create such tools; see section [Quickly create a new GATB-TOOL software](#Quickly create a new GATB-TOOL software), below.

They already exist a set of ready-to-use tools relying on GATB-CORE library: see [https://gatb.inria.fr/software/](https://gatb.inria.fr/software/)

# Project content

All the needed material of GATB-CORE is contained in the current directory in order to 
generate the wanted artifacts:  

* dynamic and static libraries holding the services component

* 120+ unit tests of the entire library

# Dependencies

The following third parties have to be already installed to compile GATB-Core:

* a **C++/11 capable compiler** (*e.g.* gcc 4.7+, clang 3.5+, Apple/clang 6.0+)
* **CMake 3.1+**

In addition, you could install these optional tools:

* **doxygen**: to compile a local copy of the GATB-Core documentation
* **cppunit**: to compile and run Unit tests

# Compile GATB-CORE


## Compile in Release mode (default)

Type:

	cd <some_directory>
	git clone https://github.com/GATB/gatb-core.git
	cd gatb-core/gatb-core
	mkdir build ; cd build ; cmake .. ; make -j8
	
## Compile in Debug mode

Type same as above, except for the CMake command:

    cmake -D CMAKE_BUILD_TYPE=Debug ..
    make -j8

## Run unit tests

* cppunit is required
* compile using the command above (Release or Debug mode)

Then type:


     # enter gatb-core build directory
     cd gatb-core/gatb-core/build
     # set verbose mode to on so that we have name of failing tests (if any)
     export CPPUNIT_VERBOSE=1
     # Copy database for unit tests
     cp -r ../test/db ./test/
     # Launch the full test suite
     cd bin
     ./gatb-core-cppunit

The gatb-core-cppunit command may also take as argument the categories of tests that show up in the verbose output, e.g. './gatb-core-cppunit TestBank'.

More about GATB-CORE code compiling instruction is available [here](http://gatb-core.gforge.inria.fr/doc/api/compilation.html).

# Work on GATB-Core code using Eclipse

Read [this documentation](https://gatb.inria.fr/use-eclipse-to-develop-gatb-core-softwares/).

# Work on GATB-Core code using Xcode

Read [this documentation](https://gatb.inria.fr/use-xcode-to-develop-gatb-core-softwares/).


# Learning GATB-Core: tutorial

You can follow [this link](https://gatb.inria.fr/gatb-programming-tutorial/) to start the GATB-Core Online Tutorial trail.

The project also contains many [code examples](https://github.com/GATB/gatb-core/tree/master/gatb-core/examples) that can be easily compiled and executed to review how to use GATB-Core APIs.

# Documentation

The complete GATB-Core documentation is available [here](http://gatb-core.gforge.inria.fr/doc/api/). It contains: API, code snippets, compile instructions, *etc*.

Nevertheless, you can create a local copy of the documentation as follows (we suppose you already compiled the c++ code, see above; requires 'doxygen'):

     cd gatb-core/gatb-core/build
     make doc

Documentation is then available in _build/html/index.html_


# kmer default sizes

By default, the library is compiled for supporting 4 ranges of kmers : 

* k1 : for kmerSize < k1  (default value 32)         
* k2 : for k1 <= kmerSize < k2 (default value 64)
* k3 : for k2 <= kmerSize < k3 (default value 96)
* k4 : for k3 <= kmerSize < k4 (default value 128)

You can customize these values through cmake, provided they rebuild the project from scratch. For instance:

    rm -Rf build; mkdir build ; cd build ; cmake -DKSIZE_LIST="64 96 128 162" ..

Tools may set a default kmer lists in their CMakeFiles.txt, as such (see for instance Minia):

    list (APPEND KSIZE_DEFAULT_LIST  32   64   96  128  160  192  224  256)


--------------------------------------------------------------------------------
# Directory content

* __README__:                  this file

* __CMakeList.txt__:           global cmake file

* __doc__:                 
    * __design__      design documentation for the component
    * __doxygen__     pages for doxygen

* __examples__:       snippets showing how to use the library                 

* __src__:            source code for the component

* __test__:           tests directory
    * __src__         source code for unit tests
    * __db__          FASTA databases for unit tests

* __thirdparty__:    third parties    

--------------------------------------------------------------------------------
# Details for 'src' directory content

It contains several sub directories, each one corresponding to one software package.

A package may be composed of sub packages; the directory hierarchy should represent
this packages tree structure.

For one atomic package (or sub package), we should have:

* directory 'api'       API of the package  

* directory 'impl'      several implementations of the API

--------------------------------------------------------------------------------
# Quickly create a new GATB-TOOL software

A GATB-TOOL is a new software relying upon GATB-CORE.

You use GATB-CORE to create a new tool project, with the following script:

    sh scripts/NewProject/NewProject.sh -d directory -n toolName

where:

    'directory' is the directory where the project will be created
    'toolName' is the name of the project.

The script will automatically creates the full path 'directory/toolName' to deploy a self-contained tool.
 
By default, the following part will be included in the project:

* a CMakeLists.txt file used for building the project
* a 'tools' directory holding a default source code using GATB-Core
* a 'scripts' directory holding a script to automatically package the tool
* an optional 'thirdparty' directory holding the gatb-core resources

The 'thirdparty' directory is only available for tool created outside the GATB-Tools repository.

Tools located within GATB-Tools rely on a common GATB-CORE sub-module already available in this repository.

The directory where the project is created has no link to any external resources. You can therefore
move it anywhere you want.

Such a project can be a start for building applications based on GATB-CORE. 

More on creating a new GATB-Core based project: [http://gatb-core.gforge.inria.fr/doc/api/new_project.html](http://gatb-core.gforge.inria.fr/doc/api/new_project.html)


#Contact

To contact a developer, request help, *etc*, use: 

    https://gatb.inria.fr/contact/
    

# License

GATB is free software; you can redistribute it and/or modify it under the [Affero GPL v3 license](http://www.gnu.org/licenses/agpl-3.0.en.html).