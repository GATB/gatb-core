--------------------------------------------------------------------------------
# Introduction

The current project provides a library providing some services packaged in software components.

All the needed material is contained in the current directory in order to 
generate the wanted artifacts:  

* dynamic and static libraries holding the services component

* unit tests of the component

* wrappers for the components in other languages (java, python)

A common way to build all the artifacts is :

	mkdir build ; cd build ; cmake .. ; make

Then, the documentation is available at _build/html/index.html_

and the unit tests can be launched by the following command line :

	build/test/unit/GatbToolsTest out.xml
	
If you don't specify an output file name, you will get no console output, but you will be able to use the exit status code of the command (useful for automatization processes).

--------------------------------------------------------------------------------
# Dependencies

The following third parties should be already installed:

* doxygen
* cppunit

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

* __external__:    third parties    

--------------------------------------------------------------------------------
# Details for 'src' directory content

It contains several sub directories, each one corresponding to one software package.

A package may be composed of sub packages; the directory hierarchy should represent
this packages tree structure.

For one atomic package (or sub package), we should have:

* directory 'api'       API of the package  

* directory 'impl'      several implementations of the API