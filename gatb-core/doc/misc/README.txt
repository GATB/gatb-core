                          ***********************************
                          *                                 *
                          *  GATB-Core binary distribution  *
                          *                                 *
                          ***********************************

==================
== INTRODUCTION ==
==================

This release provides you with a pre-compiled binary release of the GATB-Core library. 
Using it, you can easily compile home-made GATB software, as explained below.


=============
== CONTENT ==
=============

    - README.txt : the file you are reading
    - lib        : static GATB-Core library
    - include    : header files; usually, you’ll use ‘gatb/gatb_core.hpp’
    - bin        : De Bruijn graph utilities (see below)
    - test       : unit tests of the library
    - examples   : a few snippets showing how to use the library


===================
== DOCUMENTATION ==
===================

For documentation, please have a look at:

- tutorial with code snippets: https://gatb.inria.fr/gatb-programming-tutorial/
  
  This tutorial runs in action some of the code snippets available in the “examples” directory.
  
  All these snippets are also described here: 
  http://gatb-core.gforge.inria.fr/doc/api/snippets_page.html

- API documentation: http://gatb-core.gforge.inria.fr/doc/api/
  
  The complete c++ API of GATB-Core.


=====================
== COMPILING CODES ==
=====================

You can compile a code based on GATB-Core using the library as follows. Here, we simply compile a snippet 
taken from the 'examples' directory:

- on Linux:
   g++ examples/debruijn/debruijn1.cpp -Iinclude -Llib -lgatbcore -lhdf5 -ldl -lz -lpthread  -std=c++0x -O3 -o debruijn1
    
-on MacOS X:
    clang++ examples/debruijn/debruijn1.cpp -Iinclude -Llib -lgatbcore -lhdf5 -ldl -lz -lpthread  -std=c++0x -O3 -DBOOST_NO_CXX11_RVALUE_REFERENCES=1 -o debruijn1

To get started with your own code (let's say, test.cpp), here is a quick walkthrough

    wget https://github.com/GATB/gatb-core/releases/download/v1.2.2/gatb-core-1.2.2-bin-Linux.tar.gz
    tar xf gatb-core-1.2.2-bin-Linux.tar.gz
    mv gatb-core-1.2.2-bin-Linux gatb-core

Create a small example code in test.cpp:

    #include <gatb/gatb_core.hpp>
    int main (int argc, char* argv[])
    {
        // a small GATB example

        const size_t span = KMER_SPAN(1);
        typedef Kmer<span>::ModelCanonical Model;
        Model model (5);
        Model::Kmer kmer = model.codeSeed ("AAGTC", Data::ASCII);
        std::cout << "revcomp kmer: " << model.toString(kmer.revcomp())    << std::endl;
    }


Now the folder structure looks like:
  
    test.cpp
    gatb-core/
    gatb-core/include/
    gatb-core/lib/
    ...

To compile:
    
    g++ test.cpp -Igatb-core/include -Lgatb-core/lib -lgatbcore -lhdf5 -ldl -lz -lpthread  -std=c++0x -O3 -o test 

Then the program is ready to run:
    
    ./test
    
Output:
    
    revcomp kmer: GACTT

========================================
== TESTING THE LIBRARY ON YOUR SYSTEM ==
========================================

You can check the library by launching unit tests on your system, as follows:

  cd test
  export CPPUNIT_VERBOSE=1
  ./gatb-core-cppunit


=========================
== GATB-Core Utilities ==
=========================

The ‘bin’ directory contains the following tools:

- dbgh5: converts a sequence file into a De Bruin graph and stores it into a HDF5 formatted file; 
         file extension is: .h5
         (HDF5 is the data format used by GATB-Core to store compact representation of a graph)

- dbginfo: displays some information about a ‘.h5’ file created with ‘dbgh5’ tool

- gatb-h5dump: dump the content of a ‘.h5’ file created with ‘dbgh5’ tool

