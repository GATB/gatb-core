# GATB-Core example code snippets

This directory contains many ready-to-compile code snippets you can use to learn [GATB-Core c++ API](http://gatb-core.gforge.inria.fr/doc/api/).

## Dependencies

The following third parties have to be already installed to compile GATB-Core examples:

* a **C++/11 capable compiler** (*e.g.* gcc 4.7+, clang 3.5+, Apple/clang 6.0+)
* **CMake 3.1+**


## Compile GATB-CORE snippets

### Get a copy of the project

If not already done:

	cd <some_directory>
	git clone https://github.com/GATB/gatb-core.git

### Compile all snippets

	cd <some_directory>/gatb-core/gatb-core
	mkdir build
	cd build
	cmake -DGATB_CORE_INCLUDE_EXAMPLES=True -DCMAKE_BUILD_TYPE=Debug .. 
	make -j8 examples

### Compile a particular snippet

Instead of using:

	make -j8 examples

simply do:

	make -j8 bank15

to compile snippet "bank15.cpp". Apply the same recipe to compile any other code snippets.

### Simple Makefile

Alternatively, you can use the provided "makefile" script to compile a single example. It needs to be modified to point to the correct path of the GATB-core library. Maybe it will need some tweaking (try removing the "-static" flag if compilation fails).

Try it, from this folder:

   make bank/bank1

should compile the first bank1.cpp example

### Run a compiled code snippet

Have a look at the begining of each c++ source code: you'll see how to use the example programs.

For instance, taking the above example "bank15.cpp", you run it as follows:

    # from within the 'build' directory:
    ./bin/bank15 ../test/db/reads1.fa

## Documentation

Basic APIs explained:

* [Bank](http://gatb-core.gforge.inria.fr/doc/api/snippets_bank.html): read/write FastA and FastQ files (plain text and gzip)
* [Iterator](http://gatb-core.gforge.inria.fr/doc/api/snippets_iterators.html): go through a Bank by iterating over its sequences
* [k-mer](http://gatb-core.gforge.inria.fr/doc/api/snippets_kmer.html): from sequences to k-mers
* [De Bruijn Graph](http://gatb-core.gforge.inria.fr/doc/api/snippets_graph.html): from sequences to De Bruijn Graphs

Advanced APIs explained:

* [Multi-threading](http://gatb-core.gforge.inria.fr/doc/api/snippets_multithread.html): easy way to manage multi-threaded tasks on Linux and OSX
* [Storage](http://gatb-core.gforge.inria.fr/doc/api/snippets_storage.html): easy way to handle HDF5 storage

Make your own full-featured GATB-Tool:

* [Tool API](http://gatb-core.gforge.inria.fr/doc/api/snippets_tools.html): make a command-line based tool quickly

The complete GATB-Core API reference documentation is [here](http://gatb-core.gforge.inria.fr/doc/api/index.html).

## Online tutorial

Some of these code snippets are also available for direct use from the [GATB-Core online Tutorial](https://gatb.inria.fr/gatb-programming-tutorial/).

## Contact

To contact a developer, request help, *etc*, use: 

    https://gatb.inria.fr/contact/
    

## License

GATB is free software; you can redistribute it and/or modify it under the [Affero GPL v3 license](http://www.gnu.org/licenses/agpl-3.0.en.html).
