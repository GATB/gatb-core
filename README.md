

This is a branch with htslib to read sam/bam/cram files. The compilation process is not yet integrated fully with cmake. 
To compile this branch, follow the steps :


### download project

	git clone https://github.com/GATB/gatb-core.git
	cd gatb-core/gatb-core
	git checkout htsfiles


### compile htslib 

	cd thirdparty/htslib-1.9/
	autoheader
	autoconf
	./configure -prefix=$PWD/build
	make
	make install

### compile gatb core 

	cd ../../
	mkdir build; cd build; cmake  -DGATB_CORE_INCLUDE_EXAMPLES=True ..; make -j8

### test an example with a bam file

	./bin/bank2 ../test/db/ex1.bam

## License

GATB is free software; you can redistribute it and/or modify it under the [Affero GPL v3 license](http://www.gnu.org/licenses/agpl-3.0.en.html).
