#About XXX

This is a software tool relying on GATB-CORE library.

The architecture of the tool is as follows:

    * a CMakeLists.txt file used for building the project
    * a 'tools' directory holding a default source code using GATB-Core
    * a 'scripts' directory holding a script to automatically package the tool
    * a 'thirdparty' directory holding the gatb-core resources
    * a 'doc' directory
    * a 'tests' directory holding test procedures
    
The 'thirdparty' directory is only available for tool created outside the GATB-Tools repository.
Tools located within GATB-Tools rely on a common GATB-Core sub-module already available in this repository.

It is advised to use:

    * 'tests' directory to hold test procedures: scripts and/or small sized data files
    * 'scripts' directory to hold any scripts this tool relies on
    * 'doc' directory to hold tool's documentation
    * 'thirdparty' directory to hold any third-party librairies this tool relies on
    
It is worth noting that 'tools' directory is organised using sub-folders; by default, there is
at least one such sub-folder called 'XXX'. It holds the source code of the tool. However, when
considering a more complex software, it could be nice to setup several "sub-tools", each of them
making a particular feature. In that case, you can easily create several "sub-tool" folders inside
"tools", each of them having a "src" folder containing the source code, as well as a "main.cpp", for
each feature. Using this organisation has a big advantage: the provided CMakeLists.txt is aware of 
that, so you do not have to edit the CMake file when you add a new "sub-tool". As a real example, you
can have a look at the DiscoSNP software.

#License

Please not that GATB-Core is distributed under Affero-GPL license.

#Dependencies

The following third parties should be already installed:

* cmake 2.8+ (mandatory)

#Project build

For building your project, you should do the following
    
    cd [path-to-tool-home]
    mkdir build;  cd build;  cmake ..;  make -j8
    
Then, you should get a binary holding the name of the project within 'build/tools'.

Note: the first compilation should take some time since the GATB-CORE library is generated.

#Project packaging

You can prepare your tool for distribution using:
    
    ./[path-to-tool-home]/scripts/package_tool.sh -M X -m Y -p Z

With X, Y and Z being major, minor and patch release numbers, respectively.

Then, you should get two 'tar.gz' files within 'build', one containing the binary release 
and the other the source codes.

Note: the script re-builds the entire tool from its sources to ensure a clean build process.

#Examples

The project is created with a default 'main' function that dumps some information about the library.

You can find many snippets showing how to use the library. 
These snippets are located in '[path-to-gatb-core]/examples' and are split in several fields.

You can find documentation about the snippets <a href="http://gatb-core.gforge.inria.fr/doc/api/snippets_page.html">here</a>.

You can copy the content of one of the snippet file into the 'tools/XXX/src/main.cpp' file and relaunch the build.

For instance:

    cp [path-to-gatb-core]/examples/debruijn/debruijn4.cpp tools/XXX/src/main.cpp

WARNING: some examples use on purpose lambda expressions, so you will need a compiler supporting this feature for this examples.

#Binaries from gatb-core

After the project build, some gatb-core binaries are available here: 'ext/gatb-core/bin'

As gatb-core uses HDF5, you will have here some H5xxx tools from the HDF5 distribution.

You will also find two gatb-core binaries:

    * dbgh5:    builds a DeBruijn graph from a set of reads and save it as a HDF5 file
    * dbginfo:  dumps information about a DeBruinj graph build by dbgh5

