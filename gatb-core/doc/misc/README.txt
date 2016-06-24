##################################################################################################
#                                                                                                #
#          #####      #     #######  ######             #####   #######  ######   #######        # 
#         #     #    # #       #     #     #           #     #  #     #  #     #  #              #
#         #         #   #      #     #     #           #        #     #  #     #  #              #
#         #  ####  #     #     #     ######   #######  #        #     #  ######   #####          #
#         #     #  #######     #     #     #           #        #     #  #   #    #              #
#         #     #  #     #     #     #     #           #     #  #     #  #    #   #              #
#          #####   #     #     #     ######             #####   #######  #     #  #######        #
#                                                                                                #
##################################################################################################

The gatb-core project provides a library for denovo genome assembly.

The archive is made of:

    - README.txt : the file you are reading
    - lib        : a static library holding all the assembly services
    - include    : a header files directory
    - doc        : the documentation of the library
    - bin        : 'dbgh5' tool that builds a De Bruijn graph from a set of reads
    - test       : unit tests of the library
    - examples   : a few snippets showing how to use the library

For documentation, please have a look at doc/html/index.html

You can check the library by launching unit tests on your system. Just go into 'test' directory 
and run 'gatb-core-cppunit'. 

You can try the library by compiling a snippet from the 'examples' directory:

for Linux:
    g++ examples/debruijn/debruijn1.cpp -Iinclude -Llib -lgatbcore -lhdf5 -ldl -lz -lpthread  -std=c++0x -O3 -o debruijn1
    
for MacOs:
    clang++ examples/debruijn/debruijn1.cpp -Iinclude -Llib -lgatbcore -lhdf5 -ldl -lz -lpthread  -std=c++0x -O3 -DBOOST_NO_CXX11_RVALUE_REFERENCES=1 -o debruijn1
