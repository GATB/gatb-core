//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*                                                                              */
/* Graph loading from a HDF5 file, the iterate over the nodes, parallel version.*/
/*                                                                              */
/* Cmd-line: debruijn7 <h5 file>                                                */
/*                                                                              */
/* Sample: debruijn7 gatb-core/gatb-core/test/db/celegans_reads.h5              */
/*                                                                              */
/* WARNING ! THIS SNIPPET SHOWS ALSO HOW TO USE LAMBDA EXPRESSIONS, SO YOU NEED */
/* TO USE A COMPILER THAT SUPPORTS THIS FEATURE.                                */
/*                                                                              */
/* Note:                                                                        */
/*     - '.h5' file contains the HDF5 formatted representation of a de bruijn   */
/*     graph created from a set of reads.                                       */
/*     - a '.h5' file is created using dbgh5 program provided with GATB-Core.   */
/*     Basic use is as follows:                                                 */
/*        dbgh5 -in <fasta/q file> -out <h5 file>                               */
/*     You can also control kmer-size and kmer abundance, see dbgh5 help.       */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We check that the user provides at least one option (supposed to be in HDF5 format).
    if (argc < 2)
    {
        std::cerr << "You must provide a HDF5 file." << std::endl;
        return EXIT_FAILURE;
    }

    try
    {
        // We load the graph from the given graph file
        Graph graph = Graph::load (argv[1]);

        // We get an iterator for all nodes of the graph.
        GraphIterator<Node> it = graph.iterator ();

        // We choose how many cores we want to use.
        // By convention, a 0 value will use all available cores.
        size_t nbCores = 4;

        // We iterate all the nodes of the graph.

        // Note how we create a Dispatcher object (with the required number of cores) and
        // use it for dispatching the iteration in several threads.
        IDispatcher::Status status = Dispatcher(nbCores).iterate (it, [&graph] (const Node& node)
        {
            // This instruction block will be executed within the context of one of the 4 created threads.

            // We are given a current node and we can do some process on it in the current thread.
        });

        // We display some information about the dispatching
        std::cout << "we used " << status.nbCores << " cores, elapsed time is " << status.time << " msec" << std::endl;
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
