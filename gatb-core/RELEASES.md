--------------------------------------------------------------------------------
# RELEASE 1.2.0

* Assembly-inspired de Bruijn graph simplifications are available using a single command.

Here is an example:

    // removes tips, bubbles and erroneous connections, 
    // similar to SPAdes algorithm
    graph.simplify(); 

* Faster graph traversal using a single command. 

Here is an example:

    // allocates 1 byte/node to precompute adjacency for each nodes 
    // in the MPHF. 
    // Faster graph traversal (especially using neighbors()).
    graph.precomputeAdjacency(); 
       
* **Breaking API changes**

Major changes in API are:

      neighbors\<Node>(..) *becomes* neighbors(..)
      neighbors\<Edge>(..) *becomes* neighborsEdge(..)
      iterator\<Node>(..) *becomes* iterator(..)
      iterator\<BranchingNode>(..) *becomes* iteratorBranching(..) 
      node.kmer.get\<Type>() *becomes* node.template getKmer<Type>()
      successors\<Node>(..) *becomes* successors(..)  
      const Node& *becomes* Node&
           (as MPHF indices are now cached in Node objects)

      etc.. for all fonctions of the type:
      - xxx\<Node>,
      - xxx\<Edge>, 
      - xxx\<BranchingNode>,
      - xxx\<BranchingEdge>.

      
* The basic kmer type (Kmer\<>::Type) no longer has a constructor. Use [kmer].setVal(0) to set the value of the variable [kmer] to zero.

For instance, the following code:

    optimum = Kmer<span>::Type(0)

becomes:

    optimum.setVal(0);

* Graph is now a templated object (GraphTemplate\<Node\_t, Edge\_t, GraphDataVariant\_t>) behind the scenes. However this change is transparent to users of previous versions of GATB-core, as compatibility with the Graph class is preserved.
    
* bug fixes in how queries with dir=DIR_INCOMING are handled.



--------------------------------------------------------------------------------
# RELEASE 1.1.1

* Re-design to support variable number of kmer sizes
 => now, one can use the cmake variable KSIZE_LIST, for instance "cmake -DKSIZE_LIST="32 64 96" ..

* Allows "auto" value for the -abundance-min parameter

--------------------------------------------------------------------------------
# RELEASE 1.1.0

* Re-design of the SortingCountAlgorithm with introduction of interface ICountProcessor
 => it should allow development of new tools based on kmers counting

--------------------------------------------------------------------------------
# RELEASE 1.0.8

* Correction of memory alignment issue on MacOs in some cases

* Re-introduce multi-passes management in DSK

* Correction of passes number configuration with some banks inputs

* Temporary files have now unique names so dbgh5 can be launched several times in the same working directory

--------------------------------------------------------------------------------
# RELEASE 1.0.7

* Correction of scripts for new project creation and delivery process

--------------------------------------------------------------------------------
# RELEASE 1.0.6

* Speed up from x2 to x3 for kmer counting and graph construction phases (optimizations based on minimizers and improved Bloom filters). GATB's k-mer counter has been improved using techniques from KMC2, to achieve competitive running times compared to KMC2.

* Ability to store arbitrary information associated to each kmer of the graph, enabled by a minimal perfect hash function (costs only 2.61 bits/kmer of memory)

* Improved API with new possibilities (banks and kmers management)

* Many new snippets showing how to use the library.


--------------------------------------------------------------------------------
# RELEASE 1.0.5

Modifications of Kmer::Model class for kmers management
* better implementation (factorization, optimization)
* introduction of minimizers concept

WARNING ! These modifications introduced small API changes. Please read snippets kmer2 and kmer5 to see how to handle kmers now.

