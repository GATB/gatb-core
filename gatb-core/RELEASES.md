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
