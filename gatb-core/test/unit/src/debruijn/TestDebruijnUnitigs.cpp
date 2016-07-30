/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#include <CppunitCommon.hpp>

#include <gatb/system/impl/System.hpp>

#include <gatb/tools/math/LargeInt.hpp>

#include <gatb/tools/collections/impl/IteratorFile.hpp>

#include <gatb/tools/misc/api/StringsRepository.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <gatb/debruijn/impl/Graph.hpp>
#include <gatb/debruijn/impl/GraphUnitigs.hpp>
#include <gatb/debruijn/impl/Terminator.hpp>
#include <gatb/debruijn/impl/Traversal.hpp>

#include <gatb/kmer/impl/SortingCountAlgorithm.hpp>
#include <gatb/kmer/impl/BloomAlgorithm.hpp>
#include <gatb/kmer/impl/DebloomAlgorithm.hpp>

#include <gatb/bank/impl/BankStrings.hpp>
#include <gatb/bank/impl/BankSplitter.hpp>
#include <gatb/bank/impl/BankRandom.hpp>

#include <gatb/tools/storage/impl/Storage.hpp>

#include <iostream>
#include <memory>

using namespace std;

using namespace gatb::core::debruijn;
using namespace gatb::core::debruijn::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::math;
using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::storage;
using namespace gatb::core::tools::storage::impl;

extern std::string DBPATH (const string& a);

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for genomic databases management
 */
class TestDebruijnUnitigs : public Test
{
    typedef gatb::core::debruijn::impl::GraphUnitigsTemplate<32> GraphUnitigs;

    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestDebruijnUnitigs);

        CPPUNIT_TEST_GATB (debruijn_unitigs_test7_nocircular);
         //CPPUNIT_TEST_GATB (debruijn_unitigs_deletenode); // probably not appropriate, it's a weird case of a self-revcomp kmer inside a unitig and also at an extremity.
        CPPUNIT_TEST_GATB (debruijn_unitigs_test2);
        CPPUNIT_TEST_GATB (debruijn_unitigs_test4);
        CPPUNIT_TEST_GATB (debruijn_unitigs_test5);
        CPPUNIT_TEST_GATB (debruijn_unitigs_test6);
        //CPPUNIT_TEST_GATB (debruijn_unitigs_test8);// simplePathEdge tests, that api isn't implemented.
        //CPPUNIT_TEST_GATB (debruijn_unitigs_test9);
        CPPUNIT_TEST_GATB (debruijn_unitigs_test13);
        CPPUNIT_TEST_GATB (debruijn_unitigs_build);
        //CPPUNIT_TEST_GATB (debruijn_unitigs_traversal1); // would need to be fixed
        CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {}
    void tearDown ()  {}

    // SMALL VALUE NEEDED because continuous integration servers are not very powerful...
    static const u_int64_t MAX_MEMORY = 1000;

    /********************************************************************************/
    struct Info
    {
        Info() : nbNodes(0), nbNeighbors(0) {}
        Integer checksum;
        size_t nbNodes;
        size_t nbNeighbors;

        /** */
        void incNodes ()      { nbNodes++;  }
        void inc      (Integer& t)  { checksum = checksum + t; nbNeighbors++;  }
        void inc      ()  {  nbNeighbors++;  }

        /** */
        string toString ()
        {
            stringstream ss;
            ss << "[INFO  nbNodes=" << nbNodes << "  nbNeighbors=" << nbNeighbors << "  checksum=" << checksum << "]";
            return ss.str();
        }
    };

    /********************************************************************************/
    void debruijn_unitigs_test2_aux (GraphUnitigs& graph)
    {
        /** We get an iterator over all the nodes of the graph. */
        GraphIterator<NodeFast<32>> itNodes = graph.iterator();

        GraphVector<EdgeFast<32>> successors;
        Info    info;

        TimeInfo ti;
        {
            TIME_INFO (ti, "loop");

            /** We iterate all the nodes of the graph. */
            for (itNodes.first(); !itNodes.isDone(); itNodes.next())
            {
                info.incNodes();

                /** We retrieve the outcoming edges. */
                successors = graph.successorsEdge (itNodes.item());

                /** We iterate all outcoming edges. */
                //for (size_t i=0; i<successors.size(); i++)  {  info.inc (successors[i].to.kmer);  } // some type conflict.. cba.
                for (size_t i=0; i<successors.size(); i++)  {  info.inc ();  } 
            }
        }

        //cout << info.toString() <<  "  time=" << ti.getEntryByKey("loop") << endl;
    }

    /********************************************************************************/
    void debruijn_unitigs_test2_aux (StorageMode_e mode, size_t kmerSize, size_t nks, const char* seq)
    {
        size_t seqLen   = strlen (seq);

        CPPUNIT_ASSERT (seqLen >= kmerSize);

        /** We create a bank with one sequence. */
        IBank* bank = new BankStrings (seq, 0);

        /** We configure parameters for a SortingCountAlgorithm object. */
        IProperties* params = SortingCountAlgorithm<>::getDefaultProperties();
        params->setInt (STR_KMER_SIZE,          kmerSize);
        params->setInt (STR_MAX_MEMORY,         MAX_MEMORY);
        params->setInt (STR_KMER_ABUNDANCE_MIN, nks);
        params->setStr (STR_URI_OUTPUT,         "dummy");

        /** We create a DSK instance. */
        SortingCountAlgorithm<> sortingCount (bank, params);

        /** We launch DSK. */
        sortingCount.execute();

        /** We get the storage instance. */
        Storage* storage = sortingCount.getStorage();

        /** We check that the sequence has no duplicate kmers. */
        CPPUNIT_ASSERT ( (int64_t) (seqLen - kmerSize + 1) == sortingCount.getSolidCounts()->getNbItems());

        }

    /********************************************************************************/
    void debruijn_unitigs_test2 ()
    {
        const char* sequences [] =
        {
            "ACCATGTATAATTATAAGTAGGTACCTATTTTTTTATTTTAAACTGAAAT",
            "CGCTACAGCAGCTAGTTCATCATTGTTTATCAATGATAAAATATAATAAGCTAAAAGGAAACTATAAATA",
            "CGCTATTCATCATTGTTTATCAATGAGCTAAAAGGAAACTATAAATAACCATGTATAATTATAAGTAGGTACCTATTTTTTTATTTTAAACTGAAATTCAATATTATATAGGCAAAG"
        };

        size_t kmerSizes[] = { 15, 23, 31};
        size_t nks=1;

        for (size_t i=0; i<ARRAY_SIZE(sequences); i++)
        {
            for (size_t j=0; j<ARRAY_SIZE(kmerSizes); j++)
            {
                debruijn_unitigs_test2_aux (STORAGE_HDF5, kmerSizes[j], nks, sequences[i]);
            }
        }
    }
    /********************************************************************************/
    void debruijn_unitigs_test4 ()
    {
        char* seq = (char*) "ACCATGTATAATTATAAGTAGGTACCT";  // size 27
        char* rev = (char*) "AGGTACCTACTTATAATTATACATGGT";

        /** We create the graph. */
        GraphUnitigs graph = GraphUnitigs::create (new BankStrings (seq, 0), "-kmer-size 27  -abundance-min 1  -verbose 0  -max-memory %d", MAX_MEMORY);

        GraphIterator<NodeFast<32>> it = graph.iterator();  it.first();

        NodeFast<32> n1 = it.item();
        CPPUNIT_ASSERT (n1.strand == STRAND_FORWARD);
        CPPUNIT_ASSERT (graph.toString(n1).compare (seq) == 0);

        NodeFast<32> n2 = graph.reverse (n1);
        CPPUNIT_ASSERT (n2.strand == STRAND_REVCOMP);
        CPPUNIT_ASSERT (graph.toString(n2).compare (rev) == 0);
    }

    /********************************************************************************/
    void debruijn_unitigs_test5 ()
    {
        /** We create an empty graph with a given kmer size. */
        GraphUnitigs graph = GraphUnitigs::create (7);

        /** We define a string of size equal to the kmer size. */
        char* seq = (char*) "ACCAGTT";
        char* rev = (char*) "AACTGGT";

        NodeFast<32> n1 = graph.buildNode (Data(seq));
        CPPUNIT_ASSERT (graph.toString (n1) == seq);

        NodeFast<32> n2 = graph.reverse (n1);
        CPPUNIT_ASSERT (graph.toString (n2) == rev);
    }

    /********************************************************************************/
    struct debruijn_unitigs_test6_fct
    {
        debruijn_unitigs_test6_fct (const GraphUnitigs& graph) : graph(graph) {}
        const GraphUnitigs& graph;

        void operator() (NodeFast<32>& node) const
        {
            string snode = graph.toString (node);

            NodeFast<32> rev = graph.reverse (node);
            string srev = graph.toString (rev);

            /** We build a node from the reverse string. */
            NodeFast<32> rev2 = graph.buildNode (Data ((char*)srev.c_str()));
            CPPUNIT_ASSERT (graph.toString(rev2) == srev);

            /** We reverse the reversed node. */
            NodeFast<32> node2 = graph.reverse (rev2);
            CPPUNIT_ASSERT (graph.toString(node2) == snode);
        }
    };

    void debruijn_unitigs_test6 ()
    {
        char* seq = (char*) "ACCATGTATAATTATAAGTAGGTACCACGATCGATCGATCGATCGTAGCATATCGTACGATCT";

        /** We create the graph. */
        GraphUnitigs graph = GraphUnitigs::create (new BankStrings (seq, 0), "-kmer-size 27  -abundance-min 1  -verbose 0  -max-memory %d", MAX_MEMORY);

        debruijn_unitigs_test6_fct fct(graph);
        graph.iterator().iterate (fct);
    }

    /********************************************************************************/
    struct debruijn_unitigs_test7_fct
    {
        debruijn_unitigs_test7_fct (const GraphUnitigs& graph, NodeFast<32>& n1, NodeFast<32>& n2) : graph(graph), n1(n1), n2(n2) {}
        const GraphUnitigs& graph;
        NodeFast<32>& n1;
        NodeFast<32>& n2;

        void operator() (NodeFast<32>& current) const
        {
            string currentStr = graph.toString(current);
            
            if (!(currentStr==graph.toString(n1) || currentStr==graph.toString(n2) ))
                return;
            //CPPUNIT_ASSERT (currentStr==graph.toString(n1) || currentStr==graph.toString(n2) ); // restore once dummy node to address MPHF bug is resolved (see below)

            /** We get all possible edges from the current node */
            GraphVector<EdgeFast<32>> neighbors = graph.neighborsEdge(current);
            
            CPPUNIT_ASSERT (neighbors.size()>=1);

            for (size_t i=0; i<neighbors.size(); i++)
            {
                /** Shortcut. */
                EdgeFast<32>& edge = neighbors[i];

                //std::cout << "curstr "  << currentStr << std::endl;
                if (currentStr==graph.toString(n1))  // 1 neighbor
                {
                    if (neighbors.size() != 1) 
                        std::cout << "anticipation of assert fail: neighbors size of n1 " << neighbors.size() << std::endl;
                    if (graph.toString(edge.to) != "GGCGC") 
                        std::cout << "anticipation of assert fail: graph.toString(edge.to) = " << graph.toString(edge.to) << std::endl;

                    CPPUNIT_ASSERT (neighbors.size()==1);
                    CPPUNIT_ASSERT (edge.nt==NUCL_C);
                    CPPUNIT_ASSERT (edge.direction==DIR_OUTCOMING);
                    CPPUNIT_ASSERT (graph.toString(edge.from)=="AGGCG");
                    CPPUNIT_ASSERT (graph.toString(edge.to)  =="GGCGC");
                }

                if (currentStr==graph.toString(n2))  // 2 neighbors
                {
                    CPPUNIT_ASSERT (neighbors.size()==2);

                    if (graph.toString(edge.to) == "CGCCT")
                    {
                        if (graph.toString(edge.to) != "CGCCT") 
                            std::cout << "anticipation of assert fail: graph.toString(edge.to) = " << graph.toString(edge.to) << std::endl;
                        if (edge.nt != NUCL_T) 
                            std::cout << "anticipation of assert fail: edge.nt = " << (int)edge.nt << " and not T" << std::endl;
                        if (edge.direction != DIR_OUTCOMING) 
                            std::cout << "anticipation of assert fail: edge.direction = " << (int)edge.direction << std::endl;

                        CPPUNIT_ASSERT (graph.toString(edge.to)=="CGCCT");
                        CPPUNIT_ASSERT (edge.nt==NUCL_T);
                        CPPUNIT_ASSERT (edge.direction==DIR_OUTCOMING);

                    }
                    else
                    { if (graph.toString(edge.to)=="GGCGC")
                        {
                            if (graph.toString(edge.from) != "GCGCC") 
                                std::cout << "anticipation of assert fail: graph.toString(edge.from) = " << graph.toString(edge.from) << std::endl;

                            CPPUNIT_ASSERT (graph.toString(edge.from)=="GCGCC");
                            CPPUNIT_ASSERT (edge.nt==NUCL_G);
                            CPPUNIT_ASSERT (edge.direction==DIR_INCOMING);
                        }
                        else
                        {
                            std::cout << "unknown edge found" << std::endl;
                            std::cout << "anticipation of assert fail: graph.toString(edge.to) = " << graph.toString(edge.to) << " graph.toString(edge.from) = " << graph.toString(edge.from) << std::endl;
                            CPPUNIT_ASSERT (false);
                        }
                    }
                }
            }


            GraphVector<EdgeFast<32>> neighborsIncoming = graph.neighborsEdge(current, DIR_INCOMING);

            if (currentStr==graph.toString(n1))  // 0 incoming neighbor
            {
                if (neighborsIncoming.size() != 0) 
                    std::cout << "anticipation of assert fail: incoming neighbors size of n1 " << neighbors.size() << std::endl;
                CPPUNIT_ASSERT (neighborsIncoming.size()==0);
            }


            for (size_t i=0; i<neighborsIncoming.size(); i++)
            { 
                if (currentStr==graph.toString(n2))  // 1 incoming neighbor
                {
                    /** Shortcut. */
                    EdgeFast<32>& edge = neighborsIncoming[i];
                    CPPUNIT_ASSERT (neighborsIncoming.size()==1);
                    CPPUNIT_ASSERT (graph.toString(edge.to)=="GGCGC");
                    CPPUNIT_ASSERT (edge.nt==NUCL_G);
                    CPPUNIT_ASSERT (edge.direction==DIR_INCOMING);
                }
            }
        }
    };

    void debruijn_unitigs_test7 ()
    {
        /** We create the graph. */
        GraphUnitigs graph = GraphUnitigs::create (new BankStrings ("AGGCGC", "ACTGACTGACTGACTG",0),  "-kmer-size 5  -abundance-min 1  -verbose 1 -max-memory %d -out dummy -minimizer-size 3", MAX_MEMORY);

        /** We should get two kmers:
         *      - AGGCG / CGCCT
         *      - GCGCC / GGCGC
         */
        NodeFast<32> n1 = graph.buildNode ((char*)"AGGCG");
        NodeFast<32> n2 = graph.buildNode ((char*)"GCGCC");

        // We should get as neighborhood
        // GCGCC  [GCGCC --T--> CGCCT]
        // GCGCC  [GGCGC --C--> GCGCC]
        // AGGCG  [AGGCG --C--> GGCGC]

        debruijn_unitigs_test7_fct fct (graph, n1, n2);
        graph.iterator().iterate (fct);

#ifdef WITH_MPHF
        /* rerun this test with adjacency information instead of bloom */
        graph.precomputeAdjacency(1, false);
        
        graph.iterator().iterate (fct);
#endif

    }

    void debruijn_unitigs_test7_nocircular ()
    {
        // same as test7, but without the circular contig
        GraphUnitigs graph = GraphUnitigs::create (new BankStrings ("AGGCGC", "ACTGACT",0),  "-kmer-size 5  -abundance-min 1  -verbose 0 -max-memory %d -out dummy -minimizer-size 3", MAX_MEMORY);

        /** We should get two kmers:
         *      - AGGCG / CGCCT
         *      - GCGCC / GGCGC
         */
        NodeFast<32> n1 = graph.buildNode ((char*)"AGGCG");
        NodeFast<32> n2 = graph.buildNode ((char*)"GCGCC");

        // We should get as neighborhood
        // GCGCC  [GCGCC --T--> CGCCT]
        // GCGCC  [GGCGC --C--> GCGCC]
        // AGGCG  [AGGCG --C--> GGCGC]

        debruijn_unitigs_test7_fct fct (graph, n1, n2);
        graph.iterator().iterate (fct);
    }



    /********************************************************************************/
    void debruijn_unitigs_test8_aux (char* seq, size_t kmerSize)
    {
        /** We create the graph. */
        GraphUnitigs graph = GraphUnitigs::create (new BankStrings (seq, NULL),  "-kmer-size %d  -abundance-min 1  -verbose 0  -max-memory %d", kmerSize, MAX_MEMORY);

        // We get the first node.
        NodeFast<32> node = graph.buildNode (seq);

        GraphIterator<EdgeFast<32>> path = graph.simplePathEdge (node, DIR_OUTCOMING);

        for (path.first(); !path.isDone(); path.next())
        {
            /** We check that the current transition nucleotide matches the correct character in the sequence. */
            CPPUNIT_ASSERT (ascii(path.item().nt) == seq [graph.getKmerSize() + path.rank()]);
        }

        /** We check that we found the correct number of nodes. */
        CPPUNIT_ASSERT (path.rank() == strlen (seq) - kmerSize);
    }

    /** */
    void debruijn_unitigs_test8 ()
    {
        /** This sequence should not have branching nodes for kmer size big enough. */
        char* seq = (char*) "AGGCGCTAGGGTAGAGGATGATGA";

        size_t kmerSizes[] = {11, 13, 15, 17};

        for (size_t i=0; i<ARRAY_SIZE(kmerSizes); i++)  { debruijn_unitigs_test8_aux (seq, kmerSizes[i]); }
    }

    /********************************************************************************/
    void debruijn_unitigs_test9 ()
    {
        size_t kmerSize = 9;

        char* seq1 = (char*) "AGGCGCTAGGGTAGAGGATGATGA";
        char* seq2 = (char*) "AGGCGCTAGGGTATAGGATGATGA";
        //                    000000000011111111112222
        //                    012345678901234567890123
        //  difference here                ^

        /** We create the graph. */
        GraphUnitigs graph = GraphUnitigs::create (new BankStrings (seq1, seq2, NULL),  "-kmer-size %d  -abundance-min 1  -verbose 0  -max-memory %d ", kmerSize, MAX_MEMORY);

        /** We get the first node. */
        NodeFast<32> node = graph.buildNode (seq1);

        /** We get a simple path iterator starting from the beginning of the seq1. */
        GraphIterator<EdgeFast<32>> path = graph.simplePathEdge (node, DIR_OUTCOMING);

        for (path.first(); !path.isDone(); path.next())
        {
            CPPUNIT_ASSERT (graph.isSimple (path.item()));
            CPPUNIT_ASSERT (ascii(path.item().nt) == seq1 [graph.getKmerSize() + path.rank()]);
        }

        /** We check that we stopped at the first difference between the two sequences. */
        CPPUNIT_ASSERT (path.rank() == 4);   // 4 = diffOffset - kmerSize
    }

    /********************************************************************************/
    void debruijn_unitigs_test13 ()
    {
        size_t   readSize = 80;
        u_int8_t coverage = 5;
        size_t   kmerSize = 31;
        size_t   nks      = coverage;

        const char* seq =
            "GAATTCCAGGAGGACCAGGAGAACGTCAATCCCGAGAAGGCGGCGCCCGCCCAGCAGCCCCGGACCCGGGCTGGACTGGC"
            "GGTACTGAGGGCCGGAAACTCGCGGGGTCCAGCTCCCCAGAGGCCTAAGACGCGACGGGTTGCACCTCTTAAGGATCTTC"
            "CTATAAATGATGAGTATGTCCCTGTTCCTCCCTGGAAAGCAAACAATAAACAGCCTGCATTTACCATACATGTGGATGAA"
            "GCAGAAGAAATTCAAAAGAGGCCAACTGAATCTAAAAAATCAGAAAGTGAAGATGTCTTGGCCTTTAATTCAGCTGTTAC"
            "TTTACCAGGACCAAGAAAGCCACTGGCACCTCTTGATTACCCAATGGATGGTAGTTTTGAGTCTCCACATACTATGGAAA"
            "TGTCAGTTGTATTGGAAGATGAAAAGCCAGTGAGTGTTAATGAAGTACCAGACTACCATGAGGACATTCACACGTACCTT"
            "AGGGAAATGGAGGTTAAATGTAAGCCTAAAGTGGGTTACATGAAGAAACAGCCAGACATTACTAACAGTATGAGGGCTAT"
            "CCTCGTGGACTGGTTAGTTGAAGTAGGAGAAGAATATAAACTGCAGAACGAGACCCTGCATTTGGCTGTGAACTACATTG"
            "ATAGGTTTCTTTCATCCATGTCTGTGTTGAGAGGAAAACTTCAACTTGTGGGCACTGCTGCTATGCTTTTAGCCTCAAAG"
            "TTTGAAGAGATATACCCGCCAGAAGTAGCAGAGTTTGTATACATTACAGATGACACTTATACCAAGAAACAAGTTCTAAG"
            "GATGGAGCACCTAGTCTTGAAAGTCCTGGCTTTTGACTTAGCTGCACCAACAATAAATCAGTTTCTTACCCAGTACTTTT"
            "TGCATCAGCAGCCTGCAAACTGCAAAGTTGAAAGTTTAGCAATGTTTTTGGGAGAGTTAAGTTTGATAGATGCTGACCCA"
            "TATCTAAAGTATTTGCCGTCAGTTATCGCTGCAGCAGCCTTTCATTTAGCACTCTACACAGTCACAGGACAAAGCTGGCC"
            "TGAATCATTAGTACAGAAGACTGGATATACTCTGGAAACTCTAAAGCCTTGTCTCCTGGACCTTCACCAGACCTACCTCA"
            "GAGCACCACAGCACGCACAACAGTCAATAAGAGAGAAGTACAAAAATTCAAAGTATCATGGTGTTTCTCTCCTCAACCCA"
            "CCAGAGACACTAAATGTGTAACAGTGAAAAGACTGCCTTTGTTTTCTAAGACTGTAAATCACATGCAATGTATATGGTGT"
            "ACAGATTTTATCTTAGGTTTTAATTTTACAACATTTCTGAATAAGAAGAGTTATGGTCCAGTACAAATTATGGTATCTAT"
            "TACTTTTTAAATGGTTTTAATTTGTATATCTTTTGTAAATGTAACTATCTTAGATATTTGGCTAATTTTAAGTGGTTTCT";

        // We create a bank of reads from a given long sequence
        IBank* bank = new BankSplitter (new BankStrings (seq, NULL), readSize, kmerSize-1, coverage);

        // We create the graph.
        GraphUnitigs graph = GraphUnitigs::create (bank,  "-kmer-size %d  -abundance-min %d  -verbose 0  -max-memory %d", kmerSize, nks, MAX_MEMORY);

        // We check we got the correct number of solid kmers.
        CPPUNIT_ASSERT (graph.getInfo().getInt ("kmers_nb_solid") == (int) (strlen(seq) - kmerSize + 1));

    }

    /********************************************************************************/
    struct debruijn_unitigs_build_entry
    {
        debruijn_unitigs_build_entry () : nbNodes(0), nbBranchingNodes(0) {}
        size_t  nbNodes;
        Integer checksumNodes;
        size_t  nbBranchingNodes;
        Integer checksumBranchingNodes;
    };

    debruijn_unitigs_build_entry debruijn_unitigs_build_aux_aux (IBank* inputBank, bool checkNodes, bool checkBranching)
    {
        debruijn_unitigs_build_entry result;

        /** We load the graph. */
        GraphUnitigs graph = GraphUnitigs::create (inputBank,  "-kmer-size 31 -out %s -abundance-min 1  -verbose 0  -max-memory %d",                        "dummy", MAX_MEMORY);

        if (checkNodes)
        {
            GraphIterator<NodeFast<32>> iterNodes = graph.iterator();
            for (iterNodes.first(); !iterNodes.isDone(); iterNodes.next())
            { result.nbNodes++; /*result.checksumNodes = iterNodes.item().kmer + result.checksumNodes; */ /* disabled because of unsupported largeint operation apparently! */ }
        }

        if (checkBranching)
        {
        }

        return result;
    }

    /********************************************************************************/
    void debruijn_unitigs_build_aux (const char* sequences[], size_t nbSequences)
    {
        // We build the bank
        IBank* inputBank = new BankStrings (sequences, nbSequences);
        LOCAL (inputBank);

        debruijn_unitigs_build_entry r1 = debruijn_unitigs_build_aux_aux (inputBank, true,  true);
    }

    /********************************************************************************/
    void debruijn_unitigs_build ()
    {
        const char* sequences[] =
        {
            "GAATTCCAGGAGGACCAGGAGAACGTCAATCCCGAGAAGGCGGCGCCCGCCCAGCAGCCCCGGACCCGGGCTGGACTGGC",
            "GGTACTGAGGGCCGGAAACTCGCGGGGTCCAGCTCCCCAGAGGCCTAAGACGCGACGGGTTGCACCTCTTAAGGATCTTC",
            "CTATAAATGATGAGTATGTCCCTGTTCCTCCCTGGAAAGCAAACAATAAACAGCCTGCATTTACCATACATGTGGATGAA",
            "GCAGAAGAAATTCAAAAGAGGCCAACTGAATCTAAAAAATCAGAAAGTGAAGATGTCTTGGCCTTTAATTCAGCTGTTAC",
            "TTTACCAGGACCAAGAAAGCCACTGGCACCTCTTGATTACCCAATGGATGGTAGTTTTGAGTCTCCACATACTATGGAAA",
            "TGTCAGTTGTATTGGAAGATGAAAAGCCAGTGAGTGTTAATGAAGTACCAGACTACCATGAGGACATTCACACGTACCTT",
            "AGGGAAATGGAGGTTAAATGTAAGCCTAAAGTGGGTTACATGAAGAAACAGCCAGACATTACTAACAGTATGAGGGCTAT",
            "CCTCGTGGACTGGTTAGTTGAAGTAGGAGAAGAATATAAACTGCAGAACGAGACCCTGCATTTGGCTGTGAACTACATTG",
            "ATAGGTTTCTTTCATCCATGTCTGTGTTGAGAGGAAAACTTCAACTTGTGGGCACTGCTGCTATGCTTTTAGCCTCAAAG",
            "TTTGAAGAGATATACCCGCCAGAAGTAGCAGAGTTTGTATACATTACAGATGACACTTATACCAAGAAACAAGTTCTAAG",
            "GATGGAGCACCTAGTCTTGAAAGTCCTGGCTTTTGACTTAGCTGCACCAACAATAAATCAGTTTCTTACCCAGTACTTTT",
            "TGCATCAGCAGCCTGCAAACTGCAAAGTTGAAAGTTTAGCAATGTTTTTGGGAGAGTTAAGTTTGATAGATGCTGACCCA",
            "TATCTAAAGTATTTGCCGTCAGTTATCGCTGCAGCAGCCTTTCATTTAGCACTCTACACAGTCACAGGACAAAGCTGGCC",
            "TGAATCATTAGTACAGAAGACTGGATATACTCTGGAAACTCTAAAGCCTTGTCTCCTGGACCTTCACCAGACCTACCTCA",
            "GAGCACCACAGCACGCACAACAGTCAATAAGAGAGAAGTACAAAAATTCAAAGTATCATGGTGTTTCTCTCCTCAACCCA",
            "CCAGAGACACTAAATGTGTAACAGTGAAAAGACTGCCTTTGTTTTCTAAGACTGTAAATCACATGCAATGTATATGGTGT",
            "ACAGATTTTATCTTAGGTTTTAATTTTACAACATTTCTGAATAAGAAGAGTTATGGTCCAGTACAAATTATGGTATCTAT",
            "TACTTTTTAAATGGTTTTAATTTGTATATCTTTTGTAAATGTAACTATCTTAGATATTTGGCTAATTTTAAGTGGTTTCT"
        };

        debruijn_unitigs_build_aux (sequences, ARRAY_SIZE(sequences));
    }

    /********************************************************************************/

    void debruijn_unitigs_traversal1_aux_aux (bool useCopyTerminator, size_t kmerSize, const char** seqs, size_t seqsSize,
		TraversalKind traversalKind, const char* checkStr
	)
    {
        // We create a fake bank with a SNP
        IBank* bank = new BankStrings (seqs, seqsSize);

        // We load the graph
        GraphUnitigs graph = GraphUnitigs::create (bank, "-abundance-min 1  -verbose 0  -kmer-size %d  -max-memory %d",
			kmerSize, MAX_MEMORY);

        // We create a Terminator object
        MPHFTerminatorTemplate<NodeFast<32>,EdgeFast<32>,GraphUnitigs> terminator (graph);

		// We create a node from the start of the first sequence
		NodeFast<32> node = graph.buildNode (seqs[0]);

        Path_t<NodeFast<32>> path;

        // We create a Traversal instance according to the chosen traversal kind
        TraversalTemplate<NodeFast<32>,EdgeFast<32>,GraphUnitigs>* traversal = TraversalTemplate<NodeFast<32>,EdgeFast<32>,GraphUnitigs>::create (traversalKind, graph, terminator);
        LOCAL (traversal);

        traversal->traverse (node, DIR_OUTCOMING, path);


        stringstream ss;  ss << graph.toString (node) << path ;
        CPPUNIT_ASSERT (ss.str().compare(checkStr)==0);
    }

    void debruijn_unitigs_traversal1_aux (bool useCopyTerminator)
    {
        const char* seqs[] =
        {
            "CGCTACAGCAGCTAGTTCATCATTGTTTATCAATGATAAAATATAATAAGCTAAAAGGAAACTATAAATA",
            "CGCTACAGCAGCTAGTTCATCATTGTTTATCGATGATAAAATATAATAAGCTAAAAGGAAACTATAAATA"
            //      SNP HERE at pos 31      x
        };

    	debruijn_unitigs_traversal1_aux_aux (useCopyTerminator, 15, seqs, ARRAY_SIZE(seqs), TRAVERSAL_UNITIG,
			"CGCTACAGCAGCTAGTTCATCATTGTTTATC"
		);

    	debruijn_unitigs_traversal1_aux_aux (useCopyTerminator, 15, seqs, ARRAY_SIZE(seqs), TRAVERSAL_CONTIG,
			"CGCTACAGCAGCTAGTTCATCATTGTTTATCAATGATAAAATATAATAAGCTAAAAGGAAACTATAAATA"
		);
    }

    void debruijn_unitigs_traversal1 ()
    {
    	debruijn_unitigs_traversal1_aux (false);
    }




    /********************************************************************************/
    void debruijn_unitigs_deletenode_fct (const GraphUnitigs& graph) 
    {
        NodeFast<32> n1 = graph.buildNode ((char*)"AGGCG");
        //NodeFast<32> n2 = graph.buildNode ((char*)"GGCGC"); // hum it's the same as n3..
        NodeFast<32> n3 = graph.buildNode ((char*)"GCGCC");

        graph.deleteNode(n3);

        /** We get all possible edges from the current kmer (note: not from the current node). */
        GraphVector<EdgeFast<32>> neighbors = graph.neighborsEdge(n1.kmer);

        CPPUNIT_ASSERT (neighbors.size()==0);
    }

    void debruijn_unitigs_deletenode ()
    {
        // MPHF has a known bug where, when there are only like a tiny amount of elements (I tested with three), it will just return mphf(elt)=0 always.
        // so this is why I'm adding the dummy "ACTGACTGACTGACTG" sequence, to artificially increase the amount of elements in the mphf
        GraphUnitigs graph = GraphUnitigs::create (new BankStrings ("AGGCGCC", "ACTGACTGACTGACTG",0),  "-kmer-size 5  -abundance-min 1  -verbose 0  -max-memory %d", MAX_MEMORY);

        debruijn_unitigs_deletenode_fct (graph);
    }

    void debruijn_unitigs_deletenode2_fct (const GraphUnitigs& graph) 
    {
        NodeFast<32> n1 = graph.buildNode ((char*)"AGGCG");
        NodeFast<32> n2 = graph.buildNode ((char*)"GGCGA");

        graph.deleteNode(n2);

        /** We get all possible edges from the current kmer (note: not from the current node). */
        GraphVector<EdgeFast<32>> neighbors = graph.neighborsEdge(n1.kmer);

        std::cout<<"unfinished test" << std::endl;
    return; // TODO: finish it;

        CPPUNIT_ASSERT (neighbors.size()==0);

        for (size_t i=0; i<neighbors.size(); i++)
        {
            /** Shortcut. */
            EdgeFast<32>& edge = neighbors[i];

            if (neighbors.size() != 1) 
                std::cout << "anticipation of assert fail: neighbors size of n1 " << neighbors.size() << std::endl;
            if (graph.toString(edge.to) != "GGCGC") 
                std::cout << "anticipation of assert fail: graph.toString(edge.to) = " << graph.toString(edge.to) << std::endl;

            CPPUNIT_ASSERT (neighbors.size()==1);
            CPPUNIT_ASSERT (edge.nt==NUCL_C);
            CPPUNIT_ASSERT (edge.direction==DIR_OUTCOMING);
            CPPUNIT_ASSERT (graph.toString(edge.from)=="AGGCG");
            CPPUNIT_ASSERT (graph.toString(edge.to)  =="GGCGC");
        }
    }

    void debruijn_unitigs_deletenode2 ()
    {
        GraphUnitigs graph = GraphUnitigs::create (new BankStrings ("AGGCGAAGGCGT", "ACTGACTGACTGACTG",0),  "-kmer-size 5  -abundance-min 1  -verbose 0  -max-memory %d -mphf emphf", MAX_MEMORY);

        debruijn_unitigs_deletenode_fct (graph);

        /* rerun this test with adjacency information instead of bloom */
        
        GraphUnitigs graph2 = GraphUnitigs::create (new BankStrings ("AGGCGAAGGCGT", "ACTGACTGACTGACTG",0),  "-kmer-size 5  -abundance-min 1  -verbose 0  -max-memory %d -mphf emphf", MAX_MEMORY);
        graph2.precomputeAdjacency(1, false);
        
        debruijn_unitigs_deletenode2_fct (graph2);
    }



};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestDebruijnUnitigs);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestDebruijnUnitigs);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

