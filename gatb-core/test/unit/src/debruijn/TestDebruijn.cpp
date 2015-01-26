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
class TestDebruijn : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestDebruijn);

        CPPUNIT_TEST_GATB (debruijn_test2);
        CPPUNIT_TEST_GATB (debruijn_test3);
        CPPUNIT_TEST_GATB (debruijn_test4);
        CPPUNIT_TEST_GATB (debruijn_test5);
        CPPUNIT_TEST_GATB (debruijn_test6);
        CPPUNIT_TEST_GATB (debruijn_test7);
        CPPUNIT_TEST_GATB (debruijn_test8);
        CPPUNIT_TEST_GATB (debruijn_test9);
        CPPUNIT_TEST_GATB (debruijn_test10);
        CPPUNIT_TEST_GATB (debruijn_test11);
        CPPUNIT_TEST_GATB (debruijn_test12);
        CPPUNIT_TEST_GATB (debruijn_test13);
        CPPUNIT_TEST_GATB (debruijn_mutation);
        CPPUNIT_TEST_GATB (debruijn_build);
        CPPUNIT_TEST_GATB (debruijn_checksum);
        CPPUNIT_TEST_GATB (debruijn_checkbranching);
#ifdef WITH_MPHF
        CPPUNIT_TEST_GATB (debruijn_mphf);
#endif

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

        /** */
        string toString ()
        {
            stringstream ss;
            ss << "[INFO  nbNodes=" << nbNodes << "  nbNeighbors=" << nbNeighbors << "  checksum=" << checksum << "]";
            return ss.str();
        }
    };

    /********************************************************************************/
    void debruijn_test2_aux (Graph& graph)
    {
        /** We get an iterator over all the nodes of the graph. */
        Graph::Iterator<Node> itNodes = graph.iterator<Node>();

        Graph::Vector<Edge> successors;
        Info    info;

        TimeInfo ti;
        {
            TIME_INFO (ti, "loop");

            /** We iterate all the nodes of the graph. */
            for (itNodes.first(); !itNodes.isDone(); itNodes.next())
            {
                info.incNodes();

                /** We retrieve the outcoming edges. */
                successors = graph.successors<Edge> (itNodes.item());

                /** We iterate all outcoming edges. */
                for (size_t i=0; i<successors.size(); i++)  {  info.inc (successors[i].to.kmer);  }
            }
        }

        //cout << info.toString() <<  "  time=" << ti.getEntryByKey("loop") << endl;
    }

    /********************************************************************************/
    void getNearestBranchingRange (const Graph& graph, const Node& node, Node& begin, Node& end) const
    {
        begin = node;    for (Graph::Vector<Node> nodes ; (nodes = graph.predecessors<Node> (begin)).size() == 1;  begin = nodes[0])  {}
        end   = node;    for (Graph::Vector<Node> nodes ; (nodes = graph.successors  <Node> (end  )).size() == 1;  end   = nodes[0])  {}
    }

    void debruijn_check_sequence (const Graph& graph, size_t kmerSize, const char* seq)
    {
        size_t seqLen = strlen (seq);

        Graph::Iterator<Node> nodes = graph.iterator<Node> ();
        nodes.first ();

        /** We get the first node. */
        Node node = nodes.item();

        /** We compute the branching range for the node. */
        Node begin, end;     getNearestBranchingRange (graph, node, begin, end);

        /** We check that the begin kmer matches the beginning of the sequence. */
        bool check1 =
            graph.toString (begin)        == string (seq, kmerSize)  ||
            graph.toString (graph.reverse(end)) == string (seq, kmerSize);

        /** We check that the end kmer matches the end of the sequence. */
        bool check2 =
            graph.toString (end)            == string (seq + seqLen - kmerSize, kmerSize)  ||
            graph.toString (graph.reverse(begin)) == string (seq + seqLen - kmerSize, kmerSize);

        if (!check1 || !check2)
        {
            cout << "kmerSize=" << kmerSize << endl;
            cout << graph.debugString (node,  STRAND_FORWARD) << "  " << graph.debugString (node,  STRAND_REVCOMP) << endl;
            cout << graph.debugString (begin, STRAND_FORWARD) << "  " << graph.debugString (end,   STRAND_FORWARD) << endl;
            cout << graph.debugString (end,   STRAND_REVCOMP) << "  " << graph.debugString (begin, STRAND_REVCOMP) << endl;
        }

        CPPUNIT_ASSERT (check1 && check2);
    }

    /********************************************************************************/
    void debruijn_test2_aux (StorageMode_e mode, size_t kmerSize, size_t nks, const char* seq)
    {
        size_t seqLen   = strlen (seq);

        CPPUNIT_ASSERT (seqLen >= kmerSize);

        /** We create a bank with one sequence. */
        IBank* bank = new BankStrings (seq, 0);

        /** We create a storage instance. */
        Storage* storage  = StorageFactory(mode).create ("foo", true, true);
        LOCAL (storage);

        /** We create a DSK instance. */
        SortingCountAlgorithm<> sortingCount (storage, bank, kmerSize, make_pair(nks,0));

        /** We launch DSK. */
        sortingCount.execute();

        /** We check that the sequence has no duplicate kmers. */
        CPPUNIT_ASSERT ( (seqLen - kmerSize + 1) == sortingCount.getSolidCounts()->getNbItems());

        /** We create a bloom instance. */
        float nbitsPerKmer = DebloomAlgorithm<>::getNbBitsPerKmer (kmerSize, DEBLOOM_ORIGINAL);
        BloomAlgorithm<> bloom (*storage, sortingCount.getSolidCounts(), kmerSize, nbitsPerKmer, 0, BLOOM_BASIC);
        bloom.execute ();

        /** We create a debloom instance. */
        DebloomAlgorithm<> debloom (*storage, *storage, sortingCount.getSolidCounts(), kmerSize);

        /** We launch the debloom. */
        debloom.execute();
    }

    /********************************************************************************/
    void debruijn_test2 ()
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
                debruijn_test2_aux (STORAGE_HDF5, kmerSizes[j], nks, sequences[i]);
            }
        }
    }

    /********************************************************************************/
    void debruijn_test3 ()
    {
        const char* sequences [] =
        {
            "ACCATGTATAATTATAAGTAGGTACCTATTTTTTTATTTTAAACTGAAAT",
            "CGCTACAGCAGCTAGTTCATCATTGTTTATCAATGATAAAATATAATAAGCTAAAAGGAAACTATAAATA",
            "CGCTATTCATCATTGTTTATCAATGAGCTAAAAGGAAACTATAAATAACCATGTATAATTATAAGTAGGTACCTATTTTTTTATTTTAAACTGAAATTCAATATTATATAGGCAAAG"
        };

        size_t kmerSizes[] = { 13, 15, 17, 19, 21, 23, 25, 27, 29, 31};

        for (size_t i=0; i<ARRAY_SIZE(sequences); i++)
        {
            for (size_t j=0; j<ARRAY_SIZE(kmerSizes); j++)
            {
                /** We create the graph. */
                Graph graph = Graph::create (
                    new BankStrings (sequences[i], 0),
                    "-kmer-size %d  -abundance-min 1  -verbose 0  -max-memory %d", kmerSizes[j], MAX_MEMORY
                );

                debruijn_check_sequence (graph, kmerSizes[j], sequences[i]);

                /** We remove the graph. */
                graph.remove ();
            }
        }
    }

    /********************************************************************************/
    void debruijn_test4 ()
    {
        char* seq = (char*) "ACCATGTATAATTATAAGTAGGTACCT";  // size 27
        char* rev = (char*) "AGGTACCTACTTATAATTATACATGGT";

        /** We create the graph. */
        Graph graph = Graph::create (new BankStrings (seq, 0), "-kmer-size 27  -abundance-min 1  -verbose 0  -max-memory %d", MAX_MEMORY);

        Graph::Iterator<Node> it = graph.iterator<Node>();  it.first();

        Node n1 = it.item();
        CPPUNIT_ASSERT (n1.strand == STRAND_FORWARD);
        CPPUNIT_ASSERT (graph.toString(n1).compare (seq) == 0);

        Node n2 = graph.reverse (n1);
        CPPUNIT_ASSERT (n2.strand == STRAND_REVCOMP);
        CPPUNIT_ASSERT (graph.toString(n2).compare (rev) == 0);
    }

    /********************************************************************************/
    void debruijn_test5 ()
    {
        /** We create an empty graph with a given kmer size. */
        Graph graph = Graph::create (7);

        /** We define a string of size equal to the kmer size. */
        char* seq = (char*) "ACCAGTT";
        char* rev = (char*) "AACTGGT";

        Node n1 = graph.buildNode (Data(seq));
        CPPUNIT_ASSERT (graph.toString (n1) == seq);

        Node n2 = graph.reverse (n1);
        CPPUNIT_ASSERT (graph.toString (n2) == rev);
    }

    /********************************************************************************/
    struct debruijn_test6_fct
    {
        debruijn_test6_fct (const Graph& graph) : graph(graph) {}
        const Graph& graph;

        void operator() (Node& node) const
        {
            string snode = graph.toString (node);

            Node rev = graph.reverse (node);
            string srev = graph.toString (rev);

            /** We build a node from the reverse string. */
            Node rev2 = graph.buildNode (Data ((char*)srev.c_str()));
            CPPUNIT_ASSERT (graph.toString(rev2) == srev);

            /** We reverse the reversed node. */
            Node node2 = graph.reverse (rev2);
            CPPUNIT_ASSERT (graph.toString(node2) == snode);
        }
    };

    void debruijn_test6 ()
    {
        char* seq = (char*) "ACCATGTATAATTATAAGTAGGTACCACGATCGATCGATCGATCGTAGCATATCGTACGATCT";

        /** We create the graph. */
        Graph graph = Graph::create (new BankStrings (seq, 0), "-kmer-size 27  -abundance-min 1  -verbose 0  -max-memory %d", MAX_MEMORY);

        debruijn_test6_fct fct(graph);
        graph.iterator<Node>().iterate (fct);
    }

    /********************************************************************************/
    struct debruijn_test7_fct
    {
        debruijn_test7_fct (const Graph& graph, Node& n1, Node& n2) : graph(graph), n1(n1), n2(n2) {}
        const Graph& graph;
        Node& n1;
        Node& n2;

        void operator() (Node& current) const
        {
            string currentStr = graph.toString(current);

            CPPUNIT_ASSERT (currentStr==graph.toString(n1) || currentStr==graph.toString(n2) );

            /** We get all possible edges from the current kmer (note: not from the current node). */
            Graph::Vector<Edge> neighbors = graph.neighbors<Edge>(current.kmer);
            for (size_t i=0; i<neighbors.size(); i++)
            {
                /** Shortcut. */
                Edge& edge = neighbors[i];

                if (currentStr==graph.toString(n1))  // 1 neighbor
                {
                    CPPUNIT_ASSERT (neighbors.size()==1);

                    CPPUNIT_ASSERT (edge.nt==NUCL_C);
                    CPPUNIT_ASSERT (edge.direction==DIR_OUTCOMING);
                    CPPUNIT_ASSERT (graph.toString(edge.from)=="AGGCG");
                    CPPUNIT_ASSERT (graph.toString(edge.to)  =="GGCGC");
                }

                if (currentStr==graph.toString(n2))  // 2 neighbors
                {
                    CPPUNIT_ASSERT (neighbors.size()==2);

                    if (graph.toString(edge.from)=="GCGCC")
                    {
                        CPPUNIT_ASSERT (graph.toString(edge.to)=="CGCCT");
                        CPPUNIT_ASSERT (edge.nt==NUCL_T);
                        CPPUNIT_ASSERT (edge.direction==DIR_OUTCOMING);

                    }
                    else if (graph.toString(edge.from)=="GGCGC")
                    {
                        CPPUNIT_ASSERT (graph.toString(edge.to)=="GCGCC");
                        CPPUNIT_ASSERT (edge.nt==NUCL_C);
                        CPPUNIT_ASSERT (edge.direction==DIR_OUTCOMING);
                    }
                    else
                    {
                        CPPUNIT_ASSERT (false);
                    }
                }
            }
        }
    };

    void debruijn_test7 ()
    {
        /** We create the graph. */
        Graph graph = Graph::create (new BankStrings ("AGGCGC", 0),  "-kmer-size 5  -abundance-min 1  -verbose 0  -max-memory %d", MAX_MEMORY);

        /** We should get two kmers:
         *      - AGGCG / CGCCT
         *      - GCGCC / GGCGC
         */
        Node n1 = graph.buildNode ((char*)"AGGCG");
        Node n2 = graph.buildNode ((char*)"GCGCC");

        // We should get as neighborhood
        // GCGCC  [GCGCC --T--> CGCCT]
        // GCGCC  [GGCGC --C--> GCGCC]
        // AGGCG  [AGGCG --C--> GGCGC]

        debruijn_test7_fct fct (graph, n1, n2);
        graph.iterator<Node>().iterate (fct);
    }

    /********************************************************************************/
    void debruijn_test8_aux (char* seq, size_t kmerSize)
    {
        /** We create the graph. */
        Graph graph = Graph::create (new BankStrings (seq, NULL),  "-kmer-size %d  -abundance-min 1  -verbose 0  -max-memory %d", kmerSize, MAX_MEMORY);

        // We get the first node.
        Node node = graph.buildNode (seq);

        Graph::Iterator<Edge> path = graph.simplePath<Edge> (node, DIR_OUTCOMING);

        for (path.first(); !path.isDone(); path.next())
        {
            /** We check that the current transition nucleotide matches the correct character in the sequence. */
            CPPUNIT_ASSERT (ascii(path.item().nt) == seq [graph.getKmerSize() + path.rank()]);
        }

        /** We check that we found the correct number of nodes. */
        CPPUNIT_ASSERT (path.rank() == strlen (seq) - kmerSize);
    }

    /** */
    void debruijn_test8 ()
    {
        /** This sequence should not have branching nodes for kmer size big enough. */
        char* seq = (char*) "AGGCGCTAGGGTAGAGGATGATGA";

        size_t kmerSizes[] = {7, 9, 11, 13, 15, 17};

        for (size_t i=0; i<ARRAY_SIZE(kmerSizes); i++)  { debruijn_test8_aux (seq, kmerSizes[i]); }
    }

    /********************************************************************************/
    void debruijn_test9 ()
    {
        size_t kmerSize = 9;

        char* seq1 = (char*) "AGGCGCTAGGGTAGAGGATGATGA";
        char* seq2 = (char*) "AGGCGCTAGGGTATAGGATGATGA";
        //                    000000000011111111112222
        //                    012345678901234567890123
        //  difference here                ^

        /** We create the graph. */
        Graph graph = Graph::create (new BankStrings (seq1, seq2, NULL),  "-kmer-size %d  -abundance-min 1  -verbose 0  -max-memory %d ", kmerSize, MAX_MEMORY);

        /** We get the first node. */
        Node node = graph.buildNode (seq1);

        /** We get a simple path iterator starting from the beginning of the seq1. */
        Graph::Iterator<Edge> path = graph.simplePath<Edge> (node, DIR_OUTCOMING);

        for (path.first(); !path.isDone(); path.next())
        {
            CPPUNIT_ASSERT (graph.isSimple (path.item()));
            CPPUNIT_ASSERT (ascii(path.item().nt) == seq1 [graph.getKmerSize() + path.rank()]);
        }

        /** We check that we stopped at the first difference between the two sequences. */
        CPPUNIT_ASSERT (path.rank() == 4);   // 4 = diffOffset - kmerSize
    }

    /********************************************************************************/
    void debruijn_test10 ()
    {
        size_t kmerSize = 7;

        char* seq1 = (char*) "AGGCGCTAGGGAGAGGATGATGAAA";
        char* seq2 = (char*) "AGGCGCTAGGGTGAGGATGATGAAA";
        //  difference here              ^

        /** We create the graph. */
        Graph graph = Graph::create (new BankStrings (seq1, seq2, NULL),  "-kmer-size %d  -abundance-min 1  -verbose 0  -max-memory %d", kmerSize, MAX_MEMORY);

        /** We get the first node. */
        Node node = graph.buildNode (seq1);
        CPPUNIT_ASSERT (graph.toString(node) == "AGGCGCT");

        /** We retrieve the branching neighbors for the node. */
        Graph::Vector<BranchingNode> branchingNeighbors = graph.successors <BranchingNode> (node);

        /** In our example, we should have only one branching neighbor. */
        CPPUNIT_ASSERT (branchingNeighbors.size() == 1);

        /** We check this branching neighbor. */
        CPPUNIT_ASSERT (graph.toString(branchingNeighbors[0]) == "GCTAGGG");

        for (size_t i=0; i<branchingNeighbors.size(); i++)
        {
            /** We check the in/out degree to be sure that it is indeed a branching node. */
            CPPUNIT_ASSERT (graph.isBranching(branchingNeighbors[i]));
        }
    }

    /********************************************************************************/
    void debruijn_test11 ()
    {
        size_t kmerSize = 7;

        // We define some sequences used for building our test graph.
        // Note that the sequences have a difference at index==kmerSize
        const char* sequences[] =
        {
            //      x <- difference here
            "AGGCGCTAGGGAGAGGATGATGAAA",
            "AGGCGCTCGGGAGAGGATGATGAAA",
            "AGGCGCTTGGGAGAGGATGATGAAA"
        };

        // We create the graph.
        Graph graph = Graph::create (new BankStrings (sequences, ARRAY_SIZE(sequences)),  "-kmer-size %d  -abundance-min 1  -verbose 0  -max-memory %d", kmerSize, MAX_MEMORY);

        // We get the first node (should be AGGCGCT); this is a branching node.
        Node node = graph.buildNode ((char*)sequences[0]);
        CPPUNIT_ASSERT (graph.toString(node) == "AGGCGCT");

        /** We retrieve the branching neighbors for the node. */
        Graph::Vector<BranchingNode> branchingNeighbors = graph.successors <BranchingNode> (node);

        /** In our example, we should have 3 branching neighbors. */
        CPPUNIT_ASSERT (branchingNeighbors.size() == 3);

        for (size_t i=0; i<branchingNeighbors.size(); i++)
        {
            /** We should close the bubble, ie all the branching neighbors are identical. */
            CPPUNIT_ASSERT (graph.toString(branchingNeighbors[i]) == "GGGAGAG");

            /** We check the in/out degree to be sure that it is indeed a branching node. */
            CPPUNIT_ASSERT (graph.isBranching(branchingNeighbors[i]));
        }
    }

    /********************************************************************************/
    void debruijn_test12 ()
    {
        size_t kmerSize = 7;

        // We define some sequences used for building our test graph.
        // Note that the sequences have a difference at index==kmerSize
        const char* sequences[] =
        {
            //      x <- difference here
            "AGGCGCTAGGGAGAGGATGATGAAA",
            "AGGCGCTCGGGAGAGGATGATGAAA",
            "AGGCGCTTGGGAGAGGATGATGAAA"
        };

        // We create the graph.
        Graph graph = Graph::create (new BankStrings (sequences, ARRAY_SIZE(sequences)),  "-kmer-size %d  -abundance-min 1  -verbose 0  -max-memory %d", kmerSize, MAX_MEMORY);

        // We get the first node (should be AGGCGCT); this is a branching node.
        Node node = graph.buildNode ((char*)sequences[0]);
        CPPUNIT_ASSERT (graph.toString(node) == "AGGCGCT");

        /** We retrieve the branching neighbors for the node. */
        Graph::Vector<BranchingEdge> branchingNeighbors = graph.successors <BranchingEdge> (node);

        /** In our example, we should have 3 branching neighbors. */
        CPPUNIT_ASSERT (branchingNeighbors.size() == 3);

        for (size_t i=0; i<branchingNeighbors.size(); i++)
        {
            /** We should close the bubble, ie all the branching neighbors are identical. */
            CPPUNIT_ASSERT (graph.toString(branchingNeighbors[i].to) == "GGGAGAG");

            /** We check the path length. */
            CPPUNIT_ASSERT (branchingNeighbors[i].distance == 7);

            /** We check the simple path between the two branching nodes.
             *  We need first to retrieve the first (simple) neighbor from the origin. */
            Node simpleNeighbor = graph.successor<Node> (branchingNeighbors[i].from, branchingNeighbors[i].nt);

            Graph::Iterator<Edge> path = graph.simplePath<Edge> (simpleNeighbor, branchingNeighbors[i].direction);
            for (path.first(); !path.isDone(); path.next())
            {
                CPPUNIT_ASSERT (graph.toString(path->from) == string (sequences[i], path.rank()+1, kmerSize));
                CPPUNIT_ASSERT (graph.isSimple(*path));
            }
        }
    }

    /********************************************************************************/
    void debruijn_test13 ()
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
        Graph graph = Graph::create (bank,  "-kmer-size %d  -abundance-min %d  -verbose 0  -max-memory %d", kmerSize, nks, MAX_MEMORY);

        // We check we got the correct number of solid kmers.
        CPPUNIT_ASSERT (graph.getInfo().getInt ("kmers_nb_solid") == strlen(seq) - kmerSize + 1);

        // We check that we have only two branching nodes.
        CPPUNIT_ASSERT (graph.getInfo().getInt ("nb_branching") == 2);
    }

    /********************************************************************************/
    void debruijn_mutation_aux (const char* sequences[], size_t len, bool reverse)
    {
        size_t kmerSize = strlen (sequences[0]);

        // We create the graph.
        Graph graph = Graph::create (new BankStrings (sequences, len),  "-kmer-size %d  -abundance-min 1  -verbose 0  -max-memory %d", kmerSize, MAX_MEMORY);

        // We get the first node (should be AGGCGCT); this is a branching node.
        Node node = graph.buildNode ((char*)sequences[0]);
        CPPUNIT_ASSERT (graph.toString(node) == sequences[0]);

        if (reverse)  { node = graph.reverse(node); }

        Graph::Vector<Node> mutations = graph.mutate (node, kmerSize-1, 1);
        CPPUNIT_ASSERT (mutations.size() == 3);

        // We check we got the correct mutations
        for (size_t i=0; i<mutations.size(); i++)
        {
            Node item = mutations[i];

            if (reverse)  { item = graph.reverse (item); }

            CPPUNIT_ASSERT (graph.toString(item) == sequences[i+1]);
        }
    }

    /** */
    void debruijn_mutation ()
    {
        const char* sequences1[] =
        {
            "TTGCTCACATGTTCTTTCCTGCGTTATCCCA",
            "TTGCTCACATGTTCTTTCCTGCGTTATCCCC",
            "TTGCTCACATGTTCTTTCCTGCGTTATCCCT",
            "TTGCTCACATGTTCTTTCCTGCGTTATCCCG"

        };

        const char* sequences2[] =
        {
            "TTTGCTCACATGTTCTTTCCTGCGTTATCCC",
            "GTTGCTCACATGTTCTTTCCTGCGTTATCCC",
            "ATTGCTCACATGTTCTTTCCTGCGTTATCCC",
            "CTTGCTCACATGTTCTTTCCTGCGTTATCCC"
        };

        debruijn_mutation_aux (sequences1, ARRAY_SIZE(sequences1), false);
        debruijn_mutation_aux (sequences2, ARRAY_SIZE(sequences2), true);
    }

    /********************************************************************************/
    void debruijn_mphf_aux (const char* sequences[], size_t len, const int abundances[])    
    {
        size_t kmerSize = strlen (sequences[0]);

        // We create the graph.
        Graph graph = Graph::create (new BankStrings (sequences, len),  "-kmer-size %d  -abundance-min 1  -verbose 0 -mphf emphf  -max-memory %d", kmerSize, MAX_MEMORY);

        Graph::Iterator<Node> it = graph.iterator<Node> ();

        // debugging
        for (it.first(); !it.isDone(); it.next())
        {
            //std::cout << graph.toString (it.item()) << " test printing node abundance " << it.item().abundance << std::endl;
        }

        for (int i = len-1; i >=  0; i--)
        {
            // random access to nodes
            Node node = graph.buildNode ((char*)sequences[i]);
            CPPUNIT_ASSERT (graph.toString(node) == sequences[i]);
            int abundance = graph.queryAbundance(node);
            //std::cout << graph.toString(node) << " test printing node abundance " << abundance << " expected abundance:" << abundances[i] << std::endl;
            CPPUNIT_ASSERT (abundance == abundances[i]);
        }
    }

    /** */
    void debruijn_mphf ()
    {
        const char* sequences1[] =
        {
            "TTGCTCACATGTTCTTTCCTGCGTTATCCCA",
            "TTGCTCACATGTTCTTTCCTGCGTTATCCCC",
            "TTGCTCACATGTTCTTTCCTGCGTTATCCCC",
            "TTGCTCACATGTTCTTTCCTGCGTTATCCCT",
            "TTGCTCACATGTTCTTTCCTGCGTTATCCCT",
            "TTGCTCACATGTTCTTTCCTGCGTTATCCCT",
            "TTGCTCACATGTTCTTTCCTGCGTTATCCCG",
            "TTGCTCACATGTTCTTTCCTGCGTTATCCCG",
            "TTGCTCACATGTTCTTTCCTGCGTTATCCCG",
            "TTGCTCACATGTTCTTTCCTGCGTTATCCCG"
        };

        const int abundances[] = { 1, 2, 2, 3, 3, 3, 4, 4, 4, 4 };

        debruijn_mphf_aux (sequences1, ARRAY_SIZE(sequences1), abundances);
    }

    /********************************************************************************/
    struct debruijn_build_entry
    {
        debruijn_build_entry () : nbNodes(0), nbBranchingNodes(0) {}
        size_t  nbNodes;
        Integer checksumNodes;
        size_t  nbBranchingNodes;
        Integer checksumBranchingNodes;
    };

    debruijn_build_entry debruijn_build_aux_aux (const char* name, bool checkNodes, bool checkBranching)
    {
        debruijn_build_entry result;

        /** We load the graph. */
        Graph graph = Graph::load (name);

        if (checkNodes)
        {
            Graph::Iterator<Node> iterNodes = graph.iterator<Node>();
            for (iterNodes.first(); !iterNodes.isDone(); iterNodes.next())
            { result.nbNodes++; result.checksumNodes += iterNodes.item().kmer; }
        }

        if (checkBranching)
        {
            Graph::Iterator<BranchingNode> iterBranchingNodes = graph.iterator<BranchingNode>();
            for (iterBranchingNodes.first(); !iterBranchingNodes.isDone(); iterBranchingNodes.next())
            { result.nbBranchingNodes++; result.checksumBranchingNodes += iterBranchingNodes.item().kmer; }
        }

        return result;
    }

    /********************************************************************************/
    void debruijn_build_aux (const char* sequences[], size_t nbSequences)
    {
        // We build the bank
        IBank* inputBank = new BankStrings (sequences, nbSequences);
        LOCAL (inputBank);

        Graph::create (inputBank,  "-kmer-size 31 -out %s -abundance-min 1  -verbose 0  -max-memory %d",                        "g1", MAX_MEMORY);
        Graph::create (inputBank,  "-kmer-size 31 -out %s -abundance-min 1  -verbose 0 -branching-nodes none  -max-memory %d",  "g2", MAX_MEMORY);
        Graph::create (inputBank,  "-kmer-size 31 -out %s -abundance-min 1  -verbose 0 -solid-kmers-out none  -max-memory %d",  "g3", MAX_MEMORY);

        debruijn_build_entry r1 = debruijn_build_aux_aux ("g1", true,  true);
        debruijn_build_entry r2 = debruijn_build_aux_aux ("g2", true,  true);
        debruijn_build_entry r3 = debruijn_build_aux_aux ("g3", false, true);

        CPPUNIT_ASSERT (r1.nbNodes       == r2.nbNodes);
        CPPUNIT_ASSERT (r1.checksumNodes == r2.checksumNodes);

        /** Right now, we don't test Node for g3 since the Node iterator (based only on BranchingNode) is not implemented. */
        //CPPUNIT_ASSERT (r1.nbNodes       == r3.nbNodes);
        //CPPUNIT_ASSERT (r1.checksumNodes == r3.checksumNodes);

        CPPUNIT_ASSERT (r1.nbBranchingNodes       == r2.nbBranchingNodes);
        CPPUNIT_ASSERT (r1.checksumBranchingNodes == r2.checksumBranchingNodes);
        CPPUNIT_ASSERT (r1.nbBranchingNodes       == r3.nbBranchingNodes);
        CPPUNIT_ASSERT (r1.checksumBranchingNodes == r3.checksumBranchingNodes);
    }

    /********************************************************************************/
    void debruijn_build ()
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

        debruijn_build_aux (sequences, ARRAY_SIZE(sequences));
    }

    /********************************************************************************/
    void debruijn_checksum_aux2 (
        const string& readfile,
        size_t kmerSize,
        size_t integerPrecision,
        size_t nbBranching,
        const string& checksum,
        const string& bloom,
        const string& debloom,
        const string& debloomImpl
    )
    {
        size_t maxMem[]  = { 500, MAX_MEMORY };
        size_t nbCores[] = { 1, 4 };

        for (size_t i=0; i<ARRAY_SIZE(maxMem); i++)
        {
            for (size_t j=0; j<ARRAY_SIZE(nbCores); j++)
            {
                //printf ("file=%s k=%d  prec=%d  bloom=%s  debloom=%s  debloomImpl=%s  mem=%d  nbCores=%d\n",
                //    readfile.c_str(), kmerSize, integerPrecision, bloom.c_str(), debloom.c_str(), debloomImpl.c_str(), maxMem[i], nbCores[j]
                //);

                Graph graph = Graph::create (
                    "-verbose 0 -in %s -kmer-size %d -integer-precision %d  -bloom %s  -debloom %s  -debloom-impl %s  -max-memory %d  -nb-cores %d  -max-memory %d",
                    readfile.c_str(),
                    kmerSize,
                    integerPrecision,
                    bloom.c_str(),
                    debloom.c_str(),
                    debloomImpl.c_str(),
                    maxMem[i],
                    nbCores[j],
                    MAX_MEMORY
                );
                CPPUNIT_ASSERT (graph.getInfo().getStr ("checksum_branching") == checksum);
                CPPUNIT_ASSERT (graph.getInfo().getInt ("nb_branching")       == nbBranching);
            }
        }
    }

    void debruijn_checksum_aux (
        const string readfile,
        size_t kmerSize,
        size_t integerPrecision,
        size_t nbBranching,
        const string& checksum,
        bool all = false
    )
    {
        if (all)
        {
            debruijn_checksum_aux2 (readfile, kmerSize, integerPrecision, nbBranching, checksum, "basic",    "original",  "basic");
            debruijn_checksum_aux2 (readfile, kmerSize, integerPrecision, nbBranching, checksum, "cache",    "original",  "basic");
            debruijn_checksum_aux2 (readfile, kmerSize, integerPrecision, nbBranching, checksum, "neighbor", "original",  "basic");

            debruijn_checksum_aux2 (readfile, kmerSize, integerPrecision, nbBranching, checksum, "basic",    "cascading", "basic");
            debruijn_checksum_aux2 (readfile, kmerSize, integerPrecision, nbBranching, checksum, "cache",    "cascading", "basic");
            debruijn_checksum_aux2 (readfile, kmerSize, integerPrecision, nbBranching, checksum, "neighbor", "cascading", "basic");
        }

        debruijn_checksum_aux2 (readfile, kmerSize, integerPrecision, nbBranching, checksum, "neighbor", "original",  "minimizer");
        debruijn_checksum_aux2 (readfile, kmerSize, integerPrecision, nbBranching, checksum, "neighbor", "cascading", "minimizer");
    }

    void debruijn_checksum ()
    {
        /** The nbBranching and checksum values have been computed with the SVN version of minia (minia-1.6906).
         *  Now we compare to the values produced by gatb-core. */

        string filepath;

        /*****************************************/
        /**          FILE reads1.fa              */
        /*****************************************/
        filepath = DBPATH("reads1.fa");

        debruijn_checksum_aux (filepath,  31,  1, 24,   "30eb72bc69eca0d3", true);
        debruijn_checksum_aux (filepath,  31,  2, 24, "2.30eb72bc69eca0d3", true);
        debruijn_checksum_aux (filepath,  31,  3, 24, "2.30eb72bc69eca0d3", true);
        debruijn_checksum_aux (filepath,  31,  4, 24, "2.30eb72bc69eca0d3", true);

        debruijn_checksum_aux (filepath,  63,  2,  8, "92acb8443ed65990.7b4298b762ce39ff", true);
        debruijn_checksum_aux (filepath,  63,  3,  8, "92acb8443ed65990.7b4298b762ce39ff", true);
        debruijn_checksum_aux (filepath,  63,  4,  8, "92acb8443ed65990.7b4298b762ce39ff", true);

        debruijn_checksum_aux (filepath,  95,  3,  4, "71ed998e1b26a8e0.1a1d73f05438c413.1c6405c67a8fab0a", true);
        debruijn_checksum_aux (filepath,  95,  4,  4, "71ed998e1b26a8e0.1a1d73f05438c413.1c6405c67a8fab0a", true);

        debruijn_checksum_aux (filepath, 127,  4,  4, "5a5d5720302692fd.214182472e05744f.5c6b807cecb99db2.c5655b04b6a7b8fe", true);

        /*****************************************/
        /**          FILE reads3.fa.gz           */
        /*****************************************/
        filepath = DBPATH("reads3.fa.gz");

        debruijn_checksum_aux (filepath,  31,  1, 2956,    "d238698aa54e0ce2");
        debruijn_checksum_aux (filepath,  31,  2, 2956, "cc.d238698aa54e0ce2");
        debruijn_checksum_aux (filepath,  31,  3, 2956, "cc.d238698aa54e0ce2");
        debruijn_checksum_aux (filepath,  31,  4, 2956, "cc.d238698aa54e0ce2");

        debruijn_checksum_aux (filepath,  63,  2,  969,    "f0a5da085cc20ee8.90af6e8e523e8b96");
        debruijn_checksum_aux (filepath,  63,  3,  969, "40.f0a5da085cc20ee8.90af6e8e523e8b96");
        debruijn_checksum_aux (filepath,  63,  4,  969, "40.f0a5da085cc20ee8.90af6e8e523e8b96");

        debruijn_checksum_aux (filepath,  95,  3,  600,    "99817bec5ebde83b.97a71a71c72636c7.4cd2353be480a3b4");
        debruijn_checksum_aux (filepath,  95,  4,  600, "2c.99817bec5ebde83b.97a71a71c72636c7.4cd2353be480a3b4");

        debruijn_checksum_aux (filepath, 127,  4,  424, "50c7d59f28890ef3.23f38c0611dc341c.525fd5f6a6fa045b.6f1255d3695f039d");
    }

    /********************************************************************************/

    /** Check that the branching are sorted. */
    void debruijn_checkbranching ()
    {
        string filepath = DBPATH("reads3.fa.gz");

        /** We create a graph. */
        Graph graph = Graph::create ("-verbose 0 -in %s -max-memory %d", filepath.c_str(), MAX_MEMORY);

        Node::Value previous = 0;
        size_t i=0;

        /** We iterate the branching nodes. */
        Graph::Iterator<BranchingNode> it = graph.iterator<BranchingNode> ();
        for (it.first(); !it.isDone(); it.next(), i++)
        {
            if (i > 0) {  CPPUNIT_ASSERT (previous < it->kmer);  }
            previous = it->kmer;
        }
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestDebruijn);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestDebruijn);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

