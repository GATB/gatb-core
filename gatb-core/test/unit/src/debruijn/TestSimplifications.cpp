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
#include <gatb/debruijn/impl/Simplifications.hpp>
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

#define DEBUGprint(a) a

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for genomic databases management
 */
class TestSimplifications : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestSimplifications);
        CPPUNIT_TEST_GATB (debruijn_simpl_X);
        CPPUNIT_TEST_GATB (debruijn_simpl_tip);
        CPPUNIT_TEST_GATB (debruijn_simpl_bubble);
        CPPUNIT_TEST_GATB (debruijn_simpl_ec);
    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {}
    void tearDown ()  {}

    // SMALL VALUE NEEDED because continuous integration servers are not very powerful...
    static const u_int64_t MAX_MEMORY = 1000;

   /********************************************************************************/
    struct debruijn_build_entry
    {
        debruijn_build_entry () : nbNodes(0), nbNonDeletedNodes(0), nbBranchingNodes(0) {}
        size_t  nbNodes;
        size_t  nbNonDeletedNodes;
        Integer checksumNodes;
        size_t  nbBranchingNodes;
        Integer checksumBranchingNodes;
    };

    debruijn_build_entry debruijn_stats (Graph& graph, bool checkNodes, bool checkBranching)
    {
        debruijn_build_entry result;

        result.nbNodes = 0; result.nbNonDeletedNodes = 0; result.nbBranchingNodes = 0;

        if (checkNodes)
        {
            GraphIterator<Node> iterNodes = graph.iterator();
            for (iterNodes.first(); !iterNodes.isDone(); iterNodes.next())
            { result.nbNodes++; result.checksumNodes += iterNodes.item().kmer; 
            
              if (! graph.isNodeDeleted(*iterNodes)) 
              {
                    result.nbNonDeletedNodes++;
                    //std::cout << "non deleted node: " << graph.toString(*iterNodes) << std::endl;
              }
            }
        }

        if (checkBranching)
        {
            GraphIterator<BranchingNode> iterBranchingNodes = graph.iteratorBranching();
            for (iterBranchingNodes.first(); !iterBranchingNodes.isDone(); iterBranchingNodes.next())
            { result.nbBranchingNodes++; result.checksumBranchingNodes += iterBranchingNodes.item().kmer; }
        }

        return result;
    }

    /********************************************************************************/
    void debruijn_build (const char* sequences[], size_t nbSequences, int k, Graph& graph, debruijn_build_entry& r)
    {
        // We build the bank
        graph = Graph::create (new BankStrings (sequences, nbSequences),  "-kmer-size %d  -abundance-min 1  -verbose 0  -max-memory %d", k , MAX_MEMORY);

        r = debruijn_stats (graph, true,  true);
    }

    void debruijn_traversal (Graph& graph, string startSeq ,const char* checkStr)
    {
        TraversalKind traversalKind=TRAVERSAL_UNITIG; /* setting traversal=CONTIG doesn't make sense anymore, it was for miniav1 */

        // We create a Terminator object
        BranchingTerminator terminator (graph);

		// We create a node from the start of the first sequence
        string startKmer = startSeq.substr(0, graph._kmerSize);
		Node node = graph.buildNode (startKmer.c_str());

        Path path;

        // We create a Traversal instance according to the chosen traversal kind
        Traversal* traversal = Traversal::create (traversalKind, graph, terminator);
        LOCAL (traversal);

        traversal->traverse (node, DIR_OUTCOMING, path);

        stringstream ss;  ss << graph.toString (node) << path ;
        //DEBUGprint(std::cout << "result of traversal: " << ss.str() << std::endl;)
        CPPUNIT_ASSERT (ss.str().compare(checkStr)==0);
    }

    /********************************************************************************/

    /********************************************************************************/
    void debruijn_simpl_X ()
    {
        size_t kmerSize = 5;

        const char* sequences[] =
        {
            // it's a classical X for k=5 (node-centric): 1 and 2 go to 3, 3 go to 4 and 5;
            "AAAAA",
            "CCCCAAG",
            "AAACAAG",
            "CAAGA",
            "AAGAAGC",
            "AAGACCC"
        };

        // We create the graph.
        Graph graph;
        debruijn_build_entry r;
        debruijn_build(sequences, ARRAY_SIZE(sequences), kmerSize, graph, r);

        CPPUNIT_ASSERT (r.nbNodes == 14);
        if (r.nbNonDeletedNodes != 14)
            std::cout << "ancipating failed assert: nbNonDeletedNodes = " << r.nbNonDeletedNodes << " != 14" << std::endl;
        CPPUNIT_ASSERT (r.nbNonDeletedNodes == 14);

        // simplify it
        graph.simplify(1, false); // one core, no verbose

        r = debruijn_stats (graph, true,  true);
        //std::cout << "nb nodes:" << r.nbNodes << std::endl;
        CPPUNIT_ASSERT (r.nbNodes == 14);
        if (r.nbNonDeletedNodes != 2)
            std::cout << "ancipating failed assert: nbNonDeletedNodes = " << r.nbNonDeletedNodes << " != 14" << std::endl;

        CPPUNIT_ASSERT (r.nbNonDeletedNodes == 2); // only two nodes remaining: AAAAA and CAAGA
        
   }

    void debruijn_simpl_tip ()
    {
        size_t kmerSize = 21;

        const char* sequences[] =
        {
            //>works well for k=21; part of genome10K.fasta
            "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACC",
            "TGTCATCTAGTTCAACAACCAAAAAAA", //>that's the tip
            "TGTCATCTAGTTCAACAACCGTTATGCCGTCCGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACATG" //>remaining part
        };

        // We create the graph.
        Graph graph;
        debruijn_build_entry r;
        debruijn_build(sequences, ARRAY_SIZE(sequences), kmerSize, graph, r);

        //std::cout << "nb nodes:" << r.nbNodes << std::endl;
        CPPUNIT_ASSERT (r.nbNodes == 624);
        CPPUNIT_ASSERT (r.nbNonDeletedNodes == 624);

        // simplify it
        graph.simplify(1, false); // one core, no verbose


        // how many nodes left? should be as many as initially. it's a negative test: graph shouldn't be simplified
        r = debruijn_stats (graph, true,  true);
        //std::cout << "nb nodes post simp:" << r.nbNodes << std::endl;
        //std::cout << "nb non deleted nodes post simp:" << r.nbNonDeletedNodes << std::endl;
        CPPUNIT_ASSERT (r.nbNodes == 624);
        CPPUNIT_ASSERT (r.nbNonDeletedNodes == 617);

        debruijn_traversal (graph, sequences[0], "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACCGTTATGCCGTCCGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACATG");
    }


    /********************************************************************************/


    void debruijn_simpl_bubble ()
    {
        size_t kmerSize = 21;

        const char* sequences[] =
        {
            "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACC", //>works well for k=21; part of genome10K.fasta
            "TGTCATCTAGTTCAACAACCAAAATAACGACTCTTGCGCTCGGATGT", //>that's the bubble  (highly covered)
            "TGTCATCTAGTTCAACAACCAAAATAACGACTCTTGCGCTCGGATGT", //>that's the bubble
            "TGTCATCTAGTTCAACAACCAAAATAACGACTCTTGCGCTCGGATGT", //>that's the bubble
            "TGTCATCTAGTTCAACAACCAAAAAAACGACTCTTGCGCTCGGATGT", //>that's the bubble path 2, low covered
            "CGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACAAA" // >remaining part
        };

        // We create the graph.
        Graph graph;
        debruijn_build_entry r;
        debruijn_build(sequences, ARRAY_SIZE(sequences), kmerSize, graph, r);

        //DEBUGprint(std::cout << "nb nodes:" << r.nbNodes << std::endl;)
        CPPUNIT_ASSERT (r.nbNodes == 634);
        CPPUNIT_ASSERT (r.nbNonDeletedNodes == 634);

        // simplify it
        graph.simplify(1, false); // one core, no verbose


        // how many nodes left? should be as many as initially. it's a negative test: graph shouldn't be simplified
        r = debruijn_stats (graph, true,  true);
        //DEBUGprint(std::cout << "nb nodes post simp:" << r.nbNodes << std::endl;)
        //DEBUGprint(std::cout << "nb non deleted nodes post simp:" << r.nbNonDeletedNodes << std::endl;)
        CPPUNIT_ASSERT (r.nbNodes == 634);
        CPPUNIT_ASSERT (r.nbNonDeletedNodes == 613);

        debruijn_traversal (graph, sequences[0], "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACCAAAATAACGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACAAA");
    }


    void debruijn_simpl_ec ()
    {
        size_t kmerSize = 21;

        const char* sequences[] =
        {
            // from ec.fa in minia/tests
            "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACC", //>works well for k=21; part of genome10K.fasta
            "TGTCATCTAGTTCAACAACCGTTATGCCGTCCGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACATG", // >remaining part
            "TGTCATCTAGTTCAACAACCAAAAAAA", //>that's the EC
            "GGTGAACAGCACATCTTTTCGTCCTGAGGCCATATTAATTCTACTCAGATTGTCTGTAACCGGAGCTTCGGGCGTATTTTTGCGTAAGACACTGCCTAAAGGGAACATATGTGTCCAGAATAGGGTTCAACGGTGTATGAGCAAACTAGTTCAACAACCAAAAAAATTGTGTGCAAGCTACTTCTAGACCTTATTAAGTGCCCAGGAATTCCTAGGAAGGCGCGCAGCTCAAGCAATCATACATGGCGGAATGCCTGTCCACCGGGGGTTCTACTGTACCACAGTGGCCTGGATAGCTAAGCAGGTCCTGGATTGGCATGTCATCCGGAGTGATAGGCACTGCTCACGACCAGCTTGCGGACAAACGGGGTGCCCGCGCCTGCGTCCGGTAGACGAGCGATGGATTTAGACCGTTCACTGAACCCTCTAATAGGACCTCTTGCCCATCCGAGGCTTAAGC", /* >contig that is split in two by the EC, containing the last kmer of the ec (CTAGTTCAACAACCAAAAAAA) ; ---now let's just repeat some of these sequences for coverage reasons
*/
            "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACC", //>works well for k=21; part of genome10K.fasta
            "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACC", //>works well for k=21; part of genome10K.fasta
            "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACC", //>works well for k=21; part of genome10K.fasta
            "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACC", //>works well for k=21; part of genome10K.fasta
            "TGTCATCTAGTTCAACAACCGTTATGCCGTCCGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACATG", // >remaining part
            "TGTCATCTAGTTCAACAACCGTTATGCCGTCCGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACATG", // >remaining part
            "TGTCATCTAGTTCAACAACCGTTATGCCGTCCGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACATG", // >remaining part
            "TGTCATCTAGTTCAACAACCGTTATGCCGTCCGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACATG", // >remaining part
    
            "GGTGAACAGCACATCTTTTCGTCCTGAGGCCATATTAATTCTACTCAGATTGTCTGTAACCGGAGCTTCGGGCGTATTTTTGCGTAAGACACTGCCTAAAGGGAACATATGTGTCCAGAATAGGGTTCAACGGTGTATGAGCAAACTAGTTCAACAACCAAAAAAATTGTGTGCAAGCTACTTCTAGACCTTATTAAGTGCCCAGGAATTCCTAGGAAGGCGCGCAGCTCAAGCAATCATACATGGCGGAATGCCTGTCCACCGGGGGTTCTACTGTACCACAGTGGCCTGGATAGCTAAGCAGGTCCTGGATTGGCATGTCATCCGGAGTGATAGGCACTGCTCACGACCAGCTTGCGGACAAACGGGGTGCCCGCGCCTGCGTCCGGTAGACGAGCGATGGATTTAGACCGTTCACTGAACCCTCTAATAGGACCTCTTGCCCATCCGAGGCTTAAGC", // >contig that is split in two by the EC, containing the last kmer of the ec
            "GGTGAACAGCACATCTTTTCGTCCTGAGGCCATATTAATTCTACTCAGATTGTCTGTAACCGGAGCTTCGGGCGTATTTTTGCGTAAGACACTGCCTAAAGGGAACATATGTGTCCAGAATAGGGTTCAACGGTGTATGAGCAAACTAGTTCAACAACCAAAAAAATTGTGTGCAAGCTACTTCTAGACCTTATTAAGTGCCCAGGAATTCCTAGGAAGGCGCGCAGCTCAAGCAATCATACATGGCGGAATGCCTGTCCACCGGGGGTTCTACTGTACCACAGTGGCCTGGATAGCTAAGCAGGTCCTGGATTGGCATGTCATCCGGAGTGATAGGCACTGCTCACGACCAGCTTGCGGACAAACGGGGTGCCCGCGCCTGCGTCCGGTAGACGAGCGATGGATTTAGACCGTTCACTGAACCCTCTAATAGGACCTCTTGCCCATCCGAGGCTTAAGC", // >contig that is split in two by the EC, containing the last kmer of the ec
            "GGTGAACAGCACATCTTTTCGTCCTGAGGCCATATTAATTCTACTCAGATTGTCTGTAACCGGAGCTTCGGGCGTATTTTTGCGTAAGACACTGCCTAAAGGGAACATATGTGTCCAGAATAGGGTTCAACGGTGTATGAGCAAACTAGTTCAACAACCAAAAAAATTGTGTGCAAGCTACTTCTAGACCTTATTAAGTGCCCAGGAATTCCTAGGAAGGCGCGCAGCTCAAGCAATCATACATGGCGGAATGCCTGTCCACCGGGGGTTCTACTGTACCACAGTGGCCTGGATAGCTAAGCAGGTCCTGGATTGGCATGTCATCCGGAGTGATAGGCACTGCTCACGACCAGCTTGCGGACAAACGGGGTGCCCGCGCCTGCGTCCGGTAGACGAGCGATGGATTTAGACCGTTCACTGAACCCTCTAATAGGACCTCTTGCCCATCCGAGGCTTAAGC", // >contig that is split in two by the EC, containing the last kmer of the ec
            "GGTGAACAGCACATCTTTTCGTCCTGAGGCCATATTAATTCTACTCAGATTGTCTGTAACCGGAGCTTCGGGCGTATTTTTGCGTAAGACACTGCCTAAAGGGAACATATGTGTCCAGAATAGGGTTCAACGGTGTATGAGCAAACTAGTTCAACAACCAAAAAAATTGTGTGCAAGCTACTTCTAGACCTTATTAAGTGCCCAGGAATTCCTAGGAAGGCGCGCAGCTCAAGCAATCATACATGGCGGAATGCCTGTCCACCGGGGGTTCTACTGTACCACAGTGGCCTGGATAGCTAAGCAGGTCCTGGATTGGCATGTCATCCGGAGTGATAGGCACTGCTCACGACCAGCTTGCGGACAAACGGGGTGCCCGCGCCTGCGTCCGGTAGACGAGCGATGGATTTAGACCGTTCACTGAACCCTCTAATAGGACCTCTTGCCCATCCGAGGCTTAAGC" // >contig that is split in two by the EC, containing the last kmer of the ec

        };

        // We create the graph.
        Graph graph;
        debruijn_build_entry r;
        debruijn_build(sequences, ARRAY_SIZE(sequences), kmerSize, graph, r);

        //DEBUGprint(std::cout << "nb nodes:" << r.nbNodes << std::endl;)
        CPPUNIT_ASSERT (r.nbNodes == 1063);
        CPPUNIT_ASSERT (r.nbNonDeletedNodes == 1063);

        // simplify it
        graph.simplify(1, false); // one core, no verbose


        // how many nodes left? should be as many as initially. it's a negative test: graph shouldn't be simplified
        r = debruijn_stats (graph, true,  true);
        //DEBUGprint(std::cout << "nb nodes post simp:" << r.nbNodes << std::endl;)
        //DEBUGprint(std::cout << "nb non deleted nodes post simp:" << r.nbNonDeletedNodes << std::endl;)
        CPPUNIT_ASSERT (r.nbNodes == 1063);
        if (r.nbNonDeletedNodes != 1057)
            std::cout << " in anticipation of assert fail, r.nbNonDeletedNodes = " << r.nbNonDeletedNodes << std::endl;
        CPPUNIT_ASSERT (r.nbNonDeletedNodes == 1057);

        debruijn_traversal (graph, sequences[0], "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACCGTTATGCCGTCCGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACATG");
        debruijn_traversal (graph, sequences[3], "GGTGAACAGCACATCTTTTCGTCCTGAGGCCATATTAATTCTACTCAGATTGTCTGTAACCGGAGCTTCGGGCGTATTTTTGCGTAAGACACTGCCTAAAGGGAACATATGTGTCCAGAATAGGGTTCAACGGTGTATGAGCAAACTAGTTCAACAACCAAAAAAATTGTGTGCAAGCTACTTCTAGACCTTATTAAGTGCCCAGGAATTCCTAGGAAGGCGCGCAGCTCAAGCAATCATACATGGCGGAATGCCTGTCCACCGGGGGTTCTACTGTACCACAGTGGCCTGGATAGCTAAGCAGGTCCTGGATTGGCATGTCATCCGGAGTGATAGGCACTGCTCACGACCAGCTTGCGGACAAACGGGGTGCCCGCGCCTGCGTCCGGTAGACGAGCGATGGATTTAGACCGTTCACTGAACCCTCTAATAGGACCTCTTGCCCATCCGAGGCTTAAGC");
    }

};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestSimplifications);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestSimplifications);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

