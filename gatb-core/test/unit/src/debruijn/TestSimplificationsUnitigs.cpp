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

#include <gatb/debruijn/impl/GraphUnitigs.hpp>
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

#define DEBUGprint(a) //a

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for genomic databases management
 */
class TestSimplificationsUnitigs : public Test
{
    typedef gatb::core::debruijn::impl::GraphUnitigsTemplate<32> GraphUnitigs;

    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestSimplificationsUnitigs);
        CPPUNIT_TEST_GATB (debruijn_simplunitigs_ec);
        CPPUNIT_TEST_GATB (debruijn_simplunitigs_X);
        CPPUNIT_TEST_GATB (debruijn_simplunitigs_tip);
        CPPUNIT_TEST_GATB (debruijn_simplunitigs_bubble);
        CPPUNIT_TEST_GATB (debruijn_simplunitigs_bubble_snp);
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

    debruijn_build_entry debruijn_stats (GraphUnitigs& graph, bool checkNodes, bool checkBranching)
    {
        debruijn_build_entry result;

        result.nbNodes = 0; result.nbNonDeletedNodes = 0; result.nbBranchingNodes = 0;

        if (checkNodes)
        {
            GraphIterator<NodeGU> iterNodes = graph.iterator();
            for (iterNodes.first(); !iterNodes.isDone(); iterNodes.next())
            { result.nbNodes++; /*result.checksumNodes += iterNodes.item().kmer; */
            
              if (! graph.isNodeDeleted(*iterNodes)) 
              {
                    result.nbNonDeletedNodes++;
                    //std::cout << "non deleted node: " << graph.toString(*iterNodes) << std::endl;
              }
            }
        }

        return result;
    }

    /********************************************************************************/
    void debruijn_build (const char* sequences[], size_t nbSequences, int k, GraphUnitigs& graph, debruijn_build_entry& r)
    {
        // We build the bank
        graph = GraphUnitigs::create (new BankStrings (sequences, nbSequences),  "-kmer-size %d  -abundance-min 1  -verbose 0  -max-memory %d -nb-cores 1 -minimizer-size 3" /*because some of kmer sizes are 5*/, k , MAX_MEMORY);

        r = debruijn_stats (graph, true,  true);
    }

    static 
        char revcomp (char s) {
            if (s == 'A') return 'T';
            else if (s == 'C') return 'G';
            else if (s == 'G') return 'C';
            else if (s == 'T') return 'A';
            else if (s == 'a') return 't';
            else if (s == 'c') return 'g';
            else if (s == 'g') return 'c';
            else if (s == 't') return 'a';
            return 'X';
        }

    static string revcomp (const string &s) {
        string rc;
        for (signed int i = s.length() - 1; i >= 0; i--) {rc += revcomp(s[i]);}
        return rc;
    }

    void debruijn_traversal (GraphUnitigs& graph, string startSeq ,const char* checkStr, const char *checkStr2 = nullptr)
    {

        string startKmer = startSeq.substr(0, graph._kmerSize);
		NodeGU node = graph.debugBuildNode(startKmer);

        bool isolatedLeft, isolatedRight;
        float coverage = 0;
        string sequence = graph.simplePathBothDirections(node, isolatedLeft, isolatedRight, true, coverage);
        string rev_seq = revcomp(sequence);

        if ((sequence.compare(checkStr) != 0 && rev_seq.compare(checkStr) != 0 && checkStr2 == nullptr) || 
                (checkStr2 != nullptr && (sequence.compare(checkStr2) != 0 && rev_seq.compare(checkStr2) != 0)))
        {
            std::cout << "anticipation of checkStr failing, sequence: " << sequence << " checkStr: " << checkStr << std::endl;
		    graph.debugPrintAllUnitigs();
        }

        CPPUNIT_ASSERT ((sequence.compare(checkStr)==0 || rev_seq.compare(checkStr) == 0) \
           || (checkStr2 != nullptr && (sequence.compare(checkStr2) == 0 || rev_seq.compare(checkStr2) == 0)) );
    }

    /********************************************************************************/

    /********************************************************************************/
    void debruijn_simplunitigs_X ()
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
        GraphUnitigs graph;
        debruijn_build_entry r;
        debruijn_build(sequences, ARRAY_SIZE(sequences), kmerSize, graph, r);

        DEBUGprint(std::cout << "nb nodes:" << r.nbNodes << "nb nodes non deleted " << r.nbNonDeletedNodes << std::endl;)
    
        CPPUNIT_ASSERT (r.nbNodes == 10);
        if (r.nbNonDeletedNodes != 10)
            std::cout << "ancipating failed assert: nbNonDeletedNodes = " << r.nbNonDeletedNodes << " != 14" << std::endl;
        CPPUNIT_ASSERT (r.nbNonDeletedNodes == 10);

        // simplify it
        graph.simplify(1, false); // one core, no verbose

        r = debruijn_stats (graph, true,  true);
        DEBUGprint(std::cout << "post simp:" << std::endl;)
        DEBUGprint(std::cout << "nb nodes:" << r.nbNodes << "nb nodes non deleted " << r.nbNonDeletedNodes << std::endl;)
        CPPUNIT_ASSERT (r.nbNodes == 2);
        if (r.nbNonDeletedNodes != 2)
            std::cout << "ancipating failed assert: nbNonDeletedNodes = " << r.nbNonDeletedNodes << " != 14" << std::endl;

        CPPUNIT_ASSERT (r.nbNonDeletedNodes == 2); // only two nodes remaining: AAAAA and CAAGA
        
   }

    void debruijn_simplunitigs_tip ()
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
        GraphUnitigs graph;
        debruijn_build_entry r;
        debruijn_build(sequences, ARRAY_SIZE(sequences), kmerSize, graph, r);

        DEBUGprint(std::cout << "nb nodes:" << r.nbNodes << "nb nodes non deleted " << r.nbNonDeletedNodes << std::endl;)
        CPPUNIT_ASSERT (r.nbNodes == 6); // should now be identical because GraphUnitigs doesn't iterate on deleted nodes anymore
        CPPUNIT_ASSERT (r.nbNonDeletedNodes == 6);

        // simplify it
        graph.simplify(1, false); // one core, no verbose


        // how many nodes left? should be as many as initially. it's a negative test: graph shouldn't be simplified
        r = debruijn_stats (graph, true,  true);
        DEBUGprint(std::cout << "post simp:" << std::endl;)
        DEBUGprint(std::cout << "nb nodes:" << r.nbNodes << "nb nodes non deleted " << r.nbNonDeletedNodes << std::endl;)
        CPPUNIT_ASSERT (r.nbNodes == 4);
        CPPUNIT_ASSERT (r.nbNonDeletedNodes == 4);

        debruijn_traversal (graph, sequences[0], "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACCGTTATGCCGTCCGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACATG");
    }


    /********************************************************************************/


    void debruijn_simplunitigs_bubble_aux (const char* sequences[], int nb_seqs, const char* sol1, const char* sol2=nullptr)
    {
        size_t kmerSize = 21;

        // We create the graph.
        GraphUnitigs graph;
        debruijn_build_entry r;
        debruijn_build(sequences, nb_seqs, kmerSize, graph, r);

        DEBUGprint(std::cout << "nb nodes:" << r.nbNodes << "nb nodes non deleted " << r.nbNonDeletedNodes << std::endl;)
        CPPUNIT_ASSERT (r.nbNodes == 8);
        CPPUNIT_ASSERT (r.nbNonDeletedNodes == 8);

        // simplify it
        graph.simplify(1, false); // one core, no verbose


        r = debruijn_stats (graph, true,  true);
        DEBUGprint(std::cout << "post simp:" << std::endl;)
        DEBUGprint(std::cout << "nb nodes:" << r.nbNodes << "nb nodes non deleted " << r.nbNonDeletedNodes << std::endl;)
        CPPUNIT_ASSERT (r.nbNodes == 6);
        CPPUNIT_ASSERT (r.nbNonDeletedNodes == 6);

        if (0) // do a useless second pass
        {
            // simplify it
            graph.simplify(1, false); // one core, no verbose


            // how many nodes left? should be as many as initially. it's a negative test: graph shouldn't be simplified
            r = debruijn_stats (graph, true,  true);
            //
            DEBUGprint(std::cout << "post simp-again:" << std::endl;)
                DEBUGprint(std::cout << "nb nodes:" << r.nbNodes << "nb nodes non deleted " << r.nbNonDeletedNodes << std::endl;)
                CPPUNIT_ASSERT (r.nbNodes == 8);
            CPPUNIT_ASSERT (r.nbNonDeletedNodes == 6);
        }

        debruijn_traversal (graph, sequences[0], sol1, sol2);
    }

    void debruijn_simplunitigs_bubble ()
    {
 
        const char* sequences[] =
        {
            "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACC", //>works well for k=21; part of genome10K.fasta
            "TGTCATCTAGTTCAACAACCAAAATAACGACTCTTGCGCTCGGATGT", //>that's the bubble  (highly covered)
            "TGTCATCTAGTTCAACAACCAAAATAACGACTCTTGCGCTCGGATGT", //>that's the bubble
            "TGTCATCTAGTTCAACAACCAAAATAACGACTCTTGCGCTCGGATGT", //>that's the bubble
            "TGTCATCTAGTTCAACAACCAAAAAAACGACTCTTGCGCTCGGATGT", //>that's the bubble path 2, low covered
            "CGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACAAA" // >remaining part
        };
        const char* sol = "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACCAAAATAACGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACAAA";
    
        debruijn_simplunitigs_bubble_aux(sequences, ARRAY_SIZE(sequences),sol);
    }

    void debruijn_simplunitigs_bubble_snp()
    {
        const char* sequences[] =
        {
            "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACC", //>works well for k=21; part of genome10K.fasta
            "TGTCATCTAGTTCAACAACCAAAATAACGACTCTTGCGCTCGGATGT", //>that's the bubble
            "TGTCATCTAGTTCAACAACCAAAATAACGACTCTTGCGCTCGGATGT", //>that's the bubble
            "TGTCATCTAGTTCAACAACCAAAATAACGACTCTTGCGCTCGGATGT", //>that's the bubble
            "TGTCATCTAGTTCAACAACCAAAAAAACGACTCTTGCGCTCGGATGT", //>that's the bubble path 2, same coverage
            "TGTCATCTAGTTCAACAACCAAAAAAACGACTCTTGCGCTCGGATGT", //>that's the bubble path 2,
            "TGTCATCTAGTTCAACAACCAAAAAAACGACTCTTGCGCTCGGATGT", //>that's the bubble path 2, 
            "TGTCATCTAGTTCAACAACCAAAAAAACGACTCTTGCGCTCGGATGT", //>that's the bubble path 2, 
            "CGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACAAA" // >remaining part
        };
        const char* sol1 = "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACCAAAATAACGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACAAA";
        const char* sol2 = "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACCAAAAAAACGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACAAA";

        debruijn_simplunitigs_bubble_aux(sequences, ARRAY_SIZE(sequences), sol1, sol2);
 
    }
 

    void debruijn_simplunitigs_ec ()
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
        GraphUnitigs graph;
        debruijn_build_entry r;
        debruijn_build(sequences, ARRAY_SIZE(sequences), kmerSize, graph, r);

        DEBUGprint(std::cout << "nb nodes:" << r.nbNodes << "nb nodes non deleted " << r.nbNonDeletedNodes << std::endl;)
        CPPUNIT_ASSERT (r.nbNodes == 10);
        CPPUNIT_ASSERT (r.nbNonDeletedNodes == 10);

        // simplify it
        graph.simplify(1, false); // one core, no verbose


        // how many nodes left? should be as many as initially. it's a negative test: graph shouldn't be simplified
        r = debruijn_stats (graph, true,  true);
        DEBUGprint(std::cout << "post simp:" << std::endl;)
        DEBUGprint(std::cout << "nb nodes:" << r.nbNodes << "nb nodes non deleted " << r.nbNonDeletedNodes << std::endl;)
        CPPUNIT_ASSERT (r.nbNodes == 8);
        if (r.nbNonDeletedNodes != 8)
            std::cout << " in anticipation of assert fail, r.nbNonDeletedNodes = " << r.nbNonDeletedNodes << std::endl;
        CPPUNIT_ASSERT (r.nbNonDeletedNodes == 8);

        debruijn_traversal (graph, sequences[0], "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACCGTTATGCCGTCCGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACATG");
        debruijn_traversal (graph, sequences[3], "GGTGAACAGCACATCTTTTCGTCCTGAGGCCATATTAATTCTACTCAGATTGTCTGTAACCGGAGCTTCGGGCGTATTTTTGCGTAAGACACTGCCTAAAGGGAACATATGTGTCCAGAATAGGGTTCAACGGTGTATGAGCAAACTAGTTCAACAACCAAAAAAATTGTGTGCAAGCTACTTCTAGACCTTATTAAGTGCCCAGGAATTCCTAGGAAGGCGCGCAGCTCAAGCAATCATACATGGCGGAATGCCTGTCCACCGGGGGTTCTACTGTACCACAGTGGCCTGGATAGCTAAGCAGGTCCTGGATTGGCATGTCATCCGGAGTGATAGGCACTGCTCACGACCAGCTTGCGGACAAACGGGGTGCCCGCGCCTGCGTCCGGTAGACGAGCGATGGATTTAGACCGTTCACTGAACCCTCTAATAGGACCTCTTGCCCATCCGAGGCTTAAGC");
    }

};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestSimplificationsUnitigs);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestSimplificationsUnitigs);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

