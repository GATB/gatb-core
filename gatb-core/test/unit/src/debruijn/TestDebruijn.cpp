/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Copyright (c) 2013                                                      *
 *                                                                           *
 *   GATB is free software; you can redistribute it and/or modify it under   *
 *   the CECILL version 2 License, that is compatible with the GNU General   *
 *   Public License                                                          *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   CECILL version 2 License for more details.                              *
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
#include <gatb/kmer/impl/DebloomAlgorithm.hpp>

#include <gatb/bank/impl/BankStrings.hpp>

#include <gatb/tools/collections/impl/Product.hpp>
#include <gatb/tools/collections/impl/ProductFile.hpp>
#include <gatb/tools/collections/impl/ProductHDF5.hpp>

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

typedef LargeInt<1>  LocalInteger;

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for genomic databases management
 */
class TestDebruijn : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestDebruijn);

//        CPPUNIT_TEST_GATB (debruijn_test1);
//        CPPUNIT_TEST_GATB (debruijn_test2);
        CPPUNIT_TEST_GATB (debruijn_test3);
        CPPUNIT_TEST_GATB (debruijn_test4);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {}
    void tearDown ()  {}

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
    void debruijn_test1_aux (Graph& graph)
    {
        /** We get an iterator over all the nodes of the graph. */
        Graph::Iterator<Node> itNodes = graph.iterator<Node>();

        std::vector<Node> successors;
        Info    info;

        TimeInfo ti;
        {
            TIME_INFO (ti, "loop");

            /** We iterate all the nodes of the graph. */
            for (itNodes.first(); !itNodes.isDone(); itNodes.next())
            {
                info.incNodes();

                /** We retrieve the successors. */
                successors = graph.successors<Node> (itNodes.item());

                /** We iterate all the successors. */
                for (size_t i=0; i<successors.size(); i++)   {  info.inc (successors[i].kmer);  }
            }
        }

        cout << info.toString() <<  "  time=" << ti.getEntryByKey("loop") << endl;
    }


    /********************************************************************************/
    void debruijn_test2_aux (Graph& graph)
    {
        /** We get an iterator over all the nodes of the graph. */
        Graph::Iterator<Node> itNodes = graph.iterator<Node>();

        std::vector<Edge> successors;
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

        cout << info.toString() <<  "  time=" << ti.getEntryByKey("loop") << endl;
    }


    /********************************************************************************/
    void debruijn_test3_aux (const Graph& graph)
    {
        /** We get an iterator over all the nodes of the graph. */
        Graph::Iterator<Node> itNodes = graph.iterator<Node> ();

        /** We retrieve the first node. */
        itNodes.first();  //for (size_t i=1; i<=5; i++)  { itNodes.next(); }
        Node node =  itNodes.item();

        cout << "------- NODE " << graph.toString(node) << endl;

        Strand strandInit = node.strand;

        size_t i=0;
        for (std::vector<Edge> successors; (successors = graph.successors<Edge>(node)).size() > 0; i++, node = successors[0].to)
        {
        }
        cout << "nb found " << i << endl;
    }

    /********************************************************************************/
    void debruijn_test1 ()
    {
#if 0
        Graph<LocalInteger> graph = GraphFactory::createGraph <LocalInteger> (
            new Property (STR_KMER_SOLID,  "/local/users/edrezen/projects/GATB/gforge/gatb-tools/gatb-tools/tools/debloom/build/tmp.solid"),
            new Property (STR_KMER_CFP,    "/local/users/edrezen/projects/GATB/gforge/gatb-tools/gatb-tools/tools/debloom/build/tmp.debloom"),
            new Property (STR_KMER_SIZE,   "27"),
            PROP_END
        );
        //debruijn_test1_aux<LocalInteger> (graph);
        //debruijn_test2_aux<LocalInteger> (graph);
        debruijn_test3_aux<LocalInteger> (graph);
#endif
    }

    /********************************************************************************/
    void debruijn_check_sequence (const Graph& graph, size_t kmerSize, const char* seq)
    {
        size_t seqLen = strlen (seq);

        Graph::Iterator<Node> nodes = graph.iterator<Node> ();
        nodes.first ();

        /** We get the first node. */
        Node node = nodes.item();

        /** We compute the branching range for the node. */
        Node begin, end;     graph.getNearestBranchingRange (node, begin, end);

        /** We check that the begin kmer matches the beginning of the sequence. */
        bool check1 =
            graph.toString (begin)        == string (seq, kmerSize)  ||
            graph.toString (reverse(end)) == string (seq, kmerSize);

        /** We check that the end kmer matches the end of the sequence. */
        bool check2 =
            graph.toString (end)            == string (seq + seqLen - kmerSize, kmerSize)  ||
            graph.toString (reverse(begin)) == string (seq + seqLen - kmerSize, kmerSize);

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
    template<typename ProductFactory>
    void debruijn_test2_aux (size_t kmerSize, size_t nks, const char* seq)
    {
        size_t seqLen   = strlen (seq);

        CPPUNIT_ASSERT (seqLen >= kmerSize);

        /** We create a bank with one sequence. */
        IBank* bank = new BankStrings (seq, 0);

        /** We create a product instance. */
        Product<ProductFactory>* product  = ProductFactory::createProduct ("test", true, true);
        LOCAL (product);

        /** We create a DSK instance. */
        SortingCountAlgorithm<ProductFactory, LocalInteger> sortingCount (product, bank, kmerSize, nks);

        /** We launch DSK. */
        sortingCount.execute();

        /** We check that the sequence has no duplicate kmers. */
        CPPUNIT_ASSERT ( (seqLen - kmerSize + 1) == sortingCount.getSolidKmers()->getNbItems());

        /** We create a debloom instance. */
        DebloomAlgorithm<ProductFactory, LocalInteger> debloom (*product, sortingCount.getSolidKmers(), kmerSize);

        /** We launch the debloom. */
        debloom.execute();

        /** We create the graph. */
//        Graph graph = GraphFactory::createGraph (sortingCount.getSolidKmers(), debloom.getCriticalKmers(), kmerSize);
//
//        debruijn_check_sequence (graph, kmerSize, seq);
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

        size_t kmerSizes[] = { 13, 15, 17, 19, 21, 23, 25, 27, 29, 31};
        size_t nks=1;

        for (size_t i=0; i<ARRAY_SIZE(sequences); i++)
        {
            for (size_t j=0; j<ARRAY_SIZE(kmerSizes); j++)
            {
                debruijn_test2_aux <ProductHDF5Factory> (kmerSizes[j], nks, sequences[i]);
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
                /** We define options for the graph generation. */
                IProperties* props = new Properties();  LOCAL (props);
                props->add (0, STR_KMER_SIZE, "%d", kmerSizes[j]);
                props->add (0, STR_NKS,       "%d", 1);

                /** We create the graph. */
                Graph graph = Graph::create (new BankStrings (sequences[i], 0), props);

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
        Graph graph = Graph::create (new BankStrings (seq, 0));  // kmerSize=27 by default

        Graph::Iterator<Node> it = graph.iterator<Node>();  it.first();

        Node n1 = it.item();
        CPPUNIT_ASSERT (n1.strand == STRAND_FORWARD);
        CPPUNIT_ASSERT (graph.toString(n1).compare (seq) == 0);

        Node n2 = reverse (n1);
        CPPUNIT_ASSERT (n2.strand == STRAND_REVCOMP);
        CPPUNIT_ASSERT (graph.toString(n2).compare (rev) == 0);
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestDebruijn);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestDebruijn);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

