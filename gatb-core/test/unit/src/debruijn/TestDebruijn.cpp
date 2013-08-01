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

#include <gatb/tools/math/NativeInt64.hpp>
#include <gatb/tools/math/NativeInt128.hpp>
#include <gatb/tools/math/LargeInt.hpp>

#include <gatb/tools/collections/impl/IteratorFile.hpp>

#include <gatb/tools/misc/api/StringsRepository.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <gatb/debruijn/api/IGraph.hpp>
#include <gatb/debruijn/impl/GraphBasic.hpp>
#include <gatb/debruijn/impl/GraphFactory.hpp>

#include <iostream>
#include <memory>

using namespace std;

using namespace gatb::core::debruijn;
using namespace gatb::core::debruijn::impl;

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

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for genomic databases management
 */
class TestDebruijn : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestDebruijn);

        CPPUNIT_TEST_GATB (debruijn_test1);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {}
    void tearDown ()  {}

    /********************************************************************************/
    template<typename T>
    struct Info
    {
        Info() : checksum(0), nbNodes(0), nbNeighbors(0) {}
        T checksum;
        size_t nbNodes;
        size_t nbNeighbors;

        /** */
        void incNodes ()      { nbNodes++;  }
        void inc      (T& t)  { checksum = checksum + t; nbNeighbors++;  }

        /** */
        string toString ()
        {
            stringstream ss;
            ss << "[INFO  nbNodes=" << nbNodes << "  nbNeighbors=" << nbNeighbors << "  checksum=" << checksum << "]";
            return ss.str();
        }
    };

    /********************************************************************************/
    template<typename T>
    void debruijn_test1_aux (Graph<T>& graph)
    {
        /** We get an iterator over all the nodes of the graph. */
        NodeIterator<T> itNodes = graph.nodes();

        NodeSet<T> nodes (graph);
        Info<T>    info;

        TimeInfo ti;
        {
            TIME_INFO (ti, "loop");

            /** We iterate all the nodes of the graph. */
            for (itNodes.first(); !itNodes.isDone(); itNodes.next())
            {
                info.incNodes();

                /** We retrieve the successors. */
                size_t nbNodes = graph.getSuccessors (*itNodes, nodes);

                /** We iterate all the successors. */
                for (size_t i=0; i<nbNodes; i++)   {  info.inc (nodes[i].kmer);  }
            }
        }

        cout << info.toString() <<  "  time=" << ti.getEntryByKey("loop") << endl;
    }


    /********************************************************************************/
    template<typename T>
    void debruijn_test2_aux (Graph<T>& graph)
    {
        /** We get an iterator over all the nodes of the graph. */
        NodeIterator<T> itNodes = graph.nodes();

        EdgeSet<T> edges (graph);
        Info<T>    info;

        TimeInfo ti;
        {
            TIME_INFO (ti, "loop");

            /** We iterate all the nodes of the graph. */
            for (itNodes.first(); !itNodes.isDone(); itNodes.next())
            {
                info.incNodes();

                /** We retrieve the outcoming edges. */
                size_t nbEdges = graph.getOutEdges (*itNodes, edges);

                /** We iterate all outcoming edges. */
                for (size_t i=0; i<nbEdges; i++)  {  info.inc (edges[i].to.kmer);  }
            }
        }

        cout << info.toString() <<  "  time=" << ti.getEntryByKey("loop") << endl;
    }


    /********************************************************************************/
    template<typename T>
    void debruijn_test3_aux (Graph<T>& graph)
    {
        EdgeSet<T> edges (graph);

        /** We get an iterator over all the nodes of the graph. */
        NodeIterator<T> itNodes = graph.nodes();

        /** We retrieve the first node. */
        itNodes.first();  //for (size_t i=1; i<=5; i++)  { itNodes.next(); }
        Node<T> node =  itNodes.item();

        cout << "------- NODE " << node.toString() << endl;

        Strand strandInit = node.strand;

        size_t i=0;
        for (size_t nbEdges = ~0 ; nbEdges > 0; i++, node = edges[0].to)
        {
            nbEdges = graph.getOutEdges (node, edges);

            if (nbEdges > 0)  {  cout << "-------------- " << edges[0].toString (strandInit) << " --------------" << endl; }
        }
        cout << "nb found " << i << endl;
    }

    /********************************************************************************/
    void debruijn_test1 ()
    {
#if 0
        Graph<NativeInt64> graph = GraphFactory::createGraph <NativeInt64> (
            new Property (STR_KMER_SOLID,  "/local/users/edrezen/projects/GATB/gforge/gatb-tools/gatb-tools/tools/debloom/build/tmp.solid"),
            new Property (STR_KMER_CFP,    "/local/users/edrezen/projects/GATB/gforge/gatb-tools/gatb-tools/tools/debloom/build/tmp.debloom"),
            new Property (STR_KMER_SIZE,   "27"),
            PROP_END
        );
        //debruijn_test1_aux<NativeInt64> (graph);
        //debruijn_test2_aux<NativeInt64> (graph);
        debruijn_test3_aux<NativeInt64> (graph);
#endif
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION (TestDebruijn);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

