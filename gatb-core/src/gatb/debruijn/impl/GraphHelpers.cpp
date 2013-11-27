/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/debruijn/impl/GraphHelpers.hpp>

using namespace std;


#define DEBUG(a)  //printf a

extern u_int64_t incomingTable[];
extern bool hack;


/********************************************************************************/
namespace gatb {  namespace core {  namespace debruijn {  namespace impl {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
int GraphHelper::simplePathAvance (const Node& node, Direction dir) const
{
    Edge output;  return simplePathAvance (node, dir, output);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
int GraphHelper::simplePathAvance (const Node& node, Direction dir, Edge& output) const
{
    std::vector<Edge> neighbors = _graph.neighbors<Edge> (node, dir);

    /** We check we have no outbranching. */
    if (neighbors.size() == 1)
    {
        /** We check whether the neighbor has an inbranching or not. */
        if (_graph.degree (neighbors[0].to, reverse(dir)) > 1)  { return -2; }

        output = neighbors[0];

		if (hack)
		{
		    if (output.direction == DIR_INCOMING)
		    {
		        size_t span = _graph.getKmerSize();

		        if (output.to.strand == kmer::STRAND_FORWARD)
		        {
		            output.nt = ((kmer::Nucleotide) (output.to.kmer[span-1]));
		            output.nt = (kmer::Nucleotide) incomingTable[output.nt];
		        }
		        else
		        {
		            output.nt = kmer::reverse ((kmer::Nucleotide) (output.to.kmer[0]));
		            output.nt = (kmer::Nucleotide) incomingTable[output.nt];
		        }
		    }
		}

        return 1;
    }

    /** if this kmer has out-branching, don't extend it. */
    if (neighbors.size() > 1)  {  return -1;  }

    return 0;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Graph::Iterator<Node> GraphHelper::getSimplePathNodeIterator (const Node& node, Direction dir) const
{
    class LocalIterator : public tools::dp::ISmartIterator<Node>
    {
    public:

        LocalIterator (const GraphHelper& helper, const Node& node, Direction dir)
            : _helper(helper), _graph(helper.getGraph()), _dir(dir), _isDone(true)
        {
            *(this->_item) = node;
        }

        /** \copydoc  Iterator::first */
        void first()
        {
            _isDone = false;
            next ();
        }

        /** \copydoc  Iterator::next */
        void next()
        {
            Edge output;

            if (_helper.simplePathAvance (*(this->_item), _dir, output) > 0)
            {
                /** We update the currently iterated node. */
                *(this->_item) = output.to;
            }
            else
            {
                /** We can't have a non branching node in the wanted direction => iteration is finished. */
                _isDone = true;
            }
        }

        /** \copydoc  Iterator::isDone */
        bool isDone() {  return _isDone;  }

        /** \copydoc  Iterator::item */
        Node& item ()  {  return *(this->_item);  }

        /** */
        u_int64_t getNbItems () const { return 0; }

    private:
        const GraphHelper& _helper;
        const Graph&       _graph;
        Direction          _dir;
        bool               _isDone;
    };

    return Graph::Iterator<Node> (new LocalIterator(*this, node, dir));
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
