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

#include <gatb/debruijn/impl/GraphHelpers.hpp>

using namespace std;


#define DEBUG(a)  //printf a

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
int GraphHelper::simplePathAvance (const Node& node, Direction dir, kmer::Nucleotide& nt) const
{
    Edge edge;
    int res = simplePathAvance (node, dir, edge);
    nt = edge.nt;
    return res;
}

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
    Graph::Vector<Edge> neighbors = _graph.neighbors<Edge> (node, dir);

    /** We check we have no outbranching. */
    if (neighbors.size() == 1)
    {
        /** We check whether the neighbor has an inbranching or not. */
        if (_graph.degree (neighbors[0].to, reverse(dir)) > 1)  { return -2; }

        output = neighbors[0];

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
            : _helper(helper), _graph(helper.getGraph()), _dir(dir), _rank(0), _isDone(true)
        {
            *(this->_item) = node;
        }

        u_int64_t rank () const { return _rank; }

        /** \copydoc  Iterator::first */
        void first()
        {
            _rank   = 0;
            _isDone = false;
            next ();
        }

        /** \copydoc  Iterator::next */
        void next()
        {
            Edge output;

            _rank ++;

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
        u_int64_t size () const { return 0; }

    private:
        const GraphHelper& _helper;
        const Graph&       _graph;
        Direction          _dir;
        u_int64_t          _rank;
        bool               _isDone;
    };

    return Graph::Iterator<Node> (new LocalIterator(*this, node, dir));
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
