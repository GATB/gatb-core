/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
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

#ifndef _GATB_TOOLS_TERMINATOR_HPP_
#define _GATB_TOOLS_TERMINATOR_HPP_

/********************************************************************************/

#include <gatb/debruijn/impl/Graph.hpp>
#include <vector>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
namespace impl      {
/********************************************************************************/

/** \brief Interface that allows to mark certain nodes in a graph. Note: It is on the verge of becoming deprecated.
 *
 * In short: don't use it. Use the MPHF.
 *
 * This is a legacy class that was used when the graph did not have a MPHF.
 * I.e. in Minia version 1.
 * Back then, it was not possible to directly tag nodes inside the graph structure.
 * So this Terminator class uses an external data structure to mark some of 
 * the nodes (specifically, the branching nodes). It was used during graph traversal 
 * for instance.
 *
 * The Terminator interface provides such a service. It is often used
 * by the \ref Traversal class during its traversing process.
 *
 * The name Terminator is just related to termination, we needed a way to ensure that we don't traverse
 * the same contig twice, so that the algorithm terminates eventually.
 */

template <typename Node, typename Edge, typename Graph>
class TerminatorTemplate : public system::SmartPointer
{
public:

    /** Constructor
     * \param[in] graph : the graph */
    TerminatorTemplate (const Graph& graph) : _graph(graph)  {}

    /** Destructor. */
    virtual ~TerminatorTemplate ()  {}

    /** Get the graph.
     * \return the graph. */
    const Graph& getGraph() const { return _graph; }

    /** Tells whether marking nodes is allowed or not
     * \return true if marking nodes is allowed, false otherwise. */
    virtual bool isEnabled () const { return true; }

    /** Mark the provided edge
     * \param[in] edge : edge to be marked. */
    virtual void mark      (Edge& edge) = 0;

    /** Tells whether an edge is marked
     * \param[in] edge : edge to be checked
     * \return true if the edge is marked, false otherwise
     */
    virtual bool is_marked (Edge& edge)  const = 0;

    /** Mark the provided node.
     * \param[in] node : node to be marked. */
    virtual void mark      (Node& node) = 0;

    /** Tells whether a node is marked
     * \param[in] node : node to be checked
     * \return true if the node is marked, false otherwise
     */
    virtual bool is_marked (Node& node)  const = 0;

    /** Tells whether a branching node is marked
     * \param[in] node : node to be checked
     * \return true if the node is marked, false otherwise
     */
    virtual bool is_marked_branching (Node& node) const = 0;

    /** Tells whether a node is branching
     * \param[in] node : node to be checked
     * \return true if the node is branching, false otherwise.
     */
    virtual bool is_branching (Node& node) const = 0;

    /** Reset the current marked nodes/edges. */
    virtual void reset () = 0;

    /** Dump (for debug purpose). */
    virtual void dump () = 0;

protected:

    const Graph& _graph;
};

/********************************************************************************/

/** \brief Null implementation of Terminator.
 *
 * This class provides a null implementation, which means that no marks are done.
 */
template <typename Node, typename Edge, typename Graph>
class NullTerminatorTemplate :  public TerminatorTemplate<Node,Edge,Graph>
{
public:

    /** Singleton method.
     * \return a singleton instance.
     */
    static TerminatorTemplate<Node,Edge,Graph>& singleton()
    {
        static Graph dummy;
        static NullTerminatorTemplate instance(dummy); return instance;
    }

    /** \copydoc Terminator::isEnabled */
    virtual bool isEnabled () const { return false; }

    /** \copydoc Terminator::mark */
    virtual void mark      (Edge& edge) {}

    /** \copydoc Terminator::is_marked */
    virtual bool is_marked (Edge& edge)  const  { return false; };

    /** \copydoc Terminator::mark(const gatb::core::debruijn::impl::Node&)  */
    virtual void mark      (Node& node) {}

    /** \copydoc Terminator::is_marked(const gatb::core::debruijn::impl::Node&) const */
    virtual bool is_marked (Node& node)  const  { return false; }

    /** \copydoc Terminator::is_marked_branching */
    virtual bool is_marked_branching (Node& node) const { return false; }

    /** \copydoc Terminator::is_branching */
    virtual bool is_branching (Node& node) const { return false; }

    /** \copydoc Terminator::reset */
    virtual void reset () {}

    /** \copydoc Terminator::dump */
    virtual void dump () {}

//private:

    NullTerminatorTemplate (const Graph& graph) : TerminatorTemplate<Node,Edge,Graph>(graph)  {}
};

/********************************************************************************/

/** \brief Implementation of Terminator that marks branching nodes.
 */
template <typename Node, typename Edge, typename Graph>
class BranchingTerminatorTemplate :  public TerminatorTemplate<Node,Edge,Graph>
{
public:

    typedef unsigned short int  Value;

    /** Constructor
     * \param[in] graph : the graph */
    BranchingTerminatorTemplate (const Graph& graph);

    /** Copy constructor
     * \param[in] terminator: the graph */
    BranchingTerminatorTemplate (const BranchingTerminatorTemplate& terminator);

    /** Destructor. */
    ~BranchingTerminatorTemplate();

    /** \copydoc Terminator::mark */
    void mark      (Edge& edge);

    /** \copydoc Terminator::is_marked */
    bool is_marked (Edge& edge)  const;

    /** \copydoc Terminator::mark(const gatb::core::debruijn::impl::Node&)  */
    void mark      (Node& node);

    /** \copydoc Terminator::is_marked(const gatb::core::debruijn::impl::Node&) const */
    bool is_marked (Node& node) const ;

    /** \copydoc Terminator::is_marked_branching */
    bool is_marked_branching (Node& node) const ;

    /** \copydoc Terminator::is_branching */
    bool is_branching (Node& node) const ;

    /** \copydoc Terminator::reset */
    void reset();

    /** \copydoc Terminator::dump */
    void dump ();

private:

    /* Custom implementation of a map.
     * IMPORTANT : The keys set is supposed to be built only by one instance and is shared
     * with other instances through the copy constructor.
     */
    template <typename Key, typename Value>  class AssocSet
    {
    public:

        AssocSet () : isRef(false)  {  keys = new std::vector<Key>();  }

        AssocSet (const AssocSet& other) : isRef(true)
        {
        	keys = other.keys;
            liste_value.assign(this->keys->size(),0);
        }

        ~AssocSet ()  {  if (isRef==false)  { delete keys; }  }

    	void insert (const Key& elem) { keys->push_back(elem); }

        bool contains(const Key& elem)  const  {  return binary_search(keys->begin(), keys->end(), elem);  }

        int get (const Key& elem, Value& val) const
        {
            typename std::vector<Key>::const_iterator it;
            it = lower_bound(this->keys->begin(), this->keys->end(),elem);
            if (it == this->keys->end() || elem != *it) return 0;
            size_t rank = it - this->keys->begin();
            val = liste_value[rank];
            return 1;
        }

        int set (const Key& elem, const Value& val)
        {
            typename  std::vector<Key>::iterator it;
            it = lower_bound(this->keys->begin(), this->keys->end(),elem);
            if (it == this->keys->end() ||elem != *it) return 0;
            size_t rank = it - this->keys->begin();
            liste_value[rank]=val;
            return 1;
        }

        void finalize (bool doSort=true)
        {
            if (doSort) {  sort(this->keys->begin(), this->keys->end());  }
            liste_value.assign(this->keys->size(),0);
        }

        void clear()  { liste_value.assign(liste_value.size(),0); }

    private:
        std::vector<Key>*  keys;
        std::vector<Value> liste_value;
        bool isRef;
    };


    bool is_indexed (Node& node) const ;

    AssocSet<typename Node::Value, Value> branching_kmers;

    int getDelta (Edge& edge) const;
};


/** \brief MPHF implementation of Terminator.
 *
 */
template <typename Node, typename Edge, typename Graph>
class MPHFTerminatorTemplate :  public TerminatorTemplate<Node,Edge,Graph>
{
public:

    /** \copydoc Terminator::mark */
    virtual void mark      (Edge& edge)  { printf("not expecting a call to MPHFTermiantor.mark(edge)\n"); exit(1); }

    /** \copydoc Terminator::is_marked */
    virtual bool is_marked (Edge& edge)  const  { printf("not expecting a call to MPHFTermiantor.is_marked(edge)\n"); exit(1); return false; }

    /** \copydoc Terminator::mark(const gatb::core::debruijn::impl::Node&)  */
    virtual void mark      (Node& node);

    /** \copydoc Terminator::is_marked(const gatb::core::debruijn::impl::Node&) const */
    virtual bool is_marked (Node& node)  const  ;

    /** \copydoc Terminator::is_marked_branching */
    virtual bool is_marked_branching (Node& node) const { printf("not expecting a call to MPHFTermiantor.is_marked_branching\n"); exit(1); return false; }

    /** \copydoc Terminator::is_branching */
    virtual bool is_branching (Node& node) const { printf("not expecting a call to MPHFTermiantor.is_branching\n"); exit(1);return false; }

    /** \copydoc Terminator::reset */
    virtual void reset ();

    /** \copydoc Terminator::dump */
    virtual void dump () { printf("not expecting a call to MPHFTermiantor.dump\n"); exit(1); };

    MPHFTerminatorTemplate (const Graph& graph) : TerminatorTemplate<Node,Edge,Graph>(graph)  {}
};

typedef TerminatorTemplate<Node, Edge, Graph> Terminator; 
typedef MPHFTerminatorTemplate<Node, Edge, Graph> MPHFTerminator; 
typedef BranchingTerminatorTemplate<Node, Edge, Graph> BranchingTerminator; 


/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_TOOLS_TERMINATOR_HPP_ */

