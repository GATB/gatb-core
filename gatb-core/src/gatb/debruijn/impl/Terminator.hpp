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

template <typename Item> class ListSet
{
public:
    std::vector<Item> liste;
public:
    void insert (const Item& elem)  {     liste.push_back(elem); }
    void finalize() {    sort(liste.begin(), liste.end()); }
    bool contains(const Item& elem)  const {    return binary_search(liste.begin(), liste.end(), elem);  }
    uint64_t capacity() {return (u_int64_t)liste.capacity();}
    static const int bits_per_element = sizeof(Item)*8;
    ListSet(u_int64_t taille_approx)  {   liste.reserve(taille_approx);  }
    ListSet()  {}
};

/********************************************************************************/

template <typename Key, typename Value>  class AssocSet : public ListSet<Key>
{
public:
    std::vector<Value> liste_value;

public:
    int get (const Key& elem, Value& val) const
    {
        typename std::vector<Key>::const_iterator it;
        it = lower_bound(this->liste.begin(), this->liste.end(),elem);
        if (it == this->liste.end() || elem != *it) return 0;
        size_t rank = it - this->liste.begin();
        val = liste_value[rank];
        return 1;
    }

    int set (const Key& elem, const Value& val)
    {
        typename  std::vector<Key>::iterator it;
        it = lower_bound(this->liste.begin(), this->liste.end(),elem);
        if (it == this->liste.end() ||elem != *it) return 0;
        size_t rank = it - this->liste.begin();
        liste_value[rank]=val;
        return 1;
    }

    void finalize (bool doSort=true)
    {
        if (doSort) {  sort(this->liste.begin(), this->liste.end());  }
        liste_value.assign(this->liste.size(),0);
    }

    AssocSet() {}

    void clear()  { liste_value.assign(liste_value.size(),0); }


    void start_iterator()  {   iterator = this->liste.begin()-1;  }
    bool next_iterator()
    {
        iterator++;
        if (iterator==this->liste.end()) return false;
        return true;
    }

    typename std::vector<Key>::iterator iterator;
};

/********************************************************************************/

/** \brief Interface that allows to mark nodes in a graph.
 *
 * De Bruijn graphs in GATB are immutable and therefore it is not
 * possible to directly tag nodes inside the graph structure.
 *
 * Nevertheless, we can use external data structure to mark nodes,
 * during graph traversal for instance.
 *
 * The Terminator interface provides such a service. It is often used
 * by the \ref Traversal class during its traversing process.
 */
class Terminator : public system::SmartPointer
{
public:

    /** Constructor
     * \param[in] graph : the graph */
    Terminator (const Graph& graph) : _graph(graph)  {}

    /** Destructor. */
    virtual ~Terminator ()  {}

    /** Get the graph.
     * \return the graph. */
    const Graph& getGraph() const { return _graph; }

    /** Tells whether marking nodes is allowed or not
     * \return true if marking nodes is allowed, false otherwise. */
    virtual bool isEnabled () const { return true; }

    /** Mark the provided edge
     * \param[in] edge : edge to be marked. */
    virtual void mark      (const Edge& edge) = 0;

    /** Tells whether an edge is marked
     * \param[in] edge : edge to be checked
     * \return true if the edge is marked, false otherwise
     */
    virtual bool is_marked (const Edge& edge)  const = 0;

    /** Mark the provided node.
     * \param[in] node : node to be marked. */
    virtual void mark      (const Node& node) = 0;

    /** Tells whether a node is marked
     * \param[in] node : node to be checked
     * \return true if the node is marked, false otherwise
     */
    virtual bool is_marked (const Node& node)  const = 0;

    /** Tells whether a branching node is marked
     * \param[in] node : node to be checked
     * \return true if the node is marked, false otherwise
     */
    virtual bool is_marked_branching (const Node& node) const = 0;

    /** Tells whether a node is branching
     * \param[in] node : node to be checked
     * \return true if the node is branching, false otherwise.
     */
    virtual bool is_branching (const Node& node) const = 0;

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
class NullTerminator :  public Terminator
{
public:

    /** Singleton method.
     * \return a singleton instance.
     */
    static Terminator& singleton()
    {
        static Graph dummy;
        static NullTerminator instance(dummy); return instance;
    }

    /** \copydoc Terminator::isEnabled */
    virtual bool isEnabled () const { return false; }

    /** \copydoc Terminator::mark */
    virtual void mark      (const Edge& edge) {}

    /** \copydoc Terminator::is_marked */
    virtual bool is_marked (const Edge& edge)  const  { return false; };

    /** \copydoc Terminator::mark(const gatb::core::debruijn::impl::Node&)  */
    virtual void mark      (const Node& node) {}

    /** \copydoc Terminator::is_marked(const gatb::core::debruijn::impl::Node&) const */
    virtual bool is_marked (const Node& node)  const  { return false; }

    /** \copydoc Terminator::is_marked_branching */
    virtual bool is_marked_branching (const Node& node) const { return false; }

    /** \copydoc Terminator::is_branching */
    virtual bool is_branching (const Node& node) const { return false; }

    /** \copydoc Terminator::reset */
    virtual void reset () {}

    /** \copydoc Terminator::dump */
    virtual void dump () {}

//private:

    NullTerminator (const Graph& graph) : Terminator(graph)  {}
};

/********************************************************************************/

/** \brief Implementation of Terminator that marks branching nodes.
 */
class BranchingTerminator :  public Terminator
{
public:

    typedef unsigned short int  Value;

    /** Constructor
     * \param[in] graph : the graph */
    BranchingTerminator (const Graph& graph);

    /** Destructor. */
    ~BranchingTerminator();

    /** \copydoc Terminator::mark */
    void mark      (const Edge& edge);

    /** \copydoc Terminator::is_marked */
    bool is_marked (const Edge& edge)  const;

    /** \copydoc Terminator::mark(const gatb::core::debruijn::impl::Node&)  */
    void mark      (const Node& node);

    /** \copydoc Terminator::is_marked(const gatb::core::debruijn::impl::Node&) const */
    bool is_marked (const Node& node) const ;

    /** \copydoc Terminator::is_marked_branching */
    bool is_marked_branching (const Node& node) const ;

    /** \copydoc Terminator::is_branching */
    bool is_branching (const Node& node) const ;

    /** \copydoc Terminator::reset */
    void reset();

    /** \copydoc Terminator::dump */
    void dump ();

private:

    bool is_indexed (const Node& node) const ;

    AssocSet<Node::Value, Value> branching_kmers;

    int getDelta (const Edge& edge) const;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_TOOLS_TERMINATOR_HPP_ */

