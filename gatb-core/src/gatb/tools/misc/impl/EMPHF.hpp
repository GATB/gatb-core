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

/** \file EMPHF.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Minimal Perfect Hash Function
 */

#ifdef WITH_MPHF

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_EMPHF_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_EMPHF_HPP_

/********************************************************************************/

#include <gatb/tools/misc/api/MPHF.hpp>

#include <emphf/common.hpp>
#include <emphf/mphf.hpp>
#include <emphf/base_hash.hpp>
#include <emphf/mmap_memory_model.hpp>
#include <emphf/hypergraph_sorter_scan.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
/********************************************************************************/

/** \brief Minimal Perfect Hash Function
 *
 * This is a specialization of the MPHF<Key,Adaptor,exist> class for exist=true.
 * It uses EMPHF for the implementation and is most a wrapper between EMPHF and
 * GATB-CORE concepts.
 */
template<typename Key,typename Adaptator>
class MPHF<Key,Adaptator,true>
{
private:
    // adapted from compute_mphf_scan_mmap.cpp
    typedef emphf::hypergraph_sorter_scan<uint32_t, emphf::mmap_memory_model> HypergraphSorter32;
    typedef emphf::hypergraph_sorter_scan<uint64_t, emphf::mmap_memory_model> HypergraphSorter64;
    typedef emphf::jenkins64_hasher BaseHasher;
    typedef emphf::mphf<BaseHasher> mphf_t;

public:

    /** */
    static const bool enabled = true;

    /** Definition of a hash value. */
    typedef u_int64_t Value;

    /** Constructor. */
    MPHF () : isBuilt(false)  {}

    /** Constructor. */
    MPHF (tools::dp::Iterator<Key>* iterator, size_t nbElts)  : isBuilt(false)  {  build (iterator, nbElts);  }

    /** Constructor. */
    MPHF (tools::collections::Iterable<Key>* iterable)  : isBuilt(false)
    {
        tools::dp::Iterator<Key>* iterator = iterable->iterator(); LOCAL(iterator);
        build (iterator, iterable->getNbItems());
    }

    /** Build the hash function from a set of items.
     * \param[in] iterator : keys iterator
     * \param[in] nbItems : number of keys (if known) */
    void build (tools::dp::Iterator<Key>* iterator, size_t nbElts)
    {
        if (isBuilt==true) { throw system::Exception ("MFHP: built already done"); }

        iterator_wrapper kmers (iterator);

        size_t max_nodes = (size_t(std::ceil(double(nbElts) * 1.23)) + 2) / 3 * 3;
        if (max_nodes >= uint64_t(1) << 32)
        {
            HypergraphSorter64 sorter;
            mphf_t(sorter, nbElts, kmers, adaptor).swap(mphf);
        }
        else
        {
            HypergraphSorter32 sorter;
            mphf_t(sorter, nbElts, kmers, adaptor).swap(mphf);
        }

        isBuilt = true;
    }

    /** Returns the hash code for the given key. WARNING : default implementation here will
     * throw an exception.
     * \param[in] key : the key to be hashed
     * \return the hash value. */
    Value operator () (const Key& key)
    {
        return mphf.lookup (key, adaptor);
    }

    /** Returns the number of keys.
     * \return keys number */
    size_t size() const { return mphf.size(); }

private:

    mphf_t    mphf;
    Adaptator adaptor;
    bool      isBuilt;

private:
    class iterator_adaptator : public std::iterator<std::forward_iterator_tag, const Key>
    {
    public:
        iterator_adaptator()  : iterator(0), pos(0) {}

        iterator_adaptator(tools::dp::Iterator<Key>* iterator)  : iterator(iterator), pos(0)  {  iterator->first();  }

        Key const& operator*()  {  return iterator->item();  }

        iterator_adaptator& operator++()
        {
            iterator->next();
            pos++;
            if (iterator->isDone())
            {
                iterator = nullptr;
                pos = 0;
            }
            return *this;
        }

        friend bool operator==(iterator_adaptator const& lhs, iterator_adaptator const& rhs)
        {
            if (!lhs.iterator || !rhs.iterator)  {  if (!lhs.iterator && !rhs.iterator) {  return true; } else {  return false;  } }
            return rhs.pos == lhs.pos;
        }

        friend bool operator!=(iterator_adaptator const& lhs, iterator_adaptator const& rhs)  {  return !(lhs == rhs);  }

    private:
        tools::dp::Iterator<Key>* iterator;
        unsigned long pos;
    };

    class iterator_wrapper
    {
    public:
        iterator_wrapper (tools::dp::Iterator<Key>* iterator) : iterator(iterator) {}

        iterator_adaptator begin() const  {  return iterator_adaptator (iterator); }
        iterator_adaptator end  () const  {  return iterator_adaptator ();         }
        size_t        size () const  {  return 0;                        }

    private:
        // noncopyble
        iterator_wrapper(iterator_wrapper const&);
        iterator_wrapper& operator=(iterator_wrapper const&);
        tools::dp::Iterator<Key>* iterator;
    };
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_EMPHF_HPP_ */

#endif /* FOO */
