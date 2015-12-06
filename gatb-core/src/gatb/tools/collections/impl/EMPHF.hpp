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

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_EMPHF_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_EMPHF_HPP_

/********************************************************************************/

#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>

#include <emphf/common.hpp>
#include <emphf/base_hash.hpp>
#include <emphf/mmap_memory_model.hpp>

//#define USE_HEM 1 
// **warning**: HEM is buggy with k=21; it just is. I have not investigated why because we switched to boophf.
// otherwise, HEM, is faster to construct. will use more memory and might have slower queries.

#ifdef USE_HEM
#include <emphf/mphf_hem.hpp>
#else
#include <emphf/mphf.hpp>
#include <emphf/hypergraph_sorter_scan.hpp>
#endif
/********************************************************************************/
namespace gatb        {
namespace core        {
namespace tools       {
namespace collections {
namespace impl        {
/********************************************************************************/

/** \brief Minimal Perfect Hash Function
 *
 * This is a specialization of the MPHF<Key,Adaptor,exist> class for exist=true.
 * It uses EMPHF for the implementation and is most a wrapper between EMPHF and
 * GATB-CORE concepts.
 */
template<typename Key,typename Adaptator, class Progress>
class EMPHF : public system::SmartPointer
{
private:

#ifdef USE_HEM
    typedef emphf::mphf_hem<emphf::jenkins64_hasher> mphf_t; 
#else
    // adapted from compute_mphf_scan_mmap.cpp
    typedef emphf::hypergraph_sorter_scan<uint32_t, emphf::mmap_memory_model> HypergraphSorter32;
    typedef emphf::hypergraph_sorter_scan<uint64_t, emphf::mmap_memory_model> HypergraphSorter64;
    typedef emphf::jenkins64_hasher BaseHasher;
    typedef emphf::mphf<BaseHasher> mphf_t;
#endif
public:

    /** Template specialization.  */
    static const bool enabled = true;

    /** Definition of a hash value. */
    typedef u_int64_t Code;

    /** Constructor. */
    EMPHF () : isBuilt(false), nbKeys(0)  {}

    /** Build the hash function from a set of items.
     * \param[in] iterable : keys iterator
     * \param[in] progress : object that listens to the event of the algorithm */
    void build (tools::collections::Iterable<Key>* iterable, int nbThreads = 1, tools::dp::IteratorListener* progress=0)
    {
        if (isBuilt==true) { throw system::Exception ("MFHP: built already done"); }

        /** We create an iterator from the iterable. */
        tools::dp::Iterator<Key>* iter = iterable->iterator();
        LOCAL (iter);

        size_t nbElts = iterable->getNbItems();

        // a small fix, emphf for 2 nodes doesn't seem to work
        if (nbElts <= 2)  {  nbElts = 3;  }
        if (nbElts <= 3) { std::cout << "Warning: MPHF has a tiny amount of elements (" << nbElts << "), might not work correctly." << std::endl; }

        iterator_wrapper kmers (iter);

        // We may have no provided listener => use default one.
        if (progress==0)  { progress = new tools::dp::IteratorListener; }
        LOCAL (progress);

#ifdef USE_HEM
        emphf::mmap_memory_model mm;
        mphf_t(mm, nbElts, kmers, adaptor, progress).swap(mphf);
#else
        size_t max_nodes = (size_t(std::ceil(double(nbElts) * 1.23)) + 2) / 3 * 3;
        if (max_nodes >= uint64_t(1) << 32)
        {
            HypergraphSorter64 sorter;
            mphf_t(sorter, nbElts, kmers, adaptor, progress).swap(mphf);
        }
        else
        {
            HypergraphSorter32 sorter;
            mphf_t(sorter, nbElts, kmers, adaptor, progress).swap(mphf);
        }
#endif

        isBuilt = true;
        nbKeys  = iterable->getNbItems();
    }

    /** Returns the hash code for the given key. WARNING : default implementation here will
     * throw an exception.
     * \param[in] key : the key to be hashed
     * \return the hash value. */
    Code operator () (const Key& key)
    {
        return mphf.lookup (key, adaptor);
    }

    /** Returns the number of keys.
     * \return keys number */
    size_t size() const { return mphf.size(); }

    /** Load hash function from a collection*/
    size_t load (tools::storage::impl::Group& group, const std::string& name)
    {
        /** We need an input stream for the given collection given by group/name. */
        tools::storage::impl::Storage::istream is (group, name);
        /** We load the emphf object from the input stream. */
        mphf.load (is);
        /** We return the number of keys. */
        return mphf.size();
    }

    /** Save hash function to a collection
     * \return the number of bytes of the saved data. */
    size_t save (tools::storage::impl::Group& group, const std::string& name)
    {
        /** We need an output stream for the given collection given by group/name. */
        tools::storage::impl::Storage::ostream os (group, name);
        /** We save the emphf object to the output stream. */
        mphf.save (os);
        /** We set the number of keys as an attribute of the group. */
        group.addProperty ("nb_keys", misc::impl::Stringify().format("%d",nbKeys)); // FIXME: maybe overflow here
        return os.tellp();
    }

private:

    mphf_t    mphf;
    Adaptator adaptor;
    bool      isBuilt;
    size_t    nbKeys;

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
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_EMPHF_HPP_ */
