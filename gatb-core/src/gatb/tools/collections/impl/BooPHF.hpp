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

/** \file BooPHF.hpp
 *  \brief Minimal Perfect Hash Function from Guillaume
 */

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_BOOPHF_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_BOOPHF_HPP_

/********************************************************************************/

#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>

#include <BooPHF/BooPHF.h>

// let's use emphf base_hash for hashing elements here
#include <emphf/base_hash.hpp>

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
 * It uses BooPHF for the implementation and is most a wrapper between BooPHF and
 * GATB-CORE concepts.
 */
template<typename Key,typename Adaptator, class Progress>
class BooPHF : public system::SmartPointer
{
private:

    // a hash wrapper that calls emphf's hasher to produce, given an element, a single hash value for BooPHF
    class hasher_t
    {
        typedef emphf::jenkins64_hasher BaseHasher;
        BaseHasher emphf_hasher;
        Adaptator adaptor;
            
       public:
        hasher_t(){
            std::mt19937_64 rng(37); // deterministic seed
            emphf_hasher = BaseHasher::generate(rng);
        }

 

        uint64_t operator ()  (const Key& key, uint64_t seed = 0) const  {  
                if (seed != 0x33333333CCCCCCCCULL)
                    return std::get<0>(emphf_hasher(adaptor(key)));  
                return std::get<2>(emphf_hasher(adaptor(key)));   
                // this is a big hack, because I'm lazy. 
                // I wanted to return two different hashes depending on how boophf calls it
                // since I contrl BooPHF code's, I know it calls this function with 0x33333333CCCCCCCCULL as the second seed.
                }
    };

    typedef boomphf::mphf<  Key, hasher_t  > boophf_t;

public:

    /** Template specialization.  */
    static const bool enabled = true;

    /** Definition of a hash value. */
    typedef u_int64_t Code;

    /** Constructor. */
    BooPHF () : isBuilt(false), nbKeys(0)  {}

    /** Build the hash function from a set of items.
     * \param[in] iterable : keys iterator
     * \param[in] progress : object that listens to the event of the algorithm */
    void build (tools::collections::Iterable<Key>* iterable, int nbThreads, tools::dp::IteratorListener* progress=0)
    {
        if (isBuilt==true) { throw system::Exception ("MFHP: built already done"); }

        /** We create an iterator from the iterable. */
        tools::dp::Iterator<Key>* iter = iterable->iterator();
        LOCAL (iter);

        size_t nbElts = iterable->getNbItems();

        iterator_wrapper kmers (iter); // TODO use EMPHF's to prevent code duplication, or actually, put it in MPHFWrapper.

		bool withprogress = true;

		if (progress==0)
			withprogress = false;
		

        bphf =  boophf_t(nbElts, kmers, nbThreads,1.0, withprogress);

        isBuilt = true;
        nbKeys  = iterable->getNbItems();
    }

    /** Returns the hash code for the given key. WARNING : default implementation here will
     * throw an exception.
     * \param[in] key : the key to be hashed
     * \return the hash value. */
    Code operator () (const Key& key)
    {
        return bphf.lookup (key);
    }

    /** Returns the number of keys.
     * \return keys number */
    size_t size() const { return bphf.nbKeys(); }

    /** Load hash function from a collection*/
    size_t load (tools::storage::impl::Group& group, const std::string& name)
    {
        /** We need an input stream for the given collection given by group/name. */
        tools::storage::impl::Storage::istream is (group, name);
        /** We load the emphf object from the input stream. */
		bphf =  boophf_t();
        bphf.load (is);
        return size();
    }

    /** Save hash function to a collection
     * \return the number of bytes of the saved data. */
    size_t save (tools::storage::impl::Group& group, const std::string& name)
    {
        /** We need an output stream for the given collection given by group/name. */
        tools::storage::impl::Storage::ostream os (group, name);
        /** We save the emphf object to the output stream. */
        bphf.save (os);
        /** We set the number of keys as an attribute of the group. */
        group.addProperty ("nb_keys", misc::impl::Stringify().format("%d",nbKeys)); // FIXME: maybe overflow here
        return os.tellp();
    }

private:

    boophf_t  bphf;
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
        // noncopyble // FIXME: made it copyable because boophf needed it; need to see if it's correct
        //iterator_wrapper(iterator_wrapper const&);
        //iterator_wrapper& operator=(iterator_wrapper const&);
        tools::dp::Iterator<Key>* iterator;
    };
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_BOOPHF_HPP_ */
