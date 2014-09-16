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

/** \file MapMPHF.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Minimal Perfect Hash Function
 */

#ifndef _GATB_CORE_TOOLS_COLLECTION_MAP_MPHF_HPP_
#define _GATB_CORE_TOOLS_COLLECTION_MAP_MPHF_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/collections/impl/MPHF.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <vector>

/********************************************************************************/
namespace gatb        {
namespace core        {
namespace tools       {
namespace collections {
namespace impl        {
/********************************************************************************/

/** */
template <class Key, class Value, class Adaptator=AdaptatorDefault<Key> >
class MapMPHF : public system::SmartPointer
{
public:

    /** Hash type. */
    typedef MPHF<Key, Adaptator> Hash;

    /** Constant telling whether the feature is enabled or not. */
    static const bool enabled = Hash::enabled;

    /** */
    MapMPHF ()  {}

    /** Build the hash function from a set of items.
     * \param[in] iterator : keys iterator
     * \param[in] nbItems : number of keys (if known) */
    void build (tools::collections::Iterable<Key>& keys, tools::dp::IteratorListener* progress=0)
    {
        /** We build the hash function. */
        hash.build (&keys, progress);

        /** We resize the vector of Value objects. */
        data.resize (keys.getNbItems());
    }

    /** Save hash function to a collection
     * \return the number of bytes of the saved data. */
    size_t save (tools::storage::impl::Group& group, const std::string& name)  {  return hash.save (group, name);  }

    /** Load hash function from a collection*/
    void load (tools::storage::impl::Group& group, const std::string& name)
    {
        /** We load the hash function. */
        size_t nbKeys = hash.load (group, name);

        /** We resize the vector of Value objects. */
        data.resize (nbKeys);
    }

    /** Get the value for a given key
     * \param[in] key : the key
     * \return the value associated to the key. */
    Value& operator[] (const Key& key)  { return data[hash(key)]; }

    /** Get the value for a given index
     * \param[in] index : the index
     * \return the value associated to the index. */
    Value& at (typename Hash::Code code)  { return data[code]; }

    /** Get the hash code of the given key. */
    typename Hash::Code getCode (const Key& key) { return hash(key); }

    /** Get the number of keys.
     * \return keys number. */
    size_t size() const { return hash.size(); }

private:

    Hash               hash;
    std::vector<Value> data;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTION_MAP_MPHF_HPP_ */
