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

/** \file StorageTools.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds interfaces related to the Collection interface
 */

#ifndef _GATB_CORE_TOOLS_STORAGE_IMPL_STORAGE_TOOLS_HPP_
#define _GATB_CORE_TOOLS_STORAGE_IMPL_STORAGE_TOOLS_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Collection.hpp>
#include <gatb/tools/math/NativeInt8.hpp>

#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/collections/impl/Bloom.hpp>
#include <gatb/tools/collections/impl/ContainerSet.hpp>


/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace storage   {
namespace impl      {
/********************************************************************************/

class StorageTools
{
public:

    /** why is there code in headers? because else we'd have to instantiate the template, see http://stackoverflow.com/questions/8752837/undefined-reference-to-template-class-constructor */

    /** */
    static StorageTools& singleton() { static StorageTools instance; return instance; }


    /** */
    template<typename T>  void saveContainer (Group& group, const std::string& name, collections::Collection<T>* collection)
    {
        collections::Collection<T>* storageCollection = & group.getCollection<T> (name);

        tools::dp::Iterator<T>* it = collection->iterator();   LOCAL(it);
        for (it->first(); !it->isDone(); it->next())  {  storageCollection->insert (it->item());  }
        storageCollection->flush ();
    }

    /** */
    template<typename T>  collections::Container<T>*  loadContainer (Group& group, const std::string& name)
    {
        collections::Collection<T>*  storageCollection = & group.getCollection<T> (name);
        return new collections::impl::ContainerSet<T> (storageCollection->iterator());
    }

    /** */
    template<typename T>  void saveBloom (Group& group, const std::string& name, collections::impl::IBloom<T>* bloom)
    {
        collections::Collection<math::NativeInt8>* bloomCollection = & group.getCollection<math::NativeInt8> (name);

        bloomCollection->insert ((math::NativeInt8*)bloom->getArray(), bloom->getSize());

        std::stringstream ss1;  ss1 <<  bloom->getBitSize();
        std::stringstream ss2;  ss2 <<  bloom->getNbHash();

        bloomCollection->addProperty ("size",    ss1.str());
        bloomCollection->addProperty ("nb_hash", ss2.str());
        bloomCollection->addProperty ("type",    bloom->getName());
    }

    /** */
    template<typename T>  collections::impl::IBloom<T>*  loadBloom (Group& group, const std::string& name)
    {
        /** We retrieve the raw data buffer for the Bloom filter. */
        tools::collections::Collection<tools::math::NativeInt8>* bloomArray = & group.getCollection<tools::math::NativeInt8> (name);

        /** We create the Bloom fiter. */
        tools::collections::impl::IBloom<T>* bloom = tools::collections::impl::BloomFactory::singleton().createBloom<T> (
            bloomArray->getProperty("type"),
            bloomArray->getProperty("size"),
            bloomArray->getProperty("nb_hash")
        );

        /** We set the bloom with the provided array given as an iterable of NativeInt8 objects. */
        bloomArray->getItems ((tools::math::NativeInt8*&)bloom->getArray());

        /** We return the result. */
        return bloom;
    }

#ifdef WITH_MPHF
    /** */
    template<typename BaseHasher>   void saveEMPHF (Group& group, const std::string& name, void* mphf);

    /** */
    template<typename BaseHasher>     void*  loadEMPHF (Group& group, const std::string& name);
#endif

};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_STORAGE_IMPL_STORAGE_TOOLS_HPP_ */
