#ifdef WITH_MPHF
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

#include <gatb/tools/storage/impl/StorageTools.hpp>


// for mphf
#include <iostream>
#include <fstream>

#include <emphf/common.hpp>
#include <emphf/mphf.hpp>
#include <emphf/base_hash.hpp>

using namespace emphf;

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace storage   {
namespace impl      {
/********************************************************************************/


    template<typename BaseHasher>   void StorageTools::saveEMPHF (Group& group, const std::string& name, void* mphf_void)
    {
#ifdef WITH_MPHF
        // TODO: I couldn't bother converting the ostream returned by mphf to the "collection" of GATB (whatever that is),
        // so I'm hacking my way through
        
        /*collections::Collection<math::NativeInt8>* mphfCollection = & group.getCollection<math::NativeInt8> (name);
        
          mphfCollection->insert ((math::NativeInt8*)mphf->save(something like a converter from an ostream to an array of chars), mphf->size());

        std::stringstream ss1;  ss1 <<  mphf->size();
        bloomCollection->addProperty ("size",    ss1.str());
        bloomCollection->addProperty ("name",    name);*/

        logger() << "Saving mphf to disk" << std::endl;
        mphf<BaseHasher>* mphf_propercast = static_cast< emphf::mphf<BaseHasher>* >(mphf_void);
        std::ofstream os(name, std::ios::binary);
        mphf_propercast->save(os);
#endif

    }

    /** */
    template<typename BaseHasher>     void*  StorageTools::loadEMPHF (Group& group, const std::string& name)
    {
        // TODO same as function right above
        //
#ifdef WITH_MPHF
        mphf<BaseHasher> *mphf_propercast = new mphf<BaseHasher>();
        logger() << "Loading mphf from disk" << std::endl;
        std::ifstream is(name, std::ios::binary);
        mphf_propercast->load(is);
        return static_cast<void*>(mphf_propercast);
#endif
    }

    // tell compiler to instantiate those templates, else there will be undefined reference problems
    // and I don't want to put that code in the .hpp, to avoid including emphf includes in the hpp, I recall it posed problems
    template  void StorageTools::saveEMPHF<  emphf::jenkins64_hasher >(Group& group, const std::string& name, void* mphf_void);
    template  void* StorageTools::loadEMPHF< emphf::jenkins64_hasher >(Group& group, const std::string& name);


/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
#endif 

