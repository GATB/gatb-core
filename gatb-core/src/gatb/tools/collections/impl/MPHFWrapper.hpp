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

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_MPHFWRAPPER_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_MPHFWRAPPER_HPP_

/********************************************************************************/

#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>
#include <gatb/tools/misc/api/Enums.hpp> // for MPHFKind

#include <gatb/tools/collections/impl/EMPHF.hpp> 
#include <gatb/tools/collections/impl/BooPHF.hpp> 

/********************************************************************************/
namespace gatb        {
namespace core        {
namespace tools       {
namespace collections {
namespace impl        {
/********************************************************************************/

/** \brief Minimal Perfect Hash Function
 *
 * that's a wrapper for emphf or boophf
 */
template<typename Key,typename Adaptator, class Progress>
class MPHF<Key,Adaptator,Progress,true> : public system::SmartPointer
{
private:

public:

    /** Template specialization.  */
    static const bool enabled = true;

    /** Definition of a hash value. */
    typedef u_int64_t Code;

    tools::misc::MPHFKind mphfKind;
    
    /** Constructor. */
    MPHF (tools::misc::MPHFKind mphfKind) : mphfKind(mphfKind)   {}
    
    /** Constructor. */
    MPHF () : mphfKind(tools::misc::MPHF_BOOPHF /* boophf is best mphf; so it deserves default*/)   {}

    /** Build the hash function from a set of items.
     * \param[in] iterable : keys iterator
     * \param[in] progress : object that listens to the event of the algorithm */
    void build (tools::collections::Iterable<Key>* iterable, int nbThreads = 1, tools::dp::IteratorListener* progress=0)
    {
        if (mphfKind == tools::misc::MPHF_EMPHF)
            emphf.build(iterable, nbThreads, progress);
        else
        {
            if (mphfKind == tools::misc::MPHF_BOOPHF)
                boophf.build(iterable, nbThreads, progress);
            else
               std::cout << "Error: building MPHF of wrong kind (debug: " << (unsigned int)mphfKind << ")" << std::endl;
        }
    }

    /** Returns the hash code for the given key. WARNING : default implementation here will
     * throw an exception.
     * \param[in] key : the key to be hashed
     * \return the hash value. */
    Code operator () (const Key& key)
    {
        if (mphfKind == tools::misc::MPHF_EMPHF)
            return emphf(key);
        if (mphfKind == tools::misc::MPHF_BOOPHF)
            return boophf(key);
        return 0;
    }

    /** Returns the number of keys.
     * \return keys number */
    size_t size() const { 
        if (mphfKind == tools::misc::MPHF_EMPHF)
            return emphf.size(); 
        if (mphfKind == tools::misc::MPHF_BOOPHF)
            return boophf.size(); 
        
        std::cout << "Error: size of MPHF of wrong kind (debug: " << (unsigned int)mphfKind << ")" << std::endl;
        return 0;
    }

    /** Load hash function from a collection*/
    size_t load (tools::storage::impl::Group& group, const std::string& name)
    {
        if (mphfKind == tools::misc::MPHF_EMPHF)
            return emphf.load(group,name); 
        if (mphfKind == tools::misc::MPHF_BOOPHF)
            return boophf.load(group,name); 

        std::cout << "Error: loading MPHF of wrong kind (debug: " << (unsigned int)mphfKind << ")" << std::endl;
        return 0;
    }

    /** Save hash function to a collection
     * \return the number of bytes of the saved data. */
    size_t save (tools::storage::impl::Group& group, const std::string& name)
    {
        if (mphfKind == tools::misc::MPHF_EMPHF)
            return emphf.save(group,name); 
        if (mphfKind == tools::misc::MPHF_BOOPHF)
            return boophf.save(group,name); 
        
        std::cout << "Error: loading MPHF of wrong kind (debug: " << (unsigned int)mphfKind << ")" << std::endl;
        return 0;
    }

private:
    
    /** we will use alternatively one or the other;
     * just having them as variable doesn't use any memory as long as build() isn't called*/

    EMPHF<Key,Adaptator,Progress> emphf;
    BooPHF<Key,Adaptator,Progress> boophf;

};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_EMPHF_HPP_ */
