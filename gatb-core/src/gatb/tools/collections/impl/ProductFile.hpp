/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file Collection.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds interfaces related to the Collection interface
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_PRODUCT_FILE_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_PRODUCT_FILE_HPP_

/********************************************************************************/

#include <gatb/tools/collections/impl/Product.hpp>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

class ProductFileFactory
{
public:

    /** */
    static Product<ProductFileFactory>* createProduct (const std::string& name, bool autoRemove)
    {
        return new Product<ProductFileFactory> (name, autoRemove);
    }

    /** */
    static Group<ProductFileFactory>* createGroup (dp::ICell* parent, const std::string& name)
    {
        return new Group<ProductFileFactory> (parent, name);
    }

    /** */
    template<typename Type>
    static Partition<ProductFileFactory,Type>* createPartition (dp::ICell* parent, const std::string& name, size_t nb)
    {
        return new Partition<ProductFileFactory,Type> (parent, name, nb);
    }

    /** */
    template<typename Type>
    static CollectionNode<Type>* createCollection (dp::ICell* parent, const std::string& name, system::ISynchronizer* synchro)
    {
        /** We define the full qualified id of the current collection to be created. */
        std::stringstream ss;  ss << parent->getFullId().c_str() << "." << name;

        return new CollectionNode<Type> (parent, name, new CollectionFile<Type>(ss.str()));
    }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_PRODUCT_FILE_HPP_ */
