/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file ProductFile.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds interfaces related to the Collection interface
 */

#ifndef _GATB_CORE_TOOLS_STORAGE_IMPL_PRODUCT_FILE_HPP_
#define _GATB_CORE_TOOLS_STORAGE_IMPL_PRODUCT_FILE_HPP_

/********************************************************************************/

#include <cassert>
#include <gatb/tools/storage/impl/CollectionFile.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace storage   {
namespace impl      {
/********************************************************************************/

class ProductFileFactory
{
public:

    /** */
    static Product* createProduct (const std::string& name, bool deleteIfExist, bool autoRemove)
    {
        return new Product (PRODUCT_FILE, name, autoRemove);
    }

    /** */
    static Group* createGroup (ICell* parent, const std::string& name)
    {
        ICell* root = ICell::getRoot (parent);
        Product* product = dynamic_cast<Product*> (root);
        assert (product != 0);

        return new Group (product->getFactory(), parent, name);
    }

    /** */
    template<typename Type>
    static Partition<Type>* createPartition (ICell* parent, const std::string& name, size_t nb)
    {
        ICell* root = ICell::getRoot (parent);
        Product* product = dynamic_cast<Product*> (root);
        assert (product != 0);

        Partition<Type>* result = new Partition<Type> (product->getFactory(), parent, name, nb);
        return result;
    }

    /** */
    template<typename Type>
    static CollectionNode<Type>* createCollection (ICell* parent, const std::string& name, system::ISynchronizer* synchro)
    {
        /** We define the full qualified id of the current collection to be created. */
        std::string actualName = std::string("tmp.") + name;

        ICell* root = ICell::getRoot (parent);
        Product* product = dynamic_cast<Product*> (root);
        assert (product != 0);

        return new CollectionNode<Type> (product->getFactory(), parent, name, new CollectionFile<Type>(actualName));
    }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_STORAGE_IMPL_PRODUCT_FILE_HPP_ */
