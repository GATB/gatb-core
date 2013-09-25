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

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_PRODUCT_HDF5_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_PRODUCT_HDF5_HPP_

/********************************************************************************/

#include <gatb/tools/collections/impl/Product.hpp>
#include <gatb/tools/collections/impl/CollectionHDF5.hpp>
#include <gatb/system/impl/System.hpp>
#include <hdf5.h>
#include <typeinfo>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

class ProductHDF5Factory
{
public:

    /** */
    static Product<ProductHDF5Factory>* createProduct (const std::string& name, bool autoRemove)
    {
        return new ProductHDF5 (name, autoRemove);
    }

    /** */
    static Group<ProductHDF5Factory>* createGroup (dp::ICell* parent, const std::string& name)
    {
        dp::ICell* root = dp::ICell::getRoot (parent);

        Product<ProductHDF5Factory>* p1 = dynamic_cast<Product<ProductHDF5Factory>*> (root);

        ProductHDF5* product = dynamic_cast<ProductHDF5*> (p1);
        assert (product != 0);

        /** We create the instance. */
        Group<ProductHDF5Factory>* result = new Group<ProductHDF5Factory> (parent, name);

        std::string actualName = result->getFullId('/');

        /** We create the HDF5 group if needed. */
        htri_t doesExist = H5Lexists (product->getId(), actualName.c_str(), H5P_DEFAULT);
        if (doesExist <= 0)
        {
            hid_t group = H5Gcreate (product->getId(), actualName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Gclose (group);
        }

        /** We return the result. */
        return result;
    }

    /** */
    template<typename Type>
    static Partition<ProductHDF5Factory,Type>* createPartition (dp::ICell* parent, const std::string& name, size_t nb)
    {
        dp::ICell* loop = dp::ICell::getRoot (parent);
        ProductHDF5* product = dynamic_cast<ProductHDF5*> (loop);

        std::string actualName = parent->getFullId('/') + "/" + name;

        /** We create the HDF5 group if needed. */
        htri_t doesExist = H5Lexists (product->getId(), actualName.c_str(), H5P_DEFAULT);
        if (doesExist <= 0)
        {
            hid_t group = H5Gcreate (product->getId(), actualName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Gclose (group);
        }

        return new Partition<ProductHDF5Factory,Type> (parent, name, nb);
    }

    /** */
    template<typename Type>
    static CollectionNode<Type>* createCollection (dp::ICell* parent, const std::string& name, system::ISynchronizer* synchro)
    {
#if 1
        synchro = GlobalSynchro::singleton();
#endif
        dp::ICell* loop=0;
        for (loop=parent ; loop->getParent() != 0;  loop=loop->getParent())  {}
        ProductHDF5* product = dynamic_cast<ProductHDF5*> (loop);

        std::string actualName = parent->getFullId('/') + "/" + name;

        return new CollectionNode<Type> (parent, name, new CollectionHDF5<Type>(product->getId(), actualName, synchro));
    }

private:

    class GlobalSynchro
    {
    public:
        static system::ISynchronizer* singleton()
        {
            static GlobalSynchro instance;
            return instance.synchro;
        }

    private:
        GlobalSynchro ()  { synchro = system::impl::System::thread().newSynchronizer(); }
        ~GlobalSynchro () { if (synchro)  { delete synchro; } }
        system::ISynchronizer* synchro;
    };

    /** */
    class ProductHDF5 : public Product<ProductHDF5Factory>
    {
    public:
        ProductHDF5 (const std::string& name, bool autoRemove) : Product<ProductHDF5Factory> (name, autoRemove), _fileId(0), _name(name)
        {
            /** We create a new file using default properties. */
            if (system::impl::System::file().doesExist(getActualName()))
            {
                _fileId = H5Fopen (getActualName().c_str(), H5P_DEFAULT, H5P_DEFAULT);
            }
            else
            {
                _fileId = H5Fcreate (getActualName().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            }
        }

        virtual ~ProductHDF5 ()
        {
            if (_autoRemove)  { remove(); }
            H5Fclose(_fileId);
        }

        hid_t getId ()  { return _fileId; }

        void remove ()
        {
            system::impl::System::file().remove (getActualName());
        }

    private:
        hid_t       _fileId;
        std::string _name;

        std::string getActualName ()  { return _name + ".h5"; }

    };
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_PRODUCT_HDF5_HPP_ */
