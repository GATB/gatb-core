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

/** \file StorageFile.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds interfaces related to the Collection interface
 */

#ifndef _GATB_CORE_TOOLS_STORAGE_IMPL_STORAGE_FILE_HPP_
#define _GATB_CORE_TOOLS_STORAGE_IMPL_STORAGE_FILE_HPP_

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

class StorageFileFactory
{
public:

    /** */
    static Storage* createStorage (const std::string& name, bool deleteIfExist, bool autoRemove)
    {
        return new Storage (STORAGE_FILE, name, autoRemove);
    }

    /** */
    static bool exists (const std::string& name)
    {
        return false;
    }

    /** */
    static Group* createGroup (ICell* parent, const std::string& name)
    {
        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);

        return new Group (storage->getFactory(), parent, name);
    }

    /** */
    template<typename Type>
    static Partition<Type>* createPartition (ICell* parent, const std::string& name, size_t nb)
    {
        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);

        Partition<Type>* result = new Partition<Type> (storage->getFactory(), parent, name, nb);
        return result;
    }

    /** */
    template<typename Type>
    static CollectionNode<Type>* createCollection (ICell* parent, const std::string& name, system::ISynchronizer* synchro)
    {
        /** We define the full qualified id of the current collection to be created. */
        std::string actualName = std::string("tmp.") + name;

        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);

        return new CollectionNode<Type> (storage->getFactory(), parent, name, new CollectionFile<Type>(actualName));
    }
};



class StorageGzFileFactory
{
public:
    
    /** */
    static Storage* createStorage (const std::string& name, bool deleteIfExist, bool autoRemove)
    {
        return new Storage (STORAGE_GZFILE, name, autoRemove);
    }
    
    /** */
    static bool exists (const std::string& name)
    {
        return false;
    }

    /** */
    static Group* createGroup (ICell* parent, const std::string& name)
    {
        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);
        
        return new Group (storage->getFactory(), parent, name);
    }
    
    /** */
    template<typename Type>
    static Partition<Type>* createPartition (ICell* parent, const std::string& name, size_t nb)
    {
        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);
        
        Partition<Type>* result = new Partition<Type> (storage->getFactory(), parent, name, nb);
        return result;
    }
    
    /** */
    template<typename Type>
    static CollectionNode<Type>* createCollection (ICell* parent, const std::string& name, system::ISynchronizer* synchro)
    {
        /** We define the full qualified id of the current collection to be created. */
        std::string actualName = std::string("tmp.") + name;
        
        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);
        
        return new CollectionNode<Type> (storage->getFactory(), parent, name, new CollectionGzFile<Type>(actualName));
    }
};

class StorageSortedFactory
{
public:
    
    /** */
    static Storage* createStorage (const std::string& name, bool deleteIfExist, bool autoRemove)
    {
        return new Storage (STORAGE_COMPRESSED_FILE, name, autoRemove);
    }
    
    /** */
    static bool exists (const std::string& name)
    {
        return false;
    }

    /** */
    static Group* createGroup (ICell* parent, const std::string& name)
    {
        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);
        
        return new Group (storage->getFactory(), parent, name);
    }
    
    /** */
    template<typename Type>
    static Partition<Type>* createPartition (ICell* parent, const std::string& name, size_t nb)
    {
        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);
        
        Partition<Type>* result = new Partition<Type> (storage->getFactory(), parent, name, nb);
        return result;
    }
    
    /** */
    template<typename Type>
    static CollectionNode<Type>* createCollection (ICell* parent, const std::string& name, system::ISynchronizer* synchro)
    {
        /** We define the full qualified id of the current collection to be created. */
        std::string actualName = std::string("tmp.") + name;
        
        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);
        
        return new CollectionNode<Type> (storage->getFactory(), parent, name, new CollectionCountFile<Type>(actualName));
    }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_STORAGE_IMPL_STORAGE_FILE_HPP_ */
