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

/** \file StorageHDF5.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds interfaces related to the Collection interface
 */

#ifndef _GATB_CORE_TOOLS_STORAGE_IMPL_STORAGE_HDF5_HPP_
#define _GATB_CORE_TOOLS_STORAGE_IMPL_STORAGE_HDF5_HPP_

/********************************************************************************/

#include <gatb/tools/storage/impl/CollectionHDF5.hpp>
#include <gatb/system/impl/System.hpp>
#include <hdf5/hdf5.h>
#include <typeinfo>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace storage   {
namespace impl      {
/********************************************************************************/

class StorageHDF5Factory
{
public:

    /** */
    static Storage* createStorage (const std::string& name, bool deleteIfExist, bool autoRemove)
    {
        return new StorageHDF5 (STORAGE_HDF5, name, deleteIfExist, autoRemove);
    }

    /** */
    static bool exists (const std::string& name)
    {
        H5Eset_auto (0, NULL, NULL);

        bool result = false;
        {
            hid_t id = H5Fopen (name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            if (id > 0)  {  H5Fclose(id);  result = true;  }
        }
        {
            hid_t id = H5Fopen ((name+".h5").c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            if (id > 0)  {  H5Fclose(id);  result = true;  }
        }
        return result;
    }

    /** */
    static Group* createGroup (ICell* parent, const std::string& name)
    {
        StorageHDF5* storage = dynamic_cast<StorageHDF5*> (ICell::getRoot (parent));
        assert (storage != 0);

        /** We create the instance. */
        Group* result = new Group (storage->getFactory(), parent, name);

        /** We may need to create the HDF5 group. Empty name means root group, which is constructed by default. */
        if (name.empty() == false)
        {
            std::string actualName = result->getFullId('/');

            /** We create the HDF5 group if needed. */
            htri_t doesExist = H5Lexists (storage->getId(), actualName.c_str(), H5P_DEFAULT);

            if (doesExist <= 0)
            {
                hid_t group = H5Gcreate (storage->getId(), actualName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Gclose (group);
            }
        }

        /** We return the result. */
        return result;
    }

    /** */
    template<typename Type>
    static Partition<Type>* createPartition (ICell* parent, const std::string& name, size_t nb)
    {
        StorageHDF5* storage = dynamic_cast<StorageHDF5*> (ICell::getRoot (parent));
        assert (storage != 0);

        std::string actualName = parent->getFullId('/') + "/" + name;

        /** We create the HDF5 group if needed. */
        htri_t doesExist = H5Lexists (storage->getId(), actualName.c_str(), H5P_DEFAULT);
        if (doesExist <= 0)
        {
            hid_t group = H5Gcreate (storage->getId(), actualName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Gclose (group);
        }

        return new Partition<Type> (storage->getFactory(), parent, name, nb);
    }

    /** */
    template<typename Type>
    static CollectionNode<Type>* createCollection (ICell* parent, const std::string& name, system::ISynchronizer* synchro)
    {
#if 1
        synchro = GlobalSynchro::singleton();
#endif

        StorageHDF5* storage = dynamic_cast<StorageHDF5*> (ICell::getRoot (parent));
        assert (storage != 0);

        std::string actualName = parent->getFullId('/') + "/" + name;

        return new CollectionNode<Type> (storage->getFactory(), parent, name, new CollectionHDF5<Type>(storage->getId(), actualName, synchro));
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
    class StorageHDF5 : public Storage
    {
    public:
        StorageHDF5 (StorageMode_e mode, const std::string& name, bool deleteIfExist, bool autoRemove)
            : Storage (mode, name, autoRemove), _fileId(0), _name(name)
        {
            if (deleteIfExist)  {  system::impl::System::file().remove (getActualName());  }

            std::string usedName;

            /** We test the actual name. */
                 if (system::impl::System::file().doesExist(getActualName()))  { usedName = getActualName(); }
            else if (system::impl::System::file().doesExist(getName()))        { usedName = getName();       }

            if (usedName.empty()==false)
            {
                /** We open the existing file. */
                _fileId = H5Fopen (usedName.c_str(), H5P_DEFAULT, H5P_DEFAULT);
            }
            else
            {
                /** We create a new file using default properties. */
                _fileId = H5Fcreate (getActualName().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            }
        }

        virtual ~StorageHDF5 ()
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

        std::string getActualName () const  { return _name + ".h5"; }
        std::string getName       () const  { return _name;         }

    };
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_STORAGE_IMPL_STORAGE_HDF5_HPP_ */
