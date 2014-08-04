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

        return new GroupHDF5 (storage, parent, name);
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

    /************************************************************/
    class StorageHDF5 : public Storage
    {
    public:
        StorageHDF5 (StorageMode_e mode, const std::string& name, bool deleteIfExist, bool autoRemove)
            : Storage (mode, name, autoRemove), _fileId(0), _name(name)
        {
            if (deleteIfExist)  {  system::impl::System::file().remove (getActualName());  }

            /** We test the actual name exists in filesystem. */
            bool exists = system::impl::System::file().doesExist(getActualName());

            if (exists==true)
            {
                /** We open the existing file. */
                _fileId = H5Fopen (getActualName().c_str(), H5P_DEFAULT, H5P_DEFAULT);
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
        std::string _actualName;

        /** */
        std::string getActualName ()
        {
            /** We set the actual name at first call. */
            if (_actualName.empty())
            {
                _actualName = _name;
                /** We check whether the given name has a ".h5" suffix. */
                if (_name.rfind(".h5") == std::string::npos)  {  _actualName += ".h5";  }

            }
            return _actualName;
        }

        /** */
        std::string getName       () const  { return _name;         }

    };

    /************************************************************/
    class GroupHDF5 : public Group
    {
    public:
        GroupHDF5 (StorageHDF5* storage, ICell* parent, const std::string& name)
        : Group(storage->getFactory(),parent,name), _groupId(0)
        {
            /** We may need to create the HDF5 group. Empty name means root group, which is constructed by default. */
            if (name.empty() == false)
            {
                std::string actualName = this->getFullId('/');

                /** We create the HDF5 group if needed. */
                htri_t doesExist = H5Lexists (storage->getId(), actualName.c_str(), H5P_DEFAULT);

                if (doesExist <= 0)
                {
                    _groupId = H5Gcreate (storage->getId(), actualName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                }
                else
                {
                    _groupId = H5Gopen2 (storage->getId(), actualName.c_str(), H5P_DEFAULT);
                }
            }
            else
            {
                _groupId = H5Gopen2 (storage->getId(), "/", H5P_DEFAULT);
            }
        }

        /** */
        ~GroupHDF5()
        {
            /** We release the group handle. */
            H5Gclose(_groupId);
        }

        /** */
        void addProperty (const std::string& key, const std::string value)
        {
            hid_t datatype = H5Tcopy (H5T_C_S1);  H5Tset_size (datatype, H5T_VARIABLE);

            hsize_t dims = 1;
            hid_t space_id = H5Screate_simple (1, &dims, NULL);

            /** We create the attribute. */
            hid_t attrId = H5Acreate2 (_groupId, key.c_str(), datatype,  space_id, H5P_DEFAULT, H5P_DEFAULT);
            if (attrId >= 0)
            {
                /** We write the data. */
                const char* array[] = { value.c_str() };
                H5Awrite (attrId, datatype, &array);

                /** We close resources. */
                H5Aclose (attrId);
                H5Tclose (datatype);
                H5Sclose (space_id);
            }
        }

        /** */
        std::string getProperty (const std::string& key)
        {
            std::string result;
            herr_t status;

            hid_t datatype = H5Tcopy (H5T_C_S1);  H5Tset_size (datatype, H5T_VARIABLE);

            hid_t attrId = H5Aopen (_groupId, key.c_str(), H5P_DEFAULT);
            if (attrId >= 0)
            {
                hid_t space_id = H5Aget_space (attrId);

                hsize_t dims = 1;
                H5Sget_simple_extent_dims (space_id, &dims, NULL);
                char** rdata = (char **) malloc (dims * sizeof (char *));

                status = H5Aread (attrId, datatype, rdata);

                /** We set the result. */
                result.assign (rdata[0]);

                /** We release buffers. */
                status = H5Dvlen_reclaim (datatype, space_id, H5P_DEFAULT, rdata);
                free (rdata);

                /** We close resources. */
                H5Aclose (attrId);
                H5Tclose (datatype);
                H5Sclose (space_id);
            }

            return result;
        }

    private:
        hid_t _groupId;
    };
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_STORAGE_IMPL_STORAGE_HDF5_HPP_ */
