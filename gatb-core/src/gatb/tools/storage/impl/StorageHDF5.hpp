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
#include <gatb/tools/storage/impl/CollectionHDF5Patch.hpp>
#include <gatb/system/impl/System.hpp>
#include <hdf5/hdf5.h>
#include <sstream>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace storage   {
namespace impl      {
/********************************************************************************/

/** \brief Factory used for storage of kind STORAGE_HDF5
 */
class StorageHDF5Factory
{
public:

    /** Create a Storage instance.
     * \param[in] name : name of the instance to be created
     * \param[in] deleteIfExist : if the storage exits in file system, delete it if true.
     * \param[in] autoRemove : auto delete the storage from file system during Storage destructor.
     * \return the created Storage instance
     */
    static Storage* createStorage (const std::string& name, bool deleteIfExist, bool autoRemove)
    {
        return new StorageHDF5 (STORAGE_HDF5, name, deleteIfExist, autoRemove);
    }

    /** Tells whether or not a Storage exists in file system given a name
     * \param[in] name : name of the storage to be checked
     * \return true if the storage exists in file system, false otherwise.
     */
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

    /** Create a Group instance and attach it to a cell in a storage.
     * \param[in] parent : parent of the group to be created
     * \param[in] name : name of the group to be created
     * \return the created Group instance.
     */
    static Group* createGroup (ICell* parent, const std::string& name)
    {
        StorageHDF5* storage = dynamic_cast<StorageHDF5*> (ICell::getRoot (parent));
        assert (storage != 0);

        return new GroupHDF5 (storage, parent, name);
    }

    /** Create a Partition instance and attach it to a cell in a storage.
     * \param[in] parent : parent of the partition to be created
     * \param[in] name : name of the partition to be created
     * \param[in] nb : number of collections of the partition
     * \return the created Partition instance.
     */
    template<typename Type>
    static Partition<Type>* createPartition (ICell* parent, const std::string& name, size_t nb)
    {
        StorageHDF5* storage = dynamic_cast<StorageHDF5*> (ICell::getRoot (parent));
        assert (storage != 0);

        std::string actualName = parent->getFullId('/') + "/" + name;

        /** We create the HDF5 group if needed. This will just create the entry if not already existing. */
        GroupHDF5 group (storage, parent, name);

        /** If the nb of partitions is null, we try to get it from a property. */
        if (nb==0)
        {
            nb = StorageHDF5Factory::getPartitionsNb (&group);

            /** For backward compatibility only (ugly) : we look for the attribute in the parent object. */
            if (nb==0)  {  nb = StorageHDF5Factory::getPartitionsNb (dynamic_cast<Group*> (parent));  }

            /** Ok, we can't find any value. */
            if (nb==0)  {  throw system::Exception ("Partition '%s' has 0 items", name.c_str());      }
        }
        else
        {
            std::stringstream ss; ss << nb;
            group.addProperty (StorageHDF5Factory::getNbPartitionsName(), ss.str());
        }

        return new Partition<Type> (storage->getFactory(), parent, name, nb);
    }

    /** Create a Collection instance and attach it to a cell in a storage.
     * \param[in] parent : parent of the collection to be created
     * \param[in] name : name of the collection to be created
     * \param[in] synchro : synchronizer instance if needed
     * \return the created Collection instance.
     */
    template<typename Type>
    static CollectionNode<Type>* createCollection (ICell* parent, const std::string& name, system::ISynchronizer* synchro)
    {
#if 1
        synchro = GlobalSynchro::singleton();
#endif

        StorageHDF5* storage = dynamic_cast<StorageHDF5*> (ICell::getRoot (parent));
        assert (storage != 0);

        std::string actualName = parent->getFullId('/') + "/" + name;

        /** NOTE: we use here CollectionHDF5Patch and not CollectionHDF5 in order to reduce resources leaks due to HDF5.
         * (see also comments in CollectionHDF5Patch). */
        return new CollectionNode<Type> (storage->getFactory(), parent, name, new CollectionHDF5Patch<Type>(storage->getFileId(), actualName, synchro, parent->getCompressLevel()));
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

    /* */
    static const char* getNbPartitionsName()  { return "nb_partitions"; }

    /** Get the number of partitions (if any).
     * \param[in] group : check the attribute in this group
     * \return 0 if no attribute, otherwise the number of partitions */
    static size_t getPartitionsNb (Group* group)
    {
        size_t result = 0;
        if (group != 0)
        {
            std::string nbPartStr = group->getProperty (getNbPartitionsName());
            if (nbPartStr.empty()==false)  { result = atoi (nbPartStr.c_str()); }
        }
        return result;
    }

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
                _fileId = H5Fopen (getActualName().c_str(), H5F_ACC_RDWR, H5P_DEFAULT); /* FIXME: now all files are opened as read/write because I need that in Graph.cpp for opening exiting h5 files this needs better interface */
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

        hid_t getFileId ()  { return _fileId; }

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
                htri_t doesExist = H5Lexists (storage->getFileId(), actualName.c_str(), H5P_DEFAULT);

                if (doesExist <= 0)
                {
                    _groupId = H5Gcreate (storage->getFileId(), actualName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                }
                else
                {
                    _groupId = H5Gopen2 (storage->getFileId(), actualName.c_str(), H5P_DEFAULT);
                }
            }
            else
            {
                _groupId = H5Gopen2 (storage->getFileId(), "/", H5P_DEFAULT);
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

            /** We first check that the attribute exitst. */
            if ( H5Aexists(_groupId, key.c_str()) > 0)
            {
                hid_t datatype = H5Tcopy (H5T_C_S1);  H5Tset_size (datatype, H5T_VARIABLE);

                hid_t attrId = H5Aopen (_groupId, key.c_str(), H5P_DEFAULT);
                if (attrId >= 0)
                {
                    herr_t status;
                    hid_t space_id = H5Aget_space (attrId);

                    hsize_t dims = 1;
                    H5Sget_simple_extent_dims (space_id, &dims, NULL);
                    char** rdata = (char **) MALLOC (dims * sizeof (char *));

                    status = H5Aread (attrId, datatype, rdata);
                    if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Aread), status %d key", status); }

                    /** We set the result. */
                    result.assign (rdata[0]);

                    /** We release buffers. */
                    status = H5Dvlen_reclaim (datatype, space_id, H5P_DEFAULT, rdata);
                    if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Dvlen_reclaim), status %d", status); }
                    FREE (rdata);

                    /** We close resources. */
                    H5Aclose (attrId);
                    H5Tclose (datatype);
                    H5Sclose (space_id);
                }
            }

            return result;
        }
        
        void delProperty (const std::string& key)
        {
            H5Adelete(_groupId, key.c_str());
        }

        /** hack to set the attribute if it already exists: so i'm deleting and inserting again.
         * I had cleaner code based on H5Aopen / H5Awrite but I didn't know if it was the cause of 
         * a bug or not, so I opted for this failsafe solution */
        void setProperty (const std::string& key, const std::string value)
        {
            if ( H5Aexists(_groupId, key.c_str()) > 0)
            {
                delProperty(key);
            }
            addProperty(key, value);
        }

    private:
        hid_t _groupId;
    };
};

/********************************************************************************/

template<typename T1, typename T2=T1>
struct HDF5Pair : std::pair<T1,T2>
{
    HDF5Pair ()  {}

    HDF5Pair (const T1& a, const T2& b) : std::pair<T1,T2>(a,b) {}

    static hid_t hdf5 (bool& isCompound)
    {
        hid_t result = H5Tcreate (H5T_COMPOUND, sizeof(HDF5Pair));
        H5Tinsert (result, "first",   HOFFSET(HDF5Pair, first),  T1::hdf5(isCompound));
        H5Tinsert (result, "second",  HOFFSET(HDF5Pair, second), T2::hdf5(isCompound));

        isCompound = true;

        return result;
    }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_STORAGE_IMPL_STORAGE_HDF5_HPP_ */
