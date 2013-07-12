/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file Hash16.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Container implementation
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_OAHASH_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_OAHASH_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Container.hpp>
#include <gatb/tools/collections/api/Bag.hpp>
#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/system/impl/System.hpp>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

template <typename Item> class OAHash
{
    typedef std::pair<Item,u_int32_t>  element_pair;

public:

    /** */
    static int size_entry ()  {  return sizeof(element_pair); }

    /** */
    int getMaxNbItems ()  { return hash_size; }

    /** Constructor.
     * \param[in] max_memory : max memory for the hash table.*/
    OAHash (u_int64_t max_memory)
    {
        hash_size = max_memory / sizeof(element_pair);
        if (hash_size == 0)  {  throw system::Exception ("empty OAHash allocated");  }
        data = (element_pair *) calloc( hash_size, sizeof(element_pair));  //create hashtable
    }

    /** Destructor. */
    ~OAHash()  {  free (data);  }

    /** */
    void insert (const Item& graine, int value)
    {
        element_pair *element = find_slot(graine);
        if (!is_occupied(element))
            element->first = graine;
        element->second = value;
    }

    /** */
    void increment (const Item& graine)
    {
        element_pair *element = find_slot(graine);

        if (!is_occupied(element))
            element->first = graine;
        element->second = element->second + 1;
    }

    /** */
    bool get (const Item& graine, int * val=0)
    {
        element_pair *element = find_slot(graine, false);
        if (element == 0)  { return 0; }

        if (!is_occupied(element))
            return false;
        if (element->first == graine)
            if (val != NULL)
                *val = element->second;
        return true;
    }

    /** */
    bool has_key (const Item& graine)  {     return get(graine,NULL) == 1;  }

    /** */
    u_int64_t memory_usage() {     return hash_size* sizeof(element_pair); /* in bits */ }

    /** */
    float load_factor()
    {
        u_int64_t ptr = 0;
        u_int64_t nbKeys = 0;
        while (ptr < hash_size)
        {
            while ((ptr < hash_size) &&  ((data+ptr)->second == 0)  )
                ptr++;

            if (ptr == hash_size)
                break;

            nbKeys++;
            ptr++;
        }
        return (float)nbKeys/(float)hash_size;
    }

    /** */
    dp::Iterator <std::pair<Item,u_int32_t> >* iterator ()  {  return new Iterator (*this);  }


    /************************************************************/
    class Iterator : public tools::dp::Iterator <std::pair<Item,u_int32_t> >
    {
    public:

        Iterator (OAHash<Item>& aRef) : ref(aRef), iterator(0), iteratorMax(0), done(true)  {}

        /** \copydoc tools::dp::Iterator::first */
        void first()
        {
            iterator    = ref.data - 1;
            iteratorMax = ref.data + ref.hash_size;
            done        = false;

            next ();
        }

        /** \copydoc tools::dp::Iterator::next */
        void next()
        {
            while (!done)
            {
                ++iterator;
                done = (iterator >= iteratorMax);
                if (!done && iterator->second != 0)
                {
                    *this->_item = *iterator;  break;
                }
            }
        }

        /** \copydoc tools::dp::Iterator::isDone */
        bool isDone ()   {  return done; }

        /** \copydoc tools::dp::Iterator::item */
        std::pair<Item,u_int32_t>& item ()     { return *this->_item; }

    private:
        OAHash<Item>&  ref;
        element_pair*  iterator;
        element_pair*  iteratorMax;
        bool           done;
    };

protected:

    u_int64_t      hash_size;
    element_pair*  data;

    /** */
    element_pair * find_slot (const Item& key, bool exceptionOnBadKey = true)
    {
        u_int64_t ptr = oahash (key) % hash_size;
        element_pair* element = data+ptr;
        u_int64_t retries = 0;

        // search until we either find the key, or find an empty slot.
        while ( ( is_occupied(element)) && ( element->first != key ) && (retries < hash_size))
        {
            ptr = (ptr + 1) % hash_size;
            element = data+ptr;
            retries++;
        }

        if (retries == hash_size)
        {
            if (exceptionOnBadKey)  {  throw system::Exception ("OAHash: max rehashes reached: %lld (notify a developer)", hash_size);  }
            return 0;
        }

        return element;
    }


    bool is_occupied (element_pair *element)   {  return (element->second != 0); }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_OAHASH_HPP_ */
