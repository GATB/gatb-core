/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file Bloom.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Bloom implementation
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BLOOM_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BLOOM_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Container.hpp>
#include <gatb/tools/collections/api/Bag.hpp>
#include <gatb/system/impl/System.hpp>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

extern u_int8_t bit_mask [];

/********************************************************************************/
/**
 */
template <typename Item> class HashFunctors
{
public:

    /** */
    HashFunctors (size_t nbFct, u_int64_t seed=0) : _nbFct(nbFct), user_seed(seed)
    {
        generate_hash_seed ();
    }

    /** */
    u_int64_t operator ()  (const Item& key, size_t idx)  {  return hash (key, seed_tab[idx]);  }

private:

    /** */
    void generate_hash_seed ()
    {
        static const u_int64_t rbase[NSEEDSBLOOM] =
        {
            0xAAAAAAAA55555555ULL,  0x33333333CCCCCCCCULL,  0x6666666699999999ULL,  0xB5B5B5B54B4B4B4BULL,
            0xAA55AA5555335533ULL,  0x33CC33CCCC66CC66ULL,  0x6699669999B599B5ULL,  0xB54BB54B4BAA4BAAULL,
            0xAA33AA3355CC55CCULL,  0x33663366CC99CC99ULL
        };

        for (size_t i=0; i<NSEEDSBLOOM; ++i)  {  seed_tab[i] = rbase[i];  }
        for (size_t i=0; i<NSEEDSBLOOM; ++i)  {  seed_tab[i] = seed_tab[i] * seed_tab[(i+3) % NSEEDSBLOOM] + user_seed ;  }
    }

    size_t _nbFct;

    static const size_t NSEEDSBLOOM = 10;
    u_int64_t seed_tab[NSEEDSBLOOM];
    u_int64_t user_seed;
};


/********************************************************************************/
/** \brief Bloom filter implementation
 */
template <typename Item> class BloomContainer : public Container<Item>
{
public:

    /** Constructor.
     * \param[in] tai_bloom : size (in bits) of the bloom filter.
     * \param[in] nbHash : number of hash functions to use */
    BloomContainer (u_int64_t tai_bloom, size_t nbHash = 4)
        : _hash(nbHash), n_hash_func(nbHash), blooma(0), tai(tai_bloom), nchar(0)
    {
        nchar  = (1+tai/8LL);
        blooma = (unsigned char *) system::impl::System::memory().malloc (nchar*sizeof(unsigned char)); // 1 bit per elem
        system::impl::System::memory().memset (blooma, 0, nchar*sizeof(unsigned char));
    }

    /** Destructor. */
    virtual ~BloomContainer ()
    {
        system::impl::System::memory().free (blooma);
    }

    /** \copydoc Container::contains. */
    bool contains (const Item& item)
    {
        for (size_t i=0; i<n_hash_func; i++)
        {
            u_int64_t h1 = _hash (item,i) % tai;

            if ((blooma[h1 >> 3 ] & bit_mask[h1 & 7]) != bit_mask[h1 & 7])
                return false;
        }
        return true;
    }

protected:

    HashFunctors<Item> _hash;
    int n_hash_func;

    u_int8_t* blooma;
    u_int64_t tai;
    u_int64_t nchar;
};

/********************************************************************************/
/** \brief Bloom filter implementation
 */
template <typename Item> class Bloom : public BloomContainer<Item>, public Bag<Item>
{
public:

    /** \copydoc BloomContainer::BloomContainer */
    Bloom (u_int64_t tai_bloom, size_t nbHash = 4)  : BloomContainer<Item> (tai_bloom, nbHash)  {}

    /** \copydoc Bag::insert. */
    void insert (const Item& item)
    {
        for (size_t i=0; i<this->n_hash_func; i++)
        {
            u_int64_t h1 = this->_hash (item,i) % this->tai;
            this->blooma [h1 >> 3] |= bit_mask[h1 & 7];
        }
    }

    /** \copydoc Bag::flush */
    void flush ()  {}

    /** */
    void dump (const char* filename)
    {
        FILE* file = fopen(filename,"wb");
        fwrite (this->blooma, sizeof(unsigned char), this->nchar, file);
        fclose (file);
    }
};

/********************************************************************************/
/** \brief Bloom filter implementation
 */
template <typename Item> class BloomSynchronized : public Bloom<Item>
{
public:

    /** \copydoc BloomContainer::BloomContainer */
    BloomSynchronized (u_int64_t tai_bloom, size_t nbHash = 4)  : Bloom<Item> (tai_bloom, nbHash)  {}

    /** \copydoc Bag::insert. */
    void insert (const Item& item)
    {
        for (size_t i=0; i<this->n_hash_func; i++)
        {
            u_int64_t h1 = this->_hash (item,i) % this->tai;
            __sync_fetch_and_or (this->blooma + (h1 >> 3), bit_mask[h1 & 7]);
        }
    }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BLOOM_HPP_ */
