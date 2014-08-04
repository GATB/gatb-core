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
#include <gatb/tools/math/LargeInt.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/system/api/types.hpp>
#include <gatb/tools/misc/api/Enums.hpp>

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
    u_int64_t operator ()  (const Item& key, size_t idx)  {  return hash1 (key, seed_tab[idx]);  }

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

/** \brief Bloom interface
 *
 * We define a Bloom filter has a container (something that tells whether an item is here or not) and
 * a bag (something we can put items into).
 *
 * As expected, there is no Iterable interface here.
 */
template <typename Item> class IBloom : public Container<Item>, public Bag<Item>, public system::SmartPointer
{
public:
    virtual ~IBloom() {}

    virtual u_int8_t*& getArray    () = 0;
    virtual u_int64_t  getSize     () = 0;
    virtual u_int64_t  getBitSize  () = 0;
    virtual size_t     getNbHash   () const = 0;

    virtual std::string  getName   () const  = 0;
};

/********************************************************************************/

/** \brief Bloom filter implementation
 */
template <typename Item> class BloomContainer : public IBloom<Item>
{
public:

    /** Constructor.
     * \param[in] tai_bloom : size (in bits) of the bloom filter.
     * \param[in] nbHash : number of hash functions to use */
    BloomContainer (u_int64_t tai_bloom, size_t nbHash = 4)
        : _hash(nbHash), n_hash_func(nbHash), blooma(0), tai(tai_bloom), nchar(0), isSizePowOf2(false)
    {
        nchar  = (1+tai/8LL);
        blooma = (unsigned char *) system::impl::System::memory().malloc (nchar*sizeof(unsigned char)); // 1 bit per elem
        system::impl::System::memory().memset (blooma, 0, nchar*sizeof(unsigned char));

        /** We look whether the provided size is a power of 2 or not.
         *   => if we have a power of two, we can optimize the modulo operations. */
        isSizePowOf2 = (tai && !(tai & (tai - 1)));

        /** In case we have a power of 2^N, we set the size to 2^N-1 and use the classical trick:
         *     a % 2^N  <=>  a & (2^N-1)
         * */
        if (isSizePowOf2)  {  tai --;  }
    }

    /** Destructor. */
    virtual ~BloomContainer ()
    {
        system::impl::System::memory().free (blooma);
    }

    /** */
    size_t getNbHash () const { return n_hash_func; }

    /** \copydoc Container::contains. */
    bool contains (const Item& item)
    {
        if (isSizePowOf2)
        {
            for (size_t i=0; i<n_hash_func; i++)
            {
                u_int64_t h1 = _hash (item,i) & tai;
                if ((blooma[h1 >> 3 ] & bit_mask[h1 & 7]) != bit_mask[h1 & 7])  {  return false;  }
            }
        }
        else
        {
            for (size_t i=0; i<n_hash_func; i++)
            {
                u_int64_t h1 = _hash (item,i) % tai;
                if ((blooma[h1 >> 3 ] & bit_mask[h1 & 7]) != bit_mask[h1 & 7])  {  return false;  }
            }
        }
        return true;
    }

    /** */
    virtual u_int8_t*& getArray     ()  { return blooma; }
    virtual u_int64_t  getSize      ()  { return nchar;  }
    virtual u_int64_t  getBitSize   ()  { return tai;    }

    virtual std::string  getName    () const  = 0;

protected:

    HashFunctors<Item> _hash;
    int n_hash_func;

    u_int8_t* blooma;
    u_int64_t tai;
    u_int64_t nchar;
    bool      isSizePowOf2;
};

/********************************************************************************/
/** \brief Bloom filter implementation
 */
template <typename Item> class Bloom : public BloomContainer<Item>
{
public:

    /** \copydoc BloomContainer::BloomContainer */
    Bloom (u_int64_t tai_bloom, size_t nbHash = 4)  : BloomContainer<Item> (tai_bloom, nbHash)  {}

    /** \copydoc Bag::insert. */
    void insert (const Item& item)
    {
        if (this->isSizePowOf2)
        {
            for (size_t i=0; i<this->n_hash_func; i++)
            {
                u_int64_t h1 = this->_hash (item,i) & this->tai;
                this->blooma [h1 >> 3] |= bit_mask[h1 & 7];
            }
        }
        else
        {
            for (size_t i=0; i<this->n_hash_func; i++)
            {
                u_int64_t h1 = this->_hash (item,i) % this->tai;
                this->blooma [h1 >> 3] |= bit_mask[h1 & 7];
            }
        }
    }

    /** \copydoc Bag::flush */
    void flush ()  {}

    /** */
    std::string  getName () const { return "Bloom"; }

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
template <typename Item> class BloomNull : public IBloom<Item>
{
public:

    virtual ~BloomNull() {}

    u_int8_t*& getArray    () { return a; }
    u_int64_t  getSize     () { return 0; }
    u_int64_t  getBitSize  () { return 0; }
    size_t     getNbHash   () const { return 0; }

    virtual std::string  getName   () const  { return "BloomNull"; }

    bool contains (const Item& item) { return false; }
    void insert (const Item& item) {}
    void flush ()  {}
private:
    u_int8_t* a;
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
        if (this->isSizePowOf2)
        {
            for (size_t i=0; i<this->n_hash_func; i++)
            {
                u_int64_t h1 = this->_hash (item,i) & this->tai;
                __sync_fetch_and_or (this->blooma + (h1 >> 3), bit_mask[h1 & 7]);
            }
        }
        else
        {
            for (size_t i=0; i<this->n_hash_func; i++)
            {
                u_int64_t h1 = this->_hash (item,i) % this->tai;
                __sync_fetch_and_or (this->blooma + (h1 >> 3), bit_mask[h1 & 7]);
            }
        }
    }

    /** */
    std::string  getName () const { return "basic"; }
};
    

/********************************************************************************/
/** \brief Bloom filter implementation
 */
template <typename Item> class BloomCacheCoherent : public Bloom<Item>
{
public:
    
    /** Constructor.
     * \param[in] tai_bloom : size (in bits) of the bloom filter.
     * \param[in] nbHash : number of hash functions to use
     * \param[in] block_nbits : size of the block (actual 2^nbits) */
    BloomCacheCoherent (u_int64_t tai_bloom, size_t nbHash = 4,size_t block_nbits = 12)  : Bloom<Item> (tai_bloom + (1<<block_nbits), nbHash),_nbits_BlockSize(block_nbits)
    {
        _mask_block = (1<<_nbits_BlockSize) - 1;
        _reduced_tai = this->tai -  (1<<_nbits_BlockSize) ;
    }
    
    
    //for insert, no prefetch, perf is not important
     /** \copydoc Bag::insert. */
    void insert (const Item& item)
    {

        u_int64_t h0;

        h0 = this->_hash (item,0) % _reduced_tai;
        
        __sync_fetch_and_or (this->blooma + (h0 >> 3), bit_mask[h0 & 7]);

            for (size_t i=1; i<this->n_hash_func; i++)
            {
                u_int64_t h1 = h0  + (simplehash16( item, i) & _mask_block )   ;
                __sync_fetch_and_or (this->blooma + (h1 >> 3), bit_mask[h1 & 7]);
            }
        
    }
    
    /** */
    std::string  getName () const { return "cache"; }

    /** */
    u_int64_t  getBitSize   ()  { return _reduced_tai;    }
        
    /** \copydoc Container::contains. */
    bool contains (const Item& item)
    {

        u_int64_t tab_keys [20];
        u_int64_t h0;

        h0 = this->_hash (item,0) % _reduced_tai;
        __builtin_prefetch(&(this->blooma [h0 >> 3] ), 0, 3); //preparing for read

        //compute all hashes during prefetch
        for (size_t i=1; i<this->n_hash_func; i++)
        {
           tab_keys[i] =  h0  + (simplehash16( item, i) & _mask_block );// with simplest hash
            
        }
        
        
        if ((this->blooma[h0 >> 3 ] & bit_mask[h0 & 7]) != bit_mask[h0 & 7])  {  return false;  }
        
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 = tab_keys[i];
            if ((this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) != bit_mask[h1 & 7])  {  return false;  }
        }
        return true;
        
    }
    
    
private:
    u_int64_t _mask_block;
    size_t _nbits_BlockSize;
    u_int64_t _reduced_tai;
};

/********************************************************************************/

template <typename Item, size_t prec=1> class BloomGroupOld : public system::SmartPointer
{
public:

    typedef tools::math::LargeInt<prec> Result;

    /** */
    BloomGroupOld (u_int64_t size, size_t nbHash=4)
        : _hash(nbHash), _nbHash(nbHash), _size(size), _blooma(0)
    {
        _blooma = (Result*) system::impl::System::memory().malloc (_size*sizeof(Result));
        system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));
    }

    /** */
    BloomGroupOld (const std::string& uri)
        : _hash(0), _nbHash(0), _size(0), _blooma(0)
    {
        load (uri);
    }

    /** */
    ~BloomGroupOld ()  {  system::impl::System::memory().free (_blooma); }

    /** */
    std::string getName () const { return "BloomGroupOld"; }

    /** */
    void insert (const Item& item, size_t idx)
    {
        static const Result ONE (1);

        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;
#if 1
            this->_blooma[h1] |= (ONE << idx);
#else
            this->_blooma[h1].sync_fetch_and_or (ONE << idx);
#endif
        }
    }

    /** Return the size (in bytes). */
    u_int64_t getMemSize () const { return _size*sizeof(Result); }

    /** */
    void save (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "wb+");
        if (file != 0)
        {
            /** We write the nb of hash functions. */
            file->fwrite (&_nbHash, sizeof(_nbHash), 1);

            /** We write the size of the blooms. */
            file->fwrite (&_size, sizeof(_size), 1);

            /** We write the blooms info. */
            file->fwrite (_blooma, _size*sizeof(Result), 1);

#if 1
for (size_t i=0; i<10; i++)
{
    for (size_t j=0; j<prec; j++)
    {
        printf ("%8x ", _blooma[i].value[j]);
    }
    printf ("\n");
}
#endif

            delete file;
        }
    }

    /** */
    void load (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "rb+");
        if (file != 0)
        {
            /** We read the nb of hash functions. */
            file->fread (&_nbHash, sizeof(_nbHash), 1);

            /** We read the size of the blooms. */
            file->fread (&_size, sizeof(_size), 1);

            /** We allocate the array. */
            _blooma = (Result*) system::impl::System::memory().malloc (_size*sizeof(Result));
            system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));

            /** We read the blooms info. */
            file->fread (_blooma, _size*sizeof(Result), 1);

            delete file;
        }
    }

    /** */
    bool contains (const Item& item, size_t idx)
    {
        static const Result ONE (1);
        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;
            if ( (_blooma[h1] & (ONE << idx)) != (ONE << idx) )  {  return false;  }
        }
        return true;
    }

    /** */
    Result contains (const Item& item)
    {
        static const Result ZERO (0);
        Result res = ~ZERO;

        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;
            res &= _blooma [h1];
        }
        return res;
    }

private:

    HashFunctors<Item> _hash;
    size_t             _nbHash;
    u_int64_t          _size;
    Result*            _blooma;
};

/********************************************************************************/

template <typename Item, size_t prec=1> class BloomGroup : public system::SmartPointer
{
public:

    class Result
    {
    public:
        Result (u_int64_t v=0)  { memset (v); }

              u_int64_t& operator[] (size_t idx)       { return value[idx]; }
        const u_int64_t& operator[] (size_t idx) const { return value[idx]; }

        const u_int64_t* array () const { return value; }

        Result& operator&= (const Result& r)
        {
            for (size_t j=0; j<prec; j++)  { (*this)[j] &=  r[j]; }
            return *this;
        }

    private:
        u_int64_t value[prec];
        void memset (u_int64_t v)  {  system::impl::System::memory().memset (value, v, prec*sizeof(u_int64_t));  }

        friend class BloomGroup<Item,prec>;
    };

    /** */
    BloomGroup (u_int64_t size, u_int64_t maxMemory, size_t nbHash=4)
        : _hash(nbHash), _nbHash(nbHash), _size(size), _blooma(0)
    {
        printf ("BloomGroup:  size=%ld   sizeof(Result)=%d  maxMemory=%ld\n", size, sizeof(Result), maxMemory);
        if (_size*sizeof(Result) > maxMemory)
        {
            _size = maxMemory /sizeof (Result);
        }
        else
        {
            maxMemory = _size*sizeof(Result);
        }
        printf ("===> size=%ld   allocMemory=%ld\n", _size, sizeof(Result)*_size);

        _blooma = (Result*) system::impl::System::memory().malloc (_size*sizeof(Result));
        system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));
    }

    /** */
    BloomGroup (const std::string& uri)
        : _hash(0), _nbHash(0), _size(0), _blooma(0)
    {
        load (uri);
    }

    /** */
    ~BloomGroup ()  {  system::impl::System::memory().free (_blooma); }

    /** */
    std::string getName () const { return "BloomGroup"; }

    /** Return the size (in bytes). */
    u_int64_t getMemSize () const { return _size*sizeof(Result); }

    /** */
    void save (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "wb+");
        if (file != 0)
        {
            /** We write the nb of hash functions. */
            file->fwrite (&_nbHash, sizeof(_nbHash), 1);

            /** We write the size of the blooms. */
            file->fwrite (&_size, sizeof(_size), 1);

            /** We write the blooms info. */
            file->fwrite (_blooma, _size*sizeof(Result), 1);

            delete file;
        }
    }

    /** */
    void load (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "rb+");
        if (file != 0)
        {
            /** We read the nb of hash functions. */
            file->fread (&_nbHash, sizeof(_nbHash), 1);

            /** We read the size of the blooms. */
            file->fread (&_size, sizeof(_size), 1);

            /** We allocate the array. */
            _blooma = (Result*) system::impl::System::memory().malloc (_size*sizeof(Result));
            system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));

            /** We read the blooms info. */
            file->fread (_blooma, _size*sizeof(Result), 1);

            delete file;
        }
    }

    /** Insert an item in the 'idx'th Bloom filter
     * \param[in] item : item to be inserted into the filter
     * \param[in] idx : index of the filter */
    void insert (const Item& item, size_t idx)
    {
        u_int64_t q,mask;  euclidian(idx,q,mask);

        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;

#if 1
            this->_blooma[h1][q] |= mask;
#else
            __sync_fetch_and_or (this->_blooma[h1].value + q, mask);
#endif
        }
    }

    /** */
    bool contains (const Item& item, size_t idx)
    {
        u_int64_t q,mask;  euclidian(idx,q,mask);

        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;
            if ( (_blooma[h1][q] & mask) != mask )  {  return false;  }
        }
        return true;
    }

    /** */
    Result contains (const Item& item)
    {
        Result res (~0);

        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;
            res &=  _blooma [h1];
        }
        return res;
    }

private:

    HashFunctors<Item> _hash;
    size_t             _nbHash;
    u_int64_t          _size;
    Result*            _blooma;

    void euclidian (size_t idx, u_int64_t& dividend, u_int64_t& mask) const
    {
        dividend = idx / (8*sizeof(u_int64_t));
        mask     = ((u_int64_t) 1) << (idx % (8*sizeof(u_int64_t)));
    }
};


/********************************************************************************/

template <typename Item, size_t prec=1> class BloomGroupCacheCoherent : public system::SmartPointer
{
public:

    typedef tools::math::LargeInt<prec> Result;

    /** */
    BloomGroupCacheCoherent (u_int64_t size, size_t nbHash=4, size_t block_nbits=12)
        : _hash(nbHash), _nbHash(nbHash), _size(size), _blooma(0), _nbits_BlockSize(block_nbits)
    {
        _size += (1<<_nbits_BlockSize);

        _blooma = (Result*) system::impl::System::memory().malloc (_size*sizeof(Result));
        system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));

        _mask_block   = (1<<_nbits_BlockSize) - 1;
        _reduced_size = this->_size -  (1<<_nbits_BlockSize) ;
    }

    /** */
    BloomGroupCacheCoherent (const std::string& uri)
        : _hash(0), _nbHash(0), _size(0), _blooma(0)
    {
        load (uri);

        _mask_block   = (1<<_nbits_BlockSize) - 1;
        _reduced_size = this->_size -  (1<<_nbits_BlockSize) ;
    }

    /** */
    ~BloomGroupCacheCoherent ()  {  system::impl::System::memory().free (_blooma); }

    /** */
    std::string getName () const { return "BloomGroupCacheCoherent"; }

    /** Return the size (in bytes). */
    u_int64_t getMemSize () const { return _size*sizeof(Result); }

    /** */
    void save (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "wb+");
        if (file != 0)
        {
            /** We write the nb of hash functions. */
            file->fwrite (&_nbHash, sizeof(_nbHash), 1);

            /** We write the size of the blooms. */
            file->fwrite (&_size, sizeof(_size), 1);

            /** We write block size information. */
            file->fwrite (&_nbits_BlockSize, sizeof(_nbits_BlockSize), 1);

            /** We write the blooms info. */
            file->fwrite (_blooma, _size*sizeof(Result), 1);

            delete file;
        }
    }

    /** */
    void load (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "rb+");
        if (file != 0)
        {
            /** We read the nb of hash functions. */
            file->fread (&_nbHash, sizeof(_nbHash), 1);

            /** We read the size of the blooms. */
            file->fread (&_size, sizeof(_size), 1);

            /** We read block size information. */
            file->fread (&_nbits_BlockSize, sizeof(_nbits_BlockSize), 1);

            /** We allocate the array. */
            _blooma = (Result*) system::impl::System::memory().malloc (_size*sizeof(Result));
            system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));

            /** We read the blooms info. */
            file->fread (_blooma, _size*sizeof(Result), 1);

            delete file;
        }
    }

    /** */
    void insert (const Item& item, size_t idx)
    {
        static const Result ONE (1);

        /** First hash. */
        u_int64_t h0 = this->_hash (item,0) % _reduced_size;
        this->_blooma[h0] |= (ONE << idx);

        for (size_t i=1; i<this->_nbHash; i++)
        {
            /** Other hash. */
            u_int64_t h1 = h0  + (simplehash16 (item,i) & _mask_block);
            this->_blooma[h1] |= (ONE << idx);
        }
    }

    /** */
    bool contains (const Item& item, size_t idx)
    {
        static const Result ONE  (1);

        /** First hash. */
        u_int64_t h0 = this->_hash (item,0) % _reduced_size;
        if ((this->_blooma[h0] & (ONE << idx)) != (ONE << idx))  {  return false;  }

        for (size_t i=1; i<this->_nbHash; i++)
        {
            /** Other hash. */
            u_int64_t h1 = h0  + (simplehash16 (item,i) & _mask_block);
            if ((this->_blooma[h1] & (ONE << idx)) != (ONE << idx))  {  return false;  }
        }
        return true;
    }

    /** */
    Result contains (const Item& item)
    {
        static const Result ZERO (0);
        Result res = ~ZERO;

        u_int64_t h0 = this->_hash (item,0) % _reduced_size;
        res &= _blooma [h0];

        for (size_t i=1; i<this->_nbHash; i++)
        {
            u_int64_t h1 = h0 + (simplehash16(item,i) & _mask_block);
            res &= _blooma [h1];
        }
        return res;
    }

private:

    HashFunctors<Item> _hash;
    size_t             _nbHash;
    u_int64_t          _size;
    Result*            _blooma;

    u_int64_t  _mask_block;
    size_t     _nbits_BlockSize;
    u_int64_t  _reduced_size;

};

/********************************************************************************/

/** */
class BloomFactory
{
public:

    /** */
    static BloomFactory& singleton()  { static BloomFactory instance; return instance; }

    /** */
    template<typename T> IBloom<T>* createBloom (tools::misc::BloomKind kind, u_int64_t tai_bloom, size_t nbHash)
    {
        switch (kind)
        {
            case tools::misc::BLOOM_NONE:      return new BloomNull<T>          ();
            case tools::misc::BLOOM_BASIC:     return new BloomSynchronized<T>  (tai_bloom, nbHash);
            case tools::misc::BLOOM_CACHE:     return new BloomCacheCoherent<T> (tai_bloom, nbHash);
            case tools::misc::BLOOM_DEFAULT:   return new BloomCacheCoherent<T> (tai_bloom, nbHash);
            default:        throw system::Exception ("bad Bloom kind %d in createBloom", kind);
        }
    }

    /** */
    template<typename T> IBloom<T>* createBloom (std::string name, std::string sizeStr, std::string nbHashStr)
    {
        tools::misc::BloomKind kind;  parse (name, kind);
        return createBloom<T> (kind, (u_int64_t)atol (sizeStr.c_str()), (size_t)atol (nbHashStr.c_str()));
    }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BLOOM_HPP_ */
