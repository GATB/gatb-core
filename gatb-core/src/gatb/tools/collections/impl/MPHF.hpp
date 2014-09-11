#ifdef WITH_MPHF
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

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_MPHF_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_MPHF_HPP_

/********************************************************************************/

#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/api/Abundance.hpp>

#include <emphf/common.hpp>
#include <emphf/mphf.hpp>
#include <emphf/base_hash.hpp>
#include <emphf/compute_mphf_generic.hpp>
#include <emphf/mmap_memory_model.hpp>
#include <emphf/hypergraph_sorter_scan.hpp>


#include <string>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

template <size_t span, typename ValueType> class MPHF : public system::SmartPointer

{
    // OK so.. this isn't a generic MPHF, it's typed so that keys are kmers. that could be typedef'd, but I'm too lazy
    typedef typename kmer::impl::Kmer<span>::Count Count;
    typedef typename kmer::impl::Kmer<span>::Type KmerType;
    typedef typename dp::Iterator<Count> IteratorCount;

    public:

    // adapted from compute_mphf_scan_mmap.cpp
    typedef emphf::hypergraph_sorter_scan<uint32_t, emphf::mmap_memory_model> HypergraphSorter32;
    typedef emphf::hypergraph_sorter_scan<uint64_t, emphf::mmap_memory_model> HypergraphSorter64;
    typedef emphf::jenkins64_hasher BaseHasher;
    typedef emphf::mphf<BaseHasher> mphf_t;


    // gave up on writing a kmer adaptor that would read the large int raw data and used stl_string_adaptor instead
#if 0
        struct kmer_adaptor
        {
            byte_range_t operator()(typename gatb::core::kmer::impl::Kmer<span>::Count kmer_count) const
            {
                // can be optimized, that toString is ugly and slow
                const uint8_t* buf = reinterpret_cast<uint8_t const*>( /* todo */);
                const uint8_t* end = buf + /* todo */;
                return byte_range_t(buf, end);

            }
        };
#endif

    // adapted from line_iterator and file_lines from emphf:common.cpp
    // because emphf wants a proper stl iterator, not Erwan's custom classes :)

    private:
    class kmer_iterator
        : public std::iterator<std::forward_iterator_tag, const typename kmer::impl::Kmer<span>::Count> {

            public:
                kmer_iterator()
                    : solidIterable(nullptr), pos(0), itKmers(nullptr)
                {             }

                kmer_iterator(Iterable<Count>* solidIterable, int sizeKmer)
                    : solidIterable(solidIterable), pos(0), sizeKmer(sizeKmer)
                {
                    itKmers =                             solidIterable->iterator()                            ;

                    itKmers->first();
                    s.assign(itKmers->item().value.toString(sizeKmer));
                }

                std::string const& operator*() const {
                    // for some reason I can't do that instruction here:
                    //s.assign(itKmers->item().value.toString(sizeKmer));
                    // but I can in ++() and the constructor.. so I did that.
                    // bah.. c++..
                    return s;
                }

                kmer_iterator& operator++() {
                    itKmers->next();
                    pos++;
                    if (itKmers->isDone())
                    {
                        itKmers = nullptr;
                        pos = 0;

                    }
                    else
                        s.assign(itKmers->item().value.toString(sizeKmer));
                    return *this;
                }

                friend bool operator==(kmer_iterator const& lhs, kmer_iterator const& rhs)
                {

                    if (!lhs.itKmers || !rhs.itKmers) {
                        if (!lhs.itKmers && !rhs.itKmers) {
                            return true;
                        } else {
                            return false;
                        }
                    }

                    assert(lhs.itKmers == rhs.itKmers);
                    return rhs.pos == lhs.pos;
                }

                friend bool operator!=(kmer_iterator const& lhs, kmer_iterator const& rhs)
                {
                    return !(lhs == rhs);
                }

            private:
                IteratorCount* itKmers;
                Iterable<Count>* solidIterable;
                unsigned long pos;
                std::string s;
                int sizeKmer;
        };


    class kmer_iterator_wrapper
    {
        public:
        kmer_iterator_wrapper( Iterable<Count>* solidIterable, int sizeKmer) : solidIterable(solidIterable), sizeKmer(sizeKmer)
        {}

        kmer_iterator begin() const
        {
            return kmer_iterator(solidIterable, sizeKmer);
        }

        kmer_iterator end() const { return kmer_iterator(); }

        size_t size() const
        {
            printf("kmer_iterator_wrapper size called but not impemented!\n");
            return 0;
        }

        private:
        // noncopyble // blindly copies from emphf, not sure what this does
        kmer_iterator_wrapper(kmer_iterator_wrapper const&);
        kmer_iterator_wrapper& operator=(kmer_iterator_wrapper const&);

        Iterable<Count>* solidIterable;
        int sizeKmer;
    };


    public:

    /** Constructor.
     */

    // light constructor without emphf creation    
    MPHF (int kmerSize, uint64_t nbElts) : _kmerSize(kmerSize), _nbElts(nbElts)
    {
        // a small fix, emphf for 2 nodes doesn't seem to work
        if (nbElts <= 2)
        {
            nbElts = 3;
            _nbElts = 3;
        }


    }

    // full constructor
    MPHF (int kmerSize, uint64_t nbElts, Iterable<Count>* solidIterable) : _kmerSize(kmerSize), _nbElts(nbElts)
    {
        kmer_iterator_wrapper kmers(solidIterable, _kmerSize);

        fprintf(stderr, "Preparing to construct a MPHF for %ld nodes\n", (long int)nbElts);

        // a small fix, emphf for 2 nodes doesn't seem to work
        if (nbElts <= 2)
        {
            nbElts = 3;
            _nbElts = 3;
        }

        size_t max_nodes = (size_t(std::ceil(double(nbElts) * 1.23)) + 2) / 3 * 3;
        if (max_nodes >= uint64_t(1) << 32) {
            emphf::logger() << "Using 64-bit sorter" << std::endl;
            HypergraphSorter64 sorter;
            mphf_t(sorter, nbElts, kmers, adaptor).swap(mphf);
        } else {
            emphf::logger() << "Using 32-bit sorter" << std::endl;
            HypergraphSorter32 sorter;
            mphf_t(sorter, nbElts, kmers, adaptor).swap(mphf);
        }

    }

    void populateAbundances(Iterable<Count>* solidIterable)
    {
        data = new ValueType[_nbElts];

        long n = _nbElts;

        {
            emphf::logger() << "Setting abundances of " << n << " kmers." << std::endl;
            long nb_iterated = 0;
            long nb_abundances_above_precision = 0;

            IteratorCount* itKmers =                   solidIterable->iterator();

            // set counts and at the same time, test the mphf
            for (itKmers->first(); !itKmers->isDone(); itKmers->next())
            {
                uint64_t h = get_index(itKmers->item().value);

                if (h >= n) {
                    emphf::logger() << "ERROR during MPHF check: value out of bounds "
                        << h << std::endl;
                    exit(1);
                }

                int abundance = itKmers->item().abundance;
                if (abundance > 255)
                {
                    nb_abundances_above_precision++;
                    abundance = 255;
                }

                //emphf::logger() << "index of kmer " << itKmers->item().value.toString(_kmerSize) << " is " << h << " and setting abundance " << abundance << std::endl;

                set_at_index(h, abundance);
                nb_iterated ++;
            }

            if (nb_iterated != n && n > 3)
            {
                emphf::logger() << "ERROR during abundance population: itKmers iterated over " << nb_iterated << "/" << n << " k-mers only" << std::endl;
                    exit(1);

            }

            if (nb_abundances_above_precision > 0)
                emphf::logger() << "WARNING: " << nb_abundances_above_precision << " k-mer counts were clipped to 255" << std::endl;

            nb_iterated = 0;
            // you know what? let's test if the MPHF does not have collisions, it won't hurt.
            for (itKmers->first(); !itKmers->isDone(); itKmers->next())
            {
                int abundance = itKmers->item().abundance;

                unsigned char given_abundance = 0;
                get(itKmers->item().value, &given_abundance);
                int int_abundance = given_abundance;
        
                uint64_t h = get_index(itKmers->item().value);

                //emphf::logger() << "index of kmer " << itKmers->item().value.toString(_kmerSize) << " is " << h << " and expected abundance is " << abundance << " while stored abundance is " << int_abundance << std::endl;

                if (int_abundance != abundance && abundance < 255)
                {
                    emphf::logger() << "ERROR: MPHF isn't injective (abundance population failed)" << std::endl;
                    exit(1);
                }
                nb_iterated ++;
            }

            if (nb_iterated != n && n > 3)
            {
                emphf::logger() << "ERROR during abundance population check: itKmers iterated over " << nb_iterated << "/" << n << " k-mers only" << std::endl;
                    exit(1);

            }

        }
    }

    ~MPHF() { }

    /** */
    void set(const KmerType& graine, ValueType value)
    {
        uint64_t index = get_index(graine);
        data[index] = value;
    }

    /** */
    bool get(const KmerType& graine, ValueType * val=0)
    {
        uint64_t index = get_index(graine);
        *val = data[index];
        return true;
    }

    uint64_t get_index(const KmerType& graine)
    {
        std::string s = graine.toString(_kmerSize);
        return mphf.lookup(s, adaptor);
    }

    void set_at_index(uint64_t index, ValueType value)
    {
        data[index] = value;
    }

    uint64_t getSize ()
    {
        return mphf.size();
    }

    void setEMPHF(void * mphf_void)
    {
        static_cast< emphf::mphf<BaseHasher>* >(mphf_void)->swap(mphf);
    }
    
    mphf_t mphf;

protected:

    //kmer_adaptor adaptor;
    emphf::stl_string_adaptor adaptor;

    int _kmerSize;
    uint64_t _nbElts;
    ValueType *data;

};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif
#endif
