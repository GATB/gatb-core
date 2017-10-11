#ifndef OGRAPH
#define OGRAPH

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <cstdlib>
#include <stdint.h>

#include <gatb/tools/math/NativeInt64.hpp>
#include <gatb/tools/math/NativeInt128.hpp>
#include <gatb/tools/math/LargeInt.hpp>
#include <gatb/kmer/impl/Model.hpp>
using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;



namespace gatb { namespace core { namespace debruijn { namespace impl  {

enum seq_pos { SEQ_LEFT, SEQ_RIGHT } ; // indicates the location in the input sequence (regardless of orientation) 

template<size_t span>
struct kmerIndiceT{
    typedef typename Kmer<span>::Type  Type;
    //typedef __uint128_t  Type;
	uint32_t indice;
	Type kmmer;
    seq_pos position; 
};



template<size_t span>
struct comparator{bool operator()(const kmerIndiceT<span>& a , const kmerIndiceT<span>& b) { 
    if (a.kmmer != b.kmmer)
        return a.kmmer < b.kmmer; 
    // it used to be just a kmmer comparison, i'm making it a bit less undefined on identical kmers now
    return a.indice < b.indice; }};


template<size_t span>
class graph3{
	public:
        typedef typename Kmer<span>::Type  kmerType;
        //typedef __uint128_t  kmerType;
        typedef kmerIndiceT<span>  kmerIndice;
		uint k,indiceUnitigs,nbElement,minimizer,minsize,nb_pretips;
		std::string*       unitigs;
        std::vector<uint>* unitigs_abundances;
        std::vector<kmerIndice> left;
        std::vector<kmerIndice> right;
        std::vector<bool> connected_left, connected_right, indexed_left, indexed_right;
		void addvertex(std::string& str);
        void addtuple(std::tuple<std::string,uint,uint,uint>& tuple);
		void addleftmin(unsigned int mini);
		void addrightmin(unsigned int mini);
		void update_connected(kmerIndiceT<span> &ki);
		void debruijn();
        void debruijn2();
        void compaction2(uint iL, uint iR);
		void compress();
        kmerType end2int128rc(const std::string& str);
        kmerType end2int128(const std::string& str);
        kmerType beg2int128rc(const std::string& str);
        kmerType beg2int128(const std::string& str);
        kmerType rcb(kmerType min);
		void compaction(uint iR, uint iL, kmerType kmmer);
		void compact_abundances(uint i1, uint i2, bool reverse_first=false, bool reverse_second=false);
		uint size();
        bool output(uint i);
        bool clear();
        bool pre_tip_cleaning;

		graph3(uint ka, uint min,uint size, uint nb){
            indiceUnitigs=0;
            pre_tip_cleaning = false;
            nb_pretips=0;
			minsize=size;
			k=ka;
			minimizer=min;
            nbElement=nb;
            unitigs=new std::string [nbElement];
            left.reserve(nbElement);
    	    right.reserve(nbElement);
            connected_left.reserve(nbElement);
            connected_right.reserve(nbElement);
            indexed_left.reserve(nbElement);
            indexed_right.reserve(nbElement);
            unitigs_abundances=new std::vector<uint> [nbElement];
		}
};

}}}}

#endif
