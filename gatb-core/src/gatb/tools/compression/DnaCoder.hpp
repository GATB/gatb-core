/*****************************************************************************
 *   Leon: reference free compression for NGS reads
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: G.Benoit, G.Rizk, C.Lemaitre
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

#ifndef _DNACODER_HPP_
#define _DNACODER_HPP_


#include <gatb/gatb_core.hpp>
//#include "RangeCoder.hpp"
#include "Leon.hpp"
//#include "CompressionUtils.hpp"

//#define PRINT_DISTRIB

#ifdef PRINT_DISTRIB
#include <unordered_set>
#include <unordered_map>
#endif


using namespace std;
class Leon;

//====================================================================================
// ** AbstractDnaCoder
//====================================================================================
class AbstractDnaCoder
{
	public:
		AbstractDnaCoder(Leon* leon);
		
	protected:
		KmerModel _kmerModel;

		Order0Model _readTypeModel; //only 2 value in this model: with anchor or without anchor

		Order0Model _noAnchorReadModel;

		Order0Model _bifurcationModel;
		Order0Model _bifurcationBinaryModel;
		Order0Model _readAnchorRevcompModel;
	
		Leon* _leon;
		collections::impl::IBloom<kmer_type>* _bloom; // the bloom containing the solid kmers
		
		Order0Model _readSizeDeltaTypeModel;
		Order0Model _anchorPosDeltaTypeModel;
		Order0Model _anchorAddressDeltaTypeModel;
		Order0Model _NposDeltaTypeModel;
		Order0Model _errorPosDeltaTypeModel;
		

		//Sequence* _prevSequences;
		//Order0Model _isPrevReadAnchorableModel;
		//vector<Order0Model> _isPrevReadAnchorablePosModel;
		
		vector<Order0Model> _anchorAddressModel;
		
		vector<Order0Model> _anchorPosModel;
		
		vector<Order0Model> _numericModel;
		vector<Order0Model> _leftErrorModel;
		vector<Order0Model> _rightErrorModel;
		
		vector<Order0Model> _NposModel;
		vector<Order0Model> _leftErrorPosModel;
		vector<Order0Model> _rightErrorPosModel;
		
		vector<Order0Model> _readSizeValueModel;


		vector<Order0Model> _noAnchorReadSizeValueModel;
		
		size_t _kmerSize;
		unsigned int _readSize;
	
		vector<int> _leftErrorPos;
		vector<int> _rightErrorPos;
		vector<int> _Npos;
		
		void startBlock();
		void endRead();
		void codeSeedBin(KmerModel* model, kmer_type* kmer, int nt, bool right);
		void codeSeedNT(KmerModel* model, kmer_type* kmer, char nt, bool right);

		void addErrorPos(int pos, bool rightExtend);
		
		u_int64_t _seqId;
		

	
		u_int64_t _prevReadSize;
		u_int64_t _prevAnchorPos;
		u_int64_t _prevAnchorAddress;
		u_int64_t _prevNpos;
		u_int64_t _prevErrorPos;
		u_int64_t _prevNbLeftError;
	
		int _processedSequenceCount;
};

//====================================================================================
// ** DnaEncoder
//====================================================================================
class DnaEncoder : AbstractDnaCoder
{
		
	public:
		
		DnaEncoder(Leon* leon);
		DnaEncoder(const DnaEncoder& copy);
		~DnaEncoder();

		void operator()(Sequence& sequence);
		
	private:
	
	
	//pour quals
	char * _qualseq;
	int * _nb_solids;
	int _smoothing_threshold;
	unsigned int _max_read_size;
	bool _trunc_mode;
	
	void storeSolidCoverageInfo();
	void smoothQuals();
	bool apply_smoothing_at_pos(int pos);

	double char2proba(char c);
	char char2phred(char c);

	char * _bufferQuals;
	unsigned int _bufferQuals_idx;
	unsigned int _bufferQuals_size;
	
#ifdef PRINT_DISTRIB
	vector<Sequence*> _sequences;
	const int maxSequences = 100;
	vector<u_int32_t> _distrib;
	u_int64_t _outDistrib;
#endif
	
	
	
	
		RangeEncoder _rangeEncoder;
		
		#ifdef LEON_PRINT_STAT
			RangeEncoder _rangeEncoder1;
			RangeEncoder _rangeEncoder2;
			RangeEncoder _rangeEncoder3;
			RangeEncoder _rangeEncoder4;
			RangeEncoder _rangeEncoder5;
			RangeEncoder _rangeEncoder6;
		#endif
		
		//static void encodeFirstHeader();
		void writeBlock();
		void execute();
		
		void buildKmers();
		bool isReadAnchorable();
		int findExistingAnchor(u_int32_t* anchorAddress);
		
		void encodeAnchorRead(int anchorPos, u_int32_t anchorAddress);
		kmer_type buildBifurcationList(int pos, kmer_type kmer, bool rightExtend);
		//int buildBifurcationList(int pos, bool rightExtend);
		int voteMutations(int pos, int depth, bool rightExtend);
		
		void encodeNoAnchorRead();

		int getBestPath(int pos, kmer_type& kmer, bitset<4>& initRes4, bool rightExtend);
		
		Sequence* _sequence;
		char* _readseq;
		vector<kmer_type> _kmers;
		KmerModel::Iterator _itKmer;
		vector<u_int8_t> _bifurcations;
		vector<u_int8_t> _binaryBifurcations;
		vector<u_int8_t> _bifurcationTypes;
		
		//bool _isPrevReadAnchorable;
		//u_int64_t _isPrevReadAnchorablePos;

		vector<int> _solidMutaChain;
		//int _solidMutaChainPos;
		u_int64_t _totalDnaSize;
		u_int64_t _readCount;
		u_int64_t _MCtotal;
		u_int64_t _readWithoutAnchorCount;
		u_int64_t _MCuniqSolid;
		u_int64_t _MCuniqNoSolid;
		u_int64_t _MCnoAternative;
		u_int64_t _MCmultipleSolid;
	
		int _minSequenceSize;
		int _maxSequenceSize;
		//u_int64_t _MCmultipleNoSolid;
	
		int _thread_id;

	//	int _solidMutaChainStartPos;
	//	int _solidMutaChainSize;
	//	int _solidMutaChainLockTime;
};

//====================================================================================
// ** DnaDecoder
//====================================================================================
class DnaDecoder : AbstractDnaCoder
{
		
	public:
		
		DnaDecoder(Leon* leon, const string& inputFilename,tools::storage::impl::Group *  group);
		~DnaDecoder();
		
		void setup(u_int64_t blockStartPos, u_int64_t blockSize, int sequenceCount, int blockID);
		void execute();
	
		string _buffer;
		bool _finished;
		
	private:
	

	
	
		RangeDecoder _rangeDecoder;
	//	ifstream* _inputFile;
		//ofstream* _outputFile;
		u_int64_t _blockStartPos;
		u_int64_t _blockSize;
	//	int _decodedSequenceCount;
		string _currentSeq;
		ifstream* _anchorDictFile;
		
		void decodeAnchorRead();
		kmer_type extendAnchor(kmer_type kmer, int pos, bool rightExtend);
		
		void decodeNoAnchorRead();
		void endRead();
		
		int _sequenceCount;
	
	tools::storage::impl::Group *  _group;
	tools::storage::impl::Storage::istream *_inputStream;

};

class QualDecoder
{
public:
	//QualDecoder(Leon* leon, const string& inputFilename);
	QualDecoder(Leon* leon, const string& inputFilename,tools::storage::impl::Group *  group);

	~QualDecoder();
	
	void setup(u_int64_t blockStartPos, u_int64_t blockSize, int sequenceCount);
	void setup( int blockID);

	void execute();
	
	string _buffer;
	bool _finished;
	
	
private:
	Leon* _leon;

	tools::storage::impl::Group *  _group;
	
	char * _inbuffer;
	//ifstream* _inputFile;
	
	tools::storage::impl::Storage::istream *_inputStream;

	
	//ofstream* _outputFile;
	u_int64_t _blockStartPos;
	u_int64_t _blockSize;
	//int _decodedSequenceCount;
	string _currentSeq;
	int _sequenceCount;
	int _processedSequenceCount;

};
#endif /* _DNACODER_HPP_ */

