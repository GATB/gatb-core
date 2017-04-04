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

#ifndef _HEADERCODER_HPP_
#define _HEADERCODER_HPP_

/*
#include <string>
#include <vector>
#include <ctype.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>*/

//#include "RangeCoder.hpp"
#include <gatb/gatb_core.hpp>
#include "Leon.hpp"
//#include "CompressionUtils.hpp"

using namespace std;
class Leon;
//====================================================================================
// ** AbstractHeaderCoder
//====================================================================================
class AbstractHeaderCoder
{
	public:
		AbstractHeaderCoder(Leon* leon);
		
	protected:
		void addFieldColumn();
	
		enum HeaderType{HEADER_END=1, HEADER_END_MATCH, FIELD_ASCII, FIELD_NUMERIC, FIELD_DELTA, FIELD_DELTA_2, FIELD_ZERO_ONLY, FIELD_ZERO_AND_NUMERIC, HEADER_TYPE_COUNT};
		//static const int MAX_FIELD_COUNT = 200;
		
		vector<Order0Model> _typeModel;
		vector<Order0Model> _fieldIndexModel;
		vector<Order0Model> _fieldColumnModel;
		vector<Order0Model> _misSizeModel;
		vector<Order0Model> _asciiModel;
		vector< vector<Order0Model> > _numericModels;
		vector<Order0Model> _zeroModel;
		
		Order0Model _headerSizeModel;
		
		int typeOfChar(u_int8_t c, bool* isDigit);
		void splitHeader();
		void makeField();
		void endHeader();
		
		string _prevHeader;
		string _currentHeader;
		vector<unsigned int> _prevFieldPos;
		vector<unsigned int> _currentFieldPos;
		int _currentPos;
		int _fieldStartPos;
		int _prevFieldCount;
		int _fieldIndex;
		int _misIndex;
	
		vector<u_int64_t> _prevFieldValues;
		vector<u_int64_t> _currentFieldValues;
		vector<u_int64_t> _prevFieldZeroValues;
		vector<u_int64_t> _currentFieldZeroValues;
		vector<HeaderType> _prevFieldTypes;
		vector<HeaderType> _currentFieldTypes;
		
		
		bool _isCurrentFieldNumeric;
		int _currentFieldCount;
		
		Leon* _leon;
		
		void startBlock();
		
		int _processedSequenceCount;
};

//====================================================================================
// ** HeaderEncoder
//====================================================================================
class HeaderEncoder : AbstractHeaderCoder
{
		
	public:
		
		HeaderEncoder(Leon* leon);
		HeaderEncoder(const HeaderEncoder& copy);
		~HeaderEncoder();

		void operator()(Sequence& sequence);
		//int getId();
		u_int64_t _lastSequenceIndex;
		
	private:
		
		RangeEncoder _rangeEncoder;
		
		u_int64_t _totalHeaderSize;
	
		int _fieldPos;
		//int _misPrevStartPos, _misCurrentStartPos;
		int _misCurrentStartPos;
		//int _encoderFieldIndex;
		int _prevFieldSize, _currentFieldSize;
		int _lastMatchFieldIndex;
		u_int64_t _seqId;
		int _thread_id;

		//static void encodeFirstHeader();
		void writeBlock();
		
		void processNextHeader();
		void compareHeader();
		//void encode();
		//void encodeMismatch();
		void encodeNumeric();
		void encodeAscii();
		
};

//====================================================================================
// ** HeaderDecoder
//====================================================================================
class HeaderDecoder : AbstractHeaderCoder
{
		
	public:
		
		HeaderDecoder(Leon* leon, std::string & inputFilename, tools::storage::impl::Group *  group);
		~HeaderDecoder();
		
		//void processNextByte(u_int8_t byte);
		void setup(u_int64_t blockStartPos, u_int64_t blockSize, int sequenceCount, int blockID);
		void execute();
	
		string _buffer;
		bool _finished;
		
	private:
	tools::storage::impl::Group *  _group;

		RangeDecoder _rangeDecoder;
		//ifstream* _inputFile;
	tools::storage::impl::Storage::istream *_inputStream;

		//ofstream* _outputFile;
		u_int64_t _blockStartPos;
		u_int64_t _blockSize;
		
		//int _prevPos;
		void endHeader();
		//void decodeFirstHeader();
		void decodeMatch();
		void decodeAscii();
		void decodeNumeric();
		void decodeDelta();
		void decodeDelta2();
		void decodeZero();
		

		int _sequenceCount;
		
};

#endif /* _HEADERCODER_HPP_ */

