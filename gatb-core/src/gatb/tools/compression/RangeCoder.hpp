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


#ifndef _RANGECODER_HPP_
#define _RANGECODER_HPP_

#include <string>
#include <vector>
#include <ctype.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <sys/types.h>
//#include "File.hpp"

using namespace std;

const u_int64_t TOP = (u_int64_t) 1<<56;
const u_int64_t BOTTOM = (u_int64_t) 1<<48;
const u_int64_t MAX_RANGE = BOTTOM;



class Order0Model
{
	
	public:
		Order0Model(int charCount);
		~Order0Model();
		
		void clear();
		void update(uint8_t c);
		
		u_int64_t rangeLow(uint8_t c);
		u_int64_t rangeHigh(uint8_t c);
		u_int64_t totalRange();
		unsigned int charCount();
	
	private:
		vector<u_int64_t> _charRanges;
		int _charCount;
		
		void rescale();
		
};

//====================================================================================
// ** AbstractRangeCoder
//====================================================================================
class AbstractRangeCoder
{
	
	public:
		AbstractRangeCoder();
		//~AbstractRangeCoder();

	protected:
		u_int64_t _low;
		u_int64_t _range;
		
};


//====================================================================================
// ** RangeEncoder
//====================================================================================
class RangeEncoder : AbstractRangeCoder
{
	public:
		
		RangeEncoder();//(ofstream& outputFile);
		~RangeEncoder();
		
		void encode(Order0Model& model, uint8_t c);
		void flush();
		void clear();
		void clearBuffer();
		u_int8_t* getBuffer(bool reversed=false);
		u_int64_t getBufferSize();
		
		bool updateModel; //Used by leon when PRINT_STAT macro is defined
		
	private:
		vector<u_int8_t> _buffer;
		vector<u_int8_t> _reversedBuffer;
};

//====================================================================================
// ** RangeDecoder
//====================================================================================
class RangeDecoder : AbstractRangeCoder
{
	public:
		
		RangeDecoder();
		~RangeDecoder();
		
		void setInputFile(ifstream* inputFile, bool reversed=false);
		u_int8_t nextByte(Order0Model& model);
		void clear();
		
	private:
		
		ifstream* _inputFile;
		u_int64_t _code;
		bool _reversed;
		
		u_int64_t getCurrentCount(Order0Model& model);
		void removeRange(Order0Model& model, uint8_t c);
		u_int8_t getNextByte();
};



#endif /* _RANGECODER_HPP_ */

