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

#include "RangeCoder.hpp"
#include <algorithm>    // std::reverse

//#define PRINT_DEBUG_RANGE_CODER

//====================================================================================
// ** Order0Model
//====================================================================================
Order0Model::Order0Model(int charCount){
	_charCount = charCount+1;
	
	for(int i=0; i<_charCount; i++){
		_charRanges.push_back(i);
	}
}

Order0Model::~Order0Model(){
}

void Order0Model::clear(){
	for(int i=0; i<_charCount; i++){
		_charRanges[i] = i;
	}
	
}

void Order0Model::update(uint8_t c){
	//cout << rangeLow(c) <<  " " << rangeHigh(c) << " " << totalRange()<< endl;
	for(int i=c+1; i<_charCount; i++){
		_charRanges[i] += 1;
	}
	if(_charRanges[_charCount-1] >= MAX_RANGE){
		rescale();
	}
}

void Order0Model::rescale(){
	for(int i=1; i<_charCount; i++) {
		_charRanges[i] /= 2;
		if(_charRanges[i] <= _charRanges[i-1]){
			_charRanges[i] = _charRanges[i-1] + 1;
		}
	}
}

u_int64_t Order0Model::rangeLow(uint8_t c){
	return _charRanges[c];
}

u_int64_t Order0Model::rangeHigh(uint8_t c){
	return _charRanges[c+1];
}

u_int64_t Order0Model::totalRange(){
	return _charRanges[_charRanges.size()-1];
	//return _charRanges[_charCount-1];
}

unsigned int Order0Model::charCount(){
	return _charCount;
}

//====================================================================================
// ** AbstractRangeCoder
//====================================================================================
AbstractRangeCoder::AbstractRangeCoder(){
	_low = 0;
	_range = -1;
}


//====================================================================================
// ** RangeEncoder
//====================================================================================
RangeEncoder::RangeEncoder()//(ofstream& outputFile)
{	
	updateModel = true;
	//_outputFile = new OutputFile(outputFile);
}

RangeEncoder::~RangeEncoder(){
	//flush(); //ATTENTION PTET PB DE SYNCHRO ICI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//delete _outputFile;
}

void RangeEncoder::encode(Order0Model& model, uint8_t c){
	#ifdef PRINT_DEBUG_RANGE_CODER
		printf("\t\t\tencoding char: %c\n", c);
	#endif
	
	//cout << model->rangeHigh(c) -model->rangeLow(c) << endl;
	_range /= model.totalRange();
	_low += model.rangeLow(c) * _range;
	_range *= model.rangeHigh(c) - model.rangeLow(c);

	while((_low ^ _low+_range) < TOP || _range < BOTTOM ){
		if(_range < BOTTOM && ((_low ^ _low+_range)) >= TOP){
			_range = -_low & (BOTTOM-1);
		}
		
		_buffer.push_back(_low>>56);
		//cout << "RangeEncoder Output: " << (_low>>56) << endl;
		//_outputFile->writeByte(_low>>56);
		_range <<= 8;
		_low <<= 8;
	}
	
	#ifdef LEON_PRINT_STAT
		if(updateModel) model.update(c);
	#else
		model.update(c);
	#endif

}

void RangeEncoder::flush(){
	for(int i=0; i<8; i++){
		//cout << "RangeEncoder Output: " << (_low>>56) << endl;
		_buffer.push_back(_low>>56);
		//_outputFile->writeByte(_low>>56);
		_low <<= 8;
	}
	_low = 0;
	_range = -1;
}

void RangeEncoder::clear(){
	_low = 0;
	_range = -1;
	clearBuffer();
}

void RangeEncoder::clearBuffer(){
	_buffer.clear();
}

u_int8_t* RangeEncoder::getBuffer(bool reversed){
	if(reversed){
		_reversedBuffer.clear();
		_reversedBuffer.assign(_buffer.begin(), _buffer.end());
		reverse(_reversedBuffer.begin(), _reversedBuffer.end());
		return &_reversedBuffer[0];
	}
	else{
		return &_buffer[0];
	}
}

		
u_int64_t RangeEncoder::getBufferSize(){
	return _buffer.size();
}



//====================================================================================
// ** RangeDecoder
//====================================================================================
RangeDecoder::RangeDecoder()
{
}

RangeDecoder::~RangeDecoder(){
}
	
	
void RangeDecoder::setInputFile(ifstream* inputFile, bool reversed){
	_reversed = reversed;
	clear();
	_inputFile = inputFile;
	
	for(int i=0; i<8; i++){
		_code = (_code << 8) | getNextByte();
	}
}

uint8_t RangeDecoder::nextByte(Order0Model& model){
	u_int64_t count = getCurrentCount(model);
	uint8_t c;
	for(c=model.charCount()-2; model.rangeLow(c) > count; c--);

	removeRange(model, c);

	model.update(c);
	//cout << "RangeDecoder output: " << (int)c << endl;
	return c;
}
	
u_int64_t RangeDecoder::getCurrentCount(Order0Model& model){
	_range /= model.totalRange();
	return (_code-_low) / _range;
}

void RangeDecoder::removeRange(Order0Model& model, uint8_t c){
	_low += model.rangeLow(c) * _range;
	_range *= model.rangeHigh(c) - model.rangeLow(c);

	while((_low ^ _low+_range) < TOP || _range < BOTTOM ){
		if(_range < BOTTOM && ((_low ^ _low+_range)) >= TOP){
			_range = -_low & (BOTTOM-1);
		}
		
		_code = _code << 8 | getNextByte();
		_range <<= 8; 
		_low <<= 8;
	}
}

void RangeDecoder::clear(){
	_low = 0;
	_range = -1;
	_code = 0;
}

u_int8_t RangeDecoder::getNextByte(){
	u_int8_t byte;
	
	if(_reversed){
		_inputFile->seekg(-1, _inputFile->cur);
	}
	
	byte = _inputFile->get();
	
	if(_reversed){
		_inputFile->seekg(-1, _inputFile->cur);
	}
	
	return byte;
}

