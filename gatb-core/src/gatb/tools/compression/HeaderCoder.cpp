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

#include "HeaderCoder.hpp"

#include <bitset> //////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! delete
/*
#define PRINT_DEBUG_ENCODER
#define PRINT_DEBUG_DECODER
*/
		

//====================================================================================
// ** AbstractHeaderCoder
//====================================================================================
AbstractHeaderCoder::AbstractHeaderCoder(Leon* leon) :
_headerSizeModel(256)
{
	_leon = leon;
	_prevHeader = "";
	_currentHeader = "";
}
	
void AbstractHeaderCoder::addFieldColumn(){
	
	_typeModel.push_back(Order0Model(HEADER_TYPE_COUNT+1));
	_fieldIndexModel.push_back(Order0Model(256));
	_fieldColumnModel.push_back(Order0Model(256));
	_misSizeModel.push_back(Order0Model(256));
	_asciiModel.push_back(Order0Model(128));
	_zeroModel.push_back(Order0Model(256));
	
	_numericModels.push_back(vector<Order0Model>());
	for(int j=0; j<CompressionUtils::NB_MODELS_PER_NUMERIC; j++)
		_numericModels[_numericModels.size()-1].push_back( Order0Model(256) );
		
	_prevFieldPos.push_back(0);
	_currentFieldPos.push_back(0);
	_prevFieldValues.push_back(0);
	_currentFieldValues.push_back(0);
	_prevFieldTypes.push_back(FIELD_ASCII);
	_currentFieldTypes.push_back(FIELD_ASCII);
	_prevFieldZeroValues.push_back(0);
	_currentFieldZeroValues.push_back(0);
}

int AbstractHeaderCoder::typeOfChar(u_int8_t c, bool* isDigit){
	if(isdigit(c)){
		*isDigit = true;
		return 1;
	}
	else if(isalpha(c)){
		*isDigit = false;
		return 1;
	}
	else{
		*isDigit = false;
		return 2;
	}
}

void AbstractHeaderCoder::splitHeader(){
	_fieldIndex = 0;
	_fieldStartPos = 0;
	_currentPos = 0;
	_isCurrentFieldNumeric = true;
	
	u_int8_t c;
	int charType;
	bool digitOnly;
	int lastCharType = typeOfChar(_currentHeader[0], &digitOnly);
	
	for(_currentPos=0; (unsigned)_currentPos<_currentHeader.size(); _currentPos++){
		c = _currentHeader[_currentPos];
		
		digitOnly = true;
		charType = typeOfChar(c, &digitOnly);
		
		if(charType != lastCharType){
			lastCharType = charType;
			makeField();
		}
		
		if(_isCurrentFieldNumeric){
			_isCurrentFieldNumeric = digitOnly;
		}
	}
	
	makeField();
	
	_currentFieldCount = _fieldIndex;
}

void AbstractHeaderCoder::makeField(){
	if(_fieldStartPos == _currentPos) return;
	
	//Adjust the maximum number fo field column
	int currentFieldColumn = _currentFieldPos.size();
	while(currentFieldColumn <= _fieldIndex+1){
		addFieldColumn();
		currentFieldColumn = _currentFieldPos.size();
	}
	
	_currentFieldPos[_fieldIndex] = _fieldStartPos;
	_currentFieldPos[_fieldIndex+1] = _currentPos;
	
	if(_isCurrentFieldNumeric){
		string field = _currentHeader.substr(_currentFieldPos[_fieldIndex], _currentFieldPos[_fieldIndex+1]-_currentFieldPos[_fieldIndex]);
		
		int zeroCount = 0;
		if(field[0] == '0'){
			while(field[0] == '0'){
				zeroCount += 1;
				field.erase(field.begin());
			}
		}
		
		//cout << " HAHAHAHAAHAHHAHA    " << field << endl; 
		_currentFieldZeroValues[_fieldIndex] = zeroCount;
		
		u_int64_t value = strtoul(field.c_str(), NULL, 0);
		_currentFieldValues[_fieldIndex] = value;
		
		if(zeroCount > 0){
			if(value == 0){
				_currentFieldTypes[_fieldIndex] = FIELD_ZERO_ONLY;
			}
			else{
				_currentFieldTypes[_fieldIndex] = FIELD_ZERO_AND_NUMERIC;
			}
		}
		else {
			_currentFieldTypes[_fieldIndex] = FIELD_NUMERIC;
		}
		//cout << "OOOOOOOOOOOOoo   " << strtoul("15313", NULL, 0) << endl;
	}
	else{
		_currentFieldTypes[_fieldIndex] = FIELD_ASCII;
	}
	
	#ifdef PRINT_DEBUG_ENCODER
		cout << "\tField "<< _fieldIndex << ": " << _currentHeader.substr(_currentFieldPos[_fieldIndex], _currentFieldPos[_fieldIndex+1]-_currentFieldPos[_fieldIndex]) << "  Digit only? " << _isCurrentFieldNumeric << endl;
	#endif
	_fieldIndex += 1;
	_fieldStartPos = _currentPos;
	_isCurrentFieldNumeric = true;
}

void AbstractHeaderCoder::endHeader(){
	_prevFieldCount = _currentFieldCount;
	
	#ifdef PRINT_DEBUG_ENCODER
		cout << "\tField count: " << _prevFieldCount << endl;
	#endif
	
	for(int i=0; i<_prevFieldCount+1; i++){
		_prevFieldPos[i] = _currentFieldPos[i];
		_prevFieldValues[i] = _currentFieldValues[i];
		_prevFieldTypes[i] = _currentFieldTypes[i];
		_prevFieldZeroValues[i] = _currentFieldZeroValues[i];
		
		_currentFieldZeroValues[i] = 0;
	}
	_prevHeader = _currentHeader;
	_misIndex = 0;
	_fieldIndex = 0;
	
	_processedSequenceCount += 1;
}

void AbstractHeaderCoder::startBlock(){

	_currentHeader = _leon->_firstHeader;
	
	for(int i=0; (unsigned)i<_typeModel.size(); i++){
		_typeModel[i].clear();
		_fieldIndexModel[i].clear();
		_fieldColumnModel[i].clear();
		_misSizeModel[i].clear();
		_asciiModel[i].clear();
		_zeroModel[i].clear();
		
		for(int j=0; j<8; j++)
			_numericModels[i][j].clear();
			
	}
	_headerSizeModel.clear();
	
	splitHeader();
	endHeader();
	
	_processedSequenceCount = 0;
}

//====================================================================================
// ** HeaderEncoder
//====================================================================================
HeaderEncoder::HeaderEncoder(Leon* leon) :
AbstractHeaderCoder(leon) , _totalHeaderSize(0) ,_seqId(0)
{
	_thread_id = __sync_fetch_and_add (&_leon->_nb_thread_living, 1);

	
	//_firstHeader = firstHeader;
	//_rangeEncoder = new RangeEncoder();
}

HeaderEncoder::HeaderEncoder(const HeaderEncoder& copy) :
AbstractHeaderCoder(NULL), _totalHeaderSize(0),_seqId(0)
{

	
	_leon = copy._leon;
	
	_thread_id = __sync_fetch_and_add (&_leon->_nb_thread_living, 1);
	startBlock();

	//_firstHeader = copy._firstHeader;
	//_rangeEncoder = new RangeEncoder();
}

HeaderEncoder::~HeaderEncoder(){
	
	
	if( _thread_id!=0 && (_seqId+1) % _leon->getReadPerBlock() != 0 ){
		writeBlock();
	}
//	int nb_remaining =
 	__sync_fetch_and_add (&_leon->_nb_thread_living, -1);
	
	__sync_fetch_and_add(&_leon->_totalHeaderSize, _totalHeaderSize);

}


//int HeaderEncoder::getId(){
//	return ((_lastSequenceIndex / Leon::READ_PER_BLOCK) % _leon->_nb_cores);
//}

void HeaderEncoder::operator()(Sequence& sequence){
	_lastSequenceIndex = sequence.getIndex();
	_seqId = sequence.getIndex() ;


	_currentHeader = sequence.getComment();
	
	_totalHeaderSize += _currentHeader.size();
	
	processNextHeader();
	
	
	if(_processedSequenceCount >= _leon->getReadPerBlock() ){
		
		writeBlock();
		startBlock();
	}
	
}

void HeaderEncoder::writeBlock(){
	if(_rangeEncoder.getBufferSize() > 0){
		_rangeEncoder.flush();
	}
	
	int blockId = (  _seqId / _leon->getReadPerBlock())   ;

	//printf("\nheader coder writeblock   bid %i   tid %i \n",blockId, _thread_id);
	
	_leon->writeBlock(_rangeEncoder.getBuffer(), _rangeEncoder.getBufferSize(), _processedSequenceCount,blockId,true);
	_rangeEncoder.clear();
}

void HeaderEncoder::processNextHeader(){
	
	#ifdef PRINT_DEBUG_ENCODER
		cout << _prevHeader << endl;
		cout << _currentHeader << endl;
	#endif
	splitHeader();
	compareHeader();
	endHeader();
}

void HeaderEncoder::compareHeader(){
	_fieldPos = 0;
	_misCurrentStartPos = -1;
	
	for(_fieldIndex=0; _fieldIndex<_currentFieldCount; _fieldIndex++){
		
		_currentFieldSize = _currentFieldPos[_fieldIndex+1] - _currentFieldPos[_fieldIndex];
		_prevFieldSize = _prevFieldPos[_fieldIndex+1]-_prevFieldPos[_fieldIndex];
		_misCurrentStartPos = -1;
		
		HeaderType prevFieldType = _prevFieldTypes[_fieldIndex];
		HeaderType currentFieldType = _currentFieldTypes[_fieldIndex];
		
		//Comparing numeric field
		if(prevFieldType == FIELD_NUMERIC && currentFieldType == FIELD_NUMERIC){
			#ifdef PRINT_DEBUG_ENCODER
				cout << "\t\tComparing numeric fields: " <<_prevFieldValues[_fieldIndex] << " " << _currentFieldValues[_fieldIndex] << endl;
			#endif
			if(_prevFieldValues[_fieldIndex] == _currentFieldValues[_fieldIndex]){ //match
				_lastMatchFieldIndex = _fieldIndex;
				continue;
			}
			//encodeNumeric();
		}
		//Comparing field with zero only
		else if(prevFieldType == FIELD_ZERO_ONLY && currentFieldType == FIELD_ZERO_ONLY){
			#ifdef PRINT_DEBUG_ENCODER
				cout << "\t\tComparing fields with zero only: "  << endl;
			#endif
			if(_prevFieldZeroValues[_fieldIndex] == _currentFieldZeroValues[_fieldIndex]){ //match
				_lastMatchFieldIndex = _fieldIndex;
				continue;
			}
			//encodeNumeric();
		}
		//Comparing field with zero at begining and numeric
		else if(prevFieldType == FIELD_ZERO_AND_NUMERIC && currentFieldType == FIELD_ZERO_AND_NUMERIC){
			#ifdef PRINT_DEBUG_ENCODER
				cout << "\t\tComparing fields with zero at begining and numeric: "  << endl;
			#endif
			if(_prevFieldZeroValues[_fieldIndex] == _currentFieldZeroValues[_fieldIndex] && _prevFieldValues[_fieldIndex] == _currentFieldValues[_fieldIndex]){ //match
				_lastMatchFieldIndex = _fieldIndex;
				continue;
			}
			//encodeNumeric();
		}
		
		
		//Encoding numeric field
		if(currentFieldType == FIELD_NUMERIC || currentFieldType == FIELD_ZERO_ONLY || currentFieldType == FIELD_ZERO_AND_NUMERIC){
			encodeNumeric();
		}
		//Comparing ascii fields
		else{
			for(_fieldPos=0; _fieldPos<_currentFieldSize; _fieldPos++){
			
				if(_fieldIndex >= _prevFieldCount){
					_misCurrentStartPos = _fieldPos;
					break;
				}
					
				if(_fieldPos >= _prevFieldSize){
					_misCurrentStartPos = _fieldPos;
					break;
				}

				u_int8_t c = _currentHeader[_currentFieldPos[_fieldIndex]+_fieldPos];
				
				#ifdef PRINT_DEBUG_ENCODER
					cout << "\t\tComparing: " << _prevHeader[_prevFieldPos[_fieldIndex]+_fieldPos] << " " << c << "  ";
				#endif
				if(c == _prevHeader[_prevFieldPos[_fieldIndex]+_fieldPos]){
					#ifdef PRINT_DEBUG_ENCODER
						cout << "match" << endl;
					#endif
				}
				else{
					#ifdef PRINT_DEBUG_ENCODER
						cout << "mismatch" << endl;
					#endif
					_misCurrentStartPos = _fieldPos;
					break;
				}
			}
			
			if(_misCurrentStartPos != -1){
				encodeAscii();
			}
			else if(_fieldPos != _prevFieldSize){ //All the character of the current field match but there are always character in the current field of the prev header
				_misCurrentStartPos = _fieldPos;
				encodeAscii();
			}
			else{
				_lastMatchFieldIndex = _fieldIndex;
			}
			
		}
		
	}
	
	//if(_currentFieldPos[_fieldIndex]+_fieldPos == _currentHeader.size()){
	//	_misCurrentStartPos = _fieldPos;
	//	encodeMismatch();
	//}
	
	//cout << _lastMatchFieldIndex << " " << _fieldIndex << endl;
	
	//if the last field match, we have to signal to the decoder to add the last matching field of the prev header
	
	if(_lastMatchFieldIndex == _fieldIndex-1){
		_rangeEncoder.encode(_typeModel[_misIndex], HEADER_END_MATCH);
		_rangeEncoder.encode(_headerSizeModel, _currentHeader.size());
		//_misCurrentStartPos = _currentFieldSize;
		//encodeAscii();
	}
	else{
		_rangeEncoder.encode(_typeModel[_misIndex], HEADER_END);
	}
	//_misIndex += 1;
	
	//end of header
	//_rangeEncoder.encode(_typeModel[_misIndex], HEADER_END);
	//_rangeEncoder.encode(_headerSizeModel, _currentHeader.size());
	
}


void HeaderEncoder::encodeNumeric(){
	u_int64_t zeroCount = _currentFieldZeroValues[_fieldIndex];
	u_int64_t fieldValue = _currentFieldValues[_fieldIndex];
	
	HeaderType currentFieldType = _currentFieldTypes[_fieldIndex];
	
	if(currentFieldType == FIELD_ZERO_ONLY){
		#ifdef PRINT_DEBUG_ENCODER
			cout << "\t\t\tField with zero only" << endl;
			cout << "\t\t\tEnconding zero count: " << zeroCount << endl;
		#endif
		_rangeEncoder.encode(_typeModel[_misIndex], FIELD_ZERO_ONLY);
		_rangeEncoder.encode(_fieldIndexModel[_misIndex], _fieldIndex);
		_rangeEncoder.encode(_zeroModel[_misIndex], zeroCount);
		_misIndex += 1;
		return;
	}
	else if(currentFieldType == FIELD_ZERO_AND_NUMERIC){
		#ifdef PRINT_DEBUG_ENCODER
			cout << "\t\t\tField with zero and numeric" << endl;
			cout << "\t\t\tEnconding zero count: " << zeroCount << endl;
		#endif
		_rangeEncoder.encode(_typeModel[_misIndex], FIELD_ZERO_AND_NUMERIC);
		_rangeEncoder.encode(_fieldIndexModel[_misIndex], _fieldIndex);
		_rangeEncoder.encode(_zeroModel[_misIndex], zeroCount);
		_misIndex += 1;
	}

	u_int64_t value = fieldValue;
	u_int64_t prevValue = _prevFieldValues[_fieldIndex];


	//int valueByteCount = CompressionUtils::getByteCount(value);
	
	#ifdef PRINT_DEBUG_ENCODER
		cout << "\t\t\tPrev value: " << prevValue << endl;
		cout << "\t\t\tField value: " << value << "    Byte: " << valueByteCount << endl;
	#endif
	
	u_int64_t deltaValue;
	int deltaType = CompressionUtils::getDeltaValue(value, prevValue, &deltaValue);
	
	if(deltaType == 0){
		_rangeEncoder.encode(_typeModel[_misIndex], FIELD_NUMERIC);
	}
	else if(deltaType == 1){
		_rangeEncoder.encode(_typeModel[_misIndex], FIELD_DELTA);
		value = deltaValue;
	}
	else if(deltaType == 2){
		_rangeEncoder.encode(_typeModel[_misIndex], FIELD_DELTA_2);
		value = deltaValue;
	}
	
	  
	_rangeEncoder.encode(_fieldIndexModel[_misIndex], _fieldIndex);
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericModels[_misIndex], value);
	//_rangeEncoder->encode(&_fieldColumnModel[_misIndex], 0);
	//_prevFieldValues[_fieldIndex] = fieldValue;
	
	_misIndex += 1;
}

void HeaderEncoder::encodeAscii(){
	int missSize = _currentFieldSize - _misCurrentStartPos;//_currentPos - _misCurrentStartPos;
	//cout << _currentFieldSize <<  " " << _fieldPos << endl;
	_rangeEncoder.encode(_typeModel[_misIndex], FIELD_ASCII);
	_rangeEncoder.encode(_fieldIndexModel[_misIndex], _fieldIndex);
	_rangeEncoder.encode(_fieldColumnModel[_misIndex], _misCurrentStartPos);
	_rangeEncoder.encode(_misSizeModel[_misIndex], missSize);
	#ifdef PRINT_DEBUG_ENCODER
		cout << "\t\t\t<Mismatch> " << "    Type: " << "ASCII" << "    Field: " << _fieldIndex << "    Column: " << _misCurrentStartPos << "    Size: " << missSize << endl;
	#endif
	//for(int j=_misCurrentStartPos; j<_currentPos; j++){
	for(int i=_misCurrentStartPos; i < _misCurrentStartPos+missSize; i++){
		#ifdef PRINT_DEBUG_ENCODER
			cout << "\t\t\tEncoding: " << _currentHeader[_currentFieldPos[_fieldIndex]+i] << endl;
		#endif
		//cout << _currentHeader[j] << flush;
		_rangeEncoder.encode(_asciiModel[_misIndex], _currentHeader[_currentFieldPos[_fieldIndex]+i]);
	}
	//cout << endl;
	_misIndex += 1;
}











//====================================================================================
// ** HeaderDecoder
//====================================================================================
HeaderDecoder::HeaderDecoder(Leon* leon,std::string & inputFilename, tools::storage::impl::Group *  group) :
AbstractHeaderCoder(leon)
//, _rangeDecoder(inputFile)
{
	_group = group;
	_inputStream =0;
	_finished = false;

}

HeaderDecoder::~HeaderDecoder(){

	if(_inputStream !=0) delete _inputStream;

}

void HeaderDecoder::setup(u_int64_t blockStartPos, u_int64_t blockSize, int sequenceCount,int blockID){
	startBlock();
	_rangeDecoder.clear();
	
	
		if(_inputStream !=0) delete _inputStream;
	std::string datasetname = Stringify::format ("header_%i",blockID);
	
	_inputStream = new tools::storage::impl::Storage::istream  (*_group, datasetname);
	
	auto _tempcollec = & _group->getCollection<math::NativeInt8> (datasetname);
	std::string dsize = _tempcollec->getProperty ("size");
	
	_blockSize =  std::stoi(dsize); // blockSize;
	
	
	
	_rangeDecoder.setInputFile(_inputStream);

	
	
	
	_blockStartPos = blockStartPos;
	//_blockSize = blockSize;
	
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t-----------------------" << endl;
		cout << "\tDecoding block " << _blockStartPos << " - " << _blockStartPos+_blockSize << endl;
	#else
		//_leon->_progress_decode->inc(1);

	#endif

	_currentHeader.clear();
	_misIndex = 0;
	
	_sequenceCount = sequenceCount;
}

void HeaderDecoder::execute(){
	//cout << "executing" << endl;
	//decodeFirstHeader();
	
	while(_processedSequenceCount < _sequenceCount){

		u_int8_t type = _rangeDecoder.nextByte(_typeModel[_misIndex]);
		#ifdef PRINT_DEBUG_DECODER
			cout << "\t\tNext type is: " << (int)type << endl;
		#endif

		
		if(type == HEADER_END){

			endHeader();
			//i+=1;
		}
		else if(type == HEADER_END_MATCH){
			//decodeMatch();
			u_int8_t headerSize = _rangeDecoder.nextByte(_headerSizeModel);

			for(/*_fieldIndex*/; _fieldIndex < _prevFieldCount; _fieldIndex++){
				#ifdef PRINT_DEBUG_DECODER
					cout << "\t\t\tAdding from prev header: " << _prevHeader.substr(_prevFieldPos[_fieldIndex], _prevFieldPos[_fieldIndex+1]-_prevFieldPos[_fieldIndex]) << endl;
				#endif
				_currentHeader += _prevHeader.substr(_prevFieldPos[_fieldIndex], _prevFieldPos[_fieldIndex+1]-_prevFieldPos[_fieldIndex]);
				if(_currentHeader.size() >= headerSize) break;
			}
			
			endHeader();
		}
		else{
			
			decodeMatch();
			
			if(type == FIELD_ASCII){
				decodeAscii();
				_fieldIndex += 1;
				_misIndex += 1;
			}
			else if(type == FIELD_NUMERIC){
				decodeNumeric();
				_fieldIndex += 1;
				_misIndex += 1;
			}
			else if(type == FIELD_DELTA){
				decodeDelta();
				_fieldIndex += 1;
				_misIndex += 1;
			}
			else if(type == FIELD_DELTA_2){
				decodeDelta2();
				_fieldIndex += 1;
				_misIndex += 1;
			}
			else if(type == FIELD_ZERO_ONLY){
				decodeZero();
				_fieldIndex += 1;
				_misIndex += 1;
			}
			else if(type == FIELD_ZERO_AND_NUMERIC){
				decodeZero();
				_misIndex += 1;
				//decodeNumeric();
				//_fieldIndex += 1;
			}
			//_prevPos = _prevFieldPos[_fieldIndex+1];

		}
		
		//cout << "lala" << endl;
	}
	
	_finished = true;
}


void HeaderDecoder::decodeMatch(){
	u_int8_t misFieldIndex = _rangeDecoder.nextByte(_fieldIndexModel[_misIndex]);
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tMatch to field: " << (int)misFieldIndex << endl;
	#endif
	for(/*_fieldIndex*/; _fieldIndex < misFieldIndex; _fieldIndex++){
		#ifdef PRINT_DEBUG_DECODER
			cout << "\t\t\tAdding from prev header: " << _prevHeader.substr(_prevFieldPos[_fieldIndex], _prevFieldPos[_fieldIndex+1]-_prevFieldPos[_fieldIndex]) << endl;
		#endif
		_currentHeader += _prevHeader.substr(_prevFieldPos[_fieldIndex], _prevFieldPos[_fieldIndex+1]-_prevFieldPos[_fieldIndex]);
	}
}

void HeaderDecoder::decodeAscii(){
	u_int8_t misColumn = _rangeDecoder.nextByte(_fieldColumnModel[_misIndex]);
	u_int8_t misSize = _rangeDecoder.nextByte(_misSizeModel[_misIndex]);
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tDecoding   Type: ASCII     Column: " << (int)misColumn << "    Size: " << (int)misSize << endl;
	#endif
	
	if(_fieldIndex < _prevFieldCount){
		for(int fieldPos=0; fieldPos<misColumn; fieldPos++){
			_currentHeader += _prevHeader[_prevFieldPos[_fieldIndex]+fieldPos];
		}
	}
	
	for(int i=0; i<misSize; i++){
		u_int8_t c = _rangeDecoder.nextByte(_asciiModel[_misIndex]);
		
		#ifdef PRINT_DEBUG_DECODER
			cout << "\t\t\tAdding: " << c << " (" << (int)c << ")"<< endl;
		#endif
		//_currentHeader2[_currentPos] = c;
		_currentHeader += c;
		//_currentPos += 1;
	}
	
}	

void HeaderDecoder::decodeNumeric(){
	//u_int8_t misSize = _rangeDecoder.nextByte(_misSizeModel[_misIndex]);
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tDecoding   Type: NUMERIC" << endl; //"    Size: " << (int)misSize << endl;
	#endif
	
	u_int64_t value = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModels[_misIndex]);
	//_currentHeader += CompressionUtils::numberToString(value);
	
	char temp[200];
	snprintf(temp,200,"%llu",value);
	_currentHeader += string(temp);
	//_currentHeader += to_string(value); // C++11
	
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\t\tAdding: " << string(temp) << endl;
	#endif
}

void HeaderDecoder::decodeDelta(){
	//u_int8_t misSize = _rangeDecoder.nextByte(_misSizeModel[_misIndex]);
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tDecoding   Type: DELTA" << endl;//"    Size: " << (int)misSize << endl;
	#endif
	
	u_int64_t value = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModels[_misIndex]);

	value = CompressionUtils::getValueFromDelta(1, _prevFieldValues[_fieldIndex], value);
	
	char temp[200];
	snprintf(temp,200,"%llu",value);
	_currentHeader += string(temp);
	//_currentHeader += to_string(value);
	
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\t\tAdding: " << string(temp) << endl;
	#endif
}

void HeaderDecoder::decodeDelta2(){
	//u_int8_t misSize = _rangeDecoder.nextByte(_misSizeModel[_misIndex]);
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tDecoding   Type: DELTA 2" << endl;//"    Size: " << (int)misSize << endl;
	#endif
	
	u_int64_t value = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModels[_misIndex]);

	value = CompressionUtils::getValueFromDelta(2, _prevFieldValues[_fieldIndex], value);
	char temp[200];
	snprintf(temp,200,"%llu",value);
	_currentHeader += string(temp);
	//_currentHeader += to_string(value);
	
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\t\tAdding: " << string(temp) << endl;
	#endif
}

void HeaderDecoder::decodeZero(){
	u_int8_t zeroCount = _rangeDecoder.nextByte(_zeroModel[_misIndex]);
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tDecoding   Type: ZERO     Size: " << (int)zeroCount << endl;
	#endif
	
	for(int i=0; i<zeroCount; i++){
		
		#ifdef PRINT_DEBUG_DECODER
			cout << "\t\t\tAdding: 0"<< endl;
		#endif
		
		_currentHeader += '0';
	}
}

void HeaderDecoder::endHeader(){
	_buffer += _currentHeader + '\n';

	#ifdef PRINT_DEBUG_DECODER
		cout << _currentHeader << endl;
		//for(int i=0; i<_currentPos; i++){
		//	cout << _currentHeader2[i];
		//}
		//cout << endl;
	#endif
	
	
	splitHeader();
	AbstractHeaderCoder::endHeader();
	_currentHeader.clear();
	_misIndex = 0;
}


