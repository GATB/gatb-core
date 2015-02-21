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

#ifndef _COMPRESSIONUTILS_HPP_
#define _COMPRESSIONUTILS_HPP_

using namespace std;


//namespace gatb          {
//namespace core          {
//namespace tools         {
//namespace compression   {


//====================================================================================
// ** COMPRESSIONUTILS
//====================================================================================

class CompressionUtils
{
	public:
	
		//For better compression rate, each byte outputted when processing numerics is encoded using distinct models
		static const u_int8_t NB_MODELS_PER_NUMERIC = 20;
		
		
		//Take a value and return its number of byte
		static int getByteCount(u_int64_t n){
			if(n < 0x100){
				return 1;
			}
			else if(n >= 0x100 && n < 0x10000){
				return 2;
			}
			else if(n >= 0x10000 && n < 0x1000000){
				return 3;
			}
			else if(n >= 0x1000000 && n < 0x100000000){
				return 4;
			}
			else if(n >= 0x100000000 && n < 0x10000000000){
				return 5;
			}
			else if(n >= 0x10000000000 && n < 0x1000000000000){
				return 6;
			}
			else if(n >= 0x1000000000000 && n < 0x100000000000000){
				return 7;
			}
			else{
				return 8;
			}

		}

		//Variable byte enconding (VBE) http://nlp.stanford.edu/IR-book/html/htmledition/variable-byte-codes-1.html
		static void encodeNumeric(RangeEncoder& rangeEncoder, vector<Order0Model>& numericModels, u_int64_t value){

			int i=0;
			do
			{
				u_int8_t byteCode = (u_int8_t) (value & 127);
				value = value >> 7;
				if (value != 0)
				{
					byteCode |= 128;
					rangeEncoder.encode(numericModels[i], byteCode);
				}
				else
					rangeEncoder.encode(numericModels[i], byteCode);

				i += 1;
			}
			while (value != 0);


			//int valueByteCount = getByteCount(value);
			//cout << "Utils: " << (int)valueByteCount << endl;
			//rangeEncoder.encode(byteCountModel, valueByteCount);
				
			//for(int i=0; i<valueByteCount; i++){
				//cout << "Utils: " << ((value >> i*8) & 0xff) << endl;
				//rangeEncoder.encode(numericModels[i], (value >> i*8) & 0xff);
			//}
		}

		static void encodeFixedNumeric(RangeEncoder& rangeEncoder, vector<Order0Model>& numericModels, u_int64_t value, int byteCount){
			//cout << "sdf" << value << " " << byteCount << endl;
			for(int i=0; i<byteCount; i++){
				//cout << ((value >> i*8) & 0xff) << endl;
				rangeEncoder.encode(numericModels[i], (value >> i*8) & 0xff);
			}
		}

		//Variable byte enconding (VBE) http://nlp.stanford.edu/IR-book/html/htmledition/variable-byte-codes-1.html
		static u_int64_t decodeNumeric(RangeDecoder& rangeDecoder, vector<Order0Model>& numericModels){
			//u_int8_t byteCount = rangeDecoder.nextByte(byteCountModel);
			//cout << (int)byteCount << endl;
			//u_int64_t value = 0;
			//for(int i=0; i<byteCount; i++){
			//	u_int8_t byteValue = rangeDecoder.nextByte(numericModels[i]);
				//cout << "Utils: " << (byteValue << i*8) << endl;
			//	value |= ((u_int64_t)byteValue << i*8);
			//}
			int i = 0;
			u_int64_t value = 0;
			u_int64_t byteCode = 0;
			int shift = 0;
			do
			{
				byteCode = rangeDecoder.nextByte(numericModels[i]);
				value += (byteCode & 127) << shift;
				shift += 7;
				i += 1;
			}
			while (byteCode > 127);
			
			return value;
		}
		
		static u_int64_t decodeFixedNumeric(RangeDecoder& rangeDecoder, vector<Order0Model>& numericModels, int byteCount){
			
			u_int64_t value = 0;
			for(int i=0; i<byteCount; i++){
				
				u_int8_t byteValue = rangeDecoder.nextByte(numericModels[i]);
				value |= (byteValue << i*8);
			}
			
			return value;
		}
		
		/*
		static stringstream _convert;
		
		static string numberToString(u_int64_t number){
			_convert << number;
			string result(_convert.str());
			//cout << result << endl;
			return result;
		}*/
		
		
		static u_int8_t getDeltaValue(u_int64_t value, u_int64_t prevValue, u_int64_t* resultDeltaValue){
			
			bool isDelta1Valid = false;
			bool isDelta2Valid = false;

			u_int64_t deltaValue1 = value - prevValue;
			u_int64_t deltaValue2 = prevValue - value;
			
			
			if( deltaValue1 < value){ //deltaValue1 >= 0 &&
				isDelta1Valid = true;
			}
			if( deltaValue2 < value){ //deltaValue2 >= 0 &&
				isDelta2Valid = true;
			}
			
			if(isDelta1Valid && isDelta2Valid){
				if(deltaValue1 <= deltaValue2){
					*resultDeltaValue = deltaValue1;
					return 1;
				}
				else{
					*resultDeltaValue = deltaValue2;
					return 2;
				}
			}
			else if(isDelta1Valid){
				*resultDeltaValue = deltaValue1;
				return 1;
			}
			else if(isDelta2Valid){
				*resultDeltaValue = deltaValue2;
				return 2;
			}
			
			*resultDeltaValue = value;
			return 0;	
		}
	
	
	
		static u_int64_t getValueFromDelta(u_int8_t deltaType, u_int64_t prevValue, u_int64_t deltaValue){
			
			if(deltaType == 0){
				return deltaValue;
			}
			if(deltaType == 1){
				return prevValue + deltaValue;	
			}
			else if(deltaType == 2){
				return prevValue - deltaValue;
			}
			return deltaValue; //avoid warning;
		}
};

	

/********************************************************************************/
//}}}};
/********************************************************************************/
	

#endif /* _COMPRESSIONUTILS_HPP_ */

