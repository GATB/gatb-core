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

#include "DnaCoder.hpp"

/*
#define PRINT_DEBUG_EXTMUTA
#define PRINT_DEBUG_ENCODER
#define PRINT_DEBUG_DECODER
*/


char bin2NTrev[4] = {'T','G','A','C'};
//char bin2NT[4] = {'A','C','T','G'};

/*
GACGCGCCGATATAACGCGCTTTCCCGGCTTTTACCACGTCGTTGAGGGCTTCCAGCGTCTCTTCGATCGGCGTGTTGTAATCCCAGCGATGAATTTG6:2308q
			Anchor pos: 8
			Anchor: GATATAACGCGCTTTCCCGGCTTTTACCACG
			§ Anchor: 8
*/

//====================================================================================
// ** AbstractDnaCoder
//====================================================================================
AbstractDnaCoder::AbstractDnaCoder(Leon* leon) :
_kmerModel(leon->_kmerSize),
_readTypeModel(2), //only 2 value in this model: read with anchor or without anchor
//_isPrevReadAnchorableModel(2),
_noAnchorReadModel(5), _bifurcationModel(5), //5value: A, C, G, T, N
_bifurcationBinaryModel(2), //0 or 1 (lowest or highest bifurcation by alphabetical order)
_readAnchorRevcompModel(2),
_readSizeDeltaTypeModel(3),
_anchorPosDeltaTypeModel(3),
_anchorAddressDeltaTypeModel(3),
_NposDeltaTypeModel(3),
_errorPosDeltaTypeModel(3),_seqId(0)
{
	_leon = leon;
	_bloom = _leon->_bloom;
	_kmerSize = _leon->_kmerSize;
	
	
	for(int i=0; i<CompressionUtils::NB_MODELS_PER_NUMERIC; i++){
		_anchorAddressModel.push_back(Order0Model(256));
		//_isPrevReadAnchorablePosModel.push_back(Order0Model(256));
		_anchorPosModel.push_back(Order0Model(256));
		_noAnchorReadSizeValueModel.push_back(Order0Model(256));
		_readSizeValueModel.push_back(Order0Model(256));
		_NposModel.push_back(Order0Model(256));
		_leftErrorPosModel.push_back(Order0Model(256));
		//_rightErrorPosModel.push_back(Order0Model(256));
		_numericModel.push_back(Order0Model(256));
		_leftErrorModel.push_back(Order0Model(256));
		//_rightErrorModel.push_back(Order0Model(256));
	}
	
	
}

void AbstractDnaCoder::startBlock(){
	//_prevSequences = NULL;

	for(int i=0; i<CompressionUtils::NB_MODELS_PER_NUMERIC; i++){
		_anchorAddressModel[i].clear();
		//_isPrevReadAnchorablePosModel[i].clear();
		_anchorPosModel[i].clear();
		_noAnchorReadSizeValueModel[i].clear();
		_readSizeValueModel[i].clear();
		_NposModel[i].clear();
		_leftErrorPosModel[i].clear();
		//_rightErrorPosModel[i].clear();
		_numericModel[i].clear();
		_leftErrorModel[i].clear();
		//_rightErrorModel[i].clear();
	}
	_readTypeModel.clear();
	//_isPrevReadAnchorableModel.clear();
	_noAnchorReadModel.clear();
	_bifurcationModel.clear();
	_bifurcationBinaryModel.clear();
	_readAnchorRevcompModel.clear();
	_readSizeDeltaTypeModel.clear();
	_anchorPosDeltaTypeModel.clear();
	_anchorAddressDeltaTypeModel.clear();
	_NposDeltaTypeModel.clear();
	_errorPosDeltaTypeModel.clear();
	_prevReadSize = 0;
	_prevAnchorPos = 0;
	_prevAnchorAddress = 0;
	_prevNpos = 0;
	_prevErrorPos = 0;
	
	_processedSequenceCount = 0;
}

void AbstractDnaCoder::endRead(){
	_processedSequenceCount += 1;
}

void AbstractDnaCoder::codeSeedBin(KmerModel* model, kmer_type* kmer, int nt, bool right){

	
	if(right)
	{
        /** We initialize the kmer. */
        KmerModel::Kmer tmp;  tmp.set (*kmer);

		*kmer = model->codeSeedRight (tmp, nt, Data::INTEGER).value();
	}
	else
	{
        /** We initialize the canonical kmer. */
        KmerModel::Kmer tmp;  tmp.set (revcomp(*kmer, _kmerSize));

		*kmer = model->codeSeedRight (tmp, binrev[nt], Data::INTEGER).value();
		*kmer = revcomp(*kmer, _kmerSize);
	}
}

void AbstractDnaCoder::codeSeedNT(KmerModel* model, kmer_type* kmer, char nt, bool right){
	//if(nt == 'N') nt = 'A';
	return codeSeedBin(model, kmer, Leon::nt2bin(nt), right);
}


void AbstractDnaCoder::addErrorPos(int pos, bool rightExtend){

		_leftErrorPos.push_back(pos);
}

//====================================================================================
// ** DnaEncoder
//====================================================================================
DnaEncoder::DnaEncoder(Leon* leon) :
AbstractDnaCoder(leon), _itKmer(_kmerModel), _totalDnaSize(0), _readCount(0), _MCtotal(0), _readWithoutAnchorCount(0),
_MCuniqSolid (0), _MCuniqNoSolid(0), _MCnoAternative(0), _MCmultipleSolid(0)//, _MCmultipleNoSolid(0)
{
	_maxSequenceSize = 0;
	_minSequenceSize = INT_MAX;
	
	_thread_id = __sync_fetch_and_add (&_leon->_nb_thread_living, 1);

#ifdef PRINT_DISTRIB
	_distrib.resize(maxSequences);
	_outDistrib = 0;
#endif

	//pour quals
	if(! leon->_isFasta)
	{
	_max_read_size = 10000;
	_nb_solids = (int *) malloc(_max_read_size * sizeof(int) );
	_qualseq = (char *) malloc(_max_read_size*sizeof(char ));
	_bufferQuals_size = _leon->getReadPerBlock()* 200;
	_bufferQuals = (char *) malloc(_bufferQuals_size * sizeof(char ));
	_bufferQuals_idx=0;
	
	_trunc_mode = true;
	_smoothing_threshold = 2;
		
	}
}

DnaEncoder::DnaEncoder(const DnaEncoder& copy) :
AbstractDnaCoder(copy._leon), _itKmer(_kmerModel),
 _totalDnaSize(0), _readCount(0), _MCtotal(0), _readWithoutAnchorCount(0),
_MCuniqSolid (0), _MCuniqNoSolid(0), _MCnoAternative(0), _MCmultipleSolid(0)//, _MCmultipleNoSolid(0)
{
	_maxSequenceSize = 0;
	_minSequenceSize = INT_MAX;
	
#ifdef PRINT_DISTRIB
	_distrib.resize(maxSequences);
	_outDistrib = 0;
#endif

	_thread_id = __sync_fetch_and_add (&_leon->_nb_thread_living, 1);

	startBlock();

	
	//for quals
	if(! _leon->_isFasta)
	{
	_max_read_size = 10000;
	_nb_solids = (int *) malloc(_max_read_size * sizeof(int) );
	_qualseq = (char *) malloc(_max_read_size*sizeof(char ));
	_bufferQuals_size = _leon->getReadPerBlock()* 200;
	_bufferQuals = (char *) malloc(_bufferQuals_size * sizeof(char ));
	//printf("initial buffer qual size %i \n",_bufferQuals_size );

	_bufferQuals_idx =0;
	
	_trunc_mode = true;
	_smoothing_threshold = 2;
	}
	
	///
	
	#ifdef LEON_PRINT_STAT
		_rangeEncoder1.updateModel = false;
		_rangeEncoder2.updateModel = false;
		_rangeEncoder3.updateModel = false;
		_rangeEncoder4.updateModel = false;
		_rangeEncoder5.updateModel = false;
		_rangeEncoder6.updateModel = false;
	#endif

}

DnaEncoder::~DnaEncoder(){

	if(_thread_id!=0 && (_seqId+1) % _leon->getReadPerBlock() != 0 ){
		writeBlock();
	}
	//int nb_remaining =
	__sync_fetch_and_add (&_leon->_nb_thread_living, -1);

	//printf("\~ this decoder %lli seq  %lli  mctotal   %lli mltnos %p   tid %i \n",_readCount,_MCtotal,_MCmultipleNoSolid,this,_thread_id);
	__sync_fetch_and_add(&_leon->_readCount, _readCount);
	__sync_fetch_and_add(&_leon->_MCtotal, _MCtotal);
	__sync_fetch_and_add(&_leon->_readWithoutAnchorCount, _readWithoutAnchorCount);
	__sync_fetch_and_add(&_leon->_totalDnaSize, _totalDnaSize);
	__sync_fetch_and_add(&_leon->_MCuniqSolid, _MCuniqSolid);
	__sync_fetch_and_add(&_leon->_MCuniqNoSolid, _MCuniqNoSolid);
	__sync_fetch_and_add(&_leon->_MCnoAternative, _MCnoAternative);
	__sync_fetch_and_add(&_leon->_MCmultipleSolid, _MCmultipleSolid);
	
	_leon->updateMinMaxSequenceSize(_minSequenceSize,_maxSequenceSize);
	
	
	//__sync_fetch_and_add(&_leon->_MCmultipleNoSolid, _MCmultipleNoSolid);
	
	#ifdef LEON_PRINT_STAT
		__sync_fetch_and_add(&_leon->_anchorAdressSize, _rangeEncoder3.getBufferSize());
		__sync_fetch_and_add(&_leon->_anchorPosSize, _rangeEncoder2.getBufferSize());
		__sync_fetch_and_add(&_leon->_readSizeSize, _rangeEncoder1.getBufferSize());
		__sync_fetch_and_add(&_leon->_bifurcationSize, _rangeEncoder4.getBufferSize());
		__sync_fetch_and_add(&_leon->_otherSize, _rangeEncoder6.getBufferSize());
		__sync_fetch_and_add(&_leon->_noAnchorSize, _rangeEncoder5.getBufferSize());

		
		_rangeEncoder1.clear();
		_rangeEncoder2.clear();
		_rangeEncoder3.clear();
		_rangeEncoder4.clear();
		_rangeEncoder5.clear();
		_rangeEncoder6.clear();
	#endif
	
	
	//pour quals
	if(! _leon->_isFasta)
	{
		free(_nb_solids);
		free(_qualseq);
		free(_bufferQuals);
	}
}

void DnaEncoder::operator()(Sequence& sequence){

#ifdef PRINT_DISTRIB
	if(_sequences.size() > maxSequences){
		_sequences.pop_back();
	}
#endif

	_sequence = &sequence;
	//cout << _sequence->getIndex() << endl;
	_seqId = _sequence->getIndex() ;
	_readSize = _sequence->getDataSize();
	_readseq = _sequence->getDataBuffer();
		
	_totalDnaSize += _readSize ;
	
	
	_minSequenceSize = std::min(_minSequenceSize, (int) _readSize);
	_maxSequenceSize = std::max(_maxSequenceSize, (int)_readSize);
	
	
	//_lastSequenceIndex = sequence->getIndex();
	
//	if(_sequence->getIndex() % Leon::READ_PER_BLOCK == 0){

	execute();

	//_prevSequences = _sequence;
#ifdef PRINT_DISTRIB
	_sequences.insert(_sequences.begin(), _sequence);
#endif

	if(_processedSequenceCount >= _leon->getReadPerBlock() ){
		
		writeBlock();
		startBlock();
	}

}

void DnaEncoder::writeBlock(){
	if(_processedSequenceCount == 0) return;
	
	if(_rangeEncoder.getBufferSize() > 0){
		_rangeEncoder.flush();
	}
	
	int blockId = (  _seqId / _leon->getReadPerBlock())   ;
	//printf("\nTid %i  WB :  blockid %i sid %llu     size: %llu  _processedSequenceCount %i\n",_thread_id, blockId, _seqId, _rangeEncoder.getBufferSize(),_processedSequenceCount );

	//_leon->_realDnaCompressedSize += _rangeEncoder.getBufferSize();
	_leon->writeBlock(_rangeEncoder.getBuffer(), _rangeEncoder.getBufferSize(), _processedSequenceCount,blockId,false);
	_rangeEncoder.clear();
	
	if(! _leon->_isFasta)
	{
		_leon->writeBlockLena((u_int8_t*) _bufferQuals, _bufferQuals_idx ,_processedSequenceCount, blockId);
		_bufferQuals_idx = 0;
	}
	
#ifdef PRINT_DISTRIB
	cout << "----------------------------------------------------" << endl;
	for(int i=0; i<_distrib.size(); i++){
		cout << i << "    " << _distrib[i] << endl;
	}
	cout << "Adressed:    " << _outDistrib << endl;
#endif

	
}

void DnaEncoder::execute(){


	#ifdef PRINT_DEBUG_ENCODER
		cout << endl << "\tEncoding seq " << _sequence->getIndex() << endl;
		cout << "\t\t" << _readseq << endl;
	#endif
	
	//cout << _readseq << endl;
	
	_readCount +=1;
	_Npos.clear();
	
	if(_readSize < _kmerSize){
		encodeNoAnchorRead();
		if(! _leon->_isFasta)
		{
			smoothQuals();
		}
		endRead();
		return;
	}


	 
	//cout << _leon->_readCount << endl;
	//kmer_type anchorKmer = 0;
	u_int32_t anchorAddress;
	
	buildKmers(); // en profiter ici pour faire la compression des qual ?
	
	if(! _leon->_isFasta)
	{
		storeSolidCoverageInfo();
		smoothQuals();
	}
	

	//_isPrevReadAnchorable = false;
	int anchorPos = findExistingAnchor(&anchorAddress); //unsynch

	if(anchorPos == -1)
		anchorPos = _leon->findAndInsertAnchor(_kmers, &anchorAddress);  //unsynch

	//cout << anchorPos << endl;

	if(anchorPos == -1)
		encodeNoAnchorRead();
	else{
		encodeAnchorRead(anchorPos, anchorAddress);
	}
	//}
	
	endRead();

}




double DnaEncoder::char2proba(char c)
{
	int phred = c -33;
	
	double proba =  exp(-phred* log(10)/10);
	return proba;
	//Q 10 : 10% err
	//Q 20 : 1% err
	//Q 30  0.1% err ..
}


char DnaEncoder::char2phred(char c)
{
	return c -33;
}


void DnaEncoder::smoothQuals()
{
	strcpy (_qualseq, _sequence->getQuality().c_str());  // copy the qual sequence of this read in _qualseq

	if(! _leon->_lossless && _readSize >= _kmerSize)
	{
		for (unsigned int ii=0; ii< _readSize; ii++)
		{
			if ((_nb_solids[ii]>= _smoothing_threshold) || (((int) _qualseq[ii] > (int) '@') && _trunc_mode ))
			{
				apply_smoothing_at_pos (ii);
			}
		}
	}
	
	_qualseq[_readSize]='\n';
	_qualseq[_readSize+1]='\0';
	
	if( (_bufferQuals_idx+ _readSize+1 ) >= _bufferQuals_size)
	{
		//printf("_bufferQuals_size %i  _bufferQuals_idx %i  seqid %zu \n",_bufferQuals_size,_bufferQuals_idx,_sequence->getIndex()  );
		_bufferQuals_size = _bufferQuals_size * 2;
		_bufferQuals = (char *) realloc(_bufferQuals,_bufferQuals_size * sizeof(char) );
	}
	
	strcpy(_bufferQuals + _bufferQuals_idx, _qualseq);

	_bufferQuals_idx += _readSize+1 ; // with last '\n'

	//fprintf(_leon->_testQual,"%s",_qualseq); //st_qualseq.c_str()

}



bool DnaEncoder::apply_smoothing_at_pos(int pos)
{
	if(char2phred(_qualseq[pos])==0 || char2phred(_qualseq[pos])==2 )
		return false;
	
	bool ok_to_smooth= true;
	
	int diff = ('@' -  _qualseq[pos]);
	if(  diff > 10   )
	{
		if(_nb_solids[pos]>(diff-5))
			ok_to_smooth =true;
		else
			ok_to_smooth = false;
	}
	
	if(ok_to_smooth)
	{
		_qualseq[pos] = '@'; //smooth qual
		return true;
	}
	else return false;
	
}



void DnaEncoder::storeSolidCoverageInfo()
{
	kmer_type kmer, kmerMin;
	
	if(_readSize >= _max_read_size)
	{
		_max_read_size = _readSize + 1000;
		_nb_solids = (int *) realloc(_nb_solids,_max_read_size * sizeof(int) );
		_qualseq = (char *) realloc(_qualseq,_max_read_size*sizeof(char ));

	}
	memset(_nb_solids,0,_max_read_size * sizeof(int) );

	for(unsigned int ii=0; ii<_kmers.size(); ii++){
		kmer = _kmers[ii];
		kmerMin = min(kmer, revcomp(kmer, _kmerSize));

		if(_bloom->contains(kmerMin))
		{
			//increments all pos covered by the solid kmer
			for (unsigned int jj=0; jj< _kmerSize ; jj++)
			{
				_nb_solids[ii+jj] ++ ;
			}
		}
	}
}


void DnaEncoder::buildKmers(){

	
	
	for(unsigned int i=0; i<_readSize; i++){
		if(_readseq[i] == 'N'){
			_Npos.push_back(i);
			_readseq[i] = 'A';
		}
	}
	_itKmer.setData(_sequence->getData());
	
	_kmers.clear();
	for (_itKmer.first(); !_itKmer.isDone(); _itKmer.next()){
		//cout << (*_itKmer).toString(_kmerSize) << endl;
		_kmers.push_back(_itKmer->value());
	}
	

#ifdef PRINT_DISTRIB
	unordered_set<u_int64_t> H;
	for(kmer_type kmer : _kmers){
		kmer_type kmerMin = min(kmer, revcomp(kmer, _kmerSize));
		H.insert(kmerMin.getVal());
	}

	for(int i=0; i<_sequences.size(); i++){
		Sequence* sequence = _sequences[i];

		_itKmer.setData(sequence->getData());
		for (_itKmer.first(); !_itKmer.isDone(); _itKmer.next()){
			kmer_type kmerMin2 = min(_itKmer->value(), revcomp(_itKmer->value(), _kmerSize));
			if(H.find(kmerMin2.getVal()) != H.end()){
				_distrib[i] += 1;
				return;
			}
		}
	}

	_outDistrib += 1;
#endif

}

int DnaEncoder::findExistingAnchor(u_int32_t* anchorAddress){

	kmer_type kmer, kmerMin;

	
	for(unsigned int i=0; i<_kmers.size(); i++){
		kmer = _kmers[i];
		kmerMin = min(kmer, revcomp(kmer, _kmerSize));
		if(_leon->anchorExist(kmerMin, anchorAddress)){
			return i;
		}
	}
	return -1;
}

bool DnaEncoder::isReadAnchorable(){
	int nbKmerSolid = 0;
	kmer_type kmer, kmerMin;

	for(unsigned int i=0; i<_kmers.size(); i++){

		kmer = _kmers[i];
		kmerMin = min(kmer, revcomp(kmer, _kmerSize));

		if(_bloom->contains(kmerMin)){
			nbKmerSolid += 1;
			i += _kmerSize;
		}

		if(nbKmerSolid >= 2) return true;
	}

	return nbKmerSolid >= 2;

}

void DnaEncoder::encodeAnchorRead(int anchorPos, u_int32_t anchorAddress){
	#ifdef PRINT_DEBUG_ENCODER
		cout << "\t\tEncode anchor read" << endl;
	#endif
	//printf("encode  anchor read \n");

	//encode read type (0: read with anchor, 1: read without anchor)
	#ifdef LEON_PRINT_STAT
		_rangeEncoder6.encode(_readTypeModel, 0);
	#endif
	_rangeEncoder.encode(_readTypeModel, 0);
	
	//u_int64_t deltaValue;
	//u_int8_t deltaType;
	
	//Encode read size
	//deltaType = CompressionUtils::getDeltaValue(_readSize, _prevReadSize, &deltaValue);
	#ifdef LEON_PRINT_STAT
		//_rangeEncoder1.encode(_readSizeDeltaTypeModel, deltaType);
		CompressionUtils::encodeNumeric(_rangeEncoder1, _readSizeValueModel, _readSize);
	#endif
	//_rangeEncoder.encode(_readSizeDeltaTypeModel, deltaType);
	CompressionUtils::encodeNumeric(_rangeEncoder, _readSizeValueModel, _readSize);
	//_prevReadSize = _readSize;
	//printf("read size %i  deltaValue %i\n",_readSize,deltaValue);

	//Encode anchor pos
	//deltaType = CompressionUtils::getDeltaValue(anchorPos, _prevAnchorPos, &deltaValue);
	#ifdef LEON_PRINT_STAT
		//_rangeEncoder2.encode(_anchorPosDeltaTypeModel, deltaType);
		CompressionUtils::encodeNumeric(_rangeEncoder2, _anchorPosModel, anchorPos);
	#endif
	//_rangeEncoder.encode(_anchorPosDeltaTypeModel, deltaType);
	CompressionUtils::encodeNumeric(_rangeEncoder, _anchorPosModel, anchorPos);
	//_prevAnchorPos = anchorPos;
	//printf("anchor pos %i \n",anchorPos);

	//Encode anchor address
	//deltaType = CompressionUtils::getDeltaValue(anchorAddress, _prevAnchorAddress, &deltaValue);
	#ifdef LEON_PRINT_STAT
		//_rangeEncoder3.encode(_anchorAddressDeltaTypeModel, deltaType);
		CompressionUtils::encodeNumeric(_rangeEncoder3, _anchorAddressModel, anchorAddress);
	#endif
	//_rangeEncoder.encode(_anchorAddressDeltaTypeModel, deltaType);
	//if(_isPrevReadAnchorable){
		//	_rangeEncoder.encode(_isPrevReadAnchorableModel, 0);
		//	CompressionUtils::encodeNumeric(_rangeEncoder, _isPrevReadAnchorablePosModel, _isPrevReadAnchorablePos);
		//}
		//else{
		//_rangeEncoder.encode(_isPrevReadAnchorableModel, 1);
	CompressionUtils::encodeNumeric(_rangeEncoder, _anchorAddressModel, anchorAddress);
		//CompressionUtils::encodeNumeric(_rangeEncoder, _isPrevReadAnchorablePosModel, _isPrevReadAnchorablePos);
		//}
	//_prevAnchorAddress = anchorAddress;
	//printf("anchor adress %i \n",anchorAddress);

	
	kmer_type anchor = _kmers[anchorPos];
	
	//Encode a bit that says if the anchor is normal or revcomp
	if(anchor == min(anchor, revcomp(anchor, _kmerSize))){
		#ifdef LEON_PRINT_STAT
			_rangeEncoder6.encode(_readAnchorRevcompModel, 0);
		#endif
		_rangeEncoder.encode(_readAnchorRevcompModel, 0);
	}
	else{
		#ifdef LEON_PRINT_STAT
			_rangeEncoder6.encode(_readAnchorRevcompModel, 1);
		#endif
		_rangeEncoder.encode(_readAnchorRevcompModel, 1);
	}

	#ifdef PRINT_DEBUG_ENCODER
		cout << "\t\t\tAnchor pos: " << anchorPos << endl;
		cout << "\t\t\tAnchor: " << _kmers[anchorPos].toString(_kmerSize) << endl;
	#endif

	_bifurcations.clear();
	_binaryBifurcations.clear();
	_bifurcationTypes.clear();
	_leftErrorPos.clear();
	//_rightErrorPos.clear();

	
	kmer_type kmer = anchor;
	for(int i=anchorPos-1; i>=0; i--){
		kmer = buildBifurcationList(i, kmer, false);
		//i = buildBifurcationList(i, false);
		//cout << kmer.toString(_kmerSize) << endl;
	}

	kmer = anchor;
	for(unsigned int i=anchorPos+_kmerSize; i<_readSize; i++){
		//cout << "Pos: " << i << endl;
		kmer = buildBifurcationList(i, kmer, true);
		//i = buildBifurcationList(i, true);
	//for(int i=anchorPos; i<_kmers.size()-1; i++)
		//cout << kmer.toString(_kmerSize) << endl;
	}
		

	//Encode N positions
	_prevNpos = 0;
	#ifdef LEON_PRINT_STAT
		CompressionUtils::encodeNumeric(_rangeEncoder6, _numericModel, _Npos.size());
	#endif
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _Npos.size());
	for(unsigned int i=0; i<_Npos.size(); i++){
		//deltaType = CompressionUtils::getDeltaValue(_Npos[i], _prevNpos, &deltaValue);
		//_rangeEncoder.encode(_NposDeltaTypeModel, deltaType);
		#ifdef LEON_PRINT_STAT
			CompressionUtils::encodeNumeric(_rangeEncoder6, _NposModel, _Npos[i]-_prevNpos);
		#endif
		CompressionUtils::encodeNumeric(_rangeEncoder, _NposModel, _Npos[i]-_prevNpos);
		_prevNpos = _Npos[i];
	}
	
	#ifdef LEON_PRINT_STAT
		CompressionUtils::encodeNumeric(_rangeEncoder6, _leftErrorModel, _leftErrorPos.size());
	#endif
	CompressionUtils::encodeNumeric(_rangeEncoder, _leftErrorModel, _leftErrorPos.size());
	sort(_leftErrorPos.begin(), _leftErrorPos.end());
	_prevErrorPos = 0;
	for(unsigned int i=0; i<_leftErrorPos.size(); i++){
		#ifdef LEON_PRINT_STAT
			CompressionUtils::encodeNumeric(_rangeEncoder6, _leftErrorPosModel, _leftErrorPos[i]-_prevErrorPos);
		#endif
		CompressionUtils::encodeNumeric(_rangeEncoder, _leftErrorPosModel, _leftErrorPos[i]-_prevErrorPos);
		_prevErrorPos = _leftErrorPos[i];
	}

	
	u_int64_t bifType0 = 0;
	u_int64_t bifType1 = 0;
	//cout << _bifurcationTypes.size() << " " << _bifurcations.size() << " " << _binaryBifurcations.size() << endl;
	for(unsigned int i=0; i<_bifurcationTypes.size(); i++){
		u_int8_t type = _bifurcationTypes[i];
		if(type == 0){
			#ifdef LEON_PRINT_STAT
				_rangeEncoder4.encode(_bifurcationModel, _bifurcations[bifType0]);
			#endif
			//cout << Leon::nt2bin(_bifurcations[i]) << " ";
			_rangeEncoder.encode(_bifurcationModel, _bifurcations[bifType0]);
			bifType0 += 1;
		}
		else{
			#ifdef LEON_PRINT_STAT
				_rangeEncoder4.encode(_bifurcationBinaryModel, _binaryBifurcations[bifType1]);
			#endif
			//cout << Leon::nt2bin(_bifurcations[i]) << " ";
			_rangeEncoder.encode(_bifurcationBinaryModel, _binaryBifurcations[bifType1]);
			bifType1 += 1;
		}
	}

	
	
}
	
kmer_type DnaEncoder::buildBifurcationList(int pos, kmer_type kmer, bool rightExtend){
		
	char nextNt = _readseq[pos];
	int nextNtBin = Leon::nt2bin(nextNt);

	if(std::find(_Npos.begin(), _Npos.end(), pos) != _Npos.end()){
		codeSeedNT(&_kmerModel, &kmer, nextNt, rightExtend);
		return kmer;
		//return pos;
	}
	
	//kmer_type kmerMin;
	kmer_type uniqKmer;
	bool firstSolidKmer = false;
	//int uniqNt;
	//u_int8_t binNt2;
	bool isKmerSolid = false;
	
	int indexedKmerCount = 0;
	
	
	
	
	std::bitset<4> res4  = _bloom->contains4(kmer,rightExtend);
	for(int nt=0; nt<4; nt++){
		
		//mutatedKmer.printASCII(_kmerSize);
		
		if(res4[nt]){

			indexedKmerCount += 1;

			if(!firstSolidKmer){
				firstSolidKmer = true;
				uniqKmer = kmer;
				codeSeedBin(&_kmerModel, &uniqKmer, nt, rightExtend);
			}
			/*
			
			//uniqNt = nt;
			uniqKmer = mutatedKmer;
			*/
			
			if(nt == nextNtBin){
				isKmerSolid = true;
			}
		}
		
	}
	
	
	_MCtotal +=1;
	

	if(isKmerSolid){

		if(indexedKmerCount == 1){
			_MCuniqSolid += 1;
			return uniqKmer;
		}
		else if(indexedKmerCount == 2){

			char nt1 = -1;
			char nt2 = -1;

			for(int nt=0; nt<4; nt++){
				if(res4[nt]){
					//cout << "\t" << nt << endl;
					if(nt1 == -1)
						nt1 = nt;
					else if(nt2 == -1)
						nt2 = nt;
					else break;
				}
			}


			if(nt1 == nextNtBin){
				//cout << "\t0" << endl;
				_binaryBifurcations.push_back(0);
				_bifurcationTypes.push_back(1);
				_MCmultipleSolid += 1;
			}
			else if(nt2 == nextNtBin){
				//cout << "\t1" << endl;
				_binaryBifurcations.push_back(1);
				_bifurcationTypes.push_back(1);
				_MCmultipleSolid += 1;
			}
			else{

				//if(_sequence->getIndex() < 20)
				//	cout << "\tallo" << endl;
				//_MCuniqNoSolid += 1;
				//nextNt = Leon::bin2nt(nt1);
				//_bifurcations.push_back(nextNtBin);
				//_errorPos.push_back(pos);
				//_bifurcationTypes.push_back(0);
				//return uniqKmer;

				/*
				_MCuniqNoSolid += 1;
				//nextNt = Leon::bin2nt(nt1);

				_bifurcationTypes.push_back(0);
				_bifurcations.push_back(nextNtBin);
				_errorPos.push_back(pos);

				nextNtBin = getBestPath(pos, kmer, res4, rightExtend);
				//cout << (int)nextNtBin << endl;
				if(nextNtBin == -1){
					nextNtBin = nt1;
				}
				nextNt = Leon::bin2nt(nextNtBin);
				_bifurcations.push_back(nextNtBin);
				_bifurcationTypes.push_back(0);*/

			}
			//cout << "PROBLEME IN BUILD BINARY BIFURCATION (DnaEncoder - buildBifurcationList)" << endl;

			//if(_sequence->getIndex() < 10)
			//	cout << (char) Leon::bin2nt(nextNtBin) << endl;;

			codeSeedNT(&_kmerModel, &kmer, nextNt, rightExtend);
			return kmer;



		}
		else{

			_bifurcations.push_back(nextNtBin);
			_bifurcationTypes.push_back(0);
			codeSeedNT(&_kmerModel, &kmer, nextNt, rightExtend);
			return kmer;
		}

	}
	else{


		if(indexedKmerCount == 0){
			//cout << "PAF_BREAK    " << pos << endl;
			_MCnoAternative += 1;
		}
		else if(indexedKmerCount == 1){
			//cout << "PAF_UNIQ    " << pos << endl;
			_MCuniqNoSolid += 1;

			//_leon->_readWithAnchorMutationChoicesSize += 0.25;

			_bifurcations.push_back(nextNtBin);
			_bifurcationTypes.push_back(0);
			addErrorPos(pos, rightExtend);
			//_errorPos.push_back(pos);
			return uniqKmer;
		}
		else if(indexedKmerCount == 2){
			//cout << "PAF_MULTIPLE    " << pos << endl;
			_MCuniqNoSolid += 1;

		
			//encode error
			_bifurcations.push_back(nextNtBin);
			_bifurcationTypes.push_back(0);
			addErrorPos(pos, rightExtend);

			//get the first path in bufurcation
			for(int nt=0; nt<4; nt++){
				if(res4[nt]){
					nextNtBin = nt;
					nextNt = Leon::bin2nt(nt);
					break;
				}
			}

			codeSeedNT(&_kmerModel, &kmer, nextNt, rightExtend);
			return kmer;
			//return uniqKmer;
			//_bifurcationTypes.push_back(0);

			/*_bifurcations.push_back(nextNtBin);
			_bifurcationTypes.push_back(0);

			return uniqKmer;*/
		}
		else{
			_MCuniqNoSolid += 1;
			//cout << "PAF_MULTIPLE    " << pos << endl;
			/*
			_errorPos.push_back(pos);
			_bifurcations.push_back(nextNtBin);
			_bifurcationTypes.push_back(0);
			_errorPos.push_back(pos);
			return uniqKmer;*/

		}

		

		//_leon->_readWithAnchorMutationChoicesSize += 0.25;
		_bifurcations.push_back(nextNtBin);
		_bifurcationTypes.push_back(0);
		codeSeedNT(&_kmerModel, &kmer, nextNt, rightExtend);
		return kmer;

	}


	
}


int DnaEncoder::getBestPath(int pos, kmer_type& kmer, bitset<4>& initRes4, bool rightExtend){


	char ntInRead = 0;
	if(rightExtend){
		if((unsigned)(pos+1) < _readSize){
			ntInRead = _readseq[pos+1];
		}
	}
	else{
		if(pos-1 >= 0){
			ntInRead = _readseq[pos-1];
		}
	}

	int ntInReadBin = Leon::nt2bin(ntInRead);

	//int depth = 2;
	int bestNt = -1;
	bool isValid[4];
	for(int i=0; i<4; i++){
		isValid[i] = true;
	}

	//for(int j=0; j<depth; j++){

	for(int nt=0; nt<4; nt++){

		if(initRes4[nt]){

			kmer_type mutatedKmer = kmer;
			codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);

			bitset<4> res4  = _bloom->contains4(mutatedKmer, rightExtend);
			int nbSolidKmer = 0;

			for(int nt2=0; nt2<4; nt2++){

				if(res4[nt2]){
					nbSolidKmer += 1;

					if(nt2 == ntInReadBin){
						bestNt = nt;
					}
				}

			}

			if(nbSolidKmer != 1){
				isValid[nt] = false;
			}
		}
		else{
			isValid[nt] = false;
		}
	}

	//}

	int nbAlternative = 0;
	int uniqAlternative = -1;

	for(int nt=0; nt<4; nt++){
		if(isValid[nt]){
			nbAlternative += 1;
			uniqAlternative = nt;
		}
	}

	if(nbAlternative == 1){
		return uniqAlternative;
	}

	return bestNt;
}

int DnaEncoder::voteMutations(int pos, int depth, bool rightExtend){
	//kmer_type kmer;
	//int maxScore = 0;
	//int votes[4];
	
	int bestNt = 0;
	//bool isValid[4];
	
	//kmer_type mutatedKmers[4];
	//vector<int> mutations;
	//bool isMutateChoice[4];
	//kmer = _kmers[pos];
	/*
	for(int nt=0; nt<4; nt++){
		
		kmer_type mutatedKmer = kmer;
		codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
		kmer_type mutatedKmerMin = min(mutatedKmer, revcomp(mutatedKmer, _kmerSize));
		
		//mutated_kmer.printASCII(_kmerSize);
		
		if(_bloom->contains(mutatedKmerMin)){
			mutations.push_back(nt);
			mutatedKmers[nt] = mutatedKmer;
			votes[nt] = 0;
			//isMutateChoice[nt] = true;
		}
	}
	*/
	//kmer_type currentKmer;
	//cout << _readseq << endl;
	//cout << pos << ": " << kmer.toString(_kmerSize) << endl;
	
	/*
	for(int nt=0; nt<4; nt++){

		//mutatedKmer.printASCII(_kmerSize);

		if(res4[nt]){

			kmer_type mutatedKmer = kmer;
			codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);

			indexedKmerCount += 1;
			uniqNt = nt;
			uniqKmer = mutatedKmer;


			if(Leon::bin2nt(nt) == nextNt){
				isKmerSolid = true;
			}
		}

	}


	//for(int nt=0; nt<4; nt++){
		
		kmer_type mutatedKmer = kmer;
		codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
		kmer_type mutatedKmerMin = min(mutatedKmer, revcomp(mutatedKmer, _kmerSize));
		
	for(int j=0; j<depth; j++){
		//char nextNt;
		//int kmerPos;
		if(rightExtend){
			if(pos+1+_kmerSize+j >= _readSize) break;
			nextNt = _readseq[pos+1+_kmerSize+j];
		}
		else{
			if(pos-2-j < 0) break;
			nextNt = _readseq[pos-2-j];
		}

		std::bitset<4> res4  = _bloom->contains4(kmer,rightExtend);

		for(int nt=0; nt<4; nt++){

			if(res4[nt]){

			}
		}

		//cout << j << ": " << nextNt << endl;
		codeSeedNT(&_kmerModel, &mutatedKmer, nextNt, rightExtend);
		kmer_type mutatedKmerMin = min(mutatedKmer, revcomp(mutatedKmer, _kmerSize));

		if(_bloom->contains(mutatedKmerMin)){
			votes[nt] += 1;
			if(votes[nt] > maxScore){
				maxScore = votes[nt];
				bestNt = nt;
			}
		}
	}
	//}
	
	
	if(maxScore == 0){
		//cout << "No best NT" << endl;
		bestNt = -1;
	}
	//else
		//cout << "Best nt: " << bin2NT[bestNt] << endl;
	*/
	return bestNt;
	
	
}

void DnaEncoder::encodeNoAnchorRead(){
	#ifdef PRINT_DEBUG_ENCODER
		cout << "\t\tEncode no anchor read" << endl;
	#endif

	//printf("encode no anchor read \n");
	//Reinsert N because they can be encoded by the coder
	for(unsigned int i=0; i<_Npos.size(); i++){
		_readseq[_Npos[i]] = 'N';
	}
	
	#ifdef LEON_PRINT_STAT
		_rangeEncoder6.encode(_readTypeModel, 1);
	#endif
	_rangeEncoder.encode(_readTypeModel, 1);
	
	//_leon->_readWithoutAnchorSize += _readSize*0.375;
	_readWithoutAnchorCount +=1;
	
	/*
	for(int i=0; i<_readSize; i++){
		if(_readseq[i] == 'N'){
			_leon->_noAnchor_with_N_kmer_count += 1;
			break;
		}
	}
	
	bool full_N = true;
	for(int i=0; i<_readSize; i++){
		if(_readseq[i] != 'N'){
			full_N = false;
			break;
		}
	}
	if(full_N){
		_leon->_noAnchor_full_N_kmer_count += 1;
	}*/
	
	#ifdef LEON_PRINT_STAT
		CompressionUtils::encodeNumeric(_rangeEncoder5, _noAnchorReadSizeValueModel, _readSize);
	#endif
	CompressionUtils::encodeNumeric(_rangeEncoder, _noAnchorReadSizeValueModel, _readSize);
	
			
	for(unsigned int i=0; i<_readSize; i++){
		
		#ifdef LEON_PRINT_STAT
			_rangeEncoder5.encode(_noAnchorReadModel, Leon::nt2bin(_readseq[i]));
		#endif
		_rangeEncoder.encode(_noAnchorReadModel, Leon::nt2bin(_readseq[i]));
		
	}
	
}




QualDecoder::QualDecoder(Leon* leon, const string& inputFilename,tools::storage::impl::Group *  group)
//QualDecoder::QualDecoder(Leon* leon, const string& inputFilename)
{
	_group = group;
	_inputStream =0;
	_finished = false;
	_leon = leon;
	_inbuffer = NULL;
}


QualDecoder::~QualDecoder(){

	free(_inbuffer);
	if(_inputStream !=0) delete _inputStream;

}



void QualDecoder::setup( int blockID){
	
	_processedSequenceCount = 0;
	
	if(_inputStream !=0) delete _inputStream;
	
	std::string datasetname = Stringify::format ("qual_%i",blockID);

	_inputStream = new tools::storage::impl::Storage::istream  (*_group, datasetname);

	
	auto _tempcollec = & _group->getCollection<math::NativeInt8> (datasetname);
	std::string dsize = _tempcollec->getProperty ("size");
	
	_blockSize =  std::stoi(dsize); // blockSize;
	
	_inbuffer = (char * ) realloc(_inbuffer, _blockSize* sizeof(char));
	
}



void QualDecoder::setup(u_int64_t blockStartPos, u_int64_t blockSize, int sequenceCount){
	
	_processedSequenceCount = 0;

	
	_inputStream->seekg(blockStartPos, _inputStream->beg);
	
	_blockStartPos = blockStartPos;
	_blockSize = blockSize;

	


	
	
	
	_inbuffer = (char * ) realloc(_inbuffer, blockSize* sizeof(char));
	
	_sequenceCount = sequenceCount;
}



void QualDecoder::execute(){
	

//	printf("execute qual decoder _blockStartPos %llu  _blockSize %llu \n",_blockStartPos,_blockSize);

	_inputStream->read(_inbuffer,_blockSize );
	
	if(!_inputStream->good()) printf("inputstream E bad \n");

	//_inputFile->read(_inbuffer,_blockSize );
	
	//printf("----Begin decomp of Block     ----\n");

	z_stream zs;
	memset(&zs, 0, sizeof(zs));
	
	//deflateinit2 to be able to gunzip it fro mterminal
	
	//if (inflateInit2(&zs, (15+32)) != Z_OK)
	if (inflateInit (&zs) != Z_OK)
		throw    Exception ("inflate Init failed while decompressing.");
	
	zs.next_in = (Bytef*) _inbuffer ;
	zs.avail_in = _blockSize ;    // set the z_stream's input
	
	int ret;
	char outbuffer[32768];
	
	// retrieve the compressed bytes blockwise
	do {
		zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
		zs.avail_out = sizeof(outbuffer);
		
		ret = inflate (&zs, Z_SYNC_FLUSH); //ou Z_FINISH ? Z_SYNC_FLUSH
		
		if (_buffer.size() < zs.total_out) {
			// append the block to the output string
			_buffer.append(outbuffer,
							 zs.total_out - _buffer.size());
		}
		if (ret != Z_OK)
		{
			//printf("ret val %i  _blockStartPos %llu \n",ret,_blockStartPos);
			break;
		}
		else
		{
			//printf("-----block ret ok  _blockStartPos %llu  ----\n",_blockStartPos);
		}
	} while (ret == Z_OK);
	
	inflateEnd(&zs);

	_finished = true;
	
	//printf("Should be done decompressing block, in size   %llu   out size    %lu  \n",_blockSize,_buffer.size() );

}


//====================================================================================
// ** DnaDecoder
//====================================================================================
DnaDecoder::DnaDecoder(Leon* leon, const string& inputFilename,tools::storage::impl::Group *  group) :
AbstractDnaCoder(leon)
{
	_group = group;
	_inputStream =0;
	
	//_inputFile = new ifstream(inputFilename.c_str(), ios::in|ios::binary);
	_finished = false;
	
	//_anchorDictFile = new ifstream(_leon->_anchorDictFilename.c_str(), ios::in);
	
}

DnaDecoder::~DnaDecoder(){
	//delete _rangeDecoder;
	//delete _outputFile;
//	delete _inputFile;
//	delete _anchorDictFile;
	
	if(_inputStream !=0) delete _inputStream;

}

void DnaDecoder::setup(u_int64_t blockStartPos, u_int64_t blockSize, int sequenceCount,int blockID){
	startBlock();
	_rangeDecoder.clear();
	
	//_inputFile->seekg(blockStartPos, _inputFile->beg);
	//_rangeDecoder.setInputFile(_inputFile);
	
	if(_inputStream !=0) delete _inputStream;
	std::string datasetname = Stringify::format ("dna_%i",blockID);
	
	_inputStream = new tools::storage::impl::Storage::istream  (*_group, datasetname);
	
	auto _tempcollec = & _group->getCollection<math::NativeInt8> (datasetname);
	std::string dsize = _tempcollec->getProperty ("size");
	
	_blockSize =  std::stoi(dsize); // blockSize;
	_rangeDecoder.setInputFile(_inputStream);

	
	
	_blockStartPos = blockStartPos;
	_blockSize = blockSize;
	
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t-----------------------" << endl;
		cout << "\tDecoding Dna block " << _blockStartPos << " - " << _blockStartPos+_blockSize << endl;
	#else
		;//_leon->_progress_decode->inc(1);
	#endif
	
	_sequenceCount = sequenceCount;
}

void DnaDecoder::execute(){
	
	//decodeFirstHeader();
		
	while(_processedSequenceCount < _sequenceCount){
		
		//cout << "lala" << endl;
	//int i=0;
	//while(i < Leon::READ_PER_BLOCK){
	//while(_inputFile->tellg() <= _blockStartPos+_blockSize){
		//if(_leon->_readCount > 1) return;
	
	
		u_int8_t readType = _rangeDecoder.nextByte(_readTypeModel);
		//cout << "Read type: " << (int)readType << endl;

		if(readType == 0)
			decodeAnchorRead(); //ici
		else if(readType == 1)
			decodeNoAnchorRead();
			
		endRead();
		//cout << _inputFile->tellg() << " " << _blockStartPos+_blockSize << endl;
		/*
		string trueSeq = string((*_leon->_testBankIt)->getDataBuffer());
		trueSeq = trueSeq.substr(0, _readSize);
		//cout << trueSeq << endl;
		//cout << _currentSeq << endl;
		if(trueSeq != _currentSeq){
			cout << (*_leon->_testBankIt)->getIndex() << "\t\tseq different !!" << endl;
			cout << "\t\t" << trueSeq << endl;
			cout << "\t\t" << _currentSeq << endl;
			_leon->_readCount += 1;
			return;
		}
		_leon->_testBankIt->next();
		*/
		#ifdef PRINT_DEBUG_DECODER
			_readCount += 1;
			cout << _leon->_readCount << ": " << _currentSeq << endl;
		#endif
		
		//i++;
		//_leon->_readCount += 1;
		//if(i == 1) return;
		//_currentSeq.clear();
		
		//cout << (int)(_inputFile->tellg() < _blockStartPos+_blockSize) << endl;
		
	}
	
	//cout << "endooo" << endl;
	_finished = true;
	
}

void DnaDecoder::decodeAnchorRead(){
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tDecode anchor read" << endl;
	#endif

	//u_int8_t deltaType;
	//u_int64_t deltaValue;
	
	//printf("Decode anchor read \n");

	//Decode read size
	//deltaType = _rangeDecoder.nextByte(_readSizeDeltaTypeModel);
	//deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _readSizeValueModel);
//	printf("read deltaValue %llu \n",deltaValue);
	//_readSize = CompressionUtils::getValueFromDelta(deltaType, _prevReadSize, deltaValue);
	//_prevReadSize = _readSize;
	_readSize = CompressionUtils::decodeNumeric(_rangeDecoder, _readSizeValueModel);
	
//	printf("read size %i \n",_readSize);
	
	//Decode anchor pos
	//deltaType = _rangeDecoder.nextByte(_anchorPosDeltaTypeModel);
	//deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosModel);
	//int anchorPos = CompressionUtils::getValueFromDelta(deltaType, _prevAnchorPos, deltaValue);
	//_prevAnchorPos = anchorPos;
//	printf("anchor pos %i \n",anchorPos);
	int anchorPos = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosModel);

	
	//Decode anchor address
	//deltaType = _rangeDecoder.nextByte(_anchorAddressDeltaTypeModel);
	//deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorAddressModel);
	//u_int64_t anchorAddress = CompressionUtils::getValueFromDelta(deltaType, _prevAnchorAddress, deltaValue);
	//_prevAnchorAddress = anchorAddress;
	u_int64_t anchorAddress = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorAddressModel);
	
	kmer_type anchor = _leon->getAnchor(_anchorDictFile, anchorAddress); //laa
	
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\t\tRead size: " << _readSize << endl;
		cout << "\t\t\tAnchor pos: " << anchorPos << endl;
		cout << "\t\t\tAnchor adress: " << anchorAddress << endl;
		cout << "\t\t\tAnchor: " << anchor.toString(_kmerSize) << endl;
	#endif
	
	//Decode the bit that says if the anchor is revcomp or not
	if(_rangeDecoder.nextByte(_readAnchorRevcompModel) == 1)
		anchor = revcomp(anchor, _kmerSize);
		
	_currentSeq = anchor.toString(_kmerSize);	
	_leftErrorPos.clear();
	//_rightErrorPos.clear();
	_Npos.clear();
	
	//cout << _readSize << " " << anchorPos << " " << anchorAddress << endl;
	//Decode N pos
	_prevNpos = 0;
	u_int64_t NposCount = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
	for(unsigned int i=0; i<NposCount; i++){
		//deltaType = _rangeDecoder.nextByte(_NposDeltaTypeModel);
		//deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _NposModel);
		//u_int64_t nPos = CompressionUtils::getValueFromDelta(deltaType, _prevNpos, deltaValue);
		u_int64_t nPos = CompressionUtils::decodeNumeric(_rangeDecoder, _NposModel) + _prevNpos;
		_Npos.push_back(nPos);
		//cout << nPos << endl;
		_prevNpos = nPos;
		//_Npos.push_back(CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosSizeModel, _anchorPosModel));
	}
	
	//Decode error pos
	u_int64_t nbLeftError = CompressionUtils::decodeNumeric(_rangeDecoder, _leftErrorModel);
	_prevErrorPos = 0;
	for(unsigned int i=0; i<nbLeftError; i++){
		u_int64_t errorPos = CompressionUtils::decodeNumeric(_rangeDecoder, _leftErrorPosModel) + _prevErrorPos;
		addErrorPos(errorPos, true);
		_prevErrorPos = errorPos;
	}

	/*
	u_int64_t nbRightError = CompressionUtils::decodeNumeric(_rangeDecoder, _rightErrorModel);
	_prevErrorPos = anchorPos-1;
	for(int i=0; i<nbLeftError; i++){
		int deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _leftErrorPosModel);
		int errorPos = _prevErrorPos-deltaValue;
		//deltaType = _rangeDecoder.nextByte(_errorPosDeltaTypeModel);
		//deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _errorPosModel); //reprise
		//u_int64_t errorPos = CompressionUtils::getValueFromDelta(deltaType, _prevErrorPos, deltaValue);
		addErrorPos(errorPos, false);
		//_errorPos.push_back(errorPos);
		_prevErrorPos = errorPos;
		//_errorPos.push_back(CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosSizeModel, _anchorPosModel));
	}
	_prevErrorPos = anchorPos+_kmerSize;
	for(int i=0; i<nbRightError; i++){
		int deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _rightErrorPosModel);
		int errorPos = _prevErrorPos+deltaValue;
		//deltaType = _rangeDecoder.nextByte(_errorPosDeltaTypeModel);
		//deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _errorPosModel); //reprise
		//u_int64_t errorPos = CompressionUtils::getValueFromDelta(deltaType, _prevErrorPos, deltaValue);
		addErrorPos(errorPos, true);
		//_errorPos.push_back(errorPos);
		_prevErrorPos = errorPos;
		//_errorPos.push_back(CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosSizeModel, _anchorPosModel));
	}*/

	
	//Extend anchor to the left
	kmer_type kmer = anchor;
	for(int i=anchorPos-1; i>=0; i--){
		kmer = extendAnchor(kmer, i, false);
	}
	
	//Extend anchor to the right
	kmer = anchor;
	for(unsigned int i=anchorPos+_kmerSize; i<_readSize; i++){
		kmer = extendAnchor(kmer, i, true);
		//cout << "\t" << kmer.toString(_kmerSize) << endl;
	}
	
	//Inject N in the decoded read sequence
	//printf("npos s %i currseq %s \n",_Npos.size(),_currentSeq.c_str());
	for(unsigned int i=0; i<_Npos.size(); i++){
		_currentSeq[_Npos[i]] = 'N';
	}
}

kmer_type DnaDecoder::extendAnchor(kmer_type kmer, int pos, bool rightExtend){
	
	u_int8_t nextNt;
	//int nextNtBin;
	kmer_type resultKmer;
		
	if(std::find(_Npos.begin(), _Npos.end(), pos) != _Npos.end()){
		nextNt = 'A';
		if(rightExtend){
			_currentSeq += nextNt;
		}
		else{
			_currentSeq.insert(_currentSeq.begin(), nextNt);
		}
		//cout << _currentSeq << endl;
		//if(nextNt == 'N') nextN
		
		resultKmer = kmer;
		codeSeedNT(&_kmerModel, &resultKmer, nextNt, rightExtend);
		//cout << kmer.toString(_kmerSize) << endl;
		//cout << resultKmer.toString(_kmerSize) << endl;
		return resultKmer;
	}
	
	/*
	if(rightExtend){
		if(std::find(_rightErrorPos.begin(), _rightErrorPos.end(), pos) != _rightErrorPos.end()){
			nextNt = Leon::bin2nt(_rangeDecoder.nextByte(_bifurcationModel));

			if(rightExtend)
				_currentSeq += nextNt;
			else
				_currentSeq.insert(_currentSeq.begin(), nextNt);

			std::bitset<4> res4  = _bloom->contains4(kmer,rightExtend);

			for(int nt=0; nt<4; nt++){

				if(res4[nt]){
					kmer_type mutatedKmer = kmer;
					codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
					return mutatedKmer;
				}
			}
		}
	}
	else{*/
		if(std::find(_leftErrorPos.begin(), _leftErrorPos.end(), pos) != _leftErrorPos.end()){
			nextNt = Leon::bin2nt(_rangeDecoder.nextByte(_bifurcationModel));

			if(rightExtend)
				_currentSeq += nextNt;
			else
				_currentSeq.insert(_currentSeq.begin(), nextNt);

			std::bitset<4> res4  = _bloom->contains4(kmer,rightExtend);

			for(int nt=0; nt<4; nt++){

				if(res4[nt]){
					kmer_type mutatedKmer = kmer;
					codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
					return mutatedKmer;
				}
			}
		}
	//}

	
		
	//cout << kmer.toString(_kmerSize) << endl;
	kmer_type uniqKmer;
	//, mutatedSolidKmer;
	int uniqNt;
	//bool isKmerSolid = false;
	

	//kmer = _kmers[pos];

	int indexedKmerCount = 0;
	
	//cout << kmer.toString(_kmerSize) << endl;
	
	
	
	std::bitset<4> res4  = _bloom->contains4(kmer,rightExtend);
	for(int nt=0; nt<4; nt++){
		if(res4[nt]){
			kmer_type mutatedKmer = kmer;
			codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
			
			indexedKmerCount += 1;
			uniqNt = nt;
			uniqKmer = mutatedKmer;
		}
	}
		

	
	/*
	for(int nt=0; nt<4; nt++){
		//if(nt == original_nt){
		//	continue;
		//}
		
		kmer_type mutatedKmer = kmer;
		codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
		kmer_type mutatedKmerMin = min(mutatedKmer, revcomp(mutatedKmer, _kmerSize));
		
		//mutatedKmer.printASCII(_kmerSize);
		
		if(_bloom->contains(mutatedKmerMin)){
			indexedKmerCount += 1;
			uniqNt = nt;
			uniqKmer = mutatedKmer;
		}
		
	}
 */
	
	if(indexedKmerCount == 1){
		nextNt = Leon::bin2nt(uniqNt);
		//cout << "case 1         " << nextNt << endl;
		resultKmer = uniqKmer;
	}
	else if(indexedKmerCount == 2){

		char nt1 = -1;
		char nt2 = -1;

		for(int nt=0; nt<4; nt++){
			if(res4[nt]){
				if(nt1 == -1)
					nt1 = nt;
				else if(nt2 == -1)
					nt2 = nt;
				else break;
			}
		}

		u_int8_t nextBinaryNt = _rangeDecoder.nextByte(_bifurcationBinaryModel);

		//cout << (int)nextBinaryNt << endl;
		if(nextBinaryNt == 0)
			nextNt = Leon::bin2nt(nt1);
		else
			nextNt = Leon::bin2nt(nt2);

		//cout << nextNt << endl;
		resultKmer = kmer;
		codeSeedNT(&_kmerModel, &resultKmer, nextNt, rightExtend);


	}
	else{
		nextNt = Leon::bin2nt(_rangeDecoder.nextByte(_bifurcationModel));
		//cout << "case 2          "<< nextNt << endl;
		resultKmer = kmer;
		codeSeedNT(&_kmerModel, &resultKmer, nextNt, rightExtend);
	}
	
	//cout << nextNt << endl;
	
	//if(nextNt == 'N') cout << "lala" << endl;
	
	if(rightExtend){
		_currentSeq += nextNt;
	}
	else{
		_currentSeq.insert(_currentSeq.begin(), nextNt);
	}
	
	//cout << _currentSeq << endl;
	//cout << resultKmer.toString(_kmerSize) << endl;
	return resultKmer;
		
}

void DnaDecoder::decodeNoAnchorRead(){
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tDecode no anchor read" << endl;
	#endif
	
	_readSize = CompressionUtils::decodeNumeric(_rangeDecoder, _noAnchorReadSizeValueModel);
	//cout << "\tRead size: " << _readSize << endl;
	for(unsigned int i=0; i<_readSize; i++){
		_currentSeq += Leon::bin2nt(_rangeDecoder.nextByte(_noAnchorReadModel));
	}
	//endSeq();
	//cout << read << endl;
}
	
void DnaDecoder::endRead(){
	AbstractDnaCoder::endRead();
	
	_buffer += _currentSeq + '\n';
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\t\tRead: " << _currentSeq << endl;
	#endif
	//_outputFile->write(_currentSeq.c_str(), _currentSeq.size());
	_currentSeq.clear();
}




