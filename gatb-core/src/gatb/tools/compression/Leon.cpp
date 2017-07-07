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


 
#include "Leon.hpp"


using namespace std;

//#define PRINT_DEBUG
//#define PRINT_DEBUG_DECODER



//const u_int64_t ANCHOR_KMERS_HASH_SIZE = 500000000;
const char* Leon::STR_COMPRESS = "-c";
const char* Leon::STR_DECOMPRESS = "-d";
const char* Leon::STR_TEST_DECOMPRESSED_FILE = "-test-file";
const char* Leon::STR_DNA_ONLY = "-seq-only";
const char* Leon::STR_NOHEADER = "-noheader";
const char* Leon::STR_NOQUAL = "-noqual";

const char* Leon::STR_DATA_INFO = "Info";
const char* Leon::STR_INIT_ITER = "-init-iterator";

const int Leon::nt2binTab[128] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 1, 0, 0, //69
	0, 3, 0, 0, 0, 0, 0, 0, 4, 0, //79
	0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	};
const int Leon::bin2ntTab[5] = {'A', 'C', 'T', 'G', 'N'};


template<typename T>
void createDataset(tools::storage::impl::Group *  group, std::string datasetname, T data)
{
	tools::storage::impl::Storage::ostream os (*group, datasetname);
	os.write (reinterpret_cast<char const*>(&data), sizeof(data));
	os.flush();
}



template<typename T>
void readDataset(tools::storage::impl::Group *  group, std::string datasetname, T & data)
{
	tools::storage::impl::Storage::istream is (*group, datasetname);
	is.read (reinterpret_cast<char*> (&data),sizeof(data));
}





//Leon::Leon ( bool compress, bool decompress) :
Leon::Leon () :
Tool("leon"),
_progress_decode(0),_generalModel(256),_inputBank(0),_anchorDictModel(5)
{
	_MCnoAternative = _MCuniqSolid = _MCuniqNoSolid = _totalDnaSize =  _compressedSize = _readCount = _MCtotal = _nb_thread_living = 0;
	_compressed_qualSize = _anchorDictSize = _MCmultipleSolid = _anchorAdressSize = _readWithoutAnchorCount = _anchorPosSize = 0;
	_input_qualSize = _total_nb_quals_smoothed = _otherSize =  _readSizeSize =  _bifurcationSize =  _noAnchorSize = 0;
	_lossless = false;
	_storageH5file = 0;
	_bloom = 0;
	
	_isFasta = true;
	_maxSequenceSize = 0;
	_minSequenceSize = INT_MAX;
	std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);

    //_kmerSize(27)
	//_compress = compress;
	//_decompress = decompress;

	/** We don't want default options of Tool (or we want to put them at a specific location). */
	setParser (new OptionsParser ("leon"));

    getParser()->push_back (new OptionOneParam (STR_URI_FILE, "input file (e.g. FASTA/FASTQ for compress or .leon file for decompress)",   true));
    getParser()->push_back (new OptionNoParam  ("-c", "compression",   false));
    getParser()->push_back (new OptionNoParam  ("-d", "decompression", false));
    getParser()->push_back (new OptionOneParam (STR_NB_CORES, "number of cores (default is the available number of cores)", false, "0"));
    getParser()->push_back (new OptionOneParam (STR_VERBOSE,  "verbosity level",                                            false, "1", false),0,false);

    getParser()->push_back (new OptionOneParam ("-reads",     "number of reads per block (default is 50000)",               false, "50000", false),0,false);
	
	getParser()->push_back (new OptionNoParam  ("-lossless", "switch to lossless compression for qualities (default is lossy. lossy has much higher compression rate, and the loss is in fact a gain. lossy is better!)",   false));

	
    IOptionsParser* compressionParser = new OptionsParser ("compression");

    /** We add the sorting count options and hide all of them by default and display one some of them. */
    compressionParser->push_back (SortingCountAlgorithm<>::getOptionsParser(), 1, false);

    if (IOptionsParser* input = compressionParser->getParser (STR_URI_INPUT))  {  input->setName (STR_URI_FILE);  }
    if (IOptionsParser* input = compressionParser->getParser (STR_KMER_SIZE))  {  input->setVisible (true);  }

    compressionParser->push_back (new OptionOneParam(STR_KMER_ABUNDANCE, "abundance threshold for solid kmers (default inferred)", false));
	compressionParser->push_back (new OptionNoParam  (STR_INIT_ITER, "init iterator for ibank mode", false));

	compressionParser->push_back (new OptionNoParam (Leon::STR_DNA_ONLY, "store dna seq only, header and quals are discarded, will decompress to fasta (same as -noheader -noqual)", false));

	compressionParser->push_back (new OptionNoParam (Leon::STR_NOHEADER, "discard header", false));
	compressionParser->push_back (new OptionNoParam (Leon::STR_NOQUAL, "discard quality scores", false));

    IOptionsParser* decompressionParser = new OptionsParser ("decompression");
    decompressionParser->push_back (new OptionNoParam (Leon::STR_TEST_DECOMPRESSED_FILE, "check if decompressed file is the same as original file (both files must be in the same folder)", false));
	

	_subgroupInfoCollection = NULL;
	 _groupLeon = _subgroupQual = _subgroupInfo = _subgroupDict = _subgroupDNA =  _subgroupHeader = NULL;

	
    getParser()->push_back (compressionParser);
    getParser()->push_back (decompressionParser, 0, false);

	pthread_mutex_init(&findAndInsert_mutex, NULL);
	pthread_mutex_init(&writeblock_mutex, NULL);
	pthread_mutex_init(&minmax_mutex, NULL);

	
}

Leon::~Leon ()
{
	setInputBank (0);
	
	if(_storageH5file !=0)
		delete _storageH5file;
	
	if (_progress_decode)  { delete _progress_decode; }
}

void Leon::execute()
{
	_time = clock(); //Used to calculate time taken by decompression
	
	gettimeofday(&_tim, NULL);
	 _wdebut_leon = _tim.tv_sec +(_tim.tv_usec/1000000.0);
	
	_iterator_mode=false;
	if(getParser()->saw (STR_INIT_ITER))
		_iterator_mode = true;
	
	if(getParser()->saw ("-lossless"))
		_lossless = true;
		
    _compress = false;
    _decompress = false;
    if(getParser()->saw (Leon::STR_COMPRESS)) _compress = true;
    if(getParser()->saw (Leon::STR_DECOMPRESS)) _decompress = true;
	if((_compress && _decompress) || (!_compress && !_decompress)){
		cout << "Choose one option among -c (compress) or -d (decompress)" << endl << endl;
		return;
	}

	//getParser()->displayWarnings(); //pb si ici, affiche warnings apres exec dsk ,et prob option -c -d qui sont pas dans le parser car 'globales'
	
   // u_int64_t total_nb_solid_kmers_in_reads = 0;
   // int nb_threads_living;
    _nb_cores = getInput()->getInt(STR_NB_CORES);
    
    setReadPerBlock(getInput()->getInt("-reads"));
    
	//setup global
	for(int i=0; i<CompressionUtils::NB_MODELS_PER_NUMERIC; i++){
		_numericModel.push_back(Order0Model(256));
	}
	
	if(_compress){
		//#define SERIAL
		executeCompression();
	}
	else{
		executeDecompression();
	}
	

        
    
    //outputFile->flush();

    /*************************************************/
    // We gather some statistics.
    /*************************************************/
    //getInfo()->add (1, "result");
    //getInfo()->add (2, "nb solid kmers in reads", "%ld", total_nb_solid_kmers_in_reads);
    
    if(_decompress){
	//	delete _inputFile;
		if(! _iterator_mode)
		{
			delete _outputFile;
		}
		//if(! _isFasta) delete _inputFileQual;
		delete _bloom;
	}
    
}

void Leon::createBloom (){
    TIME_INFO (getTimeInfo(), "fill bloom filter");
	
	//u_int64_t solidFileSize
	
	int _auto_cutoff = 0 ;
	u_int64_t nbs = 0 ;
	u_int64_t nb_kmers_infile;
	
	//cout << _dskOutputFilename << endl;
	Storage* storage = StorageFactory(STORAGE_HDF5).load (_dskOutputFilename);
	LOCAL (storage);
	
	Partition<kmer_count> & solidCollection = storage->root().getGroup("dsk").getPartition<kmer_count> ("solid");
	
	/** We get the number of solid kmers. */
   // u_int64_t solidFileSize = solidCollection.getNbItems();
	
	nb_kmers_infile = solidCollection.getNbItems();
	//(System::file().getSize(_dskOutputFilename) / sizeof (kmer_count)); //approx total number of kmer

	if( ! getParser()->saw(STR_KMER_ABUNDANCE)){
		
		//retrieve cutoff
		
		Collection<NativeInt64>& cutoff  = storage->getGroup("histogram").getCollection<NativeInt64> ("cutoff");
		Iterator<NativeInt64>* iter = cutoff.iterator();
		LOCAL (iter);
		for (iter->first(); !iter->isDone(); iter->next())  {
			_auto_cutoff = iter->item().toInt();
		}
		//////
		
		//retrieve nb solids
		
		Collection<NativeInt64>& storagesolid  = storage->getGroup("histogram").getCollection<NativeInt64> ("nbsolidsforcutoff");
		Iterator<NativeInt64>* iter2 = storagesolid.iterator();
		LOCAL (iter2);
		for (iter2->first(); !iter2->isDone(); iter2->next())  {
			nbs = iter2->item().toInt();
		}
		//////
		
		_nks = _auto_cutoff;

		//printf("\tcutoff auto: %i  \n",_nks);

	}
	else
	{
		_auto_cutoff =0;
		nbs  = nb_kmers_infile;
	//	printf("\tcutoff user: %i (total solids %lli) \n",_nks,nbs);
	}
	
	
	


	
	
    //double lg2 = log(2);
    //float NBITS_PER_KMER = log (16*_kmerSize*(lg2*lg2))/(lg2*lg2);
    int NBITS_PER_KMER = 12;
    
	
    u_int64_t estimatedBloomSize = (u_int64_t) ((double)nbs * NBITS_PER_KMER);
    if (estimatedBloomSize ==0 ) { estimatedBloomSize = 1000; }
    
    
    //printf("raw solidFileSize %llu fsize %llu    %lu %lu \n",System::file().getSize(_solidFile), solidFileSize,sizeof (kmer_type),sizeof (kmer_count));
    
    /** We create the kmers iterator from the solid file. */
//    Iterator<kmer_count>* itKmers = createIterator<kmer_count> (
//                                                                new IteratorFile<kmer_count>(_dskOutputFilename),
//                                                                nb_kmers_infile,
//                                                                "fill bloom filter"
//                                                                );
	
	/** We create the kmers iterator from the solid file. */
    Iterator<kmer_count>* itKmers = createIterator<kmer_count> (
																solidCollection.iterator(),
																nb_kmers_infile,
																"fill bloom filter"
																);
    LOCAL (itKmers);

	
	

	
    /** We instantiate the bloom object. */
    //BloomBuilder<> builder (estimatedBloomSize, 7,tools::collections::impl::BloomFactory::CACHE,getInput()->getInt(STR_NB_CORES));
    //cout << "ESTIMATED:" << estimatedBloomSize << endl;
    //_bloomSize = estimatedBloomSize;
	
    if(_auto_cutoff){
		    getInfo()->add (0, "Abundance threshold");
        getInfo()->add (1, "cut-off (auto)", "%d", _auto_cutoff);
        getInfo()->add (1, "nb solid kmers", "%d", nbs);
    }
    else{
        getInfo()->add (0, "Abundance threshold");
        getInfo()->add (1, "cut-off", "%d", _nks);
        getInfo()->add (1, "nb solid kmers", "%d", nbs);
    }
    	
	//modif ici pour virer les kmers < auto cutoff
    BloomBuilder<> builder (estimatedBloomSize, 7,_kmerSize,tools::misc::BLOOM_NEIGHBOR,getInput()->getInt(STR_NB_CORES),_auto_cutoff);
    _bloom = builder.build (itKmers); // BLOOM_NEIGHBOR // BLOOM_CACHE


}



void Leon::executeCompression(){
	
	
	
	
	#ifdef PRINT_DEBUG
		cout << "Start compression" << endl;
	#endif

    _kmerSize      = getInput()->getInt (STR_KMER_SIZE);
    _nks      = getInput()->get(STR_KMER_ABUNDANCE) ? getInput()->getInt(STR_KMER_ABUNDANCE) : 3;
	//_nks           = getInput()->getInt (STR_KMER_ABUNDANCE);
    _inputFilename = getInput()->getStr (STR_URI_FILE);
    
	#ifdef PRINT_DEBUG
		cout << "\tInput filename: " << _inputFilename << endl;
	#endif
	
	u_int8_t infoByte = 0;
	
	
	/** We look for the beginnin of the suffix. */
	int lastindex = _inputFilename.find_last_of (".");
	
	/** We build the result. */
	string extension = _inputFilename.substr(lastindex+1);
	
	_noHeader =false;

	
	if(getParser()->saw (Leon::STR_NOHEADER))
	{
		_noHeader = true;
		infoByte |= 0x02; //no header
	}
	
	if(getParser()->saw (Leon::STR_NOQUAL))
	{
		_isFasta = true;
		infoByte |= 0x01; //fasta mode == no quals
	}
	
	
	if(getParser()->saw (Leon::STR_DNA_ONLY))
	{
		_noHeader = true;
		_isFasta = true;

		infoByte |= 0x02; //no header
		infoByte |= 0x01; //fasta mode == no quals

	}


    //_inputBank = Bank::singleton().createBank(_inputFilename);
	setInputBank (Bank::open(_inputFilename));
	
	//cout << Bank::getType(_inputFilename) << endl;
	
	
      if(_inputFilename.find(".fq") !=  string::npos || _inputFilename.find(".fastq") !=  string::npos)
	{
		getInfo()->add (0, "Input format: FastQ");
		
		if(! getParser()->saw (Leon::STR_DNA_ONLY) && ! getParser()->saw (Leon::STR_NOQUAL))
		{
			
			if (_lossless)
          getInfo()->add (0, "Quality compression: LOSSLESS mode");
       else
          getInfo()->add (0, "Quality compression: lossy mode (use '-lossless' for lossless compression)");
			
			_isFasta = false;
			
		}
		
		
	}
	//attentio a l ordre, ".fa" est aussi present dans .fastq
	else if (_inputFilename.find(".fa") !=  string::npos || _inputFilename.find(".fasta") !=  string::npos) {
		getInfo()->add (0, "Input format: FastA");
		infoByte |= 0x01;
		_isFasta = true;
		
	}
	else
	{
		getInfo()->add (0, "Input format: unknown. Input extension must be one among fasta (.fa, .fasta) or fastq (.fq, .fastq)");
		return;
	}
	

	
	std::string leonversion = Stringify::format ("%i.%i.%i", LEON_VERSION_MAJOR, LEON_VERSION_MINOR,LEON_VERSION_PATCH);
	
	
    //Redundant from dsk solid file !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    _dskOutputFilename = getInput()->get(STR_URI_OUTPUT) ?
        getInput()->getStr(STR_URI_OUTPUT) + ".h5"  :
        System::file().getBaseName (_inputFilename) + ".h5"; //_inputFilename instead of prefix GR

#if 1

    /*************************************************/
    // Sorting count part
    /*************************************************/

    {
        /** We create a DSK instance and execute it. */
        SortingCountAlgorithm<> sortingCount (_inputBank, getInput());

        sortingCount.getInput()->add (0, STR_VERBOSE, getInput()->getStr(STR_VERBOSE));

        // best compression ratio if min-abundance never below 4 (ie. each kmer of the graph is used by at least 4 reads)
        getInput()->setInt (STR_KMER_ABUNDANCE_MIN_THRESHOLD, 4);

        sortingCount.execute();
    }

#endif

    /*************************************************/
    // We create the modified file
    /*************************************************/
    
    string dir = System::file().getDirectory(_inputFilename);
    string prefix = System::file().getBaseName(_inputFilename);
    //_outputFilename = dir + "/" + System::file().getBaseName(prefix) + ".leon";
	string baseOutputname;
	if(extension.find("gz") !=string::npos)
	{
		baseOutputname = dir + "/" + System::file().getBaseName(_inputFilename) ;
	}
	else
	{
		baseOutputname = _inputFilename;
	}
	_outputFilename = baseOutputname + ".leon";

//	_outputFile = System::file().newFile(_outputFilename, "wb");
	
	_storageH5file = StorageFactory(STORAGE_HDF5).create (_outputFilename, true, false,true);

	_groupLeon    = new tools::storage::impl::Group((*_storageH5file)().getGroup      ("leon"));
	_subgroupInfo = new tools::storage::impl::Group((*_storageH5file)().getGroup   ("metadata"));
	_subgroupDict = new tools::storage::impl::Group((*_storageH5file)().getGroup   ("leon/anchors"));
	_subgroupDNA  = new tools::storage::impl::Group((*_storageH5file)().getGroup    ("leon/dna"));

	if(! _isFasta)
		_subgroupQual  = new tools::storage::impl::Group((*_storageH5file)().getGroup    ("leon/qual"));

	
	_subgroupHeader  = new tools::storage::impl::Group((*_storageH5file)().getGroup    ("leon/header"));

	
	
	
	
	_subgroupInfoCollection = & _subgroupInfo->getCollection<math::NativeInt8> ("infobyte");

	
	
	if(_isFasta)
		_subgroupInfoCollection->addProperty ("type","fasta");
	else
		_subgroupInfoCollection->addProperty ("type","fastq");
	
	
	if(_noHeader)
		_subgroupInfoCollection->addProperty ("header","false");
	else
		_subgroupInfoCollection->addProperty ("header","true");


	_subgroupInfoCollection->addProperty ("version",leonversion);

	
	//making a block here so that ostream is immediately destroyed
	//otherwise bug since the referred to _subgroupInfo is destroyed in endcommpression
	//(I do not want to destroy _subgroupInfo in leon destructor, otherwise  compression will not be flushed until leon object is destroyed, bad behavior whan compression used within code)
	{
		tools::storage::impl::Storage::ostream osInfo (*_subgroupInfo, "infobyte");
		osInfo.write (reinterpret_cast<char const*>(&infoByte), sizeof(infoByte));
		osInfo.flush();
		
		tools::storage::impl::Storage::ostream osk (*_subgroupInfo, "kmerSize");
		osk.write (reinterpret_cast<char const*>(&_kmerSize), sizeof(_kmerSize));
		osk.flush();
	}
	
	
	if(! _isFasta)
	{
		//_FileQualname = baseOutputname + ".qual";
		//_FileQual = System::file().newFile(_FileQualname, "wb");
		//_Qual_outstream  =  new tools::storage::impl::Storage::ostream  (*_groupLeon, "qualities");
	}


	

	
#ifdef PRINT_DEBUG
	cout << "\tOutput filename: " << _outputFilename << endl;
	cout << "prefix " << prefix << endl;
	cout << "dir " << dir << endl;
	cout << "dskout  " << _dskOutputFilename << endl;
#endif

	
	
    //Compression
	if(! _noHeader)
	{
		startHeaderCompression();
	}
	
	



	startDnaCompression();
	

	
	
	endCompression();
}


void Leon::endQualCompression(){


	_qualCompRate = ((double)_compressed_qualSize / _input_qualSize);

	
}


void Leon::writeBlockLena(u_int8_t* data, u_int64_t size, int encodedSequenceCount,u_int64_t blockID){


	z_stream zs;
	memset(&zs, 0, sizeof(zs));
	
	//deflateinit2 to be able to gunzip it fro mterminal
	
	//if(deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
		//			(15+16), 8, Z_DEFAULT_STRATEGY) != Z_OK)
	
	
			if (deflateInit(&zs, Z_BEST_COMPRESSION) != Z_OK)
		throw Exception ("deflateInit failed while compressing.");
	
	zs.next_in = (Bytef*) data ;
	zs.avail_in = size ;           // set the z_stream's input
	
	int ret;
	char outbuffer[32768];
	std::string outstring;
	
	// retrieve the compressed bytes blockwise
	do {
		zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
		zs.avail_out = sizeof(outbuffer);
		
		ret = deflate(&zs, Z_FINISH);
		
		if (outstring.size() < zs.total_out) {
			// append the block to the output string
			outstring.append(outbuffer,
							 zs.total_out - outstring.size());
		}
	} while (ret == Z_OK);
	
	deflateEnd(&zs);
	
	/////////////////

	
	pthread_mutex_lock(&writeblock_mutex);
	
	std::string datasetname = Stringify::format ("qual_%i",blockID);
	
	tools::storage::impl::Storage::ostream os (*_subgroupQual, datasetname);
	os.write (reinterpret_cast<char const*>( outstring.data()), outstring.size());
	os.flush();

	std::string dsize = Stringify::format ("%i",outstring.size());
	auto _tempcollec = & _subgroupQual->getCollection<math::NativeInt8> (datasetname);
	_tempcollec->addProperty ("size",dsize);

	
	_input_qualSize += size;
	_compressed_qualSize +=  outstring.size();
	
//	if ((2*(blockID+1)) > _qualBlockSizes.size() )
//	{
//		_qualBlockSizes.resize(2*(blockID+1));
//	}
//	
//	_qualBlockSizes[2*blockID] =  outstring.size();
//	_qualBlockSizes[2*blockID+1] = encodedSequenceCount;
	
	
	pthread_mutex_unlock(&writeblock_mutex);

	
	
}

void Leon::writeBlock(u_int8_t* data, u_int64_t size, int encodedSequenceCount,u_int64_t blockID, bool Header){
	if(size <= 0) return;
	
	
	//cout << "\t-----------------------" << endl;
	//printf("____ write block %i ____ \n",blockID);
	//cout << "\tWrite block " << _blockSizes.size() << endl;
	//cout << "\tSequence " << encoder->_lastSequenceIndex-READ_PER_BLOCK << " - " << encoder->_lastSequenceIndex << endl;
	//cout << "Thread id: " << thread_id << endl;
	//cout << "\tEncoded size (byte): " << size << endl;
	
	
	
	pthread_mutex_lock(&writeblock_mutex);

	if(Header)
	{
		std::string datasetname = Stringify::format ("header_%i",blockID);
		
		tools::storage::impl::Storage::ostream os (*_subgroupHeader, datasetname);
		os.write (reinterpret_cast<char const*>( data), size);
		os.flush();
		
		std::string dsize = Stringify::format ("%i",size);
		auto _tempcollec = & _subgroupHeader->getCollection<math::NativeInt8> (datasetname);
		_tempcollec->addProperty ("size",dsize);
	}
	
	
	if(!Header)
	{
		std::string datasetname = Stringify::format ("dna_%i",blockID);
		
		tools::storage::impl::Storage::ostream os (*_subgroupDNA, datasetname);
		os.write (reinterpret_cast<char const*>( data), size);
		os.flush();
		
		std::string dsize = Stringify::format ("%i",size);
		auto _tempcollec = & _subgroupDNA->getCollection<math::NativeInt8> (datasetname);
		_tempcollec->addProperty ("size",dsize);
	}
	
	_compressedSize += size;

	//int thread_id = encoder->getId();

	if ((2*(blockID+1)) > _blockSizes.size() )
	{
		_blockSizes.resize(2*(blockID+1));
	}
	
	_blockSizes[2*blockID] = size ;
	_blockSizes[2*blockID+1] = encodedSequenceCount;
	
	
	pthread_mutex_unlock(&writeblock_mutex);
	
}
		
void Leon::endCompression(){
	//_rangeEncoder.flush();
	//_outputFile->fwrite(_rangeEncoder.getBuffer(true), _rangeEncoder.getBufferSize(), 1);
	
	//tools::storage::impl::Storage::ostream osInfo (*_groupLeon, STR_DATA_INFO);

	//osInfo.write (reinterpret_cast<char const*>(_rangeEncoder.getBuffer(true)), _rangeEncoder.getBufferSize()*sizeof(char));
	//osInfo.flush();

	
	
//	printf("_rangeEncoder buffer size %i B \n",_rangeEncoder.getBufferSize());
//	_outputFile->flush();
	
	u_int64_t inputFileSize = System::file().getSize(_inputFilename.c_str());
    getInfo()->add(0, "End compression");
    getInfo()->add(1, "Input file");
    getInfo()->add(2, "name", "%s", _inputFilename.c_str());
    getInfo()->add(2, "size", "%d bytes (%ld Mb)", inputFileSize, inputFileSize/1024LL/1024LL);
	
	u_int64_t outputFileSize = System::file().getSize(_outputFilename.c_str());
	
    getInfo()->add(1, "Output file");
	getInfo()->add(2, "name", "%s", _outputFilename.c_str());
    getInfo()->add(2, "size", "%d bytes (%ld Mb)", outputFileSize, outputFileSize/1024LL/1024LL);
	
    getInfo()->add(1, "Compression");
    gettimeofday(&_tim, NULL);
    _wfin_leon  = _tim.tv_sec +(_tim.tv_usec/1000000.0);
    
    getInfo()->add(2, "Time:",  "%.2f seconds", (  _wfin_leon - _wdebut_leon) );
    getInfo()->add(2, "Speed:",  "%.2f Mb/seconds", (System::file().getSize(_inputFilename)/1000000.0) / (  _wfin_leon - _wdebut_leon) );
    getInfo()->add(2, "Rates");
    getInfo()->add(3, "overall", "%.4f (%.4f)",
                   (float)((double)outputFileSize / (double)inputFileSize),
                   (float)((double)inputFileSize / (double)outputFileSize ));
	if(! _noHeader)
	{
        getInfo()->add(3, "header only", "%.4f (%.4f)",
                       (float)_headerCompRate,
                       (float) ((double)1/_headerCompRate));
	}
	else
	{
		getInfo()->add(3, "header completely discarded in '-seq-only' mode");
	}
    getInfo()->add(3, "Sequence only", "%.4f (%.4f)",
                   (float)_dnaCompRate,
                   (float)((double)1/_dnaCompRate));
	if( ! _isFasta)
	{
        getInfo()->add(3, "Quality only", "%.4f (%.4f) [%s mode]",
                       (float)_qualCompRate,
                       (float)((double)1/_qualCompRate),
                       _lossless?"lossless":"lossy");
	}
	
	delete  _groupLeon;
	delete  _subgroupInfo;
	delete _subgroupDict;
	delete _subgroupDNA;
	if(! _isFasta)
		delete _subgroupQual;
	delete _subgroupHeader;
	
	if(_storageH5file !=0)
	{
		delete _storageH5file;
		_storageH5file =0;
	}
	

	
	
	
	//printf("\tTime: %.2fs\n", (double)(clock() - _time)/CLOCKS_PER_SEC);
	//printf("\tSpeed: %.2f mo/s\n", (System::file().getSize(_inputFilename)/1000000.0) / ((double)(clock() - _time)/CLOCKS_PER_SEC));
}
		
		




void Leon::startHeaderCompression(){
    Iterator<Sequence>* itSeq = createIterator<Sequence> (
                                                          _inputBank->iterator(),
                                                          _inputBank->estimateNbItems(),
                                                          "Compressing headers"
                                                          );
    LOCAL(itSeq);
    
    
	_totalHeaderSize = 0;
	_compressedSize = 0;
	
	#ifdef PRINT_DEBUG
		cout << endl << "Start header compression" << endl;
    #endif
    
    //write first header to file and store it in _firstHeader variable
	//ifstream inputFileTemp(getInput()->getStr(STR_URI_FILE).c_str(), ios::in);
	//getline(inputFileTemp, _firstHeader);   //should be get comment from itseq
	//inputFileTemp.close();
	itSeq->first();
	_firstHeader = itSeq->item().getComment();
	_firstHeader.erase(_firstHeader.begin());
	itSeq->reset();

	#ifdef PRINT_DEBUG
		cout << "\tFirst Header: " << _firstHeader << endl;
		cout << "\tSize: " << _firstHeader.size() << endl;
	#endif
	
	_totalHeaderSize += _firstHeader.size();
	
	//encode the size of the first header on 2 byte and the header itself
	//CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _firstHeader.size());
	
	//for(int i=0; i < _firstHeader.size(); i++){
	//	_rangeEncoder.encode(_generalModel, _firstHeader[i]);
	//}
	
	
	createDataset(_subgroupHeader,"firstheadersize",_firstHeader.size());
	
	tools::storage::impl::Storage::ostream osH (*_subgroupHeader, "firstheader");
	osH.write (reinterpret_cast<char const*>(_firstHeader.data()), _firstHeader.size());
	osH.flush();
	
	
	
	//_rangeEncoder.flush();
	//_totalHeaderCompressedSize += _rangeEncoder.getBufferSize();
	//_outputFile->fwrite(_rangeEncoder.getBuffer(), _rangeEncoder.getBufferSize(), 1);
	//_rangeEncoder.clear();
	
	//cout << "Block start pos: " << _outputFile->tell() << endl;
	
	//iterate on read sequences and compress headers
	TIME_INFO (getTimeInfo(), "header compression");

	//int nb_threads_living = 0 ;
	
	#ifdef SERIAL
		setDispatcher (new SerialDispatcher());
	#else
		setDispatcher (  new Dispatcher (_nb_cores) );
	#endif

	//getDispatcher()->iterate (itSeq,  HeaderEncoder(this, &nb_threads_living), 10000);
	getDispatcher()->iterate (itSeq,  HeaderEncoder(this), this->getReadPerBlock());
	endHeaderCompression();
}




void Leon::endHeaderCompression(){
	//u_int64_t descriptionStartPos = _outputFile->tell();
	//cout << "Description start pos: " << descriptionStartPos << endl;
	
	//CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _blockSizes.size());
	//for(int i=0; i<_blockSizes.size(); i++){
		//cout << "block size: " << _blockSizes[i] << endl;
	//	CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _blockSizes[i]);
	//}
	

	createDataset(_subgroupHeader,"nb_blocks",_blockSizes.size());

	getInfo()->add(0, "End Header compression");
  getInfo()->add(1, "# blocks", "%i", _blockSizes.size());
	tools::storage::impl::Storage::ostream os (*_subgroupHeader, "blocksizes");
	os.write (reinterpret_cast<char const*>(_blockSizes.data()), _blockSizes.size()*sizeof(u_int64_t));
	os.flush();
	

		
	_headerCompRate = ((double)_compressedSize / _totalHeaderSize);
	
	//cout << "\t\tData blocks count: " << _blockSizes.size() << endl;
	//cout << "\tBlock data size: " << _rangeEncoder.getBufferSize() << endl;
  getInfo()->add(1, "headers size", "%u", _totalHeaderSize);
  getInfo()->add(1, "headers compressed size", "%u", _compressedSize);
  getInfo()->add(1, "compression rate", "%.4f", (float)(_headerCompRate));
	//_rangeEncoder.clear();
	_blockSizes.clear();
}


void Leon::startDnaCompression(){
	#ifdef PRINT_DEBUG
		cout << endl << "Start reads compression" << endl;
    #endif
    
	

	//Create and fill bloom
    createBloom ();
   // LOCAL (_bloom); //now we need it later
    
	int64_t nbestimated = _inputBank->estimateNbItems();
	


//	_anchorKmers = new Hash16<kmer_type, u_int32_t > ( (nbestimated/10) *  sizeof(u_int32_t)  *10LL /1024LL / 1024LL ); // hmm  Hash16 would need a constructor with sizeof main entry //maybe *2 for low coverage dataset
	u_int64_t nbcreated ;
	_anchorKmers = new Hash16<kmer_type, u_int32_t > ( nbestimated/10 , &nbcreated ); //creator with nb entries given
//	printf("asked %lli entries, got %llu \n",nbestimated/10 ,nbcreated);
	
    Iterator<Sequence>* itSeq = createIterator<Sequence> (
                                                          _inputBank->iterator(),
                                                          nbestimated,
                                                          "Compressing dna"
                                                          );
    LOCAL(itSeq);
    
	//create a temporary output file to store the anchors dict
	//_dictAnchorFile = System::file().newFile(_outputFilename + ".adtemp", "wb"); 
	_dictAnchorFile = new ofstream((_outputFilename + ".adtemp").c_str(), ios::out|ios::binary);
	
	_lastAnchorValue = 0;
	_anchorAdress = 0;
	_totalDnaSize = 0;
	//_totalDnaCompressedSize = 0;
	_compressedSize = 0;
	
	//createKmerAbundanceHash();
    
	//iterate on read sequences and compress headers
	TIME_INFO (getTimeInfo(), "DNA compression");

	//int nb_threads_living = 0 ;
	
	#ifdef SERIAL
		setDispatcher (new SerialDispatcher());
	#else
		setDispatcher (  new Dispatcher (_nb_cores) );
	#endif

	//getDispatcher()->iterate (itSeq,  HeaderEncoder(this, &nb_threads_living), 10000);
	getDispatcher()->iterate (itSeq,  DnaEncoder(this), this->getReadPerBlock());
	
	endDnaCompression();
	
	if(! _isFasta)
	{
		endQualCompression();
	}
}


void Leon::endDnaCompression(){
	
	
	createDataset(_subgroupDNA,"nb_blocks",_blockSizes.size());
	
  getInfo()->add(0, "End Sequence compression");
  getInfo()->add(1, "# blocks", "%u", _blockSizes.size());
	tools::storage::impl::Storage::ostream os (*_subgroupDNA, "blocksizes");
	os.write (reinterpret_cast<char const*>(_blockSizes.data()), _blockSizes.size()*sizeof(u_int64_t));
	os.flush();

	_blockSizes.clear();
	
	writeBloom();
	writeAnchorDict();
	
	_dnaCompRate = ((double)_compressedSize / _totalDnaSize);
	
  getInfo()->add(1, "# sequences", "%u", _readCount);
  getInfo()->add(1, "# nucleotides", "%u", _totalDnaSize);

	createDataset(_subgroupInfo,"readcount",_readCount);
	createDataset(_subgroupInfo,"totalDnaSize",_totalDnaSize);
	createDataset(_subgroupInfo,"minSequenceSize",_minSequenceSize);
	createDataset(_subgroupInfo,"maxSequenceSize",_maxSequenceSize);
	
  getInfo()->add(1, "Compression rates");
  getInfo()->add(2, "overall", "%.4f (%u)", (float)_dnaCompRate, _compressedSize);

  getInfo()->add(2, "Bloom", "%.2f (%u)", ((_bloom->getSize()*100) / (double)_compressedSize), _bloom->getSize());
  getInfo()->add(2, "Anchors dict", "%.2f (%u) (%u entries)", ((_anchorDictSize*100) / (double)_compressedSize), _anchorDictSize, _anchorAdress);

  u_int64_t readsSize = _anchorAdressSize+_anchorPosSize+_readSizeSize+_bifurcationSize+_otherSize;
  getInfo()->add(2, "Reads", "%.2f (%u)", ((readsSize*100) / (double)_compressedSize), readsSize);
  getInfo()->add(3, "Anchor adress", "%.2f (%u)", ((_anchorAdressSize*100) / (double)_compressedSize), _anchorAdressSize);
  getInfo()->add(3, "Anchor pos", "%.2f (%u)", ((_anchorPosSize*100) / (double)_compressedSize), _anchorPosSize);
  getInfo()->add(3, "Read size", "%.2f (%u)", ((_readSizeSize*100) / (double)_compressedSize), _readSizeSize);
  getInfo()->add(3, "Bifurcation", "%.2f (%u)", ((_bifurcationSize*100) / (double)_compressedSize), _bifurcationSize);
  getInfo()->add(3, "Other (N, error, infoBits)", "%.2f (%u)", ((_otherSize*100) / (double)_compressedSize), _otherSize);
  getInfo()->add(2, "Read without anchor", "%.2f (%u)", ((_noAnchorSize*100) / (double)_compressedSize), _noAnchorSize);
  
  if(_anchorAdress!=0){
  getInfo()->add(1, "Reads per anchor", "%u",  _readCount / _anchorAdress);
  }
  getInfo()->add(1, "Read without anchor", "%.2f", ((double)_readWithoutAnchorCount*100) / _readCount);
  getInfo()->add(1, "De Bruijn graph");
	
	getInfo()->add(2, "Simple path", "%.2f", ((_MCuniqSolid*100)/ (double)_MCtotal));
	getInfo()->add(2, "Bifurcation", "%.2f", ((_MCmultipleSolid*100)/(double)_MCtotal));
	getInfo()->add(2, "Break", "%.2f", ((_MCnoAternative*100)/(double)_MCtotal));
	getInfo()->add(2, "Error", "%.2f", ((_MCuniqNoSolid*100)/(double)_MCtotal));
	

	delete _anchorKmers;
	System::file().remove(_dskOutputFilename);


}

void Leon::writeBloom(){
	_compressedSize += _bloom->getSize();
	

	StorageTools::singleton().saveBloom<kmer_type> (_storageH5file->getGroup(this->getName()), "bloom", _bloom, _kmerSize);

}



void Leon::writeAnchorDict(){

	_anchorRangeEncoder.flush();
	
	//todo check if the tempfile _dictAnchorFile may be avoided (with the use of hdf5 ?)
	_dictAnchorFile->write( (const char*) _anchorRangeEncoder.getBuffer(), _anchorRangeEncoder.getBufferSize());
	_dictAnchorFile->flush();
	_dictAnchorFile->close();
	_anchorRangeEncoder.clear();
	
	
	u_int64_t size = System::file().getSize(_outputFilename + ".adtemp");
	_anchorDictSize = size;
	//u_int64_t size = _anchorRangeEncoder.getBufferSize();
	_compressedSize += size;
	

	
	createDataset(_subgroupDict,"size",size);
	createDataset(_subgroupDict,"anchorAdress",_anchorAdress);

	
	//_dictAnchorFile->seeko(0, SEEK_SET);
	//_outputFile->fwrite(_dictAnchorFile, size, 1);
	ifstream tempFile((_outputFilename + ".adtemp").c_str(), ios::in|ios::binary);
	
	
	
	
	int bufsize = 4096*8;
	char * buffer = new char [bufsize];
	
	tools::storage::impl::Storage::ostream osD (*_subgroupDict, "anchorsDict");

    while (tempFile.good()) {
		tempFile.read(buffer, bufsize);
       // _outputFile->fwrite(buffer, tempFile.gcount(), 1);
		osD.write (reinterpret_cast<char const*>(buffer), tempFile.gcount());

    }
	
	osD.flush();

	tempFile.close();
	remove((_outputFilename + ".adtemp").c_str());
	delete [] buffer;

}

bool Leon::anchorExist(const kmer_type& kmer, u_int32_t* anchorAdress){
	
	if (_anchorKmers->get(kmer,anchorAdress)) //avec Hash16
	{
		return true;
	}

	return false;

}



void Leon::updateMinMaxSequenceSize(int newMin, int newMax)
{
	pthread_mutex_lock(&minmax_mutex);
	_minSequenceSize = std::min(_minSequenceSize, newMin);
	_maxSequenceSize = std::max(_maxSequenceSize, newMax);
	pthread_mutex_unlock(&minmax_mutex);

}

int Leon::findAndInsertAnchor(const vector<kmer_type>& kmers, u_int32_t* anchorAdress){
	
	pthread_mutex_lock(&findAndInsert_mutex);

		
	//cout << "\tSearching and insert anchor" << endl;
	int maxAbundance = -1;
	int bestPos;
	kmer_type bestKmer;
	//int bestPos = -1;
	

	kmer_type kmer, kmerMin;
	
	/*
	////////////
	for(int i=0; i<kmers.size(); i++){
		kmer = kmers[i];
		kmerMin = min(kmer, revcomp(kmer, _kmerSize));
		if(_bloom->contains(kmerMin)){
			encodeInsertedAnchor(kmerMin);
			_anchorKmers.insert(kmerMin, _anchorAdress);
			*anchorAdress = _anchorAdress;
			_anchorKmerCount += 1;
			_anchorAdress += 1;
			return i;
		}
	}
	return -1;
	/////////////////////*/


	//int iMin = 40;
	//int iMax = 60;
	int iMin = kmers.size()/2;
	int iMax = kmers.size()/2 + 10;
	//cout << iMin << "  " << iMax << endl;
	iMin = std::max(iMin, 0);
	iMax = std::min(iMax, (int) kmers.size());

	for(int i=iMin; i<iMax; i++){

		kmer = kmers[i];
		kmerMin = min(kmer, revcomp(kmer, _kmerSize));

		if(_bloom->contains(kmerMin)){
			maxAbundance = 0;
			bestPos = i;
			bestKmer = kmerMin;
			break;
		}
	}

	if(maxAbundance == -1){

		for(int i=0; i<iMin; i++){
			kmer = kmers[i];
			kmerMin = min(kmer, revcomp(kmer, _kmerSize));


			if(_bloom->contains(kmerMin)){
				maxAbundance = 0;
				bestPos = i;
				bestKmer = kmerMin;
				break;
			}
		}


		if(maxAbundance == -1){
			for(unsigned int i=iMax; i<kmers.size(); i++){
				kmer = kmers[i];
				kmerMin = min(kmer, revcomp(kmer, _kmerSize));

				if(_bloom->contains(kmerMin)){
					maxAbundance = 0;
					bestPos = i;
					bestKmer = kmerMin;
					break;
				}
			}
		}
	}


	/*
	for(int i=0; i<kmers.size(); i++){

		kmer = kmers[i];
		kmerMin = min(kmer, revcomp(kmer, _kmerSize));

		int abundance;
		if(_kmerAbundance->get(kmerMin, &abundance)){
			if(abundance > maxAbundance){
				maxAbundance = abundance;// + ((kmers.size()-i)*2);
				bestKmer = kmerMin;
				bestPos = i;
				//cout << maxAbundance << endl;
				//cout << bestPos << " " << abundance << " " << kmer.toString(_kmerSize) << " " << revcomp(kmer, _kmerSize).toString(_kmerSize) << endl;
			}
			//cout << abundance << endl;
		}
		else if(maxAbundance == -1 && _bloom->contains(kmerMin)){
			maxAbundance = _nks;
			bestKmer = kmerMin;
			bestPos = i;
			//cout << maxAbundance << endl;
		}
	}*/
	
	if(maxAbundance == -1)
	{
		pthread_mutex_unlock(&findAndInsert_mutex);
		return -1;
	}

	encodeInsertedAnchor(bestKmer);
	
	_anchorKmers->insert(bestKmer,_anchorAdress); //with Hash16
	//_anchorKmers[bestKmer] = _anchorAdress;
	//_anchorKmers.insert(bestKmer, _anchorAdress);
		
	
	*anchorAdress = _anchorAdress;
	//_anchorKmerCount += 1;
	_anchorAdress += 1;
	
	/*
	int val;
	for(int i=0; i<kmers.size(); i++){
		kmer = kmers[i];
		kmerMin = min(kmer, revcomp(kmer, _kmerSize));
		_kmerAbundance->remove(kmerMin, &val);
		//_kmerAbundance->insert(kmerMin, val-1);
	}*/

	pthread_mutex_unlock(&findAndInsert_mutex);
	return bestPos;
}

void Leon::encodeInsertedAnchor(const kmer_type& kmer){

	//static int i = 0;
	
	string kmerStr = kmer.toString(_kmerSize);

	for(unsigned int i=0; i<kmerStr.size(); i++){
		_anchorRangeEncoder.encode(_anchorDictModel, Leon::nt2bin(kmerStr[i]));
	}
	//i+= 1;
	//cout << i << endl;
	if(_anchorRangeEncoder.getBufferSize() >= 4096){
		_dictAnchorFile->write((const char*) _anchorRangeEncoder.getBuffer(), _anchorRangeEncoder.getBufferSize());
		_anchorRangeEncoder.clearBuffer();
	}
}











void * decoder_all_thread(void * args)
{
	
	thread_arg_decoder * targ = (thread_arg_decoder*) args;
	QualDecoder * qual_decoder = targ->qual_decoder;
	DnaDecoder * dna_decoder = targ->dna_decoder;
	HeaderDecoder * header_decoder = targ->header_decoder;
	
	if(qual_decoder!=NULL)
		qual_decoder->execute();
	
	if(header_decoder!=NULL)
		header_decoder->execute();

	dna_decoder->execute();
	
	pthread_exit(0);
}



void * decoder_dna_thread(void * args)
{
	DnaDecoder * dna_decoder = (DnaDecoder*) args;
	dna_decoder->execute();
 	pthread_exit(0);
}

void * decoder_qual_thread(void * args)
{
	QualDecoder * qual_decoder = (QualDecoder*) args;
	qual_decoder->execute();
	pthread_exit(0);
}

void Leon::executeDecompression(){
	

	_filePos = 0;
	
	#ifdef PRINT_DEBUG
    if(!_iterator_mode){
	  cout << "Start decompression" << endl;
    }
    #endif
    
    _inputFilename = getInput()->getStr(STR_URI_FILE);
	//string inputFilename = prefix + ".txt"; //".leon"
	//_outputFile = System::file().newFile(outputFilename, "wb");
    
    #ifdef PRINT_DEBUG
    if(!_iterator_mode){
	  cout << "\tInput filename: " << _inputFilename << endl;
    }
    #endif

    if (!System::file().doesExist(_inputFilename)){
        std::stringstream ss;
        ss << "File not found: " << _inputFilename;
        throw Exception (ss.str().c_str());
    }
    _storageH5file = StorageFactory(STORAGE_HDF5).create (_inputFilename,false,false,true); //open without adding extension h5
	_groupLeon = new tools::storage::impl::Group((*_storageH5file)().getGroup   ("leon"));
	_subgroupInfo = new tools::storage::impl::Group((*_storageH5file)().getGroup   ("metadata"));
	_subgroupDict = new tools::storage::impl::Group((*_storageH5file)().getGroup   ("leon/anchors"));
	_subgroupDNA = new tools::storage::impl::Group((*_storageH5file)().getGroup   ("leon/dna"));
	
	

	
	_subgroupInfoCollection = & _subgroupInfo->getCollection<math::NativeInt8> ("infobyte");

	
	tools::storage::impl::Storage::istream isInfo (*_subgroupInfo, "infobyte");
	
	
	
	
	
	string dir = System::file().getDirectory(_inputFilename);
	
	//_inputFile = new ifstream(_inputFilename.c_str(), ios::in|ios::binary);
	

	
	//remove .leon at the end :
	string inputFilename_leon_removed = System::file().getBaseName(_inputFilename);

	
	//Go to the end of the file to decode blocks informations, data are read in reversed order (from right to left in the file)
	//The first number is the number of data blocks

	
	
	//Output file
	string prefix = System::file().getBaseName(inputFilename_leon_removed); // for prefix need to remove two dots : .fastq.leon , but not all of them if another dot in the filename (it happened, mail from angry users)

	_outputFilename = dir + "/" + prefix;
	
	//printf("_outputFilename %s \n",_outputFilename.c_str());

	//Decode the first byte of the compressed file which is an info byte
	u_int8_t infoByte;
	//infoByte = _rangeDecoder.nextByte(_generalModel);
	
	
	isInfo.read (reinterpret_cast<char*> (&infoByte),sizeof(infoByte));

	



	
	//the first bit holds the file format. 0: fastq, 1: fasta
	//_isFasta = ((infoByte & 0x01) == 0x01);
	
	
	
	//Second bit : option no header
	//_noHeader = ((infoByte & 0x02) == 0x02);
	
	
	
	std::string  filetype = _subgroupInfoCollection->getProperty("type");
	if(filetype == "fasta")
		_isFasta = true;
	else
		_isFasta = false;
	
	std::string  headerinfo = _subgroupInfoCollection->getProperty("header");
	if(headerinfo == "true")
		_noHeader = false;
	else
		_noHeader = true;

	//printf("%s %s \n",filetype.c_str(),headerinfo.c_str() );

	
	if(! _isFasta)
		_subgroupQual  = new tools::storage::impl::Group((*_storageH5file)().getGroup    ("leon/qual"));
	

	_subgroupHeader  = new tools::storage::impl::Group((*_storageH5file)().getGroup    ("leon/header"));

	
	
	//printf("info byte %i  _noHeader %i  _isFasta %i \n",infoByte,_noHeader,_isFasta);
	
	if(! _isFasta)
	{
		
	//	_FileQualname =    dir + "/" +  System::file().getBaseName(inputFilename_leon_removed) + ".fastq.qual";
	//	_inputFileQual = new ifstream(_FileQualname.c_str(), ios::in|ios::binary);
	//	cout << "\tQual filename: " << _FileQualname << endl;
	}
	
  if(!_iterator_mode){
    getInfo()->add(0, "Decompression");
    getInfo()->add(1, "Input filename", "%s", _inputFilename.c_str());
  }

	if(_noHeader)
	{
		if(!_iterator_mode)
        getInfo()->add(1, "Headers were not stored, will number reads.");
	}
	
	if(_isFasta){
		if(!_iterator_mode)
        getInfo()->add(1, "Output format", "FastA");
		_outputFilename += ".fasta.d";
	}
	else{
		if(!_iterator_mode)
        getInfo()->add(1, "Output format", "FastQ");
		_outputFilename += ".fastq.d";
	}
	
	if(!_iterator_mode)
		_outputFile = System::file().newFile(_outputFilename, "wb");
	
	
	//Get kmer size
	//_kmerSize = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
	
	tools::storage::impl::Storage::istream isk (*_subgroupInfo, "kmerSize");
	isk.read (reinterpret_cast<char*> (&_kmerSize),sizeof(_kmerSize));
	
	if(!_iterator_mode)
	{
  	getInfo()->add(1, "Kmer size", "%i", _kmerSize);
	}

	std::string  leonversion =  _subgroupInfoCollection->getProperty ("version");

	if(!_iterator_mode)
	  getInfo()->add(1, "Input File compressed with Leon", "%s", leonversion.c_str());

	//cout << "\tInput File was compressed with leon version " << version_major << "."  << version_minor << "."  << version_patch  << endl;
	
	
	//	if(version_major != LEON_VERSION_MAJOR || version_minor != LEON_VERSION_MINOR  || version_patch != LEON_VERSION_PATCH )
	//	{
	//		cout << "\tWarning diff version "   << endl;
	//	}
	
	if(! _iterator_mode)
	{
		startDecompressionAllStreams();

		endDecompression();
	}

	
}

void Leon::testing_iter(){
	//printf("testing iterator \n");
	
	tools::dp::Iterator<Sequence>* iterl = new LeonIterator(*this);
	
	/*
	iterl->first();
	
	std::cout << "[" << (*iterl)->getDataSize() << "] " << (*iterl)->getComment()  << std::endl;
	std::cout << (*iterl)->toString() << std::endl;
	
	iterl->next();
	std::cout << "[" << (*iterl)->getDataSize() << "] " << (*iterl)->getComment()  << std::endl;
	std::cout << (*iterl)->toString() << std::endl;*/

	for(iterl->first();!iterl->isDone();iterl->next())
	{
		//std::cout << "[" << (*iterl)->getDataSize() << "] " << (*iterl)->getComment()  << std::endl;
		std::cout << (*iterl)->getComment() << std::endl;
		std::cout << (*iterl)->toString() << std::endl;
		std::cout << (*iterl)->getQuality() << std::endl;

	}

}

void Leon::startDecompression_setup(){
	
	_filePosHeader = 0;
	_filePosDna = 0;
	
	if(! _noHeader)
	{
		///////// header setup  /////////
		
		size_t firstHeaderSize;
		readDataset(_subgroupHeader,"firstheadersize",firstHeaderSize);
		
		
		char * tempS = (char * ) malloc(firstHeaderSize+1);
		tools::storage::impl::Storage::istream isH (*_subgroupHeader, "firstheader");
		isH.read (reinterpret_cast<char *>(tempS), firstHeaderSize);
		
		tempS[firstHeaderSize] = '\0';
		_firstHeader = std::string(tempS);
		free(tempS);
		
		
		//setup header block sizes
		size_t nb_blocks;
		readDataset(_subgroupHeader,"nb_blocks",nb_blocks);
		_headerBlockSizes.resize(nb_blocks,0);
		
		tools::storage::impl::Storage::istream is (*_subgroupHeader, "blocksizes");
		is.read (reinterpret_cast<char *>(_headerBlockSizes.data()), _headerBlockSizes.size()*sizeof(u_int64_t));
		////
		
	}
	
	
	/////// dna setup ////////////
	
	//need to init _filePosDna here
	for(unsigned int ii=0; ii<_headerBlockSizes.size(); ii+=2 )
	{
		_filePosDna += _headerBlockSizes[ii];
	}
	
	//setup dna block sizes
	size_t nb_blocks;
	readDataset(_subgroupDNA,"nb_blocks",nb_blocks);
	_blockCount = nb_blocks;
	_dnaBlockSizes.resize(nb_blocks,0);
	
	tools::storage::impl::Storage::istream is (*_subgroupDNA, "blocksizes");
	is.read (reinterpret_cast<char *>(_dnaBlockSizes.data()), _dnaBlockSizes.size()*sizeof(u_int64_t));
	////
	
	_kmerModel = new KmerModel(_kmerSize);
	
	decodeBloom();
	decodeAnchorDict();
	
	/////////// qualities setup //////////
	if(! _isFasta)
	{
		_filePosQual =0;
	}
	///////////////
	
}


void Leon::decoders_setup(){


	
	for(int i=0; i<_nb_cores; i++){
		
		if(! _isFasta)
		{
			//QualDecoder* qd = new QualDecoder(this, _FileQualname,_groupLeon);
			QualDecoder* qd = new QualDecoder(this, "qualities",_subgroupQual);
			
			_qualdecoders.push_back(qd);
		}
		
		DnaDecoder* dd = new DnaDecoder(this, _inputFilename, _subgroupDNA);
		_dnadecoders.push_back(dd);
		
		if(! _noHeader)
		{
			//HeaderDecoder* hd = new HeaderDecoder(this, _inputFilename);
			HeaderDecoder* hd = new HeaderDecoder(this, _inputFilename, _subgroupHeader);
			
			_headerdecoders.push_back(hd);
		}
	}
	
	
	_tab_threads = new pthread_t [_nb_cores];
	_targ = new thread_arg_decoder [_nb_cores];
	
}

void Leon::decoders_cleanup(){
	
	for(unsigned int i=0; i<_dnadecoders.size(); i++){
		delete _dnadecoders[i];
	}
	_dnadecoders.clear();
	
	for(unsigned int i=0; i<_headerdecoders.size(); i++){
		delete _headerdecoders[i];
	}
	_headerdecoders.clear();
	
	for(unsigned int i=0; i<_qualdecoders.size(); i++){
		delete _qualdecoders[i];
	}
	_qualdecoders.clear();
}

void Leon::decompressionDecodeBlocks(unsigned int & idx, int & livingThreadCount){


	for(int j=0; j<_nb_cores; j++){
		

		if(idx >= _dnaBlockSizes.size()) break;
		
		int blockId = idx/2 ;

		
		u_int64_t blockSize;
		int sequenceCount;
		
		livingThreadCount = j+1;
		
		QualDecoder* qdecoder;
		HeaderDecoder* hdecoder;
		DnaDecoder* ddecoder;
		
		//header decoder
		if(! _noHeader)
		{
			blockSize = _headerBlockSizes[idx];
			sequenceCount = _headerBlockSizes[idx+1];
			hdecoder = _headerdecoders[j];
			hdecoder->setup(_filePosHeader, blockSize, sequenceCount,blockId);
			_filePosHeader += blockSize;
		}
		else
		{
			hdecoder= NULL;
		}
		
		//dna decoder
		blockSize = _dnaBlockSizes[idx];
		sequenceCount = _dnaBlockSizes[idx+1];
		ddecoder = _dnadecoders[j];
		ddecoder->setup(_filePosDna, blockSize, sequenceCount,blockId);
		_filePosDna += blockSize;
		
		//qual decoder setup
		//here test if in fastq mode, put null pointer otherwise
		if(! _isFasta)
		{
			qdecoder = _qualdecoders[j];
			qdecoder->setup( blockId);
		}
		else
		{
			qdecoder= NULL;
		}
		
		
		_targ[j].qual_decoder = qdecoder;
		_targ[j].dna_decoder = ddecoder;
		_targ[j].header_decoder = hdecoder;
		
		pthread_create(&_tab_threads[j], NULL, decoder_all_thread, _targ + j);
		
		idx += 2;
		
		if(! _iterator_mode)
			this->_progress_decode->inc(1);
	}
}
void Leon::startDecompressionAllStreams(){
	

	startDecompression_setup();
	
	switch (getInput()->getInt(STR_VERBOSE))
	{
		case 0: default:    _progress_decode =  new IteratorListener ();break;
		case 1:             _progress_decode = new ProgressSynchro ( new ProgressTimer ( _blockCount/2, "Decompressing all streams"), System::thread().newSynchronizer()   );break;
			
		case 2:             _progress_decode = new ProgressSynchro ( new Progress ( _blockCount/2, "Decompressing all streams"), System::thread().newSynchronizer()   );break;
	}
	
	
	getInfo()->add(1, "Block count", "%u", _blockCount/2);


//	delete _progress_decode;
//	_progress_decode = new ProgressSynchro ( new ProgressTimer ( _blockCount/2, "Decompressing all streams"), System::thread().newSynchronizer()   );
//	_progress_decode = new ProgressSynchro ( new Progress ( _blockCount/2, "Decompressing all streams"), System::thread().newSynchronizer()   );

	_progress_decode->init();
	
	
	decoders_setup();
	

	unsigned int i = 0;
	int livingThreadCount = 0;
	
	
	
	while(i < _dnaBlockSizes.size()){
		
		//decode blocks

		decompressionDecodeBlocks(i,livingThreadCount); //this will increment i

		for(int j=0; j < livingThreadCount; j++){
			
			pthread_join(_tab_threads[j], NULL);
			
			HeaderDecoder* hdecoder = NULL;
			QualDecoder* qdecoder = NULL;
			DnaDecoder* ddecoder = _dnadecoders[j];
			

			std::istringstream  * stream_qual = NULL;
			std::istringstream  * stream_header = NULL;
			std::istringstream  * stream_dna = NULL;

			if(! _isFasta)
			{
				qdecoder = _qualdecoders[j];
				stream_qual = new std::istringstream (qdecoder->_buffer);
				qdecoder->_buffer.clear();
			}
			
			if(! _noHeader)
			{
				hdecoder = _headerdecoders[j];
				stream_header = new std::istringstream (hdecoder->_buffer);
				hdecoder->_buffer.clear();
			}
			
			stream_dna = new std::istringstream (ddecoder->_buffer);

			ddecoder->_buffer.clear();

			
			std::string line;
			std::string output_buff;

			output_buff.reserve(this->getReadPerBlock() * 500);
			
			bool reading = true;
			
			
			
			u_int64_t readid=0;
			while(reading){
				
				stringstream sint;
				sint << readid;
				
				if( ! _noHeader)
				{
					if(getline(*stream_header, line)){
						if(_isFasta)
							output_buff += ">";
						else
							output_buff += "@";
						
						output_buff +=  line + '\n';
					}
					else
						reading = false;
				}
				else
				{
					if(_isFasta)
						output_buff += "> " + sint.str() + '\n';
					else
						output_buff += "@ " + sint.str() + '\n';
					
					readid++;
				}
				 
				
				
				if(getline(*stream_dna, line)){
					output_buff +=  line + '\n';
				}
				else
					reading = false;
				
				
				if( ! _isFasta)
				{
					if(getline(*stream_qual, line)){
						output_buff += "+\n";
						output_buff +=  line + '\n';
					}
					else
						reading = false;
				}
				 
			}
			
			 
			_outputFile->fwrite(output_buff.c_str(), output_buff.size(), 1);

			if(stream_qual!= NULL) delete  stream_qual;
			if(stream_header!= NULL) delete  stream_header;
			if(stream_dna!= NULL) delete  stream_dna;
			
			
		}
		
		livingThreadCount = 0;
	}
	
	
	
	_outputFile->flush();

	for(unsigned int i=0; i<_dnadecoders.size(); i++){
		delete _dnadecoders[i];
	}
	_dnadecoders.clear();
	
	for(unsigned int i=0; i<_headerdecoders.size(); i++){
		delete _headerdecoders[i];
	}
	_headerdecoders.clear();
	
	for(unsigned int i=0; i<_qualdecoders.size(); i++){
		delete _qualdecoders[i];
	}
	_qualdecoders.clear();
	
	
	
	delete [] _tab_threads;
	delete [] _targ;
	cout << endl;
	
	
	_progress_decode->finish();
	
	delete _kmerModel;

}




void Leon::setupNextComponent( 		vector<u_int64_t>   & blockSizes    ){
	//Go to the data block position (position 0 for headers, position |headers data| for reads)
	_inputFile->seekg(_filePos, _inputFile->beg);
	
	blockSizes.clear();
	//u_int64_t size = 0;
	
	_blockCount = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
	for(unsigned int i=0; i<_blockCount; i++){
		u_int64_t blockSize = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
		blockSizes.push_back(blockSize);
		//size += blockSize;
	}
	
	

	
	//cout << "\tBlock count: " << _blockCount/2 << endl;
	/*
	for(int i=0; i<_blockSizes.size(); i++){
		cout << _blockSizes[i] << " ";
	}
	cout << endl;*/
	
}



void Leon::decodeBloom(){

	//to be removed
	////////
	//pos = tous les block sizes des header
	
	/*
	u_int64_t total_header_block_size = 0 ;
	
	for(int ii=0; ii<_headerBlockSizes.size(); ii+=2 )
	{
		total_header_block_size  += _headerBlockSizes[ii];
	}
	
	u_int64_t bloomPos =  total_header_block_size ;  

	for(int i=0; i<_dnaBlockSizes.size(); i++){
		bloomPos += _dnaBlockSizes[i];
		i += 1;
	}

	_inputFile->seekg(bloomPos, _inputFile->beg);
	*/

    
	//u_int64_t bloomBitSize = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
	//u_int64_t bloomHashCount = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
	

	//_bloom   =  new BloomNeighborCoherent<kmer_type> (bloomBitSize,_kmerSize,bloomHashCount);
	

//	_inputFile->read((char*)_bloom->getArray(), _bloom->getSize());
//////
	
	_bloom   = StorageTools::singleton().loadBloom<kmer_type>     (*_groupLeon,   "bloom");


#ifdef PRINT_DEBUG_DECODER
//		cout << "Bloom size: " << _bloom->getSize() << endl;
//		cout << "Anchor dict pos: " << _inputFile->tellg() << endl;
	#endif
	

}

void Leon::decodeAnchorDict(){
	#ifdef PRINT_DEBUG_DECODER
		cout << "\tDecode anchor dict" << endl;
	#endif
	
	
	u_int64_t anchorDictSize;
	u_int32_t anchorCount;
	readDataset(_subgroupDict,"size",anchorDictSize);
	readDataset(_subgroupDict,"anchorAdress",anchorCount);
	
	
	tools::storage::impl::Storage::istream isD (*_subgroupDict, "anchorsDict");
	
	//_anchorRangeDecoder.setInputFile(_inputFile);
	_anchorRangeDecoder.setInputFile(&isD); //seems to be working ok
	
	string anchorKmer = "";
	
//	u_int64_t dictPos = _inputFile->tellg();
	
	//KmerModel model(_kmerSize, KMER_DIRECT);

	u_int64_t currentAnchorCount = 0;
	
	while(currentAnchorCount < anchorCount){


		u_int8_t c = _anchorRangeDecoder.nextByte(_anchorDictModel);
		anchorKmer += Leon::bin2nt(c); //convert to char
		if(anchorKmer.size() == _kmerSize){
			

			//cout << anchorKmer << endl;
			//if(i<=10) cout << anchorKmer << endl;
			//cout << "1: " << anchorKmer << endl;
			
			kmer_type kmer = _kmerModel->codeSeed(anchorKmer.c_str(), Data::ASCII).value() ; //then convert to bin
			

			//could be optimized if needed
			//cout << "2: " << model.toString(kmer) << endl;
			//lala += 1;
			_vecAnchorKmers.push_back(kmer);
			
			anchorKmer.clear();
			
			currentAnchorCount += 1;

		}
	}
	
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tAnchor count: " << _vecAnchorKmers.size() << endl;
	#endif
	

}



kmer_type Leon::getAnchor(ifstream* anchorDictFile, u_int32_t adress){
	
	return _vecAnchorKmers[adress];  //icii
	
	//anchorDictFile->seekg(_kmerSize*adress);
	
	//char buffer[_kmerSize];
	
	//anchorDictFile->read(buffer, _kmerSize);
	//kmer_type kmer = model.codeSeed(anchorKmer.c_str(), Data::ASCII);
	//return _vecAnchorKmers[adress];
	//return _kmerModel->codeSeed(buffer, Data::ASCII);
}

void Leon::endDecompression(){
	
	getInfo()->add(1, "Output filename", "%s", _outputFile->getPath().c_str());

	gettimeofday(&_tim, NULL);
	_wfin_leon  = _tim.tv_sec +(_tim.tv_usec/1000000.0);
	
	getInfo()->add(1, "Time", "%.2f seconds", (  _wfin_leon - _wdebut_leon) );
	getInfo()->add(1, "Speed", "%.2f Mo/seconds", (System::file().getSize(_outputFilename)/1000000.0) / (  _wfin_leon - _wdebut_leon) );
	
	
	//Test decompressed file against original reads file (decompressed and original read file must be in the same dir)
	if(getParser()->saw (Leon::STR_TEST_DECOMPRESSED_FILE)){
		
		getInfo()->add(1, "Checking decompressed file");
		
		string dir = System::file().getDirectory(_inputFilename);

		string prefix = System::file().getBaseName(_inputFilename);;
		//while(prefix.find('.') != string::npos){
		//	int lastindex = prefix.find_last_of(".");
		//	prefix = prefix.substr(0, lastindex);
		//}
		//string prefix = System::file().getBaseName(_inputFilename);
		
		string originalFilename;
		IBank* originalBank;
		IBank* newBank;
		Iterator<Sequence>* originalBankIt;
		Iterator<Sequence>* newBankIt;
		
		if(_isFasta)
			originalFilename = dir + "/" + prefix + ".fasta";
		else
			originalFilename = dir + "/" + prefix + ".fastq";
		
		
		getInfo()->add(2, "Original file", "%s", originalFilename.c_str());
		getInfo()->add(2, "New file", "%s",  _outputFile->getPath().c_str());
		
		originalBank = Bank::open(originalFilename);
		originalBankIt = originalBank->iterator();
		originalBankIt->first();
		newBank = Bank::open(_outputFile->getPath());
		newBankIt = newBank->iterator();
		newBankIt->first();
		
		//int i=0;
		
		while(true){
			if(newBankIt->isDone()){
				if(originalBankIt->isDone())
					getInfo()->add(1, "OK");
				else
					getInfo()->add(1, "Decompressed file end but not the original file");
				break;
			}
			if(originalBankIt->isDone()){
				if(newBankIt->isDone())
					getInfo()->add(1, "OK");
				else
					getInfo()->add(1, "Original file end but not the decomrpessed file");
				break;
			}
			
			string originalHeader = (*originalBankIt)->getComment();
			string originalDna = (string((*originalBankIt)->getDataBuffer())).substr(0, (*originalBankIt)->getDataSize());
			
			
			string newHeader = (*newBankIt)->getComment();
			string newDna = (string((*newBankIt)->getDataBuffer())).substr(0, (*newBankIt)->getDataSize());
			
			if(originalHeader != newHeader){
				getInfo()->add(1, "Sequence with a different header", "%i", (*newBankIt)->getIndex());
                getInfo()->add(2, "original", "%s", originalHeader.c_str());
                getInfo()->add(2, "new", "%s", newHeader.c_str());
				break;
			}
			
			if(originalDna != newDna){
                getInfo()->add(1, "Sequence with a different DNA", "%i", (*newBankIt)->getIndex());
                getInfo()->add(2, "original", "%s", originalDna.c_str());
                getInfo()->add(2, "new", "%s", newDna.c_str());
				break;
			}
			
			originalBankIt->next();
			newBankIt->next();
			
			//i ++;
			//cout << i << endl;
			//if(i > 20) return;
		}
	}
}



Leon::LeonIterator::LeonIterator( Leon& refl)
: _leon(refl), _isDone(true) , _isInitialized(false)
{
	_stream_qual = _stream_header = _stream_dna = NULL ;

}

void Leon::LeonIterator::first()
{
	//printf("iter first\n");
	init  ();

	_idxB=0;
	_livingThreadCount=0;
	_currentTID=0;
	_isDone = false;
	_readingThreadBlock = false;
	_readid=0;
	if(_stream_qual!= NULL) delete  _stream_qual;
	if(_stream_header!= NULL) delete  _stream_header;
	if(_stream_dna!= NULL) delete  _stream_dna;
	
	next();

}

void Leon::LeonIterator::next()
{
//	printf("---------- iter next ------------\n");

	if(_livingThreadCount==0 ||
	   (( _currentTID>= _livingThreadCount) && !_readingThreadBlock )
	   )
	{
		readNextBlocks();
		if(_isDone) return;
	}
	
	if(!_readingThreadBlock)
	{
		//  assert (_currentTID < _livingThreadCount)
		readNextThreadBock();
	}
	
	assert(_readingThreadBlock);

	//try to get nex tseq from this thread block
	
	std::string line;
	stringstream sint;
	sint << _readid;
	std::string current_comment;
	std::string current_dna;
	std::string current_qual;
	
	
	if( ! _leon._noHeader)
	{
		if(getline(*_stream_header, line)){
			current_comment +=  line ;
		}
		else
			_readingThreadBlock = false;
	}
	else
	{
		current_comment += sint.str() ;
		
		_readid++;
	}
	
	
	
	if(getline(*_stream_dna, line)){
		current_dna +=  line ;
	}
	else
		_readingThreadBlock = false;
	
	
	if( ! _leon._isFasta)
	{
		if(getline(*_stream_qual, line)){
			current_qual +=  line ;
		}
		else
			_readingThreadBlock = false;
	}
	
	if(_readingThreadBlock)
	{
		Sequence *currentSeq = _item;
		currentSeq->setComment(current_comment);
		currentSeq->setQuality(current_qual);
		
		//currData.set ((char *)current_dna.c_str(), current_dna.size()  );
		
		// huum casting const char * to char *; not nice, could be fixed with strdup but want to avoid unnecessary copy,
		//the set() method *should* take a const anyway
		currentSeq->getData().set((char *)current_dna.c_str(), current_dna.size()  );
	}
	else  //reached end of current thread block, try to advance to next block
	{
		next();
	
	}

	
}


void Leon::LeonIterator::readNextBlocks()
{
//	printf("--- iter readNextBlocks  %i / %i ---\n",_idxB,_leon._dnaBlockSizes.size());

	if(_idxB >= _leon._dnaBlockSizes.size()){
		_isDone= true;
	}
	if(!_isDone)
	{
		_leon.decompressionDecodeBlocks(_idxB,_livingThreadCount); //this will update _idxB and _livingThreadCount
		_currentTID =0;
	}
	
//	printf("___ done iter readNextBlocks  %i / %i ---\n",_idxB,_leon._dnaBlockSizes.size());

	
}


//put next chunk in stream_qual,stream_header and _stream_dna
void Leon::LeonIterator::readNextThreadBock()
{
//	printf("--- iter readNextThreadBock  %i %i ---\n",_currentTID,_livingThreadCount);

	pthread_join(_leon._tab_threads[_currentTID], NULL);
	
 	_hdecoder = NULL;
	_qdecoder = NULL;
	_ddecoder = _leon._dnadecoders[_currentTID];
	
	
	if(_stream_qual!= NULL) delete  _stream_qual;
	if(_stream_header!= NULL) delete  _stream_header;
	if(_stream_dna!= NULL) delete  _stream_dna;
	
	_stream_qual = NULL;
	_stream_header = NULL;
	_stream_dna = NULL;
	
	if(! _leon._isFasta)
	{
		_qdecoder = _leon._qualdecoders[_currentTID];
		_stream_qual = new std::istringstream (_qdecoder->_buffer);
		_qdecoder->_buffer.clear();
	}
	
	if(! _leon._noHeader)
	{
		_hdecoder = _leon._headerdecoders[_currentTID];
		_stream_header = new std::istringstream (_hdecoder->_buffer);
		_hdecoder->_buffer.clear();
	}
	
	_stream_dna = new std::istringstream (_ddecoder->_buffer);
	_ddecoder->_buffer.clear();
	
	
	//std::string output_buff;
	//output_buff.reserve(READ_PER_BLOCK * 500);
	
	_readingThreadBlock = true;
	_currentTID++;

	
//	printf("___ done iter readNextThreadBock  %i %i ---\n",_currentTID,_livingThreadCount);


///	u_int64_t readid=0;
}

Leon::LeonIterator::~LeonIterator ()
{
	_leon.decoders_cleanup();
}


void Leon::LeonIterator::estimate(u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)
{
	readDataset(_leon._subgroupInfo,"readcount",number);
	readDataset(_leon._subgroupInfo,"totalDnaSize",totalSize);
	
	int maxsizei ;

	readDataset(_leon._subgroupInfo,"maxSequenceSize",maxsizei);
	maxSize = maxsizei;

}

void Leon::LeonIterator::init()
{
	if (_isInitialized == true)  { return ;}

	///printf("iter init\n");

	_leon.startDecompression_setup();
	_leon.decoders_setup();
	
	
	_isInitialized = true;
}


void Leon::LeonIterator::finalize()
{
	
	
}

//////////////////////////////////////////////////
//////////////////// BankLeon ////////////////////
//////////////////////////////////////////////////

BankLeon::BankLeon (const std::string& filename)
{
	_fname = filename;
	_leon = NULL;
	
	_leon = new Leon();

	
	std::vector<std::string> arguments;
	arguments.push_back("leon");
	arguments.push_back("-file");
	arguments.push_back(_fname);
	arguments.push_back("-d");
	arguments.push_back(Leon::STR_INIT_ITER);
	
	std::vector<char*> argv;
	for (const auto& arg : arguments)
		argv.push_back((char*)arg.data());
	argv.push_back(nullptr);
	int argc = argv.size() - 1;

	_leon->run (argc, argv.data());

}

BankLeon::~BankLeon ()
{
	if(_leon!=NULL)
		delete _leon;
	
}

u_int64_t BankLeon::getSize ()
{
	return System::file().getSize (_fname);
}

int64_t BankLeon::getNbItems () {
	u_int64_t number;
	readDataset(_leon->_subgroupInfo,"readcount",number);
	return number;
}



void BankLeon::estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)
{
	readDataset(_leon->_subgroupInfo,"readcount",number);
	readDataset(_leon->_subgroupInfo,"totalDnaSize",totalSize);
	
	int maxsizei ;
	
	readDataset(_leon->_subgroupInfo,"maxSequenceSize",maxsizei);
	maxSize = maxsizei;
	
}




/// BankLeonFactory : test hdf5  storage opening, if leon version can be found, this is a leon file

IBank* BankLeonFactory::createBank (const std::string& uri)
{
	//printf("create bank factory  Leon  %s \n",uri.c_str());
	bool isLEON = false;

	try {
		isLEON = true;
		
		auto storageH5file = StorageFactory(STORAGE_HDF5).create (uri,false,false,true); //open without adding extension h5

		
		//auto groupLeon = new tools::storage::impl::Group((*storageH5file)().getGroup   ("leon"));
		
		auto _subgroupInfo = new tools::storage::impl::Group((*storageH5file)().getGroup   ("metadata"));

		auto _subgroupInfoCollection = & _subgroupInfo->getCollection<math::NativeInt8> ("infobyte");

		std::string  leonversion =  _subgroupInfoCollection->getProperty ("version");
	//	std::cout << "leon file version : "<< leonversion << std::endl;

	} catch (system::Exception& e) {
//		 std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
		isLEON = false;
	}
	
	return (isLEON ? new BankLeon (uri) : NULL);
	
}

