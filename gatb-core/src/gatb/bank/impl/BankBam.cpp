/*****************************************************************************
 *
 *   Copyright (C) 2018 ENANCIO
*****************************************************************************/

#include <gatb/bank/impl/BankBam.hpp>
#include <gatb/bank/impl/BankComposite.hpp>

#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/api/StringsRepository.hpp>
#include <gatb/tools/misc/impl/Tokenizer.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <algorithm>
#include <string.h>
#include <errno.h>
#include <zlib.h>
using namespace std;
using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;
using namespace gatb::core::tools::misc;


/********************************************************************************/
namespace gatb {  namespace core {  namespace bank {  namespace impl {
/********************************************************************************/

	
	
	BankBam::BankBam (const std::string& filename)
	{
		_fname = filename;
	}
	
	
	BankBam::~BankBam ()
	{

	}
	
	
	u_int64_t BankBam::getSize () {
		return System::file().getSize (_fname) ;
	}
	
	void BankBam::estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize) {
		/** We create an iterator for the bank. */
		BankBam::BamIterator it (*this);
		
		/** We delegate the computation to the iterator. */
		it.estimate (number, totalSize, maxSize);
		
	}
	
	
	///////////// iterator
	
	
	BankBam::BamIterator::BamIterator( BankBam& refb)
	: _BankBam(refb), _isDone(true) , _isInitialized(false), _samfile(NULL), _bamHdr(NULL), _bam_record(NULL)
	{
		;
	}
	
	
	BankBam::BamIterator::~BamIterator ()
	{

		if(_samfile!= NULL)
		{
			hts_close(_samfile);
			_samfile = NULL;
		}
		
		if(_bam_record!=NULL)
		{
			bam_destroy1(_bam_record);
			_bam_record=NULL;
		}
	
		if(_bamHdr!=NULL)
		{
			bam_hdr_destroy(_bamHdr);
			_bamHdr =NULL;
		}
		
		
	}

	
	
	void BankBam::BamIterator::init()
	{
		//if (_isInitialized == true)  { return ;}
		
		if(_samfile!= NULL)
		{
			hts_close(_samfile);
			_samfile = NULL;
		}
		
		 _samfile =  hts_open(_BankBam._fname.c_str(),"r"); //open bam/cram/sam file

		
		if(_bam_record!=NULL)
		{
			bam_destroy1(_bam_record);
			_bam_record=NULL;
		}
		
		_bam_record = bam_init1(); //initialize a record

		
		if(_bamHdr!=NULL)
		{
			bam_hdr_destroy(_bamHdr);
			_bamHdr =NULL;
		}
		

		
		_bamHdr = sam_hdr_read(_samfile); //read header

		
		_isInitialized = true;
	}

	
	void BankBam::BamIterator::first()
	{
		//how to rewind bam file ?

		init  ();
		
		
		
		_isDone = false;
		
		next();
		
	}
	
	
	void BankBam::BamIterator::next()
	{

		
		int rsam = sam_read1(_samfile, _bamHdr, _bam_record);
		
		if(rsam >=0)
		{
			Sequence *currentSeq = _item;

		//	printf(". . . seq len %i   \n",_bam_record->core.l_qseq);
			//get dna sequence
			uint8_t * pseq = bam_get_seq(_bam_record);
			std::string outseq(_bam_record->core.l_qseq, 'N');
			for (int32_t ii = 0; ii < _bam_record->core.l_qseq; ++ii)
				outseq[ii] = seq_nt16_str[bam_seqi(pseq,ii)];
			
			
			// huum casting const char * to char *; not nice, could be fixed with strdup but want to avoid unnecessary copy,
			//the set() method *should* take a const anyway
			currentSeq->getData().set((char *)outseq.c_str(), outseq.size()  );
			
			
			//get quality sequence
			uint8_t * pqual = bam_get_qual(_bam_record);
			std::string currentQual;
			int offset = 33;
			if(pqual != NULL)
			{
				currentQual = std::string(_bam_record->core.l_qseq, ' ');
				for (int32_t ii = 0; ii < _bam_record->core.l_qseq; ++ii)
					currentQual[ii] = (char)(pqual[ii] + offset);
				
				currentSeq->setQuality(currentQual);
			}
			
			
			//get header
			std::string current_comment = std::string (bam_get_qname(_bam_record));
			currentSeq->setComment(current_comment);
			
		}
		else
		{
			_isDone= true;
			return;
		}
		
	}
	
	
	
	void BankBam::BamIterator::finalize()
	{
		
		bam_destroy1(_bam_record); _bam_record= NULL;
		sam_close(_samfile); _samfile = NULL;
		bam_hdr_destroy(_bamHdr); _bamHdr = NULL;

	}

	
	void BankBam::BamIterator::estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)
	{
		
		///returning *very* rough approximation
		number = _BankBam.getSize()/1024ULL/1024ULL * 13000;
		
		number = std::max(1000ULL ,(long long unsigned int ) number );
		
		//printf("rough estimation %llu \n",number);
		//unknown
		totalSize=number *100; //very rough too
		maxSize=0;
	}
	
	
	IBank* BankBamFactory::createBank (const std::string& uri)
	{
		
		hts_set_log_level (HTS_LOG_OFF);
		bool isBAM = false;
		
		//try to open the bam file for reading
		
		samFile * testsamfile =  hts_open(uri.c_str(),"r"); //open bam/cram/sam file
		if(testsamfile == NULL)
			return NULL;
		
		bam_hdr_t * bamHdr = sam_hdr_read(testsamfile); //read header
		if (bamHdr == NULL) {
			return NULL;
		}
		
		
		bam1_t * bam_record  = bam_init1(); //initialize a record

		//try to read first record
		int rsam = sam_read1(testsamfile, bamHdr, bam_record);
		if(rsam >=0)
			isBAM = true;
		else
			isBAM = false;
		
		bam_destroy1(bam_record);
		
		
		sam_close(testsamfile);
		bam_hdr_destroy(bamHdr);
		
		return (isBAM ? new BankBam (uri) : NULL);
	}

	
	
/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
