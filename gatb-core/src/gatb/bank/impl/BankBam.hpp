/*****************************************************************************
 *
 *   Copyright (C) 2018 ENANCIO
 *****************************************************************************/



#ifndef _GATB_CORE_BANK_IMPL_BANK_BAM_HPP_
#define _GATB_CORE_BANK_IMPL_BANK_BAM_HPP_

/********************************************************************************/
#include <zlib.h>
#include <gatb/bank/impl/AbstractBank.hpp>
#include <vector>
#include <string>


#include "sam.h"


/********************************************************************************/
namespace gatb      {
/** \brief Core package of the GATP project.
 *
 * The gatb::core package holds all the fundamental packages needed for writting
 * assembly algorithms.
 *
 * It holds some generic tools, like operating system abstraction, collections management or design patterns
 * concerns. It also holds recurrent needs as reading genomic banks, handling kmers and so on.
 */
namespace core      {
/** \brief Package for genomic databases management. */
namespace bank      {
/** \brief Implementation for genomic databases management. */
namespace impl      {
/********************************************************************************/


class BankBam : public AbstractBank
{
public:

	/** Returns the name of the bank format. */
	static const char* name()  { return "bam"; }
	
	/** Constructor.
	 * \param[in] nbSequences : number of sequences of the random bank
	 * \param[in] length : length of a sequence. */
	BankBam (const std::string& filename);
	
	/** Destructor. */
	~BankBam ();
	
	/** \copydoc IBank::getId. */
	std::string getId ()  { return _fname; }

	/** \copydoc IBank::iterator */
	tools::dp::Iterator<Sequence>* iterator ()  { return new BankBam::BamIterator (*this);}

	/** */
	int64_t getNbItems () {return -1;} ;

	/** \copydoc IBank::insert */
	void insert (const Sequence& item) {}

	/** \copydoc IBank::flush */
	void flush ()  {}

	/** \copydoc IBank::getSize */
	u_int64_t getSize () ;

	/** \copydoc IBank::estimate */
	void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize) ;

	/** \return maximum number of files. */
	static const size_t getMaxNbFiles ()  { return 0; }

	/************************************************************/
	
	class BamIterator : public tools::dp::Iterator<Sequence>
	{
	public:
		
		BamIterator (BankBam& ref);
		
		/** Destructor */
		~BamIterator ();
		
		/** \copydoc tools::dp::Iterator::first */
		void first();
		
		/** \copydoc tools::dp::Iterator::next */
		void next();
		
		/** \copydoc tools::dp::Iterator::isDone */
		bool isDone ()  { return _isDone; }
		
		/** \copydoc tools::dp::Iterator::item */
		Sequence& item ()     { return *_item; }
		
		/** Estimation of the sequences information */
		void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize) ; //todo
		
	private:
		
		/** Reference to the underlying lena instance. */
		BankBam&    _BankBam;
		
		/** Tells whether the iteration is finished or not. */
		bool _isDone;
		
		/** Tells whether the instance is initialized. */
		bool _isInitialized;
		
		
		/** Initialization method. */
		void init ();
		
		/** Finish method. */
		void finalize ();
		

		samFile * _samfile;
		bam_hdr_t * _bamHdr;
		bam1_t * _bam_record;
	};
	

protected:

	std::string _fname;
};

/********************************************************************************/

/* \brief Factory for the BankFasta class. */
class BankBamFactory : public IBankFactory
{
public:

    /** \copydoc IBankFactory::createBank */
    IBank* createBank (const std::string& uri);
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_BANK_BAM_HPP_ */
