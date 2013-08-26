/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/bank/impl/BankConverterAlgorithm.hpp>

#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>

#include <gatb/kmer/impl/Model.hpp>

// We use the required packages
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

#define DEBUG(a)  //printf a

/********************************************************************************/
namespace gatb {  namespace core { namespace bank { namespace impl {
/********************************************************************************/

static const char* progressFormat1 = "Bank: fasta to binary                ";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankConverterAlgorithm::BankConverterAlgorithm (IBank* bank,  size_t kmerSize, const std::string& outputUri)
    : Algorithm ("bankconverter"), _bankInput(0), _bankOutput(0), _outputUri(outputUri)
{
    setBankInput  (bank);
    setBankOutput (new BankBinary (outputUri, kmerSize));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankConverterAlgorithm::~BankConverterAlgorithm ()
{
    setBankInput  (0);
    setBankOutput (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankConverterAlgorithm::execute ()
{
    TIME_INFO (getTimeInfo(), "bank conversion");

    /** We may have to delete the binary file if it already exists. */
    System::file().remove (_outputUri);

    /** We get information about the FASTA bank. */
    u_int64_t number, totalSize, maxSize;
    _bankInput->estimate (number, totalSize, maxSize);

    /** We need an iterator on the input bank. */
     Iterator<Sequence>* itBank = createIterator<Sequence> (
         _bankInput->iterator(),
         number,
         progressFormat1
     );
     LOCAL (itBank);

     u_int64_t   nbSeq = 0;
     u_int64_t sizeSeq = 0;

     /** We iterate the sequences of the input bank. */
     for (itBank->first(); !itBank->isDone(); itBank->next())
     {
         nbSeq ++;
         sizeSeq += (*itBank)->getDataSize();

         /** We insert the current sequence into the output bank. */
         _bankOutput->insert (itBank->item());
     }

     /** We flush the output bank (important if it is a file output). */
     _bankOutput->flush ();

     /** We gather some statistics. */
     getInfo()->add (1, "info");
     getInfo()->add (2, "sequences number", "%ld",  nbSeq);
     getInfo()->add (2, "sequences size",   "%ld",  sizeSeq);
     getInfo()->add (2, "output size",      "%ld",  _bankOutput->getSize());
     getInfo()->add (2, "ratio",            "%.3f",  (double)sizeSeq / (double)_bankOutput->getSize());
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
