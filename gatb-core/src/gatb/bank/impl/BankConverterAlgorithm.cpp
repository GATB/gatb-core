/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
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

#include <gatb/bank/impl/BankConverterAlgorithm.hpp>

#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Banks.hpp>

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

static const char* progressFormat1 = "Bank: fasta to binary                  ";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankConverterAlgorithm::BankConverterAlgorithm (IBank* bank,  size_t kmerSize, const std::string& outputUri)
    : Algorithm ("bankconverter"), _kind(BANK_CONVERT_TMP), _bankInput(0), _bankOutput(0), _outputUri(outputUri)
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
BankConverterAlgorithm::BankConverterAlgorithm (tools::storage::impl::Storage& storage)
: Algorithm ("bankconverter"), _kind(BANK_CONVERT_NONE), _bankInput(0), _bankOutput(0)
{
    string xmlString = storage(this->getName()).getProperty ("xml");
    stringstream ss; ss << xmlString;   getInfo()->readXML (ss);
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
    /** We may have no conversion at all to do. */
    if (_kind == BANK_CONVERT_NONE)
    {
        setBankOutput (_bankInput);
        return;
    }

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

     /** We use a block for measuring the time elapsed in it. */
     {
         TIME_INFO (getTimeInfo(), "conversion");

         /** We iterate the sequences of the input bank. */
         for (itBank->first(); !itBank->isDone(); itBank->next())
         {
             nbSeq ++;
             sizeSeq += (*itBank)->getDataSize();

             /** We insert the current sequence into the output bank. */
             _bankOutput->insert (itBank->item());
         }
     }

     /** We flush the output bank (important if it is a file output). */
     _bankOutput->flush ();

     /** We gather some statistics. */
     getInfo()->add (1, "info");
     getInfo()->add (2, "input",            "%s",   _bankInput->getId().c_str());
     getInfo()->add (2, "sequences_number", "%ld",  nbSeq);
     getInfo()->add (2, "sequences_size",   "%ld",  sizeSeq);
     getInfo()->add (2, "output_size",      "%ld",  _bankOutput->getSize());
     getInfo()->add (2, "ratio",            "%.3f",  (double)sizeSeq / (double)_bankOutput->getSize());
     getInfo()->add (1, getTimeInfo().getProperties("time"));
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
