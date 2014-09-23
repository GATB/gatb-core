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
	
		g_msize =  7;
	
//	typedef Kmer<31>::ModelCanonical  Modelc;
//	typedef typename Kmer<31>::ModelCanonical::Kmer  Mmer;
//	typedef typename Kmer<31>::Type           Typem;
//	

//	freq_mmer = ( u_int64_t * )  calloc (( 1 << (2*g_msize)),sizeof(u_int64_t) ); //cpt on smaller type would be enough
//	stat_mmer = ( u_int64_t * )  calloc (( 1 << (2*g_msize)),sizeof(u_int64_t) );
//	
//	Modelc model (g_msize);
//    vector<Mmer> mmers;
	
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
			 
			 
			 //todo also in this pass :  parse the mmers of the seq itBank->item() to count mmers , to do freq based minim later
//			model.build (itBank->item().getData(), mmers);
//			 for (size_t i=0; i<mmers.size(); i++)
//			 {
//				 Typem cur =mmers[i].value();
//				 freq_mmer[cur.getVal() ] ++;
//			 }
			 
			 
         }
	
//		 max_mmer = 0;
//		 for (int ii= 0; ii< ( 1 <<(g_msize*2)); ii++)
//		 {
//			 if (freq_mmer[ii]> freq_mmer[max_mmer]) max_mmer = ii;
//			 //		typedef typename Kmer<31>::Type           Typem;
//			 //		Typem cur = ii;
//			 //		printf("%s : %lli\n",cur.toString(8).c_str(),stat_mmer[ii]);
//			 
//		 }
		 
//		 Typem mma = max_mmer;
//		 printf("Max mmer by freq is  %llu  (%s)\n",max_mmer, mma.toString(g_msize).c_str());

//		 //debug
//		 for (int ii= 0; ii< ( 1 <<(8*2)); ii++)
//		 {
//			 Typem cur = ii;
//			 printf("%s : %lli\n",cur.toString(8).c_str(),freq_mmer[ii]);
//		 }
		 //
		 
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
