/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/bank/impl/BankRandom.hpp>
#include <gatb/tools/misc/api/StringsRepository.hpp>

#include <gatb/system/impl/System.hpp>

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;
using namespace gatb::core::tools::misc;

#define DEBUG(a)  //printf a

/********************************************************************************/
namespace gatb {  namespace core {  namespace bank {  namespace impl {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankRandom::BankRandom (size_t nbSequences, size_t length)
    : _nbSequences(nbSequences), _length(length)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankRandom::~BankRandom ()
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankRandom::estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankRandom::Iterator::Iterator(const BankRandom& bank)
    : _bank(bank),_rank(0), _isDone(true), _dataRef(0)
{
    setDataRef (new Data (bank._length));

    _item->getData().setRef (_dataRef, 0, bank._length);

    srand (time(NULL));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankRandom::Iterator::~Iterator()
{
    setDataRef (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankRandom::Iterator::first()
{
    _rank = -1;
    next ();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankRandom::Iterator::next()
{
    _isDone = (++_rank >= _bank._nbSequences);
    if (!_isDone)
    {
        static char table[] = {'A', 'C', 'T', 'G' };

        char* buffer = _item->getDataBuffer();

        for (size_t i=0; i<_item->getDataSize(); i++)
        {
            buffer [i] = table[rand() % sizeof(table)/sizeof(table[0])];
        }
    }
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
