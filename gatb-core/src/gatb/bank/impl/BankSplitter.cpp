/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/bank/impl/BankSplitter.hpp>
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

typedef pair<size_t,size_t> Offset;
vector<Offset> _offsets;

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankSplitter::BankSplitter (
    IBank*   reference,
    size_t   readMeanSize,
    size_t   overlap,
    u_int8_t coverage,
    bool     random
)
    : _reference (0), _readMeanSize(readMeanSize), _coverage(coverage), _overlap(overlap), _random(random)
{
    setReference (reference);

    if (_random)  {  srand (time(NULL));  }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankSplitter::~BankSplitter ()
{
    setReference (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankSplitter::estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)
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
BankSplitter::Iterator::Iterator(const BankSplitter& bank)
    : _dataRef (0), _rank(0), _nbMax(0), _isDone(true), _itRef(0),
      _readMeanSize (bank._readMeanSize), _random(bank._random), _overlap(bank._overlap)
{
    assert (bank._readMeanSize > 0);
    assert (bank._readMeanSize > bank._overlap);

    /** We get the first sequence of the referred bank. */
    setItRef (bank._reference->iterator());
    _itRef->first();
    assert (!_itRef->isDone());

    /** We create the reference Data object (that references the provided string). */
    setDataRef (new Data ((*_itRef)->getDataBuffer()));

    _offsetMax = _dataRef->size() - _readMeanSize;

    DEBUG (("refSize=%d  _readMeanSize=%d  _overlap=%d  _offsetMax=%d\n", _reference->size(), _readMeanSize, _overlap, _offsetMax));

    _offsets.clear();
    size_t idx;
    size_t delta = _readMeanSize -_overlap;
    for (idx=0; idx<_offsetMax; idx+=delta)
    {
        _offsets.push_back (make_pair (idx,_readMeanSize));
    }
    _offsets.push_back (make_pair (idx, _dataRef->size()-idx));

    DEBUG (("FOUND %d offsets\n", _offsets.size()));
    for (size_t i=0; i<_offsets.size(); i++)
    {
        DEBUG (("   %d  %d\n", _offsets[i].first, _offsets[i].second));
    }

    _nbMax = _offsets.size() * bank._coverage;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankSplitter::Iterator::~Iterator()
{
    setDataRef(0);
    setItRef(0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankSplitter::Iterator::first()
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
void BankSplitter::Iterator::next()
{
    _isDone = (++_rank >= _nbMax);
    if (!_isDone)
    {
        //size_t offset = _random ? rand() % _offsetMax : _offsets[_rank % (_offsets.size())];
        size_t offset = _offsets[_rank % (_offsets.size())].first;
        size_t size   = _offsets[_rank % (_offsets.size())].second;

//cout << "rank=" << _rank << "  len=" << _reference->size() << "  offset=" << offset << endl;

        _item->getData().setRef (_dataRef, offset, size);
    }
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
