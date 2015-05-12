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

#ifndef _COUNT_PROCESSOR_SOLIDITY_HPP_
#define _COUNT_PROCESSOR_SOLIDITY_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/CountProcessorAbstract.hpp>
#include <gatb/tools/misc/api/Enums.hpp>
#include <gatb/tools/misc/api/Range.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

template<size_t span, class Derived>
class CountProcessorSolidityAbstract : public CountProcessorAbstract<span>
{
public:

    CountProcessorSolidityAbstract (const tools::misc::CountRange& threshold) : _threshold(threshold), _total(0), _ok(0)   {}

    virtual ~CountProcessorSolidityAbstract()  {}

    /** */
    bool process (size_t partId, const typename Kmer<span>::Type& kmer, const CountVector& count, CountNumber sum)
    {
        /** We use static polymorphism here. */
        bool result = static_cast<Derived*>(this)->check (count, sum);

        _total ++;
        if (result)  { _ok++; }
        return result;
    }

    void endPart (size_t passId, size_t partId)
    {
        CountProcessorSolidityAbstract* proto = dynamic_cast<CountProcessorSolidityAbstract*> (this->getPrototype());
        if (proto != 0)
        {
            __sync_fetch_and_add (&proto->_total, this->_total);
            __sync_fetch_and_add (&proto->_ok,    this->_ok);
        }
    }

    tools::misc::impl::Properties getProperties() const
    {
        tools::misc::impl::Properties result;
        result.add (0, "kmers");
        result.add (1, "kmers_nb_distinct",  "%ld", _total);
        result.add (1, "kmers_nb_solid",     "%ld", _ok);
        result.add (1, "kmers_nb_weak",      "%ld", _total - _ok);
        if (_total > 0)  {  result.add (1, "kmers_percent_weak", "%.1f", 100.0 - 100.0 * (double)_ok / (double)_total  );  }

        return result;
    }

protected:

    /** */
    CountProcessorAbstract<span>* doClone ()  { return new Derived (_threshold); }

    tools::misc::CountRange _threshold;

    u_int64_t _total;
    u_int64_t _ok;
};

/********************************************************************************/

template<size_t span=KMER_DEFAULT_SPAN>
class CountProcessorSoliditySum : public CountProcessorSolidityAbstract<span,CountProcessorSoliditySum<span> >
{
public:

    CountProcessorSoliditySum (const tools::misc::CountRange& threshold)
        : CountProcessorSolidityAbstract<span,CountProcessorSoliditySum<span> > (threshold)  {}

    bool check (const CountVector& count, CountNumber sum)
    {
        return this->_threshold.includes (sum);
    }

    std::string getName() const  { return std::string("sum"); }
};

/********************************************************************************/

template<size_t span=KMER_DEFAULT_SPAN>
class CountProcessorSolidityMax : public CountProcessorSolidityAbstract<span,CountProcessorSolidityMax<span> >
{
public:

    CountProcessorSolidityMax (const tools::misc::CountRange& threshold)
        : CountProcessorSolidityAbstract<span,CountProcessorSolidityMax<span> > (threshold)  {}

    bool check (const CountVector& count, CountNumber sum)
    {
        return this->_threshold.includes (*std::max_element (count.begin(),count.end()));
    }

    std::string getName() const  { return std::string("max"); }
};

/********************************************************************************/

template<size_t span=KMER_DEFAULT_SPAN>
class CountProcessorSolidityMin : public CountProcessorSolidityAbstract<span,CountProcessorSolidityMin<span> >
{
public:

    CountProcessorSolidityMin (const tools::misc::CountRange& threshold)
        : CountProcessorSolidityAbstract<span,CountProcessorSolidityMin<span> > (threshold)  {}

    bool check (const CountVector& count, CountNumber sum)
    {
        return this->_threshold.includes (*std::min_element (count.begin(),count.end()));
    }

    std::string getName() const  { return std::string("min"); }
};

/********************************************************************************/

template<size_t span=KMER_DEFAULT_SPAN>
class CountProcessorSolidityAll : public CountProcessorSolidityAbstract<span,CountProcessorSolidityAll<span> >
{
public:

    CountProcessorSolidityAll (const tools::misc::CountRange& threshold)
        : CountProcessorSolidityAbstract<span,CountProcessorSolidityAll<span> > (threshold)  {}

    bool check (const CountVector& count, CountNumber sum)
    {
        for (size_t i=0; i<count.size(); i++)  {  if (this->_threshold.includes(count[i]) == false)   { return false; }  }
        return true;
    }

    std::string getName() const  { return std::string("all"); }
};

/********************************************************************************/

template<size_t span=KMER_DEFAULT_SPAN>
class CountProcessorSolidityOne : public CountProcessorSolidityAbstract<span,CountProcessorSolidityOne<span> >
{
public:

    CountProcessorSolidityOne (const tools::misc::CountRange& threshold)
        : CountProcessorSolidityAbstract<span,CountProcessorSolidityOne<span> > (threshold)  {}

    bool check (const CountVector& count, CountNumber sum)
    {
        for (size_t i=0; i<count.size(); i++)  {  if (this->_threshold.includes(count[i]) == true)   { return true; }  }
        return false;
    }

    std::string getName() const  { return std::string("one"); }
};

/********************************************************************************/

template<size_t span>
class CountProcessorSolidityFactory
{
public:

    static ICountProcessor<span>* create (tools::misc::KmerSolidityKind kind, const tools::misc::CountRange& threshold)
    {
        switch (kind)
        {
        case tools::misc::KMER_SOLIDITY_MIN: return new CountProcessorSolidityMin<span> (threshold);
        case tools::misc::KMER_SOLIDITY_MAX: return new CountProcessorSolidityMax<span> (threshold);
        case tools::misc::KMER_SOLIDITY_ONE: return new CountProcessorSolidityOne<span> (threshold);
        case tools::misc::KMER_SOLIDITY_ALL: return new CountProcessorSolidityAll<span> (threshold);
        case tools::misc::KMER_SOLIDITY_SUM: return new CountProcessorSoliditySum<span> (threshold);
        default:  throw system::Exception ("unable to create CountProcessorSolidity instance for kind %d", kind);
        }
    }

    static ICountProcessor<span>* create (tools::misc::IProperties& props, const tools::misc::CountRange& threshold)
    {
        tools::misc::KmerSolidityKind kind;
        parse (props.getStr (STR_SOLIDITY_KIND), kind);
        return create (kind, threshold);
    }

};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _COUNT_PROCESSOR_SOLIDITY_HPP_ */
