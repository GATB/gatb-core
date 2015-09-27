/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  R.Chikhi, G.Rizk, E.Drezen
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

#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/BankKmerIterator.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/math/LargeInt.hpp>
#include <gatb/tools/math/Integer.hpp>
#include <iostream>
#include <string.h>

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
using namespace gatb::core::tools::math;
using namespace gatb::core::tools::misc::impl;

/********************************************************************************/

struct Hash  {   kmer_type operator () (kmer_type& lkmer)
{
    kmer_type kmer_hash;

    kmer_hash = lkmer ^ (lkmer >> 14);
    kmer_hash = (~kmer_hash) + (kmer_hash << 18);
    kmer_hash = kmer_hash ^ (kmer_hash >> 31);
    kmer_hash = kmer_hash * 21;
    kmer_hash = kmer_hash ^ (kmer_hash >> 11);
    kmer_hash = kmer_hash + (kmer_hash << 6);
    kmer_hash = kmer_hash ^ (kmer_hash >> 22);

    return kmer_hash;
}};

Hash theHash;

/********************************************************************************/

struct FunctorIter1
{
    void operator() (Sequence& sequence)
    {
        nbSeq++;
        dataSeq += sequence.getData().size();

        model->build (sequence.getData(), kmers);

        size_t size = kmers.size();

        for (size_t i=0; i<size; i++)
        {
            kmer_type h = theHash(kmers[i]);

            if ((h % nbPass) != pass)  { continue; }

            nbKmers ++;
            checksumKmers = checksumKmers + h;
        }
    }

    FunctorIter1 () : nbSeq(0), dataSeq(0), nbKmers(0), checksumKmers(0), pass(0), nbPass(0), model(0) {}

    vector<kmer_type> kmers;

    u_int64_t nbSeq;
    u_int64_t dataSeq;
    u_int64_t nbKmers;
    kmer_type checksumKmers;

    size_t pass;
    size_t nbPass;
    KmerModel* model;
};

template <typename Functor> void iter1 (
    IDispatcher& dispatcher,
    IBank& bank,
    KmerModel& model,
    Functor& hash,
    IteratorListener* progress=0
)
{
    // We use the provided listener if any
    LOCAL (progress);

    // We need an iterator on the FASTA bank.
    Iterator<Sequence>* itBank = bank.iterator();
    LOCAL (itBank);

    // We declare two kmer iterators for the two banks and a paired one that links them.
    KmerModel::Iterator itKmer (model);

    // We create an iterator over the paired iterator on sequences
    SubjectIterator<Sequence> itSeq (itBank, 5*1000);

    if (progress)  {  itSeq.addObserver (progress);  }

    u_int64_t total_nbKmers       = 0;
    kmer_type total_checksumKmers = 0;

    size_t nbPasses = 10;

    // We get some information about the kmers.
    u_int64_t nbSeq[nbPasses];
    u_int64_t dataSeq[nbPasses];
    u_int64_t nbKmers[nbPasses];
    kmer_type checksumKmers[nbPasses];

    for (size_t p=0; p<nbPasses; p++)
    {
        // We get current time stamp
        ITime::Value t0 = System::time().getTimeStamp();

        vector<FunctorIter1> functors (dispatcher.getExecutionUnitsNumber());

        for (size_t i=0; i<functors.size(); i++)
        {
            functors[i].pass   = p;
            functors[i].nbPass = nbPasses;
            functors[i].model  = &model;
        }

        nbSeq[p]         = 0;
        dataSeq[p]       = 0;
        nbKmers[p]       = 0;
        checksumKmers[p] = 0;

#if 0
        dispatcher.iterate (itSeq, functors);
#endif

        for (size_t i=0; i<functors.size(); i++)
        {
            nbSeq[p]         += functors[i].nbSeq;
            dataSeq[p]       += functors[i].dataSeq;
            nbKmers[p]       += functors[i].nbKmers;
            checksumKmers[p] = checksumKmers[p] + functors[i].checksumKmers;

            total_nbKmers       += functors[i].nbKmers;
            total_checksumKmers = total_checksumKmers + functors[i].checksumKmers;
        }

        // We get current time stamp
        ITime::Value t1 = System::time().getTimeStamp();

        // We dump some information about the iterated kmers;
        cout << "FOUND " << nbKmers[p] << " kmers  for " << nbSeq[p] << " sequences  (" << dataSeq[p] << " bytes)  "
            << "in " << (t1-t0) << " msec (rate " << (double)nbKmers[p] / (double) (t1>t0 ? t1-t0 : 1) << " kmers/msec),  "
            << "checksum is " << checksumKmers[p]
            << endl;
    }
    cout << endl;
    cout << "TOTAL KMERS " << total_nbKmers << "  WITH CHECKSUM " << hex << total_checksumKmers << "  WITH " << Integer::getName()  <<  endl;

    for (size_t p=0; p<nbPasses; p++)
    {
        cout << "   [" << p << "]  " << 100.0 * (double)nbKmers[p] / (double) total_nbKmers << endl;
    }
}

/********************************************************************************/

template <typename Functor> void iter2 (IBank& bank, KmerModel& model, Functor& hash, IteratorListener* progress=0)
{
#if 1
    // We use the provided listener if any
    LOCAL (progress);

    BankKmerIterator itKmerBank (bank, model);

    // We get some information about the kmers.
    u_int64_t nbKmers       = 0;
    kmer_type checksumKmers = 0;

    // We get current time stamp
    ITime::Value t0 = System::time().getTimeStamp();

    if (progress)  {  itKmerBank.addObserver (progress);  }

    for (itKmerBank.first(); !itKmerBank.isDone();  itKmerBank.next())
    {
        nbKmers      += 1;
        checksumKmers = checksumKmers + hash (*itKmerBank);
    }

    // We get current time stamp
    ITime::Value t1 = System::time().getTimeStamp();

    // We dump some information about the iterated kmers;
    cout << "FOUND " << nbKmers << " kmers "
         << "in " << (t1-t0) << " msec (rate " << (double)nbKmers / (double) (t1-t0) << " kmers/msec),  "
         << "checksum is " << checksumKmers
         << endl;
#endif
}

/********************************************************************************/

struct FunctorIter3
{
    void operator() (std::vector<kmer_type>& kmers)
    {
        Hash hash;

        size_t size = kmers.size();

        for (size_t i=0; i<size; i++)
        {
            nbKmers ++;
            checksumKmers = checksumKmers + hash (kmers[i]);
        }
    }

    u_int64_t nbKmers;
    kmer_type checksumKmers;
};

template <typename Functor> void iter3 (IBank& bank, KmerModel& model, Functor& hash, IteratorListener* progress=0)
{
    // We use the provided listener if any
    LOCAL (progress);

    BankVectorKmerIterator<kmer_type> itKmerBank (bank, model);

    // We get some information about the kmers.
    u_int64_t nbKmers       = 0;
    kmer_type checksumKmers = 0;

    // We get current time stamp
    ITime::Value t0 = System::time().getTimeStamp();

    if (progress)  {  itKmerBank.addObserver (progress);  }

    ParallelCommandDispatcher dispatcher (8);

    vector<FunctorIter3> functors (dispatcher.getExecutionUnitsNumber());

#if 0
    dispatcher.iterate (itKmerBank, functors, 1);
#endif

    for (size_t i=0; i<functors.size(); i++)
    {
        nbKmers       += functors[i].nbKmers;
        checksumKmers = checksumKmers + functors[i].checksumKmers;
    }

    // We get current time stamp
    ITime::Value t1 = System::time().getTimeStamp();

    // We dump some information about the iterated kmers;
    cout << "FOUND " << nbKmers << " kmers "
         << "in " << (t1-t0) << " msec (rate " << (double)nbKmers / (double) (t1-t0) << " kmers/msec),  "
         << "checksum is " << checksumKmers
         << endl;
}

/********************************************************************************/

int main (int argc, char* argv[])
{
    if (argc < 3)
    {
        cerr << "you must provide at least 2 arguments. Arguments are:" << endl;
        cerr << "   1) kmer size"  << endl;
        cerr << "   2) FASTA  bank" << endl;
        cerr << "   3) binary bank" << endl;
        return EXIT_FAILURE;
    }

    // We define the max size of a data line in the FASTA output file
    size_t kmerSize = atoi(argv[1]);

    // We get the URI of the FASTA bank
    string filename (argv[2]);
    string filenameBin = argc <=3 ? (filename + ".bin") : argv[3];

#if 1
    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        ParallelCommandDispatcher dispatcher (8);

        // We declare a kmer model with a given span size.
        KmerModel model (kmerSize);

        // We declare the FASTA bank
        BankFasta bank (filename);

        // We declare a binary bank
        BankBinary bankBin (filenameBin);

        if (System::file().doesExist(filenameBin) == false)
        {
            // We declare some job listener.
            Progress progress (bank.estimateNbSequences(), "FASTA to binary conversion");

            // We convert the FASTA bank in binary format
            IProperties* props = BankHelper::singleton().convert (bank, bankBin, &progress);
            LOCAL (props);
        }

        //HashNull hash;
        Hash hash;

        // TEST 1
        iter1 (dispatcher, bankBin, model, hash, new Progress (bank.estimateNbSequences(), "Iterating 1"));

        // TEST 2
        //iter2 (bankBin, model, hash, new Progress (bank.estimateNbSequences(), "Iterating 2"));

        // TEST 3
        //iter3 (bankBin, model, hash, new Progress (bank.estimateNbSequences(), "Iterating 3"));

        // We remove the binary bank
        //System::file().remove (filename + ".bin");
    }

    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }
#endif
    return EXIT_SUCCESS;
}
