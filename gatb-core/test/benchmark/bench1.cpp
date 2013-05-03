#include <gatb/system/impl/System.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>

#include <gatb/kmer/impl/Model.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/math/LargeInt.hpp>
#include <gatb/tools/math/Integer.hpp>
#include <gatb/tools/math/ttmath/ttmath.h>
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
using namespace gatb::core::tools::misc::impl;
using namespace gatb::core::tools::math;
using namespace gatb::core::tools::misc::impl;

//typedef u_int64_t kmer_type;
//typedef ttmath::UInt<2> kmer_type;
//typedef LargeInt<2> kmer_type;
typedef Integer<2> kmer_type;
//typedef __uint128_t kmer_type;


/********************************************************************************/

struct HashNull  {   kmer_type operator () (kmer_type& lkmer)
{
    return 0;
}};

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

/********************************************************************************/

class BankKmerIterator : public Iterator<kmer_type>, public AbstractSubjectIterator
{
public:

    /** Constructor.
     * \param[in] bank : the bank whose sequences are to be iterated.
     * \param[in] model : kmer model
     */
    BankKmerIterator (IBank& bank, Model<kmer_type>& model)
        : _itSeq(0), _itKmer(model), _isDone(true),  _moduloMask(1), _current(0)
    {
        /** We create an iterator over the sequences of the provided bank.
         * Note that this is a dynamic allocation, so we will have to get rid of the instance
         * in the destructor. */
        setItSeq (bank.iterator());

        /** We set the modulo mask for which we will send notification to potential listeners.
         * Such notification is done on outer loop over Sequence objects. */
        _moduloMask = (1<<10) - 1;
    }

    /** Destructor. */
    ~BankKmerIterator ()
    {
        /** We get rid of the dynamically allocated Sequence iterator. */
        setItSeq(0);
    }

    /** \copydoc Iterator::first */
    void first()
    {
        /** We begin by notifying potential listeners that the iteration is beginning. */
        notifyInit ();

        /** We reset the iteration counter. We will check when this counter is equal to our modulo;
         * in such a case, we will notify our potential listeners.  */
        _current = 0;

        /** We go to the first item of the Sequence iteration. */
        _itSeq->first ();

        /** We use a shortcut variable in order to avoid some calls to the isDone method for the Sequence iteration.
         * This is important for performance concerns since the 'isDone' kmer iterator relies on it and that we use
         * a generic Iterator<Sequence>; in other words, we have here polymorphism on Sequence iterator and we have
         * to limit such polymorphic calls when the number of calls is huge (overhead due to polymorphism).
         */
        _isDone = _itSeq->isDone();

        /** We check that we have at least one sequence to iterate. */
        if (!_isDone)
        {
            /** We configure the kmer iterator with the current sequence data. */
            _itKmer.setData ((*_itSeq)->getData());

            /** We go to the first kmer. */
            _itKmer.first ();
        }
    }

    /** \copydoc Iterator::next */
    void next()
    {
        /** We look for the next kmer. */
        _itKmer.next ();

        /** We check the case where the kmer iteration is done. */
        if (_itKmer.isDone ())
        {
            /** We have no more kmer for the current sequence, therefore we go for the next sequence. */
            _itSeq->next();

            /** We check whether we have another sequence or not. */
            _isDone = _itSeq->isDone();
            if (!_isDone)
            {

                /** We configure the kmer iterator with the current sequence data. */
                _itKmer.setData ((*_itSeq)->getData());

                /** We go to the first kmer. */
                _itKmer.first ();

                /** We may have to notify potential listeners if we looped enough items. */
                if ((_current & _moduloMask) == 0)  { notifyInc (_current);  _current=0; }

                /** We increase the iterated kmers number. */
                _current++;
            }
        }
    }

    /** \copydoc Iterator::isDone */
    bool isDone()
    {
        /** If we are done, we notify potential listeners. */
        if (_isDone) { notifyFinish(); }

        /** We return the outer loop isDone status. */
        return _isDone;
    }

    /** \copydoc Iterator::item */
    kmer_type& item ()  { return _itKmer.item(); }

private:

    /** Outer loop iterator on Sequence. */
    Iterator<Sequence>* _itSeq;
    void setItSeq (Iterator<Sequence>* itSeq)  { SP_SETATTR(itSeq); }

    /** Inner loop iterator on kmer. */
    Model<kmer_type>::Iterator _itKmer;

    /** Shortcut (for performance). */
    bool _isDone;

    u_int32_t _current;
    u_int32_t _moduloMask;
};

/********************************************************************************/

template <typename Functor> void iter1 (IBank& bank, Model<kmer_type>& model, Functor& hash, IteratorListener* progress=0)
{
    // We need an iterator on the FASTA bank.
    Iterator<Sequence>* itBank = bank.iterator();
    LOCAL (itBank);

    // We declare two kmer iterators for the two banks and a paired one that links them.
    Model<kmer_type>::Iterator itKmer (model);

    // We get some information about the kmers.
    u_int64_t nbKmers       = 0;
    kmer_type checksumKmers = 0;

#if 1
    // We create an iterator over the paired iterator on sequences
    SubjectIterator<Sequence> itSeq (*itBank, 10*1000);

    if (progress)  {  itSeq.addObserver (*progress);  }
#else
    Iterator<Sequence>& itSeq = *itBank;
#endif

    // We get current time stamp
    ITime::Value t0 = System::time().getTimeStamp();

    for (itSeq.first(); !itSeq.isDone(); itSeq.next())
    {
        // We set the data from which we want to extract kmers.
        itKmer.setData (itSeq->getData());

        // We loop the kmers for the two datas.
        for (itKmer.first(); !itKmer.isDone();  itKmer.next())
        {
            nbKmers       += 1;
            checksumKmers = checksumKmers + hash (*itKmer);
        }
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

template <typename Functor> void iter2 (IBank& bank, Model<kmer_type>& model, Functor& hash, IteratorListener* progress=0)
{
    BankKmerIterator itKmerBank (bank, model);

    // We get some information about the kmers.
    u_int64_t nbKmers       = 0;
    kmer_type checksumKmers = 0;

    // We get current time stamp
    ITime::Value t0 = System::time().getTimeStamp();

    if (progress)  {  itKmerBank.addObserver (*progress);  }

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
}

/********************************************************************************/

void BankConvert (IBank& in, IBank& out, IteratorListener* progress=0)
{
#if 1
    // We need an iterator on the FASTA bank.
    Iterator<Sequence>* itBank = in.iterator();
    LOCAL (itBank);

    SubjectIterator<Sequence> itSeq (*itBank, 100*1000);

    if (progress != 0)  {  itSeq.addObserver (*progress);  }
#else

    Bank::Iterator itBank (in);
#endif

    u_int64_t   nbSeq = 0;
    u_int64_t sizeSeq = 0;

    // We get current time stamp
    ITime::Value t0 = System::time().getTimeStamp();

    for (itSeq.first(); !itSeq.isDone(); itSeq.next())
    {
        nbSeq ++;
        sizeSeq += (itSeq)->getDataSize();

        out.insert ( itSeq.item());
    }

    out.flush ();

    // We get current time stamp
    ITime::Value t1 = System::time().getTimeStamp();

    cout << "CONVERSION IN " << (t1-t0) << " msec, "
         << nbSeq << " sequences, " << sizeSeq << " bytes, "
         << "file size " << out.getSize()
         << endl;
}


/********************************************************************************/

int main (int argc, char* argv[])
{
    if (argc < 3)
    {
        cerr << "you must provide at least 2 arguments. Arguments are:" << endl;
        cerr << "   1) kmer size"  << endl;
        cerr << "   2) FASTA bank" << endl;
        return EXIT_FAILURE;
    }

    // We define the max size of a data line in the FASTA output file
    size_t kmerSize = atoi(argv[1]);

    // We get the URI of the FASTA bank
    string filename (argv[2]);

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        // We declare a kmer model with a given span size.
        Model<kmer_type> model (kmerSize);

        // We declare the FASTA bank
        Bank bank (argc-2, argv+2);

        // We declare a binary bank
        BankBinary bankBin (filename + ".bin");

        // We convert the FASTA bank in binary format
        BankConvert (bank, bankBin, new Progress (bank.estimateNbSequences(), "FASTA to binary conversion"));

        //HashNull hash;
        Hash     hash;

        // TEST 1
        iter1 (bankBin, model, hash, new Progress (bank.estimateNbSequences(), "Iterating 1"));

        // TEST 2
        iter2 (bankBin, model, hash, new Progress (bank.estimateNbSequences(), "Iterating 2"));

        // We remove the binary bank
        //System::file().remove (filename + ".bin");
    }

    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }

    return EXIT_SUCCESS;
}
