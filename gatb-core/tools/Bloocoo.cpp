#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>

#include <gatb/kmer/impl/Model.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>
#include <gatb/tools/collections/impl/Bloom.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <iostream>
#include <map>
#include <math.h>

#include <gatb/tools/math/Integer.hpp>

#include <omptl/omptl_numeric>
#include <omptl/omptl_algorithm>

#include <omp.h>

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

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::math;

/********************************************************************************/

struct Hash  {   kmer_type operator () (kmer_type& lkmer)
{
    kmer_type kmer_hash;

    kmer_hash  = lkmer ^ (lkmer >> 14);
    kmer_hash  = (~kmer_hash) + (kmer_hash << 18);
    kmer_hash ^= (kmer_hash >> 31);
    kmer_hash  = kmer_hash * 21;
    kmer_hash ^= (kmer_hash >> 11);
    kmer_hash += (kmer_hash << 6);
    kmer_hash ^= (kmer_hash >> 22);

    return kmer_hash;
}};




/********************************************************************************/

struct Partition
{
    Partition (size_t nbPartitions) : partitions(nbPartitions), _nbPartitions(nbPartitions)
    {
        /** We delete physically the partition files. */
        for (size_t i=0; i<_nbPartitions; i++)  {  System::file().remove (getFilename(i));  }

        /** We create the partition files. */
        for (size_t i=0; i<_nbPartitions; i++)  {  partitions[i] = new BagFile<kmer_type> (getFilename(i));  }
    }

    /** */
    ~Partition ()
    {
        /** We logically delete the partition files. */
        for (size_t i=0; i<_nbPartitions; i++)  {  delete partitions[i];  }
    }

    /** */
    Bag<kmer_type>*& operator[] (size_t idx)  { return partitions[idx]; }

    /** */
    size_t size () const  { return partitions.size(); }

private:

    std::string getFilename (size_t idx)
    {
        char filename[128];  snprintf (filename, sizeof(filename), "foo_partition%d.out", idx);
        return string(filename);
    }

    vector <Bag<kmer_type>*> partitions;
    size_t _nbPartitions;
};

/********************************************************************************/

struct PartitionCache
{
    PartitionCache (Partition& partition, ISynchronizer* synchro) : _ref(partition), cache(partition.size()), _synchro(synchro)
    {
        /** We create the partition files. */
        for (size_t i=0; i<cache.size(); i++)
        {
            cache[i] = new BagCache<kmer_type> (*(partition[i]), 5*1000, synchro);
        }
    }

    PartitionCache (const PartitionCache& p) : _ref(p._ref), cache(p.cache.size()), _synchro(p._synchro)
    {
        /** We create the partition files. */
        for (size_t i=0; i<cache.size(); i++)
        {
            cache[i] = new BagCache<kmer_type> (*(p._ref[i]), 5*1000, _synchro);
        }
    }

    ~PartitionCache ()
    {
        for (size_t i=0; i<cache.size(); i++)
        {
            cache[i]->flush();
            delete cache[i];
        }
    }

    /** */
    Bag<kmer_type>*& operator[] (size_t idx)  { return cache[idx]; }

    /** */
    size_t size () const  { return cache.size(); }

private:

    Partition& _ref;
    vector <Bag<kmer_type>*> cache;
    ISynchronizer* _synchro;
};

/********************************************************************************/

class DSK
{
public:

    /** */
    struct Info
    {
        Info() : nbSeq(0), dataSeq(0), nbKmers(0), checksumKmers(0) {}

        Info (vector<Info>& infos) : nbSeq(0), dataSeq(0), nbKmers(0), checksumKmers(0)
        {
        	for (size_t i=0; i<infos.size(); i++)  {  *this += infos[i];   }
        }

        Info& operator+=  (const Info& other)
        {
            nbSeq         += other.nbSeq;
            dataSeq       += other.dataSeq;
            nbKmers       += other.nbKmers;
            checksumKmers += other.checksumKmers;
            return *this;
        }

        u_int64_t nbSeq;
        u_int64_t dataSeq;
        u_int64_t nbKmers;
        kmer_type checksumKmers;
    };

    /********************************************************************************/
    DSK (const string& filename, size_t kmerSize, ICommandDispatcher& dispatcher, size_t maxMemory)
        : _bankBinary(filename + ".bin"), _filename(filename),
          _model(kmerSize), _kmerSize(kmerSize),
          _dispatcher(dispatcher),
          _estimateSeqNb(0), _estimateSeqTotalSize(0), _estimateSeqMaxSize(0),
          _max_disk_space(0), _max_memory(maxMemory), _volume(0), _nb_passes(0), _nb_partitions(0)
    {}

    /********************************************************************************/
    void configure ()
    {
        // We create a Bank instance.
        Bank bank (_filename);

        // We get some estimations about the bank
        bank.estimate (_estimateSeqNb, _estimateSeqTotalSize, _estimateSeqMaxSize);

        // We may have to build the binary bank if not already existing.
        buildBankBinary (bank);

        // We get the available space (in MBytes) of the current directory.
        u_int64_t available_space = System::file().getAvailableSpace (System::file().getCurrentDirectory()) / 1024;

        u_int64_t bankSize = bank.getSize() / MBYTE;
        u_int64_t kmersNb  = (_estimateSeqTotalSize - _estimateSeqNb * (_kmerSize-1));

        _volume = kmersNb * sizeof(kmer_type) / MBYTE;  // in MBytes

        _max_disk_space = std::min (available_space/2, bankSize);

        if (_max_disk_space == 0)  { _max_disk_space = 10000; }

        _nb_passes = ( _volume / _max_disk_space ) + 1;

        size_t max_open_files = System::file().getMaxFilesNumber() / 2;
        u_int64_t volume_per_pass;

        do  {
            volume_per_pass = _volume / _nb_passes;
            _nb_partitions  = ( volume_per_pass / _max_memory ) + 1;

            if (_nb_partitions >= max_open_files)   { _nb_passes++;  }
            else                                    { break;         }

        } while (1);

//_nb_passes     = 4;
//_nb_partitions = 26;

        cout << "CURRENT DIRECTORY '"    << System::file().getCurrentDirectory() << "'" << endl;
        cout << "AVALAIBLE SPACE     : " << available_space << " MBytes" << endl;
        cout << "BANK SIZE           : " << bankSize << " MBytes" << endl;
        cout << "SEQUENCE NUMBER     : " << _estimateSeqNb << endl;
        cout << "SEQUENCE VOLUME     : " << _estimateSeqTotalSize / MBYTE << " MBytes" << endl;
        cout << "KMER NUMBER         : " << kmersNb << endl;
        cout << "KMERS VOLUME        : " << _volume << " MBytes" << endl;
        cout << "MAX DISK SPACE      : " << _max_disk_space << " MBytes" << endl;
        cout << "MAX MEMORY          : " << _max_memory     << " MBytes" << endl;
        cout << "NB PASSES           : " << _nb_passes << endl;
        cout << "NB PARTITIONS       : " << _nb_partitions << endl;
    }

    /********************************************************************************/
    void execute ()
    {
        /** We configure dsk. */
        configure ();

        /** We delete the solid kmers file. */
        System::file().remove ("solids.bin");

        // We use the provided listener if any
        IteratorListener* progress = new Progress (_estimateSeqNb, "DSK");
        LOCAL (progress);

        // We need an iterator on the FASTA bank.
        Iterator<Sequence>* itBank = _bankBinary.iterator();
        LOCAL (itBank);

        // We create an iterator over the paired iterator on sequences
        SubjectIterator<Sequence> itSeq (*itBank, 5*1000);

        // We add a listener to the sequences iterator.
        itSeq.addObserver (*progress);

        BagFile<kmer_type> solidKmersFile ("solids.bin");
        BagCache<kmer_type> solidKmers (solidKmersFile, 5*1000);

        /** We loop N times the bank. For each pass, we will consider a subset of the whole kmers set of the bank. */
        for (size_t pass=0; pass<_nb_passes; pass++)
        {
            /** We update the iterator listener with some information message. */
            progress->setMessage ("DSK %d/%d", pass+1, _nb_passes);

            /** 1) We fill the partition files. */
            fillPartitions (pass, itSeq);

            /** 2) We fill the kmers solid file from the partition files. */
            fillSolidKmers (solidKmers);
        }

        /** We flush the solid kmers file. */
        solidKmers.flush();
    }

    /********************************************************************************/
    void fillPartitions (size_t pass, Iterator<Sequence>& itSeq)
    {
        /** We create the partition files for the current pass. */
        Partition  partitions (_nb_partitions);

        /** We create a shared synchronizer for the partitions building. */
        ISynchronizer* synchro = System::thread().newSynchronizer();

        /** We launch the iteration of the sequences iterator with the created functors. */
        _dispatcher.iterate (itSeq, FillPartitions (&_model, _nb_passes, pass, partitions, synchro));

        /** We cleanup resources. */
        delete synchro;
    }

    /********************************************************************************/
    void fillSolidKmers (Bag<kmer_type>&  solidKmers)
    {
        TimeInfo ti;

        /** We parse each partition file. */
        for (size_t i=0; i<_nb_partitions; i++)
        {
            char filename[128];  snprintf (filename, sizeof(filename), "foo_partition%d.out", i);

            IteratorFile<kmer_type> it (filename);

            ti.start ("read");

            u_int64_t len = System::file().getSize(filename) / sizeof (kmer_type);
            vector<kmer_type> kmers (len);

            size_t k=0;  for (it.first(); !it.isDone(); it.next())  {  kmers[k++] = *it;  }

            ti.stop ("read");

            ti.start ("sort");

            omptl::sort (kmers.begin (), kmers.end ());

            ti.stop ("sort");

            u_int32_t nks = 3;
            u_int32_t max_couv  = 2147483646;
            u_int32_t abundance = 0;
            kmer_type previous_kmer = kmers.front();

            ti.start ("solid");

            for (vector<kmer_type>::iterator itKmers = kmers.begin(); itKmers != kmers.end(); ++itKmers)
            {
                kmer_type current_kmer = *itKmers;

                if (current_kmer == previous_kmer)  {   abundance++;  }
                else
                {
                    if (abundance >= nks  && abundance <= max_couv)
                    {
                        solidKmers.insert (previous_kmer);
                    }
                    abundance = 1;
                }
                previous_kmer = current_kmer;
            }

            if (abundance >= nks && abundance <= max_couv)
            {
                solidKmers.insert (previous_kmer);
            }
            ti.stop ("solid");

cout << "-> kmersNb: "  << len << "  "
     << "read: "  << ti.getEntryByKey("read")  << "  "
     << "sort: "  << ti.getEntryByKey("sort")  << "  "
     << "solid: " << ti.getEntryByKey("solid") << "  "
     << endl;

        }

    }

    /********************************************************************************/
    class FillPartitions
    {
    public:

        void operator() (Sequence& sequence)
        {
            vector<kmer_type> kmers;

            Hash hash;

            /** We build the kmers from the current sequence. */
            model->build (sequence.getData(), kmers);

            /** We loop over the kmers. */
            for (size_t i=0; i<kmers.size(); i++)
            {
                /** We hash the current kmer. */
                kmer_type h = hash (kmers[i]);

                /** We check whether this kmer has to be processed during the current pass. */
                if ((h % nbPass) != pass)  { continue; }

                kmer_type reduced_kmer = h / nbPass;

                /** We compute in which partition this kmer falls into. */
                size_t p = reduced_kmer % _cache.size();

                /** We write the kmer into the bag. */
                _cache[p]->insert (kmers[i]);
            }
        }

        FillPartitions (KmerModel* model, size_t nbPasses, size_t currentPass, Partition& partition, ISynchronizer* synchro)
            : pass(currentPass), nbPass(nbPasses), _cache(partition,synchro), model(model)
        {
        }

    private:
        size_t pass;
        size_t nbPass;
        PartitionCache _cache;
        KmerModel* model;
    };


private:

    BankBinary _bankBinary;

    string    _filename;
    KmerModel _model;
    size_t    _kmerSize;

    ICommandDispatcher& _dispatcher;

    u_int64_t _estimateSeqNb;
    u_int64_t _estimateSeqTotalSize;
    u_int64_t _estimateSeqMaxSize;

    u_int64_t _max_disk_space;
    u_int32_t _max_memory;
    u_int64_t _volume;
    u_int32_t _nb_passes;
    u_int32_t _nb_partitions;

    /** */
    void buildBankBinary (Bank& bank)
    {
        string filenameBin = _filename + ".bin";

        if (System::file().doesExist(filenameBin) == false)
        {
            // We declare a binary bank
            BankBinary bankBin (filenameBin);

            // We declare some job listener.
            Progress progress (_estimateSeqNb, "FASTA to binary conversion");

            // We convert the FASTA bank in binary format
            IProperties* props = BankHelper::singleton().convert (bank, bankBin, &progress);
            LOCAL (props);
        }
    }
};

/********************************************************************************/
/* Class Bloocoo for read correction : takes as input a bag of solid kmers (dsk result),
 insert that in a bloom filter, and correct read form it*/
/********************************************************************************/

class Bloocoo
{
private:

    size_t              _kmerSize;

    IProperties*        _props;
    ICommandDispatcher* _dispatcher;
    TimeInfo            _timeInfo;
    IteratorListener*   _progress;
    uint64_t            _nb_errors_corrected;
    uint64_t            _seq_num;

public:

    /** */
    Bloocoo (IProperties* props) : _props(props), _progress(0)
    {
        _kmerSize = props->getProperty ("kmer_size")->getInt();

        size_t nbCores = (props->getProperty ("nb_cores") ?  props->getProperty ("nb_cores")->getInt() : 0);
        _dispatcher = new ParallelCommandDispatcher (nbCores);
        _nb_errors_corrected = 0;
        _seq_num = 0;
    }

    /** */
    ~Bloocoo ()
    {
        delete _dispatcher;
        delete _progress;
    }

    /** */
    void execute ()
    {
        /** We delete the extension file. */
        System::file().remove ("extension.bin");

        double lg2 = log(2);
        float NBITS_PER_KMER = log (16*_kmerSize*(lg2*lg2))/(lg2*lg2);

        u_int64_t solidFileSize = (System::file().getSize("solids.bin") / sizeof (kmer_type));

        cout << "Solid kmers  "<< solidFileSize <<  endl;

        u_int64_t estimatedBloomSize = solidFileSize * NBITS_PER_KMER;
        cout << "estimatedBloomSize "<< estimatedBloomSize << " bits  = "<< estimatedBloomSize/1024LL/8LL/1024LL << " MB" <<  endl;

        IteratorFile<kmer_type> itSolid ("solids.bin");
         //  We create a kmer iterator that notifies listeners every 5*1000 kmers 
        SubjectIterator<kmer_type> itKmers (itSolid, 5*1000);

        _progress = new Progress (solidFileSize, "iterate solid kmers");
        itKmers.addObserver (*_progress);

        /** We create a bloom with inserted solid kmers. */
        Bloom<kmer_type>* bloom = createBloom (itKmers);

        
        
        //iterate over initial file
        
        Bank b (_props->getProperty ("input_file")->getString());
        
        // We create a sequence iterator for the bank
        Bank::Iterator itSeq (b);
        
        //  We create some listener to be notified every 1000 iterations and attach it to the iterator.
        SubjectIterator<Sequence> iter (itSeq, 10000);
        Progress fct (b.estimateNbSequences(), "Iterating and correcting sequences");    iter.addObserver (fct);
        
        char fileName[1024];
        sprintf(fileName,"corrected_%s",_props->getProperty ("input_file")->getString());

        Bank outbank (fileName);
        
        //try with a bag cahce to allow parallelization ?
        //BagCache<Sequence> outb_sync (outbank, 1000);

        
        //iterate over sequences and correct them 
        CorrectReads fct2(*bloom, outbank, *this) ; iter.iterate (fct2); //

        //CorrectReads fct2(*bloom) ; iter.iterate (fct2); //

        //_dispatcher->iterate (iter,  CorrectReads (*bloom));
       // _dispatcher->iterate (iter,  CorrectReads (*bloom,outbank));
        //_dispatcher->iterate (iter,  CorrectReads (*bloom,outb_sync));
        

        cout << "NB errors corrected " << _nb_errors_corrected << endl;
        
       // BagFile<kmer_type>  extensionFile ("extension.bin");

        //executeBuildKmerExtension (itKmers, *bloom, extensionFile);

        /** We make sure everything is put into the extension file. */
       // extensionFile.flush();

        /** We clean up resources. */
        delete bloom;
    }

    Properties getProperties ()
    {
        Properties res;

        res.add (0, "bloocoo", "");
        res.add (1, _timeInfo.getProperties("time"));

        return res;
    }

private:

    /********************************************************************************/
    Bloom<kmer_type>* createBloom (Iterator<kmer_type>& itKmers)
    {
        double lg2 = log(2);
        float NBITS_PER_KMER = log (16*_kmerSize*(lg2*lg2))/(lg2*lg2);

        u_int64_t solidFileSize = (System::file().getSize("solids.bin") / sizeof (kmer_type));

        u_int64_t estimatedBloomSize = solidFileSize * NBITS_PER_KMER;
        if (estimatedBloomSize ==0 )estimatedBloomSize = 1000;
        
        /** We instantiate the bloom object. */
        Bloom<kmer_type>* bloom = new BloomSynchronized<kmer_type> (estimatedBloomSize, (int)floorf (0.7*NBITS_PER_KMER));

        _progress->setMessage ("build bloom from solid kmers");

        _timeInfo.start ("BuildKmerBloom");
        _dispatcher->iterate (itKmers,  BuildKmerBloom (*bloom));
        _timeInfo.stop ("BuildKmerBloom");

        return bloom;
    }

    /********************************************************************************/
    void executeBuildKmerExtension (Iterator<kmer_type>& itKmers, Bloom<kmer_type>& bloom, Bag<kmer_type>& extensionBag)
    {
        KmerModel model (_kmerSize);

        ISynchronizer* synchro = System::thread().newSynchronizer();

        _progress->setMessage ("build extension from solid kmers");

        _timeInfo.start ("BuildKmerExtension");
        _dispatcher->iterate (itKmers, BuildKmerExtension (model, bloom, extensionBag, synchro));
        _timeInfo.stop ("BuildKmerExtension");

        delete synchro;
    }

    /********************************************************************************/
    struct BuildKmerBloom
    {
        void operator() (const kmer_type& kmer)  {  _bloom.insert(kmer); }
        BuildKmerBloom (Bloom<kmer_type>& bloom)  : _bloom(bloom) {}
        Bloom<kmer_type>& _bloom;
    };

    
    //functor for read correction this is the code to correct a single read
    /********************************************************************************/

    struct CorrectReads
    {
        void operator() ( Sequence& s) //no const otherwise error with tKmer.setData
        {
            
            // _bloom   is the bloom filter containing the solid kmers
            // _bloom.contains(kmer)  to ask if it contains a kmer
            
            char bin2NT[4] = {'A','C','T','G'}; 
            char bin2NTrev[4] = {'T','G','A','C'};
            char binrev[4]    = {2,3,0,1};
            
            _bloocoo._seq_num++; //counter 1 based, not thread-safe
            //printf("---- new read ----\n");
            
            char * readseq = s.getDataBuffer (); // the nucleotide sequence of the read
            size_t   sizeKmer = _bloocoo._kmerSize;
            
            
            
            KmerModel model (sizeKmer);
            KmerModel::Iterator itKmer (model);
            
            int pass;
            //multiple passes per read
            for (pass=0; pass<2; pass++)
            {
                
                kmer_type graine, graine_revcomp;
                kmer_type current_kmer ;
                kmer_type previous_kmer = 0;
                kmer_type kmer_begin = 0 ;
                kmer_type first_unindexed = 0;
                kmer_type kmer_end = 0;
                bool first_gap = true;
                
                //for each kmer in this Sequence
                
                // sets itKmer to iterate over all the kmers of the read s
                itKmer.setData (s.getData(),KMER_DIRECT);
                
                uint64_t tai_not_indexed =0;
                uint64_t tai_previous_break =0;
                uint64_t tai_indexed = 0;
                int readlen = s.getDataSize();
                bool check = false;
                
                // We iterate the kmers of this sequence
                int ii=0;
                for (itKmer.first(); !itKmer.isDone(); itKmer.next(),ii++)
                {
                    //graine est kmer qui commenece à ii,  dans le sens du read
                    //current_kmer est le min de graine avec son revcomp

                    graine = *itKmer;
                    current_kmer = std::min(  revcomp (graine, sizeKmer),graine );
                    
                    if (_bloom.contains(current_kmer)) //kmer is solid
                    {
                        tai_indexed++;
                        
                        
                        if(tai_indexed==1) //beginning of indexed zone
                        {
                            kmer_end = graine; // kmer_end should be first kmer indexed after a hole
                            
                            if(tai_not_indexed == (sizeKmer))  // this should be an isolated error, middle of the read
                            {
                                
                                //kmer_begin is the last indexed kmer, test its right extension to get the snp
                                char nt;
                                check = false;
                                
                                kmer_type temp_kmer;
                                // kmer_begin.printASCII(sizeKmer);
                                
                                for(nt=0; nt<4; nt++)
                                {
                                    temp_kmer = model.codeSeedRight (kmer_begin, nt, Data::INTEGER);
                                    
                                    //  temp_kmer.printASCII(sizeKmer);
                                    
                                    if(_bloom.contains(temp_kmer)) //kmer is indexed
                                    {
                                        //ref is the last nt of the first non indexed kmer on the reference genome
                                        //contig, pos , ref , snp
                                        //  fprintf(snp_file,"%i  %i  %c %c  \n",j,i-1,bin2NT[first_unindexed & 3],bin2NT[nt]);
                                        check = true;
                                        break;
                                    }
                                }
                                //check other kmers ?  check =false if problem
                                //for higher confidence it would also be possible to test other overlapping kmers
                                
                                if(check)
                                {
                                    readseq[ii-1]=bin2NT[nt]; //correc sequence in ram
                                    // printf("error found pos %i  read %i\n",ii, _bloocoo._seq_num);
                                    _bloocoo._nb_errors_corrected ++; //not thread safe
                                }
                                
                            }
                            else if ((tai_not_indexed < sizeKmer) && first_gap && tai_not_indexed>0)
                                // an error at the beginning of the read : tai_not_indexed is smaller even for a single snp
                                //and correct it from the right
                            {
                                check = false;
                                
                                kmer_type kmer_end_rev, tempkmer;
                                kmer_end_rev =revcomp(kmer_end, sizeKmer);
                                //kmer_end.printASCII(sizeKmer);
                                //kmer_end_rev.printASCII(sizeKmer);
                                
                                int strand =1;
                                int nt2;
                                for(nt2=0; nt2<4; nt2++)
                                {
                                    tempkmer = model.codeSeedRight (kmer_end_rev, nt2, Data::INTEGER); //to go left : go right with reverse of kmer_end
                                    if(_bloom.contains(tempkmer)) //kmer is indexed
                                    {
                                        //tempkmer.printASCII(sizeKmer);
                                        check = true;
                                        
                                        break;
                                    }
                                }
                                
                                
                                //correc error
                                if(check)
                                {
                                    //printf("error found pos %i  read %i\n",tai_not_indexed-1, _bloocoo._seq_num);
                                    //printf("%c \n",readseq[tai_not_indexed-1]);
                                    
                                    readseq[tai_not_indexed-1]=bin2NTrev[nt2]; //correc in ram
                                    _bloocoo._nb_errors_corrected ++; //not thread safe
                                    //printf("%c \n",readseq[tai_not_indexed-1]);
                                    
                                }
                            }
                            
                        }
                        
                        
                        first_gap = false;
                        
                        //   if(tai_indexed==1) tai_previous_break = tai_not_indexed; // begin of indexed zone    ( if tai_not_indexed ?)
                        if (tai_indexed > 1) // do not reset  tai_not_indexed if a single positive (might be a FP)
                            tai_not_indexed = 0; //reset previous hole size
                        if (tai_indexed==1)// begin of indexed zone
                        {
                            //printf("not indexed %lli\n",tai_not_indexed);
                            tai_not_indexed++; //  treat a single solidkmer  as an erroneous kmer
                        }
                    }
                    else //kmer contains an error
                    {
                        
                        
                        
                        if(tai_indexed==1) //begin of not indexed zone,  previous kmer was an isolated positive, probably a FP
                        {
                            
                        }
                        else if(tai_indexed > 1) // begin of not indexed zone
                        {
                            kmer_begin = previous_kmer ; //kmer_begin is the last kmer indexed before the hole
                            first_unindexed = graine;
                            //printf("indexed %lli\n",tai_indexed);
                        }
                        tai_not_indexed ++;
                        tai_indexed =0;
                        
                        
                        if(ii == (readlen-sizeKmer)) //end of the read, we should treat this  gap here
                            //correc snp with trad method
                        {
                            int nt;
                            check = false;
                            
                            kmer_type temp_kmer;
                            
                            for(nt=0; nt<4; nt++)
                            {
                                temp_kmer = model.codeSeedRight (kmer_begin, nt, Data::INTEGER);
                                if(_bloom.contains(temp_kmer)) //kmer is indexed
                                {
                                    check = true;
                                    break;
                                }
                            }
                            //verif other kmers :
                            //todo
                            
                            if(check)
                            {
                                
                                //correc error
                                //printf("%c \n",readseq[readlen-1 - tai_not_indexed +1]);
                                
                                readseq[readlen - tai_not_indexed ]=bin2NT[nt]; //correc sequence in ram
                                //printf("%c \n",readseq[readlen-1 - tai_not_indexed +1]);
                                
                                //printf("error found pos %i  read %i\n",readlen-1 - tai_not_indexed +1, _bloocoo._seq_num);
                                _bloocoo._nb_errors_corrected ++; //not thread safe
                            }
                            
                            
                            
                        }
                        
                    }
                    
                    previous_kmer = graine; // should be the kmer in its original sense
                } // end of kmers iteration over the read
                
                //                if(tai_indexed)
                //                    printf("indexed %lli\n",tai_indexed);
                //                else
                //                    printf("not indexed %lli\n",tai_not_indexed);
                
            }
            
            
            
            _outbank.insert(s); //output corrected sequence
            
            //    printf("%s\n",s.getDataBuffer());
            //    printf("%s\n",readseq);
            
        }
        
        CorrectReads (Bloom<kmer_type>& bloom, Bag<Sequence> & outbank, Bloocoo & bloocoo)  : _bloom(bloom), _outbank(outbank), _bloocoo(bloocoo) {}
      //  CorrectReads (Bloom<kmer_type>& bloom)  : _bloom(bloom)  {}

        Bloom<kmer_type>& _bloom; // the bloom containing the solid kmers
        Bag<Sequence> & _outbank; // the bloom containing the solid kmers
        Bloocoo & _bloocoo;

    };
    /********************************************************************************/
    struct BuildKmerExtension
    {
        void operator() (const kmer_type& kmer)
        {
            _itNeighbors.setSource (kmer);

            for (_itNeighbors.first(); !_itNeighbors.isDone(); _itNeighbors.next())
            {
                if (_bloom.contains (*_itNeighbors))  {  _extendBag.insert (*_itNeighbors);  }
            }
        }

        BuildKmerExtension (KmerModel& model, Bloom<kmer_type>& bloom, Bag<kmer_type>& extendBag, ISynchronizer* synchro)
            : _bloom(bloom), _extendBag(extendBag, 5*1000, synchro), _itNeighbors(model) {}

        ~BuildKmerExtension ()  {  _extendBag.flush();  }

        Bloom<kmer_type>&   _bloom;
        BagCache<kmer_type> _extendBag;
        ModelAbstract<kmer_type>::KmerNeighborIterator _itNeighbors;
    };

};


/********************************************************************************/

int main (int argc, char* argv[])
{
    
    //bloocoo command line
    if (argc < 3)
    {
        cerr << "you must provide at least 3 arguments. Arguments are:" << endl;
        cerr << "   1) kmer size"  << endl;
        cerr << "   2) FASTA  bank" << endl;
        cerr << "   3) Coverage threshold" << endl;        
        cerr << "   4) number of cores" << endl;
        cerr << "   5) memory for kmer counting (in MBytes)" << endl;
        return EXIT_FAILURE;
    }

    // We define the max size of a data line in the FASTA output file
    size_t kmerSize = atoi(argv[1]);

    // We get the URI of the FASTA bank
    string filename (argv[2]);

    
    size_t min_coverage = atoi(argv[3]);

    // We get the number of cores to be used (0 means all cores)
    size_t nbCores = argc >=5 ? atoi(argv[4]) : 0;

    // We get the memory to use (in MBytes)
    size_t maxMemory = argc >=6 ? atoi(argv[5]) : 2000;

    
    ////////
    
    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        TimeInfo ti;

        /** We create a command dispatcher. */
        ICommandDispatcher* dispatcher = 0;
        if (nbCores==1)  { dispatcher =  new SerialCommandDispatcher   ();         }
        else             { dispatcher =  new ParallelCommandDispatcher (nbCores);  }
        LOCAL (dispatcher);

        {
            /** We create an instance of DSK class. */
            DSK dsk (filename, kmerSize, *dispatcher, maxMemory);

            /** We execute dsk. */
            dsk.execute ();
        }
//
//        ///insertion of solid kmers in a bloom
//        u_int64_t solidFileSize = (System::file().getSize("solids.bin") / sizeof (kmer_type));
//
//        IteratorFile<kmer_type> itSolid ("solids.bin");
//        SubjectIterator<kmer_type> itKmers (itSolid, 5*1000);
//        
//        _progress = new Progress (solidFileSize, "iterate solid kmers");
//        itKmers.addObserver (*_progress);
//
//        
//        
//        /** We clean up resources. */
//        delete bloom;
        
        
        IProperties* props = new Properties();
        LOCAL (props);

        props->add (0, "kmer_size", argv[1]);
        props->add (0, "nb_cores",  argc >=4 ? argv[3] : "0");

        props->add (0, "solid_kmers_uri",  "solid.bin");
        props->add (0, "extend_kmers_uri", "extension.bin");
        props->add (0, "input_file", argv[2]);

        
        Bloocoo bloocoo (props);
        bloocoo.execute ();

        /** We dump some execution information. */
        RawDumpPropertiesVisitor visit;
        Properties res = bloocoo.getProperties();
        bloocoo.getProperties().accept (&visit);
    }

    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }

    cout << "finished..." << endl;

    return EXIT_SUCCESS;
}
