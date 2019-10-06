// TODO: this unit test isn't perfect. For k=121, a bug slipped through. one day, remember to code a test that covers the case of counting with large kmers.
//
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

#include <CppunitCommon.hpp>

#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Banks.hpp>
#include <gatb/bank/impl/Bank.hpp>

#include <gatb/kmer/impl/SortingCountAlgorithm.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/BankKmers.hpp>

#include <gatb/tools/misc/api/Macros.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/Histogram.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <gatb/tools/math/LargeInt.hpp>
#include <gatb/tools/math/Integer.hpp>

#include <gatb/tools/storage/impl/Storage.hpp>

#include <boost/variant.hpp>
#include <boost/mpl/for_each.hpp>

using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::storage;
using namespace gatb::core::tools::storage::impl;

using namespace gatb::core::tools::math;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

extern std::string DBPATH (const string& a);

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** Backward compatibility. */
static const int KSIZE_1 = KMER_SPAN(0);
#define KSIZE_32 (KSIZE_LIST == 32)
#if KSIZE_32
#else
static const int KSIZE_2 = KMER_SPAN(1);
static const int KSIZE_3 = KMER_SPAN(2);
#endif
struct Functor_getValue : public boost::static_visitor<Integer>    {
    template<typename T>  Integer operator() (const T& a) const  { return Integer(a.getValue());  }};

/** We define our variant according to the number of defined kmer spans. */
template<typename T>  struct ToKmerVariant  {  typedef typename Kmer<T::value>::Count type;  };
typedef boost::make_variant_over<boost::mpl::transform<IntegerList, ToKmerVariant<boost::mpl::_> >::type >::type KmerVariant;

/** \brief Test class for genomic databases management
 */
class TestDSK : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestDSK);

        CPPUNIT_TEST_GATB (DSK_check1);
        CPPUNIT_TEST_GATB (DSK_check2);
        CPPUNIT_TEST_GATB (DSK_check3);
        CPPUNIT_TEST_GATB (DSK_perBank1);
        CPPUNIT_TEST_GATB (DSK_perBank2);
        CPPUNIT_TEST_GATB (DSK_perBankKmer);
        CPPUNIT_TEST_GATB (DSK_multibank);
		 

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    static const CountNumber NKS_MAX = (CountNumber)(1<<30);  // can't use std::numeric_limits not available on old gcc

    /********************************************************************************/
    void setUp    () {}
    void tearDown () {}


        // SMALL VALUE NEEDED because continuous integration servers are not very powerful...
        static const u_int64_t MAX_MEMORY = 1000;
        
       
    /********************************************************************************/
    void DSK_check1_aux (const char* sequences[], size_t nbSequences, size_t kmerSize, size_t nks, size_t checkNbSolids)
    {
        /** We configure parameters for a SortingCountAlgorithm object. */
        IProperties* params = SortingCountAlgorithm<>::getDefaultProperties();  LOCAL (params);
        params->setInt (STR_KMER_SIZE,          kmerSize);
        params->setInt (STR_MAX_MEMORY,         MAX_MEMORY);
        params->setInt (STR_KMER_ABUNDANCE_MIN, nks);
        params->setStr (STR_URI_OUTPUT,         "foo");

        /** We create a DSK instance. */
        SortingCountAlgorithm<> dsk (new BankStrings (sequences, nbSequences), params);

        /** We launch DSK. */
        dsk.execute();

        if ((int) checkNbSolids != dsk.getInfo()->getInt("kmers_nb_solid"))
        {
            std::cout << "problem with sequences " << sequences << " kmersize " << kmerSize << " nks " << nks << " expected " << checkNbSolids << " solids, had " << dsk.getInfo()->getInt("kmers_nb_solid") <<std::endl;

        }
        CPPUNIT_ASSERT ((int) checkNbSolids == dsk.getInfo()->getInt("kmers_nb_solid"));
    }

    /********************************************************************************/
    void DSK_check1 ()
    {
        const char* s1 = "GATCCTCCCCAGGCCCCTACACCCAAT" ;

        const char* seqs1[] =  { s1 };
        DSK_check1_aux (seqs1, ARRAY_SIZE(seqs1), 27, 1, 1);
        DSK_check1_aux (seqs1, ARRAY_SIZE(seqs1), 26, 1, 2);
        DSK_check1_aux (seqs1, ARRAY_SIZE(seqs1), 27, 2, 0);
        DSK_check1_aux (seqs1, ARRAY_SIZE(seqs1), 26, 2, 0);

        const char* seqs2[] =  { s1, s1 };
        DSK_check1_aux (seqs2, ARRAY_SIZE(seqs2), 27, 1, 1);
        DSK_check1_aux (seqs2, ARRAY_SIZE(seqs2), 26, 1, 2);
        DSK_check1_aux (seqs2, ARRAY_SIZE(seqs2), 27, 2, 1);
        DSK_check1_aux (seqs2, ARRAY_SIZE(seqs2), 26, 2, 2);
        DSK_check1_aux (seqs2, ARRAY_SIZE(seqs2), 27, 3, 0);
        DSK_check1_aux (seqs2, ARRAY_SIZE(seqs2), 26, 3, 0);

        const char* seqs3[] =  { s1, s1, s1};
        DSK_check1_aux (seqs3, ARRAY_SIZE(seqs3), 27, 1, 1);
        DSK_check1_aux (seqs3, ARRAY_SIZE(seqs3), 26, 1, 2);
        DSK_check1_aux (seqs3, ARRAY_SIZE(seqs3), 27, 2, 1);
        DSK_check1_aux (seqs3, ARRAY_SIZE(seqs3), 26, 2, 2);
        DSK_check1_aux (seqs3, ARRAY_SIZE(seqs3), 27, 3, 1);
        DSK_check1_aux (seqs3, ARRAY_SIZE(seqs3), 26, 3, 2);
        DSK_check1_aux (seqs3, ARRAY_SIZE(seqs3), 27, 4, 0);
        DSK_check1_aux (seqs3, ARRAY_SIZE(seqs3), 26, 4, 0);

        const char* seqs4[] = {

            "CGCTACAGCAGCTAGTTCATCATTGTTTATCAATGATAAAATATAATAAGCTAAAAGGAAACTATAAATA"
            "ACCATGTATAATTATAAGTAGGTACCTATTTTTTTATTTTAAACTGAAATTCAATATTATATAGGCAAAG"
            "ACTTAGATGTAAGATTTCGAAGACTTGGATGTAAACAACAAATAAGATAATAACCATAAAAATAGAAATG"
            "AACGATATTAAAATTAAAAAATACGAAAAAACTAACACGTATTGTGTCCAATAAATTCGATTTGATAATT"
            "AGGTAACAATTTAACGTTAAAACCTATTCTTTTATTATCCGAAAATCCGTCGTGGAATTTGTATTAGCTT"
            "TTTTTCTACATTACCCGTTTGCGAGACAGGTGGGGTCAGACGTAGACGTAGTCTCTGGAGTCAAGACGAA"
            "ATTTTACATTTCACAATTTCCTATAGGCCGAGCAAAATTTATTAAGAACCCACAGGCATCATTACGTTTT"
            "CTTGCACAGAAGACTTCACGCTGAAGTCATTGGGCTATATTTCAACGAGACGTCTGTTGGTTTATAAAGG"
            "GCTATATTTATACAAGAATAGGAGTATGGCAGTATGCTAGGCTGGTATGTAGTACGTATACCTCCTAAGC"
            "CGAAAGGCAGTAAGTGACGATGTAATAGTTTTGAGGAAAATTACTTTTTCTGAATAATATTTTTATTTTT"
            "GTTTGCATTTTGTTAAAATTATTTACTAAATTAATGATTCTCATATGTTTTTTCATAGATTTGATGAACT"
            "ACTGTACCATCTGATTAGCGCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCA"
            "CACAACATACGAGCCGGAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTACTCACATAATTG"
            "CGTGCGCTCACTGCCCGCTTTCCAGTCGGAAACCTGTCGTGCCAGCTGCATTATTGATCCGCCACGCCCC"
            "GGGGAAAGT",

            "GTCTTCATTCAGCTGTTCTCATGATAACTAGTAATTCCTTGCTAACAATTTTTACTGAGTAGCAACCAAT"
            "TAATGTTGCCAGAATTTCATAATTGAATTTGAATTTTTTATTTTTTCCTTGATTATGCTTCAAACTCTAT"
            "GTAGTTATTTAGAGTCAATAATATTAAAGCAATCTTAATATTAACTCATTTATTTCTGATTGGCCATATT"
            "TATTTAATTCTCAACAATAATAATGATAAGTATAATAATATATTTAACTTAATAACATTTTAATCATTTT"
            "ATTTTTGTTTGTTGTGATTTTTGGACGTTGTGGTAAATAAGAAGTTTTAAGCTTATATTAATATGTTTTA"
            "CTTTTTATTTCTTAATACGAATTTAATTACCTACCCATTATATTAAGTATATGTTTTGGAATTCTTTCTG"
            "TAAAAATGTGTTTTAAATATTTTACACTTAATTATGTAGGTACCTATACATTTTTAAACTTATCGTATAA"
            "TTCTTTTAATGGTTAAATCATACAAATTAATGTGTAGAGAATAGTTTTTATAAGACTCGTTGTCAATACG"
            "TACGCATAATATAAAAAAACTGACATGTTTTAGTAAGTCGTTTTGATGCATAATAGGATTTTTACCTTTT"
            "AAAGTCTCAAGTTTTCATACAGTGGTACCTCTATATAGAACACGTTAGGCTTTACGGGGTCATTATTTCT"
            "GTTCCGATATTTTTAATGGCATAAAACTATAAACAATAACCGGTATGTATAAATGGTAC",

            "ACCTGAAGCAGTAGTTCATCATATTGCGACTGCAGAATCGATGATAAAGTGGCTTTTAGATCTAAAAGCC"
            "AATACAAAACTGAAGGAATTTGATTTGATGGATTTTAATTTTGAAAATGGATTATGATTGTCGATTGATT"
            "AACAAGTTTACTAGGTTTGAATAGAGGTGATTCTTAATATTTCAAATATTTGAATGTCATGATGAATATT"
            "ATAATTTATAATTAAAAAATATCATATTTTATTCATGGATATCAAAGCTGAAAAAATAGATATTCAAAAT"
            "CGCCTTTATAATAACCTATCATAAACTAATTAATCAATTAAATTCAGTTTTAAAAATTTAAATCCGACAA"
            "ATAAAATTCCTTCAGCTCTGTCTGGGATTTTGGTCGAAAAATTTTAAATCGAAAAAAGTTTATCTTATTC"
            "ATAATATCATTGCCAATGATATTAAAATTAATTAACAACGAATACAAATAACGTCCGACCTGTATATTGC"
            "GGGCCAACTGTTTTTATAGGAAATGTTGACCGAAAACTATTACAGATTAGATGTGTGTGTGTTTACCCTG"
            "TACAAAAATACAAGTACTATTACAACACATCATTATGTTAAATTGCCTCTATATTAATTTCTTTAAAACA"
            "CGACCAACTGCACATTAAAGTAGTTTATTTAGTACTACAGTAGATTAAATTCATTTTTGACGAAAAATTG"
            "CATTTGAAAATGGCCATTGTGTGTATAAATATTGTATACTAATATAACTCTAAATAAAGGTTTCCAGTAC"
            "CAAAGAACCAAATTTTTAATTACAACCAAAATAACTAAATCGTATTCTTTGTTAAATAGTTAAGTTTTTC"
            "GCCGATTGCTGTGCTTGACAGTCTCCTCAATTCAGAATTTCATGTAAAATAAAAATAGCGTACATATAAT"
            "GGATTGCTGTGGCATTTGGTTTGATTAATCCCAAATATTGATTCCAAATATCTATTAGCCTATTGTACCC"
            "CGGAGTACCG"
        } ;

        /** The following solid kmers numbers are computed with the original minia. */
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 9, 1, 2540);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 9, 2, 151);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 9, 3, 18);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 9, 4, 3);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 9, 5, 2);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 9, 6, 0);

        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 11, 1, 2667);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 11, 2, 41);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 11, 3, 0);

        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 13, 1, 2690);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 13, 2, 12);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 13, 3, 0);

        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 15, 1, 2691);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 15, 2, 5);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 15, 3, 0);
    }

    /********************************************************************************/
    template<size_t span>
    void DSK_check2_aux ()
    {
        /** Shortcut. */
        typedef typename Kmer<span>::Type  Type;
        typedef typename Kmer<span>::Count Count;

        size_t kmerSize = 31;
        size_t nks      = 1;

        const char* s1 = "GATCGATTCTTAGCACGTCCCCCCCTACACCCAAT" ;

        /** We configure parameters for a SortingCountAlgorithm object. */
        IProperties* params = SortingCountAlgorithm<>::getDefaultProperties();
        params->setInt (STR_KMER_SIZE,          kmerSize);
        params->setInt (STR_MAX_MEMORY,         MAX_MEMORY);
        params->setInt (STR_KMER_ABUNDANCE_MIN, nks);
        params->setStr (STR_URI_OUTPUT,         "foo");

        /** We create a DSK instance. */
        SortingCountAlgorithm<span> sortingCount (new BankStrings (s1, 0), params);

        /** We launch DSK. */
        sortingCount.execute();

        /** We iterate the solid kmers. */
        Iterator<Count>* iter = sortingCount.getSolidCounts()->iterator();
        LOCAL (iter);

        /** The following values have been computed with the original DSK.
         *  Note: we use a set and check that the solid iterator items are in this trustful set.
         *  -> we have to do this because we are not sure about the order of the iterated items.
         */
        set<Type> okValues;
        Type v1; v1.setVal( 0x1CA68D1E55561150 );
        okValues.insert (v1);
        Type v2; v2.setVal( 0x09CA68D1E5556115);
        okValues.insert (v2);
        Type v3; v3.setVal( 0x2729A34795558454);
        okValues.insert (v3);
        Type v4; v4.setVal( 0x32729A3479555845); 
        okValues.insert (v4);
        Type v5; v5.setVal( 0x0AFEE3FFF1ED8309);
        okValues.insert (v5);

        set<Type> checkValues;
        Type checksum;
        checksum.setVal(0);

        size_t idx=0;
        for (iter->first(); !iter->isDone(); iter->next(), idx++)
        {
            typename set<Type>::iterator lookup = okValues.find (iter->item().value);
            CPPUNIT_ASSERT (lookup != okValues.end());

            checkValues.insert (iter->item().value);

            checksum += iter->item().value;
        }

        CPPUNIT_ASSERT (checksum == 0x8b0c176c3b43d207);
        CPPUNIT_ASSERT (checkValues.size() == okValues.size());

        /** TO BE IMPROVED => dynamic_cast now fails because of the new storage kind of solid kmers. */
#if 0
        /** We check the result through the variant type. */
        IteratorVariant <IteratorFile, Kmer<KSIZE_1>::Count, Kmer<KSIZE_2>::Count, Kmer<KSIZE_3>::Count, Kmer<KSIZE_4>::Count > itVariant;

        IteratorFile<Count>* solid = dynamic_cast<IteratorFile<Count>*> (iter);
        CPPUNIT_ASSERT (solid != 0);

        /** We set the variant with the current T type. */
        itVariant = *solid;

        Integer checksumGeneric (Type(0));

        for (itVariant.first(); !itVariant.isDone(); itVariant.next())
        {
            boost::variant<Kmer<KSIZE_1>::Count, Kmer<KSIZE_2>::Count, Kmer<KSIZE_3>::Count, Kmer<KSIZE_4>::Count> current = itVariant.item();

            Integer val = boost::apply_visitor (Functor_getValue(),  current);

            checksumGeneric += val;
        }
        CPPUNIT_ASSERT (checksumGeneric == Integer(Type(0x8b0c176c3b43d207)));
#endif
    }

    /********************************************************************************/
    void DSK_check2 ()
    {
        DSK_check2_aux<KSIZE_1> ();
#if KSIZE_32
#else
        DSK_check2_aux<KSIZE_2> ();
        DSK_check2_aux<KSIZE_3> ();
#endif
    }

    /********************************************************************************/

    template<size_t span>
    void DSK_check3_aux (IBank* bank, size_t kmerSize, size_t nks)
    {
        /** Shortcut. */
        typedef typename Kmer<span>::Count Count;

        LOCAL (bank);

        TimeInfo ti;

        /** We configure parameters for a SortingCountAlgorithm object. */
        IProperties* params = SortingCountAlgorithm<>::getDefaultProperties();  LOCAL (params);
        params->setInt (STR_KMER_SIZE,          kmerSize);
        params->setInt (STR_MAX_MEMORY,         MAX_MEMORY);
        params->setInt (STR_KMER_ABUNDANCE_MIN, nks);
        params->setStr (STR_URI_OUTPUT,         "foo");

        /** We create a DSK instance. */
        SortingCountAlgorithm<span> sortingCount (bank, params);

        /** We launch DSK. */
        sortingCount.execute();

        /** We iterate the solid kmers. */
        Iterator<Count>* iter = sortingCount.getSolidCounts()->iterator();
        LOCAL (iter);

        // cout << "----------------------------------------------------------" << endl;

        /** TO BE IMPROVED => dynamic_cast now fails because of the new storage kind of solid kmers. */
#if 0
        /** We check the result through the variant type. */
        IteratorVariant <IteratorFile, Kmer<KSIZE_1>::Count, Kmer<KSIZE_2>::Count, Kmer<KSIZE_3>::Count, Kmer<KSIZE_4>::Count > itVar;

        IteratorFile<Count>* solid = dynamic_cast<IteratorFile<Count>*> (iter);
        CPPUNIT_ASSERT (solid != 0);

        /** We set the variant with the current T type. */
        itVar = *solid;

        Type    checksum1 = 0;
        Integer checksum2 = Integer(Type(0));

        PairedIterator <IteratorFile, Count, KmerVariant> itBoth (*solid, itVar);
        for (itBoth.first(); !itBoth.isDone(); itBoth.next())
        {
            Type    v1 = itBoth.item().first.getValue();
            Integer v2 = boost::apply_visitor (Functor_getValue(),  itBoth.item().second);

            checksum1 += (hash1 (v1,0) & 0x11111);
            checksum2 += Integer(Type((hash1 (v2,0) & 0x11111)));

             //cout << "["  << (Integer(v1) == v2 ? 'X' : ' ')  << "] "  << v1 << "  " <<  v2 << endl;
        }

        CPPUNIT_ASSERT (checksum2 == Integer(checksum1));
#endif

        /** Some performance tests. */

        size_t idx2 = 0;
        size_t idx1 = 0;

        {   TIME_INFO(ti,"1");  for (iter->first(); !iter->isDone(); iter->next()) { idx1++; }  }

        {   TIME_INFO(ti,"2");  for (iter->first(); !iter->isDone(); iter->next()) { idx2++; }  }
    }

    /** */
    void DSK_check3 ()
    {
        size_t kmerSize = 15;
        size_t nks      = 1;

        IBank* bank = new BankStrings (
            "CGCTACAGCAGCTAGTTCATCATTGTTTATCAATGATAAAATATAATAAGCTAAAAGGAAACTATAAATA"
            "ACCATGTATAATTATAAGTAGGTACCTATTTTTTTATTTTAAACTGAAATTCAATATTATATAGGCAAAG"
            "ACTTAGATGTAAGATTTCGAAGACTTGGATGTAAACAACAAATAAGATAATAACCATAAAAATAGAAATG"
            "AACGATATTAAAATTAAAAAATACGAAAAAACTAACACGTATTGTGTCCAATAAATTCGATTTGATAATT"
            "AGGTAACAATTTAACGTTAAAACCTATTCTTTTATTATCCGAAAATCCGTCGTGGAATTTGTATTAGCTT"
            "TTTTTCTACATTACCCGTTTGCGAGACAGGTGGGGTCAGACGTAGACGTAGTCTCTGGAGTCAAGACGAA"
            "ATTTTACATTTCACAATTTCCTATAGGCCGAGCAAAATTTATTAAGAACCCACAGGCATCATTACGTTTT"
            "CTTGCACAGAAGACTTCACGCTGAAGTCATTGGGCTATATTTCAACGAGACGTCTGTTGGTTTATAAAGG"
            "GCTATATTTATACAAGAATAGGAGTATGGCAGTATGCTAGGCTGGTATGTAGTACGTATACCTCCTAAGC",
            0
        );
        LOCAL (bank);

        DSK_check3_aux<KSIZE_1> (bank, kmerSize, nks);
#if KSIZE_32
#else
        DSK_check3_aux<KSIZE_2> (bank, kmerSize, nks);
        DSK_check3_aux<KSIZE_3> (bank, kmerSize, nks);
#endif
    }

    /********************************************************************************/
    template<size_t span>
    void DSK_perBank_aux (IBank* bank, size_t kmerSize, size_t nksMin, size_t nksMax, KmerSolidityKind solidityKind, size_t checkNb)
    {
        size_t maxDiskSpace = 0;
        size_t nbCores      = 1;

        /** We configure parameters for a SortingCountAlgorithm object. */
        IProperties* params = SortingCountAlgorithm<>::getDefaultProperties();
        params->setInt (STR_KMER_SIZE,          kmerSize);
        params->setInt (STR_KMER_ABUNDANCE_MIN, nksMin);
        params->setInt (STR_KMER_ABUNDANCE_MAX, nksMax);
        params->setInt (STR_MAX_MEMORY,         MAX_MEMORY);
        params->setInt (STR_MAX_DISK,           maxDiskSpace);
        params->setStr (STR_SOLIDITY_KIND,      toString(solidityKind));
        params->setStr (STR_URI_OUTPUT,         "output");

        params->add (0, STR_NB_CORES,  "%d", nbCores);

        /** We create a DSK instance. */
        SortingCountAlgorithm<span> sortingCount (bank, params);

        /** We launch DSK. */
        sortingCount.execute();

        // cout << *sortingCount.getInfo() << endl;
        // Group& sortingGroup = storage->getGroup("dsk");
        // Collection<IHistogram::Entry>& histo = sortingGroup.getCollection<IHistogram::Entry>("histogram");
        // Iterator<IHistogram::Entry>* itHisto = new TruncateIterator<IHistogram::Entry> (*histo.iterator(), 12);
        // for (itHisto->first(); !itHisto->isDone(); itHisto->next())  { cout << itHisto->item().index << " " << itHisto->item().abundance << endl; }
        // printf ("min=%ld  max=%ld  nb=%ld  check=%ld \n",
        //    nksMin, nksMax, sortingCount.getSolidCounts()->getNbItems(),checkNb
        // );
        
        if (sortingCount.getSolidCounts()->getNbItems() != (int)checkNb)
            std::cout << "counted " <<sortingCount.getSolidCounts()->getNbItems()<< " kmers, expected " << (int)checkNb << std::endl;

        CPPUNIT_ASSERT (sortingCount.getSolidCounts()->getNbItems() == (int)checkNb);
    }

    /********************************************************************************/
    void DSK_perBank1 ()
    {
        CountNumber nksMax = NKS_MAX;  // a big value here

        const char* seqs[] = {      //  KMERS ARE...
            "CGCTACAGCAGCTAGTT",    // CGCTACAGCAGCTAG  GCTACAGCAGCTAGT  CTACAGCAGCTAGTT
            "GCTACAGCAGCTAGTTA",    //                  GCTACAGCAGCTAGT  CTACAGCAGCTAGTT  TACAGCAGCTAGTTA
            "CTACAGCAGCTAGTTAC"     //                                   CTACAGCAGCTAGTT  TACAGCAGCTAGTTA  ACAGCAGCTAGTTAC
        };

        // CTACAGCAGCTAGTT is the only kmer to be present at least in each sequence
        // We have 5 different kmers

        BankComposite* album = new BankAlbum ("foo", true);  LOCAL (album);

        for (size_t i=0; i<ARRAY_SIZE(seqs); i++) {   album->addBank (new BankStrings(seqs[i],NULL)); }

        // Checks for abundance==1
        DSK_perBank_aux<KSIZE_1> (album, 15, 1, nksMax, KMER_SOLIDITY_MIN, 1);
        DSK_perBank_aux<KSIZE_1> (album, 15, 1, nksMax, KMER_SOLIDITY_MAX, 5);
        DSK_perBank_aux<KSIZE_1> (album, 15, 1, nksMax, KMER_SOLIDITY_SUM, 5);

        // Checks for abundance==2
        DSK_perBank_aux<KSIZE_1> (album, 15, 2, nksMax, KMER_SOLIDITY_MIN, 0);
        DSK_perBank_aux<KSIZE_1> (album, 15, 2, nksMax, KMER_SOLIDITY_MAX, 0);
        DSK_perBank_aux<KSIZE_1> (album, 15, 2, nksMax, KMER_SOLIDITY_SUM, 3);

        // Checks for abundance==3
        DSK_perBank_aux<KSIZE_1> (album, 15, 3, nksMax, KMER_SOLIDITY_MIN, 0);
        DSK_perBank_aux<KSIZE_1> (album, 15, 3, nksMax, KMER_SOLIDITY_MAX, 0);
        DSK_perBank_aux<KSIZE_1> (album, 15, 3, nksMax, KMER_SOLIDITY_SUM, 1);
    }

    /********************************************************************************/
    void DSK_perBank2 ()
    {
        CountNumber nksMax = NKS_MAX;

        const char* seqs[] = {
            "CGCTATCGCTA",    // CGCTA  GCTAT  CTATC  TATCG  ATCGC  TCGCT  CGCTA
            "CGCTATAGTTA",    // CGCTA  GCTAT  CTATA  TATAG  ATAGT  TAGTT  AGTTA
            "CGCTAACGCTA"     // CGCTA  GCTAA  CTAAC  TAACG  AACGC  ACGCT  CGCTA
        };

        //  KMERS ARE... (with revcomp)
        // CGCTA-TAGCG  GCTAT-ATAGC  CTATC-GATAG  TATCG-CGATA  ATCGC-GCGAT  TCGCT-AGCGA  CGCTA-TAGCG
        // CGCTA-TAGCG  GCTAT-ATAGC  CTATA-TATAG  TATAG-CTATA  ATAGT-ACTAT  TAGTT-AACTA  AGTTA-TAACT
        // CGCTA-TAGCG  GCTAA-TTAGC  CTAAC-GTTAG  TAACG-CGTTA  AACGC-GCGTT  ACGCT-AGCGT  CGCTA-TAGCG

        // CANONICAL KMERS ARE
        // CGCTA  ATAGC  CTATC  CGATA  ATCGC  AGCGA  CGCTA
        // CGCTA  ATAGC  CTATA  CTATA  ACTAT  AACTA  AGTTA
        // CGCTA  TTAGC  CTAAC  CGTTA  AACGC  ACGCT  CGCTA

        // UNIQUE KMERS ARE
        // (CGCTA,5)  (ATAGC,2)  (CTATC,1)  (CGATA,1)  (ATCGC,1)  (AGCGA,1)
        // (CTATA,2)  (ACTAT,1)  (AACTA,1)  (AGTTA,1)
        // (TTAGC,1)  (CTAAC,1)  (CGTTA,1)  (AACGC,1)  (ACGCT,1)

        // OCCURRENCES in each bank
        // CGCTA(2) ATAGC(1) CTATC(1) CGATA(1) ATCGC(1) AGCGA(1)
        // CGCTA(1) ATAGC(1)                                     CTATA(2) ACTAT(1) AACTA(1) AGTTA(1)
        // CGCTA(2)                                                                                  TTAGC(1) CTAAC(1) CGTTA(1) AACGC(1) ACGCT(1)

        BankComposite* album = new BankAlbum ("foo", true);  LOCAL (album);

        for (size_t i=0; i<ARRAY_SIZE(seqs); i++) {   album->addBank (new BankStrings(seqs[i],NULL)); }

        // Checks for abundance==1
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, nksMax, KMER_SOLIDITY_MIN, 1);
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, nksMax, KMER_SOLIDITY_MAX, 15);
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, nksMax, KMER_SOLIDITY_SUM, 15);
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, nksMax, KMER_SOLIDITY_ALL, 1);
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, nksMax, KMER_SOLIDITY_ONE, 15);

        // Checks for abundance==2
        DSK_perBank_aux<KSIZE_1> (album, 5, 2, nksMax, KMER_SOLIDITY_MIN, 0);
        DSK_perBank_aux<KSIZE_1> (album, 5, 2, nksMax, KMER_SOLIDITY_MAX, 2);
        DSK_perBank_aux<KSIZE_1> (album, 5, 2, nksMax, KMER_SOLIDITY_SUM, 3);
        DSK_perBank_aux<KSIZE_1> (album, 5, 2, nksMax, KMER_SOLIDITY_ALL, 0);
        DSK_perBank_aux<KSIZE_1> (album, 5, 2, nksMax, KMER_SOLIDITY_ONE, 2);

        // Checks for abundance==3
        DSK_perBank_aux<KSIZE_1> (album, 5, 3, nksMax, KMER_SOLIDITY_MIN, 0);
        DSK_perBank_aux<KSIZE_1> (album, 5, 3, nksMax, KMER_SOLIDITY_MAX, 0);
        DSK_perBank_aux<KSIZE_1> (album, 5, 3, nksMax, KMER_SOLIDITY_SUM, 1);
        DSK_perBank_aux<KSIZE_1> (album, 5, 3, nksMax, KMER_SOLIDITY_ALL, 0);
        DSK_perBank_aux<KSIZE_1> (album, 5, 3, nksMax, KMER_SOLIDITY_ONE, 0);

        // Checks for abundance=[1,1]  => NOTE THAT EACH SOLIDITY KIND LEADS TO A DIFFERENT VALUE...
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, 1, KMER_SOLIDITY_MIN, 1);
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, 1, KMER_SOLIDITY_MAX, 13);
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, 1, KMER_SOLIDITY_SUM, 12);
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, 1, KMER_SOLIDITY_ALL, 0);
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, 1, KMER_SOLIDITY_ONE, 14);

        // Checks for abundance=[1,2]
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, 2, KMER_SOLIDITY_MIN, 1);
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, 2, KMER_SOLIDITY_MAX, 15);
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, 2, KMER_SOLIDITY_SUM, 14);
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, 2, KMER_SOLIDITY_ALL, 1);
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, 2, KMER_SOLIDITY_ONE, 15);

        // Checks for abundance=[1,3]
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, 3, KMER_SOLIDITY_MIN, 1);
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, 3, KMER_SOLIDITY_MAX, 15);
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, 3, KMER_SOLIDITY_SUM, 14);
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, 3, KMER_SOLIDITY_ALL, 1);
        DSK_perBank_aux<KSIZE_1> (album, 5, 1, 3, KMER_SOLIDITY_ONE, 15);

        // Checks for abundance=[2,2]
        DSK_perBank_aux<KSIZE_1> (album, 5, 2, 2, KMER_SOLIDITY_MIN, 0);
        DSK_perBank_aux<KSIZE_1> (album, 5, 2, 2, KMER_SOLIDITY_MAX, 2);
        DSK_perBank_aux<KSIZE_1> (album, 5, 2, 2, KMER_SOLIDITY_SUM, 2);
        DSK_perBank_aux<KSIZE_1> (album, 5, 2, 2, KMER_SOLIDITY_ALL, 0);
        DSK_perBank_aux<KSIZE_1> (album, 5, 2, 2, KMER_SOLIDITY_ONE, 2);

        // Checks for abundance=[3,3]
        DSK_perBank_aux<KSIZE_1> (album, 5, 3, 3, KMER_SOLIDITY_MIN, 0);
        DSK_perBank_aux<KSIZE_1> (album, 5, 3, 3, KMER_SOLIDITY_MAX, 0);
        DSK_perBank_aux<KSIZE_1> (album, 5, 3, 3, KMER_SOLIDITY_SUM, 0);
        DSK_perBank_aux<KSIZE_1> (album, 5, 3, 3, KMER_SOLIDITY_ALL, 0);
        DSK_perBank_aux<KSIZE_1> (album, 5, 3, 3, KMER_SOLIDITY_ONE, 0);

        // Checks for abundance=[3,5]
        DSK_perBank_aux<KSIZE_1> (album, 5, 3, 5, KMER_SOLIDITY_MIN, 0);
        DSK_perBank_aux<KSIZE_1> (album, 5, 3, 5, KMER_SOLIDITY_MAX, 0);
        DSK_perBank_aux<KSIZE_1> (album, 5, 3, 5, KMER_SOLIDITY_SUM, 1);
        DSK_perBank_aux<KSIZE_1> (album, 5, 3, 5, KMER_SOLIDITY_ALL, 0);
        DSK_perBank_aux<KSIZE_1> (album, 5, 3, 5, KMER_SOLIDITY_ONE, 0);
    }

    /********************************************************************************/
    void DSK_perBankKmer_aux (size_t kmerSize, size_t nbBanksMax)
    {
        /** We create a bank holding all 4^k kmers =>  we will have 4^k/2 canonical kmers with abundance==2  */
        IBank* kmersBank = new BankKmers(kmerSize);  LOCAL (kmersBank);

        u_int64_t nbKmers          = (1<<(2*kmerSize));
        u_int64_t nbKmersCanonical = nbKmers / 2;

        vector<IBank*> banks;
        vector<IBank*> banksComposite;

        for (size_t i=0; i<nbBanksMax; i++)
        {
            /** We add the bank into the vector of banks. */
            banks.push_back (kmersBank);

            /** We create a new album. */
            stringstream ss;  ss << "foo" << i;

            BankComposite* compo = new BankAlbum (ss.str(), true);
            for (size_t j=0; j<banks.size(); j++)  { compo->addBank (kmersBank); }
            banksComposite.push_back (compo);
        }

        for (size_t i=0; i<banksComposite.size(); i++)  { banksComposite[i]->use();  }

        size_t abundMax = 100000;

        /** Now, we should have N albums Ai (i in [1..N]), with Ai holding i times the kmers bank. */
        for (size_t i=0; i<banksComposite.size(); i++)
        {
            /** Shortcut. */
            IBank* current = banksComposite[i];

            /** We check the number of items. */
            CPPUNIT_ASSERT ((u_int64_t)current->getNbItems() == (i+1)*nbKmers);

            for (size_t abundMin=1; abundMin<=nbBanksMax; abundMin++)
            {
                /** For mode "min" : each bank can have at most coverage==2 for each kmer. */
                size_t checkMin = abundMin > 2 ? 0 : nbKmersCanonical;

                /** For mode "max" : all bank have coverage==2 for each kmer. */
                size_t checkMax = abundMin > 2 ? 0 : nbKmersCanonical;

                /** For mode "sum" : N times composite banks have 2N coverage for each kmer. */
                size_t checkSum = abundMin > 2*(i+1) ? 0 : nbKmersCanonical;

                DSK_perBank_aux<KSIZE_1> (current, kmerSize, abundMin, abundMax, KMER_SOLIDITY_MIN, checkMin);
                DSK_perBank_aux<KSIZE_1> (current, kmerSize, abundMin, abundMax, KMER_SOLIDITY_MAX, checkMax);
                DSK_perBank_aux<KSIZE_1> (current, kmerSize, abundMin, abundMax, KMER_SOLIDITY_SUM, checkSum);
            }
        }

        for (size_t i=0; i<banksComposite.size(); i++)  { banksComposite[i]->forget();  }
    }

    /********************************************************************************/
    void DSK_perBankKmer ()
    {
        /** We can't use even kmer sizes because of palindroms (could have some kmers with coverage 1 instead of 2). */
        DSK_perBankKmer_aux ( 9, 4);
        DSK_perBankKmer_aux (11, 3);
    }

    /********************************************************************************/
    struct DSK_multibank_aux  {  template<typename U> void operator() (U)
    {
        size_t kmerSize = U::value-1;
        size_t nks      = 2;
        string filepath = DBPATH("album.txt");

        /** We configure parameters for a SortingCountAlgorithm object. */
        IProperties* params = SortingCountAlgorithm<>::getDefaultProperties();
        params->setInt (STR_KMER_SIZE,          kmerSize);
        params->setInt (STR_MAX_MEMORY,         MAX_MEMORY);
        params->setInt (STR_KMER_ABUNDANCE_MIN, nks);
        params->setStr (STR_URI_OUTPUT,         "foo");

        IBank* bank = Bank::open(filepath);

        /** We create a DSK instance. */
        SortingCountAlgorithm<U::value> dsk (bank, params);

        /** We launch DSK. */
        dsk.execute();
    }};

    void DSK_multibank ()
    {
        // used to crash on macos because non aligned uint128 in dsk 2nd part

        boost::mpl::for_each<gatb::core::tools::math::IntegerList>(DSK_multibank_aux());
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestDSK);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestDSK);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

