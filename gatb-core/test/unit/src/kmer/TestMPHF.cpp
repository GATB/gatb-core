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

#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankStrings.hpp>

#include <gatb/kmer/impl/SortingCountAlgorithm.hpp>
#include <gatb/kmer/impl/MPHFAlgorithm.hpp>

#include <gatb/tools/misc/api/Macros.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <gatb/tools/math/LargeInt.hpp>

#include <gatb/tools/storage/impl/Storage.hpp>

#include <gatb/tools/collections/impl/BooPHF.hpp>

using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::dp;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::storage;
using namespace gatb::core::tools::storage::impl;

using namespace gatb::core::tools::math;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for genomic databases management
 */
class TestMPHF : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestMPHF);

        CPPUNIT_TEST_GATB (MPHF_check1);
        CPPUNIT_TEST_GATB (MPHF_check2);

        // no mphf1 anymore
        CPPUNIT_TEST_GATB (test_mphf2);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    () {}
    void tearDown () {}

        // SMALL VALUE NEEDED because continuous integration servers are not very powerful...
        static const u_int64_t MAX_MEMORY = 1000;


    /** Shortcuts. */
    typedef Kmer<32>::Count Count;
    typedef Kmer<32>::Type  Type;

    /********************************************************************************/
    void MPHF_check1 ()
    {
        size_t kmerSize = 11;
        size_t nks      = 1;

        const char* seqs[] = {
            "CGCTACAGCAGCTAGTTCATCATTGTTTATCAATGATAAAATATAATAAGCTAAAAGGAAACTATAAATA"
            "ACCATGTATAATTATAAGTAGGTACCTATTTTTTTATTTTAAACTGAAATTCAATATTATATAGGCAAAG"
        } ;

        /** We configure parameters for a SortingCountAlgorithm object. */
        IProperties* params = SortingCountAlgorithm<>::getDefaultProperties();
        params->setInt (STR_KMER_SIZE,          kmerSize);
        params->setInt (STR_MAX_MEMORY,         MAX_MEMORY);
        params->setInt (STR_KMER_ABUNDANCE_MIN, nks);
        params->setStr (STR_URI_OUTPUT,         "foo");

        /** We create a DSK instance. */
        SortingCountAlgorithm<> sortingCount (new BankStrings (seqs, ARRAY_SIZE(seqs)), params);

        /** We launch DSK. */
        sortingCount.execute();

        if (sortingCount.getSolidCounts()->getNbItems() != (int)(strlen(seqs[0]) - kmerSize + 1))
            std::cout << "problem with sortingcount nb items: " << sortingCount.getSolidCounts()->getNbItems()  << " != " << (int)(strlen(seqs[0]) - kmerSize + 1) << std::endl;


        CPPUNIT_ASSERT (sortingCount.getSolidCounts()->getNbItems() == (int)(strlen(seqs[0]) - kmerSize + 1) );

        /** We get the storage instance. */
        Storage* storage = sortingCount.getStorage();


        /** We create a mphf instance. */
        MPHFAlgorithm<> mphf (storage->getGroup("dsk"), "mphf", sortingCount.getSolidCounts(), sortingCount.getSolidKmers(), 1, true);

        /** We actually execute the mphf construction. */
        mphf.execute();
            

        if (mphf.getAbundanceMap() == 0)
            std::cout << "could not get abundance map" << std::endl;

        CPPUNIT_ASSERT (mphf.getAbundanceMap() != 0);

        MPHFAlgorithm<>::AbundanceMap& theMap = * mphf.getAbundanceMap();

        // below are quick tests

        if (theMap.size() != 130)
            std::cout << "incorrect map size:" << theMap.size() << " != 130" << std::endl;

        CPPUNIT_ASSERT (theMap.size() == 130);

        typedef /*typename*/ Kmer<32>::ModelCanonical       Model;
        typedef /*typename*/ Kmer<32>::ModelCanonical::Kmer Kmer;
        Model model (11);

        Kmer kmer = model.codeSeed ("ACCATGTATAA", Data::ASCII);

		
        theMap.at(kmer.value()) = 4;

        if (theMap.at(kmer.value()) != 4)
            std::cout << "bad map value " << theMap[kmer.value()]  << " != 4"  << std::endl;
        CPPUNIT_ASSERT (theMap.at(kmer.value()) == 4);
    }

    /********************************************************************************/
    void MPHF_check2 ()
    {
        /** We define our MPHF type for kmers. */
        typedef BooPHF<Type> MPHF;

        size_t kmerSize = 11;
        size_t nks      = 1;

        const char* seqs[] = {
            "CGCTACAGCAGCTAGTTCATCATTGTTTATCAATGATAAAATATAATAAGCTAAAAGGAAACTATAAATA"
            "ACCATGTATAATTATAAGTAGGTACCTATTTTTTTATTTTAAACTGAAATTCAATATTATATAGGCAAAG"
        } ;

        //////////////////////////////////////////////////
        // PART 1 : we get solid kmers from SortingCount
        //////////////////////////////////////////////////

        /** We configure parameters for a SortingCountAlgorithm object. */
        IProperties* params = SortingCountAlgorithm<>::getDefaultProperties();
        params->setInt (STR_KMER_SIZE,          kmerSize);
        params->setInt (STR_KMER_ABUNDANCE_MIN, nks);
        params->setInt (STR_MAX_MEMORY,         MAX_MEMORY);
        params->setStr (STR_URI_OUTPUT,         "foo");

        IBank* bank = new BankStrings (seqs, ARRAY_SIZE(seqs));
        LOCAL (bank);

        /** We create a DSK instance. */
        SortingCountAlgorithm<> sortingCount (bank, params);

        /** We launch DSK. */
        sortingCount.execute();

        // Shortcut
        Iterable<Type>* solids = sortingCount.getSolidKmers();
        LOCAL (solids);

        size_t nbSolids = solids->getNbItems();

        //////////////////////////////////////////////////
        // PART 2 : we build a hash and save it
        //////////////////////////////////////////////////
        TimeInfo ti;

        /** We create a hash for the result of DSK. */
        MPHF hash1;

        /** We build the hash function. */
        ti.start("build");
        hash1.build (solids);
        ti.stop("build");

        CPPUNIT_ASSERT (hash1.size() == nbSolids);

        // We save the hash object in the dedicated storage group.
        ti.start("save");
        hash1.save (sortingCount.getStorage()->getGroup("dsk"), "mphf");
        ti.stop("save");

        // We need a vector to check codes existence
        vector<bool> check (nbSolids);
        CPPUNIT_ASSERT (check.size()==nbSolids);

        // We get an iterator on the solid kmers.
        Iterator<Type>* itSolids = solids->iterator();
        LOCAL (itSolids);

        // We loop over the solid kmers
        for (itSolids->first(); !itSolids->isDone(); itSolids->next())
        {
            // We get the hash code for the current key
            MPHF::Code code = hash1 (itSolids->item());

            // cout << "KEY=" << itSolids->item().value << "  VALUE=" << code << endl;

            // We check we never saw that code before
            CPPUNIT_ASSERT (check[code]==false);
            check[code]=true;
            CPPUNIT_ASSERT (check[code]==true);
        }

        // We check that all codes have been seen
        for (size_t i=0; i<check.size(); i++)  { CPPUNIT_ASSERT(check[i]==true); }

        typedef Kmer<32>::ModelCanonical       Model;
        typedef Kmer<32>::ModelCanonical::Kmer Kmer;

        // We create a kmer
        Model model (kmerSize);
        Kmer kmer = model.codeSeed ("ACCATGTATAA", Data::ASCII);

        // We check that it is known by the hash function
        hash1 (kmer.value());

        //////////////////////////////////////////////////
        // PART 3 : we read the hash from a storage
        //////////////////////////////////////////////////
        MPHF hash2;
        hash2.load (sortingCount.getStorage()->getGroup("dsk"), "mphf");

        CPPUNIT_ASSERT (hash1.size() == hash2.size());

        // We check that both hash instances have the same content.
        for (itSolids->first(); !itSolids->isDone(); itSolids->next())
        {
            // We get the hash code for the current key
           CPPUNIT_ASSERT (hash1 (itSolids->item()) == hash2 (itSolids->item()));
        }
    }


    /********************************************************************************/
    void test_mphf2 (void)
    {
        /** Shortcuts. */
        typedef int Key;
        typedef BooPHF<Key>   Hash;
        typedef Hash::Code HashValue;

            // We create a list of keys.
            Key values[] = {1,2,3,5,8,13,21,34,55,89};
            std::list<Key> l (values, values + sizeof(values)/sizeof(values[0]) );

            // We create a bag from this list
            BagFile<Key> bagFile ("./keys");
            for (list<Key>::iterator it=l.begin(); it!=l.end(); ++it)  {  bagFile.insert (*it);  }
            bagFile.flush();

            // We create an iterable from the file
            IterableFile<Key> iterableFile ("./keys");

            CPPUNIT_ASSERT (iterableFile.getNbItems() == (int)l.size());

            // We create our hash function object
            Hash hash;

            // We build the hash function with the keys in the iterator
            hash.build (&iterableFile);

            CPPUNIT_ASSERT (hash.size() == l.size());

            // We need a vector to check codes existence
            vector<bool> check (l.size());

            // We dump each hash code for each key
            Iterator<Key>* itKeys = iterableFile.iterator();  LOCAL (itKeys);

            for (itKeys->first(); !itKeys->isDone(); itKeys->next())
            {
                // We get the hash code for the current key
                HashValue code = hash (itKeys->item());

                // We check we never saw that code before
                CPPUNIT_ASSERT (check[code]==false);
                check[code]=true;
                CPPUNIT_ASSERT (check[code]==true);
            }

            // We check that all codes have been seen
            for (size_t i=0; i<check.size(); i++)  { CPPUNIT_ASSERT(check[i]==true); }
        }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestMPHF);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestMPHF);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/
