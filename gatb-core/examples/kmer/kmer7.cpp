//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

#include <vector>
#include <memory>
using namespace std;

/********************************************************************************/
/*                              Kmer statistics                                 */
/*                                                                              */
/* This snippet is a little bit more complex and uses different features.       */
/*                                                                              */
/* Its purpose is to get the kmers distribution of a bank, for a given kmer     */
/* size K. Each sequence of size N is split in pieces of size K, for instance:  */
/*     ATCTCGGGCTAGCTCTCGATAAGC => for K=3, ATC TCG GGC TAG CTC TCG ATA AGC     */
/* Then we count the number of occurrences of these pieces, and sort the result.*/
/* Finally, this distribution is saved in a Storage object (HDF5 format here).  */
/* The resulting data can be shown with HDF5 tools (provided with GATB) by:     */
/*      h5dump output.h5                                                        */
/* You can directly get a gnuplot output with (4 here is the provided kmer size)*/
/*                                                                              */
/*      h5dump -y -w 1 -d distrib/4 -o out output.h5 > /dev/null                */
/*      gnuplot -p -e 'plot [][0:] "out" with lines'                            */
/*                                                                              */
/* NOTE: don't use too big kmer size because of the potential huge memory usage.*/
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("KmerStats");
    parser.push_back (new OptionOneParam (STR_URI_INPUT,  "bank reference", true));
    parser.push_back (new OptionOneParam (STR_KMER_SIZE,  "kmer size",      true));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        size_t kmerSize = options->getInt (STR_KMER_SIZE);
        size_t nbKmers  = 1 << (2*kmerSize);

        // We create a kmer model.
        Kmer<>::ModelDirect model (kmerSize);

        // We open the bank.
        IBank* bank = Bank::open (options->getStr(STR_URI_INPUT));
        LOCAL (bank);

        vector<NativeInt64> distrib (nbKmers);
        size_t totalKmers = 0;

        // We define an iterator that encapsulates the sequences iterator with progress feedback
        ProgressIterator<Sequence> iter (*bank, "iterate bank");

        // We loop the bank
        for (iter.first(); !iter.isDone(); iter.next())
        {
            // Shortcut
            Sequence& seq = iter.item();

            size_t len  = seq.getDataSize() / kmerSize;
            char*  data = seq.getDataBuffer();

            // We iterate the sequence data by block of size kmerSize
            for (size_t i=0; i<len; i++, data += kmerSize)
            {
                // We get the kmer value of the current block
                Kmer<>::ModelDirect::Kmer kmer = model.codeSeed (data, seq.getDataEncoding());

                // We update the occurrences number for this kmer value
                distrib [kmer.value().toInt()] += 1;
            }
        }

        // We sort the distribution
        std::sort (distrib.begin(), distrib.end());

        // We create the output storage object
        auto_ptr<Storage> storage (StorageFactory(STORAGE_HDF5).create ("output", true, false));

        // We create a data set in our storage object
        Collection<NativeInt64>& distribCollect = storage->getGroup("distrib").getCollection<NativeInt64> (options->getStr (STR_KMER_SIZE));

        // We tag our data set with the user parameters
        distribCollect.addProperty ("name",      options->getStr (STR_URI_INPUT));
        distribCollect.addProperty ("kmer_size", options->getStr (STR_KMER_SIZE));

        // We insert the distribution into our storage object
        distribCollect.insert (distrib);

        // We can get a gnuplot output with the following (for kmerSize=5)
        //      h5dump -y -w 1 -d distrib/5 -o out output.h5 > /dev/null ; gnuplot -p -e 'plot [][0:] "out" with lines'
    }
    catch (OptionFailure& e)
    {
        return e.displayErrors (cout);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
