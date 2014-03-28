//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

#include <vector>
#include <memory>
using namespace std;

int main (int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cerr << "You must provide a bank URL and a kmer size" << std::endl;
        return EXIT_FAILURE;
    }

    try
    {
        size_t kmerSize = atoi(argv[2]);
        size_t nbKmers  = 1 << (2*kmerSize);

        // We create a kmer model.
        Kmer<>::Model model (kmerSize, KMER_DIRECT);

        // We open the bank.
        BankFasta bank (argv[1]);

        vector<NativeInt64> distrib (nbKmers);
        size_t totalKmers = 0;

        ProgressIterator<Sequence> iter (bank, "iterate bank");
        iter.iterate ([&] (Sequence& seq)
        {
            size_t len  = seq.getDataSize() / kmerSize;
            char*  data = seq.getDataBuffer();

            for (size_t i=0; i<len; i++, data += kmerSize)
            {
                Kmer<>::Type kmer = model.codeSeed (data, seq.getDataEncoding());

                distrib [kmer.toInt()] += 1;
            }
        });

        std::sort (distrib.begin(), distrib.end());

        auto_ptr<Storage> storage (StorageFactory(STORAGE_HDF5).create ("output", false, false));

        Collection<NativeInt64>& distribCollect = storage->getGroup("distrib").getCollection<NativeInt64> (argv[2]);
        distribCollect.addProperty ("name",      argv[1]);
        distribCollect.addProperty ("kmer_size", argv[2]);

        distribCollect.insert (distrib);

        // h5dump -y -w 1 -d distrib/5 -o out output.h5 > /dev/null ; gnuplot -p -e 'plot [][0:] "out" with lines'
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
