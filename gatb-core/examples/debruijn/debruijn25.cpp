//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

// We use the required packages
using namespace std;

/********************************************************************************/
/********************************************************************************/
int main (int argc, char* argv[])
{
    const size_t SPAN = KMER_SPAN(1);

    /** Shortcuts. */
    typedef Kmer<SPAN>::Type  Type;
    typedef Kmer<SPAN>::Count Count;
    typedef Kmer<SPAN>::ModelCanonical ModelCanon;
    typedef Kmer<SPAN>::ModelMinimizer <ModelCanon> Model;

    size_t kmerSize = 33;
    size_t mmerSize = 11;

    /** We create a command line parser. */
    OptionsParser parser ("GraphStats");
    parser.push_back (new OptionOneParam (STR_URI_INPUT, "bank input",  true));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        string filename = options->getStr (STR_URI_INPUT);

        /** We create the solid kmers. */
        Graph graph = Graph::create ("-in %s -kmer-size %d  -bloom none -out toto.h5  -abundance-min 1", filename.c_str(), kmerSize);

        /** We get the information of the solid kmers from the HDF5 file. */
        Storage* storage = StorageFactory(STORAGE_HDF5).load ("toto.h5");   LOCAL (storage);
        Group& dskGroup = storage->getGroup("dsk");

        /** We get the solid kmers partition. */
        Partition<Count>& partition = dskGroup.getPartition<Count> ("solid");

        /** We create two kmers models. */
        Model  model   (kmerSize,   mmerSize);
        Model  modelK1 (kmerSize-1, mmerSize);

        // We declare an output Bank
        BankBinary outputBank (System::file().getBaseName(filename) + ".bin");

        /** We create a sequence with BINARY encoding. */
        Sequence seq (Data::ASCII);

        /** We get an iterator over the [kmer,abundance] of solid kmers. */
        Iterator<Count>* it = partition.iterator();   LOCAL (it);

        /** We iterate the solid kmers. */
        for (it->first(); !it->isDone(); it->next())
        {
            Type current = it->item().value;

            cout << "kmer=" << it->item().value << "  minimizer=" << model.getMinimizerValue(current)
                << "  abundance=" << it->item().abundance << endl;

            /** We interpret the kmer value as a Data object. */
            Data data (Data::BINARY);
            data.setRef ((char*) &current, model.getKmerSize());

            modelK1.iterate (data, [&] (const Model::Kmer& k, size_t idx)
            {
                /** Shortcut. */
                Type miniminizerCurrent = k.minimizer().value();

                cout << "-> " << k.value()
                     << " minimizer=" << miniminizerCurrent << " "
                     << modelK1.getMmersModel().toString (miniminizerCurrent)
                     << endl;

                string tmp = modelK1.getMmersModel().toString (miniminizerCurrent);

                /** We interpret the minimizer value as a Data object. */
                seq.getData().setRef ((char*)tmp.c_str(), modelK1.getMmersModel().getKmerSize());

                /** We insert the sequence into the binary bank. */
                outputBank.insert (seq);
            });
        }

        /** We flush the output bank. */
        outputBank.flush();

        /** We iterate the output bank. */
        outputBank.iterate ([&] (const Sequence& s)
        {
            /** We get the kmer corresponding to the current sequence. */
            ModelCanon::Kmer mini = modelK1.getMmersModel().codeSeed (s.getDataBuffer(), Data::BINARY);

            cout << "mini=" << mini.value() << "  " << modelK1.getMmersModel().toString (mini.value()) << endl;
        });
    }
    catch (OptionFailure& e)
    {
        return e.displayErrors (std::cout);
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }
}
//! [snippet1]
