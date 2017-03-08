/* 
 * given N datasets
 * outputs a matrix containing the following information:
 *
 * [kmer] [abundance in dataset 1] [abundance in dataset 2] etc..
 *
 * for all input kmers in the union of the datasets
 *
 * Cmd-line: kmer16 -in <fasta/q file_1>
 *
 * Sample: kmer16 -in gatb-core/gatb-core/test/db/reads1.fa
 *
 */
#include <gatb/gatb_core.hpp>
using namespace std;

/********************************************************************************/
template<size_t span>
class CountProcessorCustom : public CountProcessorAbstract<span>
{
public:

    // We need the kmer size to dump kmers values as nucl strings
    CountProcessorCustom (size_t kmerSize, ISynchronizer* synchro) : kmerSize(kmerSize), synchro(synchro) {}

    virtual ~CountProcessorCustom () {}

    CountProcessorAbstract<span>* clone ()  { return new CountProcessorCustom<span>(kmerSize, synchro); }

    virtual bool process (size_t partId, const typename Kmer<span>::Type& kmer, const CountVector& count, CountNumber sum)
    {
        // get a mutex
        LocalSynchronizer sync (synchro);
        // output kmer and counts
        cout << kmer.toString(kmerSize) << " ";
        for (size_t i=0; i<count.size(); i++)  {  cout << count[i] << " ";  }
        cout  << endl;
        return true;
    }

private:
    size_t kmerSize;
    ISynchronizer *synchro;
};

/********************************************************************************/
template<size_t span>  struct MainLoop  {  void operator () (IProperties* options)
{
    // We force the solidity kind (otherwise default value "sum" will consider N banks as a single one)
    options->setStr(STR_SOLIDITY_KIND, "all");

    // We create a SortingCountAlgorithm instance.
    SortingCountAlgorithm<span> algo (options);

    // global synchronization
    ISynchronizer* synchro = System::thread().newSynchronizer();

    // We create a custom count processor and give it to the sorting count algorithm
    algo.addProcessor (new CountProcessorCustom<span> (options->getInt(STR_KMER_SIZE), synchro));

    // We launch the algorithm
    algo.execute();
}};

/********************************************************************************/
int main (int argc, char* argv[])
{
    // We create a command line parser for the sorting count algorithm
    IOptionsParser* parser = SortingCountAlgorithm<>::getOptionsParser ();
    parser->push_back (new OptionOneParam (STR_NB_CORES, "nb cores",  false, "1"));
    parser->push_back (new OptionOneParam (STR_VERBOSE,  "verbosity", false, "1"));

    // We launch our functor
    return Algorithm::mainloop <MainLoop> (parser, argc, argv);
}


