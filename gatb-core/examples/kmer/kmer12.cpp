//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <numeric>
using namespace std;

/********************************************************************************/
/*                              Sorting count                                   */
/*                                                                              */
/* This snippet shows how to count the common kmers between several input banks */
/*                                                                              */
/* It shows how to inherit from the CountProcessor class, i.e. the class that   */
/* receives kmer counts notifications from the SortingCountAlgorithm.           */
/*                                                                              */
/* Cmd-line: kmer12 -in <fasta/q file_1>[,<fasta/q file_2>,...]                 */
/*                                                                              */
/* Sample: kmer12 -in gatb-core/gatb-core/test/db/reads1.fa,gatb-core/gatb-core/test/db/reads2.fa  */
/*                                                                              */
/********************************************************************************/

// We define our own CountProcessor class. Here, we want to count the number of
// common kmers between each pair of input banks.
//
// We refines the 'process' method : from the vector of counts (one count per bank)
// for the current kmer, we check for all possible pairs that the two counts are greater
// than some threshold. In such a case, we increment the number of common kmers for
// the current pair.
//
// If we have N input banks, we will have N.(N+1)/2 possible pairs. For storing the number
// of common kmers for each pair, we use a vector of size N.(N+1)/2 and use a 'diagonal'
// representation with the 'offset' method. For instance, for N=4, we have:
//      0   4   7   9
//          1   5   8       Here, row and column are the indexes of the banks.
//              2   6       For instance, offset(1,2)=5
//                  3       The main diagonal corresponds to the number of kmers in each bank
//
// At the end of the sorting count, one can get the number of common kmers between the
// ith and jth banks with the method 'getCount(i,j)'
//
template<size_t span>
class CountProcessorCustom : public CountProcessorAbstract<span>
{
public:

    // Constructor.
    CountProcessorCustom (size_t nbBanks, const CountRange& range)  : _nbBanks(nbBanks), _range(range)
    {
        // We configure the vector for the N.(N+1)/2 possible pairs
        _countTotal.resize (_nbBanks*(_nbBanks+1)/2);
    }

    // Destructor
    virtual ~CountProcessorCustom () {}

    /********************************************************************/
    /*   METHODS CALLED ON THE PROTOTYPE INSTANCE (in the main thread). */
    /********************************************************************/

    // We need to implement the cloning method for our custom count processor
    CountProcessorAbstract<span>* clone ()  {  return new CountProcessorCustom (_nbBanks, _range);  }

    // When several cloned CountProcessor instances have finished to process their partition,
    // we need to update the main CountProcessor instance with the collected information from the clones.
    void finishClones (std::vector<ICountProcessor<span>*>& clones)
    {
        for (size_t i=0; i<clones.size(); i++)
        {
            /** We have to recover type information. */
            if (CountProcessorCustom* clone = dynamic_cast<CountProcessorCustom*> (clones[i]))
            {
                for (size_t i=0; i<this->_countTotal.size(); i++)  { this->_countTotal[i] += clone->_countTotal[i];  }
            }
        }
    }

    /********************************************************************/
    /*   METHODS CALLED ON ONE CLONED INSTANCE (in a separate thread).  */
    /********************************************************************/

    // We refine the 'process' method
    virtual bool process (size_t partId, const typename Kmer<span>::Type& kmer, const CountVector& count, CountNumber sum)
    {
        // We try each possible pair (i,j) through the 'i' and 'j' loops
        for (size_t i=0; i<count.size(); i++)
        {
            if (_range.includes(count[i]) == true)
            {
                for (size_t j=i; j<count.size(); j++)
                {
                    if (_range.includes(count[j]) == true)
                    {
                        // Both 'i' and 'j' banks hold the current kmer with an abundance in the correct range
                        // so we increase the number of common kmers for this pair (i,j)
                        _countTotal [offset(i,j)] ++;
                    }
                }
            }
        }

        return true;
    }

    /*****************************************************************/
    /*                          MISCELLANEOUS.                       */
    /*****************************************************************/

    // We can return the total number of common kmers in each pair of banks
    u_int64_t getNbItems ()  {  return std::accumulate (_countTotal.begin(), _countTotal.end(), 0);  }

    // After execution, this method returns the number of common kmers between banks 'i' and 'j'
    size_t getCount (size_t i, size_t j) const { return _countTotal[offset(i,j)]; }

protected:

    // Get a unique integer from the pair (i,j)
    size_t offset (size_t i, size_t j) const { return _nbBanks*(j-i) - ((j-i)*((j-i)-1))/2 + i; }

private:

    size_t         _nbBanks;
    vector<size_t> _countTotal;
    CountRange     _range;
};

/********************************************************************************/

template<size_t span>  struct MainLoop  {  void operator () (IProperties* options)
{
    // We get the min abundance threshold
    size_t nks = options->getInt(STR_KMER_ABUNDANCE_MIN);

    // We get the number of input banks.
    // For instance, if the uri is "file1.fa,file2.fa", this number will be 2
    size_t nbSources = Bank::getCompositionNb (options->getStr(STR_URI_INPUT));

    // We create a SortingCountAlgorithm instance.
    SortingCountAlgorithm<span> algo (options);

    // We create a custom count processor and give it to the sorting count algorithm
    CountProcessorCustom<span>* processor = new CountProcessorCustom<span> (nbSources, CountRange(nks,1<<30));
    algo.addProcessor(processor);

    // We launch the algorithm
    algo.execute();

    // We first dump the number of kmers in each bank.
    printf ("\n");
    for (size_t i=0; i<nbSources; i++)  {  printf ("[%2d] %9d \n", i, processor->getCount(i,i));  }
    printf ("\n");

    // We dump statistics about each pair (i,j)
    for (size_t i=0; i<nbSources; i++)
    {
        for (size_t j=i+1; j<nbSources; j++)
        {
            printf ("[%2d,%2d] %9d  %5.1f  %5.1f \n", i, j,
                processor->getCount(i,j),
                100.0*(double)processor->getCount(i,j) / (double) processor->getCount(i,i),
                100.0*(double)processor->getCount(i,j) / (double) processor->getCount(j,j)
            );
        }
    }
}};

/********************************************************************************/
/*                              main function                                   */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We create a command line parser for the sorting count algorithm
    IOptionsParser* parser = SortingCountAlgorithm<>::getOptionsParser ();

    parser->push_back (new OptionOneParam (STR_NB_CORES, "nb cores",  false, "0"));
    parser->push_back (new OptionOneParam (STR_VERBOSE,  "verbosity", false, "1"));

    // We launch our functor
    return Algorithm::mainloop <MainLoop> (parser, argc, argv);
}
//! [snippet1]
