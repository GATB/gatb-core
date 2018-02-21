/* 
 * KAT-like analysis of a genome assembly
 *
 * given 2 datasets
 * outputs the 2D histogram of kmers
 *
 * must give 2 datasets, first readset, then genome
 * Sample command: histo2D -in readset.fasta,genome.fasta -kmer-size 25 > resultfile
 *
 */
#include <gatb/gatb_core.hpp>
using namespace std;

int dim1 = 10000; // max kmer occurence in read set
int dim2 = 10;    // max kmer occurence in genome

#define IDX(i,j) ((i) + (j)*dim1)

/********************************************************************************/
template<size_t span>
class CountProcessorCustom : public CountProcessorAbstract<span>
{
public:

    // We need the kmer size to dump kmers values as nucl strings
    CountProcessorCustom (size_t kmerSize, ISynchronizer* synchro, u_int64_t * histo_table_2D) : kmerSize(kmerSize), synchro(synchro), _histo_table_2D(histo_table_2D) {}

    virtual ~CountProcessorCustom () {}

    CountProcessorAbstract<span>* clone ()  { return new CountProcessorCustom<span>(kmerSize, synchro,_histo_table_2D); }

    virtual bool process (size_t partId, const typename Kmer<span>::Type& kmer, const CountVector& count, CountNumber sum)
    {
		if( count[0]  < dim1  &&   count[1]  < dim2)
			__sync_fetch_and_add ( & _histo_table_2D[IDX( count[0] , count[1] )], 1);
        return true;
    }

private:
    size_t kmerSize;
    ISynchronizer *synchro;
	u_int64_t * _histo_table_2D;
	
};

/********************************************************************************/
template<size_t span>  struct MainLoop  {  void operator () (IProperties* options)
{

	//check that the user gave 2 input banks
	std::string banknames = options->getStr(STR_URI_INPUT);
	size_t n = std::count(banknames.begin(), banknames.end(), ',');
	if( (n+1) != 2)
	{
		printf("There must be 2 input banks.\n");
		exit(1);
	}
	
	u_int64_t * histo_table_2D = (u_int64_t *) calloc( dim1 *  dim2 , sizeof(u_int64_t *));
	
    // We force the solidity kind (otherwise default value "sum" will consider N banks as a single one)
    options->setStr(STR_SOLIDITY_KIND, "all");

    // We create a SortingCountAlgorithm instance.
    SortingCountAlgorithm<span> algo (options);

    // global synchronization
    ISynchronizer* synchro = System::thread().newSynchronizer();

    // We create a custom count processor and give it to the sorting count algorithm
    algo.addProcessor (new CountProcessorCustom<span> (options->getInt(STR_KMER_SIZE), synchro,histo_table_2D));

    // We launch the algorithm
    algo.execute();

	//outptut the 2D histogram
	for(int ii=0; ii< dim1; ii++)
	{
		printf("%5i:\t",ii);
		for(int jj=0; jj< dim2; jj++)
		{
			printf("\t%6lli",histo_table_2D[IDX(ii,jj)]);
		}
		printf("\n");
	}
	free(histo_table_2D);
	
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
