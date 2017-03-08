//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>
#include <fstream>
#include <string.h>

/********************************************************************************/
/*                              Minimizers                                      */
/*                                                                              */
/* This snippet computes the total number of minimizers with a parallel implem  */
/*          with the use of the dispatcher                                      */
/*                                                                              */
/* Cmd-line: kmer17 -in <fasta/q file> -kmer-size <value> ...                   */
/*                  -minimizer-size <size>                                      */
/*                                                                              */
/* Sample: kmer17 -in gatb-core/gatb-core/test/db/reads1.fa -kmer-size 11 ...   */
/*               -minimizer-size 11                                             */
/*                                                                              */
/********************************************************************************/

// We use the required packages
using namespace std;

/** Kmer span definition. */
const size_t span = KMER_SPAN(0);

/** Some shortcuts. */
typedef Kmer<span>::ModelDirect    ModelDirect;
typedef Kmer<span>::ModelCanonical    ModelCanonical;

//typedef Kmer<span>::ModelMinimizer<ModelDirect> ModelMinimizer;
typedef Kmer<span>::ModelMinimizer<ModelCanonical> ModelMinimizer;

class seqProcessor
{
public:
	seqProcessor(int kmerSize,int mmerSize,u_int64_t * nb_kmers_global,u_int64_t * nb_minim_global) : _nb_kmers_global(nb_kmers_global), _nb_minim_global(nb_minim_global), _kmerSize(kmerSize), _mmerSize(mmerSize)
	{
		// We declare a kmer model and a minimizer model
		_nb_minim_total=0;
		_nb_kmers_total=0;
		_model = new ModelMinimizer(_kmerSize, _mmerSize);
	}
	
	seqProcessor(seqProcessor const& ref)
	{
		_nb_minim_total=0;
		_nb_kmers_total=0;
		_nb_minim_global=ref._nb_minim_global;
		_nb_kmers_global=ref._nb_kmers_global;
		_mmerSize=ref._mmerSize;
		_kmerSize=ref._kmerSize;
		
		_model = new ModelMinimizer(_kmerSize, _mmerSize);
	}
	
	void operator() (Sequence& seq)
	{
		
		// We iterate the kmers (and minimizers) of the current sequence.
		_model->iterate (seq.getData(), [&] (const ModelMinimizer::Kmer& kmer, size_t idx)
						 {
							 
							 // We may have to update the number of different minimizer in the current sequence
							 if (kmer.hasChanged()==true)  {  _nb_minim_total++ ;   }
							 
							 _nb_kmers_total ++;
						 });
	}
	
	// increment the global counter in a thread safe way in the destructor
	~seqProcessor()
	{
		__sync_fetch_and_add(_nb_kmers_global,_nb_kmers_total);
		__sync_fetch_and_add(_nb_minim_global,_nb_minim_total);
	}
	
private:
	ModelMinimizer *_model;
	u_int64_t _nb_minim_total;
	u_int64_t * _nb_kmers_global;
	u_int64_t * _nb_minim_global;
	int _kmerSize, _mmerSize;
	
	u_int64_t _nb_kmers_total;
};

/********************************************************************************/
int main (int argc, char* argv[])
{
	static const char* STR_BATCH_SIZE = "-batch";
	
	/** We create a command line parser. */
	OptionsParser parser ("KmerTest");
	parser.push_back (new OptionOneParam (STR_URI_INPUT,      "bank input",     true));
	parser.push_back (new OptionOneParam (STR_KMER_SIZE,      "kmer size",      true));
	parser.push_back (new OptionOneParam (STR_MINIMIZER_SIZE, "minimizer size", true));
	parser.push_back (new OptionOneParam (STR_NB_CORES,   "number of cores",     false, "0"));
	parser.push_back (new OptionOneParam (STR_BATCH_SIZE,   "batch size",     false, "1000"));
	
	try
	{
		/** We parse the user options. */
		IProperties* options = parser.parse (argc, argv);
		
		int nbcores = options->getInt(STR_NB_CORES);
		int batch_size = options->getInt(STR_BATCH_SIZE);
		
		IDispatcher* _dispatcher  = new Dispatcher (nbcores);
		
		string bankFilename = options->getStr(STR_URI_INPUT);
		
		// We get the kmer and minimizer sizes.
		size_t kmerSize = options->getInt(STR_KMER_SIZE);
		size_t mmerSize = options->getInt(STR_MINIMIZER_SIZE);
		
		u_int64_t nbKmers       = 0;
		u_int64_t nbMinimizers  = 0;
		
		// We declare a bank instance defined by a list of filenames
		IBank* bank = Bank::open (bankFilename);
		
		// We create an iterator over this bank.
		ProgressIterator<Sequence> itSeq (*bank);
		
		_dispatcher->iterate(itSeq,seqProcessor(kmerSize,mmerSize,&nbKmers,&nbMinimizers),batch_size);
		
		printf("nb cores %i batch size %i \n",nbcores,batch_size);
		printf("nb minimizer total : %llu\n",nbMinimizers);
		printf("nb kmers total : %llu\n",nbKmers);
		
	}
	
	catch (OptionFailure& e)
	{
		return e.displayErrors (std::cout);
	}
	catch (Exception& e)
	{
		std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
	}
	
	return EXIT_SUCCESS;
}
//! [snippet1]
