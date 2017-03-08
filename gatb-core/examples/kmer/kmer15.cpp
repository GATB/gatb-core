//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <numeric>
using namespace std;

/********************************************************************************/
/*                              Sorting count                                   */
/*                                                                              */
/* This snippet counts specific kmers (core and bound) for each solid partition */
/*                                                                              */
/* a core kmer  := all its 8 neighbors are in the same partition                */
/* a bound kmer := a kmer that is not core.                                     */
/*                                                                              */
/* Cmd-line: kmer15 -in <fasta/q file_1>                                        */
/*                                                                              */
/* Sample: kmer15 -in gatb-core/gatb-core/test/db/reads1.fa                     */
/*                                                                              */
/********************************************************************************/

template<size_t span>
class CountProcessorCustom : public CountProcessorAbstract<span>
{
public:

    typedef Kmer<>::Type                        Type;
    typedef Kmer<>::ModelCanonical              ModelCanon;
    typedef Kmer<>::ModelMinimizer<ModelCanon>  Model;

    // Constructor for the prototype instance
    CountProcessorCustom (SortingCountAlgorithm<span>* algo)
        : _algo(algo), _repart(0), _model(0),_nbIsolated(0)
    {
    }

    // Constructor for the clone instance
    CountProcessorCustom (Repartitor* repart, Model* model)
        : _algo(0), _repart(0), _model(0), _nbIsolated(0)
    {
        setRepart (repart);
        setModel  (model);
    }

    // Destructor
    virtual ~CountProcessorCustom ()
    {
        setRepart (0);
        setModel  (0);
    }

    /********************************************************************/
    /*   METHODS CALLED ON THE PROTOTYPE INSTANCE (in the main thread). */
    /********************************************************************/

    // We need to implement the cloning method for our custom count processor
    CountProcessorAbstract<span>* clone ()  {  return new CountProcessorCustom (_repart, _model);  }

    //
    void begin (const Configuration& config)
    {
        setRepart (_algo->getRepartitor());
        setModel  (new Model (config._kmerSize, config._minim_size));
    }

    void end ()
    {
        printf ("FINAL : nbIsolated = %ld\n", _nbIsolated);
    }

    //
    void finishClones (std::vector<ICountProcessor<span>*>& clones)
    {
        for (size_t i=0; i<clones.size(); i++)
        {
            if (CountProcessorCustom* clone = dynamic_cast<CountProcessorCustom*> (clones[i]))
            {
                _nbIsolated += clone->_nbIsolated;
            }
        }
    }

    /********************************************************************/
    /*   METHODS CALLED ON ONE CLONED INSTANCE (in a separate thread).  */
    /********************************************************************/

    // We refine the 'process' method
    virtual bool process (size_t partId, const Type& kmer, const CountVector& count, CountNumber sum)
    {
if (count[0] < 3)  { return false; }

        bool inCore = true;

        // We iterate the neighbors of the current kmer
        _model->iterateNeighbors (kmer, [&] (const Type& neighbor)
        {
            // We get the minimizer value of the current kmer.
            u_int64_t mini = _model->getMinimizerValue(neighbor);

            // We get the partition index of the neighbor from its minimizer value and the minimizer repartition table.
            u_int64_t neighborPart = (*_repart) (mini);

            // We look whether or not the neighbor is in the same partition as the current one
            if (neighborPart != partId)  { inCore = false; }
        });

        if (inCore)  { _core.insert  (kmer);  }
        else         { _bound.insert (kmer);  }

        return true;
    }

    void endPart (size_t passId, size_t partId)
    {
        for (set<Type>::iterator itCore = _core.begin(); itCore != _core.end(); ++itCore)
        {
            size_t nbNeighbors = 0;

            // We iterate the neighbors of the current kmer
            _model->iterateNeighbors (*itCore, [&] (const Type& neighbor)
            {
                if (_core.find(neighbor) != _core.end())
                {
                    nbNeighbors ++;
                }
                else if (_bound.find(neighbor) != _bound.end())
                {
                    nbNeighbors ++;
                }
            });

            if (nbNeighbors == 0)  { _nbIsolated++; }
        }
    }


    /*****************************************************************/
    /*                          MISCELLANEOUS.                       */
    /*****************************************************************/

private:

    SortingCountAlgorithm<span>* _algo;

    Repartitor* _repart;
    void setRepart (Repartitor* repart)  { SP_SETATTR(repart); }

    Model*      _model;
    void setModel (Model* model) { SP_SETATTR(model); }

    set<Type> _core;
    set<Type> _bound;

    u_int64_t _nbIsolated;
};

/********************************************************************************/

template<size_t span>  struct MainLoop  {  void operator () (IProperties* options)
{
    // We create a SortingCountAlgorithm instance.
    SortingCountAlgorithm<span> algo (options);

    // We create a custom count processor and give it to the sorting count algorithm
    CountProcessorCustom<span>* processor = new CountProcessorCustom<span> (&algo);
    algo.addProcessor(processor);

    // We launch the algorithm
    algo.execute();

//    double total = processor->getNbKmersTotal();
//    printf ("total : %9ld\n", processor->getNbKmersTotal());
//    printf ("core  : %9ld  %5.1f \n", processor->getNbKmersCore(),  100.0 * (double) processor->getNbKmersCore()  / total);
//    printf ("bound : %9ld  %5.1f \n", processor->getNbKmersBound(), 100.0 * (double) processor->getNbKmersBound() / total);
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
