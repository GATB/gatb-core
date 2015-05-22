//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
using namespace std;

/********************************************************************************/
/*                              Sorting count                                   */
/*                                                                              */
/* This snippet shows how to count the common kmers between several input banks */
/*                                                                              */
/* It shows how to inherit from the CountProcessor class, i.e. the class that   */
/* receives kmer counts notifications from the SortingCountAlgorithm.           */
/*                                                                              */
/********************************************************************************/

template<size_t span>
class CountProcessorCustom : public CountProcessorAbstract<span>
{
public:

    // Constructor.
    CountProcessorCustom (
        const vector <CountProcessorHistogram<span>* >& histogramProcessors,
        Partition < HDF5Pair<NativeInt64> >& values
    )
        : _histogramProcessors(histogramProcessors), _values (values) {}

    // Destructor
    virtual ~CountProcessorCustom ()
    {
        for (size_t i=0; i<_histogramProcessors.size(); i++)  {  delete _histogramProcessors[i];  }
    }

    /********************************************************************/
    /*   METHODS CALLED ON THE PROTOTYPE INSTANCE (in the main thread). */
    /********************************************************************/

    // End of the algorithm : we 'end' each histogram
    void end ()   {  for (size_t i=0; i<_histogramProcessors.size(); i++)  {  _histogramProcessors[i]->end();  }  }

    // We need to implement the cloning method for our custom count processor
    CountProcessorAbstract<span>* clone ()
    {
        // We clone each CountProcessorHistogram
        vector <CountProcessorHistogram<span>* > clones;
        for (size_t i=0; i<_histogramProcessors.size(); i++)  {  clones.push_back ((CountProcessorHistogram<span>*)_histogramProcessors[i]->clone());  }

        // We return an instance of our custom processor with the CountProcessorHistogram clones
        return new CountProcessorCustom (clones, _values);
    }

    /********************************************************************/
    /*   METHODS CALLED ON ONE CLONED INSTANCE (in a separate thread).  */
    /********************************************************************/

    // At the end of one partition processing, we get some stats on the N cloned histograms
    void endPart (size_t passId, size_t partId)
    {
        int min_auto_threshold = 3;

        for (size_t i=0; i<_histogramProcessors.size(); i++)
        {
            IHistogram* histo = _histogramProcessors[i]->getHistogram();

            // We compute the threshold
            histo->compute_threshold (min_auto_threshold);

            // We insert the computed solid cutoff and first peak values into the current ith partition
            HDF5Pair<NativeInt64> p (histo->get_solid_cutoff(), histo->get_first_peak());
            _values[i].insert (p);
        }
    }

    // We forward the information to each clone histogram
    virtual bool process (size_t partId, const typename Kmer<span>::Type& kmer, const CountVector& count, CountNumber sum)
    {
        // Note that we provide the 'sum' argument as being the count of the ith bank
        for (size_t i=0; i<_histogramProcessors.size(); i++)  {  _histogramProcessors[i]->process (partId, kmer, count, count[i]);  }
        return true;
    }

protected:

    vector <CountProcessorHistogram<span>* > _histogramProcessors;

    Partition<HDF5Pair<NativeInt64> >& _values;
};

/********************************************************************************/

template<size_t span>  struct MainLoop  {  void operator () (IProperties* options)
{
    // We create a SortingCountAlgorithm instance.
    SortingCountAlgorithm<span> algo (options);

    // We get the number of input banks. For instance, if the uri is "file1.fa,file2.fa", this number will be 2
    size_t nbSources = Bank::getCompositionNb (options->getStr(STR_URI_INPUT));

    // We create a storage for histograms.
    Storage* histogramStorage = StorageFactory(STORAGE_HDF5).create("histograms", true, false);
    LOCAL (histogramStorage);

    Group&  histoGroup = histogramStorage->getGroup("histograms");

    // We create a vector of histogram processors, one for each bank
    vector <CountProcessorHistogram<span>* > histogramProcessors;
    for (size_t i=0; i<nbSources; i++)
    {
        histogramProcessors.push_back (new CountProcessorHistogram<span> (& histoGroup.getGroup(Stringify::format("%d", i))));
    }

    // We create a custom count processor that encapsulates N CountProcessorHistogram instances
    ICountProcessor<span>* processor = new CountProcessorCustom<span> (
        histogramProcessors,
        histogramStorage->root().getPartition<HDF5Pair<NativeInt64> >("values", nbSources)
    );

    // We configure the sorting count algorithm with this custom processor instance
    algo.addProcessor (processor);

    // We launch the SortingCountAlgorithm instance
    algo.execute();

    // We dump some information
    cout << algo.getConfig().getProperties() << endl;

    // We dump some information collected by each CountProcessorHistogram
    for (size_t i=0; i<histogramProcessors.size(); i++)
    {
        cout << histogramProcessors[i]->getProperties() << endl;
    }

    // It is now possible to dump the successive estimations of cutoff and first peak for the Xth bank  (ie replace X by an integer) :
    //  h5dump -y -d values/X histograms.h5 | grep [0-9] | grep -v [A-Z].* | paste - -
}};

/********************************************************************************/
/*                              main function                                   */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We launch our functor
    return Algorithm::mainloop <MainLoop> (SortingCountAlgorithm<>::getOptionsParser(), argc, argv);
}
//! [snippet1]
