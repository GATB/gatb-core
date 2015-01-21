//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <list>
#include <iostream>

/********************************************************************************/
/*                Command line parsing through IOptionParser                    */
/********************************************************************************/
int main (int argc, char* argv[])
{
    IOptionsParser* root = new OptionsParser ("root");

    root->push_back (SortingCountAlgorithm<>::getOptionsParser(false));
    root->push_back (DebloomAlgorithm<>::getOptionsParser());

    IOptionsParser* child1 = new OptionsParser ("child1");

    root->push_back (child1);

    if (IOptionsParser* bloomParser = root->getParser ("bloom"))  {  bloomParser->setVisible(false);  }

    child1->push_back (new OptionOneParam (STR_NB_CORES,          "nb cores (0 for all)",   false, "0"));
    child1->push_back (new OptionOneParam (STR_MAX_MEMORY,        "max memory (in MBytes)", false, "2000"));
    child1->push_back (new OptionOneParam (STR_MAX_DISK,          "max disk   (in MBytes)", false, "0"));

    child1->push_back (BranchingAlgorithm<>::getOptionsParser());

    try
    {
        root->parse (argc, argv);
    }
    catch (OptionFailure& e)
    {
        e.displayErrors(std::cout);
    }

    printf ("-------------------------------------------------------------\n");
    RawDumpPropertiesVisitor v;
    root->getProperties()->accept (&v);

    delete root;
}
//! [snippet1]
