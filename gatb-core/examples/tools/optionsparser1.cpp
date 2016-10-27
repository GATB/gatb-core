//! [snippet1]

/********************************************************************************/
/*                                                                              */
/*      Command line parsing using OptionsParser and inheritance                */
/*      of arguments coming from other GATB-Core components.                    */          
/*                                                                              */
/********************************************************************************/

// We include GATB-Core
#include <gatb/gatb_core.hpp>

#include <list>
#include <iostream>

/********************************************************************************/
// START Application
int main (int argc, char* argv[])
{
  // We create a new parser ang give it a name
  IOptionsParser* root = new OptionsParser ("root");

  // inherit arguments from SortingCountAlgorithm and we request that none of
  // its argument are mandatory
  root->push_back (SortingCountAlgorithm<>::getOptionsParser(false));

  // Now, we setup our own parser
  IOptionsParser* child1 = new OptionsParser ("child1");
  
  // and we add that specific parser to the main parser that also
  // contains arguments from Graph
  root->push_back (child1);

  // we declare are own arguments
  child1->push_back (
    new OptionOneParam (
      STR_NB_CORES,
      "nb cores (0 for all)",
      false,
      "0"
    ));
  child1->push_back (
    new OptionOneParam (
      STR_MAX_MEMORY, 
      "max memory (in MBytes)", 
      false, 
      "2000"
    ));
  child1->push_back (
    new OptionOneParam (
      STR_MAX_DISK, 
      "max disk (in MBytes)",
      false,
      "0"
    ));

  try
  {
    root->parse (argc, argv);
  }
  catch (OptionFailure& e)
  {
    e.displayErrors(std::cout);
  }

  // We can dump the many arguments and their values
  printf ("-------------------------------------------------------------\n");
  RawDumpPropertiesVisitor v;
  root->getProperties()->accept (&v);

  // So, from this point we can use argument's values in the software.

  // .../...

  // when we do not need anymore the parser... delete it!
  delete root;
}
//! [snippet1]
