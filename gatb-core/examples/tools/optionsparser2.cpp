//! [snippet1]
// GATB-Core Online Tutorial 

/********************************************************************************/
/*                                                                              */
/*                Command line parsing through OptionsParser                    */
/*                                                                              */
/********************************************************************************/

// We include GATB-Core
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
// START Application
int main (int argc, char* argv[])
{
  // We create a new parser and give it a name
  OptionsParser parser ("CmdLineSnippet");

  // first argument of our program
  parser.push_back (
    // this is the class to use to setup an argument that requires a value
    new OptionOneParam (
      "-g",          // argument to use on the cmd-line
      "graph input", // help message
      true           // is it mandatory?
    ));

  //second argument
  parser.push_back (
    new OptionOneParam (
      "-v",   
      "verbosity (0:no display, 1: display kmers, 2: display distrib",  
      false, 
      "0" // default value
    ));

  //third argument; no value expected
  parser.push_back (
    // this is the class to use to setup an argument that 
    // DOES NOT require a value
    new OptionNoParam (
      "-h", "help", false
    ));

  // Now let's start the business of our application
  try
  {
    // First to do: analyse the software arguments using our Parser.
    // Here, if a mandatory argument is missing (-g in this example), then
    // an exception is raised and the software automatically dump a full help
    // message on the command-line... thanks to GATB-Core OptionsParser API.
    IProperties* options = parser.parse (argc, argv);

    // Here is how to retrieve the value for a particular option
    std::cout << "graph input file is: =" << options->getStr("-g") << std::endl;

    // Here is how to check if we have a particular no-value argument
    std::cout << "help requested? " << parser.saw("-h") << std::endl;
  }
  catch (OptionFailure& e)
  {
    e.displayErrors(std::cout);
  }
}
//! [snippet1]
