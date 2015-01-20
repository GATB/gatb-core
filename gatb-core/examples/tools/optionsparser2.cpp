//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                Command line parsing through IOptionParser                    */
/********************************************************************************/
int main (int argc, char* argv[])
{
    IOptionsParser* root = new OptionsParser ("root");
    LOCAL (root);

    root->push_back (new OptionNoParam (STR_HELP, "help",   false));

    try
    {
        root->parse (argc, argv);
    }
    catch (OptionFailure& e)
    {
        e.displayErrors(std::cout);
    }

    std::cout << "status=" << root->saw(STR_HELP) << std::endl;
}
//! [snippet1]
