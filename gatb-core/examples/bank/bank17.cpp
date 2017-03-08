//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                              Bank album                                      */
/*                                                                              */
/* A BankAlbum file is a bank defined by a list of URI of other banks.          */
/*                                                                              */
/* This code produces on std::out the total bank size, then each sequence       */
/* description.                                                                 */
/*                                                                              */
/* Cmd-line: bank17 -in <album file>                                            */
/*                                                                              */
/* Sample: bank17 -in gatb-core/gatb-core/test/db/album.txt                     */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("BankStats");
    parser.push_back (new OptionOneParam (STR_URI_INPUT, "bank input",   true));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        std::string filename = options->getStr(STR_URI_INPUT);

        //! [snippet17_album]
        // We declare a Bank instance.
        BankAlbum bank (filename);

        // We dump some information about the bank
        std::cout << "cummulated files sizes : " << bank.getSize() << std::endl;

        // We create an iterator on the bank
        Iterator<Sequence>* it = bank.iterator();
        LOCAL (it);

        // We iterate the sequences of the bank
        for (it->first(); !it->isDone(); it->next())
        {
            // We dump some information about the sequence.
            std::cout << "comment " << it->item().getComment() << std::endl;
        }
        //! [snippet17_album]
    }
    catch (OptionFailure& e)
    {
        return e.displayErrors (std::cout);
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }
}
//! [snippet1]
