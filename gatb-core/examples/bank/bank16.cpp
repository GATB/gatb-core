//! [snippet16]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                              Bank opening                                    */
/*                                                                              */
/* This snippet shows how to open a bank with a URI and without specifying      */
/* the actual type of the bank. The correct bank format is found by analyzing   */
/* the URI and/or the content of the resource given by the URI                  */
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

        //! [snippet16_bank]
        // We get an instance of IBank from the URI.
        IBank* bank = Bank::open (filename);

        //! [snippet16_seq]
        // We create an iterator on the bank
        Iterator<Sequence>* it = bank->iterator();

        // We iterate the sequences of the bank
        for (it->first(); !it->isDone(); it->next())
        {
            // We get a shortcut on the current sequence and its data
            Sequence& seq  = it->item();
            Data&     data = seq.getData();

            // We dump some information about the sequence.
            std::cout << "comment " << seq.getComment() << std::endl;

            // We dump each nucleotide. NOTE: the output depends on the data encoding
            for (size_t i=0; i<data.size(); i++)  {  std::cout << data[i];  }  std::cout << std::endl;
        }

        //! [snippet16_seq]
        // The bank and the iterator have been allocated on the heap, so we have to delete them
        delete it;
        delete bank;
        //! [snippet16_bank]
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
//! [snippet16]
