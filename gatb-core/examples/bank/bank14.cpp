//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                              Bank management                                 */
/*                                                                              */
/* This snippet shows how to open a FASTA bank and iterate its sequences.       */
/* Some attributes of the iterated Sequence objects are used.                   */
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

        // We get information about the bank.
        u_int64_t nbSequences=0, dataSize=0, seqMaxSize=0, seqMinSize=~0;

        // We declare an input Bank and use it locally
        IBank* inputBank = Bank::open (options->getStr(STR_URI_INPUT));
        LOCAL (inputBank);

        ProgressIterator<Sequence> it (*inputBank, "iterate");
        for (it.first(); !it.isDone(); it.next())
        {
            Data& data = it.item().getData();

            nbSequences ++;
            if (data.size() > seqMaxSize)  { seqMaxSize = data.size(); }
            if (data.size() < seqMinSize)  { seqMinSize = data.size(); }
            dataSize += data.size ();
        }

        std::cout << "data size         : " << dataSize     << std::endl;
        std::cout << "sequence number   : " << nbSequences  << std::endl;
        std::cout << "sequence max size : " << seqMaxSize   << std::endl;
        std::cout << "sequence min size : " << seqMinSize   << std::endl;
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
