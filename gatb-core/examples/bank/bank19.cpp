//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                              Bank checksum                                   */
/*                                                                              */
/* This snippet shows a checkum value on the genomic data.                      */
/*                                                                              */
/* Cmd-line: bank19 -in <fasta/q file>                                          */
/*                                                                              */
/* Sample: bank19 -in gatb-core/gatb-core/test/db/reads1.fa                     */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("BankCheckum");
    parser.push_back (new OptionOneParam (STR_URI_INPUT, "bank input",   true));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        // We get information about the bank.
        u_int64_t checksum1=0;
        u_int64_t checksum2=0;

        // We declare an input Bank and use it locally
        IBank* inputBank = Bank::open (options->getStr(STR_URI_INPUT));
        LOCAL (inputBank);

        ProgressIterator<Sequence,ProgressNone> it (*inputBank, "iterate");
        for (it.first(); !it.isDone(); it.next())
        {
            Data& data = it.item().getData();

            for (size_t i=0; i<data.size(); i++)
            {
                checksum1 += data[i];
                checksum2 += data[i] * it.item().getIndex();
            }
        }

        std::cout << "checksum1=" << std::hex << checksum1 << "  checksum2=" << checksum2 << std::endl;
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
