#include <gatb/system/impl/System.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <iostream>

// We use the required packages
using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;
using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;
using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;


/********************************************************************************/

int main (int argc, char* argv[])
{
    if (argc < 2)
    {
        cerr << "you must provide at least one argument:" << endl;
        cerr << "   1) URI of input FASTA bank" << endl;
        cerr << "   2) URI of output binary bank (use input + suffix if not set)" << endl;
        cerr << "you must provide at least one argument:" << endl;
        return EXIT_FAILURE;
    }

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        // We get the URI of the FASTA bank
        string inputFilename (argv[1]);

        // We get the URI of the binary bank
        string outputFilename (argc >= 3 ? argv[2] : inputFilename + ".bin");

        // We declare the FASTA bank
        Bank bank (inputFilename);

        // We declare a binary bank
        BankBinary bankBin (outputFilename);

        // We declare some job listener.
        Progress progress (bank.estimateNbSequences(), "FASTA to binary conversion");

        // We convert the FASTA bank in binary format
        IProperties* props = BankHelper::singleton().convert (bank, bankBin, &progress);
        LOCAL (props);

        // We dump information about the conversion
        RawDumpPropertiesVisitor v;
        props->accept (&v);
    }

    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }

    return EXIT_SUCCESS;
}
