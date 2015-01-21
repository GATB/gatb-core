//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>
#include <fstream>
#include <string>

// We use the required packages
using namespace std;

/********************************************************************************/
/*                                Bank copy                                     */
/*                                                                              */
/* This snippet shows how to copy a bank into another one, with some possible   */
/* modifications. Note that this snippet can be really useful for re-formatting */
/* FASTA banks                                                                  */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 2)
    {
        cerr << "you must provide at least the FASTA file path. Arguments are:" << endl;
        cerr << "   1) FASTA file path" << endl;
        cerr << "   2) data lines size (60 by default, 0 means only one line for sequence data)" << endl;
        cerr << "   3) type of comments: 0 for new comments, 1 for only ids, 2 (default) for full comments" << endl;
        cerr << "   4) number of sequences to keep (all by default)" << endl;
        return EXIT_FAILURE;
    }

    // We get the file name from the user arguments
    const char* filename = argv[1];

    // We define the max size of a data line in the FASTA output file
    size_t dataLineSize = argc >= 3 ?  atoi(argv[2]) : 60;

    // By convention, setting data line size to 0 will dump data on a single line
    if (dataLineSize==0)  { dataLineSize = ~0; }

    // We define the kind of output comments
    BankFasta::Iterator::CommentMode_e mode = BankFasta::Iterator::FULL;
    if (argc >= 4)
    {
        int m = atoi (argv[3]);
             if (m==0)  { mode = BankFasta::Iterator::NONE;   }
        else if (m==1)  { mode = BankFasta::Iterator::IDONLY; }
    }

    // We set the number of sequences to be kept.
    u_int64_t nbSequences =  (argc >= 5 ? atol (argv[4]) : ~0);

    try  {
        // We declare a Bank instance.
        BankFasta b (filename);

        // We create an iterator over this bank.
        // Note : here, we must use specifically BankFasta::Iterator in order to set the mode
        BankFasta::Iterator itSeq (b, mode);

        // We encapsulate it with a truncation iterator
        TruncateIterator<Sequence> it (itSeq, nbSequences);

        size_t idxSeq = 1;

        // We loop over sequences.
        for (it.first(); !it.isDone(); it.next(), idxSeq++)
        {
            // We dump the comment into the file according to the user mode
            if (!it->getComment().empty())   {  cout << ">" << it->getComment() << endl;  }
            else                             {  cout << ">seq=" << idxSeq << " len=" << it.item().getDataSize() << endl;  }

            // shortcut
            size_t len = it->getDataSize();

            // We dump the data with fixed sized columns
            for (size_t i=0; i<len; )
            {
                for (size_t j=0; j<dataLineSize && i<len; j++, i++)
                {
                    cout << (char)  (it->getData() [i]);
                }
                cout << endl;
            }
        }
    }
    catch (Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
