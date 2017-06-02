//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/* !!!!!  WARNING !!!!!
 * DO NOT EDIT: snippet used to test Leon!
 */

/********************************************************************************/
/*                              Bank management                                 */
/*                                                                              */
/* Same as bank27 but provide results in a convenient way for script handling.  */
/*                                                                              */
/* Cmd-line: bank28 -in <fasta/q file>                                          */
/*                                                                              */
/* Sample: bank28 -in gatb-core/gatb-core/test/db/reads1.fa                     */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("BankStats");
    parser.push_back (new OptionOneParam (STR_URI_INPUT, "bank input",   true));
    parser.push_back (new OptionOneParam (STR_KMER_SIZE, "k-mer size",   false));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        // We get information about the bank.
        u_int64_t nbSequences=0, dataSize=0, seqMaxSize=0, seqMinSize=~0, nbSmallSequences=0;

        u_int64_t kmerSize = options->get(STR_KMER_SIZE) ? 
          options->getInt(STR_KMER_SIZE) : 31;
        // We declare an input Bank and use it locally
        IBank* inputBank = Bank::open (options->getStr(STR_URI_INPUT));
        LOCAL (inputBank);

        Iterator<Sequence>* it = inputBank->iterator();
        for (it->first(); !it->isDone(); it->next())
        {
            Data& data = it->item().getData();

            nbSequences ++;
            if (data.size() > seqMaxSize)  { seqMaxSize = data.size(); }
            if (data.size() < seqMinSize)  { seqMinSize = data.size(); }
            if (data.size() < kmerSize)    { nbSmallSequences++; }
            dataSize += data.size ();
        }
        printf("%u %u %u %u %u", dataSize, nbSequences, seqMaxSize, seqMinSize, nbSmallSequences);
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
