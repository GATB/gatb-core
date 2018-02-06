//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                          Read two banks                                      */
/*                                                                              */
/* This snippet shows how to read two banks at the same time. Iterated items    */
/* are pairs of two sequences. This may be useful to read pair ends banks for   */
/* instance.                                                                    */
/*                                                                              */
/* This code produces on std::out the sequence description of each sequence     */
/* contained in the two sequence files.                                         */
/*                                                                              */
/* See also bank17.cpp where we do the same using a single album file.          */
/*                                                                              */
/* Cmd-line: bank18 -one <fasta/q file> -two <fasta/q file>                     */
/*                                                                              */
/* Sample: bank18 -one gatb-core/gatb-core/test/db/reads1.fa                    */
/*                -two gatb-core/gatb-core/test/db/reads2.fa                    */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    static const char* STR_BANK1 = "-one";
    static const char* STR_BANK2 = "-two";

    /** We create a command line parser. */
    OptionsParser parser ("TwoBanks");
    parser.push_back (new OptionOneParam (STR_BANK1, "bank one",   true));
    parser.push_back (new OptionOneParam (STR_BANK2, "bank two",   true));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        // We open the two banks
        IBank* bank1 = Bank::open (options->getStr(STR_BANK1));  LOCAL (bank1);
        IBank* bank2 = Bank::open (options->getStr(STR_BANK2));  LOCAL (bank2);

        // We iterate the two banks. Note how we provide two iterators from the two banks.
        PairedIterator<Sequence> *  itPair  = new  PairedIterator<Sequence> (bank1->iterator(), bank2->iterator());
		LOCAL(itPair);

		ProgressIterator< std::pair <Sequence, Sequence>> progress_iter (itPair, "coucou la paire", bank1->estimateNbItems());
		

        for (progress_iter.first(); !progress_iter.isDone(); progress_iter.next())
        {
            Sequence& s1 = itPair->item().first;
            Sequence& s2 = itPair->item().second;

            std::cout << "seq1.len=" << s1.getDataSize() << " seq2.len=" << s2.getDataSize() << std::endl;
        }
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
