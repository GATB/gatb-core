//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

// We use the required packages
using namespace std;

static const char* STR_FILTER_RATIO = "-filter-ratio";

/********************************************************************************/
/*                  Filter sequences having too many bad letters.               */
/*                                                                              */
/* This snippet shows how to keep sequences that have enough "good" data.       */
/*                                                                              */
/* A sequence is retained only if nb(A,C,G,T)/seq_size>=filter-ratio .          */
/*                                                                              */
/* Cmd-line: bank13 -in <fasta/q file> -filter-ratio <[0.0..1.0]>               */
/*                                                                              */
/* Sample: bank13 -in gatb-core/gatb-core/test/db/reads1.fa -filter-ratio 0.8   */
/*         (output bank will be file:  reads1.fa_filtered)                      */
/*                                                                              */
/* WARNING ! THIS SNIPPET SHOWS ALSO HOW TO USE LAMBDA EXPRESSIONS, SO YOU NEED */
/* TO USE A COMPILER THAT SUPPORTS THIS FEATURE.                                */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("BankFilter");
    parser.push_back (new OptionOneParam (STR_URI_INPUT,     "bank input",   true));
    parser.push_back (new OptionOneParam (STR_FILTER_RATIO,  "skip a sequence if 'good letters number / seq.len > X'",   false, "0.8"));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        /** Shortcuts. */
        double percentThreshold = options->getDouble(STR_FILTER_RATIO);

        /** We open the input bank. */
        IBank* inBank = Bank::open (options->getStr(STR_URI_INPUT));
        LOCAL (inBank);

        /** We create the output inBank. */
        IBank* outBank = new BankFasta (options->getStr(STR_URI_INPUT) + "_filtered");
        LOCAL (outBank);

        /** We iterate the inBank. NOTE: WE USE A LAMBDA EXPRESSION HERE. */
        inBank->iterate ([&] (Sequence& s)
        {
            /** Shortcut. */
            char* data = s.getDataBuffer();

            size_t nbOK = 0;
            for (size_t i=0; i<s.getDataSize(); i++)
            {
                if (data[i]=='A' || data[i]=='C' || data[i]=='G' || data[i]=='T')  { nbOK++; }
            }

            if ((double)nbOK / (double)s.getDataSize() > percentThreshold)  {  outBank->insert (s);  }
        });

        /** We flush the output bank. */
        outBank->flush();
    }

    catch (OptionFailure& e)
    {
        return e.displayErrors (cout);
    }
    catch (Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }
}
//! [snippet1]
