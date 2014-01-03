/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/gatb_core.hpp>

#include <sstream>

/********************************************************************************/
void SaveAsFasta (IBank* bank, const std::string& uri)
{
    BankFasta output (uri);

    Iterator<Sequence>* it = bank->iterator();
    LOCAL (it);

    std::stringstream ss;

    size_t count=0;
    for (it->first(); !it->isDone(); it->next())
    {
        // Shortcut.
        Sequence& seq = it->item();

        // We define some decent header
        ss.str ("");
        ss << count++ << "__len__" << seq.getDataSize();
        seq.setComment (ss.str());

        output.insert (seq);
    }
}

/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("bankgen");

    const char* OUTPUT_PREFIX = "-out";
    const char* SEQ_LEN       = "-seq-len";
    const char* READ_LEN      = "-read-len";
    const char* OVERLAP_LEN   = "-overlap-len";
    const char* COVERAGE      = "-coverage";

    parser.push_back (new OptionOneParam (OUTPUT_PREFIX,  "output prefix",               true));
    parser.push_back (new OptionOneParam (SEQ_LEN,        "sequence length",             false,  "1000"));
    parser.push_back (new OptionOneParam (READ_LEN,       "read length",                 false,  "100" ));
    parser.push_back (new OptionOneParam (OVERLAP_LEN,    "overlap between two reads",   false,  "50" ));
    parser.push_back (new OptionOneParam (COVERAGE,       "coverage",                    false,  "3" ));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        /** We create the random sequence. */
        IBank* randomBank = new BankRandom (1, options->getInt(SEQ_LEN));
        LOCAL (randomBank);

        /** We create the reads bank. */
        IBank* readsBank = new BankSplitter (
            randomBank,
            options->getInt(READ_LEN),
            options->getInt(OVERLAP_LEN),
            options->getInt(COVERAGE)
        );
        LOCAL (readsBank);

        /** We save the random bank. */
        SaveAsFasta (randomBank, options->getStr(OUTPUT_PREFIX) + "_random.fa");

        /** We save the reads bank. */
        SaveAsFasta (readsBank, options->getStr(OUTPUT_PREFIX) + "_reads.fa");
    }
    catch (OptionFailure& e)
    {
        e.getParser().displayErrors (stdout);
        e.getParser().displayHelp   (stdout);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
