//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>
#include <sstream>

// We use the required packages
using namespace std;

static const char* STR_BANKS_NB       = "-split";
static const char* STR_MAX_INPUT_SIZE = "-max-size";
static const char* STR_OUTPUT_FASTQ   = "-fastq";
static const char* STR_OUTPUT_GZ      = "-gz";

/********************************************************************************/
/*                         Bank split                                           */
/*                                                                              */
/* This snippet shows how to split a bank into smaller banks and how to create  */
/* an album bank (ie a list of URL of banks). Such an album bank could be used  */
/* as bank input by other tools.                                                */
/* Note: all the generated files are put in a directory created by the snippet. */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("BankSplitter");
    parser.push_back (new OptionOneParam (STR_URI_INPUT,      "bank reference",            true));
    parser.push_back (new OptionOneParam (STR_MAX_INPUT_SIZE, "average db size per split", true));
    parser.push_back (new OptionOneParam (STR_URI_OUTPUT_DIR, "output directory",          false, "."));
    parser.push_back (new OptionNoParam  (STR_OUTPUT_FASTQ,   "fastq output",              false));
    parser.push_back (new OptionNoParam  (STR_OUTPUT_GZ,      "gzip output",               false));

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        /** Shortcuts. */
        u_int64_t maxDbSize = options->getInt(STR_MAX_INPUT_SIZE);

        // We declare an input Bank
        IBank* inputBank = Bank::open (options->getStr(STR_URI_INPUT));
        LOCAL (inputBank);

        // We get the basename of the input bank.
        string inputBasename = System::file().getBaseName (options->getStr(STR_URI_INPUT));

        /** We set the name of the output directory. */
        stringstream ss;  ss << inputBasename << "_S" << maxDbSize;
        string outputDirName = ss.str();

        /** We create the output directory. */
        string outputDir = options->getStr(STR_URI_OUTPUT_DIR) + "/" + outputDirName;
        System::file().mkdir (outputDir, S_IRWXU);

        // We create the album bank.
        BankAlbum album (outputDir + "/album.txt");

        /** We get estimations about the bank. */
        u_int64_t number, totalSize, maxSize;
        inputBank->estimate (number, totalSize, maxSize);

        u_int64_t estimationNbSeqToIterate = number;

        // We create an iterator over the input bank
        ProgressIterator<Sequence> itSeq (*inputBank, "split");

        // We loop over sequences to get the exact number of sequences.
          int64_t nbBanksOutput = -1;
        u_int64_t nbSequences   =  0;
        u_int64_t dbSize        = ~0;

        bool isFastq   = options->get(STR_OUTPUT_FASTQ) != 0;
        bool isGzipped = options->get(STR_OUTPUT_GZ)    != 0;

        IBank* currentBank = 0;

        for (itSeq.first(); !itSeq.isDone(); itSeq.next())
        {
            if (dbSize > maxDbSize)
            {
                if (currentBank != 0)  { currentBank->flush();  currentBank->finalize(); }

                nbBanksOutput ++;

                /** We build the uri of the current bank. */
                stringstream ss;  ss << inputBasename << "_" << nbBanksOutput << (isFastq ? ".fastq" : ".fasta");
                if (isGzipped) { ss << ".gz"; }

                /** We create a new bank and put it in the album. */
                currentBank = album.addBank (outputDir, ss.str(), isFastq, isGzipped);

                /** We reinit the db size counter. */
                dbSize = 0;
            }

            dbSize += itSeq->getDataSize();

            /** We insert the sequence into the current output bank. */
            currentBank->insert (*itSeq);
        }

        if (currentBank != 0)  { currentBank->flush(); }
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
