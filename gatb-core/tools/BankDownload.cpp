//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

// We use the required packages
using namespace std;

static const char* STR_BANK_URL = "-url";

/********************************************************************************/
/********************************************************************************/

/** This program loads banks from an GET, then iterates the bank to get some
 * statistics.
 *
 * The fine part is that these statistics are saved as metadata of the inode of the
 * bank file (supported by some file systems). Then, it is possible to get the
 * statistics with the shell 'getfattr' command.
 *
 * NOTE:
 *  1) if the downloaded file is gzipped, then the final bank is gunzipped
 *
 *  2) if the downloaded file is SRA format, then the final bank is converted in FASTA
 *     (provided the fact that the fastq-dump binary is available on the system)
 *
 * EXAMPLE:
 *      BankDownload -url ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/proteomes/HUMAN.fasta.gz
 *
 *      then, you can type "getfattr -d HUMAN.fasta" and get the following information:
 *
 *          # file: HUMAN.fasta
 *          user.bank_name="HUMAN"
 *          user.bank_type="prot"
 *          user.data_size="35208089"
 *          user.download_date="20140612_175614"
 *          user.download_url="ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/proteomes/HUMAN.fasta.gz"
 *          user.file_size="44489858"
 *          user.file_url="/local/databases/ref/proteomes/HUMAN.fasta"
 *          user.seq_max_size="35991"
 *          user.seq_min_size="11"
 *          user.seq_number="88949"
 */

int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("BankDownload");
    parser.push_back (new OptionOneParam (STR_BANK_URL,       "bank url",   true));
    parser.push_back (new OptionOneParam (STR_URI_OUTPUT_DIR, "output dir", false));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        /** Shortcuts */
        string url     = options->getStr(STR_BANK_URL);
        string dateStr = System::time().getDateString();

        /** We get the bank name from the url. */
        string bankFile = url;
        const size_t last_slash_idx = bankFile.find_last_of("\\/");
        if (std::string::npos != last_slash_idx)  {  bankFile.erase (0, last_slash_idx + 1);  }

        /** We erase the file if already exists. */
        System::file().remove (bankFile);

        /** We download the bank with wget. */
        stringstream ss1;     ss1 << "wget '" << url << "'";
        int res1 = system (ss1.str().c_str());
        if (res1!=EXIT_SUCCESS)  { throw Exception ("Error with wget command"); }

        /** We test whether we got the file or not. */
        if (System::file().doesExist(bankFile) == false)  {  throw Exception("Unable to download '%s'", url.c_str());  }

        /** We test whether the file is gzipped, SRA, ... */
        bool isGzipped = false;
        bool isSRA     = false;
        const size_t period_idx = bankFile.rfind('.');
        string finalBankFile = bankFile;
        if (std::string::npos != period_idx)
        {
            if (bankFile.substr(period_idx).compare(".gz")==0)
            {
                isGzipped = true;
                finalBankFile.erase (period_idx);
            }

            if (bankFile.substr(period_idx).compare(".sra")==0)
            {
                isSRA = true;
                finalBankFile.erase (period_idx);
                finalBankFile += ".fasta";
            }
        }

        if (isGzipped==true)
        {
            /** We erase the file if already exists. */
            System::file().remove (finalBankFile);

            /** We unzip the file. */
            stringstream ss2;     ss2 << "gunzip " << bankFile;
            int res2 = system (ss2.str().c_str());
            if (res2!=EXIT_SUCCESS)  { throw Exception ("Error with gunzip command"); }
        }

        if (isSRA==true)
        {
            /** We convert the SRA file to FASTA format. */
            stringstream ss3;     ss3 << "fastq-dump --fasta 60 " << bankFile;
            int res3 = system (ss3.str().c_str());
            if (res3!=EXIT_SUCCESS)  { throw Exception ("Error with fastq-dump command"); }
        }

        // We open the bank
        BankFasta bank (finalBankFile);

        /** We get information about the bank. */
        u_int64_t nbSequences=0, dataSize=0, seqMaxSize=0, seqMinSize=~0;
        u_int64_t nbLetters=0, nbNucleotides=0;

        ProgressIterator<Sequence,Progress> it (bank, "iterate", 1000);
        for (it.first(); !it.isDone(); it.next())
        {
            // Shortcuts
            Sequence& seq = it.item();
            Data&     data = seq.getData();

            nbSequences ++;
            if (data.size() > seqMaxSize)  { seqMaxSize = data.size(); }
            if (data.size() < seqMinSize)  { seqMinSize = data.size(); }
            dataSize += data.size ();

            for (size_t i=0; i<seq.getDataSize(); i++, nbLetters++)
            {
                char c = tolower(seq.getDataBuffer()[i]);
                if (c=='a' || c=='c' || c=='g' || c=='t') { nbNucleotides ++; }
            }
        }

        string dataKind = "??";
        if (nbLetters > 0)  {  dataKind = (float)nbNucleotides / (float)nbLetters > 0.9 ? "nucl" : "prot"; }

        string bankName = System::file().getBaseName(finalBankFile);

        /** We get the file size. */
        u_int64_t fileSize = System::file().getSize(finalBankFile);

        /** We move the file if needed. */
        string finalBankUrl = finalBankFile;
        if (options->get(STR_URI_OUTPUT_DIR) != 0)
        {
            finalBankUrl = options->getStr(STR_URI_OUTPUT_DIR) + "/" + finalBankUrl;
            system (Stringify::format ("mv %s %s/", finalBankFile.c_str(), options->getStr(STR_URI_OUTPUT_DIR).c_str()).c_str());
        }

        /** We tag the file with information. */
        System::file().setAttribute (finalBankUrl, "bank_name",     "%s",  bankName.c_str());
        System::file().setAttribute (finalBankUrl, "bank_type",     "%s",  dataKind.c_str());
        System::file().setAttribute (finalBankUrl, "data_size",     "%ld", dataSize);
        System::file().setAttribute (finalBankUrl, "download_date", "%s",  dateStr.c_str());
        System::file().setAttribute (finalBankUrl, "download_url",  "%s",  url.c_str());
        System::file().setAttribute (finalBankUrl, "seq_number",    "%ld", nbSequences);
        System::file().setAttribute (finalBankUrl, "seq_max_size",  "%ld", seqMaxSize);
        System::file().setAttribute (finalBankUrl, "seq_min_size",  "%ld", seqMinSize);
        System::file().setAttribute (finalBankUrl, "file_url",      "%s",  finalBankUrl.c_str());
        System::file().setAttribute (finalBankUrl, "file_size",     "%ld", fileSize);

        /** Cleanup. */
        if (isSRA == true)
        {
            System::file().remove (bankFile);
        }
    }

    catch (OptionFailure& e)
    {
        e.getParser().displayErrors (stdout);
        e.getParser().displayHelp   (stdout);
        return EXIT_FAILURE;
    }
    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }
}
//! [snippet1]
