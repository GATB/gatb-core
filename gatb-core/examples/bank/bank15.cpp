//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                              Bank management                                 */
/*                                                                              */
/* This snippet shows how to open a FASTA bank and iterate its sequences        */
/* to provide some stats: data size, nb. sequences, etc.                        */
/*                                                                              */
/* Cmd-line: bank15 <fasta/q file>                                              */
/*                                                                              */
/* Sample: bank15 gatb-core/gatb-core/test/db/reads1.fa                         */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We get the file name from the user arguments
    const char* filename = argc >= 2 ? argv[1] : "";

    // We get information about the bank.
    u_int64_t nbSequences=0, dataSize=0, seqMaxSize=0, seqMinSize=~0;

    // We declare a Bank instance.
    IBank* bank = Bank::open (filename);
    LOCAL (bank);
    Iterator<Sequence>* it = bank->iterator();
    LOCAL (it);
    
    // We loop over sequences.
    for (it->first(); !it->isDone(); it->next())
    {
        Sequence& seq = it->item();

        Data& data = seq.getData();

        nbSequences ++;
        if (data.size() > seqMaxSize)  { seqMaxSize = data.size(); }
        if (data.size() < seqMinSize)  { seqMinSize = data.size(); }
        dataSize += data.size ();
    }

    std::cout << "# letters         : " << dataSize     << std::endl;
    std::cout << "# sequences       : " << nbSequences  << std::endl;
    std::cout << "sequence max size : " << seqMaxSize   << std::endl;
    std::cout << "sequence min size : " << seqMinSize   << std::endl;
}
//! [snippet1]
