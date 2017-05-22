//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                              Bank management                                 */
/*                                                                              */
/* This snippet shows how to open a FastQ file and iterate its sequences        */
/* using a filter discarding all reads having a Phred score >= 30.              */
/*                                                                              */
/* Cmd-line: bank24 <fastq file>                                                */
/*                                                                              */
/* Sample: bank23 gatb-core/gatb-core/test/db/sample.fastq                      */
/*                                                                              */
/********************************************************************************/
int threshold = 30;

int computeMeanPhredScore(const std::string& quality){
    int score=0;
    // quality information is supposed to be FastQ-Sanger-encoded: 
    // each letter in quality string is an ASCII code ranging from 33 to 126.
    for(char c : quality){
        score += (c-33);
    }
    return score / quality.size();
}

struct QualityFilter { 
    bool operator () (Sequence& seq) const  { 
        return computeMeanPhredScore(seq.getQuality()) >= threshold;
    } 
};

/********************************************************************************/
// START Application
int main (int argc, char* argv[])
{
  // We check that the user provides at least one option: a Fasta/FastQ file.
  // Online GATB-Tutorial: this argument is automatically filled in with an 
  // appropriate file.
  if (argc < 2)
  {
    std::cerr << "Please, provide a sequence file." << std::endl;
    return EXIT_FAILURE;
  }

  // We define a try/catch block in case some method fails (bad filename for instance)
  try
  {
    // We declare an input Bank and use it locally
    IBank* inputBank = Bank::open (argv[1]);
    LOCAL (inputBank);

    // We create an iterator over this bank using some filtering system
    FilterIterator<Sequence,QualityFilter> it (inputBank->iterator(), QualityFilter());
    
    // We loop over sequences.
    for (it.first(); !it.isDone(); it.next())
    {
      // Shortcut
      Sequence& seq = it.item();

      // We dump the sequence quality
      std::cout << "[" << seq.getQuality() << "] " << computeMeanPhredScore(seq.getQuality()) << std::endl;

    }
  }
  catch (Exception& e)
  {
    std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
  }
}

//! [snippet1]

