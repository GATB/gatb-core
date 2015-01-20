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
/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 3)  {  std::cout << "You must provide two banks" << std::endl;    return EXIT_FAILURE; }

    // We get the file names from the user arguments
    const char* file1 = argv[1];
    const char* file2 = argv[2];

    // We open the two banks
    IBank* bank1 = Bank::open (file1);  LOCAL (bank1);
    IBank* bank2 = Bank::open (file2);  LOCAL (bank2);

    // We get one iterator for each bank
    Iterator<Sequence>* it1 = bank1->iterator();  LOCAL (it1);
    Iterator<Sequence>* it2 = bank2->iterator();  LOCAL (it2);

    // We iterate the two banks
    PairedIterator<Sequence> itPair (it1, it2);
    for (itPair.first(); !itPair.isDone(); itPair.next())
    {
        Sequence& s1 = itPair.item().first;
        Sequence& s2 = itPair.item().second;

        std::cout << "seq1.len=" << s1.getDataSize() << " seq2.len=" << s2.getDataSize() << std::endl;
    }
}
//! [snippet1]
