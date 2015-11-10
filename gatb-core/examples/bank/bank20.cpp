#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
void dump (IBank* bank, size_t depth=0)
{
    LOCAL (bank);

    for (size_t i=0; i<depth; i++)  { cout << "  "; }  cout << bank->getId() << endl;

    BankComposite* composite = dynamic_cast<BankComposite*> (bank);
    if (composite != 0)
    {
        const std::vector<IBank*> children = composite->getBanks();
        for (size_t i=0; i<children.size(); i++)  {  dump (children[i], depth+1); }
    }
}

/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("BankDump");
    parser.push_back (new OptionOneParam (STR_URI_INPUT, "bank input",   true));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        /** We dump the bank hierarchy. */
        dump (Bank::open (options->getStr(STR_URI_INPUT)));
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
