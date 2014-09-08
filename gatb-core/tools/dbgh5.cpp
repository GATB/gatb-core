/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#include <gatb/gatb_core.hpp>

/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser = Graph::getOptionsParser();
    parser.setName ("dbgh5");

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

         /** We create the graph with the provided options. */
         Graph graph = Graph::create (options);

         /** We dump some information about the graph. */
         if (options->getInt(STR_VERBOSE) > 0)  {  std::cout << graph.getInfo() << std::endl;  }
    }
    catch (OptionFailure& e)
    {
        e.getParser().displayErrors   (stdout);
        e.getParser().displayWarnings (stdout);
        e.getParser().displayHelp     (stdout);
        e.getParser().displayVersion  (stdout);
        return EXIT_FAILURE;
    }
    catch (Exception& e)
    {
        std::cout << std::endl << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
    catch (std::string& msg)
    {
        std::cout << std::endl << "EXCEPTION: " << msg << std::endl;
        return EXIT_FAILURE;
    }
    catch (const char* msg)
    {
        std::cout << std::endl << "EXCEPTION: " << msg << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
