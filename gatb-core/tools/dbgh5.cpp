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
#include <sstream>

/********************************************************************************/
static const char* STR_EMAIL        = "-email";
static const char* STR_EMAIL_FORMAT = "-email-fmt";

int  manageException (IOptionsParser* parser, IProperties* options, const std::string& message);
void sendEmail       (IProperties* options, IProperties* graphInfo);

/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    IOptionsParser* parser = Graph::getOptionsParser();  LOCAL (parser);

    OptionsParser* parserGeneral = dynamic_cast<OptionsParser*> (parser->getParser ("general"));
    if (parserGeneral != 0)
    {
        /** We add an option to send the statistics by email. */
        parserGeneral->push_back (new OptionOneParam (STR_EMAIL,        "send statistics to the given email address", false));
        parserGeneral->push_back (new OptionOneParam (STR_EMAIL_FORMAT, "'raw' or 'xml'",                             false, "raw"));
    }

    try
    {
        /** We parse the user options. */
        IProperties* options = parser->parse (argc, argv);

        if (options->get(STR_HELP))     {  parser->displayHelp    (std::cout);  return EXIT_SUCCESS; }
        if (options->get(STR_VERSION))  {  parser->displayVersion (std::cout);  return EXIT_SUCCESS; }

        /** We create the graph with the provided options. */
        Graph graph = Graph::create (options);

        /** We dump some information about the graph. */
        if (options->getInt(STR_VERBOSE) > 0)  {  std::cout << graph.getInfo() << std::endl;  }

        /** We may have to send statistics by email. */
        if (options->get(STR_EMAIL))  {  sendEmail (options, &graph.getInfo());  }
    }
    catch (OptionFailure& e)
    {
        e.getParser().displayErrors   (std::cout);
        e.getParser().displayWarnings (std::cout);
        e.getParser().displayHelp     (std::cout);
        return EXIT_FAILURE;
    }
    catch (Exception& e)
    {
        return manageException (parser, parser->getProperties(), e.getMessage());
    }
    catch (std::string& msg)
    {
        return manageException (parser, parser->getProperties(), msg);
    }
    catch (const char* msg)
    {
        return manageException (parser, parser->getProperties(), msg);
    }

    return EXIT_SUCCESS;
}

/********************************************************************************/
void sendEmail (IProperties* options, IProperties* graphInfo)
{
    std::string outfmt = options->getStr(STR_EMAIL_FORMAT);

    /** We build the mail content. */
    std::stringstream output;
    if (outfmt == "raw")
    {
        RawDumpPropertiesVisitor visit (output);
        graphInfo->accept (&visit);
    }
    else if (outfmt == "xml")
    {
        XmlDumpPropertiesVisitor visit (output);
        graphInfo->accept (&visit);
    }
    else
    {
        /** Unknown format. */
        throw Exception ("Unable to send email because of unknown format '%s'", outfmt.c_str());
    }

    /** We build the mail command. */
    std::stringstream cmd;
    cmd << "echo \"" << output.str()
        << "\" | mail -s \"[dbgh5] " << System::file().getBaseName(options->getStr(STR_URI_INPUT)) << "\" "
        << options->getStr(STR_EMAIL);

    /** We execute the command. */
   ::system (cmd.str().c_str());
}

/********************************************************************************/
int manageException (IOptionsParser* parser, IProperties* options, const std::string& message)
{
    std::cout << std::endl << "EXCEPTION: " << message << std::endl;

    if (parser != 0)
    {
        parser->displayErrors   (std::cout);
        parser->displayWarnings (std::cout);
        parser->displayHelp     (std::cout);
    }

    if (options && options->get(STR_EMAIL))
    {
        Properties props;  props.add(0, "error", message.c_str());
        sendEmail (options, &props);
    }

    return EXIT_FAILURE;
}

