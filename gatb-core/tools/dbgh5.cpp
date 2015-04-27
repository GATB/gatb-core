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
static const char* STR_CHECK        = "-check";
static const char* STR_CHECK_DUMP   = "-check-dump";

int    manageException (IProperties* options, const std::string& message);
void   sendEmail       (IProperties* options, IProperties* graphInfo);
size_t checkResult     (const Graph& graph, IProperties* inputProps);

/********************************************************************************/
int main (int argc, char* argv[])
{
    size_t nbErrors = 0;

    /** We create a command line parser. */
    IOptionsParser* parser = Graph::getOptionsParser();  LOCAL (parser);

    IOptionsParser* parserGeneral = parser->getParser ("general");
    if (parserGeneral != 0)
    {
        /** We add an option to send the statistics by email. */
        parserGeneral->push_back (new OptionOneParam (STR_EMAIL,        "send statistics to the given email address", false));
        parserGeneral->push_back (new OptionOneParam (STR_EMAIL_FORMAT, "'raw' or 'xml'",                             false, "raw"));
        parserGeneral->push_back (new OptionOneParam (STR_CHECK,        "check result with previous result",          false));
        parserGeneral->push_back (new OptionOneParam (STR_CHECK_DUMP,   "dump some properties of the created graph into a file", false));
    }

    try
    {
        /** We first check that we have at least one argument. Otherwise, we print some information and leave. */
        if (argc == 1)
        {
            std::cerr << LibraryInfo::getInfo();
            OptionsHelpVisitor visitor (std::cerr);
            parser->accept (visitor, 0);
            return EXIT_SUCCESS;
        }

        /** We parse the user options. */
        IProperties* props = parser->parse (argc, argv);

        /** We create the graph with the provided options. */
        Graph graph = Graph::create (props);

        /** We may have to check the result. */
        if (props->get (STR_CHECK) > 0)  {  nbErrors = checkResult (graph, props);  }

        /** We dump some information about the graph. */
        if (props->getInt(STR_VERBOSE) > 0)  {  std::cout << graph.getInfo() << std::endl;  }

        /** We may have to send statistics by email. */
        if (props->get(STR_EMAIL))  {  sendEmail (props, &graph.getInfo());  }
    }
    catch (OptionFailure& e)
    {
        return e.displayErrors (std::cout);
    }
    catch (Exception& e)
    {
        return manageException (parser->getProperties(), e.getMessage());
    }
    catch (std::string& msg)
    {
        return manageException (parser->getProperties(), msg);
    }
    catch (const char* msg)
    {
        return manageException (parser->getProperties(), msg);
    }

    return nbErrors == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
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
int manageException (IProperties* options, const std::string& message)
{
    std::cout << std::endl << "EXCEPTION: " << message << std::endl;

    if (options && options->get(STR_EMAIL))
    {
        Properties props;  props.add(0, "error", message.c_str());
        sendEmail (options, &props);
    }

    return EXIT_FAILURE;
}

/********************************************************************************/
size_t checkResult (const Graph& graph, IProperties* inputProps)
{
    size_t nbErrors = 0;

    string filecheck = inputProps->getStr(STR_CHECK);

    /** We get the graph properties. */
    IProperties& graphProps = graph.getInfo();

    /** We read the check file. */
    Properties checkProps;  checkProps.readFile (filecheck);

    /** We need properties for output file if needed. */
    Properties outputProps;
    Properties tmpProps;

    /** We get the keys to be checked. */
    set<string> keys = checkProps.getKeys();
    for (set<string>::iterator it = keys.begin(); it != keys.end(); ++it)
    {
        if (graphProps.get(*it)==0)
        {
            tmpProps.add (0, "unknown", "%s", (*it).c_str());
            continue;
        }

        string v1 = checkProps.getStr(*it);
        string v2 = graphProps.getStr(*it);

        if (v1 != v2)
        {
            tmpProps.add (0, "diff", "%s", (*it).c_str());
            tmpProps.add (1, "val",  "%s", v1.c_str());
            tmpProps.add (1, "val",  "%s", v2.c_str());
            nbErrors++;
        }

        /** We put the graph property into the output props. */
        outputProps.add (0, *it, v2);
    }

    graphProps.add (1, "check", "%d/%d", keys.size()-nbErrors, keys.size());
    graphProps.add (2, tmpProps);

    if (inputProps->get(STR_CHECK_DUMP))
    {
        ofstream outputFile (inputProps->getStr(STR_CHECK_DUMP).c_str());
        if (outputFile.is_open())
        {
            RawDumpPropertiesVisitor visitor (outputFile, 40, ' ');
            outputProps.accept (&visitor);
        }
    }

    return nbErrors;
}
