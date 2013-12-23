/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file StringsRepository.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Pool of strings providing information to end users
 */

#ifndef _GATB_CORE_TOOLS_MISC_STRINGS_REPOSITORY_HPP_
#define _GATB_CORE_TOOLS_MISC_STRINGS_REPOSITORY_HPP_

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Tools package */
namespace tools     {
/** \brief Misc interfaces */
namespace misc      {
/********************************************************************************/

/** \brief Pool of strings
 *
 * This class provides constant strings used throughout the code. It may be interesting
 * in order to have a central point for constant strings management:
 *      - ease translation in different languages
 *      - entry point for strings obfuscation if needed
 *
 * It could also be possible to read the strings from a configuration file.
 *
 * The class defines one (static) method per constant string to be used. Note that we
 * also (see below) define a macro definition that eases the use of such a facility.
 */

class MessageRepository
{
public:

    /** \brief Singleton method.
     *
     * This method could return different types for the string repository, for translation
     * for instance.
     *
     * \return the singleton instance.
     */
    static MessageRepository& singleton()  { static MessageRepository instance; return instance; }

    const char* BANK_bad_file_number    () { return "bank files number is %d but should be in [1..%d]"; }
    const char* BANK_bad_file_path      () { return "unable to find file '%s'"; }
    const char* BANK_unable_open_file   () { return "error opening file: %s"; }
    const char* BANK_unable_write_file  () { return "unable to write into file"; }
};

/********************************************************************************/

/** \brief Pool of strings
 */
class StringRepository
{
public:
    static StringRepository& singleton()  { static StringRepository instance; return instance; }

    const char* db             ()  { return "-db";             }
    const char* kmer_size      ()  { return "-kmer-size";      }
    const char* nks            ()  { return "-nks";            }
    const char* max_memory     ()  { return "-max-memory";     }
    const char* max_disk       ()  { return "-max-disk";       }
    const char* kmer_solid     ()  { return "-kmer-solid";     }
    const char* kmer_cFP       ()  { return "-kmer-cFP";       }
    const char* prefix         ()  { return "-prefix";         }
    const char* progress_bar   ()  { return "-bargraph";       }
    const char* nb_cores       ()  { return "-nb-cores";       }
    const char* partition_type ()  { return "-partition-type"; }
    const char* uri_histogram  ()  { return "-histo";          }
    const char* uri_debloom    ()  { return "-debloom";        }
    const char* uri_input      ()  { return "-in";             }
    const char* uri_output     ()  { return "-out";            }
    const char* uri_output_dir ()  { return "-out-dir";        }
    const char* verbose        ()  { return "-verbose";        }
    const char* help           ()  { return "-help";           }
};

/********************************************************************************/

/** Shortcuts. */
#define STR_URI_DB              gatb::core::tools::misc::StringRepository::singleton().db ()
#define STR_KMER_SIZE           gatb::core::tools::misc::StringRepository::singleton().kmer_size ()
#define STR_NKS                 gatb::core::tools::misc::StringRepository::singleton().nks ()
#define STR_MAX_MEMORY          gatb::core::tools::misc::StringRepository::singleton().max_memory ()
#define STR_MAX_DISK            gatb::core::tools::misc::StringRepository::singleton().max_disk ()
#define STR_KMER_SOLID          gatb::core::tools::misc::StringRepository::singleton().kmer_solid ()
#define STR_KMER_CFP            gatb::core::tools::misc::StringRepository::singleton().kmer_cFP ()
#define STR_PREFIX              gatb::core::tools::misc::StringRepository::singleton().prefix ()
#define STR_PROGRESS_BAR        gatb::core::tools::misc::StringRepository::singleton().progress_bar ()
#define STR_NB_CORES            gatb::core::tools::misc::StringRepository::singleton().nb_cores ()
#define STR_PARTITION_TYPE      gatb::core::tools::misc::StringRepository::singleton().partition_type ()
#define STR_URI_HISTOGRAM       gatb::core::tools::misc::StringRepository::singleton().uri_histogram ()
#define STR_URI_DEBLOOM         gatb::core::tools::misc::StringRepository::singleton().uri_debloom ()
#define STR_URI_INPUT           gatb::core::tools::misc::StringRepository::singleton().uri_input ()
#define STR_URI_OUTPUT          gatb::core::tools::misc::StringRepository::singleton().uri_output ()
#define STR_URI_OUTPUT_DIR      gatb::core::tools::misc::StringRepository::singleton().uri_output_dir ()
#define STR_VERBOSE             gatb::core::tools::misc::StringRepository::singleton().verbose ()
#define STR_HELP                gatb::core::tools::misc::StringRepository::singleton().help ()

/********************************************************************************/

/** Shortcuts. */
#define STR_BANK_bad_file_number    gatb::core::tools::misc::MessageRepository::singleton().BANK_bad_file_number ()
#define STR_BANK_bad_file_path      gatb::core::tools::misc::MessageRepository::singleton().BANK_bad_file_path ()
#define STR_BANK_unable_open_file   gatb::core::tools::misc::MessageRepository::singleton().BANK_unable_open_file ()
#define STR_BANK_unable_write_file  gatb::core::tools::misc::MessageRepository::singleton().BANK_unable_write_file ()

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_STRINGS_REPOSITORY_HPP_ */
