/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Copyright (c) 2013                                                      *
 *                                                                           *
 *   GATB is free software; you can redistribute it and/or modify it under   *
 *   the CECILL version 2 License, that is compatible with the GNU General   *
 *   Public License                                                          *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   CECILL version 2 License for more details.                              *
 *****************************************************************************/

/** \file Bank.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Interface definition for genomic databases management
 */

#ifndef _GATB_CORE_BANK__IMPL_BANK_HPP_
#define _GATB_CORE_BANK__IMPL_BANK_HPP_

/********************************************************************************/

#include <gatb/bank/api/IBank.hpp>

#include <vector>
#include <string>

/********************************************************************************/
namespace gatb      {
/** \brief Core package of the GATP project.
 *
 * The gatb::core package holds all the fundamental packages needed for writting
 * assembly algorithms.
 *
 * It holds some generic tools, like operating system abstraction, collections management or design patterns
 * concerns. It also holds recurrent needs as reading genomic banks, handling kmers and so on.
 */
namespace core      {
/** \brief Package for genomic databases management. */
namespace bank      {
/** \brief Implementation for genomic databases management. */
namespace impl      {
/********************************************************************************/

/** \brief Interface for reading genomic databases.
 *
 * Sample of use:
 * \snippet bank1.cpp  snippet1_bank
 */
class Bank : public IBank
{
public:

    /** Constructor.
     * \param[in] filenames : uri list of the banks. */
    Bank (const std::vector<std::string>& filenames);

    /** Constructor.
     * \param[in] argc : number of filenames
     * \param[in] argv : filenames
     */
    Bank (int argc, char* argv[]);

    /** Constructor.
     * \param[in] filename : uri of the bank. */
    Bank (const std::string& filename);

    /** Destructor. */
    ~Bank ();

    /** \copydoc IBank::iterator */
    tools::dp::Iterator<Sequence>* iterator ()  { return new Iterator (*this); }

    /** \copydoc IBank::insert */
    void insert (const Sequence& item);

    /** \copydoc IBank::flush */
    void flush ()  {}

    /** \copydoc IBank::getSize */
    u_int64_t getSize ()  { return filesizes; }

    /** \copydoc IBank::estimateNbSequences */
    u_int64_t estimateNbSequences ();

    /** \copydoc IBank::estimateNbSequences */
    size_t estimateMaxSequenceLength ();

    /**
     * \return maximum number of files. */
    static const size_t getMaxNbFiles ()  { return 30; }

    /************************************************************/

    /** \brief Specific Iterator impl for Bank class
     *
     * This implementation relies on the initial code from Minia. It wraps the
     * Iterator interface with the Minia code.
     *
     * Note that we made some effort not to put implementation code
     * here in the header; see in particular buffered_file and buffered_strings
     * attributes whose type is void* (and not the implementation type defined in
     * the cpp file).
     *
     * Note the small improvement compared to Minia: it is possible to create an
     * Iterator that provides (or not) sequence comments, according to the corresponding
     * parameter given to the Iterator constructor.
     *
     *  <b>IMPROVEMENTS</b>:
     *  - in case we have several banks to read, we could have at one time only one stream opened on the currently
     *  iterated file. The current implementation opens all streams, which may be avoided.
     */
    class Iterator : public tools::dp::Iterator<Sequence>
    {
    public:

        /** Define what kind of comment we want to retrieve. Such comments can be retrieved through
         * gatb::core::bank::Sequence
         */
        enum CommentMode_e
        {
            /** Empty comments are provided to clients. */
            NONE,
            /** Comments with only the FASTA ID are provided to clients. \n
             *  Ex: 'ENSTTRP00000009639pep:novel'*/
            IDONLY,
            /** Full comments are provided to clients. \n
             *  Ex: 'ENSTTRP00000001236pep:novel genescaffold:turTru1:GeneScaffold_3311:182575:189152:1'*/
            FULL
        };

        /** Constructor.
         * \param[in] ref : the associated iterable instance.
         * \param[in] commentMode : kind of comments we want to retrieve
         */
        Iterator (Bank& ref, CommentMode_e commentMode = FULL);

        /** Destructor */
        ~Iterator ();

        /** \copydoc tools::dp::Iterator::first */
        void first();

        /** \copydoc tools::dp::Iterator::next */
        void next();

        /** \copydoc tools::dp::Iterator::isDone */
        bool isDone ()  { return _isDone; }

        /** \copydoc tools::dp::Iterator::item */
        Sequence& item ()     { return *_item; }

        /** Estimation of the number of sequences. Used (by delegation) by the Bank class.
         * \return the sequences number estimation. */
        u_int64_t estimateNbSequences ();

    private:

        /** Reference to the underlying Iterable instance. */
        Bank&    _ref;

        /** Tells what kind of comments we want as a client of the iterator. */
        CommentMode_e  _commentsMode;

        /** Tells whether the iteration is finished or not. */
        bool _isDone;

        /** Initialization method. */
        void init ();

        /** Finish method. */
        void finalize ();

        int   index_file; // index of current file

        void** buffered_file;
        void*  buffered_strings;   // variable_string_t *read, *dummy, *header;

        bool get_next_seq           (tools::misc::Vector<char>& data);
        bool get_next_seq           (tools::misc::Vector<char>& data, char **cheader, int *hlen, CommentMode_e mode);

        bool get_next_seq_from_file (tools::misc::Vector<char>& data, int file_id);
        bool get_next_seq_from_file (tools::misc::Vector<char>& data, char **cheader, int *hlen, int file_id, CommentMode_e mode);
    };

protected:

    friend class Iterator;

    /** List of URI of the banks. */
    std::vector<std::string> _filenames;

    u_int64_t filesizes;  // estimate of total size for all files
    int       nb_files;   // total nb of files

    /** File handle for inserting sequences into the bank. */
    FILE* _insertHandle;

    /** Initialization method (compute the file sizes). */
    void init ();
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK__IMPL_BANK_HPP_ */
