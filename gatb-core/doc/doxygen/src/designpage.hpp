/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  R.Chikhi, G.Rizk, E.Drezen
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

/**\page design_page Design of the library
 *
 ****************************************************************************************************
 ****************************************************************************************************
 * \section coredesign Global design of gatb::core package
 *
 ****************************************************************************************************
 * \subsection coredesign_intro Introduction
 *
 * As its name suggests, the gatb::core package holds all the fundamental packages needed for writting
 * assembly algorithms.
 *
 * It holds some generic services, like:
 *      - operating system abstraction,
 *      - collections management
 *      - design patterns
 *      - genomic banks read/write
 *      - kmers management
 *      - De Bruijn graph management
 *
 * \image html CoreGlobal.png "Global view of the gatb::core package"
 *
 *
 ****************************************************************************************************
 * \subsection coredesign_system  Package gatb::core::system
 *
 * The gatb::core::system package holds operations related to operating system.
 *
 * Making such an abstraction is interesting when considering the following points:
 * - client code should rely on interfaces defined here and not on specific implementations (for linux, macos, ...)
 * - client code should not use compilation flags for one specific system (like \#ifdef __LINUX__)
 * - improved operating system operations should be available to clients without modification in their code
 *
 * The idea is so to list all necessary operating system operations that are likely to be used, to make several
 * interfaces that gather some operations set and to provide specific implementations for several system (linux,
 * macos, windows).
 *
 * \image html SystemGlobal.png "Global view of the gatb::core::system package"
 *
 ****************************************************************************************************
 * \subsection coredesign_tools  Package gatb::core::tools
 *
 * The gatb::core::tools package holds generic operations used throughout the project. They are not specific
 * to genome nor assembly concerns and are defined here in a separate package.
 *
 * Roughly speaking, we may find the following operations:
 *      - design patterns tools (iterator, observer, smart pointer, etc...)
 *      - collections (container, bag, iterable, etc...)
 *      - misc
 *
 * \image html ToolsGlobal.png "Global view of the gatb::core::tools package"
 *
 *
 *****************************************************************************************************
 * \subsection coredesign_bank Package gatb::core::bank
 *
 * The gatb::core::bank package holds operations related to genomic databases management.
 *
 * In particular, we define a gatb::core::bank::IBank interface that is an iterable over gatb::core::bank::Sequence
 * instances. This interface has a lot of implementations, including FASTA format.
 *
 *
 ******************************************************************************************************
 * \subsection coredesign_kmer Package gatb::core::kmer
 *
 * The gatb::core::kmer package holds operations related to kmer management.
 *
 * A typical case is to create a IBank instance and then an iterator on it; such an iterator can feed the
 * kmer model for iterating all the kmer from the sequences.
 *
 *
 *******************************************************************************************************
 * \subsection coredesign_debruijn Package gatb::core::debruijn
 *
 * The gatb::core::debruijn package holds operations related to De Bruijn graph.
 *
 * The main API is the gatb::core::debruijn::impl::Graph class. It holds several features:
 *      - create graph
 *      - load graph
 *      - iterate nodes of the graph
 *      - get neighbors of a node
 *
 * The implementation of the class relies on a probabilistic data structure named
 * <a href="http://en.wikipedia.org/wiki/Bloom_filter">Bloom filter</a>. Its usage allows to
 * have a very low memory footprint, even for big graphs.
 */
