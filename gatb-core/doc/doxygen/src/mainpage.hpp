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

/** \mainpage GATB Core Documentation
 *
 ****************************************************************************************************
 ****************************************************************************************************
 * \section intro What is GATB ?
 *
 *  GATB means "Genome Assembly Tool Box".
 *
 *  It can be seen as software components providing services for assembling and analyzing genomes of various kinds.
 *
 * The current documentation is for release 1.0.6. You can download it <a href="http://gatb.inria.fr/binaries-url">here</a>.
 *
 ****************************************************************************************************
 ****************************************************************************************************
 * \section purpose Purpose of the GATB core library
 *
 *
 * First, gatb::core holds a C++ library that provides an easy way to create bioinformatics programs.
 *
 * It supports the following operations natively:
 *  - FASTA/FASTQ parsing
 *  - K-mer counting
 *  - Minimizer computation of k-mers, partitioning of datasets by minimizers
 *  - de Bruijn graph construction
 *  - de Bruijn graph traversal operations (contigs, unitigs)
 *
 *  One structure is central to the GATB project: the <a href="http://en.wikipedia.org/wiki/De_Bruijn_graph">De Bruijn graph</a>.
 *  This sort of data structure is today widely used in NGS software (like assembly software).
 *
 * First, gatb::core holds a C++ library that provides an easy way to handle a central structure in the
 * GATB project: the <a href="http://en.wikipedia.org/wiki/De_Bruijn_graph">De Bruijn graph</a>.
 * This kind of data structure is today widely used in
 * <a href="http://en.wikibooks.org/wiki/Next_Generation_Sequencing_(NGS)">NGS</a>
 * software (like assembly software).
 *
 * More precisely, the library provides means to build and use De Bruijn graphs with a low memory footprint,
 * which comes initially from the <a href="http://minia.genouest.org">minia</a> assembly tool.
 *
 * The documentation you are reading is the official documentation of the gatb::core reference API. The
 * audience is therefore developers interested in creating bio-informatics software.
 *
 * The GATB project has been published in <a href="http://bioinformatics.oxfordjournals.org/content/early/2014/06/30/bioinformatics.btu406.abstract?keytype=ref&ijkey=koXTzqbf4Nt1kVO">BioInformatics in 2014</a>.
 *
 ****************************************************************************************************
 ****************************************************************************************************
 * \section struct Services provided by the GATB core library
 *
 * From the client point of view, the gatb::core package provides:
 *   - libraries that offer low level genomic operations, up to the De Bruign graph creation
 *   - tests of the libraries
 *   - snippets showing how to use the library
 *   - specific binaries that rely on the libraries
 *   - wrappers of the libraries services for several langages (java, python, ...)
 *
 * You will find here the code documentation for namespaces, classes, methods of the different
 * components that compose the <b>gatb::core</b> design.
 *
 ****************************************************************************************************
 ****************************************************************************************************
 * \section howto How can I use the GATB core library ?
 *
 * As a starting point, it is strongly recommended to have a look at \ref snippets_page. You will find
 * there information about the compilation process and how to create a new project based on gatb::core.
 *
 * You will find also a lot of snippets showing gatb::core
 * in action.
 *
 ****************************************************************************************************
 ****************************************************************************************************
 * \section contact Contact
 *
 * You can get support on the BioStars forum <a href="http://www.biostars.org/p/101393">here</a>.
 *
 * You can also have general information about the <a href="http://gatb.inria.fr">GATB project</a>. You will
 * find <a href="http://gatb.inria.fr/tutorials">here</a> high level tutorials about GATB.
 *
 ****************************************************************************************************
 ****************************************************************************************************
 * \section other Other material
 *
 * You can also read the related pages:
 *      - \ref download_page
 *      - \ref design_page
 *      - \ref dbgh5_page
 */
