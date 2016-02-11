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
 * GATB means "Genome Analysis Toolbox with de-Bruijn graph".
 * 
 * The GATB-CORE project provides a set of highly efficient algorithms to analyse NGS data sets. These 
 * methods enable the analysis of data sets of any size on multi-core desktop computers, including very 
 * huge amount of reads data coming from any kind of organisms such as bacteria, plants, animals and 
 * even complex samples (e.g. metagenomes). More: https://project.inria.fr/gatb/.
 *
 * GATB is made two main parts:
 *
 * - the GATB-CORE library: for development purpose, GATB-CORE enables the creation of new software tools
 * - the GATB-Tools: contains ready-to-use softwares relying on GATB-CORE. More <a href="https://project.inria.fr/gatb/software/">here</a>.
 *
 * The GATB project has been published in <a href="http://bioinformatics.oxfordjournals.org/content/30/20/2959">BioInformatics in 2014</a>.
 * There are also several publications about GATB use cases and tools available <a href="https://project.inria.fr/gatb/publications/">here</a>.
 * 
 ****************************************************************************************************
 ****************************************************************************************************
 * \section purpose Purpose of the GATB core library
 *
 * gatb::core is a high-performance and low memory footprint C++ library.
 *
 * It supports the following operations natively:
 *  - FASTA/FASTQ parsing
 *  - K-mer counting
 *  - Minimizer computation of k-mers, partitioning of datasets by minimizers
 *  - de Bruijn graph construction
 *  - de Bruijn graph traversal operations (contigs, unitigs)
 *
 * One structure is central to the GATB project: the <a href="http://en.wikipedia.org/wiki/De_Bruijn_graph">De Bruijn graph</a>.
 * This sort of data structure is today widely used in NGS software (like assembly softwares).
 *
 * So, one can say that GATB-CORE library provides means to build and use De Bruijn graphs with a low memory footprint,
 * which comes initially from the <a href="https://project.inria.fr/gatb/software/minia/">minia</a> assembly tool.
 *
 * The documentation you are reading is the official documentation of the gatb::core reference API. The
 * audience is therefore developers interested in creating bioinformatics softwares.
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
 * \section howto How can I make a new software using GATB core library ?
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
 * You can also have general information about the <a href="https://project.inria.fr/gatb/">GATB project</a>. You will
 * find <a href="https://project.inria.fr/gatb/tutorials">here</a> high level tutorials about GATB.
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
