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
 *  It can be seen as software components providing services for assembling genomes of various kinds.
 *
 ****************************************************************************************************
 ****************************************************************************************************
 * \section purpose Purpose of the GATB core library
 *
 * The documentation you are reading is about the GATB core library.
 *
 * One purpose of this library is to provide an easy way to handle a central structure in the
 * GATB project: the <a href="http://en.wikipedia.org/wiki/De_Bruijn_graph">De Bruijn graph</a>
 *
 * More precisely, the intent is to provide means to build De Bruijn graphs with a low memory footprint,
 * which comes initially from the <a href="http://minia.genouest.org">minia</a> assembly tool.
 *
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
 * As a starting point, it is strongly recommended to have a look at \ref snippets_page
 *
 ****************************************************************************************************
 ****************************************************************************************************
 * \section contact Contact
 *
 * You can get support on the following <a href="mailto:gatb-core-support@lists.gforge.inria.fr">mailing list</a>
 *
 * You can also have general information about the <a href="http://gatb.inria.fr">GATB project</a>
 *
 ****************************************************************************************************
 ****************************************************************************************************
 * \section other Other material
 *
 * You can also read the related pages:
 *      - \ref download_page
 *      - \ref design_page
 *      - \ref tests_page
 */
