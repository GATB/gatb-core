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

/**\page dbgh5_page Software dbgh5
 *
 ****************************************************************************************************
 ****************************************************************************************************
 * \section dbgh5_intro Purpose of dbgh5
 *
 * The dbgh5 binary creates a De Bruijn graph from an input bank (likely reads).
 *
 * The output graph is saved in a <a href="http://www.hdfgroup.org/HDF5/">HDF5</a> file.  Such a file can be
 * read through the gatb::core API by other tools. A classical workflow is the following:
 *
 * \image html Transfo.png "From reads to HDF5"
 *
 *****************************************************************************************************
 ****************************************************************************************************
 * \section dbgh5_available Availability
 *
 * dbgh5 is provided in the gatb::core distribution (see download page). You can get directly a binary
 * or you can also compile it by yourself.
 *
 ******************************************************************************************************
 ****************************************************************************************************
 * \section dbgh5_algorithm General algorithm
 *
 * dbgh5 does all or part of the following steps:
 *      - counts the kmers of the input bank
 *      - put these kmers into a Bloom filter
 *      - compute an additional structure called cFP (for critical false positives)
 *      - locate some specific nodes in the De Bruijn graph
 *
 * After these steps, one has a HDF5 file that holds the De Bruijn graph information. It can be then
 * read through the gatb::core API.
 *
 * ******************************************************************************************************
 ****************************************************************************************************
 * \section dbgh5_options Options
 *
 * dbgh5 is highly configurable; the minimal requirement is the input set of reads, so you can just type for instance:
 * \code
 * dbgh5 -in /somewhere/myreads.fasta
 * \endcode
 *
 * you can just type 'dbgh5' and you will get all the available options.
 * \code
[graph options]

   [kmer count options]
          -in              (1 arg) :    reads file
          -kmer-size       (1 arg) :    size of a kmer  [default '31']
          -abundance-min   (1 arg) :    min abundance threshold for solid kmers  [default '3']
          -abundance-max   (1 arg) :    max abundance threshold for solid kmers  [default '4294967295']
          -histo-max       (1 arg) :    max number of values in kmers histogram  [default '10000']
          -solidity-kind   (1 arg) :    way to compute solids (sum, min or max)  [default 'sum']
          -solid-kmers-out (1 arg) :    output file for solid kmers  [default '']
          -out             (1 arg) :    output file  [default '']
          -out-dir         (1 arg) :    output directory  [default '.']
          -minimizer-type  (1 arg) :    minimizer type (0=lexi, 1=freq)  [default '0']
          -minimizer-size  (1 arg) :    size of a minimizer  [default '8']

   [bloom options]
          -bloom        (1 arg) :    bloom type ('basic', 'cache', 'neighbor')  [default 'neighbor']
          -debloom      (1 arg) :    debloom type ('none', 'original' or 'cascading')  [default 'cascading']
          -debloom-impl (1 arg) :    debloom impl ('basic', 'minimizer')  [default 'minimizer']

   [branching options]
          -branching-nodes (1 arg) :    branching type ('none' or 'stored')  [default 'stored']
          -topology-stats  (1 arg) :    topological information level (0 for none)  [default '0']

   [general options]
          -nb-cores          (1 arg) :    number of cores  [default '0']
          -max-disk          (1 arg) :    max disk   (in MBytes)  [default '0']
          -max-memory        (1 arg) :    max memory (in MBytes)  [default '2000']
          -verbose           (1 arg) :    verbosity level  [default '1']
          -email             (1 arg) :    send statistics to the given email address  [default '']
          -email-fmt         (1 arg) :    'raw' or 'xml'  [default 'raw']
 * \endcode
 *
 *
 */
