emphf
=====

`emphf` is a minimal perfect hashing library for large-scale key sets focused on
speed and low memory usage. A minimal perfect hash function (MPHF) is a data
structure that maps injectively a static set of strings of size `n` to the
integer set `[0, n-1]`. The overall space usage of the MPHFs generated with
`emphf` is `2.61 n` bits plus a small constant factor.

The algorithms used in this library are described in the paper
[*Cache-Oblivious Peeling of Random Hypergraphs*](http://arxiv.org/abs/1312.0526)
by Djamal Belazzougui, Paolo Boldi, Giuseppe Ottaviano, Rossano Venturini, and
Sebastiano Vigna.

All the algorithms implemented here construct a MPHF for a given key set using a
standard scheme based on the Majewski–Wormald–Havas–Czech (MWHC) construction,
that is based on finding a peeling order of a random 3-hypergraph.

* `compute_mphf_seq` computes the peeling order using the standard in-memory
  linear-time peeling algorithm, implemented with the xor-trick described in the
  paper. This is the fastest MPHF construction implementation that we know of.

* `compute_mphf_scan_mmap` computes the peeling order using the new algorithm
  proposed in the paper, using an mmapped temporary file for its data structure,
  thus in an external-memory setting.

* `compute_mphf_hem` uses a technique described in F. C. Botelho, R. Pagh, and
  N. Ziviani, “Practical perfect hashing in nearly optimal space”, that splits
  the key sets in bucket that are small enough that the perfect hash function is
  easy to compute in memory for each bucket. While the construction is faster
  than by using `compute_mphf_scan_mmap`, the resulting data structure takes
  slightly more space and it is slightly slower.

`compute_mphf_seq` and `compute_mphf_scan_mmap` construct the same data
structure, but the second should be used when the space needed to construct the
MPHF does not fit in main memory. The obtained MPHF data structures can perform
lookups significantly faster than other implementations using the same algorithm
(see the paper for details).

The script `test_all.py` executes all the algorithms on a given file (or the
standard UNIX dictionary if none is provided) and checks that the generated hash
function is indeed minimal and perfect.

The project uses CMake. To build it on Unix systems it should be sufficient to
do the following:

    $ cmake .
    $ make
