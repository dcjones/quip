![Quip](http://cs.washington.edu/homes/dcjones/quip/logo.png)\

Quip compresses next-generation sequencing data with extreme prejudice. It
supports input and output in the FASTQ and SAM/BAM formats, compressing large
datasets to as little as 15% of their original size.

Installation
============

Binaries
--------

In the future will will provide binaries for a number of operating systems. In
the mean time, you must install from source.


Source from GitHub
------------------

To install from github, you will need a C compiler, as well as a relatively recent
version of [automake](http://www.gnu.org/software/automake/) and
[autoconf](http://www.gnu.org/software/autoconf/).

First clone the git repository with,
  
    git clone git://github.com/dcjones/quip.git

Enter the quip directory
  
    cd quip

Generate the configure script using

    autoreconf -i

Then configure, compile, and install

    ./configure && make install



Source from a tarball
---------------------

You will need only a C compiler to install from a source tarball.

Extract the source tarball

    tar xzf quip-1.1.0.tar.gz

Enter the quip directory

    cd quip-1.1.0

Then configure, compile, and install

    ./configure && make install



Usage
=====

Quip mostly works the same as `gzip` or `bzip2`. For example, if you do something like,

    quip reads.fastq

You will get a file called `reads.fastq.qp` that is significantly smaller than
the original. Note that unlike `gzip` or `bzip2`, this will not delete the
original file.

For more details, see `man quip` after installing.


Algorithms
==========

Quip implements an ensemble of compression techniques specially built to
compress sequencing reads as much as possible. The basis of the algorithm is
statistical compression of read identifiers, quality scores, and nucleotide
sequences using arithmetic coding. In addition, we implement reference-based
compression in which aligned reads are stored a positions within a genome.
And, when no reference genome is available, an extremely efficient de novo
assembly algorithm can transparently construct one.

In the following sections we give specific details. 


Read IDs
--------

The only requirement of read identifiers is that they uniquely identify the
read. A single integer would do, but typically each read comes with a complex
string containing the instrument name, run identifier, flow cell identifier,
and tile coordinates. Much of this information is the same for every read and
is simply repeated, inflating the file size.

To remove this redundancy, we use a form of delta encoding. A  parser
tokenizes the ID into separate fields which are then compared to the previous
ID. Tokens that remain the same from read to read (e.g. instrument name)
can be compressed to a negligible amount of space --- arithmetic coding produces
codes of less than 1 bit in such cases. Numerical tokens are recognized and
stored efficiently, either directly or as an offset from the token in the same
position in previous read. Otherwise non-identical tokens are encoded by
matching as much of the prefix as possible to the previous read's token before
directly encoding the non-matching suffix.

The end result is that read IDs, which are often 50 bytes or longer, are
typically stored in 2-4 bytes. Notably, in reads produced from Illumina
instruments, most parts of the ID can be compressed to consume almost no
space; the remaining few bytes are accounted for by tile coordinates. These
coordinates are almost never needed in downstream analysis, so removing them
as a preprocessing step would shrink file sizes even further. The parser used
is suitably general so that no change to the compression algorithm would be
needed.

Nucleotide Sequences
--------------------

To compress nucleotide sequences, we adopt a very simple model based on high-order
Markov chains. The nucleotide at a given position in a read is predicted
using the preceding twelve positions. This model uses more memory than
traditional general-purpose compression algorithms (4^13 = 67,108,864
parameters are needed, each represented in 32 bits)  but it is simple and
extremely efficient (very little computation is required and run time is
limited primarily by memory latency, as lookups in such a large table result
in frequent cache misses).

An order-12 Markov chain also requires a very large amount of
data to train, but there is no shortage with the datasets we wish to
compress. Though less adept at compressing extremely short files,
after compressing several million reads the parameters are tightly
fit to the nucleotide composition of the dataset so that the remaining reads
will be highly compressed. Compressing larger files only results in a tighter
fit and higher compression.


Quality Scores
--------------

It has been previously noted that the quality score at a given position is
highly correlated with the score at the preceding position (Kozanitis 2011).
This makes a Markov chain a natural model, but unlike
nucleotides, quality scores are over a much larger alphabet (typically 41--46
distinct scores). This limits the order of the Markov chain: long chains will
require a great deal of space and take a unrealistic amount of data to train.

To reduce the number of parameters, we use an order-3 Markov chain, but
coarsely bin the distal two positions. In addition to the preceding three
positions, we condition on the position within the read and a running count of
the number large jumps in quality scores between adjacent positions (where
a ``large jump'' is defined as `|q_{i} - q_{i-1}| > 1`), which allows reads
with highly variable quality scores to be encoded using separate models. Both
of these variables are binned to control the number of parameters.


Reference-based Compression
---------------------------

We have also implemented lossless reference-based compression. Given aligned
reads in SAM or BAM format, and the reference sequence to which they are
aligned (in FASTA format), the reads are compressed preserving all information
in the SAM/BAM file, including the header, read IDs, alignment information,
and all optional fields allowed by the SAM format. Unaligned reads are
retained and compressed using the Markov chain model.

Assembly-based Compression
--------------------------

To complement the reference-based approach, we developed an assembly-based
approach which offers some of the advantages of reference-based compression,
but requires no external sequence database and produces files which are
entirely self-contained. We use the first (by default) 2.5 million reads to
assemble contigs which are then used in place of a reference sequence
database to encode aligned reads compactly as positions.

Once contigs are assembled, read sequences are aligned using a simple ``seed
and extend'' method: 12-mer seeds are matched using a hash table, and
candidate alignments are evaluated using Hamming distance. The best
(lowest Hamming distance) alignment is chosen, assuming it falls below a given
cutoff, and the read is encoded as a position within the contig set. Roughly,
this can be thought of as a variation on the Lempel-Ziv algorithm: as
sequences are read, they are matched to previously observed data, or in this
case, contigs assembled from previously observed data. These contigs are
not explicitly stored, but rather reassembled during decompression.


Traditionally, de novo assembly is extremely computationally intensive. The
most commonly used technique involves constructing a de Bruijn graph, a
directed graph in which each vertex represents a nucleotide k-mer present in
the data for some fixed k (e.g., k = 25 is a common choice). A directed
edge from a k-mer u to v occurs if and only if the (k - 1)-mer suffix
of u is also the prefix of v. In principle, given such a graph, an
assembly can be produced by finding an Eulerian path, i.e., a path that
follows each edge in the graph exactly once (Pevzner 2001). In practice,
since NGS data has a non-negligible error rate, assemblers augment each vertex
with the number of observed occurrences of the k-mer and leverage these
counts using a variety of heuristics to filter out spurious paths.

A significant bottleneck of the de Bruijn graph approach is building an
implicit representation of the graph by counting and storing k-mer
occurrences in a hash table. The assembler implemented in Quip 
overcomes this bottleneck to a large extent by using a data structure based on
the Bloom filter to count k-mers. The Bloom filter (Bloom 1970) is a
probabilistic data structure that represents a set of elements extremely
compactly, at the cost of elements occasionally colliding and incorrectly
being reported as present in the set. It is probabilistic in the sense that
these collisions occur pseudo-randomly, determined by the size of the table
and the hash functions chosen, but generally with low probability.

The Bloom filter is generalized in the counting Bloom filter, in which an
arbitrary count can be associated with each element (Fan 2000), and
further refined in the d-left counting Bloom filter (dlCBF)
(Bonomi 2006), which requires significantly less space than the already
quite space efficient counting Bloom filter. We base our assembler on a
realization of the dlCBF. Because we use a probabilistic data structure,
k-mers are occasionally reported to have incorrect (inflated) counts. The
assembly can be made less accurate by these incorrect counts, however a poor
assembly only results in slightly increasing the compression ratio.
Compression remains lossless regardless of the assembly quality, and in
practice collisions in the dlCBF occur at a very low rate (this is explored in
the results section).

Given a probabilistic de Bruijn graph, we assemble contigs using a very simple
greedy approach. A read sequence is used as a seed and extended on both ends
one nucleotide at a time by repeatedly finding the most abundant k-mer that
overlaps the end of the contig by k-1 bases. More sophisticated heuristics
have been developed, but we choose to focus on efficiency, sacrificing a
degree of accuracy.

Counting k-mers efficiently with the help of Bloom filters was previously
explored by (Melsted 2011), who use it in addition, rather than
in place, of a hash table. The Bloom filter is used as a ``staging area'' to
store k-mers occurring only once, reducing the memory required by the hash
table. Concurrently with our work, (Pell 2011) have also developed a
probabilistic de Bruijn graph assembler, but using a non-counting Bloom
filter. While they demonstrate a significant reduction in memory use, unlike
other de Bruijn graph assemblers, only the presence or absence of a k-mer is
stored, not its abundance, which is essential information when the goal is
producing accurate contigs.



Additional Features
-------------------

In designing the file format used by Quip, we included several useful features
to protect data integrity. First, output is divided into blocks of several
magabytes each. In each block a separate 64 bit checksum is computed for read
identifiers, nucleotide sequences, and quality scores. When the archive is
decompressed, these checksums are recomputed on the decompressed data and
compared to the stored checksums, verifying the correctness of the output. The
integrity of an archived dataset can also be checked with the `quip --test`
command.

Apart from data corruption, reference-based compression creates the
possibility of data loss if the reference used for compression is lost, or an
incorrect reference is used. To protect against this, Quip files store a 64
bit hash of the reference sequence, ensuring that the same sequence is used
for decompression. To assist in locating the correct reference, the file name,
and the lengths and names of the sequences used in compression are also stored
and accessible without decompression.

Additionally, block headers store the number of reads and bases compressed in
the block, allowing summary statistics of a dataset to be listed without
decompression using the `quip --list` command.


References
==========

Christos Kozanitis, Chris Saunders, Semyon Kruglyak, Vineet Bafna, and George
Varghese. Compressing genomic sequence fragments using SlimGene. Journal of
Computational Biology : a Journal of Computational Molecular Cell Biology,
18(3):401–13, March 2011. ISSN 1557-8666. doi: 10.1089/cmb.2010.0253.

P.A.Pevzner, H.Tang, and M.S.Waterman. An Eulerian path approach to DNA fragment
assembly. Proceedings of the National Academy of Sciences of the United States
of America, 98(17):9748–53, August 2001. ISSN 0027-8424. doi:
10.1073/pnas.171285098.

Li Fan, Pei Cao, J. Almeida, and A.Z. Broder. Summary cache: a scalable wide-
area Web cache sharing protocol. IEEE/ACM Transactions on Networking,
8(3):281–293, June 2000. ISSN 10636692. doi: 10.1109/90.851975.

Flavio Bonomi, Michael Mitzenmacher, and Rina Panigrahy. An improved
construction for counting Bloom filters. 14th Annual European Symposium on
Algorithms, LNCS 4168, pages 684–695, 2006.

Pall Melsted and Jonathan K Pritchard. Efficient counting of k-mers in DNA
sequences using a Bloom Filter. BMC Bioinformatics, 12(1):333, 2011. ISSN
1471-2105. doi: 10.1186/1471-2105-12-333.

Jason Pell, Arend Hintze, R Canino-Koning, and Adina Howe. Scaling metagenome
sequence assembly with probabilistic de Bruijn graphs. Arxiv preprint arXiv:,
I(1):1–11, 2011.
