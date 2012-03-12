
http://sourceforge.net/projects/quip-compressor/

Algorithms used in Quip
=======================

Primarily Quip works by building statistical models of the read ids, nucleotide
sequences, and quality scores which are then compressed using arithmetic coding.
In addition, we use a highly efficient de novo assembler, to assemble contiguous
regions and store nucleotide sequences as positions within these contigs.

In the following sections we give specific details. 


Read IDs
--------

Read IDs are compressed as string deltas. Delimited fields and numbers are
parsed from IDs and compared against the previous ID. The algorithm looks for
numbers that increase in sequential manner, and attempts to match as many parts
of the ID as possible.


Quality Scores
--------------

The model for quality scores is informed by three things:
     
    1. The previous three quality scores in the read.
    3. The position of the quality score within the read.
    2. A running delta which gives a rough trajectory of the qualities.


Nucleotide Sequences
--------------------

Nucleotide sequences are modeled primarily using a high-order Markov-chain.

We use the additionally strategy of removing redundancy by performing de-novo
assembly of the first 100,000,000 bases and aligning every read against the
generated contigs. For reads with good alignments, we can save substantial space
by storing the alignment rather than the sequence. This approach is most
effective with sequencing data from small genomes or from a small subset of a
large genome (e.g. RNA-Seq or Chip-Seq).

With DNA sequencing of large genomes we are able to obtain some savings by
assembly long repeat regions. Yet, in most RNA-Seq data we achieve a huge
savings by assembling many full length mRNA transcripts.

With this approach we attempt to combine the benefits of reference-based
compression with those of a stand-alone compression program. Compression and
decompression are very fast, an no database is needed to perform either. This is
of particular importance for sequencing data in which no reference genome
exists: metagenomics, or transcriptome sequencing of species without genome are
obvious targets.

To undertake this, we implemented our own specially tuned short read aligner and
de novo assembler. The latter required new techniques that to our knowledge have
not been previously explored. 


Probabilistic De Bruijn Graphs
------------------------------

Quip uses an exceedingly efficient (time and space) variation of a De Bruijn
graph assembler. In the assembler, the De Bruijn graph is a represented
internally by a table keeping track of the number occurrences of nucleotide
k-mers in the short reads (where k = 25 is typical, but many multiple values
might be used). In nearly every commonly used assembler this is done with a hash
table. Hash tables can be implemented very efficiently, but with a few
sacrifices we can do much better.
 
The assembly algorithm used by Quip works using a probabilistic data structure
based on the Bloom filter, called a "d-left counting Bloom filter", first
described by Bonomi, et al, in:

      Bonomi, F., Mitzenmacher, M., Panigrahy, R., Singh, S., & Varghese, G.
      (2006). An improved construction for counting bloom filters. 14th Annual
      European Symposium on Algorithms, LNCS 4168 (pp. 684–695). Springer.

The data structure works very much like a hash table, but with two exceptions.
First, it uses vastly less memory, and secondly it has a small but non-zero
false positive rate: occasionally an incorrect count will be associated with a
k-mer. Because the only function of the assembler is compression, a small
false-positive rate is acceptable. With the Bloom filter assembly becomes a
practical option for compression: we can assemble one million 100nt reads in
seconds with a few hundred of megabyte, rather than several gigabytes that would
otherwise be needed.


Additional Features
-------------------

To ensure data integrity, quip attaches independent CRC64 checksums to each block
of IDs, sequences, and quality scores. This way, any data corruption can be
not just detected, but localized. The integrity a compressed FASTQ file
can be checked with the 'quip --test' command.

Quip also stores summary information about the compressed reads. With the
'quip --list' command, the number of compressed reads and bases is listed
without needing to decompress the file or count lines.

