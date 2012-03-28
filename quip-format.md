
The Quip File Format
====================

Quip compresses FASTQ files with extreme prejudice. This document gives a brief
overview of the file format that it outputs. This is not information that any
user needs or would want to know, but will be useful for anyone intending to
hack on the quip source code.

Conventions


The Quip format consists of a file header, zero or more blocks of compressed
reads, and a footer. Separating the reads into blocks may seem a bit old
fashioned, but it allows us to do fancy things to subsets of the data (e.g., try
to assembly the reads).


Quip Header
-----------

    +---+---+---+---+---+---+---+
    |      Magic Number     | V |
    +---+---+---+---+---+---+---+


The header of a Quip file in seven bytes.  The first six bytes of a Quip file
are a magic number used to quickly identify that the file is of the correct
format. The magic number defined as,

    const uint8_t quip_header_magic[6] = 
        {0xff, 'Q', 'U', 'I', 'P', 0x00};

Which in hexidecimal is,

          FF 51 55 49 50 00

The seventh and last byte of the header is a number used to identify the version
of the file format.


Block Header
------------

    Compressed reads are split into blocks, each block beginning with the
    following header.

    +---+---+---+---+---+---+---+---+
    |  Read Count   |   Base Count  |
    +---+---+---+---+---+---+---+---+

    +---+---+---+---+---+---+---+---+
    |    Read Len.  |    Run Len.   |  ...
    +---+---+---+---+---+---+---+---+

    +---+---+---+---+---+---+
    | Q | U |    Run Len.   |  ...  (Where, Q is the base quality score, and 
    +---+---+---+---+---+---+        U is quality score used to encode 'N's)

    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    | Uncomp. Bytes |  Comp. Bytes  |        CRC64 Checksum         |   (ID chunk description)
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+

    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    | Uncomp. Bytes |  Comp. Bytes  |        CRC64 Checksum         |   (Sequence chunk description)
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+

    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    | Uncomp. Bytes |  Comp. Bytes  |        CRC64 Checksum         |   (Quality chunk description)
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+

    A block header consists of a 4-byte integer giving the number of reads in
    the block, followed by a 4-byte integer giving the number of bases.

    Read-lengths are listed using run-length encoding.

    Guesses at the quality score scheme are encoded using run length encoding.
    
    Following this, 4-byte uncompressed and compressed byte counts and 8-byte
    checksums are given for read IDs, sequences, and quality scores,
    respectively.



Block
-----

    +---+ ... +---+---+ ... +---+---+ ... +---+
    |  Comp. Ids  | Comp. Seqs. | Comp. Quals.|
    +---+ ... +---+---+ ... +---+---+ ... +---+

A block is simply compressed chunks of read ids, sequences, and finally
qualities.


Compressed Id Chunk
-------------------

    +---+ ... +---+
    |    Data     |
    +---+ ... +---+

Compressed ids are pure arithmetic coded data.


Compressed Sequence Chunk
-------------------------

    +---+---+---+---+---+ ... +---+
    |  Contig Len   |    Data     |
    +---+---+---+---+---+ ... +---+

The compressed sequence chunk is lead by a 4-byte unsigned integer giving
the length of the supercontig---a single contig formed by concatenating
each of the seperately assembled contigs.

This number may be zero, indicating that assembly based compression is not
being used.


Compressed Quality Chunk
------------------------

    +---+ ... +---+
    |    Data     |
    +---+ ... +---+

Compressed qualities are also pure arithmetic coded data.


Quip Footer
-----------

The sequence of blocks is terminated by four zero bytes, indicating
an empty block and the end of the stream.

    +---+---+---+---+
    |       0       |
    +---+---+---+---+


