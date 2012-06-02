
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

    +---+---+---+---+---+---+---+---+
    |      Magic Number     | V | F |
    +---+---+---+---+---+---+---+---+


The header of a Quip file in seven bytes.  The first six bytes of a Quip file
are a magic number used to quickly identify that the file is of the correct
format. The magic number defined as,

    const uint8_t quip_header_magic[6] = 
        {0xff, 'Q', 'U', 'I', 'P', 0x00};

Which in hexidecimal is,

          FF 51 55 49 50 00

The seventh, penultimate byte of the header `V` is a number used to identify
the version of the file format.

Lastly, the F byte gives flags is a bitwise OR of flags determining how the
file was compressed:
    0:   whether the compression is reference-based
    1:   whether de novo assembly of unaligned reads was used
    2-7: reserved for future use

If reference-based compression was used, the next 8-bytes gives a hash of the
reference sequence, to prevent the incorrect reference sequence being used in
decompression.

    +---+---+---+---+---+---+---+---+
    |  Reference Sequence Checksum  |
    +---+---+---+---+---+---+---+---+

This is followed by the file-name of the reference sequence, which is preceded by its length.

    +---+---+---+---+---+---+---+---+---+- ... -+
    | Ref. Name Len | Reference File Name       |
    +---+---+---+---+---+---+---+---+---+- ... -+

We then list an index of all the sequences present in the reference along with
their lengths.

    +---+---+---+---+
    |   Num. Seqs.  |
    +---+---+---+---+

    +---+---+---+---+---+---+---+---+- ... -+---+---+---+---+---+---+---+---+
    | Seq. Name Len.| Sequence Name         |       Sequence Length         |
    +---+---+---+---+---+---+---+---+- ... -+---+---+---+---+---+---+---+---+

    ...

None of this information is strictly necessary for decompression, but very useful
when one needs to locate the correct reference sequence.


Auxiliary Data
--------------

This is followed by axillary data used in the compression. In particular, the
header of a compressed SAM/BAM file.

    +---+---+----+---+---+---+---+---+---+---+ ... +---+
    | T |        Aux. Data Length        |  Aux. Data  |
    +---+---+----+---+---+---+---+---+---+---+ ... +---+

Where `T` is a code giving the source format of the auxiliary data.



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

    +---+---+---+---+---+
    | Q |    Run Len.   |  ...  (Where Q is assumed the base quality score.)
    +---+---+---+---+---+        

    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    | Uncomp. Bytes |  Comp. Bytes  |        CRC64 Checksum         |   (ID chunk description)
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+

    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    | Uncomp. Bytes |  Comp. Bytes  |        CRC64 Checksum         |   (Aux chunk description)
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


Compressed Aux Chunk
--------------------

    +---+ ... +---+
    |    Data     |
    +---+ ... +---+

Compressed ids are pure arithmetic coded data.


Compressed Sequence Chunk
-------------------------

    +---+ ... +---+
    |    Data     |
    +---+ ... +---+

Compressed sequences are pure arithmetic coded data.


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


