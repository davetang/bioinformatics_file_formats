## Table of Contents

- [Bioinformatics file formats](#bioinformatics-file-formats)
  - [GFF/GTF File Format](#gffgtf-file-format)
  - [FASTA File Format](#fasta-file-format)
  - [FASTQ File Format](#fastq-file-format)
  - [SAM File Format](#sam-file-format)
  - [BED File Format](#bed-file-format)
  - [VCF File Format](#vcf-file-format)

# Bioinformatics file formats

![https://xkcd.com/927/](assets/standards.png)

I lose track of all the bioinformatics file formats, so finally I'll keep track of them here.

## GFF/GTF File Format

<https://www.ensembl.org/info/website/upload/gff.html?redirect=no>

Fields must be **tab-separated**. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'.

* `seqname` - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
* `source` - name of the program that generated this feature, or the data source (database or project name)
* `feature` - feature type name, e.g. Gene, Variation, Similarity
* `start` - Start position* of the feature, with sequence numbering starting at 1.
* `end` - End position* of the feature, with sequence numbering starting at 1.
* `score` - A floating point value.
* `strand` - defined as + (forward) or - (reverse).
* `frame` - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
* `attribute` - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

## FASTA File Format

<https://www.ncbi.nlm.nih.gov/genbank/fastaformat/>

A text-based format for representing nucleotide or amino acid sequences. Each sequence begins with a description line (defline) followed by lines of sequence data.

* `>` - The description line starts with a greater-than symbol
* `sequence identifier` - A unique identifier for the sequence immediately following the '>'
* `description` - Optional description text after the identifier, separated by a space
* `sequence lines` - One or more lines of sequence data (typically 60-80 characters per line)

Sequences can contain standard IUPAC nucleotide codes (A, C, G, T, U, N) or amino acid codes. Gaps are represented by '-'.

## FASTQ File Format

<https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#fastq-files>

A text-based format for storing both biological sequences and their corresponding quality scores. Each record consists of four lines.

* `line 1` - Begins with '@' followed by a sequence identifier and optional description
* `line 2` - The raw sequence letters (A, C, G, T, N)
* `line 3` - Begins with '+' optionally followed by the same sequence identifier
* `line 4` - Quality scores encoded as ASCII characters (same length as line 2)

Quality encoding varies by platform. Illumina 1.8+ uses Phred+33 encoding where ASCII character 33 ('!') represents quality 0.

## SAM File Format

<https://samtools.github.io/hts-specs/SAMv1.pdf>

Sequence Alignment/Map format for storing read alignments against a reference sequence. Fields are **tab-separated**.

* `QNAME` - Query template name (read name)
* `FLAG` - Bitwise flag indicating alignment properties (e.g., 0x1=paired, 0x4=unmapped, 0x10=reverse strand)
* `RNAME` - Reference sequence name (chromosome)
* `POS` - 1-based leftmost mapping position
* `MAPQ` - Mapping quality (Phred-scaled probability that the mapping position is wrong)
* `CIGAR` - CIGAR string describing alignment (M=match, I=insertion, D=deletion, N=skipped, S=soft clip, H=hard clip)
* `RNEXT` - Reference name of the mate/next read ('=' if same as RNAME, '*' if unavailable)
* `PNEXT` - Position of the mate/next read
* `TLEN` - Observed template length (insert size)
* `SEQ` - Segment sequence
* `QUAL` - ASCII-encoded base qualities (Phred+33)

Optional fields follow in TAG:TYPE:VALUE format (e.g., NM:i:2 for edit distance).

## BED File Format

<https://genome.ucsc.edu/FAQ/FAQformat.html#format1>

Browser Extensible Data format for describing genomic intervals. Fields are **tab-separated**. The first three fields are required.

* `chrom` - Chromosome name (e.g., chr1, chrX)
* `chromStart` - 0-based start position of the feature
* `chromEnd` - End position of the feature (not included in the display)

Optional fields (must be included in order):

* `name` - Name of the BED line
* `score` - Score between 0 and 1000
* `strand` - '+' or '-'
* `thickStart` - Starting position where the feature is drawn thickly
* `thickEnd` - Ending position where the feature is drawn thickly
* `itemRgb` - RGB color value (e.g., 255,0,0)
* `blockCount` - Number of blocks (exons)
* `blockSizes` - Comma-separated list of block sizes
* `blockStarts` - Comma-separated list of block starts relative to chromStart

Note: BED uses 0-based, half-open coordinates (start is 0-based, end is exclusive).

## VCF File Format

<https://samtools.github.io/hts-specs/VCFv4.3.pdf>

Variant Call Format for storing gene sequence variations. Fields are **tab-separated**. The file contains meta-information lines (starting with '##'), a header line (starting with '#'), and data lines.

* `CHROM` - Chromosome name
* `POS` - 1-based position of the variant
* `ID` - Semicolon-separated list of unique identifiers (e.g., dbSNP rs number) or '.' if unknown
* `REF` - Reference base(s) at the position
* `ALT` - Comma-separated list of alternate alleles
* `QUAL` - Phred-scaled quality score for the assertion made in ALT
* `FILTER` - 'PASS' if the position has passed all filters, otherwise a semicolon-separated list of failed filters
* `INFO` - Semicolon-separated list of additional information (key=value pairs)
* `FORMAT` - Colon-separated list of genotype fields (only present if genotype data exists)
* `SAMPLE(s)` - Genotype information for each sample, corresponding to FORMAT fields
