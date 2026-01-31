## Table of Contents

- [Bioinformatics file formats](#bioinformatics-file-formats)
  - [Coordinate Systems](#coordinate-systems)
  - [GFF/GTF File Format](#gffgtf-file-format)
  - [FASTA File Format](#fasta-file-format)
  - [FASTQ File Format](#fastq-file-format)
  - [SAM File Format](#sam-file-format)
  - [BED File Format](#bed-file-format)
  - [VCF File Format](#vcf-file-format)
  - [AnnData (h5ad) File Format](#anndata-h5ad-file-format)

# Bioinformatics file formats

![https://xkcd.com/927/](assets/standards.png)

I lose track of all the bioinformatics file formats, so finally I'll keep track of them here.

## Coordinate Systems

Different bioinformatics file formats use different coordinate systems, which is a common source of off-by-one errors. Understanding these systems is essential when converting between formats.

### 0-based vs 1-based Coordinates

* **1-based**: The first base of a sequence is position 1. This is the "natural" counting system used by biologists.
* **0-based**: The first base of a sequence is position 0. This is common in programming languages and some file formats.

### Closed vs Half-open Intervals

* **Closed interval**: Both the start and end positions are included in the interval.
* **Half-open interval**: The start position is included, but the end position is excluded.

### Mathematical Interval Notation

Intervals can be expressed using bracket notation where square brackets $[\ ]$ indicate inclusion and parentheses $(\ )$ indicate exclusion:

| Notation | Name | Description |
|----------|------|-------------|
| $[a, b]$ | Closed | Both endpoints included: $a \leq x \leq b$ |
| $[a, b)$ | Half-open (left-closed, right-open) | Start included, end excluded: $a \leq x < b$ |
| $(a, b]$ | Half-open (left-open, right-closed) | Start excluded, end included: $a < x \leq b$ |
| $(a, b)$ | Open | Both endpoints excluded: $a < x < b$ |

### Coordinate Systems by File Format

| Format | Base | Interval | Notation | Example for bases 1-3 |
|--------|------|----------|----------|----------------------|
| GFF/GTF | 1-based | Closed | $[\text{start}, \text{end}]$ | start=1, end=3 |
| SAM | 1-based | Closed | $[\text{start}, \text{end}]$ | POS=1 (+ CIGAR 3M) |
| VCF | 1-based | Closed | $[\text{start}, \text{end}]$ | POS=1 |
| BED | 0-based | Half-open | $[\text{start}, \text{end})$ | start=0, end=3 |
| BAM (internal) | 0-based | Half-open | $[\text{start}, \text{end})$ | pos=0, end=3 |
| UCSC database | 0-based | Half-open | $[\text{start}, \text{end})$ | chromStart=0, chromEnd=3 |

### Conversion Example

Consider a 3-nucleotide feature (e.g., `ATG`) at the very beginning of a chromosome:

```
Sequence:    A  T  G  C  A  A  ...
1-based:     1  2  3  4  5  6  ...
0-based:     0  1  2  3  4  5  ...
```

To represent the `ATG` region:

* **GFF/GTF** (1-based, closed): `start=1, end=3`
* **BED** (0-based, half-open): `start=0, end=3`
* **SAM** (1-based): `POS=1` with CIGAR `3M`

The length of an interval can be calculated as:
* Closed interval: $\text{length} = \text{end} - \text{start} + 1$
* Half-open interval: $\text{length} = \text{end} - \text{start}$

This is why half-open coordinates are convenient for programming—the length calculation is simpler, and adjacent features share the same boundary coordinate without overlap.

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

## AnnData (h5ad) File Format

<https://anndata.readthedocs.io/en/latest/fileformat-prose.html>

AnnData (Annotated Data) is a Python data structure from the `anndata` package, commonly used with scanpy for single-cell genomics analysis. AnnData objects are saved to disk as `.h5ad` files, which are based on the HDF5 format.

### Core Components

* `X` - The primary data matrix (cells × genes). Can be dense (numpy array) or sparse (scipy sparse matrix). Typically contains expression counts or normalized values.
* `obs` - Observation (cell) metadata as a pandas DataFrame. Rows correspond to cells, columns are annotations (e.g., cell type, batch, cluster labels). Index contains cell barcodes/IDs.
* `var` - Variable (gene/feature) metadata as a pandas DataFrame. Rows correspond to genes, columns are annotations (e.g., gene symbols, highly variable flags). Index contains gene IDs.
* `uns` - Unstructured annotation as a dictionary. Stores miscellaneous data that doesn't fit elsewhere (e.g., color palettes, analysis parameters, UMAP settings).

### Additional Layers and Embeddings

* `layers` - Alternative matrix representations with the same dimensions as `X`. Common uses include storing raw counts in `layers['counts']` while `X` holds normalized data.
* `obsm` - Multi-dimensional observation annotations (cells × dimensions). Stores embeddings like PCA (`obsm['X_pca']`), UMAP (`obsm['X_umap']`), or t-SNE coordinates.
* `varm` - Multi-dimensional variable annotations (genes × dimensions). Stores gene-level embeddings like PCA loadings (`varm['PCs']`).
* `obsp` - Pairwise observation annotations (cells × cells). Stores cell-cell relationships like neighborhood graphs (`obsp['connectivities']`, `obsp['distances']`).
* `varp` - Pairwise variable annotations (genes × genes). Stores gene-gene relationships.

### File Structure

The h5ad format stores AnnData as an HDF5 file with the following structure:

```
file.h5ad
|-- X                    # Main data matrix
|-- obs                  # Cell metadata (DataFrame stored as HDF5 group)
|-- var                  # Gene metadata (DataFrame stored as HDF5 group)
|-- uns/                 # Unstructured data (nested dictionaries/arrays)
|-- obsm/                # Cell embeddings (e.g., X_pca, X_umap)
|-- varm/                # Gene embeddings
|-- obsp/                # Cell-cell graphs
|-- varp/                # Gene-gene relationships
+-- layers/              # Alternative matrices
```

### Reading and Writing

```python
import anndata as ad
import scanpy as sc

# Read h5ad file
adata = ad.read_h5ad('file.h5ad')
# or
adata = sc.read_h5ad('file.h5ad')

# Write h5ad file
adata.write_h5ad('output.h5ad')

# Write with compression
adata.write_h5ad('output.h5ad', compression='gzip')
```

### Sparse Matrix Support

The `X` matrix and layers often use sparse formats to efficiently store single-cell data (which is typically >90% zeros). Common formats include CSR (Compressed Sparse Row) and CSC (Compressed Sparse Column).
