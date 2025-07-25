## Table of Contents

- [Bioinformatics file formats](#bioinformatics-file-formats)
  - [GFF/GTF File Format](#gffgtf-file-format)

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
