# WAAFLE

**WAAFLE** (a **W**orkflow to **A**nnotate **A**ssemblies and **F**ind **L**GT Events) is a method for finding novel LGT (Lateral Gene Transfer) events in assembled metagenomes, including those from human microbiomes. Please direct questions to the [WAAFLE "topic" in the bioBakery Support Forum](https://forum.biobakery.org/c/Microbial-community-profiling/WAAFLE).

_Large data sets described in the WAAFLE publication, including synthetic datasets and WAAFLE profiles of HMP1-II and HMP2 metagenomes, are [available for download here](https://huttenhower.sph.harvard.edu/waafle). Additional derived data will be made available as supporting information alongside the WAAFLE publication itself (see citation below, and check back for a final version)._

## Citation

The WAAFLE manuscript has been accepted!

> Tiffany Y. Hsu*, Etienne Nzabarushimana*, Dennis Wong, Chengwei Luo, Robert G. Beiko, Morgan Langille, Curtis Huttenhower, Long H. Nguyen**, Eric A. Franzosa**. _Profiling lateral gene transfer events in the human microbiome using WAAFLE_. Nature Microbiology. (In press.)
> 
> \* = co-lead; \*\* = co-supervised


In the meantime, if you use WAAFLE in your work, please cite the WAAFLE repository on GitHub: https://github.com/biobakery/waafle.

## Software requirements

* Python 3+ or 2.7+
* Python `numpy` (tested with v1.13.3)
* [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (tested with v2.6.0)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (for performing contig-level QC; tested with v2.2.3)

## Installation

### From Anaconda

WAAFLE is available as a `conda` recipe from our `conda` channel, `biobakery`. Installing WAAFLE in this way will automatically satisfy the requirements listed above.

If you haven't already, set channel priorities as follows:

```
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
$ conda config --add channels biobakery
```

Then install WAAFLE with:

```
$ conda install waafle -c biobakery
```

### From PyPI

Install the WAAFLE software and Python dependencies with `pip`:

```
$ pip install waafle
```

This will _not_ satisfy the non-Python dependencies listed above (e.g. BLAST).

### From Source

You can also clone the WAAFLE package from GitHub:

```
$ git clone https://github.com/biobakery/waafle.git
```

## Database requirements

WAAFLE requires two input databases: 1) a BLAST-formatted **nucleotide sequence database** and 2) a corresponding **taxonomy file**. The versions used in the WAAFLE publication, which are based on the ChocoPhlAn 2.0 pangenome catalog, are available for download here:

* [chocophlan2.tar.gz](http://huttenhower.sph.harvard.edu/waafle_data/chocophlan2.tar.gz) (4.3 GB)
* [chocophlan_taxonomy.tsv](http://huttenhower.sph.harvard.edu/waafle_data/chocophlan2_taxonomy.tsv) (<1 MB)

The BLAST database must be unpacked before it can be used. You can do this from the command line with:

 ```
 $ tar xzfv chocophlan2.tar.gz
 ```
 
 Which will create a folder of BLAST database files that can be referenced as `chocophlan2/chocophlan2`.

See [Advanced Topics](#advanced-topics) at the bottom of this document for guidance on creating your own WAAFLE-compatible databases.

## Input data

An individual WAAFLE run requires one or two non-fixed inputs: 1) a file containing **metagenomic contigs** and (optionally) 2) a GFF file describing **gene/ORF coordinates** along those contigs.

### Input contigs (required)

Contigs should be provided as nucleotide sequences in FASTA format. Contigs are expected to have unique, BLAST-compatible headers. WAAFLE is optimized for working with fragmentary contigs from partially assembled metagenomes (spanning 2-20 genes, or roughly 1-20 kb). WAAFLE is not optimized to work with extremely long contigs (100s of kbs), scaffolds, or closed genomes. The WAAFLE developers recommend [MEGAHIT](https://github.com/voutcn/megahit) as a general-purpose metagenomic assembler.

* [A sample contigs input file](https://raw.githubusercontent.com/biobakery/waafle/refs/heads/main/demo/input/demo_contigs.fna)

### Input ORF calls (optional)

The optional GFF file, if provided, should conform to the [GFF format]([https://useast.ensembl.org/info/website/upload/gff.html). The WAAFLE developers recommend [Prodigal](https://github.com/hyattpd/Prodigal) as a general-purpose ORF caller with GFF output. WAAFLE can alternatively generate a GFF file as part of its normal operation ([as described below](#step-2-gene-calling-with-waafle_genecaller)).

* [A sample GFF file produced by WAAFLE](https://github.com/biobakery/waafle/blob/master/demo/output/demo_contigs.gff)
* [A sample GFF file produced by Prodigal](https://github.com/biobakery/waafle/blob/master/demo/output_prodigal/demo_contigs.prodigal.gff)

## Performing a WAAFLE analysis

Analyzing a set of contigs with WAAFLE requires performing four steps: 1) subjecting the contigs to homology-based search, 2) identifying genes / open reading frames (ORFs) along the contigs, 3) combining the results of Steps 1 and 2 to identify candidate LGT events, and 4) performing contig-level quality control (QC) to weed out misassembled contigs. All of these steps can be performed using WAAFLE utilities. Steps 2 and 4 can optionally be performed outside of WAAFLE using other methods.

The rate-limiting step in the WAAFLE workflow (by far) is the upstream `blastn` search of the user's contigs against the ChocoPhlAn 2 database, which varies closely with the total input size (as measured by number of contigs or total assembly size). Using a typical laptop-sized computing instance, we observed a rate of ~1 CPU minute per megabase of input sequence for the upstream search and ~1 CPU second per megabase of input for both the internal gene-calling and LGT-scoring steps.

### Step 1: Homology-based search with  `waafle_search`

`waafle_search` is a light wrapper around `blastn` to help guide the nucleotide-level search of your metagenomic contigs against a WAAFLE-formatted database (for example, it ensures that all of the non-default BLAST output fields required for downstream processing are generated).

A sample call to `waafle_search` with input contigs `contigs.fna` and a blast database located in `chocophlan2` would be:

```
$ waafle_search contigs.fna chocophlan2/chocophlan2
```

By default, this produces an output file `contigs.blastout` in the same location as the input contigs. See the `--help` menu for additional configuration options.

### Step 2: Gene calling with `waafle_genecaller`

WAAFLE can identify gene coordinates of interest directly from the BLAST output produced in the previous step:

`$ waafle_genecaller contigs.blastout`

This produces a file in GFF format called  `contigs.gff` for use in the next and last WAAFLE step. See the `--help` menu for additional configuration options.

[As described above](#input-orf-calls-optional), a user can alternatively provide a GFF file for their contigs that was generated outside of WAAFLE using another method, e.g. [Prodigal](https://github.com/hyattpd/Prodigal).

### Step 3: Identify candidate LGT events with `waafle_orgscorer`

The next and most critical step of a WAAFLE analysis is combining the BLAST output generated in Step 1 with a GFF file (generated in Step 2 or with an external method) to 1) taxonomically score genes along the length of the input contigs and then 2) identify those contigs as derived from a single clade or a pair of clades (i.e. putative LGT). Assuming you have run Steps 1 and 2 as described above, a sample call to `waafle_orgscorer` would be:

```
$ waafle_orgscorer \
  contigs.fna \
  contigs.blastout \
  contigs.gff \
  chocophlan2_taxonomy.tsv
```

This will produce three output files which divide and describe your contigs as putative LGT contigs, single-clade (no-LGT) contigs, and unclassified contigs (e.g. those containing no genes):

* `contigs.lgt.tsv`
* `contigs.no_lgt.tsv`
* `contigs.unclassified.tsv`

These files and their formats are described in more detailed below (see "WAAFLE outputs").

`waafle_orgscorer` offers many options for fine-tuning your analysis. The various analysis parameters have been pre-optimized for maximum specificity on both short contigs (containing as little as two partial genes) and longer contigs (10s of genes). These options are detailed in the `--help` menu:

### Step 4: Filter out misassembled contigs

While WAAFLE's method of LGT detection has been optimized to distinguish LGT from other biological events (e.g. gene deletion), it does so assuming that the input contigs are biologically valid. This is notable as misassembled contigs, in particular chimeric assembly of genetic material from 2+ organisms, can spuriously manifest as LGT. It is therefore critical to filter out low-quality contigs before performing further downstream analyses on WAAFLE outputs.

WAAFLE provides a set of utilities that can aid in the process of performing contig-level QC / filtering using the outputs of Steps 1-3. This method is [described below](#contig-level-quality-control). Alternatively, a user can perform contig-level QC (and remove/correct erroneous contigs) before starting their WAAFLE analysis using an external method, e.g. [metaMIC](https://github.com/ZhaoXM-Lab/metaMIC).

Notably, all contigs included in the demo files linked from this manual passed the QC checks imposed by the WAAFLE QC utilities referenced above.

## WAAFLE outputs

The `contigs.lgt.tsv` output file lists the details of putative LGT contigs. Its fields are a superset of the types of fields included in the other output files. The following represents the first two lines/rows of a `contigs.lgt.tsv` file *transposed* such that first line (column headers) is shown as the first column and the details of the first LGT contig (second row) are shown as the second column:

```
CONTIG_NAME          12571
CALL                 lgt
CONTIG_LENGTH        9250
MIN_MAX_SCORE        0.856
AVG_MAX_SCORE        0.965
SYNTENY              AABAAAA
DIRECTION            B>A
CLADE_A              s__Ruminococcus_bromii
CLADE_B              s__Faecalibacterium_prausnitzii
LCA                  f__Ruminococcaceae
MELDED_A             --
MELDED_B             --
TAXONOMY_A           r__Root|k__Bacteria|p__Firmicutes|...|g__Ruminococcus|s__Ruminococcus_bromii
TAXONOMY_B           r__Root|k__Bacteria|p__Firmicutes|...|g__Faecalibacterium|s__Faecalibacterium_prausnitzii
LOCI                 252:668:-|792:1367:-|1557:2360:-|2724:3596:-|4540:5592:+|5608:7977:+|8180:8425:+
ANNOTATIONS:UNIPROT  R5E4K6|D4L7I2|D4JXM0|D4L7I1|D4L7I0|None|D4L7H8
```

The fields in detail:

* **`CONTIG_NAME`**: the name of the contig from the input FASTA file.
* **`CALL`**: indicates that this was an LGT contig.
* **`CONTIG_LENGTH`**: the length of the contig in nucleotides.
* **`MIN_MAX_SCORE`**: the minimum score for the pair of clades explaining the contig along the length of the contig. (i.e. the score for identifying this as a putative LGT contig, with a default threshold of 0.8.)
* **`AVG_MAX_SCORE`**: the average score for the pair of clades explaining the contig (used for ranking multiple potential explanations of the contig).
* **`SYNTENY`**: the pattern of genes assigned to the A or B clades along the contig. `*` indicates that the gene could be contributed by either clade; `~` indicates an ignored gene; `!` indicates a problem (should not occur).
* **`DIRECTION`**: indicates this as a putative B-to-A transfer event, as determined from synteny (A genes flank the inserted B gene). `A?B` indicates that the direction was uncertain.
* **`CLADE_A`**: the name of clade A.
* **`CLADE_B`**: the name of clade B.
* **`LCA`**: the lowest common ancestor of clades A and B. A higher-level LCA indicates a more remote LGT event.
* **`MELDED_A`**: if using a meld reporting option, the individual melded clades will be listed here. For example, if a contig could be explained by a transfer from *Genus_1 species_1.1* to either *Genus_2 species_2.1* or *Genus_2 species_2.2*, this field would list `species_2.1; species 2.2` and *Genus 2* would be given as `CLADE_A`.
* **`MELDED_B`**: *see above*.
* **`TAXONOMY_A`**: the full taxonomy of `CLADE_A`.
* **`TAXONOMY_B`**: the full taxonomy of `CLADE_B`.
* **`LOCI`**: Ccordinates of the loci (genes) that were considered for this contig in format `START:STOP:STRAND`.
* **`ANNOTATIONS:UNIPROT`**: indicates that UniProt annotations were provided for the genes in the input sequence database (in format `UNIPROT=IDENTIFIER`). The best-scoring UniProt annotation for each gene is given here. (Additional annotations would appear as additional, similarly-formatted columns in the output.)

## Contig-level quality control

WAAFLE bundles two utilities, `waafle_junctions` and `waafle_qc`, that can be used to aid in the identification and removal of low-quality contigs.

### Quantifying junction support with `waafle_junctions`

`waafle_junctions` maps reads to contigs to quantify support for individiual gene-gene junctions. Specifically, two forms of support are considered/computed:

1. Sequencing fragments (paired reads) that span the gene-gene junction. These are less common for junctions that are larger than typical sequencing insert sizes (~300 nts).

2. Relative junction coverage compared to the mean coverage of the two flanking genes. If the junction is not well covered relative to its flanking genes, it may represent a non-biological overlap.

A sample call to `waafle_junctions` looks like:

```
$ waafle_junctions \
  contigs.fna \
  contigs.gff \
  --reads1 contigs_reads.1.fq \
  --reads2 contigs_reads.2.fq \
```

With this call, `waafle_junctions` will use [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to index the contigs and then align the input reads (pairwise) against the index to produce a SAM file. (`waafle_junctions` can also interpret a mapping from an existing SAM file.) The alignment results are then interpreted to score individual junctions, producing an output file for each.

* `contigs.junctions.tsv`

A sample report for an individual junction looks like:

```
CONTIG             SRS011086_k119_10006
GENE1              1:1363:-
GENE2              1451:2382:-
LEN_GENE1          1363
LEN_GENE2          932
GAP                87
JUNCTION_HITS      5
COVERAGE_GENE1     8.1079
COVERAGE_GENE2     11.3573
COVERAGE_JUNCTION  10.9101
RATIO              1.1210
```

This report indicates that the junction between genes 1 and 2 (which may or may not be an LGT junction) was well supported: it was spanned by 5 mate-pairs (`JUNCTION_HITS=5`) and had coverable coverage (`RATIO=1.12`) to the mean of its flanking genes (8.11 and 11.4).

`waafle_junctions` can be tuned to produce additional gene- and nucleotide-level quality reports. Consult the `--help` menu for a full list of options.

### Using junction data for contig QC with `waafle_qc`

The `waafle_qc` script interprets the output of `waafle_junctions` to remove contigs with weak read support at one or more junctions. Currently, the script focuses on the junctions flanking LGT'ed genes among putative LGT-containing contigs.

A sample call to `waafle_qc` looks like:

```
$ waafle_qc contigs.lgt.tsv contigs.junctions.tsv
```

Where `contigs.junctions.tsv` is the output of `waafle_junctions` on this set of contigs and its underlying reads. This produces a file `contigs.lgt.tsv.qc_pass`: a subset of the original LGT calls that were supported by read-level evidence.

By default, a junction is supported if it was contained in 2+ mate-pairs *or* had >0.5x the average coverage of its two flanking genes. These thresholds are tunable with the `--min-junction-hits` and `--min-junction-ratio` parameters of `waafle_qc`, respectively. Consult the `--help` menu for a full list of options.

## Advanced topics

### Formatting a sequence database for WAAFLE

WAAFLE performs a nucleotide-level search of metagenomic contigs against a collection of taxonomically annotated protein-coding genes (*not* complete genomes). A common way to build such a database is to combine a collection of microbial pangenomes of interest. The protein-coding genes should be organized in FASTA format and then indexed for use with `blastn`. For example, the FASTA file `waafledb.fnt` would be indexed as:

```
$ makeblastdb -in waafledb.fnt -dbtype nucl
```

WAAFLE expects the input FASTA sequence headers to follow a specific format. At a minimum, the headers must contain a unique sequence ID (`GENE_123` below) followed by `|` (pipe) followed by a taxon name or taxonomic ID (`SPECIES_456` below):

```
>GENE_123|SPECIES_456
```

Headers can contain additional `|`-delimited fields corresponding to functional annotations of the gene. These fields have the format `SYSTEM=IDENTIFIER` and multiple such fields can be included per header, as in:

```
>GENE_123|SPECIES_456|PFAM=P00789|EC=1.2.3.4
```
Headers are allowed to contain different sets of functional annotation fields. WAAFLE currently expects at most one annotation per annotation system per gene; this will be improved in future versions. (We currently recommend annotating genes against the [UniRef90 and UniRef50](http://www.uniprot.org/help/uniref) databases to enable link-outs to more detailed functional annotations in downstream analyses.)

WAAFLE assumes that the taxa listed in the sequence database file are all at the same taxonomic level (for example, all genera or all species or all strains).

### Formatting a WAAFLE taxonomy file

WAAFLE requires a taxonomy file to understand the taxonomic relationships among the taxa whose genes are included in the sequence database. The taxonomy file is a tab delimited list of child-parent relationships, for example:

```
k__Animalia      r__Root
p__Chordata      k__Animalia
c__Mammalia      p__Chordata
o__Primates      c__Mammalia
f__Hominidae     o__Primates
g__Homo          f__hominidae
s__Homo sapiens  g__Homo
```

While the format of this file is relatively simple, it has a number of critical structural constraints that must be obeyed:

* All taxon names/IDs used in the sequence database must appear is the taxonomy file.
 
* The file must contain a root taxon from which all other taxa descend (the root taxon should be named `r__Root`, as above).

* All taxon names/IDs used in the sequence database must be the same distance from the root.

The following properties of the taxonomy file are optional:

* Additional taxa *below* the level of the taxa in the sequence file can be included in the taxonomy file. For example, a species-level sequence database could contain isolate genomes as children of the species-level clades in the taxonomy file. (WAAFLE can use such entries as "leaf support" for LGT events.)

* We recommend prefixing taxonomic clades to indicate their level. For example, `g__Homo` identifies *Homo* as a genus above.

## Contributions ##
This work was supported by the National Institutes of Health grants U54DE023798 (CH), R24DK110499 (CH), and K23DK125838 (LHN), the American Gastroenterological Association Research Foundation’s Research Scholars Award (LHN), and the Crohn’s and Colitis Foundation Career Development Award (LHN). The content is solely the responsibility of the authors and does not necessarily represent the official views of the NIH. We thank April Pawluk, Kelsey N. Thompson, Kristen Curry, and Todd Treangen for comments on the manuscript and helpful discussions. We also acknowledge Monia Michaud, Casey Dulong, and Yan Yan for their help with validation experiments. Computational work was conducted on the FASRC Cannon cluster supported by the FAS Division of Science Research Computing Group at Harvard University.
