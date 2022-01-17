# EDIR -- Exome Database of Interspersed Repeats


EDIR is dataset for interspersed repeats in coding regions of the human genome. These are defined by two short, identical, repeated sequences (7-20 bp) separated by a maximum of 1000 bp from one another, and filtered by those with at least one repeat within an exon, or both in introns flanking an exon. The database can be queried using the function and files in the `Query database` directory.


---
## Query database
#### Download and usage 
* Download the contents of "Query database" directory (1.54 GB). Make sure all files are in the **same directory**!
* `repeat_USE.R` -- contains init for querying the database and usage example
   1. Set the `path` variable (line 25) with the directory
   2. Run `source()` (line 28)
   3. See bottom of script for usage examples
 
---

### Gene Lookup Function
> gene_lookup (gene, length, mindist, maxdist, summary = F)

#### Description
   Queries EDIR database for a gene of interest. Results can be refined by defining the optional parameters of repeat sequence length, minimum/maximum distance.

#### Parameters
   | Parameters   | Type   | Description                                                         | Default |
   | ------------ |:------:| :-------------------------------------------------------------------|:-------:|
   | gene         | str    | Gene name to look up (case insensitive)                             | NA      |
   | length       | int    | Length of repeat sequence. If not defined, will output all results. | NA      |
   | mindist      | int    | Minimum distance between repeat sequences.                          | 0       |
   | maxdist      | int    | Maximum distance between repeat sequences.                          | 1000    |
   | summary      | bool   | Whether to return output summary.                                   | FALSE   |

#### Function returns
If `summary` is FALSE
   * Results summary printed to console with an overview of the results. Total runtime is also included.
   * The function returns the detailed results of the query.

If `summary` is TRUE
   * Results summary printed to console with an overview of the results. Total runtime is also included.
   * The function returns a tibble with:
      * `$summary` identical to the overview printed in the console.
      * `$results` containing the detailed results of the query.

---
## Generate database
This directory contains the main scripts and input used to generate the EDIR database. These are divided into 5 steps in the pipeline:
1. Input setup
   * Ensembl reference *Homo sapiens* genome (GRCh38) from [Ensembl](ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/)
   * [UCSC Known canonical transcripts](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/)
   * Gene and transcript information from [BiomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
3. Scan genome for 7bp sequences
5. Annotate 7bp output with gene info
6. Expand to 8-20bp and annotate
7. Filter data

---
## Contact
DT Laura Vo Ngoc doan.vongoc@vub.be
