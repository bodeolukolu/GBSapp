
<p align="right">
<img src="https://github.com/bodeolukolu/GBSapp/blob/master/misc/GBSapp_logo.PNG" width="173" height="111">
</p>

# Introduction
GBSapp v2.2 is an automated pipeline for variant calling and filtering. The pipeline intuitively integrates existing/novel best practices, some of which can be controlled by user-defined parameters. It optimizes memory and speed at various steps of the pipeline, for example, a novel approach performs compression and decompression of unique reads before and after read alignment, respectively. Summary reports and visualizations allow for QC at each step of the pipeline.

For questions, bugs, and suggestions, please contact bolukolu@utk.edu.

## Features
- Easy use and designed for biologist.
- Dosage-based variant calling and filtering.
- Fully-automated: “walk-away” and “walk-through” mode.
- Allows for use of haplomes (haplotype-resolved), subgenomes (haploid), and pan-genomes (haploid or haplotype-resolved) references.
- Minimizes excess heterozygosity and allele dropout.
- Variant calling implemented for up to ploidy of 8.
- Input data: shotgun WGS, reduced representation sequence (e.g., OmeSeq-qRRS, GBS, ddRADseq), RNAseq and multiplexed-PCR data.
- Can subsample shotgun whole genome data for variant calling, i.e. in silico reduced representation sequencing (RRS).
- Parallelization of job on multiple compute cluster nodes (spark cluster infrastructure not required).
- Splice-aware aligner (STAR) allows for RNAseq data as input (recommended only for haploid or diploid genomes).
- Generates variant sequence context (useful for applications such as oligo/primer design & sequenced-based phylogenetic analysis).
- Variant calling delineates SNP from uniquely mapped and paralogs.
- Increase speed of variant calling based on dynamic downsampling (avoids allele dropout due to biased downsampling).
- Fast alignment due to joint-alignment method.
- Visualizations for report and QC.

- Functions under-development:
  - calling microhaplotypes
  - estimating ploidy level and aneuploidy



# Contents
- [Installation](#installation)
- [Usage](#usage)
  - [Basic usage](#Basic_usage)
  - [Project directory setup](#Project_directory_setup)
  - [Overview of workflow](#Overview_of_workflow)
  - [Configuration](#Configuration)
- [Related Software](#Related_Software)
- [Select Article Referencing GBSapp](#Select_Article_Referencing_GBSapp)
- [Acknowledgment](#Acknowledgment)
- [Troubleshooting](#Troubleshooting)
- [Versioning](#Versioning)
- [License](#License)

## Installation
- Currently, GBSapp is only available for unix-based systems.
- Clone or download the Git repository to your desired folder.
```bash
git clone https://github.com/bodeolukolu/GBSapp.git
```
- Installation occurs automatically the first time you run the pipeline.
- To install dependencies without running a job:  
```bash
GBSapp_dir/GBSapp install
```
- With the exception of R, and Python (v2.6 or greater) all dependencies are installed to a local directory within GBSapp.


**Dependencies:**<br />
```
Installed on first run of pipeline:
-----------------------------------
NextGenMap (ngm), minnimap2, STAR, bwa, samtools, picard, bcftools, GATK, java, R-ggplot2, CMplot, and R-AGHmatrix


Pre-install before running GBSapp:
----------------------------------
- R
- python v2.6 or greater
```


## Usage
### Basic Usage
The project directory should contain the following files and directories:
- **config file**: specifies run parameters (for details: GBSapp_vx.x/examples/config).
- **samples directory**: contains “se” and/or “pe” (“paired” and “single” name format acceptable) fastq file(s). Paired-end (pe) sample filenames might require formatting so that they end in “_R1.fastq” or  “.R1.fastq” and “_R2.fastq” or “.R2.fastq” (for details: GBSapp_vx.x/misc/format_fastq_filenames.txt).
- **refgenomes directory**: contains fasta file(s) of the reference genome.<br />* <u>Genomes with subgenome assemblies in single fasta file</u>: such as allopolyploids and segmental allopolyploids might require formatting to split fasta file into multiple file containing each subgenome (for details: GBSapp_vx.x/misc/split_subgenomes_format_fasta_headers.txt).<br />* <u>Haploids, diploids and autopolyploids with single reference genomes</u>: splitting fasta file not required.<br />* <u>Note</u>:formatting of fasta headers to contain minimal text (e.g. >Chr05) might be required (for details: GBSapp_v0.1/misc/format_fasta_headers.txt)
- **for help**: --help or -h
- **for version**: --version or -v

From command line, run GBSapp with options shown below (absolute or relative path)
```
$ bash	<path-to-GBSapp-directory/GBSapp>	<path-to-project-directory>
```

### Project directory setup
A project directory should contain the following sub-directories:
- **samples folder**: this contains your quality filtered sequence data.
- **refgenomes folder**: this contains your reference genome fasta file.
- **config.sh file**: a template of the configuration file is provided in the examples folder of the GBSapp download.

<img src="https://github.com/bodeolukolu/GBSapp/blob/master/misc/project_dir_setup.PNG" width="782" height="308">


### Overview of workflow
The figure below outlines the order of steps in the GBSapp pipeline
- In Progress

### Configuration
Using a text editor, save a file containing any of the following variables as 'config.sh' file and include it in your project directory.

**General parameters**

|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|threads|na|number of cores/processors|integer|Optional|
|walkaway|true|run in walk-away or walk-through mode|true or false|Optional|
|cluster|false|run on compute cluster node (default: slurm) or workstation|true or false|Optional|
|nodes|1|number of nodes|integer|Optional|
|gap_split_align|true|minimap2 if true (large gaps and SVs) and NGM if false(small gaps)|true or false|Optional|
|RNA|false|RNA-seq reads as input (STAR aligner)|true or false|Optional|
|variant_caller|gatk|gatk (recommended) or bcftools(diploid only, might be better for some genomes)|string|Optional|
|samples_alt_dir|false|links samples in separate directory to project directory|true or false|Optional|
|lib_type|RRS|RRS (reduced representation sequence e.g. GBS, ddRADseq, qRRS) or WGS (shotgun whole genome sequence)|string|Optional|
|subsample_WGS_in_silico_qRRS|false|Fast alternative to variant calling on whole genome data|true or false|Optional|



**Variant calling parameters**

|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|ploidy|na|value = 1,2,4,6, or 8|integer|Required|
|haplome_number|1|specify single value or range (comma delimited) up to maximum haplome i.e. ploidy level|integer|Optional|
|ref1|na|reference subgenome as .fasta file. Anchor-genome when other pangenomes/subgenomes are specified |integer|Optional|
|ref2|na|2nd reference genome as .fasta file|integer|Optional|
|ref3|na|3rd reference genome as .fasta file|integer|Optional|
|ploidy_ref1|na|ploidy-level for subgenome 1|integer|Optional|
|ploidy_ref2|na|ploidy-level for subgenome 2, only specify for subgenome specific variants|integer|Optional|
|ploidy_ref3|na|ploidy-level for subgenome 3, only specify for subgenome specific variants|integer|Optional|
|hap_ref|na|haplotype-resolved reference genome (# of haplomes typically = ploidy level) |integer|Optional|
|Get_Chromosome|na|variant calling on specific chromosomes, scaffolds,and contigs|comma delimited string(s)|optional|
|Exclude_Chromosome|na|variant calling to exclude specific chromosomes, scaffolds,and contigs|comma delimited string(s)|optional|

**note: haploid assemblies of pangenomes and subgenomes should be in individual fasta files (up to 3 fasta files)*
**note: short prefix for pangenome/subgenome pseudomolecules should be unique (i.e. >TF_Chr01 and >TL_Chr01 fasta sequence header for Ipomoea trifida and I. triloba, respectively)*
**note: for polyploids with variable ploidy of subgenomes specify ploidy of subgenomes (e.g. hexapploid: ploidy_ref1=4 and ploidy_ref2=2)*
**note: reference haplotype-resolved assembly (hap_ref) should be a single fasta file and should be chromosome-level assembly*
**note: designate chromosomes of haplotype-resolved/subgenome assemblies with single character suffix (alphabets: A-Z and a-z) e.g. Chr01A*


**Variant filtering parameters**

|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|p1|na|maternal parent (specified only for biparental populations)|string|Optional|
|p2|na|paternal parent (specified only for biparental populations)|string|Optional|
|biallelic|false|filter to output only biallelic variants|string|Optional|
|genotype_missingness|1|maximum proportion of missing genotypes allowed per sample|comma delimited decimal number(s)|Optional|
|sample_missingness|1|maximum proportion of missing samples allowed per variant|comma delimited decimal number(s)|Optional|
|exclude_samples|na|sample IDs to be exclude from filtered variant data set|comma delimited string(s)|Optional|
|select_samples|na|limit variant filtering to samples IDs in file delimited by newline|filename|Optional|
|minRD_1x|2|minimum read depth threshold|integer|Optional|
|minRD_2x|6|minimum read depth threshold|integer|Optional|
|minRD_4x|25|minimum read depth threshold|integer|Optional|
|minRD_6x|45|minimum read depth threshold|integer|Optional|
|minRD_8x|100|minimum read depth threshold|integer|Optional|
|pseg|0.001|p-value threshold for chi-square test of segregation distortion|decimal number|Optional|
|maf|0.02|minor allele frequency threshold|decimal number|Optional|
|filtered_vcf|false|generate filtered vcf file|string|Optional|



**Advanced parameters**
|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|max_pseudoMol|1000|maximum # of pseudomolecules (scaffold/contig) before stitching into non-contiguous pseudo-chromosomes|integer|Optional|
|uniquely_mapped|true|include uniquely mapped for variant calling |string|Optional|
|paralogs|true|include paralogs for variant calling |string|Optional|
|minmapq|20|minimum mapping quality|integer|Optional|
|downsample_2x|50|value for unbiased downsampling for 2x ploidy|integer|Optional|
|downsample_4x|100|value for unbiased downsampling for 4x ploidy|integer|Optional|
|downsample_6x|150|value for unbiased downsampling for 6x ploidy|integer|Optional|
|downsample_8x|200|value for unbiased downsampling for 8x ploidy|integer|Optional|
|maxHaplotype|128|maximum number of haplotypes per haploid genome across population(increase for polyploids/high heterozygosity/high background mutational load)|integer|Optional|
|use_softclip|false|use soft-clipped bases for variant calling|string|Optional|
|joint_calling|false| cohort calling will be performed if set to false|string|Optional|
|keep_gVCF|false|keep sample gVCF files, if additional samples will be included for future joint calling)|string|Optional|
|RE1|NA|sequence motif at start of R1 reads|string|Optional|
|RE2|NA|sequence motif at start of R2 reads|string|Optional|
|filter_ExcHet|false|test and filter for excess heterozygosity|string|Optional|
|genomecov_est|false|compute genome coverage for each sample|string|Optional|

**note: na indicates that variable is user-defined or hard-coded/computed intuitively, as well as a function of ploidy.*

Below is an example of a configuration file:

**config.sh**
```
# General_parameters
###################################################
threads=16
walkaway=true
cluster=true
nodes=1
gap_split_align=true
RNA=false
variant_caller=gatk
samples_alt_dir=false
lib_type=RRS
subsample_WGS_in_silico_qRRS=false

# Variant calling
###################################################
ploidy=6
haplome_number=6
# Variant calling with haploid subgenome(s)
# Anchored to ref1 for loci conserved across all subgenomes
ref1=TF.fasta
ref2=TL.fasta
ploidy_ref1=4
ploidy_ref2=2
# Variant calling with haplotype-resolved reference genome or pangenomes
# Anchored to haplome_ref1 for loci conserved across all haplomes
hap_ref=Ib.fasta
# exclue or limit variant calling to specific chromosomes
Get_Chromosome=TF_Chr01,TF_Chr02
Exclude_Chromosome=TF_Chr00,TL_Chr00
genomecov_est=false

# SNP-filtering:
####################################################
p1=M9
p2=M19
biallelic=false
genotype_missingness=1
sample_missingness=1
exclude_samples=S1,S2,S3
select_samples=pop.txt
minRD_2x=6
minRD_4x=25
minRD_6x=45
pseg=0.001
maf=0.05
filtered_vcf=true

# Advanced_parameters
###################################################
max_pseudoMol=1000
uniquely_mapped=true
paralogs=true
minmapq=20
downsample_2x=50
downsample_4x=100
downsample_6x=150
downsample_8x=200
variant_intervals=false
interval_list=variant_intervals.list
interval_list_ref1=variant_intervals_TF.list
interval_list_ref2=variant_intervals_TL.list
maxHaplotype=128
use_softclip=false
joint_calling=false
keep_gVCF=false
RE1=TGCAT
RE2=CATG
filter_ExcHet=false
```

Alternatively, a configuration file (outlined below) specifying only the ploidy level is sufficient to run GBSapp.

**config.sh**
```
### Variant calling
###################################################
ploidy=2
```
Since most of the parameters are hard-coded in an intuitive manner, by specifying only the ploidy levels, the pipelines determines the other parameters as stated below:<br />
- threads: computes available number of cores (n) and and uses n-2 threads
- defaults: refer to parameters above

## Related Software
- [ngsComposer: Empirical Base-call error-filtering and read preprocessing pipeline.](https://github.com/bodeolukolu/ngsComposer)
- [Qmatey: Quantitative Metagenomic Alignment and Taxonomic Exact-matching.](https://github.com/bodeolukolu/Qmatey)



## Select Article Referencing GBSapp
1. ngsComposer: an automated pipeline for empirically based NGS data quality filtering. [Kuster et al. 2021](https://doi.org/10.1093/bib/bbab092)
2. Genome-wide association study identified candidate genes controlling continuous storage root formation and bulking in hexaploid sweetpotato. [Bararyenya et al. 2020](https://doi.org/10.1186/s12870-019-2217-9)
3. Sequencing depth and genotype quality: accuracy and breeding operation considerations for genomic selection applications in autopolyploid crops [Gemenet et al. 2020](https://doi.org/10.1007/s00122-020-03673-2)
4. Genetic Diversity and Population Structure of the USDA Sweetpotato (Ipomoea batatas) Germplasm Collections Using GBSpoly [Wadl et al. 2019](https://doi:10.3389/fpls.2018.01166)

## Acknowledgment
This package has been developed as part of the [Genomic Tools for Sweetpotato Improvement project](https://sweetpotatogenomics.cals.ncsu.edu/) (GT4SP) and [SweetGAINS](https://cgspace.cgiar.org/handle/10568/106838), both funded by [Bill & Melinda Gates Foundation](https://www.gatesfoundation.org/).

## Troubleshooting
**Pre-Installation of R and Python:**<br />
```
    - To view R and python version (or to check if installed), from the terminal type:
          $ R --version
          $ python --version

    - For Ubuntu, install R and python using apt:
          $ sudo apt update && sudo apt upgrade
          $ sudo apt install r-base
          $ sudo apt install python3

    - For macOS, install using homebrew:
          brew install r
          brew install python3
```
**If GATK can't find python:**<br />
```
- Make sure python v2.6 or greater is installed and then type the command below in terminal
- $ sudo ln -sf /usr/bin/python2 /usr/bin/python
- or
- If you are using python3
- $ sudo apt update
- $ sudo apt install python-is-python3
```
**If samtools and bcftools doesn't install properly:**<br />
```
While the installation of samtools and bcftools are automated, the installation requires some dependencies:
  $ sudo apt-get update
  $ sudo apt-get install gcc
  $ sudo apt-get install make
  $ sudo apt-get install libbz2-dev
  $ sudo apt-get install zlib1g-dev
  $ sudo apt-get install libncurses5-dev
  $ sudo apt-get install libncursesw5-dev
  $ sudo apt-get install liblzma-dev
  $ sudo apt-get install libcurl4-gnutls-dev
  $ sudo apt-get install libssl-dev
```
**If NextGenMap doesn't install properly:**<br />
```
While the installation of NextGenMap (ngm) is automated, the installation requires some dependencies:
  $ sudo apt install cmake
```
**Problem with amount of memory and/or processors/cores specified:**<br />
```
- This might be due to specifing values greater than available resources
- Re-submit job with appropriate values or modify header of the GBSapp_run.sh batch file.
- If using compute cluster managers other than SLURM, header of GBSapp_run.sh batch can also be modified to fit the syntax of the cluster manager been used.
```
## Versioning
Versioning will follow major.minor.patch <a href="https://semver.org">semantic versioning format</a>.

## License
<a href="https://github.com/bodeolukolu/GBSapp/blob/master/LICENSE">Apache License Version 2.0</a>
