
<p align="right">
<img src="https://github.com/bodeolukolu/GBSapp/blob/master/misc/GBSapp_logo.PNG" width="173" height="111">
</p>

# Introduction
GBSapp (v. 0.2.4) is an automated pipeline for variant calling and filtering (including microhaplotypes or multi-SNP markers). The pipeline intuitively integrates existing and novel best practices, some of which can be controlled by user-defined parameters. It optimizes memory and speed at various points in the pipeline, for example, a novel approach performs compression and decompression of unique reads before and after pre-variant read alignment, respectively. Summary reports and visualizations allow for QC at each step of the pipeline.

For questions, bugs, and suggestions, please contact bolukolu@utk.edu.

## Features
- Fully-automated: “walk-away” and “walk-through” mode.
- Dosage-based variant/haplotype calling and filtering.
- Haploid (1x), Diploid (2x), Tetraploid (4x), Hexaploid (6x), and Octoploid (8x).
- Captures and codes variable dosage/copy number in paleopolyploids.
- Can restrict variants call to reads to uniquely mapped
- Multi-locus variants (paralogs-derived) can also be produced (annotates variant with reference genome copy number).
- Additional haplotype-based filtering (useful for high-fidelity targeted-genotyping).
- Generates microhaplotypes (multi-SNP) markers
- Produces sequence context of variants
- Dynamic unbiased down-sampling on locus-by-locus basis (preserving allelic ratios).
- Easy to learn, and use.

## Contents
- [Installation](#installation)
- [Usage](#usage)
  - [Basic usage](#Basic_usage)
  - [Overview of workflow](#Overview_of_workflow)
  - [Configuration](#Configuration)
- [Related Software](#Related_Software)
- [Select Article Referencing GBSapp](#Select_Article_Referencing_GBSapp)
- [Acknowledgment](#Acknowledgment)
- [Troubleshooting](#Troubleshooting)
- [Versioning](#Versioning)
- [License](#License)

## Installation
- Currently, GBSapp is only available for unix-based systems (i.e. macOS and Linux/UNIX).
- Clone or download the Git repository to your desired folder.
```bash
git clone https://github.com/bodeolukolu/GBSapp.git
```
- Installation occurs automatically the first time you run the pipeline.
- With the exception of R and Python (v2.6 or greater) all dependencies are installed to a local directory within GBSapp.


**Dependencies:**<br />
```
Installed on first run of pipeline:
-----------------------------------
bwa, picard, samtools, bcftools, GATK, python, R, R-ggplot2, and R-AGHmatrix


Pre-install before running GBSapp:
----------------------------------
- R
- python v2.6 or greater
```


## Usage
### Basic Usage
The project directory should contain the following files adn directories:
- **(1) config file**: specifies run parameters (for details: GBSapp_vx.x/examples/config).
- **(2) samples directory**: contains “se” and/or “pe” (“paired” and “single” name format acceptable) fastq file(s). Paired-end (pe) sample filenames might require formatting so that they end in “_R1.fastq” or  “.R1.fastq” and “_R2.fastq” or “.R2.fastq” (for details: GBSapp_vx.x/misc/format_fastq_filenames.txt).
- **(3)	refgenomes directory**: contains fasta file(s) of the reference genome.<br />* <u>Genomes with subgenome assemblies in single fasta file</u>: such as allopolyploids and segmental allopolyploids might require formatting to split fasta file into multiple file containing each subgenome (for details: GBSapp_vx.x/misc/split_subgenomes_format_fasta_headers.txt).<br />* <u>Haploids, diploids and autopolyploids with single reference genomes</u>: splitting fasta file not required.<br />* <u>Note</u>:formatting of fasta headers to contain minimal text (e.g. >Chr05) might be required (for details: GBSapp_v0.1/misc/format_fasta_headers.txt)

![alt text](https://github.com/bodeolukolu/GBSapp/blob/master/misc/project_dir_setup.PNG?raw=true)

From command line, run GBSapp as shown below (absolute or relative path)
```bash
$ bash	<path-to-GBSapp-directory/GBSapp.sh>	<path-to-project-directory>
```

### Overview of workflow
The figure below outlines the order of steps in the GBSapp pipeline
- In Progress

### configuration
Using a text editor, save a file containing any of the following variables as 'config' file (no file extension) and include it in your project directory.

**General parameters**

|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|threads|2|number of cores/processors|integer|Optional|
|walkaway|true|run in walk-away or walk-through mode|true or false|Optional|
|cluster|false|run on compute cluster node (default: slurm) or workstation|true or false|Optional|

**Variant calling parameters**

|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|ploidy|na|ploid level = 1,2,4,6, or 8 (ploidy>8 allowed for autopolyploid). Used to determine maximum number of haplotypes per sample|integer|Required|
|ref1|na|reference genome as .fasta file (haploids, diploids, autopolyploids, or 1st subgenome in allopolyploids)|integer|Optional|
|ref2|na|2nd reference genome as fasta file|integer|Optional|
|ref3|na|3rd reference genome as fasta file|integer|Optional|
|ref4|na|4th reference genome as fasta file|integer|Optional|
|ploidy_ref1|na|ploid level for subgenome 1|integer|Required|
|ploidy_ref2|na|ploid level for subgenome 1|integer|Required|
|ploidy_ref3|na|ploid level for subgenome 1|integer|Required|
|ploidy_ref4|na|ploid level for subgenome 1|integer|Required|
|multilocus|true|use reads that align to single or multiple loci|true or false|Optional|
|paleopolyploid|false|capture/code variable dosage/copy number (i.e. 2x,4x,6x, and 8x)|true or false|Optional|



**Note: na indicates that variable is user-defined or hard-coded/computed intuitively, as well as a function of ploidy.*

**Variant filtering parameters**

|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|p1|na|maternal parent (specified only for biparental populations)|string|Optional|
|p2|na|paternal parent (specified only for biparental populations)|string|Optional|
|genotype_missingness|0.2|maximum proportion of missing genotypes allowed per sample (multiple comma delimited thresholds permitted)|decimal number|Optional|
|sample_missingness|0.2|maximum proportion of missing samples allowed per variant (multiple comma delimited thresholds permitted)|decimal number|Optional|
|exclude_samples|na|sample IDs to be exclude from filtered variant data set|comma delimited strings|Optional|
|minRD_1x|2|minimum read depth threshold|integer|Optional|
|minRD_2x|6|minimum read depth threshold|integer|Optional|
|minRD_4x|25|minimum read depth threshold|integer|Optional|
|minRD_6x|45|minimum read depth threshold|integer|Optional|
|minRD_8x|100|minimum read depth threshold|integer|Optional|
|pseg|0.001|p-value threshold for chi-square test of segregation distortion|decimal number|Optional|
|maf|0.02|minor allele frequency threshold|decimal number|Optional|
|snpformats|false|variant data set with alleles formatted as base (A,C,G,T) and/or degenerate notation|true or false|Optional|


**Advanced parameters**
|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|haplome_number|1|number of haplome(s) represented by reference genome|integer|Optional|
|maxHaplotype|6| maximum number of haplotypes per haploid genome across population|integer|Optional|
|unbiased_downsample|25| maximum read number per alignment start position (per haploid genome). 0 = no downsample|integer|Optional|
|biased_downsample|0| maximum read number per alignment start position (per haploid genome). 0 = no downsample|integer|Optional|
|softclip|true| do not use soft Clipped bases (recommended)|string|Optional|



Below is an example of a configuration file:

**config**
```
# General parameters
###################################################
threads=24
walkaway=true
cluster=true

# Variant calling
###################################################
ploidy=6
ref1=TF.fasta
ref2=TL.fasta
ploidy_ref1=4
ploidy_ref2=2
paleopolyploid=false
multilocus=true

# SNP-filtering:
####################################################
p1=Beauregard
p2=Tanzania
genotype_missingness=0.1,0.2,0.3
sample_missingness=0.1,0.2,0.3
exclude_samples=S1,S2,S3
minRD_2x=6
minRD_4x=25
minRD_6x=45
pseg=0.001
maf=0.02
snpformats=false

# Advanced parameters
###################################################
haplome_number=1
maxHaplotype=6
unbiased_downsample=25
biased_downsample=0
softclip=true
```

Alternatively, a configuration file (outlined below) specifying only the pploidy level is sufficient to run GBSapp.

**config**
```
# Variant calling
###################################################
ploidy=2
```
Since most if the parameters are hard-coded in an intuitive manner, by specifying only the ploidy levels, the pipelines determines the all the other parameters as stated below:<br />
- threads: computes available number of cores (n) and n-2 correspond
- defaults: walkaway=true, luster=false, reference genome ID is generated if a reference genome fasta files is provided, copy=1, paleopolyploid=false, genotype and sample missingness rate are 0.2, 2x ploidy uses minRD=6, maf=0.02, no samples are exclude and samples are considered to be unrelated.

## Related Software
- Next-Generation Sequence (NGS) data filtering
    - [ngsComposer: Empirical Base-call error-filtering and read preprocessing pipeline.](https://github.com/ryandkuster/ngsComposer)

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
- $ sudo ln -sf /usr/bin/python3 /usr/bin/python
```
**If samtools and bcftools doesn't install properly:**<br />
```
While the installation of samtools and bcftools are automated, the installation requires some dependencies. Consider typing the commands below in terminal, delete samtools and bcftools (from within ./GBSapp/tools/), and then re-run GBSapp in terminal:
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
**Problem with amount of memory and/or processors/cores specified:**<br />
```
- This might be due to specifing values greater than available resources
- Re-submit job with appropriate values or modify header of the GBSapprun.sh batch file.
- If using compute cluster managers other than SLURM, header of GBSapprun.sh batch can also be modified to fit the syntax of the cluster manager been used.
```
## Versioning
Versioning will follow major.minor.patch <a href="https://semver.org">semantic versioning format</a>.

## License
<a href="https://github.com/ryandkuster/composer/blob/master/LICENSE">Apache License Version 2.0</a>
