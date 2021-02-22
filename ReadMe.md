
# Introduction
GBSapp (v. 0.2.1) is an automated pipeline variant and haplotype calling/filtering. The pipeline integrates existing and novel best practices, some of which can be controlled by user-defined parameters. It ensures accurate dosage-based variant/haplotype calling, fine-scale filtering based on each variant data point (rather than averaging across samples/variants). It optimizes memory and speed at various points in the pipeline such as a novel approach that performs sequence read compression/decompression independently on each unique read before and after pre-processing. Intermediate summary reports and visualizations allow for QC at each step of the pipeline.

<p align="center">
<img src="https://github.com/bodeolukolu/GBSapp/blob/master/misc/GBSapp_logo.PNG" width="690" height="444">
</p>

For questions, bugs, and suggestions, please contact bolukolu@utk.edu.

## Features
- Fully-automated: “walk-away” and “walk-through” mode.
- Dosage-based variant/haplotype calling and filtering.
- Haploid (1x), Diploid (2x), Tetraploid (4x), Hexaploid (6x), and Octoploid (8x).
- Captures and codes variable dosage/copy number in paleopolyploids.
- Can restrict variants call to single copy sequences (i.e. multi-locus variant/paralog test and filtering).
- Additional haplotype-based filtering (useful for targeted sequencing of single locus variants).
- Generates variant/haplotype calls and their sequence context.
- Easy to learn, and use.
- Designed by biologists (please don't run away!)

## Contents
- [Installation](#installation)
- [Usage](#usage)
  - [Basic usage](#Basic usage)
  - [Overview of workflow](#Overview of workflow)
  - [Configuration](#Configuration)
- [Related Software](#Related Software)
- [Select Article Referencing GBSapp](#Select Article Referencing GBSapp)
- [Acknowledgment](#Acknowledgment)
- [Troubleshooting](#Troubleshooting)
- [Versioning](#Versioning)
- [License](#License)

## Installation
- Currently, GBSapp is only available for unix-based systems (i.e. macOS and linux).
- Clone or download the Git repository to your desired folder.
```bash
$ git clone https://github.com/bodeolukolu/GBSapp.git
```
- Installation occurs automatically the first time you run the pipeline.
- With the exception of R, all dependencies are intalled to a local directory within GBSapp.


**Dependencies:**<br />
```
bwa, picard, samtools, bcftools, GATK, java, R, R-ggplot2, R-AGHmatrix
```
Install R before running GBSapp. R installation:
```
    - To view R version (or to check if installed), from the terminal type:
          $ R --version

    - For Ubuntu, Install R  using apt:
          $ sudo apt update
          $ sudo apt install r-base

    - For macOS, install using homebrew:
          brew install r
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
- DESIGN IN PROGRESS

### configuration
Using a text editor, save a file containing any of the following variables as 'config' file (no file extension) and include it in your project directory.

**General parameters**

|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|threads|cores-2|numer of cores/processors|integer|Optional|
|walkaway|true|run in walk-away or walk-through mode|true or false|Optional|
|cluster|false|run on compute cluster node or workstation|true or false|Optional|

**Variant calling parameters**

|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|ploidy|na|ploid level = 1,2,4,6, or 8 (ploidy>8 allowed for autopolyploid)|integer|Required|
|ref1|na|reference genome as .fasta file (haploids, diploids, autopolyploids, or 1st subgenome in allopolyploids)|integer|Optional|
|ref2|na|2nd reference genome as .fasta file|integer|Optional|
|ref3|na|3rd reference genome as .fasta file|integer|Optional|
|ref4|na|4th reference genome as .fasta file|integer|Optional|
|copy|1|maximum copy number of sequence in genome/subgenome|integer|Optional|
|paleopolyploid|false|capture/code variable dosage/copy number (i.e. 2x,4x,6x, and 8x)|true or false||Optional|
|maxn_popallele|500|maximum number of alleles at a locus across individuals in a population|integer|Optional|
|ncohorts|1|number of cohorts for Joint-genotyping. ncohorts=no (very slow) indicates multi-sample variant calling without generating single-sample gVCF files |integer, yes, or no|Optional|
|scaleRD|1|downsampling approach that ensures all unique reads are represented in exact proportions |integer|Optional|
|maxRD|.|threshold for reads with excessive coverage (likely derived from paralogs/repetitive sequences)|integer|Optional|
|minRD|1|threshold for number of times a unique occurs. Attempts to eliminates haplotypes derived from bad base calls)|integer|Optional|

**Note: na indicates that variable is user-defined or hard-coded as a function of ploidy.*

**Variant filtering parameters**

|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|p1|na|maternal parent (specified only for biparental populations)|string|Optional|
|p2|na|paternal parent (specified only for biparental populations)|string|Optional|
|genotype_missingness|0.2|maximum proportion of missing genotypes allowed per sample (multiple comma delimited thresholds permitted)|decimal number|Optional|
|sample_missingness|0.2|maximum proportion of missing samples allowed per variant (multiple comma delimited thresholds permitted)|decimal number|Optional|
|minRD_1x|2|minimum read depth threshold|integer|Optional|
|minRD_2x|6|minimum read depth threshold|integer|Optional|
|minRD_4x|25|minimum read depth threshold|integer|Optional|
|minRD_6x|45|minimum read depth threshold|integer|Optional|
|minRD_8x|100|minimum read depth threshold|integer|Optional|
|maf|0.02|minor allele frequency threshold|decimal number|Optional|
|snpformats|false|variant data set with alleles formatted as base (A,C,G,T) and/or degenerate notation|true or false|Optional|
|exclude_samples|na|sample IDs to be exclude from filtered variant data set|comma delimited strings|Optional|

Below is an example of a configuration file:

**config**
```
# General parameters
####################################################
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
copy=1
paleopolyploid=false

# SNP-filtering:
####################################################
p1=Beauregard
p2=Tanzania
genotype_missingness=0.1,0.2,0.3
sample_missingness=0.1,0.2,0.3
minRD_2x=6
minRD_4x=25
minRD_6x=45
maf=0.02
snpformats=false
exclude_samples=S1,S1,S1
```

Alternatively, a configuration file (outlined below) may only need to include the parameters for a run.

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
1. Genome-wide association study identified candidate genes controlling continuous storage root formation and bulking in hexaploid sweetpotato. [Bararyenya et al. 2020](https://doi.org/10.1186/s12870-019-2217-9)
2. Sequencing depth and genotype quality: accuracy and breeding operation considerations for genomic selection applications in autopolyploid crops [Gemenet et al. 2020](https://doi.org/10.1007/s00122-020-03673-2)
3. Genetic Diversity and Population Structure of the USDA Sweetpotato (Ipomoea batatas) Germplasm Collections Using GBSpoly [Wadl et al. 2019](https://doi:10.3389/fpls.2018.01166)

## Acknowledgment
This package has been developed as part of the [Genomic Tools for Sweetpotato Improvement project](https://sweetpotatogenomics.cals.ncsu.edu/) (GT4SP) and [SweetGAINS](https://cgspace.cgiar.org/handle/10568/106838), both funded by [Bill & Melinda Gates Foundation](https://www.gatesfoundation.org/).
ss
## Troubleshooting
- IN PROGRESS

## Versioning
Versioning will follow major.minor.patch <a href="https://semver.org">semantic versioning format</a>.

## License
<a href="https://github.com/ryandkuster/composer/blob/master/LICENSE">Apache License Version 2.0</a>
