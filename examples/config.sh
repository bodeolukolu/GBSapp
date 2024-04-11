# General_parameters
###################################################
threads=16
walkaway=true
cluster=true
nodes=1
samples_alt_dir=false
lib_type=RRS
subsample_WGS_in_silico_qRRS=false

# Variant calling
###################################################
ploidy=6
haplome_number=1,2,3,4,5,6
# Variant calling with haploid subgenome(s)
# Anchored to ref1 for loci conserved across all subgenomes/pangenomes
ref1=TF.fasta
ref2=TL.fasta
ploidy_ref1=4
ploidy_ref2=2
# Variant calling with haplotype-resolved reference genome or pangenomes
# Anchored to haplome_ref1 for loci conserved across all haplomes
genomes_ref=Ib.fasta
# exclue or limit variant calling to specific chromosomes
Get_Chromosome=TF_Chr01,TF_Chr02
Exclude_Chromosome=TF_Chr00,TL_Chr00

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
genomecov_est=false
