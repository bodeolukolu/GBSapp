# General_parameters
###################################################
threads=16
walkaway=true
cluster=false
nodes=1
RNA=false
variant_caller=gatk
samples_alt_dir=false
lib_type=RRS
subsample_WGS_in_silico_qRRS=false

# Variant calling
###################################################
ploidy=6
ref1=TF.fasta
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
minRD_6x=45
pseg=0.001
maf=0.05
filtered_vcf=true

# Advanced_parameters
###################################################
max_pseudoMol=1000
uniquely_mapped=true
paralogs=true
downsample_6x=150
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
