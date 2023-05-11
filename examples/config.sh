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
ref1=TF.fasta
ref2=TL.fasta
ploidy_ref1=4
ploidy_ref2=2
Get_Chromosome=TF_Chr01,TF_Chr02
Exclude_Chromosome=TF_Chr00,TL_Chr00

# SNP-filtering:
####################################################
p1=M9
p2=M19
biallelic=false
genotype_missingness=0.1,0.2,0.3
sample_missingness=0.1,0.2,0.3
exclude_samples=S1,S2,S3
minRD_2x=6
minRD_4x=25
minRD_6x=45
pseg=0.001
maf=0.05

# Advanced_parameters
###################################################
uniquely_mapped=true
paralogs=false
minmapq=20
downsample_2x=50
downsample_4x=100
downsample_6x=150
downsample_8x=200
variant_intervals=true
interval_list=variant_intervals.list
interval_list_ref1=variant_intervals_TF.list
interval_list_ref2=variant_intervals_TL.list
maxHaplotype=128
joint_calling=false
keep_gVCF=false
RE1=TGCAT
RE2=CATG
