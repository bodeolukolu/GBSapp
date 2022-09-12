# General_parameters
###################################################
threads=24
walkaway=true
cluster=true
nodes=1
samples_alt_dir=false
lib_type=RRS

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
multilocus=false
minmmapq=20
maxHaplotype=128
haplome_number=1
softclip=false
joint_calling=false
keep_gVCF=false
RE1=TGCAT
RE2=CATG
