# General_parameters
###################################################
threads=24
walkaway=true
cluster=true
nodes=1
samples_alt_dir=false

# Variant calling
###################################################
ploidy=6
ref1=TF.fasta
ref2=TL.fasta
ploidy_ref1=4
ploidy_ref2=2


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

# Advanced_parameters
###################################################
maxHaplotype=128
haplome_number=1
mhap_freq=1
altpos=false
softclip=true
ncohorts=1
keep_gVCF=false
maxindel=100
PEdist=250