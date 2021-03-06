
red=$'\e[1;31m'
green=$'\e[1;32m'
yellow=$'\e[1;33m'
blue=$'\e[1;34m'
magenta=$'\e[1;35m'
cyan=$'\e[1;36m'
white=$'\e[0m'

######################################################################################################################################################
export GBSapp_dir=${GBSapp_dir%/*}
export projdir=${projdir%/*}

pop=${projdir%/}
pop=${pop##*/}


######################################################################################################################################################
# Software defined parameters
mkdir -p ${GBSapp_dir}/tools
cd ${GBSapp_dir}/tools
bash ../scripts/install.sh

######################################################################################################################################################
# tools
export bwa=${GBSapp_dir}/tools/bwa*/bwa && bwa=${bwa//'//'/'/'}
export samtools=${GBSapp_dir}/tools/samtools*/samtools && samtools=${samtools//'//'/'/'}
export bcftools=${GBSapp_dir}/tools/bcftools*/bcftools && bcftools=${bcftools//'//'/'/'}
picard=${GBSapp_dir}/tools/picard.jar && picard=${picard//'//'/'/'}
GATK=${GBSapp_dir}/tools/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar && GATK=${GATK//'//'/'/'}
java=${GBSapp_dir}/tools/jdk8*/bin/java && java=${java//'//'/'/'}
