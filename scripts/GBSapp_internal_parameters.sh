
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

cluster="$(grep cluster config.sh)"
cluster=${cluster//*=}

######################################################################################################################################################
# Software defined parameters

Rout=$(R --version | head -n 3)
if [ -z "$Rout" ];then
  echo -e "${white}- R not available ${white}"

	if [[ "$cluster" == true ]]; then
    module add R
  fi
	R --version | head -n 3
fi

Rout=$(R --version | head -n 3)
if [ -z "$Rout" ];then
  echo -e "${white}- install R before proceeding ${white}"
  echo -e "${white}- dependencies for R in linux: <sudo apt install libcurl4-openssl-dev> and <sudo apt install libssl-dev>"
fi

pythonout=$(python2 --version | head -n 3)
if [ -z "$pythonout" ];then
	module add python2
	python --version | head -n 3
fi
pythonout=$(python2 --version | head -n 3)
if [ -z "$pythonout" ];then
  echo -e "${white}- install Python2 before proceeding ${white}"
fi

######################################################################################################################################################
# tools
export bbmap=${GBSapp_dir}/tools/bbmap/bbmap.sh
export bbwrap=${GBSapp_dir}/tools/bbmap/bbwrap.sh
export samtools=${GBSapp_dir}/tools/samtools*/samtools && samtools=${samtools//'//'/'/'}
export bcftools=${GBSapp_dir}/tools/bcftools*/bcftools && bcftools=${bcftools//'//'/'/'}
export picard=${GBSapp_dir}/tools/picard.jar && picard=${picard//'//'/'/'}
export GATK=${GBSapp_dir}/tools/gatk-4.2.2.0/gatk && GATK=${GATK//'//'/'/'}
export java=${GBSapp_dir}/tools/jdk8*/bin/java && java=${java//'//'/'/'}

alias python=python2

if command -v pigz &>/dev/null; then
  export gzip=pigz
  export gunzip=unpigz
  export zcat="unpigz -c"
else
  export gzip=gzip
  export gunzip=gunzip
  export zcat=zcat
fi


samtoolsout=$($samtools --version | head -n 3)
if [ -z "$samtoolsout" ];then
  echo -e "${white}- samtools installation within GBSapp is probably missing a dependency on host system ${white}"
  echo -e "${white}- GBSapp will use host system samtools installation ${white}"
  if [[ "$cluster" == true ]]; then
    module add samtools
  fi
  export samtools=samtools
  $samtools --version | head -n 3
else
  $samtools --version | head -n 3
fi

bcftoolsout=$($bcftools --version | head -n 3)
if [ -z "$bcftoolsout" ];then
  echo -e "${white}- bcftools installation within GBSapp is probably missing a dependency on host system ${white}"
  echo -e "${white}- GBSapp will use host system bcftools installation ${white}"
  if [[ "$cluster" == true ]]; then
    module add bcftools
  fi
  export bcftools=bcftools
  $bcftools --version | head -n 3
else
  $bcftools --version | head -n 3
fi
