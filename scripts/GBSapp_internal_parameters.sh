
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
export pop=${pop##*/}

cluster="$(grep cluster config.sh)"
export cluster=${cluster//*=}

######################################################################################################################################################
# Software defined parameters

export slurm_module=$((module --version) 2>&1 | head -n1)

if [[ "$slurm_module" =~ "Module" ]]; then
	module unload python
	module add python
	pythonversion=$((python --version) 2>&1)
	if [[ "$pythonversion" =~ "Python 3" ]]; then
		echo -e "${white}\n- Using $pythonversion\n ${white}"
	fi
	module unload R 2> /dev/null
  module add R 2> /dev/null
	module add r 2> /dev/null
  Rversion=$((R --version) 2>&1)
  if [[ "$Rversion" =~ "R version" ]]; then
    echo -e "${white}\n- Using $Rversion\n ${white}"
  fi
fi

Rversion=$((R --version) 2>&1)
if [[ "$Rversion" =~ "R version" ]]; then
  echo -e "${white}\n- Using $Rversion\n ${white}"
else
  echo -e "${white}- install R before proceeding ${white}"
  echo -e "${white}- dependencies for R in linux: <sudo apt install libcurl4-openssl-dev> and <sudo apt install libssl-dev>"
fi

pythonversion=$((python --version) 2>&1)
if [[ "$pythonversion" =~ "Python 3" ]]; then
  echo -e "${white}\n- Using $pythonversion\n ${white}"
else
	echo -e "${white}- install python before proceeding ${white}"
fi

javaversion=$((update-alternatives --list java | grep 'java-17-openjdk') 2>&1)
mkdir -p ~/bin
PATH=~/bin:$PATH
rm ~/bin/java
ln -s $javaversion ~/bin/java
javaversion=$((java -version) 2>&1)
if [[ "$javaversion" =~ "17.0.14+7" ]]; then
	echo -e "${white}\n- Using $javaversion\n ${white}"
else
	mkdir -p ~/bin
	PATH=~/bin:$PATH
	rm ~/bin/java
	ln -s ${GBSapp_dir}/tools/jdk*/bin/java ~/bin/java
	javaversion=$((java -version) 2>&1)
	if [[ "$javaversion" =~ "17.0.14+7" ]]; then
		echo -e "${white}\n- Using $javaversion\n ${white}"
	else
		echo -e "${white}- install java version build 17.0.14+7 before proceeding ${white}"
	fi
fi
mkdir -p ~/bin
PATH=~/bin:$PATH
rm ~/bin/java
ln -s ${GBSapp_dir}/tools/*jdk*/bin/java ~/bin/java




######################################################################################################################################################
# tools
export genmap=${GBSapp_dir}/tools/genmap-build/bin/genmap
export ngm=${GBSapp_dir}/tools/*NextGenMap*/bin/ngm*/ngm
export star=${GBSapp_dir}/tools/STAR*/source/STAR
export minimap2=${GBSapp_dir}/tools/minimap2*/minimap2
export paftools=${GBSapp_dir}/tools/minimap2*/minimap2/misc/paftools.js
export wig2bed=${GBSapp_dir}/tools/wig2bed_merged.py
export bwa=${GBSapp_dir}/tools/bwa/bwa
export samtools=${GBSapp_dir}/tools/samtools*/samtools && samtools=${samtools//'//'/'/'}
export bcftools=${GBSapp_dir}/tools/bcftools*/bcftools && bcftools=${bcftools//'//'/'/'}
export bedtools=${GBSapp_dir}/tools/bedtools2/bin/bedtools && bedtools=${bedtools//'//'/'/'}
export picard=${GBSapp_dir}/tools/picard.jar && picard=${picard//'//'/'/'}
export GATK=${GBSapp_dir}/tools/gatk*/gatk && GATK=${GATK//'//'/'/'}
export java=${GBSapp_dir}/tools/*jdk*/bin/java && java=${java//'//'/'/'}
export consambig=${GBSapp_dir}/tools/EMBOSS-6.6.0/emboss/consambig
export mafft=${GBSapp_dir}/tools/mafft/mafft
if command -v pigz &>/dev/null; then
  export gzip=pigz
else
	export gzip=gzip
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

bedtoolsout=$($bedtools --version | head -n 3)
if [ -z "$bedtoolsout" ];then
  echo -e "${white}- bedtools installation within GBSapp is probably missing a dependency on host system ${white}"
  echo -e "${white}- GBSapp will use host system bcftools installation ${white}"
  if [[ "$cluster" == true ]]; then
    module add bedtools
  fi
  export bedtools=bedtools
  $bedtools --version | head -n 3
else
  $bedtools --version | head -n 3
fi

cd $GBSapp_dir/tools/R
Rscript ../../scripts/R/check_R_tools.R "${GBSapp_dir}/tools/R"
cd $projdir
