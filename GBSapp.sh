#!/bin/bash

magenta=$'\e[1;35m'
white=$'\e[0m'
yellow=$'\e[1;33m'
blue=$'\e[1;34m'
# GBSapp_dir=$(dirname "$0")
# projdir=$1
threads=""
# relpath=$(pwd)

GBSapp_dir="$( cd -- "$(dirname "$0 ")" >/dev/null 2>&1 ; pwd -P )/"
projdir="$( cd -- "$(dirname "$1 ")" >/dev/null 2>&1 ; pwd -P )/"


# if [ "${GBSapp_dir:0:1}" = "." ]; then
#   if [ "${GBSapp_dir:0:15}" = "../../../../../" ]; then
#     GBSapp_dir="${relpath%/*/*/*/*/*}${GBSapp_dir//*..}"
#   fi
#   if [ "${GBSapp_dir:0:12}" = "../../../../" ]; then
#     GBSapp_dir="${relpath%/*/*/*/*}${GBSapp_dir//*..}"
#   fi
#   if [ "${GBSapp_dir:0:9}" = "../../../" ]; then
#     GBSapp_dir="${relpath%/*/*/*}${GBSapp_dir//*..}"
#   fi
#   if [ "${GBSapp_dir:0:6}" = "../../" ]; then
#     GBSapp_dir="${relpath%/*/*}${GBSapp_dir//*..}"
#   fi
#   if [ "${GBSapp_dir:0:3}" = "../" ]; then
#     GBSapp_dir="${relpath%/*}${GBSapp_dir//*..}"
#   fi
#   if [ "${GBSapp_dir:0:2}" != ".." ]; then
#     if [ "${GBSapp_dir:0:1}" = "." ]; then
#       GBSapp_dir="${relpath}${GBSapp_dir:1}"
#     fi
#   fi
# fi
# if [ "${GBSapp_dir: -1}" != "/" ]; then
#   GBSapp_dir="${GBSapp_dir}/"
# fi
#
# if [ "${projdir:0:1}" = "." ]; then
#   if [ "${projdir:0:15}" = "../../../../../" ]; then
#     projdir="${relpath%/*/*/*/*/*}${projdir//*..}"
#   fi
#   if [ "${projdir:0:12}" = "../../../../" ]; then
#     projdir="${relpath%/*/*/*/*}${projdir//*..}"
#   fi
#   if [ "${projdir:0:9}" = "../../../" ]; then
#     projdir="${relpath%/*/*/*}${projdir//*..}"
#   fi
#   if [ "${projdir:0:6}" = "../../" ]; then
#     projdir="{relpath%/*/*}${projdir//*..}"
#   fi
#   if [ "${projdir:0:3}" = "../" ]; then
#     projdir="${relpath%/*}${projdir//*..}"
#   fi
#   if [ "${projdir:0:2}" = "./" ]; then
#     projdir="${relpath}${projdir:1}"
#   fi
#   if [ "${projdir}" = . ]; then
#     projdir="${relpath}"
#   fi
# fi
# if [ "${projdir:-1}" != "/" ]; then
#   projdir="${projdir}/"
# fi

####################################################################################################################
####################################################################################################################
####################################################################################################################

cd $projdir
echo -e "${white}\n##################################################################################\n"
echo -e "${yellow}- Program:	GBSapp"
echo -e "${yellow}- Version:	0.2.5"
echo -e "${yellow}- Description:	Automated Pipeline for Variant/Haplotype Calling and Filtering"
echo -e "${yellow}- Contact:	Bode Olukolu <bolukolu@utk.edu> ${white}"
echo -e "${white}\n##################################################################################\n"

d=($(find . -maxdepth 1 -type d | wc -l))
if [ $d == 3 ]; then
  echo ""
else
  echo -e "${magenta}- Expecting only 2 folders (i.e. refgenomes and samples) ${white}"
  echo -e "${magenta}- Do you want to continue running GBSapp? ${white}"
  read -p "- y(YES) or n(NO)? " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    exit 1
  else
      echo -e "${white}\n##############################################################################"
  fi
fi

f=($(find . -maxdepth 1 -type f | wc -l))
if [ $f == 1 ]; then
  echo ""
else
  echo -e "${magenta}- Expecting only 1 file (i.e. config) ${white}"
  echo -e "${magenta}- Do you want to continue running GBSapp? ${white}"
  read -p "- y(YES) or n(NO)? " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    exit 1
  else
      echo -e "${white}\n\n##############################################################################"
  fi
fi

####################################################################################################################
####################################################################################################################
####################################################################################################################

cd ${GBSapp_dir}
if [[ ! -d "./tools" ]]; then
  mkdir tools
fi

cd $projdir
cd samples
if [ -d "paired" ]; then
 mv paired pe
fi
if [ -d "single" ]; then
 mv single se
fi

# Initial questions before running walkaway
cd $projdir

free=$(df . | awk 'NR==2{print $4}')
required=$(du -s . | awk '{print $1}')
required=$((required * 120/100))
if [[ "$free" -lt "$required" ]]; then
	echo -e "${magenta}- You might not have enough disk space. Free: $((free/1000000))G; Required(approx. 1.2x the size of fastq files): $((required/1000000))G  ${white}\n"
	echo -e "${magenta}- Do you want to continue running GBSapp? ${white}"
	read -p "- Do you want to continue running GBSapp? " -n 1 -r
	if [[ ! $REPLY =~ ^[Yy]$ ]]; then
		printf '\n'
		echo -e "${magenta}- Exiting GBSapp ${white}\n"
		sleep 5 && exit 1
	fi
fi

if [ -d "./samples/pe" ]; then
  echo -e "${magenta}- GBSapp will modify your sample fastq files by concatenating <pe-R1-reads> and <se-reads> ${white}"
  echo -e "${magenta}- Do you want create a copy/backup of your sample fastq file? ${white}"
  read -p "- y(YES) or n(NO)? " -n 1 -r
  if [[ $REPLY = ^[Yy]$ ]]; then
    printf "\n"
    echo -e "${magenta}- A copy of your sample fastq files will be stored in user-defined directory or in ${projdir}/samples_backup/ ${white}"
    echo -e "${magenta}- Specify alternate directory (Absolute/Full path) or Press Enter? ${white}"
    read alt_dir
    if [[ -z $alt_dir ]]; then
    	cp -r samples samples_backup
    else
	cp -r samples $alt_dir
    fi
  else
    echo -e "${white}\n##############################################################################"
  fi
else
  :
fi


####################################################################################################################
####################################################################################################################
####################################################################################################################

main() {
string="${GBSapp_dir}/scripts/GBSapp_internal_parameters.sh"
string2=${string//'//'/'/'}

cd $projdir
awk '{ sub("\r$",""); print}' config > unixformat.sh
mv unixformat.sh config
string="$(awk '/ref/&&/=/&&/.fa/ {print}' config | wc -l)"
if [[ "$string" -eq 0 ]]; then
string=1
fi
string3=${GBSapp_dir}/scripts/subgenome_ref_files_${string}.sh

cd $projdir
printf "#""!""/bin/bash \n\n" > header.txt
printf "GBSapp_dir=${GBSapp_dir}\n" > fetchdir.txt
printf "projdir=${projdir}" >> fetchdir.txt

string="$(grep walkaway config)"
string=${string//*=}
cluster="$(grep cluster config)"
cluster=${cluster//*=}
thread_node="$(grep ^threads config)"
thread_node=${thread_node//*=}

if [ -z "$string" ]; then
 string=true
fi
if [ "$cluster" == false ]; then
 unset cluster
fi





if [[ "$string" == false ]]; then
  echo -e "${magenta}- GBSapp will run in walk-through mode\n\n############################################################################## ${white}"
  if [ -z "$thread_node" ]; then
    echo -e "${magenta}- GBSapp will use (total processors/cores)-2 ${white}"
    echo -e "${magenta}- Quit (ctrl + c), else GBSapp will continue in 10 seconds ${white}"
    sleep 5
  fi

  cat header.txt config fetchdir.txt $string2 $string3 | awk '{ sub("\r$",""); print}' > GBSapp_run.sh
  rm header.txt fetchdir.txt
  bash ${projdir}GBSapp_run.sh
fi





if [[ "$string" == true ]]; then
if [ -z $cluster ]; then
  echo -e "${magenta}- GBSapp will run in walkaway mode\n ${white}"
  echo -e "${white}##############################################################################\n ${white}"
  if [ -z "$thread_node" ]; then
    echo -e "${magenta}- GBSapp will use (total processors/cores)-2 ${white}"
    echo -e "${magenta}- Quit (ctrl + c), else GBSapp will continue in 10 seconds ${white}"
    sleep 5
  fi

  echo -e "${magenta}- Do you want to perform read alignments and alignment post-processing? ${white}"
  read -p "- y(YES) or n(NO) " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  	printf '\n'
	printf "alignments=0 \n" > steps.txt
  	echo -e "${magenta}- skipping read_alignments and alignment post-processing ${white}"
    echo -e "${white}\n##############################################################################\n"
  else
    printf '\n'
    printf "alignments=1 \n" > steps.txt
    echo -e "${white}\n##############################################################################\n"
  fi

  echo -e "${magenta}- Do you want to perform SNP/variant calling? ${white}"
  read -p "- y(YES) or n(NO) " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  	printf '\n'
	printf "snp_calling=0 \n" >> steps.txt
  	echo -e "${magenta}- skipping SNP/variant calling ${white}"
    echo -e "${white}\n##############################################################################\n"
  else
    printf '\n'
    printf "snp_calling=1 \n" >> steps.txt
    echo -e "${white}\n##############################################################################\n"
  fi

  echo -e "${magenta}- Do you want to perform SNP/variant filtering? ${white}"
  read -p "- y(YES) or n(NO) " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  	printf '\n'
	printf "snp_filtering=0 \n" >> steps.txt
  	echo -e "${magenta}- skipping SNP/variant filtering ${white}"
    echo -e "${white}\n##############################################################################\n"
  else
    printf '\n'
    printf "snp_filtering=1 \n" >> steps.txt
    echo -e "${white}\n##############################################################################\n"
  fi

  echo -e "${magenta}- Do you want to generate visualizations for genotype accuracy and ploidy estimation? ${white}"
  read -p "- y(YES) or n(NO) " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  	printf '\n'
	printf "SNPaccuracy_ReadDepth=0 \n" >> steps.txt
  	echo -e "${magenta}- skipping visualizations for genotype accuracy and ploidy estimation ${white}"
    echo -e "${white}\n##############################################################################\n"
  else
    printf '\n'
    printf "SNPaccuracy_ReadDepth=1 \n" >> steps.txt
    echo -e "${white}\n##############################################################################\n"
  fi

  echo -e "${magenta}- Do you want to perform sequence-based haplotyping and filtering? ${white}"
  read -p "- y(YES) or n(NO) " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  	printf '\n'
	printf "paralogfiltering_haplotyping=0 \n" >> steps.txt
  	echo -e "${magenta}- skipping sequence-based haplotyping and filtering ${white}"
    echo -e "${white}\n##############################################################################\n"
  else
    printf '\n'
    printf "paralogfiltering_haplotyping=1 \n" >> steps.txt
    echo -e "${white}\n##############################################################################\n"
  fi

  printf "#""!""/bin/bash \n" > header.txt
  printf "walkaway=true \n\n" > walkaway.txt
  printf "GBSapp_dir=${GBSapp_dir}\n" > fetchdir.txt
  printf "projdir=${projdir}" >> fetchdir.txt
  cat header.txt steps.txt config walkaway.txt fetchdir.txt $string2 $string3 | awk '{ sub("\r$",""); print}' > GBSapp_run.sh
  rm header.txt walkaway.txt fetchdir.txt steps.txt
  nohup bash ${projdir}GBSapp_run.sh > terminal.out 2>&1 &
fi
fi





if [[ "$string" == true ]]; then
if [[ "$cluster" == true ]]; then
  if [ -z "$thread_node" ]; then
    echo -e "${magenta}- Please provide number of threads config (for cluster node) ${white}"
    echo -e "${magenta}- GBSapp will quit ${white}"
    sleep 2
    exit 1
  fi
  echo -e "${magenta}- GBSapp will run in walkaway mode\n ${white}"
  echo -e "${white}##############################################################################\n ${white}"

  echo -e "${magenta}- Do you want to perform read alignments and alignment post-processing? ${white}"
  read -p "- y(YES) or n(NO) " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  	printf '\n'
	printf "alignments=0 \n" > steps.txt
  	echo -e "${magenta}- skipping read_alignments and alignment post-processing ${white}"
    echo -e "${white}\n##############################################################################\n"
  else
    printf '\n'
    printf "alignments=1 \n" > steps.txt
    echo -e "${white}\n##############################################################################\n"
  fi

  echo -e "${magenta}- Do you want to perform SNP/variant calling? ${white}"
  read -p "- y(YES) or n(NO) " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  	printf '\n'
	printf "snp_calling=0 \n" >> steps.txt
  	echo -e "${magenta}- skipping SNP/variant calling ${white}"
    echo -e "${white}\n##############################################################################\n"
  else
    printf '\n'
    printf "snp_calling=1 \n" >> steps.txt
    echo -e "${white}\n##############################################################################\n"
  fi

  echo -e "${magenta}- Do you want to perform SNP/variant filtering? ${white}"
  read -p "- y(YES) or n(NO) " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  	printf '\n'
	printf "snp_filtering=0 \n" >> steps.txt
  	echo -e "${magenta}- skipping SNP/variant filtering ${white}"
    echo -e "${white}\n##############################################################################\n"
  else
    printf '\n'
    printf "snp_filtering=1 \n" >> steps.txt
    echo -e "${white}\n##############################################################################\n"
  fi

  echo -e "${magenta}- Do you want to generate visualizations for genotype accuracy and ploidy estimation? ${white}"
  read -p "- y(YES) or n(NO) " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  	printf '\n'
	printf "SNPaccuracy_ReadDepth=0 \n" >> steps.txt
  	echo -e "${magenta}- skipping visualizations for genotype accuracy and ploidy estimation ${white}"
    echo -e "${white}\n##############################################################################\n"
  else
    printf '\n'
    printf "SNPaccuracy_ReadDepth=1 \n" >> steps.txt
    echo -e "${white}\n##############################################################################\n"
  fi

  echo -e "${magenta}- Do you want to perform sequence-based haplotyping and filtering? ${white}"
  read -p "- y(YES) or n(NO) " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  	printf '\n'
	printf "paralogfiltering_haplotyping=0 \n" >> steps.txt
  	echo -e "${magenta}- skipping sequence-based haplotyping and filtering ${white}"
    echo -e "${white}\n##############################################################################\n"
  else
    printf '\n'
    printf "paralogfiltering_haplotyping=1 \n" >> steps.txt
    echo -e "${white}\n##############################################################################\n"
  fi

  printf "#""!""/bin/bash \n#SBATCH -c ${thread_node} \n\n" > cluster_header.sh
  printf "walkaway=true \n\n" > walkaway.txt
  cat cluster_header.sh steps.txt config walkaway.txt fetchdir.txt $string2 $string3 | awk '{ sub("\r$",""); print}' > GBSapp_run.sh
  rm fetchdir.txt cluster_header.sh header.txt steps.txt walkaway.txt
  slurm_check=$(sbatch --version)
  if [[ $slurm_check == slurm* ]]; then
    sbatch ${projdir}GBSapp_run.sh
  else
    echo -e "${magenta}- Cluster manager not SLURM, revise batch file header in ${projdir}GBSapp_run.sh to for cluster manager syntax ${white}"
    echo -e "${magenta}- Revise batch file header in ${projdir}GBSapp_run.sh for your cluster manager ${white}"
  fi
fi
fi
}
cd $projdir
time main

echo -e "${magenta}- For jobs running in background, monitor progress in terminal.out or slurm-xxxxxx.out (slurm cluster manager) ${white}"
echo -e "${magenta}- The log.out file can help with troubleshooting and when reporting a bug ${white}"
