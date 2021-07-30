#!/bin/bash
red=$'\e[1;31m'
green=$'\e[1;32m'
yellow=$'\e[1;33m'
blue=$'\e[1;34m'
magenta=$'\e[1;35m'
cyan=$'\e[1;36m'
white=$'\e[0m'


######################################################################################################################################################
main () {
  echo -e "${blue}\n############################################## ${yellow}\n- downloading and installing BWA ${blue}\n##############################################${white}"
  wget https://sourceforge.net/projects/bio-bwa/files/latest/download
  tar -vxjf download*; rm download*; cd bwa*; wait; make; wait; cd ..
}
dirtool=bwa*
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (bwa) ${white}"
  main &>> ./log.out
fi


main () {
  echo -e "${blue}\n############################################## ${yellow}\n- downloading and installing samtools ${blue}\n##############################################${white}"
  wget https://sourceforge.net/projects/samtools/files/latest/download &&
  tar -vxjf download*; rm download*; cd samtools*; wait; make; wait; cd ..
}
dirtool=samtools*
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (samtools) ${white}"
  main &>> ./log.out
fi

main () {
  echo -e "${blue}\n############################################## ${yellow}\n- downloading and installing bcftools ${blue}\n##############################################${white}"
  wget https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2
  tar -vxjf bcftools-1.11.tar.bz2; rm bcftools-1.11.tar.bz2
  for i in $(ls -d bcftools*); do
	   cd bcftools*; ./configure --prefix=${GBSapp_dir}/tools/${i}/; wait; make; wait; make install; wait; cd ..
  done

}
dirtool=bcftools*
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (bcftools) ${white}"
  main &>> ./log.out
fi


main () {
  echo -e "${blue}\n############################################## ${yellow}\n- downloading PICARD tools ${blue}\n##############################################${white}"
  wget https://github.com/broadinstitute/picard/releases/download/2.25.6/picard.jar
}
dirtool=picard*
if [ -f $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (picard tools) ${white}"
  main &>> ./log.out
fi


main () {
  echo -e "${blue}\n############################################## \n- installiing java 1.8. ${blue}\n##############################################${white}"
  wget https://github.com/AdoptOpenJDK/openjdk8-binaries/releases/download/jdk8u222-b10/OpenJDK8U-jre_x64_linux_hotspot_8u222b10.tar.gz
  tar -xvf OpenJDK8U-jre_x64_linux_hotspot_8u222b10.tar.gz; rm *tar.gz
}
dirtool=jdk*
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (java 1.8)${white}"
  main &>> ./log.out
fi

dirtool=R
if [ -d $dirtool ]; then
  :
else
  mkdir ./R
fi


main () {
  echo -e "${green}\n############################################## \n- downloading GATK \n##############################################${white}"
  wget -O GATK4.1.9.0.zip "https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip"
  unzip GATK4.1.9.0.zip
  rm GATK4.1.9.0.zip
  # javaloc=${GBSapp_dir}/jdk8u222-b10-jre/bin/java
  # cd gatk*
  # awk -v javaloc=$javaloc 'NR==208{gsub(/java/,javaloc)}1' gatk > temp
  # mv temp gatk
  # chmod 777 *
  # cd ../

}
dirtool=gatk*
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (GATK) ${white}"
  main &>> ./log.out
fi




main () {
echo -e "${blue}\n############################################## \n- installing R-package: ggplot2  ${blue}\n##############################################${white}"
  cd ./R
  R -e 'install.packages("ggplot2", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
}
dirtool=./R/ggplot2
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of R-package: ggplot2 ${white}"
  main &>> ./log.out
fi
