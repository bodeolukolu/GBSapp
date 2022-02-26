#!/bin/bash

red=$'\e[1;31m'
green=$'\e[1;32m'
yellow=$'\e[1;33m'
blue=$'\e[1;34m'
magenta=$'\e[1;35m'
cyan=$'\e[1;36m'
white=$'\e[0m'


main () {
  echo -e "${blue}\n############################################## ${yellow}\n- downloading and installing NGM ${blue}\n##############################################${white}"
  wget https://github.com/Cibiv/NextGenMap/tarball/master -O NGM.tar.gz
  tar xvfz NGM.tar.gz; rm NGM.tar.gz
  cd Cibiv-NextGenMap*
  mkdir -p build/
  cd build/
  cmake ..
  make
}
dirtool=bbmap
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (BBmap) ${white}"
  main &>> ./log.out
fi


main () {
  echo -e "${blue}\n############################################## ${yellow}\n- downloading and installing bcftools ${blue}\n##############################################${white}"
  wget https://github.com/samtools/bcftools/releases/download/1.13/bcftools-1.13.tar.bz2 &&
  tar -xvf bcftools-1.13.tar.bz2;  cd bcftools-1.13; make
  rm ../bcftools-1.13.tar.bz2
  cd ..
}
dirtool=bcftools-1.13
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (bcftools) ${white}"
  main &>> ./log.out
fi


main () {
  echo -e "${blue}\n############################################## ${yellow}\n- downloading and installing samtools ${blue}\n##############################################${white}"
  wget https://sourceforge.net/projects/samtools/files/latest/download
  tar -xvjf download*; rm download*; cd samtools*
  make
  cd ..
}
dirtool=samtools*
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (samtools) ${white}"
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
  wget https://github.com/adoptium/temurin8-binaries/releases/download/jdk8u322-b06/OpenJDK8U-jdk_x64_linux_hotspot_8u322b06.tar.gz
  tar -xvf OpenJDK8U-jdk_x64_linux_hotspot_8u322b06.tar.gz; rm *tar.gz
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
  wget
  wget -O GATK4.2.2.0.zip "https://github.com/broadinstitute/gatk/releases/download/4.2.2.0/gatk-4.2.2.0.zip"
  unzip GATK4.2.2.0.zip
  rm GATK4.2.2.0.zip
}
dirtool=gatk*
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (GATK) ${white}"
  main &>> ./log.out
fi


main () {
  echo -e "${white}\n############################################## ${orange}\n- check for R installation ${white}\n##############################################${white}"
  if R --version; then
    :
  else
    module add R
    if R --version; then
      :
    else
      echo -e "${white}- install R before proceeding ${white}"
      echo -e "${white}- dependencies for R in linux: <sudo apt install libcurl4-openssl-dev> and <sudo apt install libssl-dev>"
    fi
  fi
}
main &>> ./log.out


main () {
echo -e "${blue}\n############################################## \n- installing R-package: ggplot2  ${blue}\n##############################################${white}"
  mkdir -p R
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
