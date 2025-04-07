#!/bin/bash

red=$'\e[1;31m'
green=$'\e[1;32m'
yellow=$'\e[1;33m'
blue=$'\e[1;34m'
magenta=$'\e[1;35m'
cyan=$'\e[1;36m'
white=$'\e[0m'

tools_dir=$(pwd)

main0 () {
  echo -e "${blue}\n############################################## ${yellow}\n- downloading and installing bamtools ${blue}\n##############################################${white}"
  wget https://cmake.org/files/v2.8/cmake-2.8.0-Linux-x86_64.tar.gz
  tar -zxvf cmake-2.8.0-Linux-x86_64.tar.gz
  rm cmake-2.8.0-Linux-x86_64.tar.gz
  export PATH="$tools_dir/cmake-2.8.0-Linux-i386/bin:$PATH"
  git clone https://github.com/pezmaster31/bamtools.git
  wget https://github.com/pezmaster31/bamtools/archive/refs/heads/master.zip &&
  unzip master.zip; rm master.zip; mv bamtools-master bamtools
  cd bamtools; mkdir build
  cd build
  cmake ..
  make -j$(nproc)
  cd $tools_dir
}
dirtool=bamtools*
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (bamtools) ${white}"
  main0 &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} bamtools did not install properly ${white}"
  fi
fi

main1 () {
  echo -e "${blue}\n############################################## ${yellow}\n- downloading and installing NGM ${blue}\n##############################################${white}"
  wget https://cmake.org/files/v2.8/cmake-2.8.0-Linux-x86_64.tar.gz
  tar -zxvf cmake-2.8.0-Linux-x86_64.tar.gz
  rm cmake-2.8.0-Linux-x86_64.tar.gz
  export PATH="$tools_dir/cmake-2.8.0-Linux-i386/bin:$PATH"
  wget https://github.com/Cibiv/NextGenMap/tarball/master -O NGM.tar.gz
  tar -xvzf NGM.tar.gz
  rm NGM.tar.gz
  cd Cibiv-NextGenMap*/ || exit 1
  mkdir -p build/release
  cd build/release || exit 1
  cmake ../../ -DCMAKE_PREFIX_PATH="../../../bamtools/build/src/" -DCMAKE_BUILD_TYPE=Release
  make -j$(nproc)
  cd "$tools_dir"
  rm -rf "$tools_dir/cmake-2.8.0-Linux-x86_64/"
  cd $tools_dir
}
dirtool=*NextGenMap*
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (NextGenMap aligner: NGM) ${white}"
  main1 &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} NGM did not install properly ${white}"
  fi
fi


main2 () {
  echo -e "${blue}\n############################################## ${yellow}\n- downloading and installing bcftools ${blue}\n##############################################${white}"
  wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 &&
  tar -xvf bcftools-1.17.tar.bz2;  cd bcftools-1.17; make
  rm ../bcftools-1.17.tar.bz2
  cd $tools_dir
}
dirtool=bcftools-1.17
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (bcftools) ${white}"
  main2 &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} bcftools did not install properly ${white}"
  fi
fi


main3 () {
  echo -e "${blue}\n############################################## ${yellow}\n- downloading and installing bedtools ${blue}\n##############################################${white}"
  wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz &&
  tar -zxvf bedtools-2.29.1.tar.gz; cd bedtools2; make
  rm ../bedtools-2.29.1.tar.gz
  cd $tools_dir
}
dirtool=bedtools*
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (bedtools) ${white}"
  main3 &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} bedtools did not install properly ${white}"
  fi
fi


main4 () {
  echo -e "${blue}\n############################################## ${yellow}\n- downloading and installing samtools ${blue}\n##############################################${white}"
  wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 &&
  tar -xvf samtools-1.17.tar.bz2;  cd samtools-1.17; make
  rm ../samtools-1.17.tar.bz2
  cd $tools_dir
}
dirtool=samtools*
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (samtools) ${white}"
  main4 &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} samtools did not install properly ${white}"
  fi
fi


main5 () {
  echo -e "${blue}\n############################################## ${yellow}\n- downloading PICARD tools ${blue}\n##############################################${white}"
  wget https://github.com/broadinstitute/picard/releases/download/2.25.6/picard.jar
  cd $tools_dir
}
dirtool=picard*
if [ -f $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (picard tools) ${white}"
  main5 &>> ./log.out
  if [ ! -f $dirtool ]; then
      echo -e "${magenta} PICARD tools did not install properly ${white}"
  fi
fi


main () {
  echo -e "${blue}\n############################################## \n- installiing java 17.0.14+7 ${blue}\n##############################################${white}"
  wget https://github.com/adoptium/temurin17-binaries/releases/download/jdk-17.0.14%2B7/OpenJDK17U-jdk_x64_linux_hotspot_17.0.14_7.tar.gz
  tar -xvf OpenJDK17U-jdk_x64_linux_hotspot_17.0.14_7.tar.gz; rm *tar.gz
  cd $tools_dir
}
dirtool=*jdk*
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (java 1.8)${white}"
  main &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} java did not install properly ${white}"
  fi
fi

main6 () {
  echo -e "${blue}\n############################################## \n- installiing Emboss ${blue}\n##############################################${white}"
  wget -m 'ftp://emboss.open-bio.org/pub/EMBOSS/'
  mv emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz ./
  rm -rf emboss.open-bio.org/
  gunzip EMBOSS-6.6.0.tar.gz
  tar xvf EMBOSS-6.6.0.tar
  cd EMBOSS-6.6.0/
  ./configure --without-x
  make
  cd ../
  rm EMBOSS-6.6.0.tar
  cd $tools_dir
}
dirtool=EMBOSS*
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (EMBOSS)${white}"
  main6 &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} Emboss did not install properly ${white}"
  fi
fi


main7 () {
  echo -e "${blue}\n############################################## \n- installing mafft ${blue}\n##############################################${white}"
  wget https://mafft.cbrc.jp/alignment/software/mafft-7.505-without-extensions-src.tgz
  tar zxvf mafft-7.505-without-extensions-src.tgz
  rm mafft-7.505-without-extensions-src.tgz
  cd mafft-7.505-without-extensions/core/
  awk 'NR!=3{print $0}' Makefile | awk 'NR>1{print $0}' | cat <(printf "BINDIR = ${GBSapp_dir}tools/mafft\n") - | \
  cat <(printf "PREFIX = ${GBSapp_dir}tools/mafft\n") - > makefile.tmp
  mv makefile.tmp Makefile
  make clean
  make
  make install
  cd $tools_dir
}
dirtool=mafft-7.505-without-extensions
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (mafft)${white}"
  main7 &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} mafft did not install properly ${white}"
  fi
fi

main8 () {
  echo -e "${green}\n############################################## \n- downloading GATK \n##############################################${white}"
  wget -O GATK4.6.1.0.zip "https://github.com/broadinstitute/gatk/releases/download/4.6.1.0/gatk-4.6.1.0.zip"
  unzip GATK4.6.1.0.zip
  rm GATK4.6.1.0.zip
  cd $tools_dir
}
dirtool=gatk*
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (GATK) ${white}"
  main8 &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} GATK did not install properly ${white}"
  fi
fi


main0 () {
  if which R; then
    :
  else
    module add R
    if which R; then
      :
    else
      echo -e "${magenta}- install R before proceeding ${white}" > /dev/tty
      echo -e "${magenta}- dependencies for R in linux: <sudo apt install libcurl4-openssl-dev> and <sudo apt install libssl-dev> ${white}" > /dev/tty
      exit 1
    fi
  fi
}
echo -e "${white}\n############################################## ${orange}\n- check for R installation ${white}\n##############################################${white}" &>> ./log.out
main0 &>> ./log.out


dirtool=R
if [ -d $dirtool ]; then
  cd $dirtool
  R_dir=$(pwd)
else
  mkdir ./$dirtool
  cd $dirtool
  R_dir=$(pwd)
fi


main9 () {
  echo -e "${blue}\n############################################## \n- installing R-package: reshape2  ${blue}\n##############################################${white}"
  R -e 'install.packages("reshape2", dependencies = TRUE, lib="./")'
  cd $R_dir
}
dirtool=./reshape2
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of R-package: reshape2 ${white}"
  main9 &>> ./log.out
  cd $R_dir
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} R-package: reshape2 did not install properly ${white}"
  fi
fi


main10 () {
  echo -e "${blue}\n############################################## \n- installing R-package: ggplot2  ${blue}\n##############################################${white}"
  R -e 'install.packages("ggplot2", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  cd ../
}
dirtool=./ggplot2
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of R-package: ggplot2 ${white}"
  main10 &>> ./log.out
  cd $R_dir
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} R-package: ggplot2 did not install properly ${white}"
  fi
fi


main11 () {
  echo -e "${blue}\n############################################## \n- installing R-package: CMplot  ${blue}\n##############################################${white}"
  R -e 'install.packages("CMplot", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  cd ../
}
dirtool=./CMplot
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of R-package: CMplot ${white}"
  main11 &>> ./log.out
  cd $R_dir
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} R-package: CMplot did not install properly ${white}"
  fi
fi
