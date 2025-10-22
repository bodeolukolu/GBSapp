#!/bin/bash

red=$'\e[1;31m'
green=$'\e[1;32m'
yellow=$'\e[1;33m'
blue=$'\e[1;34m'
magenta=$'\e[1;35m'
cyan=$'\e[1;36m'
white=$'\e[0m'

tools_dir=$(pwd)

main_ngm () {
  echo -e "${blue}\n############################################## ${yellow}\n- downloading and installing NGM ${blue}\n##############################################${white}"
wget https://cmake.org/files/v3.28/cmake-3.28.3-linux-x86_64.tar.gz
tar -xzf cmake-3.28.3-linux-x86_64.tar.gz
rm cmake-3.28.3-linux-x86_64.tar.gz
cd cmake-3.28.3-linux-x86_64
export PATH="$tools_dir/cmake-3.28.3-linux-x86_64/bin:$PATH"
cd $tools_dir
git clone https://github.com/pezmaster31/bamtools.git
cd bamtools && mkdir build && cd build && cmake .. -DCMAKE_CXX_STANDARD=98 && make
cd $tools_dir
wget https://github.com/Cibiv/NextGenMap/tarball/master -O NGM.tar.gz
tar -xvzf NGM.tar.gz
rm NGM.tar.gz
cd Cibiv-NextGenMap*/ || exit 1
mkdir -p build/release
cd build/release || exit 1
cmake ../../ \
  -DCMAKE_PREFIX_PATH="../../../bamtools/build/src/" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_STANDARD=98
make -j$(nproc)
cd "$tools_dir"
rm -rf "$tools_dir/cmake-3.28.3-linux-x86_64/"
cd $tools_dir
}
dirtool=*NextGenMap*
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (NextGenMap aligner: NGM) ${white}"
  main_ngm &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} NGM did not install properly ${white}"
  fi
fi


main_minimap2 () {
  echo -e "${blue}\n############################################## ${yellow}\n- downloading and installing minimap2 ${blue}\n##############################################${white}"
  git clone https://github.com/lh3/minimap2
  cd minimap2 && make
  cd $tools_dir
}
dirtool=minimap2
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (minimap2) ${white}"
  main_minimap2 &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} NGM did not install properly ${white}"
  fi
fi


main_bcftools () {
  echo -e "${blue}\n############################################## ${yellow}\n- downloading and installing bcftools ${blue}\n##############################################${white}"
  wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 &&
  tar -xvf bcftools-1.22.tar.bz2;  cd bcftools-1.22; make
  rm ../bcftools-1.22.tar.bz2
  cd $tools_dir
}
dirtool=bcftools-1.22
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (bcftools) ${white}"
  main_bcftools &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} bcftools did not install properly ${white}"
  fi
fi


main_bedtools () {
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
  main_bedtools &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} bedtools did not install properly ${white}"
  fi
fi


main_samtools () {
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
  main_samtools &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} samtools did not install properly ${white}"
  fi
fi


main_picard () {
  echo -e "${blue}\n############################################## ${yellow}\n- downloading PICARD tools ${blue}\n##############################################${white}"
  wget https://github.com/broadinstitute/picard/releases/download/2.25.6/picard.jar
  cd $tools_dir
}
dirtool=picard*
if [ -f $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of dependency (picard tools) ${white}"
  main_picard &>> ./log.out
  if [ ! -f $dirtool ]; then
      echo -e "${magenta} PICARD tools did not install properly ${white}"
  fi
fi


main_java () {
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
  main_java &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} java did not install properly ${white}"
  fi
fi

main_emboss () {
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
  main_emboss &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} Emboss did not install properly ${white}"
  fi
fi


main_mafft () {
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
  main_mafft &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} mafft did not install properly ${white}"
  fi
fi

main_gatk () {
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
  main_gatk &>> ./log.out
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} GATK did not install properly ${white}"
  fi
fi

main_bwa () {
  echo -e "${white}\n############################################## ${orange}\n- downloading and installing BWA aligner ${white}\n##############################################${white}"
  git clone https://github.com/lh3/bwa.git
  cd bwa; make
  cd $tools_dir
}
dirtool=bwa*
if [ -d $dirtool ]; then
  :
else
  echo -e "${white}- Performing installation of dependency (bwa) ${white}"
  main_bwa &>> ./log.out
fi

main_star () {
  echo -e "${white}\n############################################## ${orange}\n- downloading and installing STAR aligner ${white}\n##############################################${white}"
  wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz
  tar -xzf 2.7.11b.tar.gz
  cd STAR-2.7.11b
  cd STAR/source
  make STAR
  cd $tools_dir*
  rm 2.7.11b.tar.gz
}
dirtool=STAR
if [ -d $dirtool ]; then
  :
else
  echo -e "${white}- Performing installation of dependency (STAR) ${white}"
  main_star &>> ./log.out
fi

main_R () {
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
main_R &>> ./log.out
dirtool=R
if [ -d $dirtool ]; then
  cd $dirtool
  R_dir=$(pwd)
else
  mkdir ./$dirtool
  cd $dirtool
  R_dir=$(pwd)
fi


main_reshape2 () {
  echo -e "${blue}\n############################################## \n- installing R-package: reshape2  ${blue}\n##############################################${white}"
  R -e 'install.packages("reshape2", dependencies = TRUE, lib="./")'
  cd $R_dir
}
dirtool=./reshape2
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of R-package: reshape2 ${white}"
  main_reshape2 &>> ./log.out
  cd $R_dir
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} R-package: reshape2 did not install properly ${white}"
  fi
fi


main_ggplot2 () {
  echo -e "${blue}\n############################################## \n- installing R-package: ggplot2  ${blue}\n##############################################${white}"
  R -e 'install.packages("ggplot2", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  cd ../
}
dirtool=./ggplot2
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of R-package: ggplot2 ${white}"
  main_ggplot2 &>> ./log.out
  cd $R_dir
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} R-package: ggplot2 did not install properly ${white}"
  fi
fi


main_CMplot () {
  echo -e "${blue}\n############################################## \n- installing R-package: CMplot  ${blue}\n##############################################${white}"
  R -e 'install.packages("CMplot", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  cd ../
}
dirtool=./CMplot
if [ -d $dirtool ]; then
  :
else
  echo -e "${magenta}- Performing installation of R-package: CMplot ${white}"
  main_CMplot &>> ./log.out
  cd $R_dir
  if [ ! -d $dirtool ]; then
      echo -e "${magenta} R-package: CMplot did not install properly ${white}"
  fi
fi
