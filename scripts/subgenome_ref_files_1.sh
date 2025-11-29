if [ -z "$slurm_module" ]; then
 export slurm_module=true
fi
if [ -z "$threads" ]; then
	export threads=$(nproc --all)
	if [[ "$threads" -ge 4 ]]; then
		export threads=$((threads-2))
	fi
fi
if [ -z "$aligner" ]; then
  export aligner=minimap2
fi
if [ -z "$RNA" ]; then
 export RNA=false
fi
if [ -z "$variant_caller" ]; then
 export variant_caller=gatk
fi
if [ "$variant_caller" == "GATK" ]; then
 export variant_caller=gatk
fi
if [ "$variant_caller" == "BCFtools" ] || [ "$variant_caller" == "BCFTOOLS" ]; then
 export variant_caller=bcftools
fi
if [ "$variant_caller" == "bcftools" ]; then
 export ploidy=2
 export ploidy_ref1=2
fi
if [ -z "$lib_type" ]; then
 export lib_type=RRS
fi
if [ -z "$subsample_WGS_in_silico_qRRS" ]; then
 export subsample_WGS_in_silico_qRRS=false
fi
if [ -z "$nodes" ]; then
 export nodes=1
fi
if [ -z "$biallelic" ]; then
	export biallelic=false
fi
if [ -z "$paralogs" ]; then
	export paralogs=false
fi
if [ -z "$max_pseudoMol" ]; then
	export max_pseudoMol=1000
fi
if [ -z "$uniquely_mapped" ]; then
	export uniquely_mapped=true
fi
if [ -z "$minmapq" ]; then
	export minmapq=20
fi
if [ "$RNA" == "true" ]; then
	export minmapq=1
fi
if [ -z "$maxHaplotype" ]; then
	export maxHaplotype=128
fi
if [ -z "$haplome_number" ]; then
	export haplome_number=1
fi
if [[ "$hap_ref" ]]; then
  export hap_ref=$genome_ref
fi
if [[ "$hap_ref" ]]; then
  export ref1=$hap_ref
fi
if [ -z "$p2" ]; then
  if   [ "$p1" ]; then
	 export p2=$p1
 fi
fi
if [ -z "$use_softclip" ]; then
	export use_softclip=false
fi
if [ "$use_softclip" == "false" ]; then
	export dont_use_softclip=true
fi
if [ "$use_softclip" == "true" ]; then
	export dont_use_softclip=false
fi

if [ -z "$joint_calling" ]; then
	export joint_calling=false
fi
if [ -z "$keep_gVCF" ]; then
	export keep_gVCF=false
fi
if [ -z "$filter_ExcHet" ]; then
  filter_ExcHet=false
fi
if [ -z "$filtered_vcf" ]; then
  filtered_vcf=true
fi
if [ -z "$genomecov_est" ]; then
  genomecov_est=false
fi
mkdir -p "${projdir}"/tmp
export TMPDIR="${projdir}"/tmp

if [[ "${projdir##*/}" =~ ['!@#$%^&*()_+'] ]]; then
  echo -e "${magenta}- the project directory is: ${projdir##*/}" > "${projdir}"/job_killed_error.txt
  echo -e "${magenta}- project directory name should be alphanumeric" >> job_killed_error.txt
  echo -e "${magenta}- rename project directory without special characters and resubmit job" >> "${projdir}"/job_killed_error.txt
  echo -e "${magenta}- GBSapp will quit in 5 seconds ${white}\n" >> "${projdir}"/job_killed_error.txt
  sleep 10
  exit 0
fi
wait

main () {
	cd $projdir
	cd samples
	export nfiles=$(wc -l ../$samples_list | awk '{print $1}')
	export totalk=$(awk '/^MemTotal:/{print $2}' /proc/meminfo)
	export loopthreads=2
	if [[ "$threads" -gt 1 ]]; then
	  export N=$((threads/2))
	  export ram1=$(($totalk/$N))
	else
	  export N=1 && export loopthreads=$threads
    export ram1=$totalk
	fi
  export ngmthreads=$(echo "$totalk*0.000001" | bc)
  export ngmthreads=$(printf "%.0f" $ngmthreads)
  export ngmthreads=$((ngmthreads/10))
  if [[ "$ngmthreads" -gt "$threads" ]]; then
    export ngmthreads=$threads
  fi

  export ram1=$((ram1/1000000))
	export Xmx1=-Xmx${ram1}G
	export ram2=$(echo "$totalk*0.0000009" | bc)
	export ram2=${ram2%.*}
	export Xmx2=-Xmx${ram2}G
	if [[ "$nfiles" -lt "$N" ]]; then
	  export N=$nfiles && export loopthreads=$threads
	fi

	if [[ "$threads" -le 6 ]]; then
		export prepthreads=$threads
		export Xmxp=$Xmx2
		export prepN=1
	else
		export prepthreads=6
		export prepN=$(( threads / prepthreads ))
		export ramprep=$(( ram2 / prepN ))
		export Xmxp=-Xmx${ramprep}G
	fi

	if [[ "$threads" -le 4 ]]; then
		export gthreads=$threads
		export Xmxg=$Xmx2
		export gN=1
	else
    export ramg=20
		export Xmxg=-Xmx${ramg}G
		export gN=$(($ram2/$ramg))
    if [[ "$gN" -lt 1 ]]; then export gN=1; fi
		export gthreads=$(($threads /$gN ))
    if [[ "$gN" -gt "$((threads/4))" ]]; then
      gN=$((threads/4))
    fi
    if [[ "$gthreads" -lt 4 ]]; then
      export gthreads=4
      export gN=$(($threads / $gthreads ))
      export ramg=$(( $ram2 / $gN ))
      export Xmxg=-Xmx${ramg}G
    fi
	fi
}
cd $projdir
main &>> log.out


cd $projdir
echo -e "${blue}\n############################################################################## ${yellow}\n- Index Reference Genome \n${blue}##############################################################################${white}\n"
main () {
cd $projdir
cd refgenomes
for f in $(ls *.FASTA 2> /dev/null); do if [[ "$f" == *.FASTA ]]; then echo -e "${magenta}- changing upper case FASTA extension of reference genome(s) to lower case"; mv $f ${f%.FASTA}.fasta  2> /dev/null; fi; done
for f in $(ls *.FA 2> /dev/null); do if [[ "$f" == *.FA ]]; then echo -e "${magenta}- changing upper case FA extension of reference genome(s) to lower case"; mv $f ${f%.FA}.fa  2> /dev/null; fi; done
for f in $(ls *.FNA 2> /dev/null); do if [[ "$f" == *.FNA ]]; then echo -e "${magenta}- changing upper case FNA extension of reference genome(s) to lower case"; mv $f ${f%.FNA}.fna  2> /dev/null; fi; done
for f in $(ls *.FAS 2> /dev/null); do if [[ "$f" == *.FAS ]]; then echo -e "${magenta}- changing upper case FAS extension of reference genome(s) to lower case"; mv $f ${f%.FAS}.fas  2> /dev/null; fi; done
for f in $(ls *.FFN 2> /dev/null); do if [[ "$f" == *.FFN ]]; then echo -e "${magenta}- changing upper case FFN extension of reference genome(s) to lower case"; mv $f ${f%.FFN}.ffn  2> /dev/null; fi; done
for f in $(ls *.FAA 2> /dev/null); do if [[ "$f" == *.FAA ]]; then echo -e "${magenta}- changing upper case FAA extension of reference genome(s) to lower case"; mv $f ${f%.FAA}.faa  2> /dev/null; fi; done
wait
for i in *.gz; do
	sleep $((RANDOM % 2))
  gunzip $i >/dev/null 2>&1
done

cd refgenomes
mv ${ref1%.f*}.nohardmasked.fasta $ref1
mv $ref1 ../
rm -rf ./*
mv ../$ref1 ./

if ls ./*.dict 1> /dev/null 2>&1; then
	:
else
	checknfastafiles=$(ls *.f* | grep -v .fai | grep -v .ngm | grep -v _original.fasta | wc -l)
	if [[ $checknfastafiles -gt 1 ]]; then
		echo -e "${magenta}- expecting only 1 fasta file for reference genome ${white}\n"
		echo -e "${magenta}- GBSapp will quit in 5 seconds ${white}\n"
		sleep 5 && exit 1
	fi
	if [ -z "$ref1" ]; then
		for ref in *.f*; do
			sleep $((RANDOM % 2))
      ref1=${ref%.fa*}.fasta
		done
	fi
	export ncontigscaffold=$(grep '>' $ref1 | wc -l)
	if [[ $ncontigscaffold -gt $max_pseudoMol ]]; then
		nfakechr=$((threads/2))
		cat $ref1 | awk '/^>/{if(N>0) printf("\n"); ++N; printf("%s\t",$0);next;} {printf("%s",$0);}END{printf("\n");}' > panref0.txt
		awk 'BEGIN{srand() }
		{ lines[++d]=$0 }
			END{
			while (1){
				if (e==d) {break}
					RANDOM = int(1 + rand() * d)
					if ( RANDOM in lines  ){
					print lines[RANDOM]
					delete lines[RANDOM]
					++e
				}
			}
		}' panref0.txt > panref0.fasta
		flength=$(wc -l panref0.fasta | awk '{print $1}'); nsplit=$(( flength / nfakechr ))
		split -a 2 -d -l $nsplit panref0.fasta Chr
		rm panref0*
		for i in Chr*; do
			sleep $((RANDOM % 2))
      Nstitch=$(printf "N%.0s" $(seq 24))

			awk -v Nstitch=$Nstitch '{print $1"\t"$2}' $i | awk '{gsub(/\t/,"\n"); print $0}' > ${i}.fasta
			grep -v '>'  ${i}.fasta |  awk '/^>/{if(N>0) printf("\n"); ++N; printf("%s\t",$0);next;} {printf("%s",$0);}END{printf("\n");}' | fold -w 100 > ${i}.txt
			#generate index of concatenated contigs/scaffold
			awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' ${i}.fasta | awk '{gsub(/>/,""); print $0}' |\
			awk '{ ORS = NR % 2 ? "\t" : "\n" } 1' | awk '{print $1"\t"$2+100}' | awk -F "\t" '{sum+=$2;print($1,"\t",sum-$2);}' |\
			awk -v chrid=${i} -v refgenome=${ref1%.f*} '{print refgenome"_"chrid"\t"$1"\t"$2}' >> contigscaffold_index.txt
		done
		cp $ref1 ${ref1%.f*}_original.fasta
		cat /dev/null > $ref1
		for filename in Chr*.txt; do
			sleep $((RANDOM % 2))
      echo ">""${filename%.txt}" >> $ref1
			cat "$filename" >> $ref1
		done
		rm Chr*
	fi
  wait


	mkdir split
	awk '/^>/{close("file.fasta"f);f++}{print $0 > "./split/file"f}' "${ref1}" & PID=$!
	wait $PID
	cd split
	export checksplit=$(wc -c file* | grep -v 'total' | awk '($1 > 500000000 )' | wc -l)
	if [[ "$checksplit" -gt 0 ]]; then
		for i in file*; do
			sleep $((RANDOM % 2))
      name=$(awk 'NR==1{print $1}' "$i" | awk '{gsub(/>/,""); print}')
			awk 'NR>1{print $0}' $i > ${name}.fasta
			count=$(wc -c ${name}.fasta | awk '{print $1}')
			tr -d '\n' < "${name}".fasta | fold -w 100 |\
			split -d -l 5000000  & PID=$!
		  wait $PID
			for outfile in x*; do
				sleep $((RANDOM % 2))
        awk -v name="$name" -v outfile="$outfile" -i inplace 'BEGINFILE{print ">"name"_"outfile}{print}' $outfile
				mv $outfile ${name}_${outfile}.txt
			done
		done
		cd ../
		sleep $((RANDOM % 2))
    mv $ref1 ./old_"${ref1%.f*}_fasta.txt"
		cat ./split/*Chr*.txt > $ref1
		wait
		rm -rf split
	else
		cd ${projdir}/refgenomes
		rm -rf split
	fi
  wait
fi

cd $projdir
cd refgenomes
if [[ "$aligner" == "ngm" ]]; then
  if [[ $(ls ./*.ngm 1> /dev/null 2>&1 | wc -l) -gt 0 ]]; then
  	echo -e "${magenta}- indexed genome available ${white}\n"
  	if [ -z "$ref1" ]; then
  		for ref in *.f*; do
  			sleep $((RANDOM % 2))
        ref1=${ref%%.f*}.fasta
  		done
  	fi
  else
  	echo -e "${magenta}- indexing single reference subgenome ${white}\n"
  	if [ -z "$ref1" ]; then
  		for ref in *.f*; do
  			sleep $((RANDOM % 2))
        ref1=${ref%%.f*}.fasta
  		done
  	fi
  	awk '{ sub("\r$",""); print}' $ref1 | awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' | \
    awk '/>/{gsub(/a$/,"A");gsub(/b$/,"B");gsub(/c$/,"C");gsub(/d$/,"D");gsub(/e$/,"E");gsub(/f$/,"F");gsub(/g$/,"G");gsub(/h$/,"H");}1' > ref.txt
  	n=">${ref1%.f*}_"
  	awk '{ sub("\r$",""); print}' ref.txt | awk -v n="$n" '{gsub(n,">"); print}' | awk -v n="$n" '{gsub(/>/,n); print}' > $ref1
  	rm ref.txt
  	$samtools faidx $ref1
  	$java -jar $picard CreateSequenceDictionary REFERENCE=$ref1 OUTPUT=${ref1%.f*}.dict
    $ngm -r $ref1
    $genmap index -F $ref1 -I ./genmap_out >/dev/null 2>&1
    $genmap map -K 100 -E 2 -I ./genmap_out -O ./genmap_out/ -t -w -bg --threads $threads >/dev/null 2>&1
    if [[ "$paralogs" == "false" ]]; then
      threshold=1
    else
      threshold=8
    fi
    python "$wig2bed" "./genmap_out/${ref1%.f*}.genmap.wig" "$threshold" "./genmap_out/lowmap_merged.bed" > /dev/null 2>&1
    $bedtools maskfasta -fi $ref1 -bed ./genmap_out/lowmap_merged.bed -fo ${ref1%.f*}.hardmasked.fasta >/dev/null 2>&1
    mv $ref1 ${ref1%.f*}.nohardmasked.fasta
    mv ${ref1%.f*}.hardmasked.fasta $ref1
    $ngm -r $ref1
  fi
fi
wait
if [[ "$aligner" == "minimap2" ]]; then
  if [[ $(ls ./*.mmi 1> /dev/null 2>&1 | wc -l) -gt 0 ]]; then
  	echo -e "${magenta}- indexed genome available ${white}\n"
  	if [ -z "$ref1" ]; then
  		for ref in *.f*; do
  			sleep $((RANDOM % 2))
        ref1=${ref%%.f*}.fasta
  		done
  	fi
  else
  	echo -e "${magenta}- indexing single reference subgenome ${white}\n"
  	if [ -z "$ref1" ]; then
  		for ref in *.f*; do
  			sleep $((RANDOM % 2))
        ref1=${ref%%.f*}.fasta
  		done
  	fi
  	awk '{ sub("\r$",""); print}' $ref1 | awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' | \
    awk '/>/{gsub(/a$/,"A");gsub(/b$/,"B");gsub(/c$/,"C");gsub(/d$/,"D");gsub(/e$/,"E");gsub(/f$/,"F");gsub(/g$/,"G");gsub(/h$/,"H");}1' > ref.txt
  	n=">${ref1%.f*}_"
  	awk '{ sub("\r$",""); print}' ref.txt | awk -v n="$n" '{gsub(n,">"); print}' | awk -v n="$n" '{gsub(/>/,n); print}' > $ref1
  	rm ref.txt

    $samtools faidx $ref1
    $java -jar $picard CreateSequenceDictionary REFERENCE=$ref1 OUTPUT=${ref1%.f*}.dict
    $minimap2 -d ${ref1%.f*}.mmi $ref1
    $genmap index -F $ref1 -I ./genmap_out >/dev/null 2>&1
    $genmap map -K 100 -E 2 -I ./genmap_out -O ./genmap_out/ -t -w -bg --threads $threads >/dev/null 2>&1
    if [[ "$paralogs" == "false" ]]; then
      threshold=1
    else
      threshold=8
    fi
    python "$wig2bed" "./genmap_out/${ref1%.f*}.genmap.wig" "$threshold" "./genmap_out/lowmap_merged.bed" > /dev/null 2>&1
    $bedtools maskfasta -fi $ref1 -bed ./genmap_out/lowmap_merged.bed -fo ${ref1%.f*}.hardmasked.fasta >/dev/null 2>&1
    mv $ref1 ${ref1%.f*}.nohardmasked.fasta
    mv ${ref1%.f*}.hardmasked.fasta $ref1
    $minimap2 -d ${ref1%.f*}.mmi $ref1
  fi
fi
wait
if [[ "$RNA" == "true" ]]; then
  mkdir -p star_index
  if [[ "$(ls ./star_index/* 1> /dev/null 2>&1 | wc -l)" -eq 0 ]]; then
    gff=${ref1%.fasta}.gff*
    "$star" --runThreadN $threads \
         --runMode genomeGenerate \
         --genomeDir star_index \
         --genomeFastaFiles $ref1 \
         --sjdbGTFfile $gff \
         --sjdbOverhang 136
  fi
fi

declare -a arr=("${ref1%.f*}.dict" "${ref1}" "${ref1}.fai")
for file in "${arr[@]}"; do
  sleep $((RANDOM % 2))
	if test -f $file; then
		:
	else
		echo -e "${magenta}- reference genome not indexed: missing $file ${white}\n"
		echo -e "${magenta}- make sure reference genome is indexed ${white}\n"
		echo -e "${magenta}- GBSapp will quit in 5 seconds ${white}\n"
		sleep 5 && exit 1
	fi
done
}
cd $projdir
if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
	time main 2>> ${projdir}/log.out
fi


echo -e "${blue}\n############################################################################## ${yellow}\n- Organizing sample fastq files \n${blue}##############################################################################${white}\n"
main () {
	cd $projdir
	cd samples
  if [[ "$(ls -A *.f* 2> /dev/null | wc -l)" -eq 0 ]] || [[ "$(ls -A ./preprocess/alignment/*sam* 2> /dev/null | wc -l)" -gt 0 ]] ; then
    if [[ ! -d "pe" ]] && [[ ! -d "se" ]]; then
      :> filename_reformatted.txt
      :> flushed_reads.txt
    fi
  fi
	if test ! -f filename_reformatted.txt; then
		if [ -d "se" ]; then
			:
		else
			mkdir -p se
		fi
		if [ -d "pe" ]; then
			:
		else
			mkdir -p pe
		fi
		cd ${projdir}/samples/se
		if [ -z "$(ls -A ../pe 2> /dev/null)" ]; then
			if [ -z "$(ls -A ../se 2> /dev/null)" ]; then
				cd ../
				for i in *.f*; do (
          sleep $((RANDOM % 2))
          if [[ "$i" == *"R2.f"* ]]; then
            :
          else
  					if [[ "$i" == *.R1* ]]; then
  						mv $i ${i/.R1/}
  					elif [[ "$i" == *_R1* ]]; then
  						mv $i ${i/_R1/}
  					else
  						:
  					fi
          fi ) &
          if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
          	wait
          fi
				done
        wait
			fi
		fi

		cd ${projdir}/samples/se
		if [ "$(ls -A ../se 2> /dev/null)" ]; then
      if [ -z "$(ls -A ../pe 2> /dev/null)" ]; then
				echo -e "${magenta}- only single-end reads available in se-folder ${white}\n"
				for i in *.f*; do (
					sleep $((RANDOM % 2))
          if [[ "$i" == *".R1"* ]]; then
						mv "$i" ../${i/.R1/} && rm "$i" 2> /dev/null
					elif [[ "$i" == *_R1* ]]; then
						mv "$i" ../${i/_R1/} && rm "$i" 2> /dev/null
					else
						mv "$i" ../$i && rm "$i" 2> /dev/null
					fi ) &
          if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
          	wait
          fi
				done
        wait
			fi
		fi

		cd ${projdir}/samples/pe
		if [ "$(ls -A ../pe 2> /dev/null)" ]; then
      if [ -z "$(ls -A ../se 2> /dev/null)" ]; then
				echo -e "${magenta}- only paired-end reads available in pe-folder ${white}\n"
				for i in *R1.f*; do (
					sleep $((RANDOM % 2))
          mv "${i%R1.f*}"R2.f* ../ 2> /dev/null && rm ${i%R1.f*}R2.f* 2> /dev/null
					if [[ "$i" == *.R1* ]]; then
						mv "$i" ../${i/.R1/} 2> /dev/null && rm "$i" 2> /dev/null
					elif [[ "$i" == *_R1* ]]; then
						mv "$i" ../${i/_R1/} 2> /dev/null && rm "$i" 2> /dev/null
					else
						echo -e "${magenta}- check paired-end filenames for proper filename format, i.e. .R1 or _R1 and .R2 or _R2  ${white}\n"
						echo -e "${magenta}- Do you want to continue running GBSapp? ${white}\n"
						read -p "- y(YES) or n(NO) " -n 1 -r
						if [[ ! $REPLY =~ ^[Yy]$ ]]; then
							printf '\n'
							exit 1
						fi
					fi ) &
          if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
          	wait
          fi
				done
        wait
			fi
		fi

		cd ${projdir}/samples/pe
		if [ "$(ls -A ../se 2> /dev/null)" ]; then
			if [ "$(ls -A ../pe 2> /dev/null)" ]; then
				for i in *R1.f*; do (
					sleep $((RANDOM % 2))
          mv ${i%R1.f*}R2.f* ../ &&
					if [[ "$i" == *.R1* ]]; then
						cat $i ../se/${i%.R1.f*}.R2.f* ../se/${i%.R1.f*}.R1.f* > ../${i} 2> /dev/null && rm "$i" ../se/${i%.R1.f*}.R1.f* ../se/${i%.R1.f*}.R2.f* 2> /dev/null
						mv ../"${i}" ../${i/.R1/} 2> /dev/null && rm ../${i} 2> /dev/null
					elif [[ "$i" == *_R1* ]]; then
						cat "$i" ../se/${i%_R1.f*}_R1.f* ../se/${i%_R1.f*}_R2.f* > ../${i} 2> /dev/null && rm "$i" ../se/${i%.R1.f*}_R1.f* ../se/${i%.R1.f*}_R2.f* 2> /dev/null
						mv ../"${i}" ../${i/_R1/} 2> /dev/null && rm ../"${i}" 2> /dev/null
					else
						echo -e "${magenta}- check paired-end filenames for proper filename format, i.e. .R1 or _R1 and .R2 or _R2 ${white}\n"
						echo -e "${magenta}- Do you want to continue running GBSapp? ${white}\n"
						read -p "- y(YES) or n(NO) " -n 1 -r
						if [[ ! $REPLY =~ ^[Yy]$ ]]; then
							printf '\n'
							exit 1
						fi
					fi ) &
          if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
          	wait
          fi
				done
        wait
			fi
		fi

		cd ${projdir}/samples/
    find . -type d -empty -delete
		sampno=$(ls -1 | wc -l)
		if [[ "$sampno" == "0" ]]; then
			echo -e "${magenta}- \n- samples folder is empty, exiting pipeline ${white}\n"
			exit 1
    fi
		echo filename_reformatted > filename_reformatted.txt
	fi

	if test ! -f flushed_reads.txt; then
		if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
			for i in *.f*; do (
        if [[ "$i" == *"_uniq.fasta"* || "$i" == *"fq.gz"* ]]; then
          :
        else
  				if [[ $(file $i 2> /dev/null) =~ gzip ]]; then
  					fa_fq=$(zcat $projdir/samples/$i 2> /dev/null | head -n1 | cut -c1-1)
  				else
  					fa_fq=$(cat $projdir/samples/$i | head -n1 | cut -c1-1)
  				fi

  				if [[ $(file $i 2> /dev/null) =~ gzip ]]; then
  					if [[ "${fa_fq}" == "@" ]]; then
              awk 'NR%4 == 1 {print} NR%4 == 2 {print}' <(zcat $i) | \
              awk '{gsub(/AAAAAAAAAA/,"\t"); gsub(/CCCCCCCCCC/,"\t"); gsub(/GGGGGGGGGG/,"\t"); gsub(/TTTTTTTTTT/,"\t"); print $1}' | \
              awk '/^@/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
              awk '{print $1"\t"$2"\t+\t"$2}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$4); print}' | awk '{gsub(/\t/,"\n")}1' | gzip > nopolymers_${i} &&
              mv nopolymers_${i} $i &&
  						awk 'NR%2==0' <(zcat $i) | awk 'NR%2==1' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=1000; i++){x=int(rand()*NR) + 1; print a[x];}}' > ${i}_length_distribution.txt
  					fi
  					if [[ "${fa_fq}" == ">" ]]; then
              awk '{gsub(/AAAAAAAAAA/,"\t"); gsub(/CCCCCCCCCC/,"\t"); gsub(/GGGGGGGGGG/,"\t"); gsub(/TTTTTTTTTT/,"\t"); print $1}' <(zcat $i) | gzip > nopolymers_${i} &&
              mv nopolymers_${i} $i &&
  						awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' <(zcat $i) | \
  						awk 'NR%2==0' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=1000; i++){x=int(rand()*NR) + 1; print a[x];}}' > ${i}_length_distribution.txt
  					fi
  				else
  					if [[ "${fa_fq}" == "@" ]]; then
              awk 'NR%4 == 1 {print} NR%4 == 2 {print}' $i | \
              awk '{gsub(/AAAAAAAAAA/,"\t"); gsub(/CCCCCCCCCC/,"\t"); gsub(/GGGGGGGGGG/,"\t"); gsub(/TTTTTTTTTT/,"\t"); print $1}' | \
              awk '/^@/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
              awk '{print $1"\t"$2"\t+\t"$2}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$4); print}' | awk '{gsub(/\t/,"\n")}1' > nopolymers_${i} &&
              mv nopolymers_${i} $i &&
  						awk 'NR%2==0' $i | awk 'NR%2==1' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=1000; i++){x=int(rand()*NR) + 1; print a[x];}}' > ${i}_length_distribution.txt
  					fi
  					if [[ "${fa_fq}" == ">" ]]; then
              awk '{gsub(/AAAAAAAAAA/,"\t"); gsub(/CCCCCCCCCC/,"\t"); gsub(/GGGGGGGGGG/,"\t"); gsub(/TTTTTTTTTT/,"\t"); print $1}' $i > nopolymers_${i} &&
              mv nopolymers_${i} $i &&
  						awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' $i | \
  						awk 'NR%2==0' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=1000; i++){x=int(rand()*NR) + 1; print a[x];}}' > ${i}_length_distribution.txt
  					fi
  				fi
        fi ) &
        if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
          wait
        fi
			done
      wait

      :> length_distribution.txt
      for lenfile in $(find . -name "*_length_distribution.txt"); do
        cat $lenfile >> length_distribution.txt &&
        wait
      done
  		wait
      rm *_length_distribution.txt 2> /dev/null
      awk '{print length($0)}' length_distribution.txt | awk '$1>=64' | sort -n > tmp.txt && :> length_distribution.txt &&
      mv tmp.txt length_distribution.txt &&
			export max_seqread_len=$(awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}' length_distribution.txt)
			rm length_distribution.txt
      wait

      for i in *.f*; do (
        if [[ "$i" == *"_uniq.fasta"* || "$i" == *"fq.gz"* ]]; then
          :
        else
  				if [[ $(file $i 2> /dev/null) =~ gzip ]]; then
  					fa_fq=$(zcat $projdir/samples/$i 2> /dev/null | head -n1 | cut -c1-1)
  				else
  					fa_fq=$(cat $projdir/samples/$i | head -n1 | cut -c1-1)
  				fi

  				if [[ $(file $i 2> /dev/null) =~ gzip ]]; then
  					if [[ "${fa_fq}" == "@" ]]; then
              awk 'NR%2==0' <(zcat $i) | awk 'NR%2==1' | awk -v max=$max_seqread_len '{print substr($0,1,max)}' | awk 'NF' | \
              awk '{print "@seq"NR"\t"$1"\t"$1}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | \
              awk '{print $1"\n"$2"\n+\n"$3}' | gzip > "${i%.f*}"_tmp.fa.gz &&
              mv "${i%.f*}_tmp.fa.gz" "${i%.f*}.fastq.gz" &&
              wait
  					fi
  					if [[ "${fa_fq}" == ">" ]]; then
              grep -v '^>' <(zcat $i) | awk -v max=$max_seqread_len '{print substr($0,1,max)}' | awk 'NF' | \
              awk '{print "@seq"NR"\n"$1}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | \
              awk '{print $1"\n"$2"\n+\n"$3}' | gzip > "${i%.f*}"_tmp.fa.gz &&
              mv "${i%.f*}_tmp.fa.gz" "${i%.f*}.fastq.gz" &&
              wait
  					fi
  				else
  					if [[ "${fa_fq}" == "@" ]]; then
              awk 'NR%2==0' $i | awk 'NR%2==1' | awk -v max=$max_seqread_len '{print substr($0,1,max)}' | awk 'NF' | \
  						awk '{print "@seq"NR"\t"$1"\t"$1}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | \
              awk '{print $1"\n"$2"\n+\n"$3}' | gzip > "${i%.f*}"_tmp.fa.gz &&
              mv "${i%.f*}_tmp.fa.gz" "${i%.f*}.fastq.gz" &&
              wait
  					fi
  					if [[ "${fa_fq}" == ">" ]]; then
              grep -v '^>' $i | awk -v max=$max_seqread_len '{print substr($0,1,max)}' | awk 'NF' | \
              awk '{print "@seq"NR"\n"$1}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | \
              awk '{print $1"\n"$2"\n+\n"$3}' | gzip > "${i%.f*}"_tmp.fa.gz &&
              mv "${i%.f*}_tmp.fa.gz" "${i%.f*}.fastq.gz" &&
              wait
  					fi
  				fi
  		  fi ) &
        if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
          wait
        fi
      done
      wait
		fi
		find . -type d -empty -delete
		printf "Improvement in flushed reads already implemented""\n" > flushed_reads.txt
	fi
	if [[ "$lib_type" =~ "WGS" ]] || [[ "$lib_type" =~ "wgs" ]]; then
    if [[ "$subsample_WGS_in_silico_qRRS" == false ]]; then
      printf "Improvement in flushed reads not required for shotgun WGS data""\n" > flushed_reads.txt
    fi
    if test ! -f flushed_reads.txt || [[ "$(ls ${projdir}/preprocess/alignment/*_redun.sam.gz 2> /dev/null | wc -l)" -eq 0 ]]; then
      if [[ "$subsample_WGS_in_silico_qRRS" == true ]]; then
        for i in *.f*; do (
        if [[ "$i" == *"_uniq.fasta"* || "$i" == *"_R1_uniq.fasta"* || "$i" == *"_R2_uniq.fasta"* || "$i" == *"_uniq.hold.fasta"* || "$i" == *"_R1_uniq.hold.fasta"* || "$i" == *"_R2_uniq.hold.fasta"* || "$i" == *"fq.gz"* ]]; then
          :
        else
            if [[ $(file $i 2> /dev/null) =~ gzip ]]; then
              fa_fq=$(zcat ${projdir}/samples/$i 2> /dev/null | head -n1 | cut -c1-1)
            else
              fa_fq=$(cat ${projdir}/samples/$i | head -n1 | cut -c1-1)
            fi
            wait

            if [[ $(file $i 2> /dev/null) =~ gzip ]]; then
              if [[ "${fa_fq}" == "@" ]]; then
                awk 'NR%2==0' <(zcat $i) | awk 'NR%2==1' | awk '{gsub(/ATGCAT/,"ATGCAT\nATGCAT");}1' | awk 'length >= 64 && length <= 600' | \
                grep '^ATGCAT.*ATGCAT$' | awk '{print ">frag"NR"\n"$0}' | $gzip > ${i%.f*}_RE1.fasta.gz &&
                awk 'NR%2==0' <(zcat $i) | awk 'NR%2==1' | awk '{gsub(/CATG/,"CATG\nCATG");}1' | awk 'length >= 64 && length <= 600' | \
                grep '^CATG.*CATG$' | awk '{print ">frag"NR"\n"$0}' | $gzip > ${i%.f*}_RE2.fasta.gz &&
                awk 'NR%2==0' <(zcat $i) | awk 'NR%2==1' | awk '{gsub(/ATGCAT/,"ATGCAT\nATGCAT");gsub(/CATG/,"CATG\nCATG");}1' | awk 'length >= 64 && length <= 600' | \
                grep '^ATGCAT.*CATG$\|^CATG.*ATGCAT$' | awk '{print ">frag"NR"\n"$0}' | $gzip > ${i%.f*}_RE1RE2.fasta.gz &&
                cat ${i%.f*}_RE1.fasta.gz ${i%.f*}_RE2.fasta.gz ${i%.f*}_RE1RE2.fasta.gz > ${i%.f*}.fasta.gz &&
                rm "${i%.f*}"_RE1.fasta.gz "${i%.f*}"_RE2.fasta.gz "${i%.f*}"_RE1RE2.fasta.gz &&
                wait
              fi
              if [[ "${fa_fq}" == ">" ]]; then
                awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' <(zcat $i) | awk 'NR%2==0' | awk '{gsub(/ATGCAT/,"ATGCAT\nATGCAT");}1' | \
                awk 'length >= 64 && length <= 600' | grep '^ATGCAT.*ATGCAT$' | awk '{print ">frag"NR"\n"$0}' | $gzip > ${i%.f*}.tmp1.gz &&
                awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' <(zcat $i) | awk 'NR%2==0' | awk '{gsub(/CATG/,"CATG\nCATG");}1' | \
                awk 'length >= 64 && length <= 600' | grep '^CATG.*CATG$' | awk '{print ">frag"NR"\n"$0}' | $gzip > ${i%.f*}.tmp2.gz &&
                awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' <(zcat $i) | awk 'NR%2==0' | awk '{gsub(/ATGCAT/,"ATGCAT\nATGCAT");gsub(/CATG/,"CATG\nCATG");}1' | \
                awk 'length >= 64 && length <= 600' | grep '^ATGCAT.*CATG$\|^CATG.*ATGCAT$' | awk '{print ">frag"NR"\n"$0}' | $gzip > ${i%.f*}.tmp3.gz &&
                cat ${i%.f*}.tmp1.gz ${i%.f*}.tmp2.gz ${i%.f*}.tmp3.gz > ${i%.f*}.tmp.gz &&
                rm "${i%.f*}".fasta.gz &&
                mv "${i%.f*}".tmp.gz "${i%.f*}".fasta.gz &&
                rm "${i%.f*}".tmp1.gz "${i%.f*}".tmp2.gz "${i%.f*}".tmp3.gz "${i%.f*}".tmp.gz &&
                wait
              fi
            else
              if [[ "${fa_fq}" == "@" ]]; then
                awk 'NR%2==0' $i | awk 'NR%2==1' | awk '{gsub(/ATGCAT/,"ATGCAT\nATGCAT")}1' | awk 'length >= 64 && length <= 600' | \
                grep '^ATGCAT.*ATGCAT$' | awk '{print ">frag"NR"\n"$0}' | $gzip > ${i%.f*}_RE1.fasta.gz &&
                awk 'NR%2==0' $i | awk 'NR%2==1' | awk '{gsub(/CATG/,"CATG\nCATG");}1' | awk 'length >= 64 && length <= 600' | \
                grep '^CATG.*CATG$' | awk '{print ">frag"NR"\n"$0}' | $gzip > ${i%.f*}_RE2.fasta.gz &&
                awk 'NR%2==0' $i | awk 'NR%2==1' | awk '{gsub(/ATGCAT/,"ATGCAT\nATGCAT");gsub(/CATG/,"CATG\nCATG");}1' | awk 'length >= 64 && length <= 600' | \
                grep '^ATGCAT.*CATG$\|^CATG.*ATGCAT$' | awk '{print ">frag"NR"\n"$0}' | $gzip > ${i%.f*}_RE1RE2.fasta.gz &&
                cat "${i%.f*}"_RE1.fasta.gz "${i%.f*}"_RE2.fasta.gz "${i%.f*}"_RE1RE2.fasta.gz > "${i%.f*}".fasta.gz &&
                rm "${i%.f*}"_RE1.fasta.gz "${i%.f*}"_RE2.fasta.gz "${i%.f*}"_RE1RE2.fasta.gz &&
                wait
              fi
              if [[ "${fa_fq}" == ">" ]]; then
                awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $i | awk 'NR%2==0' | awk '{gsub(/ATGCAT/,"ATGCAT\nATGCAT");}1' | \
                awk 'length >= 64 && length <= 600' | grep '^ATGCAT.*ATGCAT$\|^CATG.*CATG$' | awk '{print ">frag"NR"\n"$0}' | $gzip > ${i%.f*}.tmp1.gz &&
                awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $i | awk 'NR%2==0' | awk '{gsub(/CATG/,"CATG\nCATG");}1' | \
                awk 'length >= 64 && length <= 600' | grep '^CATG.*CATG$' | awk '{print ">frag"NR"\n"$0}' | $gzip > ${i%.f*}.tmp2.gz &&
                awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $i | awk 'NR%2==0' | awk '{gsub(/ATGCAT/,"ATGCAT\nATGCAT");gsub(/CATG/,"CATG\nCATG");}1' | \
                awk 'length >= 64 && length <= 600' | grep '^ATGCAT.*CATG$\|^CATG.*ATGCAT$' | awk '{print ">frag"NR"\n"$0}' | $gzip > ${i%.f*}.tmp3.gz &&
                cat ${i%.f*}.tmp1.gz ${i%.f*}.tmp2.gz ${i%.f*}.tmp3.gz > ${i%.f*}.tmp.gz &&
                rm "${i%.f*}".fasta.gz &&
                mv "${i%.f*}".tmp.gz "${i%.f*}".fasta.gz &&
                rm "${i%.f*}".tmp1.gz "${i%.f*}".tmp2.gz "${i%.f*}".tmp3.gz "${i%.f*}".tmp.gz &&
                wait
              fi
            fi
            # grep -v '^>' <(zcat ${i%.f*}.fasta.gz) | awk '{print substr($0,1,64)}' | awk '{print "@B"NR"\t"$1"\t"$1}' | \
            # awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | awk '{print $1"\n"$2"\n+\n"$3}' | gzip > "${i%.f*}"_tmp1.fq.gz &&
            # grep -v '^>' <(zcat ${i%.f*}.fasta.gz) | awk -v max=$max_seqread_len 'length == max' | awk -v max=$max_seqread_len '{print substr($0,65,max)}' | \
            # awk '{print "@E"NR"\t"$1"\t"$1}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | awk '{print $1"\n"$2"\n+\n"$3}' | gzip > "${i%.f*}"_tmp2.fq.gz &&
            # rm ${i%.f*}.fasta.gz &&
            # cat "${i%.f*}"_tmp1.fq.gz "${i%.f*}"_tmp2.fq.gz > ${i%.f*}.fastq.gz &&
            # rm "${i%.f*}"_tmp1.fq.gz "${i%.f*}"_tmp2.fq.gz &&
            wait
          fi
          ) &
          if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
            wait
          fi
        done
        wait
        printf "Improvement in flushed reads already implemented""\n" > flushed_reads.txt
      fi
    fi
	fi
}
cd $projdir
cd samples
if [ -d "pe" ]; then
	fqpass=$(find ./pe -maxdepth 1 -name '*_R1.f*' -o -name '*_R2.f*' -o -name '*.R1.f*' -o -name '*.R2.f*' | wc -l)
	fqfail=$(ls ./pe/* | wc -l)
	fqfail=$((fqfail-fqpass))
	if [[ "$fqfail" -lt 1 ]]; then
		if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
			if test ! -f flushed_reads.txt; then
        time main &>> ${projdir}/log.out
      fi
		fi
	else
		echo -e "${magenta}- samples' PE fastq filenames requires formatting (i.e. needs to end in "_R1.fastq" or ".R1.fastq" and "_R2.fastq" or ".R2.fastq") ${white}\n"
		sleep 5 && exit 1
	fi
else
	if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
		if test ! -f flushed_reads.txt; then
      time main &>> ${projdir}/log.out
    fi
	fi
fi


main () {
	cd ${projdir}
	mkdir -p samples/tmp
  mkdir -p preprocess/tmp
  mkdir -p preprocess/alignment
	mkdir -p snpcall/tmp
	mkdir -p alignment_summaries
	mkdir -p ./alignment_summaries/copy_number

	if [[ "$samples_list" == "samples_list_node_1.txt" ]] && [[ "$(ls -A ./samples/*.f* 2> /dev/null | wc -l)" -gt 0 ]]; then
		for i in samples_list_node_*.txt; do
			:> ${i%.txt}_hold.txt
			while IFS="" read -r line; do
        ls -l ./samples/$line | awk '{print $5"\t"$9}' >> ${i%.txt}_hold.txt
			done < <(grep -v '_tmp.fa' $i | grep -v _R2.f | grep -v _uniq.fasta | grep -v _uniq_R1.fasta | grep -v _uniq_R2.fasta | grep -v _uniq.hold.fasta | grep -v _uniq_R1.hold.fasta | grep -v _uniq_R2.hold.fasta | grep -v fq.gz)
			sort -nr -k1 ${i%.txt}_hold.txt | awk '{gsub(/.\/samples\//,""); print $2}' > $i
			rm "${i%.txt}"_hold.txt
		done
	fi
}
if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
	time main &>> ${projdir}/log.out
fi
if [[ "$(cat ${projdir}/samples_list_node_* | wc -l)" -lt 5 ]]; then
  export p1=""
  export p2=""
fi
if [[ "$(cat ${projdir}/samples_list_node_* | wc -l)" -lt 50 ]] && [[ -z $maf ]]; then
  maf=0
fi


######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${yellow}\n- Performing Read Alignments & Alignment Post-Processing\n${blue}##############################################################################${white}\n"

cd $projdir
if [[ $nodes -gt 1 ]] && [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
	rm -rf /tmp/"${samples_list%.txt}" 2> /dev/null &&
	mkdir -p /tmp/${samples_list%.txt}/refgenomes /tmp/${samples_list%.txt}/samples &&
  mkdir -p /tmp/${samples_list%.txt}/preprocess/tmp /tmp/${samples_list%.txt}/preprocess/alignment /tmp/${samples_list%.txt}/snpcall/tmp &&
	cp -r ${projdir}/refgenomes/* /tmp/${samples_list%.txt}/refgenomes/ &&
  wait
	while IFS="" read -r i || [ -n "$i" ]; do
		sleep $((RANDOM % 2))
    cp ${i%.f*}_R1_uniq.fasta.gz /tmp/${samples_list%.txt}/samples/ 2> /dev/null &&
		cp ${projdir}/preprocess/alignment/${i%.f*}_redun.sam.gz /tmp/${samples_list%.txt}/preprocess/alignment/ 2> /dev/null &&
		cp ${projdir}/preprocess/${i%.f*}_*_precall.bam* /tmp/${samples_list%.txt}/preprocess/ 2> /dev/null &&
		wait
	done < <(cat ${projdir}/${samples_list})
fi

main () {

	cd $projdir
	cd samples

	if [[ "$samples_list" == "samples_list_node_1.txt" ]] && test ! -f ${projdir}/organize_files_done.txt; then
    while IFS="" read -r i || [ -n "$i" ]; do (
      sleep $((RANDOM % 2))
      if test ! -f ${projdir}/preprocess/alignment/${i%.f*}_redun.sam.gz && test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai; then
        if test ! -f ${i%.f*}_uniq.fasta.gz; then
          zcat $i | awk 'NR%2==0' | awk 'NR%2' | gzip > ${i%.f*}_R1_uniq.txt.gz 2> /dev/null &&
          wait
          if test -f ${i%.f*}_R2.fastq.gz; then
            zcat ${i%.f*}_R2.fastq.gz | awk 'NR%2==0' | awk 'NR%2' | gzip > ${i%.f*}_R2_uniq.txt.gz 2> /dev/null &&
            wait
          fi
          if test -f ${i%.f*}.R2.fastq.gz; then
            zcat ${i%.f*}.R2.fastq.gz | awk 'NR%2==0' | awk 'NR%2' | gzip > ${i%.f*}_R2_uniq.txt.gz 2> /dev/null &&
            wait
          fi
          if test ! -f ${i%.f*}_R2.fastq.gz && test ! -f ${i%.f*}.R2.fastq.gz; then
            touch ${i%.f*}_R2_uniq.txt && gzip ${i%.f*}_R2_uniq.txt &&
            wait
          fi
        fi
      fi) &
      if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
        wait
      fi
    done < <(cat ${projdir}/samples_list_node_* )
    wait

    # compress reads
    while IFS="" read -r i || [ -n "$i" ]; do (
      sleep $((RANDOM % 2))
      if test ! -f ${projdir}/preprocess/alignment/${i%.f*}_redun.sam.gz && test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai; then
        if test ! -f ${i%.f*}_uniq.fasta.gz; then
          cat ${i%.f*}_R1_uniq.txt.gz ${i%.f*}_R2_uniq.txt.gz | zcat | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{$1=$1};1' | gzip > ${i%.f*}_repeatn.txt.gz &&
          repeatn=$(zcat ${i%.f*}_repeatn.txt.gz | awk '{print $1}' | sort -n | awk '{all[NR] = $1} END{print all[int(NR*0.999 - 0.5)]}')
          zcat ${i%.f*}_repeatn.txt.gz | awk -v pat="$repeatn" '$1 < pat' | \
          awk '{gsub(/AAAAAAAAAA$/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");}1' | \
  				awk '{gsub(/CCCCCCCCCC$/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");}1' | \
  				awk '{gsub(/GGGGGGGGGG$/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");}1' | \
  				awk '{gsub(/TTTTTTTTTT$/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");}1' | \
  				awk '{gsub(/^AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");}1' | \
  				awk '{gsub(/^CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");}1' | \
  				awk '{gsub(/^GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");}1' | \
  				awk '{gsub(/^TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");}1' | \
  				awk '{gsub(/AAAAA~/,"~");gsub(/CCCCC~/,"~");gsub(/GGGGG~/,"~");gsub(/TTTTT~/,"~");gsub(/~AAAAA/,"~");gsub(/~CCCCC/,"~");gsub(/~GGGGG/,"~");gsub(/~TTTTT/,"~");}1' | \
  				awk '{gsub(/AA~/,"~");gsub(/CC~/,"~");gsub(/GG~/,"~");gsub(/TT~/,"~");gsub(/~AA/,"~");gsub(/~CC/,"~");gsub(/~GG/,"~");gsub(/~TT/,"~");gsub(/AA~/,"~");gsub(/CC~/,"~");gsub(/GG~/,"~");gsub(/TT~/,"~");gsub(/~AA/,"~");gsub(/~CC/,"~");gsub(/~GG/,"~");gsub(/~TT/,"~");}1' | \
  				awk '{gsub(/A~/,"~");gsub(/C~/,"~");gsub(/G~/,"~");gsub(/T~/,"~");gsub(/~A/,"~");gsub(/~C/,"~");gsub(/~G/,"~");gsub(/~T/,"~");gsub(/~/,"");}1' | \
          awk '{print ">seq"NR"_pe-"$1"\n"$2}' | gzip > ${i%.f*}_uniq.fasta.gz &&
          rm "${i%.f*}"_R1_uniq.txt.gz "${i%.f*}"_R2_uniq.txt.gz ${i%.f*}_repeatn.txt.gz 2> /dev/null &&
          wait
        fi
      fi ) &
      if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
        wait
      fi
    done < <(cat ${projdir}/samples_list_node_* )
    wait

		cd $projdir/samples
		find . -size 0 -delete  2> /dev/null
	fi
	wait
	touch ${projdir}/organize_files_done.txt


  if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
    touch ${projdir}/compress_done.txt
  fi

  cd ${projdir}/samples

  if [[ $nodes -gt 1 ]]; then
    if [[ "$samples_list" != "samples_list_node_1.txt" ]]; then
      rm -rf /tmp/"${samples_list%.txt}" 2> /dev/null &&
      wait
    fi
    mkdir -p /tmp/${samples_list%.txt}/refgenomes /tmp/${samples_list%.txt}/samples &&
    mkdir -p /tmp/${samples_list%.txt}/preprocess/tmp /tmp/${samples_list%.txt}/preprocess/alignment /tmp/${samples_list%.txt}/snpcall/tmp &&
    :> ${projdir}/queue_move_${samples_list%.txt}.txt &&
    queue_move=$(ls ${projdir}/queue_move_samples_list_node_* | wc -l)
    while [[ "$queue_move" -gt 1 ]]; do
      rm "${projdir}"/queue_move_${samples_list%.txt}.txt; sleep $[ ( $RANDOM % 120 )  + 30 ]s &&
      :> ${projdir}/queue_move_${samples_list%.txt}.txt &&
      queue_move=$(ls ${projdir}/queue_move_samples_list_node_* | wc -l) &&
      wait
    done
    cp -rn ${projdir}/refgenomes/* /tmp/${samples_list%.txt}/refgenomes/ &&
    wait
    while IFS="" read -r i || [ -n "$i" ]; do
      cp -rn ${i%.f*}_uniq.fasta.gz /tmp/${samples_list%.txt}/samples/ 2> /dev/null &&
      cp -rn ${projdir}/preprocess/alignment/${i%.f*}_redun.sam.gz /tmp/${samples_list%.txt}/preprocess/alignment/ 2> /dev/null &&
      cp -rn ${projdir}/preprocess/${i%.f*}_*_precall.bam* /tmp/${samples_list%.txt}/preprocess/ 2> /dev/null &&
      wait
    done < <(cat ${projdir}/${samples_list})
    wait
    rm "${projdir}"/queue_move_${samples_list%.txt}.txt &&
    wait
  fi

  # Perform read alignments
  cd ${projdir}
  if [[ "$haplome_number" ]]; then
    max_value=0
    for maxhap in ${haplome_number//,/ }; do
      if [[ "$maxhap" -gt "$max_value" ]]; then max_value=$maxhap; fi
    done
    export tophap=$((max_value*8))
  else
    export tophap=8
  fi
  if [[ $nodes -eq 1 ]]; then cd ${projdir}/samples/ ; fi
  if [[ $nodes -gt 1 ]] && test -f ${projdir}/GBSapp_run_node_1.sh; then cd /tmp/${samples_list%.txt}/samples/ ; fi
  if [[ "$RNA" == "false" ]]; then
    if [[ "$aligner" == "ngm" ]]; then
      if test ! -f ${projdir}/precall_done.txt && test ! -f ${projdir}/alignment_done; then
        while IFS="" read -r alignfq || [ -n "$alignfq" ]; do
          sleep $((RANDOM % 2))
          if test ! -f ../preprocess/alignment/${alignfq%.f*}_redun.sam.gz; then
            $ngm -r ../refgenomes/$ref1 --qry ${alignfq%.f*}_uniq.fasta.gz -o ../preprocess/alignment/${alignfq%.f*}_redun.sam -R 20 -t $threads --topn $tophap --strata 8 --affine &&
            awk '/@HD/ || /@SQ/{print}' ../preprocess/alignment/${alignfq%.f*}_redun.sam 2> /dev/null > ../preprocess/alignment/${alignfq%.f*}_redun_head.sam
            grep -v '^@' ../preprocess/alignment/${alignfq%.f*}_redun.sam 2> /dev/null | awk 'BEGIN{FS=OFS="\t"} !($10 == "*" && $6 !~ /^\*$/) {print}' | \
            cat ../preprocess/alignment/${alignfq%.f*}_redun_head.sam - | gzip  > ../preprocess/alignment/${alignfq%.f*}_redun.sam.gz &&
            rm ../preprocess/alignment/"${alignfq%.f*}"_redun.sam ../preprocess/alignment/"${alignfq%.f*}"_redun_head.sam &&
            cp -rn ../preprocess/alignment/${alignfq%.f*}_redun.sam.gz ${projdir}/preprocess/alignment/ &&
            rm "${alignfq%.f*}"_uniq.fasta.gz &&
            wait
          fi
        done < <( cat ${projdir}/${samples_list} )
      fi
    fi
    wait
    if [[ "$aligner" == "minimap2" ]]; then
      if test ! -f ${projdir}/precall_done.txt && test ! -f ${projdir}/alignment_done; then
        while IFS="" read -r alignfq || [ -n "$alignfq" ]; do
          sleep $((RANDOM % 2))
          if test ! -f ../preprocess/alignment/${alignfq%.f*}_redun.sam.gz; then
            $minimap2 -x sr --secondary=no --min-occ-floor=1000 ../refgenomes/${ref1%.f*}.mmi <(zcat ${alignfq%.f*}_uniq.fasta.gz) | awk '$10 <= 10 {print $1}' | \
            awk 'NR==FNR{keep[$1]=1;next}
                   (NR%2==1){id=substr($1,2);p=keep[id]}
                   {if(p)print}' - <(zcat ${alignfq%.f*}_uniq.fasta.gz) | \
            $minimap2 -t $threads -ax splice --secondary=no -f 0.0005 -N 8 -n 2 -m 25 ../refgenomes/${ref1%.f*}.mmi - > ${alignfq%.f*}_all.sam &&
            grep '^@' ${alignfq%.f*}_all.sam > ${alignfq%.f*}_header.sam &&
            grep -v '^@' ${alignfq%.f*}_all.sam | awk '$6 ~ /N/' | awk 'BEGIN{FS=OFS="\t"} !($10 == "*" && $6 !~ /^\*$/) {print}' > ${alignfq%.f*}_spliced_reads.sam &&
            grep -v '^@' ${alignfq%.f*}_all.sam | awk '$6 !~ /N/' | awk 'BEGIN{FS=OFS="\t"} !($10 == "*" && $6 !~ /^\*$/) {print}' > ${alignfq%.f*}_unspliced_reads.sam &&
            cat ${alignfq%.f*}_header.sam ${alignfq%.f*}_spliced_reads.sam > ${alignfq%.f*}_rna.sam &&
            cat ${alignfq%.f*}_header.sam ${alignfq%.f*}_unspliced_reads.sam > ${alignfq%.f*}_dna.sam &&
            $samtools view -bS ${alignfq%.f*}_rna.sam > ${alignfq%.f*}_rna.bam &&
            $samtools view -bS ${alignfq%.f*}_dna.sam > ${alignfq%.f*}_dna.bam &&
            $samtools fasta ${alignfq%.f*}_rna.bam > ${alignfq%.f*}_rna_reads.fasta &&
            $samtools fasta ${alignfq%.f*}_dna.bam > ${alignfq%.f*}_dna_reads.fasta &&

            $minimap2 -t $threads -ax sr --secondary=no -f 0.0005 -N 8 -n 2 -m 25 ../refgenomes/${ref1%.f*}.mmi ${alignfq%.f*}_dna_reads.fasta | $samtools sort -o ${alignfq%.f*}_dna_final.bam &&
            $samtools index ${alignfq%.f*}_dna_final.bam &&

            $minimap2 -t $threads -N 8 -ax splice -uf -k14 --secondary=no -f 0.0005 -N 8 -n 2 -m 25 ../refgenomes/${ref1%.f*}.mmi ${alignfq%.f*}_rna_reads.fasta | $samtools sort -o ${alignfq%.f*}_rna_final.bam &&

            $samtools index ${alignfq%.f*}_rna_final.bam &&
            $GATK SplitNCigarReads -R ../refgenomes/$ref1 -I ${alignfq%.f*}_rna_final.bam -O ${alignfq%.f*}_rna_final_spliced.bam &&
            $samtools merge ${alignfq%.f*}_redun.bam ${alignfq%.f*}_dna_final.bam ${alignfq%.f*}_rna_final_spliced.bam &&
            $samtools sort -o ${alignfq%.f*}_sorted.bam ${alignfq%.f*}_redun.bam &&
            $samtools view -h ${alignfq%.f*}_sorted.bam > ../preprocess/alignment/${alignfq%.f*}_redun.sam &&
            rm ${alignfq%.f*}_*sam ${alignfq%.f*}_*bam ${alignfq%.f*}_*bai ${alignfq%.f*}_rna_reads* ${alignfq%.f*}_dna_reads* &&
            awk '/@HD/ || /@SQ/{print}' ../preprocess/alignment/${alignfq%.f*}_redun.sam 2> /dev/null > ../preprocess/alignment/${alignfq%.f*}_redun_head.sam
            grep -v '^@' ../preprocess/alignment/${alignfq%.f*}_redun.sam 2> /dev/null | awk 'BEGIN{FS=OFS="\t"} !($10 == "*" && $6 !~ /^\*$/) {print}' | \
            cat ../preprocess/alignment/${alignfq%.f*}_redun_head.sam - | gzip  > ../preprocess/alignment/${alignfq%.f*}_redun.sam.gz &&
            rm ../preprocess/alignment/"${alignfq%.f*}"_redun.sam ../preprocess/alignment/"${alignfq%.f*}"_redun_head.sam &&
            cp -rn ../preprocess/alignment/${alignfq%.f*}_redun.sam.gz ${projdir}/preprocess/alignment/ &&
            rm "${alignfq%.f*}"_uniq.fasta.gz &&
            wait
          fi
        done < <( cat ${projdir}/${samples_list} )
      fi
    fi
  fi
  wait
  if [[ "$RNA" == "true" ]]; then
    mkdir -p ../preprocess/staralign
    if test ! -f ${projdir}/precall_done.txt && test ! -f ${projdir}/alignment_done; then
      while IFS="" read -r alignfq || [ -n "$alignfq" ]; do
        sleep $((RANDOM % 2))
        if test ! -f ../preprocess/alignment/${alignfq%.f*}_redun.sam.gz; then
          awk 'BEGIN {RS=">"; ORS=""} NR>1 {header=substr($0, 1, index($0, "\n") - 1); seq=substr($0, index($0, "\n") + 1); gsub(/\n/, "", seq); if (length(seq) >= 32) print ">" header "\n" seq "\n"}' <(zcat ${alignfq%.f*}_uniq.fasta.gz) | \
          awk 'BEGIN {OFS = "\n"} /^>/ {header = substr($0, 2); next} {print "@" header, $0, "+", gensub(/./, "I", "g", $0)}' | gzip > ${alignfq%.f*}_uniq.fastq.gz
          $star --runThreadN $threads \
            --genomeDir ${projdir}/refgenomes/star_index \
            --readFilesIn ${alignfq%.f*}_uniq.fastq.gz \
            --readFilesCommand zcat \
            --outFileNamePrefix ../preprocess/staralign/${alignfq%.f*}_redun_ \
            --outSAMtype SAM \
            --outFilterMultimapNmax 8 \
            --twopassMode Basic
          wait
          awk '/@HD/ || /@SQ/{print}' ../preprocess/staralign/${alignfq%.f*}_redun_Aligned.out.sam 2> /dev/null > ../preprocess/alignment/${alignfq%.f*}_redun_head.sam
          grep -v '^@' ../preprocess/staralign/${alignfq%.f*}_redun_Aligned.out.sam 2> /dev/null | awk 'BEGIN{FS=OFS="\t"} !($10 == "*" && $6 !~ /^\*$/) {print}' | \
          cat ../preprocess/alignment/${alignfq%.f*}_redun_head.sam - > ../preprocess/alignment/${alignfq%.f*}_redun.sam &&
          $samtools view -S -b ../preprocess/alignment/${alignfq%.f*}_redun.sam > ../preprocess/alignment/${alignfq%.f*}_redun.bam &&
          $GATK SplitNCigarReads -R ${projdir}/refgenomes/$ref1 -I ../preprocess/alignment/${alignfq%.f*}_redun.bam -O ../preprocess/alignment/${alignfq%.f*}_redunsplit.bam &&
          $samtools view -h ../preprocess/alignment/${alignfq%.f*}_redunsplit.bam | gzip > ../preprocess/alignment/${alignfq%.f*}_redun.sam.gz &&
          rm ../preprocess/alignment/${alignfq%.f*}_redun.sam ../preprocess/alignment/${alignfq%.f*}_redun.bam ../preprocess/alignment/${alignfq%.f*}_redunsplit* &&
          rm ../preprocess/alignment/"${alignfq%.f*}"_redun_head.sam &&
          cp -rn ../preprocess/alignment/${alignfq%.f*}_redun.sam.gz ${projdir}/preprocess/alignment/ &&
          rm ${alignfq%.f*}_uniq.fasta.gz ${alignfq%.f*}_uniq.fastq.gz &&
          rm -rf ../preprocess/staralign/*
          wait
        fi
      done < <( cat ${projdir}/${samples_list} )
    fi
  fi
  wait

  touch ${projdir}/alignment_done_${samples_list}

  cd ${projdir}/preprocess
  if [[ "$samples_list" == "samples_list_node_1.txt" ]] && test ! -f ${projdir}/alignment_done; then
    while [[ "$(ls ${projdir}/alignment_done_samples_list_node_* | wc -l)" -lt "$nodes" ]]; do
      sleep 300
    done
    touch ${projdir}/alignment_done
  fi

  while [[ ! -f ${projdir}/alignment_done ]]; do
    sleep 300
  done
  wait
  if [[ $nodes -eq 1 ]]; then cd ${projdir}/preprocess/; fi
  if [[ $nodes -gt 1 ]] && test -f ${projdir}/GBSapp_run_node_1.sh; then cd /tmp/${samples_list%.txt}/preprocess/; fi

  if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
    touch ${projdir}/compress_done.txt
  fi
  if [[ "$ploidy" -eq 1 ]]; then
    if [[ -z "$downsample_1x" ]]; then
      downsample=50
    else
      downsample="$downsample_1x"
    fi
  fi
  if [[ "$ploidy" -eq 2 ]]; then
    if [[ -z "$downsample_2x" ]]; then
      downsample=100
    else
      downsample="$downsample_2x"
    fi
  fi
  if [[ "$ploidy" -eq 4 ]]; then
    if [[ -z "$downsample_4x" ]]; then
      downsample=200
    else
      downsample="$downsample_4x"
    fi
  fi
  if [[ "$ploidy" -eq 6 ]]; then
    if [[ -z "$downsample_6x" ]]; then
      downsample=300
    else
      downsample=$downsample_6x
    fi
  fi
  if [[ "$ploidy" -eq 8 ]]; then
    if [[ -z "$downsample_8x" ]]; then
      downsample=400
    else
      downsample=$downsample_8x
    fi
  fi

  while IFS="" read -r i || [ -n "$i" ]; do (
    # awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) > ${i%.f*}_redun_head.sam
    # grep -v '^@' <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | awk -F"\t" 'BEGIN{FS=OFS="\t"} {$11="*"; print $0}' | cat ${i%.f*}_redun_head.sam - > ./alignment/${i%.f*}_redun.sam &&
    # rm ./alignment/${i%.f*}_redun.sam.gz && $gzip ./alignment/${i%.f*}_redun.sam &&
    # rm ${i%.f*}_redun_head.sam &&
    printf '\n###---'${i%.f*}'---###\n' > ${projdir}/alignment_summaries/${i%.f*}_summ.txt &&
    zcat ./alignment/${i%.f*}_redun.sam.gz | grep -v '^@' | awk '{gsub(/_/,"\t",$1);}1' | awk '{gsub(/se-/,"",$2);gsub(/pe-/,"",$2);}1' | tr ' ' '\t' > ${i%.f*}_full.sam &&
    split -l 10000 ${i%.f*}_full.sam ${i%.f*}_chunk_ && rm ${i%.f*}_full.sam &&
    find . -name '${i%.f*}_chunk_*' -print0 | xargs -0 -P "$prepthreads" -I{} bash -c 'awk '\''{for(i=0;i<=$2-1;i++) print $0}'\'' "$1" > "$1.out"' _ {} &&
    cat ${i%.f*}_chunk_* | awk '!($2="")1' | awk '{$1=$1"_"NR}1' | awk '{gsub(/ /,"\t");}1' > ${i%.f*}_full.sam &&
    cat <(zcat ./alignment/${i%.f*}_redun.sam.gz | grep '^@') ${i%.f*}_full.sam | gzip > ${i%.f*}_full.sam.gz &&
    rm ${i%.f*}_chunk_* ${i%.f*}_full.sam &&
    $samtools flagstat ${i%.f*}_full.sam.gz >> ${projdir}/alignment_summaries/${i%.f*}_summ.txt &&
    printf '########################################################################################################\n\n' >> ${projdir}/alignment_summaries/${i%.f*}_summ.txt &&
    printf 'copy_number\tFrequency\tPercentage\n' > ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_Read_histogram.txt &&
    $samtools view -F4 <(zcat ${i%.f*}_full.sam.gz 2> /dev/null) | awk '{print $1}' | awk '{gsub(/_pe-/,"\t");gsub(/seq/,"");}1' | \
    awk '{while ($2-- > 0) print $1}' | awk '{!seen[$0]++}END{for (i in seen) print seen[i]}' | awk '{!seen[$0]++}END{for (i in seen) print i, seen[i]}' > ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number.txt  &&
    awk 'NR==FNR{sum+= $2; next;} {printf("%s\t%s\t%3.3f%%\t%3.0f\n",$1,$2,100*$2/sum,100*$2/sum)}' ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number.txt ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number.txt | awk '$4 > 0' > ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot.txt &&
    unset IFS; printf "%s\t%s\t%s\t%*s\n" $(sed 's/$/ |/' ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot.txt) | tr ' ' '|' | sort -T ./tmp/ -k1,1 -n >> ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_Read_histogram.txt &&
    rm ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number.txt ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot.txt ${i%.f*}_full.sam.gz &&
    wait

    if [[ "$hap_ref" ]]; then
    max_hap=0
    for maxhap in ${haplome_number//,/ }; do if [[ "$maxhap" -gt "$max_hap" ]]; then max_hap=$maxhap; fi; done
    min_hap=999
    for minhap in ${haplome_number//,/ }; do if [[ "$minhap" -lt "$min_hap" ]]; then min_hap=$minhap; fi; done
    export max_hap=$maxhap
    export min_hap=$minhap


      if test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai; then
        if [[ "$paralogs" == false ]] && [[ "$uniquely_mapped" == true ]]; then
          awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) > ${i%.f*}_heading.sam &&
          $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
          awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0 || $5 >= min) {print $0}}' | awk -F '\t' 'BEGIN{OFS="\t"} {print $1"_"substr($3,length($3),1)"\t"$0}' | \
          awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | awk '$1==1{print $0}' | \
          awk '{$1=$2=""}1' | awk '$1=$1' | tr " " "\t" | awk -F '\t' 'BEGIN{OFS="\t"} {print $1"_"substr($3,1,length($3)-1)"\t"$0}' > ${i%.f*}_uniq1.sam &&
          awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' ${i%.f*}_uniq1.sam | tr " " "\t" | tr '\r' '\n' | \
          awk -v min=$min_hap -v max=$max_hap '$1>=min || $1<=max {print $0}' | awk '{$1=$2=""}1' | awk '$1=$1' | tr " " "\t" | cat ${i%.f*}_heading.sam - > ${i%.f*}_uniq2.sam &&
          $samtools sort -@ $gthreads -T ${projdir}/preprocess/tmp ${i%.f*}_uniq2.sam > ${i%.f*}_uniq2.bam &&
          rm ${i%.f*}_uniq1.sam ${i%.f*}_uniq2.sam &&
          $samtools view ${i%.f*}_uniq2.bam -o ${i%.f*}_uniqsorted.sam &&
          awk 'BEGIN{OFS="\t"} {print $1"_"$3"\t"$0}' ${i%.f*}_uniqsorted.sam | awk -F '\t' 'BEGIN{OFS="\t"}{print substr($1,1,length($1)-1)"\t"$0}' | awk '{$2=""}1' | \
          awk '!seen[$1]++' | awk '{$1=""}1' | awk '$1=$1' | awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | \
          tr -s " " | tr " " "\t" > ${i%.f*}_uniq.sam &&
          rm ${i%.f*}_uniq2.bam ${i%.f*}_uniqsorted.sam &&
          wait
        fi

        if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == false ]]; then
          awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) > ${i%.f*}_heading.sam &&
          $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
          awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0 || $5 >= min) {print $0}}' | awk -F '\t' 'BEGIN{OFS="\t"} {print $1"_"substr($3,length($3),1)"\t"$0}' | \
          awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | awk '$1<=6{print $0}' | awk '{$1=$2=""}1' | \
          awk '$1=$1' | tr " " "\t" | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {print $0}}' 2> /dev/null | \
          awk -F '\t' 'BEGIN{OFS="\t"} {print $1"_"$3"\t"$0}' | awk -F '\t' 'BEGIN{OFS="\t"} {$1=substr($1, 1, length($1)-1); print}' > ${i%.f*}_uniqeq1.sam &&
          awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' ${i%.f*}_uniqeq1.sam | tr " " "\t" | tr '\r' '\n' | \
          awk -v min=$min_hap -v max=$((max_hap*6)) '$1>=min || $1<=max{print $0}' | awk '{$1=$2=""}1' | awk '$1=$1' | tr " " "\t" > ${i%.f*}_uniqeq2.sam &&
          $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
          awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0 || $5 >= min) {print $0}}' | awk -F '\t' 'BEGIN{OFS="\t"} {print $1"_"substr($3,length($3),1)"\t"$0}' | \
          awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | awk '$1<=6{print $0}' | awk '{$1=$2=""}1' | \
          awk '$1=$1' | tr " " "\t" | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 > 0) {print $0}}' 2> /dev/null | awk -F '\t' 'BEGIN{OFS="\t"} {print $1"_"$3"\t"$0}' | \
          awk -F '\t' 'BEGIN{OFS="\t"} {$1=substr($1, 1, length($1)-1); print}' > ${i%.f*}_uniq1.sam &&
          awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' ${i%.f*}_uniq1.sam | tr " " "\t" | tr '\r' '\n' | \
          awk -v min=$min_hap -v max=$((max_hap*6)) '$1>=min || $1<=max{print $0}' | awk '{$1=$2=""}1' | awk '$1=$1' | tr " " "\t"  | cat - ${i%.f*}_uniqeq2.sam | \
          cat ${i%.f*}_heading.sam - > ${i%.f*}_uniq2.sam &&
          $samtools sort -@ $gthreads -T ${projdir}/preprocess/tmp ${i%.f*}_uniq2.sam > ${i%.f*}_uniq2.bam &&
          rm ${i%.f*}_uniqeq1.sam ${i%.f*}_uniqeq2.sam ${i%.f*}_uniq1.sam ${i%.f*}_uniq2.sam &&
          $samtools view ${i%.f*}_uniq2.bam -o ${i%.f*}_uniqsorted.sam &&
          awk -F '\t' 'BEGIN{OFS="\t"} {print $1"_"substr($3,length($3),1)"\t"$0}' ${i%.f*}_uniqsorted.sam | awk -F '\t' 'BEGIN{OFS="\t"} {print +a[$1] "\t" $0; a[$1]++}' | \
          awk '{if($1==0){print $0}}' | awk -F '\t' 'BEGIN{OFS="\t"}{print substr($2,1,length($2)-1)"\t"$0}' | awk '{$2=$3=""}1' | \
          awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | awk -v min=$min_hap -v max=$max_hap '$1>=min || $1<=max {print $3}' > ${i%.f*}_uniq2.list &&
          awk 'NR==FNR {a[$1]++; next} $1 in a' ${i%.f*}_uniq2.list ${i%.f*}_uniqsorted.sam | awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | \
          awk 'BEGIN{OFS="\t"} {print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" > ${i%.f*}_uniq.sam &&
          rm ${i%.f*}_uniq2.bam ${i%.f*}_uniqsorted.sam ${i%.f*}_uniq2.list &&
          wait
        fi
        if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == true ]]; then
          awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) > ${i%.f*}_heading.sam &&
          $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
          awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0 || $5 >= min) {print $0}}' | awk -F '\t' 'BEGIN{OFS="\t"} {print $1"_"substr($3,length($3),1)"\t"$0}' | \
          awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | awk '$1==1{print $0}' | \
          awk '{$1=$2=""}1' | awk '$1=$1' | tr " " "\t" | awk -F '\t' 'BEGIN{OFS="\t"} {print $1"_"substr($3,1,length($3)-1)"\t"$0}' > ${i%.f*}_uniq1.sam &&
          awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' ${i%.f*}_uniq1.sam | tr " " "\t" | tr '\r' '\n' | \
          awk -v min=$min_hap -v max=$max_hap '$1>=min || $1<=max {print $0}' | awk '{$1=$2=""}1' | awk '$1=$1' | tr " " "\t" | cat ${i%.f*}_heading.sam - > ${i%.f*}_uniq2.sam &&
          $samtools sort -@ $gthreads -T ${projdir}/preprocess/tmp ${i%.f*}_uniq2.sam > ${i%.f*}_uniq2.bam &&
          rm ${i%.f*}_uniq1.sam ${i%.f*}_uniq2.sam &&
          $samtools view ${i%.f*}_uniq2.bam -o ${i%.f*}_uniqsorted.sam &&
          awk 'BEGIN{OFS="\t"} {print $1"_"$3"\t"$0}' ${i%.f*}_uniqsorted.sam | awk -F '\t' 'BEGIN{OFS="\t"}{print substr($1,1,length($1)-1)"\t"$0}' | awk '{$2=""}1' | \
          awk '!seen[$1]++' | awk '{$1=""}1' | awk '$1=$1' | awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | \
          tr -s " " | tr " " "\t" > ${i%.f*}_uniqpart.sam &&
          rm ${i%.f*}_uniq2.bam ${i%.f*}_uniqsorted.sam &&
          wait
          awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) > ${i%.f*}_heading.sam &&
          $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
          awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0 || $5 >= min) {print $0}}' | awk -F '\t' 'BEGIN{OFS="\t"} {print $1"_"substr($3,length($3),1)"\t"$0}' | \
          awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | awk '$1<=6{print $0}' | awk '{$1=$2=""}1' | \
          awk '$1=$1' | tr " " "\t" | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {print $0}}' 2> /dev/null | \
          awk -F '\t' 'BEGIN{OFS="\t"} {print $1"_"$3"\t"$0}' | awk -F '\t' 'BEGIN{OFS="\t"} {$1=substr($1, 1, length($1)-1); print}' > ${i%.f*}_uniqeq1.sam &&
          awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' ${i%.f*}_uniqeq1.sam | tr " " "\t" | tr '\r' '\n' | \
          awk -v min=$min_hap -v max=$((max_hap*6)) '$1>=min || $1<=max{print $0}' | awk '{$1=$2=""}1' | awk '$1=$1' | tr " " "\t" > ${i%.f*}_uniqeq2.sam &&
          $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
          awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0 || $5 >= min) {print $0}}' | awk -F '\t' 'BEGIN{OFS="\t"} {print $1"_"substr($3,length($3),1)"\t"$0}' | \
          awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | awk '$1<=6{print $0}' | awk '{$1=$2=""}1' | \
          awk '$1=$1' | tr " " "\t" | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 > 0) {print $0}}' 2> /dev/null | awk -F '\t' 'BEGIN{OFS="\t"} {print $1"_"$3"\t"$0}' | \
          awk -F '\t' 'BEGIN{OFS="\t"} {$1=substr($1, 1, length($1)-1); print}' > ${i%.f*}_uniq1.sam &&
          awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' ${i%.f*}_uniq1.sam | tr " " "\t" | tr '\r' '\n' | \
          awk -v min=$min_hap -v max=$((max_hap*6)) '$1>=min || $1<=max{print $0}' | awk '{$1=$2=""}1' | awk '$1=$1' | tr " " "\t"  | cat - ${i%.f*}_uniqeq2.sam | \
          cat ${i%.f*}_heading.sam - > ${i%.f*}_uniq2.sam &&
          $samtools sort -@ $gthreads -T ${projdir}/preprocess/tmp ${i%.f*}_uniq2.sam > ${i%.f*}_uniq2.bam &&
          rm ${i%.f*}_uniqeq1.sam ${i%.f*}_uniqeq2.sam ${i%.f*}_uniq1.sam ${i%.f*}_uniq2.sam &&
          $samtools view ${i%.f*}_uniq2.bam -o ${i%.f*}_uniqsorted.sam &&
          awk -F '\t' 'BEGIN{OFS="\t"} {print $1"_"substr($3,length($3),1)"\t"$0}' ${i%.f*}_uniqsorted.sam | awk -F '\t' 'BEGIN{OFS="\t"} {print +a[$1] "\t" $0; a[$1]++}' | \
          awk '{if($1==0){print $0}}' | awk -F '\t' 'BEGIN{OFS="\t"}{print substr($2,1,length($2)-1)"\t"$0}' | awk '{$2=$3=""}1' | \
          awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | awk -v min=$min_hap -v max=$max_hap '$1>=min || $1<=max {print $3}' > ${i%.f*}_uniq2.list
          awk 'NR==FNR {a[$1]++; next} $1 in a' ${i%.f*}_uniq2.list ${i%.f*}_uniqsorted.sam | awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | \
          awk 'BEGIN{OFS="\t"} {print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" > ${i%.f*}_uniqAll.sam &&
          rm ${i%.f*}_uniq2.bam ${i%.f*}_uniqsorted.sam ${i%.f*}_uniq2.list &&
          wait
          cat ${i%.f*}_uniqpart.sam ${i%.f*}_uniqAll.sam | awk '!a[$0]++' > ${i%.f*}_uniq.sam &&
          rm ${i%.f*}_uniqpart.sam ${i%.f*}_uniqAll.sam
        fi
        wait

        awk '{while ($1-- > 0) print $0}' ${i%.f*}_uniq.sam | awk '{print "seq"NR"_"$0}' | tr -s ' ' | tr ' ' '\t' > ${i%.f*}_${ref1%.f*}.sam &&
        rm ${i%.f*}_uniq.sam &&
        awk '{print $4}' ${i%.f*}_${ref1%.f*}.sam | awk '{printf "%d00\n", $1/100}' | paste - ${i%.f*}_${ref1%.f*}.sam | \
        shuf | awk -v max="$downsample" 'a[$1$4]++ < max' | awk '{$1=""}1' | awk '{$1=$1};1' | tr -s " " | tr " " "\t" | awk '{$11=$10}1' | \
        awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$11); print}' | cat ${i%.f*}_heading.sam - > ${i%.f*}_${ref1%.f*}_downsample.sam &&
        :> ${i%.f*}_${ref1%.f*}.sam &&
        mv ${i%.f*}_${ref1%.f*}_downsample.sam ${i%.f*}_${ref1%.f*}.sam &&
        rm ${i%.f*}_heading.sam
        wait

        j="${i%.f*}_${ref1%.f*}.sam"
        cat ${j} | grep -vE 'I[[:digit:]]+I' | grep -vE 'D[[:digit:]]+D' | grep -vE 'D[[:digit:]]+I' | grep -vE 'I[[:digit:]]+D'  > cleaned_${j};  mv cleaned_${j} ${j} &&
        $java $Xmx2 -XX:ParallelGCThreads=$gthreads -Djava.io.tmpdir=${projdir}/preprocess/tmp -jar $picard SortSam I=$j O=${j%.sam*}.bam  SORT_ORDER=coordinate  VALIDATION_STRINGENCY=LENIENT TMP_DIR=${projdir}/preprocess/tmp && \
        $java $Xmx2 -XX:ParallelGCThreads=$gthreads -jar $picard BuildBamIndex INPUT=${j%.sam*}.bam VALIDATION_STRINGENCY=LENIENT && \
        $java $Xmx2 -XX:ParallelGCThreads=$gthreads -jar $picard AddOrReplaceReadGroups I=${j%.sam*}.bam O=${j%.sam*}_precall.bam RGLB=${i%.f*} RGPL=illumina RGPU=run RGSM=${i%.f*} VALIDATION_STRINGENCY=LENIENT && \
        $samtools index ${j%.sam*}_precall.bam &&
        rm $j ${j%.sam*}.bam ${j%.sam*}.bai
        wait
        if [[ $nodes -gt 1 ]]; then cp /tmp/${samples_list%.txt}/preprocess/${j%.sam*}_precall.bam* ${projdir}/preprocess/; fi
        wait
      fi
    else
      if test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai; then
        if [[ "$paralogs" == false ]] && [[ "$uniquely_mapped" == true ]]; then
          awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) > ${i%.f*}_heading.sam &&
          $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
          awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0 || $5 >= min) {print $0}}' | awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | \
          tr " " "\t" | tr '\r' '\n' | awk '$1==1{print $0}' | awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | \
          awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" > ${i%.f*}_uniq.sam &&
          wait
        fi
        if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == true ]]; then
          awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) > ${i%.f*}_heading.sam &&
          $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
          awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0 || $5 >= min) {print $0}}' | awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | \
          tr " " "\t" | tr '\r' '\n' | awk '$1==1 || $1<=6{print $0}' | awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | \
          awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {print $0}}' 2> /dev/null > ${i%.f*}_uniqeq.sam
          $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
          awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | awk '$1==1 || $1<=6{print $0}' | awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | \
          awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 > 0) {print $0}}' 2> /dev/null | cat - ${i%.f*}_uniqeq.sam | \
          awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" > ${i%.f*}_uniq.sam &&
          rm ${i%.f*}_uniqeq.sam &&
          wait
        fi
        if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == false ]]; then
          awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) > ${i%.f*}_heading.sam &&
          $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
          awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0 || $5 >= min) {print $0}}' | awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | \
          tr " " "\t" | tr '\r' '\n' | awk '$1==1{print $0}' | awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | \
          awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" > ${i%.f*}_uniq.sampart &&
          wait
          $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
          awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0 || $5 >= min) {print $0}}' | awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | \
          tr " " "\t" | tr '\r' '\n' | awk '$1==1 || $1<=6{print $0}' | awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | \
          awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {print $0}}' 2> /dev/null > ${i%.f*}_uniqeq.sam
          $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
          awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | awk '$1==1 || $1<=6{print $0}' | awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | \
          awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 > 0) {print $0}}' 2> /dev/null | cat - ${i%.f*}_uniqeq.sam | \
          awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" > ${i%.f*}_uniqAll.sam &&
          rm ${i%.f*}_uniqeq.sam &&
          wait
          awk 'NR==FNR{a[$0]=1;next}!a[$0]' ${i%.f*}_uniqpart.sam ${i%.f*}_uniqAll.sam > ${i%.f*}_uniq.sam &&
          ${i%.f*}_uniqpart.sam ${i%.f*}_uniqAll.sam
        fi
        wait


        awk '{while ($1-- > 0) print $0}' ${i%.f*}_uniq.sam | awk '{print "seq"NR"_"$0}' | tr -s ' ' | tr ' ' '\t' > ${i%.f*}_${ref1%.f*}.sam &&
        rm ${i%.f*}_uniq.sam &&
        awk '{print $4}' ${i%.f*}_${ref1%.f*}.sam | awk '{printf "%d00\n", $1/100}' | paste - ${i%.f*}_${ref1%.f*}.sam | \
        shuf | awk -v max="$downsample" 'a[$1$4]++ < max' | awk '{$1=""}1' | awk '{$1=$1};1' | tr -s " " | tr " " "\t" | awk '{$11=$10}1' | \
        awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$11); print}' | cat ${i%.f*}_heading.sam - > ${i%.f*}_${ref1%.f*}_downsample.sam &&
        :> ${i%.f*}_${ref1%.f*}.sam &&
        mv ${i%.f*}_${ref1%.f*}_downsample.sam ${i%.f*}_${ref1%.f*}.sam &&
        rm ${i%.f*}_heading.sam
        wait

        j="${i%.f*}_${ref1%.f*}.sam"
        cat ${j} | grep -vE 'I[[:digit:]]+I' | grep -vE 'D[[:digit:]]+D' | grep -vE 'D[[:digit:]]+I' | grep -vE 'I[[:digit:]]+D'  > cleaned_${j};  mv cleaned_${j} ${j} &&
        $java $Xmx2 -XX:ParallelGCThreads=$gthreads -Djava.io.tmpdir=${projdir}/preprocess/tmp -jar $picard SortSam I=$j O=${j%.sam*}.bam  SORT_ORDER=coordinate  VALIDATION_STRINGENCY=LENIENT TMP_DIR=${projdir}/preprocess/tmp && \
        $java $Xmx2 -XX:ParallelGCThreads=$gthreads -jar $picard BuildBamIndex INPUT=${j%.sam*}.bam VALIDATION_STRINGENCY=LENIENT && \
        $java $Xmx2 -XX:ParallelGCThreads=$gthreads -jar $picard AddOrReplaceReadGroups I=${j%.sam*}.bam O=${j%.sam*}_precall.bam RGLB=${i%.f*} RGPL=illumina RGPU=run RGSM=${i%.f*} VALIDATION_STRINGENCY=LENIENT && \
        $samtools index ${j%.sam*}_precall.bam &&
        rm $j ${j%.sam*}.bam ${j%.sam*}.bai
        wait
        if [[ $nodes -gt 1 ]]; then cp /tmp/${samples_list%.txt}/preprocess/${j%.sam*}_precall.bam* ${projdir}/preprocess/; fi
        wait
      fi
    fi
    ) &
    if [[ $(jobs -r -p | wc -l) -ge $prepN ]]; then
      wait
    fi
  done < <(cat ${projdir}/${samples_list})
  wait && touch ${projdir}/precall_done_${samples_list}
  wait
  find . -maxdepth 1 -type f -name "*" | grep -v 'precall' | xargs rm 2> /dev/null
  wait


	cd ${projdir}/samples
	if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
		precall=$(ls ${projdir}/precall_done_samples_list_node_* | wc -l)
		while [[ "$precall" -lt $nodes ]]; do sleep 300; precall=$(ls ${projdir}/precall_done_samples_list_node_* | wc -l); done
		sleep $((RANDOM % 2))
    if [[ $precall == $nodes ]] && test ! -f ${projdir}/alignment_summaries/refgenome_paralogs.txt; then
			cd ${projdir}/alignment_summaries
			cat *_summ.txt > alignment_summaries_reads.txt &&
      rm *_summ.txt &&
			# Total number of reads per samples
			awk '/###---/ || /QC-passed/{print}' alignment_summaries_reads.txt | cut -d\+ -f1 | tr -d '\n' | \
			awk  'gsub(/---###/, "\t", $0)' | awk  'gsub(/###---/, "", $0)' | tr ' ' '\n' > total_reads.txt &&
			# Total number of mapped reads per samples
			cat alignment_summaries_reads.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /0_mapped/{print}' |\
			tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
			awk 'gsub("\\+0mapped", "\t", $0)' | cut -d\: -f1 > total_reads_mapped.txt &&
			# Total number of mapped paired reads per samples
			cat alignment_summaries_reads.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /properly_paired/{print}' |\
			tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
			awk 'gsub("\\+0properlypaired", "\t", $0)' | cut -d\: -f1 > total_reads_paired.txt &&
			echo -e "Samples\tTotal\tMapped\tPerc_Mapped\tPE_Mapped\t%_PE_Mapped" > summary_precall.txt &&
			awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' total_reads_mapped.txt  total_reads.txt  | \
			awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' total_reads_paired.txt - | awk '{gsub(/ /,"\t"); print $0}' | \
			cat summary_precall.txt - | awk '{print $1"\t"$2"\t"$3"\t"$4}' > Tabulated_Alignment_Read_Summaries.txt &&
			rm total_* summary_precall.txt &> /dev/null &&
			rm ${projdir}/samples/metrics.txt ${projdir}/preprocess/metrics.txt &> /dev/null &&
      wait
		fi
    if [[ "$genomecov_est" == true ]] && [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
      cd ${projdir}/alignment_summaries
      printf "Sample\tGenome_Coverage(percentage)\n" > summary_genomecov.txt
      genome_size=$(awk '{print $3}' ../refgenomes/${ref1%.f*}.dict | awk '{gsub(/LN:/,"");}1' | awk '{s+=$1}END{print s}')
      for i in ../preprocess/alignment/*_redun.sam.gz; do
        $samtools view -bS <(zcat $i 2> /dev/null) | $samtools sort - > ${i%*_redun.sam.gz}.bam
        cov=$($bedtools genomecov -ibam ${i%*_redun.sam.gz}.bam -bga | awk '{print ($4>1)? 1 : $4}' | awk -v pat=$genome_size '{s+=$1}END{print (s/pat)*100}')
        printf "${i%*_redun.sam.gz}\t$cov\n"  | awk '{gsub(/..\/preprocess\/alignment\//,"");}1' >> summary_genomecov.txt
        rm ${i%*_redun.sam.gz}.bam 2> /dev/null
      done
      # for i in ../preprocess/*bam; do
      #   cov=$($bedtools genomecov -ibam $i -bga | awk '{print ($4>1)? 1 : $4}' | awk -v pat=$genome_size '{s+=$1}END{print (s/pat)*100}')
      #   printf "${i%_*_*}\t$cov\n"  | awk '{gsub(/..\/preprocess\//,"");}1' >> summary_genomecov.txt
      # done
    fi
		wait && touch ${projdir}/alignment_summary_done.txt
	fi

}
cd $projdir
if [ "$walkaway" == false ]; then
	echo -e "${magenta}- Do you want to perform read alignments and alignment post-processing? ${white}\n"
	read -p "- y(YES) or n(NO) " -t 36000 -n 1 -r
	if [[ ! $REPLY =~ ^[Yy]$ ]]; then
		printf '\n'
		echo -e "${magenta}- skipping read alignments and alignment post-processing ${white}\n"
	else
		printf '\n'
		if test -f ${projdir}/alignment_done.txt; then
			echo -e "${magenta}- read alignments and alignment post-processing already performed ${white}\n"
		else
			echo -e "${magenta}- performing read alignments and alignment post-processing ${white}\n"
			time main &>> log.out
		fi
	fi
fi
if [ "$walkaway" == true ]; then
	if [ "$alignments" == 1 ]; then
		if test -f ${projdir}/alignment_done.txt; then
			echo -e "${magenta}- read alignments and alignment post-processing already performed ${white}\n"
		else
			echo -e "${magenta}- performing read alignments and alignment post-processing ${white}\n"
			time main &>> log.out
		fi
	else
		echo -e "${magenta}- skipping read alignments and alignment post-processing ${white}\n"
	fi
fi


######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${yellow}\n- Performing Variant Calling with GATK HaplotypeCaller\n${blue}##############################################################################${white}\n"
main () {

cd ${projdir}
:> ${projdir}/compress_done.txt
:> ${projdir}/precall_done_${samples_list}

if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
  mkdir -p ${projdir}/preprocess/processed &&
	mv "${projdir}"/preprocess/processed/*_precall.bam* "${projdir}"/preprocess/ 2> /dev/null &&
	wait
fi
wait
:> call0

while [[ ! -f "${projdir}/call0" ]]; do sleep 5; done
wait
if [[ $nodes -gt 1 ]]; then
  mkdir -p /tmp/${samples_list%.txt}/refgenomes /tmp/${samples_list%.txt}/samples &&
  mkdir -p /tmp/${samples_list%.txt}/preprocess/tmp /tmp/${samples_list%.txt}/preprocess/alignment /tmp/${samples_list%.txt}/snpcall/tmp &&
  cp -rn ${projdir}/refgenomes/* /tmp/${samples_list%.txt}/refgenomes/ &&
  wait
  while IFS="" read -r i || [ -n "$i" ]; do
    sleep $((RANDOM % 2))
    cp -rn ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam* /tmp/${samples_list%.txt}/preprocess/ &&
    wait
  done < <(cat ${projdir}/${samples_list})
fi


cd $projdir
echo -e "${magenta}- performing SNP calling ${white}\n"
cd $projdir
cd preprocess
mkdir -p processed
if [[ ! -s $interval_list ]]; then export interval_list=""; fi

if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
	if [[ "$joint_calling" == true ]] && [[ "$variant_caller" == "gatk" ]]; then
		j=-I; input=""; k=""
		for i in *_precall.bam; do
			k="${j} ${i}"; input="${input} ${k}"
		done
		Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/${ref1%.*}.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' | awk -v pat=${ref1%.f*} '$0 ~ pat')
		if [[ -n "$Exclude_Chromosome" ]]; then
			for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
				Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
			done
		fi
		if [[ -z "$Get_Chromosome" ]]; then
      if [[ -z "$interval_list" ]]; then
  			for selchr in $Get2_Chromosome; do (
  				if [[ "$(ls ${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf* 2> /dev/null | wc -l)" -eq 0 ]]; then
  					$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/$ref1 -L ${selchr} ${input} -ploidy $ploidy -O ${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf.gz --max-reads-per-alignment-start 0 --minimum-mapping-quality $minmapq --dont-use-soft-clipped-bases $dont_use_softclip --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) &&
  					gunzip ${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf.gz &&
  					wait
  				fi
  				) &
  				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
  					wait
  				fi
  			done
      else
        for selchr in $Get2_Chromosome; do (
          cat ${projdir}/${interval_list} | grep $selchr > ${projdir}/variant_intervals_${selchr}.list &&
          if [[ "$(ls ${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf* 2> /dev/null | wc -l)" -eq 0 ]]; then
            $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/$ref1 -L ${projdir}/variant_intervals_${selchr}.list ${input} -ploidy $ploidy -O ${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf.gz --max-reads-per-alignment-start 0 --minimum-mapping-quality $minmapq --dont-use-soft-clipped-bases $dont_use_softclip --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) &&
            gunzip ${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf.gz &&
            rm ${projdir}/variant_intervals_${selchr}.list &&
            wait
          fi
          ) &
          if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
            wait
          fi
        done
      fi
			wait
			cd ../snpcall
			grep -h '^#' ${pop}_${ploidy}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt &&
			cat ${pop}_${ploidy}x_*_raw.vcf | awk '!/^#/' > all.vcf &&
			cat vcf_header.txt all.vcf > ${pop}_${ref1%.f*}_${ploidy}x_raw.vcf
			rm vcf_header.txt all.vcf *.vcf.gz.tb* ${pop}_${ploidy}x_*_raw.vcf 2> /dev/null &&
      wait
		else
			if [[ "$(ls ${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf* 2> /dev/null | wc -l)" -eq 0 ]]; then
				echo $Get_Chromosome | tr ',' '\n' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ${projdir}/refgenomes/${ref1%.fasta}.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2":1-"$3}' > ${projdir}/refgenomes/${ref1%.fasta}.list &&
				$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/$ref1 -L ${projdir}/refgenomes/${ref1%.fasta}.list ${input} -ploidy $ploidy -O ${projdir}/snpcall/${pop}_${ploidy}x_raw.vcf.gz --max-reads-per-alignment-start 0 --minimum-mapping-quality $minmapq --dont-use-soft-clipped-bases $dont_use_softclip --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) &&
				cd ../snpcall
				gunzip ${pop}_${ploidy}x_raw.vcf.gz &&
        mv ${pop}_${ploidy}x_raw.vcf ${pop}_${ref1%.f*}_${ploidy}x_raw.vcf
				wait
			fi
		fi
	fi
fi

if [[ "$joint_calling" == false ]]&& [[ "$variant_caller" == "gatk" ]]; then

  if [[ $nodes -eq 1 ]]; then cd ${projdir}/preprocess/; fi
  if [[ $nodes -gt 1 ]] && test -f ${projdir}/GBSapp_run_node_1.sh; then cd /tmp/${samples_list%.txt}/preprocess/; fi
  if [[ ! -d "${projdir}"/snpcall/cohorts_1 ]]; then
  	while IFS="" read -r i || [ -n "$i" ]; do (
  		sleep $((RANDOM % 2))
        if [[ "$(ls ${projdir}/snpcall/${pop}_${ref1%.f*}_${ploidy}x_raw.vcf* 2> /dev/null | wc -l)" -eq 0 ]]; then
    			if test ! -f "${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf"; then
    					if [[ -z "$Get_Chromosome" ]]; then
                if [[ -z "$interval_list" ]]; then
      						$GATK --java-options "$Xmxg -Djava.io.tmpdir=../snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ../refgenomes/$ref1 -I ${i%.f*}_${ref1%.f*}_precall.bam -ploidy $ploidy -O ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz -ERC GVCF --max-reads-per-alignment-start 0 --minimum-mapping-quality $minmapq --dont-use-soft-clipped-bases $dont_use_softclip --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) &&
                  wait
                else
                  $GATK --java-options "$Xmxg -Djava.io.tmpdir=../snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ../refgenomes/$ref1 -L ${projdir}/${interval_list} -I ${i%.f*}_${ref1%.f*}_precall.bam -ploidy $ploidy -O ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz -ERC GVCF --max-reads-per-alignment-start 0 --minimum-mapping-quality $minmapq --dont-use-soft-clipped-bases $dont_use_softclip --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) &&
      						wait
                fi
    					else
    						echo $Get_Chromosome | tr ',' '\n' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ../refgenomes/${ref1%.fasta}.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2":1-"$3}' > ../refgenomes/${ref1%.fasta}.list
    						$GATK --java-options "$Xmxg -Djava.io.tmpdir=../snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ../refgenomes/$ref1 -L ../refgenomes/${ref1%.fasta}.list -I ${i%.f*}_${ref1%.f*}_precall.bam -ploidy $ploidy -O ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz -ERC GVCF --max-reads-per-alignment-start 0 --minimum-mapping-quality $minmapq --dont-use-soft-clipped-bases $dont_use_softclip --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) &&
    						wait
    					fi
              wait
    					mv ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf.gz &&
              rm ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz.tb* 2> /dev/null &&
              gunzip ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf.gz &&
              $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" IndexFeatureFile -I ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf &&
    					wait
    			fi
    		fi
    		mv ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam* ${projdir}/preprocess/processed/ 2> /dev/null
        rm ${i%.f*}_${ref1%.f*}_precall.bam* 2> /dev/null
        wait ) &
  		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
  			wait
  		fi
  	done < <(cat ${projdir}/${samples_list})
  fi
	wait
  cd ${projdir}/preprocess/
	printf "variant calling completed on ${samples_list}" > ${projdir}/call1_${samples_list}

	if [[ "$samples_list" == "samples_list_node_1.txt" ]] && [[ "$(ls ${projdir}/snpcall/${pop}_${ref1%.f*}_${ploidy}x_raw.vcf* 2> /dev/null | wc -l)" -eq 0 ]]; then
		call1=$(ls ${projdir}/call1_samples_list_node_* | wc -l)
		while [[ "$call1" -lt $nodes ]]; do sleep 300; call1=$(ls ${projdir}/call1_samples_list_node_* | wc -l); done
		sleep $((RANDOM % 2))
    if [[ $call1 == $nodes ]]; then
			cd ${projdir}/snpcall
      if [[ ! -d cohorts_1 ]]; then
        mkdir -p cohorts_1
        mv *.g.vcf* ./cohorts_1/
      fi

			for dir in cohorts*/; do
				cd $dir
				j=--variant; input=""; k=""
				for i in *.g.vcf; do
					k="${j} ${i}"; input="${input} ${k}"
				done
				if [[ -z "$Get_Chromosome" ]]; then
					Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/${ref1%.*}.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' | awk -v pat=${ref1%.f*} '$0 ~ pat')
				else
					Get2_Chromosome=$(echo $Get_Chromosome | tr ',' '\n')
				fi
				if [[ -n "$Exclude_Chromosome" ]]; then
					for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
						Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
					done
				fi

				for selchr in $Get2_Chromosome; do (
          if [[ ! -f "${pop}_${ploidy}x_${selchr}_raw.vcf.gz.tbi" ]]; then
            rm -rf ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz 2> /dev/null
            rm -rf ${pop}_${ploidy}x_${selchr}_raw 2> /dev/null
            export TILEDB_DISABLE_FILE_LOCKING=1
            if [[ -z "$interval_list" ]]; then
              $GATK --java-options "$Xmx1 -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$loopthreads" GenomicsDBImport ${input} -L ${selchr} --genomicsdb-workspace-path ${pop}_${ploidy}x_${selchr}_raw --genomicsdb-shared-posixfs-optimizations true --batch-size 50 --merge-input-intervals &&
              wait
            else
              cat ${projdir}/${interval_list} | grep $selchr > ${projdir}/variant_intervals_${selchr}.list &&
              $GATK --java-options "$Xmx1 -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$loopthreads" GenomicsDBImport ${input} -L ${projdir}/variant_intervals_${selchr}.list --genomicsdb-workspace-path ${pop}_${ploidy}x_${selchr}_raw --genomicsdb-shared-posixfs-optimizations true --batch-size 50 --merge-input-intervals &&
              rm ${projdir}/variant_intervals_${selchr}.list &&
              wait
            fi
          fi
          wait ) &
      		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      			wait
      		fi
				done
        wait


				for selchr in $Get2_Chromosome; do (
          export TILEDB_DISABLE_FILE_LOCKING=1
					if [[ ! -f "${pop}_${ploidy}x_${selchr}_raw.vcf.gz.tbi" ]]; then
            if [[ -z "$interval_list" ]]; then
  						$GATK --java-options "$Xmx1 -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$loopthreads" GenotypeGVCFs -R ${projdir}/refgenomes/$ref1 -L ${selchr} -V gendb://${pop}_${ploidy}x_${selchr}_raw -O ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz &&
  						rm -r ${pop}_${ploidy}x_${selchr}_raw &&
  						mv ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz ${pop}_${ploidy}x_${selchr}_raw.vcf.gz &&
  						mv ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz.tbi ${pop}_${ploidy}x_${selchr}_raw.vcf.gz.tbi &&
  						wait
            else
              cat ${projdir}/${interval_list} | grep $selchr > ${projdir}/variant_intervals_${selchr}.list &&
              $GATK --java-options "$Xmx1 -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$loopthreads" GenotypeGVCFs -R ${projdir}/refgenomes/$ref1 -L ${projdir}/variant_intervals_${selchr}.list -V gendb://${pop}_${ploidy}x_${selchr}_raw -O ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz &&
              rm ${projdir}/variant_intervals_${selchr}.list &&
              rm -r ${pop}_${ploidy}x_${selchr}_raw &&
              mv ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz ${pop}_${ploidy}x_${selchr}_raw.vcf.gz &&
              mv ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz.tbi ${pop}_${ploidy}x_${selchr}_raw.vcf.gz.tbi &&
              wait
            fi
          fi
					if LC_ALL=C gzip -l ${pop}_${ploidy}x_${selchr}_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
						:
					else
						rm ../cohorts*/${pop}_${ploidy}x_*_raw.vcf.gz* 2> /dev/null
						rm ../cohorts*/${pop}_${ploidy}x_raw.vcf* 2> /dev/null
						rm ../${pop}_${ploidy}x_raw_cohorts*.vcf* 2> /dev/null
						echo -e "${magenta}- \n- SNP calling failed probably due to insufficient memory ${white}\n"
						echo -e "${magenta}- \n- Exiting pipeline in 5 seconds ${white}\n"
						sleep 5 && exit 1
					fi
          wait ) &
      		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      			wait
      		fi
				done
				wait
				for g in ${pop}_${ploidy}x_*_raw.vcf.gz; do
					gunzip $g
				done
				grep -h '^#' ${pop}_${ploidy}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt &&
				cat ${pop}_${ploidy}x_*_raw.vcf | awk '!/^#/' > all.vcf &&
				cat vcf_header.txt all.vcf > ${pop}_${ploidy}x_raw.vcf &&
				rm vcf_header.txt all.vcf &&
				rm ${pop}_${ploidy}x_*_raw.vcf ${pop}_${ploidy}x_*_raw.vcf.gz.tb* 2> /dev/null &&
				$bcftools view -I ${pop}_${ploidy}x_raw.vcf -O z -o ${pop}_${ploidy}x_raw.vcf.gz &&
				$bcftools index ${pop}_${ploidy}x_raw.vcf.gz &&

				$bcftools annotate -x FORMAT/PL ${pop}_${ploidy}x_raw.vcf.gz > ../${pop}_${ploidy}x_raw_${dir%/}.vcf &&
				cd ../
				$bcftools view -I ${pop}_${ploidy}x_raw_${dir%/}.vcf -O z -o ${pop}_${ploidy}x_raw_${dir%/}.vcf.gz &&
				$bcftools index ${pop}_${ploidy}x_raw_${dir%/}.vcf.gz &&
				wait
			done
			wait
			if [[ `ls -1 *cohorts*.vcf.gz 2> /dev/null | wc -l` -gt 1 ]]; then
				$bcftools merge *cohorts*.vcf.gz --force-samples -m all > ${pop}_${ref1%.f*}_${ploidy}x_raw.vcf &&
        wait
			else
				cp *cohorts*.vcf.gz ${pop}_${ref1%.f*}_${ploidy}x_raw.vcf.gz &&
				gunzip ${pop}_${ref1%.f*}_${ploidy}x_raw.vcf.gz &&
        wait
			fi

			if [[ "$keep_gVCF" == true ]]; then
				mkdir -p keep_gVCF &&
				mv ./cohorts*/*.g.vcf ./keep_gVCF &&
				rm -r *cohorts*
			else
				rm -r cohorts*
				rm -r *cohorts*
			fi
			cd ${projdir}/preprocess
		fi
	fi
	wait

	######################

	if [[ $nodes -gt 1 ]] && test -f ${projdir}/GBSapp_run_node_1.sh; then
		touch ${projdir}/queue_move_${samples_list%.txt}.txt
		queue_move=$(ls ${projdir}/queue_move_samples_list_node_* | wc -l)
		while [[ "$queue_move" -gt 1 ]]; do
			sleep $((RANDOM % 2))
      rm ${projdir}/queue_move_${samples_list%.txt}.txt; sleep $[ ( $RANDOM % 120 )  + 30 ]s
			:> ${projdir}/queue_move_${samples_list%.txt}.txt
			queue_move=$(ls ${projdir}/queue_move_samples_list_node_* | wc -l)
		done
		cd /tmp/${samples_list%.txt}/preprocess/
		mv * ${projdir}/preprocess/ && cd ${projdir}
		rm -rf /tmp/${samples_list%.txt} 2> /dev/null
		rm ${projdir}/queue_move_${samples_list%.txt}.txt
	fi
fi

if [[ "$variant_caller" == "bcftools" ]]; then
  cd ${projdir}/preprocess
  ls *.bam > bam_list.txt
  $bcftools mpileup -Ou -f ../refgenomes/$ref1 -b bam_list.txt | \
  $bcftools call -mv -Ov -o ${projdir}/snpcall/${pop}_${ref1%.f*}_${ploidy}x_raw.vcf &&
  rm bam_list.txt && cd ${projdir}/snpcall/ &&
  $bcftools view -I ${pop}_${ref1%.f*}_${ploidy}x_raw.vcf -O z -o ${pop}_${ref1%.f*}_${ploidy}x_raw.vcf.gz &&
  $bcftools index ${pop}_${ref1%.f*}_${ploidy}x_raw.vcf.gz &&
  mv ${projdir}/preprocess/*_${ref1%.f*}_precall.bam* ${projdir}/preprocess/processed/ 2> /dev/null
  wait
fi

}
cd $projdir
if [ "$walkaway" == false ]; then
	echo -e "${magenta}- Do you want to perform SNP/variant calling? ${white}\n"
	read -p "- y(YES) or n(NO) " -t 36000 -n 1 -r
	if [[ ! $REPLY =~ ^[Yy]$ ]]; then
		printf '\n'
		echo -e "${magenta}- skipping SNP/variant calling ${white}\n"
	else
		printf '\n'
		echo -e "${magenta}- performing SNP/variant calling ${white}\n"
		time main &>> log.out
	fi
fi
if [ "$walkaway" == true ]; then
	if [ "$snp_calling" == 1 ]; then
		echo -e "${magenta}- performing SNP/variant calling ${white}\n"
		time main &>> log.out
	else
		echo -e "${magenta}- skipping SNP/variant calling ${white}\n"
	fi
fi


######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${yellow}\n- Performing SNP Filtering\n${blue}##############################################################################${white}\n"
main () {
cd $projdir

if test -d snpfilter; then
	number_snpfilter=$( ls -d snpfilter* | wc -l )
	mv snpfilter snpfilter_"${number_snpfilter}"
fi


mkdir snpfilter
cd snpfilter

# large_numerous_chrom () {
#   #retrieve the original SNP positions from unsplit chromosomes/pseudomolecules.
#   cd $projdir
#   cd snpcall
#   if [[ "$checksplit" -gt 0 ]]; then
#   	for ln in *_raw0.vcf; do
#   		cp $ln ${ln%.vcf}_positions_split.vcf
#   		awk '/^#CHROM/{close("file.vcf"f);f++}{print $0 > "file"f}' $ln
#   		awk 'NR==1{print}' file1 | cat file - > file0.vcf
#   		for j in $(seq 1 10); do
#   			awk 'NR>1{print}' file1 | awk '{print $1,"\t",$2,"\t",$0}' | awk '{gsub(/_x0/,"\t"); print}' | \
#   			awk -v splits=$j -F '\t' 'BEGIN{OFS="\t"} $2 ~ splits {$3=$3+500000000}1' | awk 'BEGIN{OFS="\t"} !($2="")' | awk 'BEGIN{OFS="\t"} !($3="")' | \
#   			awk 'BEGIN{OFS="\t"} !($3="")' | awk 'BEGIN{OFS="\t"} !($3="")' | awk '{gsub(/\t\t/,"\t"); print }' | \
#   			sort -k1,1 -k2n,2 | awk NF > file1.vcf
#   		done
#   		cat file0.vcf file1.vcf > $ln
#   		rm file file1 file0.vcf file1.vcf
#   	done
#
#   	cd ${projdir}/alignment_summaries
#   	lnr=refgenome_paralogs.txt
#   	cp $lnr ${lnr%.txt}_positions_split.txt
#   	awk '/^CHROM/{close("file.txt"f);f++}{print $0 > "file"f}' $lnr
#   	awk 'NR==1{print}' file1 | cat file - > file0.txt
#   	for j in $(seq 1 10); do
#   		awk 'NR>1{print}' file1 | awk '{print $1,"\t",$2,"\t",$0}' | awk '{gsub(/_x0/,"\t"); print}' | \
#   		awk -v splits=$j -F '\t' 'BEGIN{OFS="\t"} $2 ~ splits {$3=$3+500000000}1' | awk 'BEGIN{OFS="\t"} !($2="")' | awk 'BEGIN{OFS="\t"} !($3="")' | \
#   		awk 'BEGIN{OFS="\t"} !($3="")' | awk 'BEGIN{OFS="\t"} !($3="")' | awk '{gsub(/\t\t/,"\t"); print }' | \
#   		sort -k1,1 -k2n,2 | awk NF > file1.txt
#   	done
#   	cat file0.txt file1.txt > $lnr
#   	rm file file1 file0.txt file1.txt
#   	cd ${projdir}/snpcall
#   else
#   	:
#   fi
#
#   #retrieve the original SNP positions from conigs/scaffolds before concatenation
#   cd $projdir
#   cd snpcall
#   export ncontigscaffold=$(grep '>' ${projdir}/refgenomes/${ref1%.fasta}_original.fasta &> /dev/null | wc -l)
#   if [[ ! -f ./index_code.txt ]]; then
#   	if [[ $ncontigscaffold -gt $max_pseudoMol ]]; then
#   		echo -e "${magenta}- retrieving SNP positions based on contigs/scaffold annotation ${white}\n"
#   		for nc in *_raw0.vcf; do
#   			if [[ "${nc}" =~ "vcf" ]]; then
#   				awk '/^#CHROM/{close("file.vcf"f);f++}{print $0 > "file"f}' $nc
#   				awk 'NR==1{print}' file1 | cat file - > file0.txt; rm file
#   				awk 'NR>1{print $1,"\t",$2}' file1 > file_snp_stitch.txt
#   				awk 'NR>1{print $0}' file1 > file1.txt; rm file1
#   			fi
#   			awk 'NR%10000==1{x="splitindex"++i;}{print > x}'  file_snp_stitch.txt & PID=$!
# 			  wait $PID
#   			for s in splitindex* ; do
#   				touch ${s}_index_code.txt
#   				while IFS="" read -r p || [ -n "$p" ]; do
#   					sleep $((RANDOM % 2))
#             fchr=$(echo $p | awk '{print $1}'); fpos=$(echo $p | awk '{print $2}')
#   					scf=$(awk -v fchr=$fchr '$1 == fchr' ${projdir}/refgenomes/contigscaffold_index.txt | \
#   					awk -v t=$fpos 'BEGIN {var=t; highest=0}{ j = $NF;if ( j <= var && j > highest ) { highest=j} } END {print highest}' )
#   					awk -v fchr=$fchr '$1 == fchr' ${projdir}/refgenomes/contigscaffold_index.txt | awk -v scf=$scf -v fchr=$fchr -v fpos=$fpos '$3 == scf {print fchr"\t"fpos"\t"$0}' >> ${s}_index_code.txt
#   				done < $s
#   			done
#   			wait
#   			cat *_index_code.txt > index_code.txt && rm splitindex*
#   			awk 'BEGIN {FS=OFS="\t"}{ $5 = $2 - $5 } 1' index_code.txt > file2.txt
#   			awk 'BEGIN {FS=OFS="\t"}!($3="")' file2.txt | awk '{gsub(/\t\t/,"\t"); print}' | awk '{gsub(/ /,""); print}' > index_code.txt
#   			paste -d "\t" index_code.txt file1.txt | awk -F"\t" 'BEGIN {FS=OFS="\t"}{$1=$2=$5=$6=""; print $0}' | tr -s '\t' | sed -e 's/^[ \t]*//' > file2.txt
#   			cat file0.txt file2.txt > $nc
#   			rm file0.txt file1.txt file2.txt index_code.txt
#   		done
#   		wait
#   		rm file_snp_stitch.txt
#
#
#   		cd ${projdir}/alignment_summaries
#   		awk 'NR==1{print}' refgenome_paralogs.txt > file0.txt
#   		awk 'NR>1{print $1,"\t",$2}' refgenome_paralogs.txt > file_snp_stitch.txt
#   		awk 'NR>1{print $0}' refgenome_paralogs.txt > file1.txt
#   		awk 'NR%10000==1{x="splitindex"++i;}{print > x}'  file_snp_stitch.txt & PID=$!
# 		  wait $PID
#   		for s in splitindex*; do
#   			while IFS="" read -r p || [ -n "$p" ]; do
#   			sleep $((RANDOM % 2))
#         fchr=$(echo $p | awk '{print $1}'); fpos=$(echo $p | awk '{print $2}')
#   			scf=$(awk -v fchr=$fchr '$1 == fchr' ${projdir}/refgenomes/contigscaffold_index.txt | \
#   			awk -v t=$fpos 'BEGIN {var=t; highest=0}{ j = $NF;if ( j <= var && j > highest ) { highest=j} } END {print highest}' )
#   			awk -v fchr=$fchr '$1 == fchr' ${projdir}/refgenomes/contigscaffold_index.txt | awk -v scf=$scf -v fchr=$fchr -v fpos=$fpos '$3 == scf {print fchr"\t"fpos"\t"$0}' >> ${s}_index_code.txt
#   			done < $s
#   		done
#   		wait
#   		cat *_index_code.txt > index_code.txt && rm splitindex*
#   		awk 'BEGIN {FS=OFS="\t"}{ $5 = $2 - $5 } 1' index_code.txt > file2.txt
#   		awk 'BEGIN {FS=OFS="\t"}!($3="")' file2.txt | awk '{gsub(/\t\t/,"\t"); print}' | awk '{gsub(/ /,""); print}' > index_code.txt
#   		paste -d "\t" index_code.txt file1.txt | awk -F"\t" 'BEGIN {FS=OFS="\t"}{$1=$2=$5=$6=""; print $0}' | tr -s '\t' | sed -e 's/^[ \t]*//' > file2.txt
#   		cat file0.txt file2.txt > refgenome_paralogs.txt
#   		rm file0.txt file1.txt file2.txt file_snp_stitch.txt
#   		cd ${projdir}/snpcall
#   	fi
#   fi
# }

if [[ "$ploidy" -eq 1 ]]; then
	mkdir 1x
fi
if [[ "$ploidy" -eq 2 ]]; then
	mkdir 2x
fi
if [[ "$ploidy" -eq 4 ]]; then
	mkdir 4x
fi
if [[ "$ploidy" -eq 6 ]]; then
	mkdir 6x
fi
if [[ "$ploidy" -eq 8 ]]; then
	mkdir 8x
fi

if [ -z $genotype_missingness ]; then
	genotype_missingness=1
else
	genotype_missingness=$( echo $genotype_missingness | awk '{gsub(/,/," ")}1' )
fi
if [ -z $sample_missingness ]; then
	sample_missingness=1
else
	sample_missingness=$( echo $sample_missingness | awk '{gsub(/,/," ")}1' )
fi
if [ -z $exclude_samples ]; then
	exclude_samples=NULL
fi
if [ -z $select_samples ]; then
	select_samples=NULL
fi
if [[ -z "$select_samples" ]]; then export select_samples=NULL; fi
if [[ ! -z "${projdir}/$select_samples" ]]; then
  awk 'BEGIN {FS = "\t";OFS = ""}{print $0,".fasta.gz"}' "${projdir}/$select_samples" | \
  awk '{gsub(/.fasta.gz.fasta.gz/,".fasta.gz");}1' > "${projdir}"/fetch_samples_seq.txt
fi
if [[ -z "${projdir}/$select_samples" ]]; then cat "${projdir}"/samples_list_node_* > "${projdir}"/fetch_samples_seq.txt; fi
if [[ -z "${projdir}/$select_samples" ]] && [[ "$exclude_samples" ]]; then
  echo $exclude_samples | tr ',' '\n' | awk 'BEGIN {FS = "\t";OFS = ""}{print $0,".fasta.gz"}' > "${projdir}"/fetch_excluded.txt
  grep -vFf "${projdir}"/fetch_excluded.txt "${projdir}"/fetch_samples_seq.txt > "${projdir}"/fetch_samples_seq0.txt
  mv "${projdir}"/fetch_samples_seq0.txt "${projdir}"/fetch_samples_seq.txt
fi


if [ -z $minRD_1x ]; then
	minRD_1x=2
fi
if [ -z $minRD_2x ]; then
	minRD_2x=6
fi
if [ -z $minRD_4x ]; then
	minRD_4x=25
fi
if [ -z $minRD_6x ]; then
	minRD_6x=45
fi
if [ -z $minRD_8x ]; then
	minRD_8x=100
fi
if [ -z $pseg ]; then
	pseg=0.001
fi
if [ -z $maf ]; then
	maf=0.02
fi

cd ${projdir}/snpcall
if [ -z "$(ls -A *_DP_GT.txt 2> /dev/null)" ]; then
	if [ -z "$(ls -A *_x.vcf 2> /dev/null)" ]; then
		if [ -z "$(ls -A *_raw.vcf 2> /dev/null)" ]; then
			for g in *_raw.vcf.gz; do gunzip $g;	done
		fi
	fi
fi
wait



file1xG=$( if [ "$(ls -A *_1x_DP_GT.txt 2> /dev/null)" ]; then ls *_1x_DP_GT.txt | wc -l;  else echo 0; fi )
file1xV=$( if [ "$(ls -A *_1x_raw.vcf 2> /dev/null)" ]; then ls *_1x_raw.vcf | head -n 1 | wc -l; else echo 0; fi )
if [[ "${file1xG}" -lt 1 ]]; then
	if [[ "${file1xV}" -gt 0 ]]; then
    if test ! -f ${projdir}/vcf1x_trimmed.txt; then
  		for i in *_1x_raw.vcf; do
        refg=${i%_1x_raw.vcf} && refg=${refg#*_} && refg=${refg%%_*}.fasta
        samz1=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz1=$((samz1*98)) && samz1=$((samz1/100)) &&
        samz2=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz2=$((samz2*80)) && samz2=$((samz2/100)) &&
        wait
        awk -v pat="0/0:0,0:0" -v samz1="$samz1" 'gsub(pat,pat) < samz1' $i | awk -v pat="./.:0,0:0"  -v samz2="$samz2" 'gsub(pat,pat) < samz2' | \
        awk -v pat="./.:.:."  -v samz2="$samz2" 'gsub(pat,pat) < samz2' > ${i}.tmp &&
        mv ${i}.tmp ${i} &&

        wait
        $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" IndexFeatureFile -I $i &&
  			$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$refg -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac &&
  			wait
  			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf &&
  			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%10000==2{x=file"1x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
  		  wait $PID
  			rm "${i%.vcf}"0.vcf* "${i%.vcf}"trim.vcf* &&
        wait
  		done
      :> ${projdir}/vcf1x_trimmed.txt
    fi
    samz1=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz1=$((samz1*98)) && samz1=$((samz1/100)) &&
    samz2=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz2=$((samz2*80)) && samz2=$((samz2/100)) &&
    wait
    for ptrimvcf in *rawSPLIT*.vcf; do
      awk -v pat="0:0,0:0" -v samz1="$samz1" 'gsub(pat,pat) < samz1' $ptrimvcf | awk -v pat=".:0,0:0"  -v samz2="$samz2" 'gsub(pat,pat) < samz2' | \
      awk -v pat="./.:.:."  -v samz2="$samz2" 'gsub(pat,pat) < samz2' > ${ptrimvcf}.tmp &&
      mv "${ptrimvcf}".tmp "${ptrimvcf}" &&
      wait
    done
    wait
		Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 1x "${GBSapp_dir}/tools/R" "1" "$filter_ExcHet" 2> /dev/null &&
    rm "${projdir}"/vcf1x_trimmed.txt 2> /dev/null &&
    wait
	fi
fi
wait
file2xG=$( if [ "$(ls -A *_2x_DP_GT.txt 2> /dev/null)" ]; then ls *_2x_DP_GT.txt | wc -l;  else echo 0; fi )
file2xV=$( if [ "$(ls -A *_2x_raw.vcf 2> /dev/null)" ]; then ls *_2x_raw.vcf | head -n 1 | wc -l; else echo 0; fi )
if [[ "${file2xG}" -lt 1 ]]; then
	if [[ "${file2xV}" -gt 0 ]]; then
    if test ! -f ${projdir}/vcf2x_trimmed.txt; then
  		for i in *_2x_raw.vcf; do
        refg=${i%_2x_raw.vcf} && refg=${refg#*_} && refg=${refg%%_*}.fasta
        if [[ "$(cat ${projdir}/samples_list_node_* | wc -l)" -ge 100 ]]; then
          samz1=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz1=$((samz1*98)) && samz1=$((samz1/100)) &&
          samz2=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz2=$((samz2*80)) && samz2=$((samz2/100)) &&
          wait
          awk -v pat="0/0:0,0:0" -v samz1="$samz1" 'gsub(pat,pat) < samz1' $i | awk -v pat="./.:0,0:0"  -v samz2="$samz2" 'gsub(pat,pat) < samz2' | \
          awk -v pat="./.:.:."  -v samz2="$samz2" 'gsub(pat,pat) < samz2' > ${i}.tmp &&
          mv ${i}.tmp ${i}  &&
          wait
        fi
        wait
        $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" IndexFeatureFile -I $i &&
  			$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$refg -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac &&
  			wait
  			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf &&
  			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%10000==2{x=file"2x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
  			wait $PID
  			rm "${i%.vcf}"0.vcf* "${i%.vcf}"trim.vcf* &&
        wait
  		done
      :> ${projdir}/vcf2x_trimmed.txt
    fi
    if [[ "$(cat ${projdir}/samples_list_node_* | wc -l)" -ge 100 ]]; then
      samz1=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz1=$((samz1*98)) && samz1=$((samz1/100)) &&
      samz2=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz2=$((samz2*80)) && samz2=$((samz2/100)) &&
      wait
      for ptrimvcf in *rawSPLIT*.vcf; do
        awk -v pat="0/0:0,0:0" -v samz1="$samz1" 'gsub(pat,pat) < samz1' $ptrimvcf | awk -v pat="./.:0,0:0"  -v samz2="$samz2" 'gsub(pat,pat) < samz2' | \
        awk -v pat="./.:.:."  -v samz2="$samz2" 'gsub(pat,pat) < samz2' > ${ptrimvcf}.tmp &&
        mv "${ptrimvcf}".tmp "${ptrimvcf}" &&
        wait
      done
      wait
    fi
    wait
		Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 2x "${GBSapp_dir}/tools/R" "1" "$filter_ExcHet" 2> /dev/null &&
    rm ${projdir}/vcf2x_trimmed.txt 2> /dev/null &&
    wait
	fi
fi
wait
file4xG=$( if [ "$(ls -A *_4x_DP_GT.txt 2> /dev/null)" ]; then ls *_4x_DP_GT.txt | wc -l;  else echo 0; fi )
file4xV=$( if [ "$(ls -A *_4x_raw.vcf 2> /dev/null)" ]; then ls *_4x_raw.vcf | head -n 1 | wc -l; else echo 0; fi )
if [[ "${file4xG}" -lt 1 ]]; then
	if [[ "${file4xV}" -gt 0 ]]; then
    if test ! -f ${projdir}/vcf4x_trimmed.txt; then
  		for i in *_4x_raw.vcf; do
        refg=${i%_4x_raw.vcf} && refg=${refg#*_} && refg=${refg%%_*}.fasta
        if [[ "$(cat ${projdir}/samples_list_node_* | wc -l)" -ge 100 ]]; then
          samz1=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz1=$((samz1*98)) && samz1=$((samz1/100)) &&
          samz2=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz2=$((samz2*80)) && samz2=$((samz2/100)) &&
          wait
          awk -v pat="0/0/0/0:0,0:0" -v samz1="$samz1" 'gsub(pat,pat) < samz1' $i | awk -v pat="./././.:0,0:0"  -v samz2="$samz2" 'gsub(pat,pat) < samz2' | \
          awk -v pat="./././.:.:."  -v samz2="$samz2" 'gsub(pat,pat) < samz2' > ${i}.tmp &&
          mv ${i}.tmp ${i} &&
          wait
        fi
        wait
        $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" IndexFeatureFile -I $i &&
  			$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$refg -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac &&
  			wait
  			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf &&
  			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%10000==2{x=file"4x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
  			wait $PID
  			rm "${i%.vcf}"0.vcf* "${i%.vcf}"trim.vcf* &&
        wait
  		done
      :> ${projdir}/vcf4x_trimmed.txt
    fi
    if [[ "$(cat ${projdir}/samples_list_node_* | wc -l)" -ge 100 ]]; then
      samz1=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz1=$((samz1*98)) && samz1=$((samz1/100)) &&
      samz2=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz2=$((samz2*80)) && samz2=$((samz2/100)) &&
      wait
      for ptrimvcf in *rawSPLIT*.vcf; do
        awk -v pat="0/0/0/0:0,0:0" -v samz1="$samz1" 'gsub(pat,pat) < samz1' $ptrimvcf | awk -v pat="./././.:0,0:0"  -v samz2="$samz2" 'gsub(pat,pat) < samz2' | \
        awk -v pat="./././.:.:."  -v samz2="$samz2" 'gsub(pat,pat) < samz2' > ${ptrimvcf}.tmp &&
        mv "${ptrimvcf}".tmp "${ptrimvcf}" &&
        wait
      done
      wait
    fi
    wait
		Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 4x "${GBSapp_dir}/tools/R" "1" "$filter_ExcHet" 2> /dev/null &&
    rm ${projdir}/vcf4x_trimmed.txt 2> /dev/null &&
    wait
	fi
fi
wait
file6xG=$( if [ "$(ls -A *_6x_DP_GT.txt 2> /dev/null)" ]; then ls *_6x_DP_GT.txt | wc -l;  else echo 0; fi )
file6xV=$( if [ "$(ls -A *_6x_raw.vcf 2> /dev/null)" ]; then ls *_6x_raw.vcf | head -n 1 | wc -l; else echo 0; fi )
if [[ "${file6xG}" -lt 1 ]]; then
	if [[ "${file6xV}" -gt 0 ]]; then
    if test ! -f ${projdir}/vcf6x_trimmed.txt; then
  		for i in *_6x_raw.vcf; do
        refg=${i%_6x_raw.vcf} && refg=${refg#*_} && refg=${refg%%_*}.fasta
        if [[ "$(cat ${projdir}/samples_list_node_* | wc -l)" -ge 100 ]]; then
          samz1=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz1=$((samz1*98)) && samz1=$((samz1/100)) &&
          samz2=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz2=$((samz2*80)) && samz2=$((samz2/100)) &&
          wait
          awk -v pat="0/0/0/0/0/0:0,0:0" -v samz1="$samz1" 'gsub(pat,pat) < samz1' $i | awk -v pat="./././././.:0,0:0"  -v samz2="$samz2" 'gsub(pat,pat) < samz2' | \
          awk -v pat="./././././.:.:."  -v samz2="$samz2" 'gsub(pat,pat) < samz2' > ${i}.tmp &&
          mv ${i}.tmp ${i} &&
          wait
        fi
        wait
        $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" IndexFeatureFile -I $i &&
  			$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$refg -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac &&
        wait
  			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf &&
  			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%10000==2{x=file"6x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
  			wait $PID
        rm "${i%.vcf}"0.vcf* "${i%.vcf}"trim.vcf* &&
        wait
  		done
      :> ${projdir}/vcf6x_trimmed.txt
    fi
    if [[ "$(cat ${projdir}/samples_list_node_* | wc -l)" -ge 100 ]]; then
      samz1=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz1=$((samz1*98)) && samz1=$((samz1/100)) &&
      samz2=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz2=$((samz2*80)) && samz2=$((samz2/100)) &&
      wait
      for ptrimvcf in *rawSPLIT*.vcf; do
        awk -v pat="0/0/0/0/0/0:0,0:0" -v samz1="$samz1" 'gsub(pat,pat) < samz1' $ptrimvcf | awk -v pat="./././././.:0,0:0"  -v samz2="$samz2" 'gsub(pat,pat) < samz2' | \
        awk -v pat="./././././.:.:."  -v samz2="$samz2" 'gsub(pat,pat) < samz2' > ${ptrimvcf}.tmp &&
        mv "${ptrimvcf}".tmp "${ptrimvcf}" &&
        wait
      done
      wait
    fi
    wait
		Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 6x "${GBSapp_dir}/tools/R" "1" "$filter_ExcHet" 2> /dev/null &&
    rm ${projdir}/vcf6x_trimmed.txt 2> /dev/null &&
    wait
	fi
fi
wait
file8xG=$( if [ "$(ls -A *_8x_DP_GT.txt 2> /dev/null)" ]; then ls *_8x_DP_GT.txt | wc -l;  else echo 0; fi )
file8xV=$( if [ "$(ls -A *_8x_raw.vcf 2> /dev/null)" ]; then ls *_8x_raw.vcf | head -n 1 | wc -l; else echo 0; fi )
if [[ "${file8xG}" -lt 1 ]]; then
	if [[ "${file8xV}" -gt 0 ]]; then
    if test ! -f ${projdir}/vcf8x_trimmed.txt; then
  		for i in *_8x_raw.vcf; do
        refg=${i%_8x_raw.vcf} && refg=${refg#*_} && refg=${refg%%_*}.fasta
        if [[ "$(cat ${projdir}/samples_list_node_* | wc -l)" -ge 100 ]]; then
          samz1=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz1=$((samz1*98)) && samz1=$((samz1/100)) &&
          samz2=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz2=$((samz2*80)) && samz2=$((samz2/100)) &&
          wait
          awk -v pat="0/0/0/0/0/0/0/0:0,0:0" -v samz1="$samz1" 'gsub(pat,pat) < samz1' $i | awk -v pat="./././././././.:0,0:0"  -v samz2="$samz2" 'gsub(pat,pat) < samz2' | \
          awk -v pat="./././././././.:.:."  -v samz2="$samz2" 'gsub(pat,pat) < samz2' > ${i}.tmp &&
          mv ${i}.tmp ${i} &&
          wait
        fi
        wait
        $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" IndexFeatureFile -I $i &&
  			$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$refg -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac &&
        wait
  			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf &&
  			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%10000==2{x=file"8x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
        wait $PID
        rm "${i%.vcf}"0.vcf* "${i%.vcf}"trim.vcf* &&
        wait
  		done
      :> ${projdir}/vcf8x_trimmed.txt
    fi
    if [[ "$(cat ${projdir}/samples_list_node_* | wc -l)" -ge 100 ]]; then
      samz1=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz1=$((samz1*98)) && samz1=$((samz1/100)) &&
      samz2=$(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') && samz2=$((samz2*80)) && samz2=$((samz2/100)) &&
      wait
      for ptrimvcf in *rawSPLIT*.vcf; do
        awk -v pat="0/0/0/0/0/0/0/0:0,0:0" -v samz1="$samz1" 'gsub(pat,pat) < samz1' $ptrimvcf | awk -v pat="./././././././.:0,0:0"  -v samz2="$samz2" 'gsub(pat,pat) < samz2' | \
        awk -v pat="./././././././.:.:."  -v samz2="$samz2" 'gsub(pat,pat) < samz2' > ${ptrimvcf}.tmp &&
        mv "${ptrimvcf}".tmp "${ptrimvcf}" &&
        wait
      done
      wait
    fi
    wait
		Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 8x "${GBSapp_dir}/tools/R" "1" "$filter_ExcHet" 2> /dev/null &&
    rm ${projdir}/vcf8x_trimmed.txt 2> /dev/null &&
    wait
	fi
fi
wait


filetest=*x.vcf*
if [ -z "$(ls -A *x.vcf* 2> /dev/null)" ]; then
	for v in *_DP_GT.txt; do
		vcfdose=${v%_DP*}; vcfdose=${vcfdose#*_}; out=$(ls *${vcfdose}_raw.vcf)
		for raw in $out; do
			grep '^#' $raw  > ${raw%_raw.vcf}.vcf &&
			awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$1,$2]){print b}}' $raw <(awk '{print $1"\t"$2}' $v | awk '!seen[$0] {print} {++seen[$0]}') >> ${raw%_raw.vcf}.vcf &&
      $bcftools sort ${raw%_raw.vcf}.vcf | gzip > ${raw%_raw.vcf}.vcf.gz &&
      rm ${raw%_raw.vcf}.vcf &&
      wait
		done
	done
	wait
fi

if [ "$(ls -A *.vcf 2> /dev/null)" ]; then
	for v in *.vcf; do gzip $v; done
fi
wait


cd ${projdir}/snpcall
if [ `ls -1 *biallelic* 2>/dev/null | wc -l` -eq 0 ]; then
  for bial in *_DP_GT.txt; do
    bialv=${bial%_DP_GT.txt} &&
    grep -v '^#' <(zcat *${bialv#*_}_raw.vcf.gz) | awk '{print $1"\t"$2"\t"$5}' | grep -v ',' | awk '{print $1"\t"$2}' | \
    awk 'NR==FNR{a[$1,$2]=$0;next}(($1,$2) in a)' - <(awk '{gsub(/ /,"\t");}1' $bial) | cat <(head -n1  $bial | awk '{gsub(/ /,"\t");}1') - > ${bial%.txt}_biallelic.txt
  done
fi


cd $projdir
cd samples
# window1=$(ls -S | head -1 | xargs zcat -fq | awk '{ print length }' | sort -n | tail -1)
# window=$((window1 + 20))
cd ${projdir}


for smiss in ${sample_missingness//,/ }; do
  for gmiss in ${genotype_missingness//,/ }; do
    if [[ -z "$p1" ]]; then
    	if [ -d "${projdir}/snpfilter/1x" ]; then
    		cd ${projdir}/snpfilter &&
    		cp -r 1x 1x_diversity_gmiss"${gmiss}"_smiss"${smiss}" &&
    		cd ./1x_diversity_gmiss"${gmiss}"_smiss"${smiss}" &&
    		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_1x.R $pop $gmiss $smiss $minRD_1x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome_number $biallelic $select_samples 2> /dev/null &&
    		rm "${pop}"_1x_rawRD"${minRD_1x}"_DP_GT.txt "${pop}"_1x_DP_GT.txt "${pop}"_1x_rd"${minRD_1x}".txt &&
    		mkdir visualizations && mv ./*.tiff ./visualizations/ 2> /dev/null &&
        wait

        if [[ "$variant_intervals" == true ]]; then
          mkdir -p variant_intervals &&
          cd variant_intervals &&
          wait
          while IFS="" read -r p || [ -n "$p" ]; do
            sleep $((RANDOM % 2))
            cat ../${pop}_1x_rd${minRD_1x}_maf${maf}_dose.txt | grep $p | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
            paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
            awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
            rm variant_intervals_${p}.txt &&
            :> variant_intervals_${p}.txt &&
            grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
            awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
            chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep $p | awk '{gsub(/LN:/,""); print $3}') &&
            awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
            rm vbreak_break.txt variant_intervals_${p}.tmp &&
            wait
            for vbreak in vbreak_*; do
              awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
              awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt &&
              wait
            done
            sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
            mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
            rm vbreak_*
            wait
          done < <(awk 'NR>1{print $2}' ../${pop}_1x_rd${minRD_1x}_maf${maf}_dose.txt | sort | uniq)
          cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
          rm variant_intervals_*.txt && cd ../
          rm -rf variant_intervals &&
          wait
        fi

        # Extract sequence context of variants
        mv ${projdir}/preprocess/processed/*_precall.bam ${projdir}/preprocess/ 2> /dev/null
        if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            awk '{print $1"\t"$2"\t"$3}' ${pop}_1x_rd${minRD_1x}_maf${maf}_dose.txt > snplist.txt &&
            awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
            awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
            for wint in {100,200}; do
              awk -v wint="$wint" '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
              awk -v wint="$wint" '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
              wait
            done
            awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
            mv snplist_haps.tmp snplist_haps.txt &&
            rm snplist_round.txt &&
            wait
          fi
          if [[ -n "$RE1" ]]; then
            mkdir -p seq_context_"${RE1}" &&
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
              rm ./seq_context_"${RE1}"/${i%.f*}_seqcontext.tmp 2> /dev/null &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
          fi
          if [[ -n "$RE2" ]]; then
            mkdir -p seq_context_"${RE2}" &&
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
              rm ./seq_context_"${RE2}"/${i%.f*}_seqcontext.tmp 2> /dev/null &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
          fi
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            mkdir consensus_seq_context
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              cat ./seq_context_"${RE1}"/"${i%.f*}"* ./seq_context_"${RE2}"/"${i%.f*}"* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
              maxlen=$(awk '{print $2}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}') &&
              awk -v maxlen="$maxlen" '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
              mv ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta 2> /dev/null &&
              while IFS="" read -r aln || [ -n "$aln" ]; do
    						awk -v aln="$aln" '$1 == aln {print}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
                bash $mafft --quiet --localpair --maxiterate 1000 ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp > ./consensus_seq_context/"${i%.f*}"_seqcontext.msf &&
    						bash $consambig -sequence ./consensus_seq_context/"${i%.f*}"_seqcontext.msf -outseq ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp 2> /dev/null &&
    						awk -v aln="${i%.f*}~~~${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp >> ./consensus_seq_context/"${i%.f*}"_seqcontext.cons &&
    						wait
    					done < <(awk '{print $1}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | sort | uniq )
              rm ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf 2> /dev/null &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
            wait
            for nn in ./consensus_seq_context/*.cons; do
              for run in $(seq 1 10); do
                awk '{gsub(/n$|N$/,"");}1' "$nn" 2> /dev/null > "${nn}".tmp && mv "${nn}".tmp "${nn}" 2> /dev/null &&
                wait
              done
            done
            wait
            seqid=$(cat ./consensus_seq_context/*.cons 2> /dev/null | grep '>' | awk '{gsub(/>/,""); gsub(/~~~Chr/,"\t~~~Chr");}1' | awk '{print $2}' | sort | uniq)
            wait
            mkdir -p ./consensus_seq_context/seqid_combine_samples
            wait
            for splitseqid in $seqid; do
              filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}" && filename="${filename/~~~}" &&
              cat ./consensus_seq_context/*.cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
              awk -v pat="$splitseqid" '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/"${filename}" &&
              minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{print length}' | sort -n | head -n 1) &&
              awk -v minlen="$minlen" '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{gsub(/n/,"N",$2); gsub(/~~~/,"\t");}1'| awk '{print $1"\n"$3}' > ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp &&
              mv ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp ./consensus_seq_context/seqid_combine_samples/"${filename}" 2> /dev/null &&
              wait
            done
            wait
            mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/ 2> /dev/null &&
            mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*.cons ./consensus_seq_context/sample_consensus_seqs 2> /dev/null &&
            mv seq_context_TGCAT ./consensus_seq_context/ 2> /dev/null &&
            mv seq_context_CATG ./consensus_seq_context/ 2> /dev/null &&
            wait
          fi
        fi
    		find . -type f -empty -delete
    	fi
    	wait
    	if [ -d "${projdir}/snpfilter/2x" ]; then
    		cd ${projdir}/snpfilter &&
    		cp -r 2x 2x_diversity_gmiss"${gmiss}"_smiss"${smiss}" &&
    		cd ./2x_diversity_gmiss"${gmiss}"_smiss"${smiss}" &&
    		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_2x.R $pop $gmiss $smiss $minRD_2x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome_number $biallelic $select_samples 2> /dev/null &&
    		rm "${pop}"_2x_rawRD"${minRD_2x}"_DP_GT.txt "${pop}"_2x_DP_GT.txt "${pop}"_2x_rd"${minRD_2x}".txt &&
    		mkdir visualizations && mv ./*.tiff ./visualizations/ 2> /dev/null &&
        wait

        if [[ "$variant_intervals" == true ]]; then
          mkdir -p variant_intervals &&
          cd variant_intervals &&
          wait
          while IFS="" read -r p || [ -n "$p" ]; do
            sleep $((RANDOM % 2))
            cat ../${pop}_2x_rd${minRD_2x}_maf${maf}_dose.txt | grep $p | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
            paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
            awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
            rm variant_intervals_${p}.txt &&
            :> variant_intervals_${p}.txt &&
            grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
            awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
            chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep $p | awk '{gsub(/LN:/,""); print $3}') &&
            awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
            rm vbreak_break.txt variant_intervals_${p}.tmp &&
            wait
            for vbreak in vbreak_*; do
              awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
              awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt &&
              wait
            done
            sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
            mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
            rm vbreak_* &&
            wait
          done < <(awk 'NR>1{print $2}' ../${pop}_2x_rd${minRD_2x}_maf${maf}_dose.txt | sort | uniq)
          cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
          rm variant_intervals_*.txt && cd ../
          rm -rf variant_intervals &&
          wait
        fi

        # Extract sequence context of variants
        mv ${projdir}/preprocess/processed/*_precall.bam ${projdir}/preprocess/ 2> /dev/null
        if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            awk '{print $1"\t"$2"\t"$3}' ${pop}_2x_rd${minRD_2x}_maf${maf}_dose.txt > snplist.txt &&
            awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
            awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
            wait
            for wint in {100,200}; do
              awk -v wint="$wint" '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
              awk -v wint="$wint" '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
              wait
            done
            awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
            mv snplist_haps.tmp snplist_haps.txt &&
            rm snplist_round.txt &&
            wait
          fi
          if [[ -n "$RE1" ]]; then
            mkdir -p seq_context_"${RE1}" &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
              rm ./seq_context_"${RE1}"/${i%.f*}_seqcontext.tmp 2> /dev/null &&
              wait
            done < <(awk 'NR==1{print $0}' ${pop}_2x_rd${minRD_2x}_maf${maf}_dose.txt | awk '{gsub(/SNP\tCHROM\tPOS\tREF\tALT\t/,""); gsub(/\t/,".fastq.gz\n");}1')
          fi
          if [[ -n "$RE2" ]]; then
            mkdir -p seq_context_"${RE2}" &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
              rm ./seq_context_"${RE2}"/${i%.f*}_seqcontext.tmp 2> /dev/null &&
              wait
            done < <(awk 'NR==1{print $0}' ${pop}_2x_rd${minRD_2x}_maf${maf}_dose.txt | awk '{gsub(/SNP\tCHROM\tPOS\tREF\tALT\t/,""); gsub(/\t/,".fastq.gz\n");}1')
          fi
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            mkdir consensus_seq_context &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              cat ./seq_context_"${RE1}"/"${i%.f*}"* ./seq_context_"${RE2}"/"${i%.f*}"* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
              maxlen=$(awk '{print $2}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}') &&
              awk -v maxlen="$maxlen" '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
              mv ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta 2> /dev/null &&
              while IFS="" read -r aln || [ -n "$aln" ]; do
    						awk -v aln="$aln" '$1 == aln {print}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
                bash $mafft --quiet --localpair --maxiterate 1000 ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp > ./consensus_seq_context/"${i%.f*}"_seqcontext.msf &&
    						bash $consambig -sequence ./consensus_seq_context/"${i%.f*}"_seqcontext.msf -outseq ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp 2> /dev/null &&
    						awk -v aln="${i%.f*}~~~${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp >> ./consensus_seq_context/"${i%.f*}"_seqcontext.cons &&
    						wait
    					done < <(awk '{print $1}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | sort | uniq )
              rm ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf 2> /dev/null &&
              wait
            done < <(awk 'NR==1{print $0}' ${pop}_2x_rd${minRD_2x}_maf${maf}_dose.txt | awk '{gsub(/SNP\tCHROM\tPOS\tREF\tALT\t/,""); gsub(/\t/,".fastq.gz\n");}1')
            wait
            for nn in ./consensus_seq_context/*.cons; do
              for run in $(seq 1 10); do
                awk '{gsub(/n$|N$/,"");}1' "$nn" 2> /dev/null > "${nn}".tmp && mv "${nn}".tmp "${nn}" 2> /dev/null &&
                wait
              done
            done
            wait
            seqid=$(cat ./consensus_seq_context/*.cons 2> /dev/null | grep '>' | awk '{gsub(/>/,""); gsub(/~~~Chr/,"\t~~~Chr");}1' | awk '{print $2}' | sort | uniq)
            wait
            mkdir -p ./consensus_seq_context/seqid_combine_samples
            wait
            for splitseqid in $seqid; do
              filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}" && filename="${filename/~~~}" &&
              cat ./consensus_seq_context/*.cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
              awk -v pat="$splitseqid" '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/"${filename}" &&
              minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{print length}' | sort -n | head -n 1) &&
              awk -v minlen="$minlen" '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{gsub(/n/,"N",$2); gsub(/~~~/,"\t");}1'| awk '{print $1"\n"$3}' > ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp &&
              mv ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp ./consensus_seq_context/seqid_combine_samples/"${filename}" 2> /dev/null &&
              wait
            done
            wait
            mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/ 2> /dev/null &&
            mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*.cons ./consensus_seq_context/sample_consensus_seqs 2> /dev/null &&
            mv seq_context_TGCAT ./consensus_seq_context/ 2> /dev/null &&
            mv seq_context_CATG ./consensus_seq_context/ 2> /dev/null &&
            wait
          fi
        fi
    		find . -type f -empty -delete
    	fi
    	wait
    	if [ -d "${projdir}/snpfilter/4x" ]; then
    		cd ${projdir}/snpfilter &&
    		cp -r 4x 4x_diversity_gmiss"${gmiss}"_smiss"${smiss}" &&
    		cd ./4x_diversity_gmiss"${gmiss}"_smiss"${smiss}" &&
    		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_4x.R $pop $gmiss $smiss $minRD_4x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome_number $biallelic $select_samples 2> /dev/null &&
    		rm ${pop}_4x_rawRD${minRD_4x}_DP_GT.txt ${pop}_4x_DP_GT.txt ${pop}_4x_rd${minRD_4x}.txt 2> /dev/null &&
    		mkdir visualizations && mv ./*.tiff ./visualizations/ 2> /dev/null &&
        wait

        if [[ "$variant_intervals" == true ]]; then
          mkdir -p variant_intervals &&
          cd variant_intervals &&
          wait
          while IFS="" read -r p || [ -n "$p" ]; do
            sleep $((RANDOM % 2))
            cat ../${pop}_4x_rd${minRD_4x}_maf${maf}_dose.txt | grep $p | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
            paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
            awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
            rm variant_intervals_${p}.txt &&
            :> variant_intervals_${p}.txt &&
            grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
            awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
            chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep $p | awk '{gsub(/LN:/,""); print $3}') &&
            awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
            rm vbreak_break.txt variant_intervals_${p}.tmp &&
            for vbreak in vbreak_*; do
              awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
              awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt &&
              wait
            done
            sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
            mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
            rm vbreak_* &&
            wait
          done < <(awk 'NR>1{print $2}' ../${pop}_4x_rd${minRD_4x}_maf${maf}_dose.txt | sort | uniq)
          cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
          rm variant_intervals_*.txt && cd ../
          rm -rf variant_intervals &&
          wait
        fi

        # Extract sequence context of variants
        mv ${projdir}/preprocess/processed/*_precall.bam ${projdir}/preprocess/ 2> /dev/null
        if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            awk '{print $1"\t"$2"\t"$3}' ${pop}_4x_rd${minRD_4x}_maf${maf}_dose.txt > snplist.txt &&
            awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
            awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
            for wint in {100,200}; do
              awk -v wint="$wint" '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
              awk -v wint="$wint" '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
              wait
            done
            awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
            mv snplist_haps.tmp snplist_haps.txt &&
            rm snplist_round.txt &&
            wait
          fi
          if [[ -n "$RE1" ]]; then
            mkdir -p seq_context_${RE1} &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
              rm ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
          fi
          if [[ -n "$RE2" ]]; then
            mkdir -p seq_context_${RE2} &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
              rm ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
          fi
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            mkdir consensus_seq_context &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              cat ./seq_context_"${RE1}"/"${i%.f*}"* ./seq_context_"${RE2}"/"${i%.f*}"* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
              maxlen=$(awk '{print $2}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}') &&
              awk -v maxlen="$maxlen" '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
              mv ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta 2> /dev/null &&
              while IFS="" read -r aln || [ -n "$aln" ]; do
    						awk -v aln="$aln" '$1 == aln {print}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
                bash $mafft --quiet --localpair --maxiterate 1000 ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp > ./consensus_seq_context/"${i%.f*}"_seqcontext.msf &&
    						bash $consambig -sequence ./consensus_seq_context/"${i%.f*}"_seqcontext.msf -outseq ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp 2> /dev/null &&
    						awk -v aln="${i%.f*}~~~${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp >> ./consensus_seq_context/"${i%.f*}"_seqcontext.cons &&
    						wait
    					done < <(awk '{print $1}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | sort | uniq )
              rm ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf 2> /dev/null &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
            wait
            for nn in ./consensus_seq_context/*.cons; do
              for run in $(seq 1 10); do
                awk '{gsub(/n$|N$/,"");}1' "$nn" 2> /dev/null > "${nn}".tmp && mv "${nn}".tmp "${nn}" 2> /dev/null &&
                wait
              done
            done
            seqid=$(cat ./consensus_seq_context/*.cons 2> /dev/null | grep '>' | awk '{gsub(/>/,""); gsub(/~~~Chr/,"\t~~~Chr");}1' | awk '{print $2}' | sort | uniq) &&
            mkdir -p ./consensus_seq_context/seqid_combine_samples
            wait
            for splitseqid in $seqid; do
              filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}" && filename="${filename/~~~}" &&
              cat ./consensus_seq_context/*.cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
              awk -v pat="$splitseqid" '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/"${filename}" &&
              minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{print length}' | sort -n | head -n 1) &&
              awk -v minlen="$minlen" '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{gsub(/n/,"N",$2); gsub(/~~~/,"\t");}1'| awk '{print $1"\n"$3}' > ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp &&
              mv ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp ./consensus_seq_context/seqid_combine_samples/"${filename}" 2> /dev/null &&
              wait
            done
            mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/ 2> /dev/null &&
            mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*.cons ./consensus_seq_context/sample_consensus_seqs 2> /dev/null &&
            mv seq_context_TGCAT ./consensus_seq_context/ 2> /dev/null &&
            mv seq_context_CATG ./consensus_seq_context/ 2> /dev/null &&
            wait
          fi
        fi
    		find . -type f -empty -delete
    	fi
    	wait
    	if [ -d "${projdir}/snpfilter/6x" ]; then
    		cd ${projdir}/snpfilter &&
    		cp -r 6x 6x_diversity_gmiss"${gmiss}"_smiss"${smiss}" &&
    		cd ./6x_diversity_gmiss"${gmiss}"_smiss"${smiss}" &&
    		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_6x.R $pop $gmiss $smiss $minRD_6x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome_number $biallelic $select_samples 2> /dev/null &&
    		rm ${pop}_6x_rawRD${minRD_6x}_DP_GT.txt ${pop}_6x_DP_GT.txt ${pop}_6x_rd${minRD_6x}.txt 2> /dev/null &&
    		mkdir visualizations && mv ./*.tiff ./visualizations/ 2> /dev/null &&
        wait

        if [[ "$variant_intervals" == true ]]; then
          mkdir -p variant_intervals
          cd variant_intervals
          while IFS="" read -r p || [ -n "$p" ]; do
            sleep $((RANDOM % 2))
            cat ../${pop}_6x_rd${minRD_6x}_maf${maf}_dose.txt | grep $p | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
            paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
            awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
            rm variant_intervals_${p}.txt &&
            :> variant_intervals_${p}.txt &&
            grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
            awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
            chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep $p | awk '{gsub(/LN:/,""); print $3}') &&
            awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
            rm vbreak_break.txt variant_intervals_${p}.tmp &&
            wait
            for vbreak in vbreak_*; do
              awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
              awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt
              wait
            done
            sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
            mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
            rm vbreak_* &&
            wait
          done < <(awk 'NR>1{print $2}' ../${pop}_6x_rd${minRD_6x}_maf${maf}_dose.txt | sort | uniq)
          cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
          rm variant_intervals_*.txt && cd ../
          rm -rf variant_intervals &&
          wait
        fi

        # Extract sequence context of variants
        mv ${projdir}/preprocess/processed/*_precall.bam ${projdir}/preprocess/ 2> /dev/null
        if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            awk '{print $1"\t"$2"\t"$3}' ${pop}_6x_rd${minRD_6x}_maf${maf}_dose.txt > snplist.txt &&
            awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
            awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
            for wint in {100,200}; do
              awk -v wint="$wint" '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
              awk -v wint="$wint" '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
              wait
            done
            awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
            mv snplist_haps.tmp snplist_haps.txt &&
            rm snplist_round.txt &&
            wait
          fi
          if [[ -n "$RE1" ]]; then
            mkdir -p seq_context_${RE1} &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' | \
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
          fi
          if [[ -n "$RE2" ]]; then
            mkdir -p seq_context_${RE2}
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' | \
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
          fi
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            mkdir consensus_seq_context &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              cat ./seq_context_"${RE1}"/"${i%.f*}"* ./seq_context_"${RE2}"/"${i%.f*}"* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
              maxlen=$(awk '{print $2}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}') &&
              awk -v maxlen="$maxlen" '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
              mv ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta 2> /dev/null &&
              while IFS="" read -r aln || [ -n "$aln" ]; do
    						awk -v aln="$aln" '$1 == aln {print}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
                bash $mafft --quiet --localpair --maxiterate 1000 ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp > ./consensus_seq_context/"${i%.f*}"_seqcontext.msf &&
    						bash $consambig -sequence ./consensus_seq_context/"${i%.f*}"_seqcontext.msf -outseq ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp 2> /dev/null &&
    						awk -v aln="${i%.f*}~~~${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp >> ./consensus_seq_context/"${i%.f*}"_seqcontext.cons &&
    						wait
    					done < <(awk '{print $1}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | sort | uniq )
              rm ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf 2> /dev/null &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
            wait
            for nn in ./consensus_seq_context/*.cons; do
              for run in $(seq 1 10); do
                awk '{gsub(/n$|N$/,"");}1' "$nn" 2> /dev/null > "${nn}".tmp && mv "${nn}".tmp "${nn}" 2> /dev/null &&
                wait
              done
            done
            wait
            seqid=$(cat ./consensus_seq_context/*.cons 2> /dev/null | grep '>' | awk '{gsub(/>/,""); gsub(/~~~Chr/,"\t~~~Chr");}1' | awk '{print $2}' | sort | uniq)
            wait
            mkdir -p ./consensus_seq_context/seqid_combine_samples
            wait
            for splitseqid in $seqid; do
              filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}" && filename="${filename/~~~}" &&
              cat ./consensus_seq_context/*.cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
              awk -v pat="$splitseqid" '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/"${filename}" &&
              minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{print length}' | sort -n | head -n 1) &&
              awk -v minlen="$minlen" '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{gsub(/n/,"N",$2); gsub(/~~~/,"\t");}1'| awk '{print $1"\n"$3}' > ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp &&
              mv ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp ./consensus_seq_context/seqid_combine_samples/"${filename}" 2> /dev/null &&
              wait
            done
            wait
            mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/ 2> /dev/null &&
            mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*.cons ./consensus_seq_context/sample_consensus_seqs 2> /dev/null &&
            mv seq_context_TGCAT ./consensus_seq_context/ 2> /dev/null &&
            mv seq_context_CATG ./consensus_seq_context/ 2> /dev/null &&
            wait
          fi
        fi
    		find . -type f -empty -delete
    	fi
    	wait
    	if [ -d "${projdir}/snpfilter/8x" ]; then
    		cd ${projdir}/snpfilter &&
    		cp -r 8x 8x_diversity_gmiss"${gmiss}"_smiss"${smiss}" &&
    		cd ./8x_diversity_gmiss"${gmiss}"_smiss"${smiss}" &&
    		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_8x.R $pop $gmiss $smiss $minRD_8x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome_number $biallelic $select_samples 2> /dev/null &&
    		rm ${pop}_8x_rawRD${minRD_8x}_DP_GT.txt ${pop}_8x_DP_GT.txt ${pop}_8x_rd${minRD_8x}.txt 2> /dev/null &&
    		mkdir visualizations && mv ./*.tiff ./visualizations/ 2> /dev/null &&
        wait

        if [[ "$variant_intervals" == true ]]; then
          mkdir -p variant_intervals &&
          cd variant_intervals &&
          wait
          while IFS="" read -r p || [ -n "$p" ]; do
            sleep $((RANDOM % 2))
            cat ../${pop}_8x_rd${minRD_8x}_maf${maf}_dose.txt | grep $p | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
            paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
            awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
            rm variant_intervals_${p}.txt &&
            :> variant_intervals_${p}.txt &&
            grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
            awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
            chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep $p | awk '{gsub(/LN:/,""); print $3}') &&
            awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
            rm vbreak_break.txt variant_intervals_${p}.tmp &&
            wait
            for vbreak in vbreak_*; do
              awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
              awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt &&
              wait
            done
            sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
            mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
            rm vbreak_* &&
            wait
          done < <(awk 'NR>1{print $2}' ../${pop}_8x_rd${minRD_8x}_maf${maf}_dose.txt | sort | uniq)
          cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
          rm variant_intervals_*.txt && cd ../
          rm -rf variant_intervals &&
          wait
        fi

        # Extract sequence context of variants
        mv ${projdir}/preprocess/processed/*_precall.bam ${projdir}/preprocess/ 2> /dev/null
        if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            awk '{print $1"\t"$2"\t"$3}' ${pop}_8x_rd${minRD_8x}_maf${maf}_dose.txt > snplist.txt &&
            awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
            awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
            wait
            for wint in {100,200}; do
              awk -v wint="$wint" '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
              awk -v wint="$wint" '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
              wait
            done
            awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
            mv snplist_haps.tmp snplist_haps.txt &&
            rm snplist_round.txt &&
            wait
          fi
          if [[ -n "$RE1" ]]; then
            mkdir -p seq_context_${RE1} &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
              rm ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
          fi
          if [[ -n "$RE2" ]]; then
            mkdir -p seq_context_${RE2} &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
              rm ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
          fi
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            mkdir consensus_seq_context &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              cat ./seq_context_"${RE1}"/"${i%.f*}"* ./seq_context_"${RE2}"/"${i%.f*}"* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
              maxlen=$(awk '{print $2}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}') &&
              awk -v maxlen="$maxlen" '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
              mv ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta 2> /dev/null &&
              while IFS="" read -r aln || [ -n "$aln" ]; do
    						awk -v aln="$aln" '$1 == aln {print}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
                bash $mafft --quiet --localpair --maxiterate 1000 ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp > ./consensus_seq_context/"${i%.f*}"_seqcontext.msf &&
    						bash $consambig -sequence ./consensus_seq_context/"${i%.f*}"_seqcontext.msf -outseq ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp 2> /dev/null &&
    						awk -v aln="${i%.f*}~~~${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp >> ./consensus_seq_context/"${i%.f*}"_seqcontext.cons &&
    						wait
    					done < <(awk '{print $1}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | sort | uniq )
              rm ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf 2> /dev/null &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
            wait
            for nn in ./consensus_seq_context/*.cons; do
              for run in $(seq 1 10); do
                awk '{gsub(/n$|N$/,"");}1' "$nn" 2> /dev/null > "${nn}".tmp && mv "${nn}".tmp "${nn}" 2> /dev/null &&
                wait
              done
            done
            wait
            seqid=$(cat ./consensus_seq_context/*.cons 2> /dev/null | grep '>' | awk '{gsub(/>/,""); gsub(/~~~Chr/,"\t~~~Chr");}1' | awk '{print $2}' | sort | uniq)
            wait
            mkdir -p ./consensus_seq_context/seqid_combine_samples
            wait
            for splitseqid in $seqid; do
              filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}" && filename="${filename/~~~}" &&
              cat ./consensus_seq_context/*.cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
              awk -v pat="$splitseqid" '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/"${filename}" &&
              minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{print length}' | sort -n | head -n 1) &&
              awk -v minlen="$minlen" '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{gsub(/n/,"N",$2); gsub(/~~~/,"\t");}1'| awk '{print $1"\n"$3}' > ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp &&
              mv ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp ./consensus_seq_context/seqid_combine_samples/"${filename}" 2> /dev/null &&
              wait
            done
            wait
            mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/ 2> /dev/null &&
            mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*.cons ./consensus_seq_context/sample_consensus_seqs 2> /dev/null &&
            mv seq_context_TGCAT ./consensus_seq_context/ 2> /dev/null &&
            mv seq_context_CATG ./consensus_seq_context/ 2> /dev/null &&
            wait
          fi
        fi
    		find . -type f -empty -delete
    	fi
    	wait
    fi

    if [[ "$p1" ]]; then
    	if [ -d "${projdir}/snpfilter/2x" ]; then
    		cd ${projdir}/snpfilter &&
    		cp -r 2x 2x_biparental_gmiss"${gmiss}"_smiss"${smiss}" &&
    		cd ./2x_biparental_gmiss"${gmiss}"_smiss"${smiss}" &&
    		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_2x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_2x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$pseg" "$haplome_number" $biallelic $select_samples 2> /dev/null &&
    		rm "${pop}"_2x_rawRD"${minRD_2x}"_DP_GT.txt "${pop}"_2x_DP_GT.txt "${pop}"_2x_rd"${minRD_2x}".txt 2> /dev/null &&
    		mkdir visualizations && mv ./*.tiff ./visualizations/ 2> /dev/null &&
        wait

        if [[ "$variant_intervals" == true ]]; then
          mkdir -p variant_intervals &&
          cd variant_intervals &&
          wait
          while IFS="" read -r p || [ -n "$p" ]; do
            sleep $((RANDOM % 2))
            cat ../${pop}_2x_rd${minRD_2x}_noSDdose.txt | grep $p | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
            paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
            awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
            rm variant_intervals_${p}.txt &&
            :> variant_intervals_${p}.txt &&
            grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
            awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
            chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep $p | awk '{gsub(/LN:/,""); print $3}') &&
            awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
            rm vbreak_break.txt variant_intervals_${p}.tmp &&
            wait
            for vbreak in vbreak_*; do
              awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
              awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt &&
              wait
            done
            sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
            mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
            rm vbreak_* &&
            wait
          done < <(awk 'NR>1{print $2}' ../${pop}_2x_rd${minRD_2x}_noSDdose.txt | sort | uniq)
          cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
          rm variant_intervals_*.txt && cd ../
          rm -rf variant_intervals &&
          wait
        fi

        # Extract sequence context of variants
        mv ${projdir}/preprocess/processed/*_precall.bam ${projdir}/preprocess/ 2> /dev/null
        if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            awk '{print $1"\t"$2"\t"$3}' ${pop}_2x_rd${minRD_2x}_noSDdose.txt > snplist.txt &&
            awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
            awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
            for wint in {100,200}; do
              awk -v wint="$wint" '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
              awk -v wint="$wint" '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
              wait
            done
            awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
            mv snplist_haps.tmp snplist_haps.txt &&
            rm snplist_round.txt &&
            wait
          fi
          if [[ -n "$RE1" ]]; then
            mkdir -p seq_context_"${RE1}" &&
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
              rm ./seq_context_"${RE1}"/${i%.f*}_seqcontext.tmp 2> /dev/null &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
          fi
          if [[ -n "$RE2" ]]; then
            mkdir -p seq_context_"${RE2}" &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
              rm ./seq_context_"${RE2}"/${i%.f*}_seqcontext.tmp 2> /dev/null &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
          fi
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            mkdir consensus_seq_context &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              cat ./seq_context_"${RE1}"/"${i%.f*}"* ./seq_context_"${RE2}"/"${i%.f*}"* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
              maxlen=$(awk '{print $2}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}') &&
              awk -v maxlen="$maxlen" '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
              mv ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta 2> /dev/null &&
              while IFS="" read -r aln || [ -n "$aln" ]; do
    						awk -v aln="$aln" '$1 == aln {print}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
                bash $mafft --quiet --localpair --maxiterate 1000 ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp > ./consensus_seq_context/"${i%.f*}"_seqcontext.msf &&
    						bash $consambig -sequence ./consensus_seq_context/"${i%.f*}"_seqcontext.msf -outseq ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp 2> /dev/null &&
    						awk -v aln="${i%.f*}~~~${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp >> ./consensus_seq_context/"${i%.f*}"_seqcontext.cons &&
    						wait
    					done < <(awk '{print $1}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | sort | uniq )
              rm ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf 2> /dev/null &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
            wait
            for nn in ./consensus_seq_context/*.cons; do
              for run in $(seq 1 10); do
                awk '{gsub(/n$|N$/,"");}1' "$nn" 2> /dev/null > "${nn}".tmp && mv "${nn}".tmp "${nn}" 2> /dev/null &&
                wait
              done
            done
            wait
            seqid=$(cat ./consensus_seq_context/*.cons 2> /dev/null | grep '>' | awk '{gsub(/>/,""); gsub(/~~~Chr/,"\t~~~Chr");}1' | awk '{print $2}' | sort | uniq)
            wait
            mkdir -p ./consensus_seq_context/seqid_combine_samples
            wait
            for splitseqid in $seqid; do
              filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}" && filename="${filename/~~~}" &&
              cat ./consensus_seq_context/*.cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
              awk -v pat="$splitseqid" '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/"${filename}" &&
              minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{print length}' | sort -n | head -n 1) &&
              awk -v minlen="$minlen" '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{gsub(/n/,"N",$2); gsub(/~~~/,"\t");}1'| awk '{print $1"\n"$3}' > ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp &&
              mv ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp ./consensus_seq_context/seqid_combine_samples/"${filename}" 2> /dev/null &&
              wait
            done
            wait
            mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/ 2> /dev/null &&
            mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*.cons ./consensus_seq_context/sample_consensus_seqs 2> /dev/null &&
            mv seq_context_TGCAT ./consensus_seq_context/ 2> /dev/null &&
            mv seq_context_CATG ./consensus_seq_context/ 2> /dev/null &&
            wait
          fi
        fi
    		find . -type f -empty -delete
    	fi
    	wait
    	if [ -d "${projdir}/snpfilter/4x" ]; then
    		cd ${projdir}/snpfilter &&
    		cp -r 4x 4x_biparental_gmiss"${gmiss}"_smiss"${smiss}" &&
    		cd ./4x_biparental_gmiss"${gmiss}"_smiss"${smiss}" &&
    		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_4x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_4x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$pseg" "$haplome_number" $biallelic $select_samples 2> /dev/null &&
    		rm ${pop}_4x_rawRD${minRD_4x}_DP_GT.txt ${pop}_4x_DP_GT.txt ${pop}_4x_rd${minRD_4x}.txt 2> /dev/null &&
    		mkdir visualizations && mv ./*.tiff ./visualizations/ 2> /dev/null &&
        wait

        if [[ "$variant_intervals" == true ]]; then
          mkdir -p variant_intervals &&
          cd variant_intervals &&
          wait
          while IFS="" read -r p || [ -n "$p" ]; do
            sleep $((RANDOM % 2))
            cat ../${pop}_4x_rd${minRD_4x}_noSDdose.txt | grep $p | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
            paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
            awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
            rm variant_intervals_${p}.txt &&
            :> variant_intervals_${p}.txt &&
            grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
            awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
            chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep $p | awk '{gsub(/LN:/,""); print $3}') &&
            awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
            rm vbreak_break.txt variant_intervals_${p}.tmp &&
            wait
            for vbreak in vbreak_*; do
              awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
              awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt &&
              wait
            done
            sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
            mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
            rm vbreak_* &&
            wait
          done < <(awk 'NR>1{print $2}' ../${pop}_4x_rd${minRD_4x}_noSDdose.txt | sort | uniq)
          cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
          rm variant_intervals_*.txt && cd ../
          rm -rf variant_intervals &&
          wait
        fi

        # Extract sequence context of variants
        mv ${projdir}/preprocess/processed/*_precall.bam ${projdir}/preprocess/ 2> /dev/null
        if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            awk '{print $1"\t"$2"\t"$3}' ${pop}_4x_rd${minRD_4x}_noSDdose.txt > snplist.txt &&
            awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
            awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
            for wint in {100,200}; do
              awk -v wint="$wint" '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
              awk -v wint="$wint" '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
              wait
            done
            awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
            mv snplist_haps.tmp snplist_haps.txt &&
            rm snplist_round.txt &&
            wait
          fi
          if [[ -n "$RE1" ]]; then
            mkdir -p seq_context_${RE1} &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
              rm ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
          fi
          if [[ -n "$RE2" ]]; then
            mkdir -p seq_context_${RE2} &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
              rm ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
          fi
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            mkdir consensus_seq_context &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              cat ./seq_context_"${RE1}"/"${i%.f*}"* ./seq_context_"${RE2}"/"${i%.f*}"* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
              maxlen=$(awk '{print $2}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}') &&
              awk -v maxlen="$maxlen" '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
              mv ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta 2> /dev/null &&
              while IFS="" read -r aln || [ -n "$aln" ]; do
    						awk -v aln="$aln" '$1 == aln {print}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
                bash $mafft --quiet --localpair --maxiterate 1000 ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp > ./consensus_seq_context/"${i%.f*}"_seqcontext.msf &&
    						bash $consambig -sequence ./consensus_seq_context/"${i%.f*}"_seqcontext.msf -outseq ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp 2> /dev/null &&
    						awk -v aln="${i%.f*}~~~${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp >> ./consensus_seq_context/"${i%.f*}"_seqcontext.cons &&
    						wait
    					done < <(awk '{print $1}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | sort | uniq )
              rm ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf 2> /dev/null &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
            wait
            for nn in ./consensus_seq_context/*.cons; do
              for run in $(seq 1 10); do
                awk '{gsub(/n$|N$/,"");}1' "$nn" 2> /dev/null > "${nn}".tmp && mv "${nn}".tmp "${nn}" 2> /dev/null &&
                wait
              done
            done
            wait
            seqid=$(cat ./consensus_seq_context/*.cons 2> /dev/null | grep '>' | awk '{gsub(/>/,""); gsub(/~~~Chr/,"\t~~~Chr");}1' | awk '{print $2}' | sort | uniq)
            wait
            mkdir -p ./consensus_seq_context/seqid_combine_samples
            wait
            for splitseqid in $seqid; do
              filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}" && filename="${filename/~~~}" &&
              cat ./consensus_seq_context/*.cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
              awk -v pat="$splitseqid" '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/"${filename}" &&
              minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{print length}' | sort -n | head -n 1) &&
              awk -v minlen="$minlen" '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{gsub(/n/,"N",$2); gsub(/~~~/,"\t");}1'| awk '{print $1"\n"$3}' > ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp &&
              mv ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp ./consensus_seq_context/seqid_combine_samples/"${filename}" 2> /dev/null &&
              wait
            done
            wait
            mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/ 2> /dev/null &&
            mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*.cons ./consensus_seq_context/sample_consensus_seqs 2> /dev/null &&
            mv seq_context_TGCAT ./consensus_seq_context/ 2> /dev/null &&
            mv seq_context_CATG ./consensus_seq_context/ 2> /dev/null &&
            wait
          fi
        fi
    		find . -type f -empty -delete
    	fi
    	wait
    	if [ -d "${projdir}/snpfilter/6x" ]; then
    		cd ${projdir}/snpfilter &&
    		cp -r 6x 6x_biparental_gmiss"${gmiss}"_smiss"${smiss}" &&
    		cd ./6x_biparental_gmiss"${gmiss}"_smiss"${smiss}" &&
    		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_6x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_6x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$pseg" "$haplome_number" $biallelic $select_samples 2> /dev/null &&
    		rm ${pop}_6x_rawRD${minRD_6x}_DP_GT.txt ${pop}_6x_DP_GT.txt ${pop}_6x_rd${minRD_6x}.txt 2> /dev/null &&
    		mkdir visualizations && mv ./*.tiff ./visualizations/ 2> /dev/null &&
        wait

        if [[ "$variant_intervals" == true ]]; then
          mkdir -p variant_intervals &&
          cd variant_intervals &&
          wait
          while IFS="" read -r p || [ -n "$p" ]; do
            sleep $((RANDOM % 2))
            grep $p | ../${pop}_6x_rd${minRD_6x}_noSDdose.txt | grep $p | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
            paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
            awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
            rm variant_intervals_${p}.txt &&
            :> variant_intervals_${p}.txt &&
            grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
            awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
            chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep $p | awk '{gsub(/LN:/,""); print $3}') &&
            awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
            rm vbreak_break.txt variant_intervals_${p}.tmp &&
            for vbreak in vbreak_*; do
              awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
              awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt &&
              wait
            done
            sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
            mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
            rm vbreak_* &&
            wait
          done < <(awk 'NR>1{print $2}' ../${pop}_6x_rd${minRD_6x}_noSDdose.txt | sort | uniq)
          cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
          rm variant_intervals_*.txt && cd ../
          rm -rf variant_intervals &&
          wait
        fi

        # Extract sequence context of variants
        mv ${projdir}/preprocess/processed/*_precall.bam ${projdir}/preprocess/ 2> /dev/null
        if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            awk '{print $1"\t"$2"\t"$3}' ${pop}_6x_rd${minRD_6x}_noSDdose.txt > snplist.txt &&
            awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
            awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
            wait
            for wint in {100,200}; do
              awk -v wint="$wint" '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
              awk -v wint="$wint" '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
              wait
            done
            awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
            mv snplist_haps.tmp snplist_haps.txt &&
            rm snplist_round.txt &&
            wait
          fi
          if [[ -n "$RE1" ]]; then
            mkdir -p seq_context_${RE1}
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
              rm ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
          fi
          if [[ -n "$RE2" ]]; then
            mkdir -p seq_context_${RE2} &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
              rm ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
          fi
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            mkdir consensus_seq_context &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              cat ./seq_context_"${RE1}"/"${i%.f*}"* ./seq_context_"${RE2}"/"${i%.f*}"* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
              maxlen=$(awk '{print $2}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}') &&
              awk -v maxlen="$maxlen" '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
              mv ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta 2> /dev/null &&
              while IFS="" read -r aln || [ -n "$aln" ]; do
    						awk -v aln="$aln" '$1 == aln {print}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
                bash $mafft --quiet --localpair --maxiterate 1000 ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp > ./consensus_seq_context/"${i%.f*}"_seqcontext.msf &&
    						bash $consambig -sequence ./consensus_seq_context/"${i%.f*}"_seqcontext.msf -outseq ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp 2> /dev/null &&
    						awk -v aln="${i%.f*}~~~${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp >> ./consensus_seq_context/"${i%.f*}"_seqcontext.cons &&
    						wait
    					done < <(awk '{print $1}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | sort | uniq )
              rm ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf 2> /dev/null &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
            wait
            for nn in ./consensus_seq_context/*.cons; do
              for run in $(seq 1 10); do
                awk '{gsub(/n$|N$/,"");}1' "$nn" 2> /dev/null > "${nn}".tmp && mv "${nn}".tmp "${nn}" 2> /dev/null &&
                wait
              done
            done
            wait
            seqid=$(cat ./consensus_seq_context/*.cons 2> /dev/null | grep '>' | awk '{gsub(/>/,""); gsub(/~~~Chr/,"\t~~~Chr");}1' | awk '{print $2}' | sort | uniq)
            wait
            mkdir -p ./consensus_seq_context/seqid_combine_samples
            wait
            for splitseqid in $seqid; do
              filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}" && filename="${filename/~~~}" &&
              cat ./consensus_seq_context/*.cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
              awk -v pat="$splitseqid" '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/"${filename}" &&
              minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{print length}' | sort -n | head -n 1) &&
              awk -v minlen="$minlen" '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{gsub(/n/,"N",$2); gsub(/~~~/,"\t");}1'| awk '{print $1"\n"$3}' > ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp &&
              mv ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp ./consensus_seq_context/seqid_combine_samples/"${filename}" 2> /dev/null &&
              wait
            done
            wait
            mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/ 2> /dev/null &&
            mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*.cons ./consensus_seq_context/sample_consensus_seqs 2> /dev/null &&
            mv seq_context_TGCAT ./consensus_seq_context/ 2> /dev/null &&
            mv seq_context_CATG ./consensus_seq_context/ 2> /dev/null &&
            wait
          fi
        fi
    		find . -type f -empty -delete
    	fi
    	wait
    	if [ -d "${projdir}/snpfilter/8x" ]; then
    		cd ${projdir}/snpfilter &&
    		cp -r 8x 8x_biparental_gmiss"${gmiss}"_smiss"${smiss}" &&
    		cd ./8x_biparental_gmiss"${gmiss}"_smiss"${smiss}" &&
    		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_8x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_8x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$pseg" "$haplome_number" $biallelic $select_samples 2> /dev/null &&
    		rm ${pop}_8x_rawRD${minRD_8x}_DP_GT.txt ${pop}_8x_DP_GT.txt ${pop}_8x_rd${minRD_8x}.txt 2> /dev/null &&
    		mkdir visualizations && mv ./*.tiff ./visualizations/ 2> /dev/null &&
        wait

        if [[ "$variant_intervals" == true ]]; then
          mkdir -p variant_intervals &&
          cd variant_intervals &&
          wait
          while IFS="" read -r p || [ -n "$p" ]; do
            sleep $((RANDOM % 2))
            cat ../${pop}_8x_rd${minRD_8x}_noSDdose.txt | grep $p | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
            paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
            awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
            rm variant_intervals_${p}.txt &&
            :> variant_intervals_${p}.txt &&
            grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
            awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
            chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep $p | awk '{gsub(/LN:/,""); print $3}') &&
            awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
            rm vbreak_break.txt variant_intervals_${p}.tmp &&
            wait
            for vbreak in vbreak_*; do
              awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
              awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt &&
              wait
            done
            sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
            mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
            rm vbreak_* &&
            wait
          done < <(awk 'NR>1{print $2}' ../${pop}_8x_rd${minRD_8x}_noSDdose.txt | sort | uniq)
          cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
          rm variant_intervals_*.txt && cd ../
          rm -rf variant_intervals &&
          wait
        fi

        # Extract sequence context of variants
        mv ${projdir}/preprocess/processed/*_precall.bam ${projdir}/preprocess/ 2> /dev/null
        if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            awk '{print $1"\t"$2"\t"$3}' ${pop}_8x_rd${minRD_8x}_noSDdose.txt > snplist.txt &&
            awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
            awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
            wait
            for wint in {100,200}; do
              awk -v wint="$wint" '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
              awk -v wint="$wint" '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
              wait
            done
            awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
            mv snplist_haps.tmp snplist_haps.txt &&
            rm snplist_round.txt &&
            wait
          fi
          if [[ -n "$RE1" ]]; then
            mkdir -p seq_context_${RE1} &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
              rm ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
          fi
          if [[ -n "$RE2" ]]; then
            mkdir -p seq_context_${RE2} &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
              awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
              rm ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
          fi
          if [[ -n "$RE1" ]] || [[ -n "$RE2" ]]; then
            mkdir consensus_seq_context &&
            wait
            while IFS="" read -r i || [ -n "$i" ]; do
              sleep $((RANDOM % 2))
              cat ./seq_context_"${RE1}"/"${i%.f*}"* ./seq_context_"${RE2}"/"${i%.f*}"* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
              maxlen=$(awk '{print $2}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}') &&
              awk -v maxlen="$maxlen" '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
              mv ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta 2> /dev/null &&
              while IFS="" read -r aln || [ -n "$aln" ]; do
    						awk -v aln="$aln" '$1 == aln {print}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp &&
                bash $mafft --quiet --localpair --maxiterate 1000 ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp > ./consensus_seq_context/"${i%.f*}"_seqcontext.msf &&
    						bash $consambig -sequence ./consensus_seq_context/"${i%.f*}"_seqcontext.msf -outseq ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp 2> /dev/null &&
    						awk -v aln="${i%.f*}~~~${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp >> ./consensus_seq_context/"${i%.f*}"_seqcontext.cons &&
    						wait
    					done < <(awk '{print $1}' ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta | sort | uniq )
              rm ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf 2> /dev/null &&
              wait
            done < <(cat ${projdir}/fetch_samples_seq.txt)
            wait
            for nn in ./consensus_seq_context/*.cons; do
              for run in $(seq 1 10); do
                awk '{gsub(/n$|N$/,"");}1' "$nn" 2> /dev/null > "${nn}".tmp && mv "${nn}".tmp "${nn}" 2> /dev/null &&
                wait
              done
            done
            wait
            seqid=$(cat ./consensus_seq_context/*.cons 2> /dev/null | grep '>' | awk '{gsub(/>/,""); gsub(/~~~Chr/,"\t~~~Chr");}1' | awk '{print $2}' | sort | uniq)
            wait
            mkdir -p ./consensus_seq_context/seqid_combine_samples
            wait
            for splitseqid in $seqid; do
              filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}" && filename="${filename/~~~}" &&
              cat ./consensus_seq_context/*.cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
              awk -v pat="$splitseqid" '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/"${filename}" &&
              minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{print length}' | sort -n | head -n 1) &&
              awk -v minlen="$minlen" '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/"${filename}" | awk '{gsub(/n/,"N",$2); gsub(/~~~/,"\t");}1'| awk '{print $1"\n"$3}' > ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp &&
              mv ./consensus_seq_context/seqid_combine_samples/"${filename}".tmp ./consensus_seq_context/seqid_combine_samples/"${filename}" 2> /dev/null &&
              wait
            done
            wait
            mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/ 2> /dev/null &&
            mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*.cons ./consensus_seq_context/sample_consensus_seqs 2> /dev/null &&
            mv seq_context_TGCAT ./consensus_seq_context/ 2> /dev/null &&
            mv seq_context_CATG ./consensus_seq_context/ 2> /dev/null &&
            wait
          fi
        fi
    		find . -type f -empty -delete
    	fi
    fi
  done
  wait
done
wait
mv ${projdir}/preprocess/*_precall.bam ${projdir}/preprocess/processed/ 2> /dev/null

cd ${projdir}/snpfilter/
rm ${projdir}/fetch_samples_seq.txt
find . -type f -empty -delete
find . -type d -empty -delete
for snpfilter_dir in */; do
	cd $snpfilter_dir &&
  mkdir -p visualizations &&
  for tiff in ./*.tiff; do if [ -f "$tiff" ]; then mv ./*.tiff ./visualizations/; fi; done
	smmiss_thresh=${snpfilter_dir#*smiss} &&
	smmiss_thresh=${smmiss_thresh%*/} &&
	smmiss_thresh=$(echo "$smmiss_thresh * 100" | bc) &&
	awk -v smisst=$smmiss_thresh '(NR>1) && ($2 <= smisst)' sample_missing_rate* > retained_samples.txt 2> /dev/null &&
	cd ../ &&
  wait
done
wait
wc -l *gmiss*/*dose.txt | awk '{print $2"\t"$1-1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > gmiss_smiss_titration.txt &&
wc -l *gmiss*/eliminated* | awk '{print $2"\t"$1-1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > eliminated_samples.txt &&
wc -l *gmiss*/retained_samples.txt | awk '{print $2"\t"$1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > retained_samples.txt &&
echo -e "gmiss_smiss_thresholds\t#_of_retained_samples\t#_of_SNPs\t#_of_eliminated_samples\n----------------------\t-----------------------\t---------\t-----------------------" > summary_precall.txt &&
awk 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,"\t",a[$1],"\t",$2}' gmiss_smiss_titration.txt eliminated_samples.txt | \
awk 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,"\t",a[$1],"\t",$2,"\t",$3}' retained_samples.txt - | \
cat summary_precall.txt - > gmiss_smiss.txt &&
rm gmiss_smiss_titration.txt eliminated_samples.txt retained_samples.txt summary_precall.txt 2> /dev/null &&
wait
for snpfilter_dir in */; do
if [[ "$(ls ${snpfilter_dir}/*dose.txt 2> /dev/null | wc -l)" -eq 0 ]]; then
rm -rf ${snpfilter_dir} &&
wait
fi
done
wait
ls ./*/*maf*.txt 2> /dev/null | grep -v 'maf0.txt' | grep -v 'dose' | grep -v 'binary' | xargs rm 2> /dev/null &&
ls ./*/*_plusSD.txt 2> /dev/null | xargs rm 2> /dev/null &&
ls ./*/*SD_1_G*G*.txt 2> /dev/null | xargs rm 2> /dev/null &&
wait



cd "$projdir"/snpfilter
export n="${ref1%.f*}"
{
  for snpfilter_dir in */; do
    if [ -d "${projdir}/snpfilter/${snpfilter_dir}" ]; then
      cd "${projdir}/snpfilter/${snpfilter_dir}" &&
      export ploidydir=${snpfilter_dir:0:1} &&

  		for i in $(ls *dose.txt); do
        if [[ "$filtered_vcf" == "true" ]]; then
    			export ARselect=$( echo $i | awk -v FS='rd' '{print $1}') &&
    			export ARfile=$(ls ../../snpcall/${ARselect}*AR.txt 2> /dev/null) &&
          export arr=$(cat ${projdir}/samples_list_node_* | awk '{gsub(/.fastq/,"\t.fastq");gsub(/.fq/,"\t.fq");}1' | awk '{print $1}' | tr '\n' ',' | awk '{gsub(/\t/,",");gsub(/ /,",");gsub(/^,/,"");gsub(/,$/,"");}1') &&
          export arr2=$(grep "CHROM" $i | awk '{$1=$2=$3=$4=$5=""}1' | tr -s ' ' | awk '{gsub(/ pvalue/,"");}1' | awk '{gsub(/\t/,",");gsub(/ /,",");gsub(/^,/,"");gsub(/,$/,"");}1') &&
          export darr=$(echo ${arr[@]},${arr2[@]} | tr ',' '\n' | sort | uniq -u | tr '\n' ',' | awk '{gsub(/\t/,",");gsub(/ /,",");gsub(/^,/,"");gsub(/,$/,"");}1') &&
          export darr2=$(echo ${arr[@]},${arr2[@]} | tr ',' '\n' | sort | uniq -u | awk '!/^$/' | awk '{print $1"_AR"}' | tr '\n' ',' | awk '{gsub(/\t/,",");gsub(/ /,",");gsub(/^,/,"");gsub(/,$/,"");}1') &&
          echo $darr | tr ',' '\n' > darr.txt &&
          echo $darr2 | tr ',' '\n' > darr2.txt &&
          if [[ -s darr.txt ]]; then printf "999_999_999\n" > darr.txt; fi
          wait
          if [[ -s darr2.txt ]]; then printf "999_999_999\n" > darr2.txt; fi
          wait

    			Rscript "${GBSapp_dir}"/scripts/R/heterozygote_vs_allele_ratio.R "$i" "$ARfile" "${ploidydir}x" "1" "darr2.txt" "${GBSapp_dir}/tools/R" 2> /dev/null &&
    			awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$2,$3]){print b}}' $ARfile $i | awk '{gsub(/NA/,"na"); print $1"_"$2"\t"$0}' | \
    			awk -v pat1="${n}_Chr" -v pat2="${n}_chr" '{gsub(pat1,"Chr");gsub(pat2,"Chr");gsub(/CHROM_POS/,"SNP");}1' > ${i%.txt}_AR_metric.txt &&
          wait

          if test ! -f $projdir/split_done.txt; then
            export vcfdose=${i%_rd*}; vcfdose=${vcfdose#*_} &&
            for gunv in $projdir/snpcall/*x.vcf.gz; do
              gunzip $gunv &&
              wait
            done
            wait
            for split in $projdir/snpcall/*${vcfdose}.vcf; do
              refg=${split%_*} && refg=${refg##*/} && refg=${refg#*_} && refg=${refg%%_*}.fasta
            	$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" IndexFeatureFile -I $split &&
            	$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$refg -V $split -O ${split%.vcf}_split.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac &&
            	gzip $split &&
              $bcftools view -I ${split%.vcf}_split.vcf -O z -o ${split%.vcf}_split.vcf.gz &&
              $bcftools index ${split%.vcf}_split.vcf.gz &&
              :> $projdir/split_done.txt &&
              wait
            done
            wait
          fi
          if [[ "$(ls $projdir/snpcall/*_split.vcf.gz  2> /dev/null | wc -l)" -gt 1 ]]; then
            $bcftools merge $projdir/snpcall/*_split.vcf.gz --force-samples -m all > ${i%rd*}split.vcf &&
            gzip ${i%rd*}split.vcf &&
            wait
          else
            for k in $projdir/snpcall/*_split.vcf.gz; do cp $k ./${i%rd*}split.vcf.gz; wait; done
          fi
          wait

          for i in *dose*; do
            awk -v pat1="${n}_Chr" -v pat2="${n}_chr" '{gsub(pat1,"Chr");gsub(pat2,"Chr"); print $0}' $i > ${i%.txt}_hold.txt &&
            mv ${i%.txt}_hold.txt $i &&
            wait
          done
          wait
          zcat ${i%rd*}split.vcf.gz | grep '^#' | awk -v pat1="${n}_Chr" -v pat2="${n}_chr" '{gsub(pat1,"Chr");gsub(pat2,"Chr");gsub(/chr/,"Chr");}1' > ${i%.txt}_header.vcf &&
          awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$2,$3]){print b}}' <(zcat ${i%rd*}split.vcf.gz | grep -v '^#' | awk -v pat1="${n}_Chr" -v pat2="${n}_chr" '{gsub(pat1,"Chr");gsub(pat2,"Chr");gsub(/chr/,"Chr");}1') <(awk '!seen[$0] {print} {++seen[$0]}' $i) | \
          sort -Vk1,1 -Vk2,2 | cat ${i%.txt}_header.vcf - > ${i%dose.txt}.vcf &&
          rm ${i%.txt}_header.vcf 2> /dev/null &&
          rm ${i%rd*}split.vcf.gz 2> /dev/null &&

          # grep '^##' ${i%dose.txt}.vcf > ${i%.txt}_header.vcf &&
          # grep -v '^##' ${i%dose.txt}.vcf | awk '{gsub(/#CHROM/,"CHROM");}1' > ${i%dose.txt}_tmp.vcf &&
          # rm ${i%dose.txt}.vcf &&
          # {
          #   Rscript "${GBSapp_dir}"/scripts/R/recode_vcf.R "${i%dose.txt}_tmp.vcf" "$i" "${i%.txt}_AR_metric.txt" "${ploidy}x" "darr.txt" "$filter_ExcHet" "${GBSapp_dir}/tools/R"
          # } & PID_AR=$!
          # wait $PID_AR
          # cat <(grep '^##' ${i%.txt}_header.vcf) <(awk '{gsub(/CHROM/,"#CHROM");}1' dose_temp.vcf) | \
          # awk '{gsub(/0,0:.:.:./,"0,0:.:."); gsub(/0,0,0,0:.:.:./,"0,0:.:."); gsub(/0,0,0,0,0,0:.:.:./,"0,0:.:."); gsub(/0,0,0,0,0,0,0,0:.:.:./,"0,0:.:."); }1' 2> /dev/null > ${i%dose.txt}.vcf &&
          # mv AR_temp.txt ${i%.txt}_AR_metric.txt &&
          # rm ${i%.txt}_header.vcf 2> /dev/null &&
          # rm ${i%dose.txt}_tmp.vcf dose_temp.vcf 2> /dev/null &&
          # wait
          #
    			# grep -v 'CHROM' ${i%.txt}_AR_metric.txt | awk -F'\t' '{$1=$2=$3=$4=$5=""}1' | awk -F'\t' '{gsub(/na/,"0");}1' | \
    			# awk '{gsub(/\t/," ");}1' | awk '{gsub(/-/,"");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/0,/,",");}1' | awk '{gsub(/,0$/,",");}1' | \
    			# awk -F',' -v OFS=',' -v OFMT='%0.3g' '{s=0; numFields=0; for(i=2; i<=NF;i++){if(length($i)){s+=$i; numFields++}} print (numFields ? s/numFields : 0)}' | \
    			# cat <(printf "Allele_ratio_mean\n") - | paste <(awk '{print $1"\t"$2"\t"$3}' ${i%.txt}_AR_metric.txt) - > ${i%.txt}_AR_mean.txt &&
          gzip "${i%dose.txt}".vcf 2> /dev/null &&
          if ls *_.vcf.gz 1> /dev/null 2>&1; then mv ${i%dose.txt}.vcf.gz ${i%_dose.txt}.vcf.gz; fi
          wait

        else
          for i in *dose*; do
            awk -v pat1="${n}_Chr" -v pat2="${n}_chr" '{gsub(pat1,"Chr");gsub(pat2,"Chr"); print $0}' $i > ${i%.txt}_hold.txt &&
            mv ${i%.txt}_hold.txt $i &&
            wait
          done
        fi
        wait

  			if [[ "$ploidy" -le 2 ]]; then
          Rscript "${GBSapp_dir}"/scripts/R/hapmap_format.R "$i" "${GBSapp_dir}/tools/R" 2> /dev/null &&
          mv outfile.hmp.txt "${i%dose.txt}.hmp.txt" 2> /dev/null &&
          wait
        fi

        mv ${i%.txt}_AR_metric.txt ${i%dose.txt}_AR_metric.txt 2> /dev/null &&
        mv ${i%.txt}_AR_mean.txt ${i%dose.txt}_AR_mean.txt 2> /dev/null &&
        find . -type f -empty -delete
        rm *AR*.txt 2> /dev/null &&
        wait
  		done
  		wait

      # for i in ./*/*vcf.gz; do zcat $i | awk 'BEGIN {OFS = FS} {sub(/,.*/,"",$5); print}' | gzip > ${i%.vcf.gz}_1alt_allele.vcf.gz; done
  		# cd unique_mapped
  		# for i in $(ls *dose_unique_mapped.txt 2> /dev/null); do
  		# 	ARselect=${i%rd*}
  		# 	ARfile=$(ls ../../../snpcall/${ARselect}*AR.txt 2> /dev/null)
  		# 	Rscript "${GBSapp_dir}"/scripts/R/heterozygote_vs_allele_ratio_uniqfiltered.R "$i" "$ARfile" "${ploidy}x" "1" "${GBSapp_dir}/tools/R"
  		# 	awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$2,$3]){print b}}' $ARfile $i | awk '{gsub(/NA/,"na"); print $1"_"$2"\t"$0}' | \
  		# 	awk -v n="$n" '{gsub(n,""); gsub(/CHROM_POS/,"SNP");}1' > ${i%.txt}_AR_metric.txt
      #
  		# 	vcfdose=${i%_rd*}; vcfdose=${vcfdose#*_}
  		# 	zcat ../../../snpcall/*${vcfdose}.vcf.gz | grep '^#' > ${i%.txt}.vcf
  		# 	awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$2,$3]){print b}}' <(zcat ../../../snpcall/*${vcfdose}.vcf.gz) $i >> ${i%.txt}.vcf
  		# 	arr=$(grep "CHROM" $i | awk '{$1=$2=$3=$4=$5=""}1' | tr -s ' ' | awk '{gsub(/ pvalue/,"");}1' | awk '{gsub(/\t/,",");gsub(/ /,",");gsub(/^,/,"");}1')
  		# 	$bcftools view -s "$arr2" ${i%.txt}.vcf > tmp.vcf && mv tmp.vcf ${i%.txt}.vcf
      #
  		# 	grep -v '^##' ${i%.txt}.vcf | awk '{gsub(/#CHROM/,"CHROM");}1' > ${i%.txt}_tmp.vcf
  		# 	Rscript "${GBSapp_dir}"/scripts/R/recode_vcf.R "${i%.txt}_tmp.vcf" "$i" "${i%.txt}_AR_metric.txt" "${ploidy}x" "${GBSapp_dir}/tools/R"
  		# 	grep '^##' ${i%.txt}.vcf | cat - <(awk '{gsub(/CHROM/,"#CHROM");}1' dose_temp.vcf) > ${i%.txt}.vcf
  		# 	mv AR_temp.txt ${i%.txt}_AR_metric.txt
  		# 	rm ${i%.txt}_tmp.vcf dose_temp.vcf
  		# 	gzip ${i%.txt}.vcf
      #
  		# 	grep -v 'CHROM' ${i%.txt}_AR_metric.txt | awk -F'\t' '{$1=$2=$3=$4=$5=""}1' | awk -F'\t' '{gsub(/na/,"0");}1' | \
  		# 	awk '{gsub(/\t/," ");}1' | awk '{gsub(/-/,"");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/0,/,",");}1' | awk '{gsub(/,0$/,",");}1' | \
  		# 	awk -F',' -v OFS=',' -v OFMT='%0.3g' '{s=0; numFields=0; for(i=2; i<=NF;i++){if(length($i)){s+=$i; numFields++}} print (numFields ? s/numFields : 0)}' | \
  		# 	cat <(printf "Allele_ratio_mean\n") - | paste <(awk '{print $1"\t"$2"\t"$3}' ${i%.txt}_AR_metric.txt) - > ${i%.txt}_AR_mean.txt
      #
  		# 	Rscript "${GBSapp_dir}"/scripts/R/hapmap_format.R "$i" "${GBSapp_dir}/tools/R"
      #   mv outfile.hmp.txt "${i%.txt}.hmp.txt"
      #   wait
  		# done
  		# wait
      #
  		# for i in $(ls *dose_multi_mapped.txt 2> /dev/null); do
  		# 	ARselect=${i%rd*}
  		# 	ARfile=$(ls ../../../snpcall/${ARselect}*AR.txt 2> /dev/null)
  		# 	Rscript "${GBSapp_dir}"/scripts/R/heterozygote_vs_allele_ratio_multifiltered.R "$i" "$ARfile" "${ploidy}x" "1" "${GBSapp_dir}/tools/R"
  		# 	awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$2,$3]){print b}}' $ARfile $i | awk '{gsub(/NA/,"na"); print $1"_"$2"\t"$0}' | \
  		# 	awk -v n="$n" '{gsub(n,""); gsub(/CHROM_POS/,"SNP");}1' > ${i%.txt}_AR_metric.txt
      #
  		# 	vcfdose=${i%_rd*}; vcfdose=${vcfdose#*_}
  		# 	zcat ../../../snpcall/*${vcfdose}.vcf.gz | grep '^#' > ${i%.txt}.vcf
  		# 	awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$2,$3]){print b}}' <(zcat ../../../snpcall/*${vcfdose}.vcf.gz) $i >> ${i%.txt}.vcf
  		# 	arr=$(grep "CHROM" $i | awk '{$1=$2=$3=$4=$5=""}1' | tr -s ' ' | awk '{gsub(/ pvalue/,"");}1' | awk '{gsub(/\t/,",");gsub(/ /,",");gsub(/^,/,"");}1')
  		# 	$bcftools view -s "$arr2" ${i%.txt}.vcf > tmp.vcf && mv tmp.vcf ${i%.txt}.vcf
      #
  		# 	grep -v '^##' ${i%.txt}.vcf | awk '{gsub(/#CHROM/,"CHROM");}1' > ${i%.txt}_tmp.vcf
  		# 	Rscript "${GBSapp_dir}"/scripts/R/recode_vcf.R "${i%.txt}_tmp.vcf" "$i" "${i%.txt}_AR_metric.txt" "${ploidy}x" "${GBSapp_dir}/tools/R"
  		# 	grep '^##' ${i%.txt}.vcf | cat - <(awk '{gsub(/CHROM/,"#CHROM");}1' dose_temp.vcf) > ${i%.txt}.vcf
  		# 	mv AR_temp.txt ${i%.txt}_AR_metric.txt
  		# 	rm ${i%.txt}_tmp.vcf dose_temp.vcf
  		# 	gzip ${i%.txt}.vcf
      #
  		# 	grep -v 'CHROM' ${i%.txt}_AR_metric.txt | awk -F'\t' '{$1=$2=$3=$4=$5=""}1' | awk -F'\t' '{gsub(/na/,"0");}1' | \
  		# 	awk '{gsub(/\t/," ");}1' | awk '{gsub(/-/,"");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/0,/,",");}1' | awk '{gsub(/,0$/,",");}1' | \
  		# 	awk -F',' -v OFS=',' -v OFMT='%0.3g' '{s=0; numFields=0; for(i=2; i<=NF;i++){if(length($i)){s+=$i; numFields++}} print (numFields ? s/numFields : 0)}' | \
  		# 	cat <(printf "Allele_ratio_mean\n") - | paste <(awk '{print $1"\t"$2"\t"$3}' ${i%.txt}_AR_metric.txt) - > ${i%.txt}_AR_mean.txt
      #
      #   Rscript "${GBSapp_dir}"/scripts/R/hapmap_format.R "$i" "${GBSapp_dir}/tools/R"
      #   mv outfile.hmp.txt "${i%.txt}.hmp.txt"
      #   wait
  		# done
  		# wait
      #
  		# for i in $(ls *dose* 2> /dev/null); do
  		# 	awk -v n="$n" '{gsub(n,""); print $0}' $i > ${i%.txt}_hold.txt
  		# 	mv ${i%.txt}_hold.txt $i
  		# done

  		rm *tmp.vcf *header.vcf darr.txt darr2.txt 2> /dev/null &&
      wait
  	fi
    wait
  done
  wait
} & PID_genvcf=$!
wait $PID_genvcf
wait
rm ${projdir}/snpcall/*_split.vcf* 2> /dev/null &&
rm $projdir/split_done.txt &&
wait


}
cd $projdir
if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
	if [ "$walkaway" == false ]; then
		echo -e "${magenta}- Do you want to perform SNP/variant filtering? ${white}\n"
		read -p "- y(YES) or n(NO) " -t 36000 -n 1 -r
		if [[ ! $REPLY =~ ^[Yy]$ ]]; then
			printf '\n'
			echo -e "${magenta}- skipping SNP/variant filtering ${white}\n"
		else
			printf '\n'
			echo -e "${magenta}- performing SNP/variant filtering ${white}\n"
			time main &>> log.out
		fi
	fi
	if [ "$walkaway" == true ]; then
		if [ "$snp_filtering" == 1 ]; then
			echo -e "${magenta}- performing SNP/variant filtering ${white}\n"
			time main &>> log.out
		else
			echo -e "${magenta}- skipping SNP/variant filtering ${white}\n"
		fi
	fi
fi


######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${yellow}\n- Visualizations for Genotype Accuracy and ploidy estimation\n${blue}##############################################################################${white}\n"
main () {
cd $projdir
cd snpfilter
for snpfilter_dir in */; do
	cd "$snpfilter_dir"
	dose=${snpfilter_dir%_gmiss*} && \
	dose=${dose%_*} && \
	awk 'BEGIN{OFS="\t"}{print $2,$3}' *dose_RefPrefix.txt > CHROM_POS.txt && \
	cat "$projdir"/snpcall/*"${dose}"*.vcf | grep -Fwf CHROM_POS.txt - | awk '{gsub(/#/,""); print $0}' > snp_allele_depth.txt && \
	cd "$projdir"/snpfilter
done
wait
for snpfilter_dir in */; do (
	cd "$snpfilter_dir"
	dose=${snpfilter_dir%_gmiss*} && \
	dose=${dose%_*} && \
	poptype=${snpfilter_dir%_*_*} && \
	poptype=${poptype//*_} && \
	getrd=$( echo *dose.txt | awk  'BEGIN{OFS=FS="_rd"};{print $2}' | awk  'BEGIN{OFS=FS="_"};{print $1}' ) && \
	if [[ $poptype == biparental ]]; then
		cd ../"$snpfilter_dir"
		Rscript "${GBSapp_dir}"/scripts/R/SNPaccuracy_biparental_ReadDepth.R $p1 $p2 $getrd "${GBSapp_dir}"/tools/R 2> /dev/null &&
		wait
	fi
	if [[ $poptype == diversity ]]; then
		cd ../"$snpfilter_dir"
		Rscript "${GBSapp_dir}"/scripts/R/SNPaccuracy_diversity_ReadDepth.R $getrd "${GBSapp_dir}"/tools/R 2> /dev/null &&
		wait
	fi
	wait
	rm CHROM_POS.txt snp_allele_depth.txt && \
	cd "$projdir"/snpfilter ) &
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		wait
	fi
done
wait
cd "$projdir"/snpfilter
for snpfilter_dir in */; do
	cd "$snpfilter_dir"  && \
	cd genotype_accuracy  && \
	i=0  && \
	for f in *; do
	    d=Variants_Set_$(printf %04d $((i/1000+1)))  && \
	    mkdir -p "$d"  && \
	    mv "$f" "$d"  && \
	    let i++;
	done
	wait
	cd "$projdir"/snpfilter
done
wait
find . -type f -empty -delete
find . -type d -empty -delete
}
cd $projdir
if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
	if [ "$walkaway" == false ]; then
		echo -e "${magenta}- Do you want to generate visualizations for genotype accuracy and ploidy estimation? ${white}\n"
		read -p "- y(YES) or n(NO) " -t 36000 -n 1 -r
		if [[ ! $REPLY =~ ^[Yy]$ ]]; then
			printf '\n'
			echo -e "${magenta}- skipping visualizations for genotype accuracy and ploidy estimation ${white}\n"
		else
			printf '\n'
			echo -e "${magenta}- generating visualizations for genotype accuracy and ploidy estimation ${white}\n"
			#time main &
		fi
	fi
	if [ "$walkaway" == true ]; then
		if [ "$SNPaccuracy_ReadDepth" == 1 ]; then
			echo -e "${magenta}- generating visualizations for genotype accuracy and ploidy estimation ${white}\n"
			#time main &
		else
			echo -e "${magenta}- skipping visualizations for genotype accuracy and ploidy estimation ${white}\n"
		fi
	fi
fi


######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${yellow}\n- Perform Sequence-based Haplotyping and Filtering\n${blue}##############################################################################${white}\n"
main () {
cd $projdir
cd refgenomes
file=${ref1%.f*}.dict
if test -f "$file"; then
	echo -e "${magenta}- indexed subgenome already exist ${white}\n"
else
	$bwa index -a bwtsw ${ref1}
	$samtools faidx ${ref1}
	$java -jar $picard CreateSequenceDictionary REFERENCE= ${ref1} OUTPUT=${ref1%.f*}.dict
fi

cd $projdir
cd snpfilter
for snpfilter_dir in */; do
	cd $snpfilter_dir
	if [ -d subref ]; then rm -r subref; fi
	if [ -d subsamples ]; then rm -r subsamples; fi
	if [ -d paralog_haplo_filter ]; then rm -r paralog_haplo_filter; fi
	if test -f "$file"; then rm -r haplo.txt; fi
	mkdir subref
	mkdir subsamples
	mkdir paralog_haplo_filter
	cd subref
	export ref1=${ref1}
	awk '{print "$samtools faidx ../../../refgenomes/${ref1} " $2 " >> refpos.fasta"}' ../snplist_nonredun.txt > snpseq_context.sh
	bash snpseq_context.sh
	wait
	rm *.sh
	$bwa index -a bwtsw refpos.fasta
	$samtools faidx refpos.fasta
	$java -jar $picard CreateSequenceDictionary REFERENCE= refpos.fasta OUTPUT= refpos.dict

	cd ../../../preprocess
	for sample in *_precall.bam; do
		outfile=${sample%_${ref1%.fa*}*}
		$samtools bam2fq $sample | awk 'NR%2==0' | awk 'NR%2==1' | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{$1=$1};1' | \
		awk -v sname=${outfile} '{print ">"sname"_seq"NR"-"$1"\n"$2}' | $gzip > ../snpfilter/${snpfilter_dir}/subsamples/${outfile}_precall.fasta.gz
		cd ../snpfilter/${snpfilter_dir}/subref
		$bwa mem -t $loopthreads refpos.fasta ../subsamples/${outfile}_precall.fasta.gz > ../subsamples/${outfile}.sam
		cd ../subsamples/
		$samtools view -b -F 4 ${outfile}.sam | $samtools bam2fq - | awk 'NR%4==1' > ${outfile}_uniqh.fasta
		$samtools view -b -F 4 ${outfile}.sam | $samtools bam2fq - | awk 'NR%2==0' | awk 'NR%2==1' > ${outfile}_uniqs.fasta
		paste -d '\t' ${outfile}_uniqh.fasta ${outfile}_uniqs.fasta > ${outfile}_uniq.fasta
		linelen=$(awk '{print $2}' ${outfile}_uniq.fasta | awk '{print length}' | sort -nr | uniq -c | awk -v sum=$(wc -l ${outfile}_uniq.fasta | awk '{print $1}') '{print $1/sum,$2}' | awk '$1>=0.05' | sort -n -k1,1 | awk 'NR==1{print $2}')
		awk -v len=$linelen 'length($2)>=len' ${outfile}_uniq.fasta | awk -v len=$linelen '{print $1,substr($2,1,len)}' | \
		awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t"); print}' | awk '{print $2"\t"$3}'> ${outfile}.fasta
		for j in $(LC_ALL=C; sort -n -k1,1 ${outfile}.fasta | awk '{print $1}' | uniq); do
			awk -v n="$j" '$0~n{print $0}' ${outfile}.fasta | awk -v n=$j '{for(i=0;i<n;i++) print}' | awk -F' ' '{print $2}' >> ${outfile}_FR.fasta
		done
		awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}'  ${outfile}_FR.fasta | awk '{$1=$1};1' | \
		awk -v sname=${outfile} '{print ">"sname"_seq"NR"-"$1"\n"$2}' | $gzip > ${outfile}.fasta.gz
		rm ${outfile}_precall.fasta.gz ${outfile}_precall.fasta.gz ${outfile}.sam ${outfile}_uniq* ${outfile}.fasta ${outfile}_FR.fasta
		cd ../../../preprocess
	done

	cd ../snpfilter/$snpfilter_dir
	awk '{print "$samtools faidx ../../refgenomes/${ref1} " $2 " >> snpseq_context.fasta"}' snplist_nonredun.txt > snpseq_context.sh
	bash snpseq_context.sh
	awk 1 ORS=',' snpseq_context.fasta | awk '{gsub(/>/,"\n>"); print}' | awk 'NF > 0' | sort -u -t, -k2,2 | tr ',' '\n' | awk 'NF > 0' > temp && mv temp snpseq_context.fasta
	$bwa index -a bwtsw snpseq_context.fasta
	$samtools faidx snpseq_context.fasta
	$java  -jar $picard CreateSequenceDictionary REFERENCE= snpseq_context.fasta OUTPUT=snpseq_context.dict
	for sample in ./subsamples/*.fasta.gz; do
		$bwa mem -t $threads snpseq_context.fasta $sample | $samtools view -F 4 > ${sample%.fasta.gz}_align.txt
		awk '{print $1"\t"$3}' ${sample%.fasta.gz}_align.txt | awk 'BEGIN{OFS="\t"} {gsub(/^.*-/,"",$1); print}' > ${sample%.fasta.gz}_alignh.txt
		$samtools bam2fq ${sample%.fasta.gz}_align.txt | awk 'NR%2==0' | awk 'NR%2==1' > ${sample%.fasta.gz}_aligns.txt
		paste -d '\t' ${sample%.fasta.gz}_alignh.txt ${sample%.fasta.gz}_aligns.txt > ${sample%.fasta.gz}_align.txt
		rm $sample ${sample%.fasta.gz}_alignh.txt ${sample%.fasta.gz}_aligns.txt
	done
	rm snpseq_context.sh snpseq_context*
	rm -rf subref

	for sample in ./subsamples/*_align.txt; do
		outfile=${sample#*/*/} && outfile=${outfile%_align.txt}
		rdvalues=$(awk '{print $2}' $sample | sort -n | uniq)
		for i in $rdvalues; do
			awk -v n=${i} '$0~n{print}'  $sample | sort -s -nrk1,1 > ${sample%_align*}_snp.txt
			rdthreshold=$(awk 'NR == 1 {print $1}' ${sample%_align*}_snp.txt)
			rdthreshold=$(awk -v rdthreshold=$rdthreshold -v ploidy=$ploidy 'BEGIN{printf "%0.6f\n", (rdthreshold*(0.1/ploidy))}')
			awk -v rd=$rdthreshold '{if ($1 >= rd) {print $1"\t"$2"\t"$3}}' ${sample%_align.txt}_snp.txt > ${sample%_align.txt}_haplotypes1.txt
			haplotypes=$(wc -l ${sample%_align.txt}_haplotypes1.txt | awk '{print $1}')
			if [[ "$haplotypes" -le "$ploidy" ]]; then
				haplo_RD=$(awk '{print $1"/"}' ${sample%_align.txt}_haplotypes1.txt | tr -d '\n' )
				haplo=$(awk '{print $3"/"}' ${sample%_align.txt}_haplotypes1.txt | tr -d '\n' )
				printf "$i\t$haplo_RD\t$haplo\t$i" | awk '{gsub("/\t","\t"); print}' | awk '{print $1"\t"$2"\t"$3}' >> ./paralog_haplo_filter/${outfile}_haplotypes.txt
			fi
		done
		rm ${sample%_align.txt}_snp.txt ${sample%_align.txt}_haplotypes1.txt ${sample}
	done

	rmdir subsamples
	cat ./paralog_haplo_filter/*_haplotypes.txt | awk -F '\t' 'BEGIN{OFS="\t"}{print $1}' | sort -u > ./paralog_haplo_filter/SNP_files.txt
	for i in ./paralog_haplo_filter/*_haplotypes.txt; do
		happrecall=${i#*/*/}; happrecall=${happrecall%_haplotypes.txt}
		printf "SNP\t${happrecall}_haplosRD\t${happrecall}_haplos\n" > ${i%_haplotypes.txt}_hapmissingSNP.txt
		awk -F '\t' 'BEGIN{OFS="\t"}{print $1}' $i | cat - ./paralog_haplo_filter/SNP_files.txt | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{print $0"\tNA\tNA"}' | cat $i - | sort -u | sort -k1,1 >> ${i%_haplotypes.txt}_hapmissingSNP.txt
	done
	rm ./paralog_haplo_filter/SNP_files.txt

	touch ./paralog_haplo_filter/consensus_haplos.txt
	if [ -z "$p1" ]; then
		if [ -z "$p2" ]; then
			for i in ./paralog_haplo_filter/*_hapmissingSNP.txt; do
				awk '{print $3}' $i  | awk '{gsub("/","\t"); print}' | paste -d '\t' ./paralog_haplo_filter/consensus_haplos.txt - > consensus_haplos1.txt
				cat consensus_haplos1.txt | while read -r line; do
					sleep $((RANDOM % 2))
          compress=$(echo $line | tr '\t' '\n' | sort | uniq | awk 'NF > 0' | tr '\n' '\t')
					printf '..%s..' "$compress\n" >> consensus_haplos2.txt
				done
				cat consensus_haplos2.txt > ./paralog_haplo_filter/consensus_haplos.txt
				rm consensus_haplos1.txt && rm consensus_haplos2.txt
			done
		else
			for i in ./paralog_haplo_filter/${p1}*_hapmissingSNP.txt ./paralog_haplo_filter/${p2}*_hapmissingSNP.txt; do
				awk '{print $3}' $i  | awk '{gsub("/","\t"); print}' | paste -d '\t' ./paralog_haplo_filter/consensus_haplos.txt - > consensus_haplos1.txt
				cat consensus_haplos1.txt | while read -r line; do
					sleep $((RANDOM % 2))
          compress=$(echo $line | tr '\t' '\n' | sort | uniq | awk 'NF > 0' | tr '\n' '\t')
					printf '..%s..' "$compress\n" >> consensus_haplos2.txt
				done
				cat consensus_haplos2.txt > ./paralog_haplo_filter/consensus_haplos.txt
				rm consensus_haplos1.txt && rm consensus_haplos2.txt
			done
		fi
	fi
	cat ./paralog_haplo_filter/consensus_haplos.txt | while read -r line; do
		sleep $((RANDOM % 2))
    compress=$(echo $line | tr ' ' '\t' | tr '\t\t' '\t' | tr '\t' '\n' | sort | uniq | awk 'NF > 0' | awk '$1 != "NA"' | tr '\n' '\t')
		nHap=$(printf '..%s..' "$compress" | tr '\t' '\n' | awk 'NF > 0' | wc -l)
		printf '..%s..' "$nHap\t$compress\n" >> ./paralog_haplo_filter/consensus_haplo2.txt
	done
	cat ./paralog_haplo_filter/consensus_haplo2.txt > ./paralog_haplo_filter/consensus_haplos.txt
	rm ./paralog_haplo_filter/consensus_haplo2.txt

	#source: http://www.dayofthenewdan.com/2012/12/26/AWK_Linear_Regression.html
	fitline () {
		awk 'BEGIN { FS = "[ ,\t]+" }
		NF == 2 { x_sum += $1
		y_sum += $2
		xy_sum += $1*$2
		x2_sum += $1*$1
		num += 1
		x[NR] = $1
		y[NR] = $2
		}
		END { mean_x = x_sum / num
		mean_y = y_sum / num
		mean_xy = xy_sum / num
		mean_x2 = x2_sum / num
		slope = (mean_xy - (mean_x*mean_y)) / (mean_x2 - (mean_x*mean_x))
		inter = mean_y - slope * mean_x
		for (i = num; i > 0; i--) {
		ss_total += (y[i] - mean_y)**2
		ss_residual += (y[i] - (slope * x[i] + inter))**2
		}
		r2 = 1 - (ss_residual / ss_total)
		printf("Slope      :  %g\n", slope)
		printf("Intercept  :  %g\n", inter)
		printf("R-Squared  :  %g\n", r2)
		}'
	}

	ploidy_level=${snpfilter_dir%x_*}
	if [[ "$ploidy_level" -eq 2 ]]; then
		for i in $( ls -S ./paralog_haplo_filter/*_hapmissingSNP.txt); do
			awk 'NR==1{print $3}' $i > ${i%_hapmissingSNP.txt}_hapcoded.txt
		done
		for j in $(seq 2 1 $(wc -l ./paralog_haplo_filter/consensus_haplos.txt | awk '{print $1}' )); do
			awk -v j=$j 'NR == j' ./paralog_haplo_filter/consensus_haplos.txt | tr '\t' '\n' | awk 'NR>1' | awk 'NF > 0' | awk '{print NR"\t"$1}' > ./paralog_haplo_filter/consensus_haplo_line.txt
			for i in ./paralog_haplo_filter/*_hapmissingSNP.txt; do
				nhaplo=$(awk -v j=$j 'NR==j{print $3}' $i | tr '/' '\n' | awk 'NF > 0' | wc -l)
				if [[ $nhaplo -eq 1 ]]; then
					awk -v j=$j 'NR==j{print $3}' $i | tr '/' '\n' | awk '{print ".\t"$0"\n.\t"$0}' | cat ./paralog_haplo_filter/consensus_haplo_line.txt - | sort -k2,2 | uniq -f1 -d | \
					awk '{gsub(/.\tNA/,"NA\tNA"); print}' | awk '{print $1" "$1}' | tr ' ' '/' | awk '{gsub("NA/NA","NA"); print}' >> ${i%_hapmissingSNP.txt}_hapcoded.txt
				fi
				if [[ $nhaplo -eq 2 ]]; then
					awk -v j=$j 'NR==j{print $3}' $i | tr '/' '\n' | awk '{print ".\t"$0"\n.\t"$0}' | cat ./paralog_haplo_filter/consensus_haplo_line.txt -  | sort -k2,2 | uniq -f1 -d | \
					awk '{gsub(/.\tNA/,"NA\tNA"); print}' | awk '{print $1}' | tr '\n' ' ' | awk '{$1=$1};1' | tr ' ' '/' | awk '{gsub("NA/NA","NA"); print}' >> ${i%_hapmissingSNP.txt}_hapcoded.txt
				fi
			done; wait
		done
	fi
	if [[ "$ploidy_level" -eq 4 ]]; then
		for i in $( ls -S ./paralog_haplo_filter/*_hapmissingSNP.txt); do
			awk 'NR==1{print $3}' $i > ${i%_hapmissingSNP.txt}_hapcoded.txt
		done
		for j in $(seq 2 1 $(wc -l ./paralog_haplo_filter/consensus_haplos.txt | awk '{print $1}' )); do
			awk -v j=$j 'NR == j' ./paralog_haplo_filter/consensus_haplos.txt | tr '\t' '\n' | awk 'NR>1' | awk 'NF > 0' | awk '{print NR"\t"$1}' > ./paralog_haplo_filter/consensus_haplo_line.txt
			for i in ./paralog_haplo_filter/*_hapmissingSNP.txt; do
				nhaplo=$(awk -v j=$j 'NR==j{print $3}' $i | tr '/' '\n' | awk 'NF > 0' | wc -l)
				exp=$(awk -v j=$j 'NR==j{print $3}' $i | tr '/' ' ' | awk 'NF > 0')
				if [[ $nhaplo -eq 1 ]]; then
					awk -v j=$j 'NR==j{print $3}' $i | tr '/' '\n' | awk '{print ".\t"$0"\n.\t"$0}' | cat ./paralog_haplo_filter/consensus_haplo_line.txt - | sort -k2,2 | uniq -f1 -d | \
					awk '{gsub(/.\tNA/,"NA\tNA"); print}' | awk '{print $1" "$1" "$1" "$1}' | tr ' ' '/' | awk '{gsub("NA/NA","NA"); print}' >> ${i%_hapmissingSNP.txt}_hapcoded.txt
				fi
if [[ $nhaplo -eq 2 ]]; then
obs1=$(printf "0 1 3")
obs2=$(printf "0 2 2")
exp1=$(( paste <(printf '..%s..' "$obs1\n" | tr ' ' '\n') <(printf '..%s..' "$exp\n" | tr ' ' '\n') | fitline | awk 'NR==3{print $3}' ))
exp2=$(( paste <(printf '..%s..' "$obs2\n" | tr ' ' '\n') <(printf '..%s..' "$exp\n" | tr ' ' '\n') | fitline | awk 'NR==3{print $3}' ))
	if ("$exp1" > "$exp2"); then echo yes; fi
		awk -v j=$j 'NR==j{print $3}' $i | tr '/' '\n' | awk '{print ".\t"$0"\n.\t"$0}' | cat ./paralog_haplo_filter/consensus_haplo_line.txt - | sort -k2,2 | uniq -f1 -d | \
		awk '{gsub(/.\tNA/,"NA\tNA"); print}' | awk '{print $1" "$1}' | tr ' ' '/' | awk '{gsub("NA/NA","NA"); print}' >> ${i%_hapmissingSNP.txt}_hapcoded.txt
fi
if [[ $nhaplo -eq 3 ]]; then
	awk -v j=$j 'NR==j{print $3}' $i | tr '/' '\n' | awk '{print ".\t"$0"\n.\t"$0}' | cat ./paralog_haplo_filter/consensus_haplo_line.txt - | sort -k2,2 | uniq -f1 -d | \
	awk '{gsub(/.\tNA/,"NA\tNA"); print}' | awk '{print $1" "$1" "$1}' | tr ' ' '/' | awk '{gsub("NA/NA","NA"); print}' >> ${i%_hapmissingSNP.txt}_hapcoded.txt
fi
				if [[ $nhaplo -eq 4 ]]; then
					awk -v j=$j 'NR==j{print $3}' $i | tr '/' '\n' | awk '{print ".\t"$0"\n.\t"$0}' | cat ./paralog_haplo_filter/consensus_haplo_line.txt -  | sort -k2,2 | uniq -f1 -d | \
					awk '{gsub(/.\tNA/,"NA\tNA"); print}' | awk '{print $1}' | tr '\n' ' ' | awk '{$1=$1};1' | tr ' ' '/' | awk '{gsub("NA/NA","NA"); print}' >> ${i%_hapmissingSNP.txt}_hapcoded.txt
				fi
			done; wait
		done
	fi
	if [[ "$ploidy_level" -eq 6 ]]; then
		for i in $( ls -S ./paralog_haplo_filter/*_hapmissingSNP.txt); do
			awk 'NR==1{print $3}' $i > ${i%_hapmissingSNP.txt}_hapcoded.txt
		done
		for j in $(seq 2 1 $(wc -l ./paralog_haplo_filter/consensus_haplos.txt | awk '{print $1}' )); do
			awk -v j=$j 'NR == j' ./paralog_haplo_filter/consensus_haplos.txt | tr '\t' '\n' | awk 'NR>1' | awk 'NF > 0' | awk '{print NR"\t"$1}' > ./paralog_haplo_filter/consensus_haplo_line.txt
			for i in ./paralog_haplo_filter/*_hapmissingSNP.txt; do
				nhaplo=$(awk -v j=$j 'NR==j{print $3}' $i | tr '/' '\n' | awk 'NF > 0' | wc -l)
				if [[ $nhaplo -eq 1 ]]; then
					awk -v j=$j 'NR==j{print $3}' $i | tr '/' '\n' | awk '{print ".\t"$0"\n.\t"$0}' | cat ./paralog_haplo_filter/consensus_haplo_line.txt - | sort -k2,2 | uniq -f1 -d | \
					awk '{gsub(/.\tNA/,"NA\tNA"); print}' | awk '{print $1" "$1" "$1" "$1" "$1" "$1}' | tr ' ' '/' | awk '{gsub("NA/NA","NA"); print}' >> ${i%_hapmissingSNP.txt}_hapcoded.txt
				fi
				if [[ $nhaplo -eq 6 ]]; then
					awk -v j=$j 'NR==j{print $3}' $i | tr '/' '\n' | awk '{print ".\t"$0"\n.\t"$0}' | cat ./paralog_haplo_filter/consensus_haplo_line.txt -  | sort -k2,2 | uniq -f1 -d | \
					awk '{gsub(/.\tNA/,"NA\tNA"); print}' | awk '{print $1}' | tr '\n' ' ' | awk '{$1=$1};1' | tr ' ' '/' | awk '{gsub("NA/NA","NA"); print}' >> ${i%_hapmissingSNP.txt}_hapcoded.txt
				fi
			done; wait
		done
	fi
	if [[ "$ploidy_level" -eq 8 ]]; then
		for i in $( ls -S ./paralog_haplo_filter/*_hapmissingSNP.txt); do
			awk 'NR==1{print $3}' $i > ${i%_hapmissingSNP.txt}_hapcoded.txt
		done
		for j in $(seq 2 1 $(wc -l ./paralog_haplo_filter/consensus_haplos.txt | awk '{print $1}' )); do
			awk -v j=$j 'NR == j' ./paralog_haplo_filter/consensus_haplos.txt | tr '\t' '\n' | awk 'NR>1' | awk 'NF > 0' | awk '{print NR"\t"$1}' > ./paralog_haplo_filter/consensus_haplo_line.txt
			for i in ./paralog_haplo_filter/*_hapmissingSNP.txt; do
				nhaplo=$(awk -v j=$j 'NR==j{print $3}' $i | tr '/' '\n' | awk 'NF > 0' | wc -l)
				if [[ $nhaplo -eq 1 ]]; then
					awk -v j=$j 'NR==j{print $3}' $i | tr '/' '\n' | awk '{print ".\t"$0"\n.\t"$0}' | cat ./paralog_haplo_filter/consensus_haplo_line.txt - | sort -k2,2 | uniq -f1 -d | \
					awk '{gsub(/.\tNA/,"NA\tNA"); print}' | awk '{print $1" "$1" "$1" "$1" "$1" "$1" "$1" "$1}' | tr ' ' '/' | awk '{gsub("NA/NA","NA"); print}' >> ${i%_hapmissingSNP.txt}_hapcoded.txt
				fi
				if [[ $nhaplo -eq 8 ]]; then
					awk -v j=$j 'NR==j{print $3}' $i | tr '/' '\n' | awk '{print ".\t"$0"\n.\t"$0}' | cat ./paralog_haplo_filter/consensus_haplo_line.txt -  | sort -k2,2 | uniq -f1 -d | \
					awk '{gsub(/.\tNA/,"NA\tNA"); print}' | awk '{print $1}' | tr '\n' ' ' | awk '{$1=$1};1' | tr ' ' '/' | awk '{gsub("NA/NA","NA"); print}' >> ${i%_hapmissingSNP.txt}_hapcoded.txt
				fi
			done; wait
		done
	fi

	rm  ./paralog_haplo_filter/consensus_haplo_line.txt
	awk 'NR>1{$1=""; print $0}' ./paralog_haplo_filter/consensus_haplos.txt | awk '{$1=$1};1' | awk '{gsub(/ /,""); print}' | awk 'BEGIN{print "Haplotype_sequences"}1' > ./paralog_haplo_filter/Haplotype_Seq_Context.txt
	awk '{print $1}' $(ls ./paralog_haplo_filter/*_hapmissingSNP.txt | head -n 1 ) | paste -d '\t' - ./paralog_haplo_filter/*_hapcoded.txt ./paralog_haplo_filter/Haplotype_Seq_Context.txt | \
	awk -v nsample=$(ls ./paralog_haplo_filter/*_hapcoded.txt | wc -l) '{print (gsub(/NA/,"NA"))/nsample"\t"$0}' | awk -v n=$genotype_missingness '$1<=n {print}' | awk '{$1=""}1' | awk 'BEGIN{OFS="\t"}{$1=$1};1' > ./paralog_haplo_filter/${pop}_rd${ploidy}_population_haplotype.txt

	cat ./paralog_haplo_filter/${pop}_rd${ploidy}_population_haplotype.txt | while read -r line; do
	sleep $((RANDOM % 2))
  hap=$( echo $line | awk '{print $1}' )
	awk -v hap=$hap '$0~hap {print $1}' snplist.txt >> ./paralog_haplo_filter/haplo_filtered_snp.txt
	done

	minRD=${snpfilter_dir%x*}
	for i in $(ls ./paralog_haplo_filter/${pop}_${ploidy}x_rd${minRD}_maf${maf}_* | grep -v population_haplotype ) ; do
	awk 'BEGIN{OFS="\t"}NR==FNR {h[$1] = $0; next} {print $0,h[$1]}' haplo_filtered_snp.txt $i | awk  -F '\t' 'BEGIN{OFS="\t"}$NF!=""' | awk 'BEGIN{OFS="\t"}NF{NF-=1};1' > ./paralog_haplo_filter/${i%.txt}_hapfiltered.txt
	done
done
}
cd $projdir
if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
	if [ "$walkaway" == false ]; then
		echo -e "${magenta}- Do you want to perform sequence-based haplotyping and filtering? ${white}\n"
		read -p "- y(YES) or n(NO) " -t 36000 -n 1 -r
		if [[ ! $REPLY =~ ^[Yy]$ ]]; then
			printf '\n'
			echo -e "${magenta}- skipping sequence-based haplotyping and filtering ${white}\n"
		else
			printf '\n'
			echo -e "${magenta}- performing sequence-based haplotyping and filtering ${white}\n"
			#time main &>> log.out
			echo -e "${magenta}- Under Development ${white}\n"
		fi
	fi
	if [ "$walkaway" == true ]; then
		if [ "$paralogfiltering_haplotyping" == 1 ]; then
			echo -e "${magenta}- performing sequence-based haplotyping and filtering ${white}\n"
			#time main &>> log.out
			echo -e "${magenta}- Under Development ${white}\n"
		else
			echo -e "${magenta}- skipping sequence-based haplotyping and filtering ${white}\n"
		fi
	fi
fi


######################################################################################################################################################
cd ${projdir}
if [[ "$samples_list" == "samples_list_node_1.txt" ]] && [[ -d "snpfilter" ]]; then
  find ../ -size 0 -delete >/dev/null 2>&1 &&
  rm call*  2> /dev/null &&
  ls ./snpfilter/*/*_plusSD.txt 2> /dev/null | xargs rm 2> /dev/null &&
  ls ./snpfilter/*/*SD_1_G*G*.txt 2> /dev/null | xargs rm 2> /dev/null &&
  ls *node*.txt | grep -v 'samples_list' | xargs rm 2> /dev/null &&
	rm steps.txt 2> /dev/null &&
	mv ${projdir}/GBSapp_run_node_1.sh ${projdir}/GBSapp_run_node_1_done.sh 2> /dev/null &&
	mv ${projdir}/GBSapp_run_node.sh ${projdir}/GBSapp_run_node_done.sh 2> /dev/null &&
  wait
  if [[ "$biallelic" == true ]]; then mv snpfilter snpfilter_biallelic; fi
  touch ../Analysis_Complete
  wait
else
	touch ../Analysis_Complete_${samples_list}
fi
wait
echo -e "${magenta}- Run Complete. ${white}\n"
