
if [ -z "$slurm_module" ]; then
 export slurm_module=true
fi
if [ -z "$threads" ]; then
	export threads=$(nproc --all)
	if [[ "$threads" -ge 4 ]]; then
		export threads=$((threads-2))
	fi
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
if [ -z "$uniquely_mapped" ]; then
	export uniquely_mapped=true
fi
if [ -z "$minmapq" ]; then
	export minmapq=20
fi
if [ -z "$maxHaplotype" ]; then
	export maxHaplotype=128
fi
if [ -z "$haplome_number" ]; then
	export haplome_number=1
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
mkdir -p "${projdir}"/tmp
export TMPDIR="${projdir}"/tmp


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
	  export N=1 && export loopthreads=threads
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
for i in *.gz; do
	gunzip $i >/dev/null 2>&1
done

if ls ./*.ngm 1> /dev/null 2>&1; then
	:
else
	checknfastafiles=$(ls *.f* | grep -v .fai | grep -v .ngm | grep -v _original.fasta | wc -l)
	if [[ $checknfastafiles -gt 3 ]]; then
		echo -e "${magenta}- expecting only 3 fasta files for reference genome ${white}\n"
		echo -e "${magenta}- GBSapp will quit in 5 seconds ${white}\n"
		sleep 5 && exit 1
	fi

	reffilenames=($ref1 $ref2 $ref3)
	for refg in "${reffilenames[@]}"; do
		export ncontigscaffold=$(grep '>' $refg | wc -l)
		if [[ $ncontigscaffold -gt 1000 ]]; then
			nfakechr=$((threads/2))
			cat $refg | awk '/^>/{if(N>0) printf("\n"); ++N; printf("%s\t",$0);next;} {printf("%s",$0);}END{printf("\n");}' > panref0.txt
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
				awk -v chrid=${i} -v refgenome=${refg%.f*} '{print refgenome"_"chrid"\t"$1"\t"$2}' >> contigscaffold_index.txt
			done
			cp $refg ${refg%.f*}_original.fasta
			cat /dev/null > $refg
			for filename in Chr*.txt; do
				echo ">""${filename%.txt}" >> $refg
				cat "$filename" >> $refg
			done
			rm Chr*
		fi
	done

	export ncontigscaffold=$(grep '>' $ref1 | wc -l)
	if [[ $ncontigscaffold -gt 1000 ]]; then
		mkdir split
		for reffile in *.f*; do
			awk -v reffile=${reffile%.f*} '/^>/{close("file.fasta"f);f++}{print $0 > "./split/file_"reffile"_"f}' "$reffile" & PID=$!
		  wait $PID
		done
		cd split
		export checksplit=$( wc -c file_* | head -n -1 | awk '($1 > 500000000 )' | wc -l )
		if [[ "$checksplit" -gt 0 ]]; then
			for i in file*; do
				name=$(awk 'NR==1{print $1}' "$i" | awk '{gsub(/>/,""); print}')
				awk 'NR>1{print $0}' $i > ${name}.fasta
				count=$(wc -c ${name}.fasta | awk '{print $1}')
				tr -d '\n' < "${name}".fasta | fold -w 100 |\
				split -d -l 5000000
				for outfile in x*; do
					awk -v name="$name" -v outfile="$outfile" -i inplace 'BEGINFILE{print ">"name"_"outfile}{print}' $outfile
					subgenome=${i%.fasta_*}; subgenome=${subgenome##file_}
					mv $outfile ${subgenome}_${name}_${outfile}.txt
				done
			done
			cd ../
			for i in *.f*; do
				mv $i ./old_"${i%.f*}_fasta.txt"
				cat ./split/${i%.f*}_Chr* > $i
			done
			wait
			rm -r split
		else
			cd ${projdir}/refgenomes
			rm -r split
		fi
	fi
fi


cd $projdir
cd refgenomes
if ls ./*.ngm 1> /dev/null 2>&1; then
	echo -e "${magenta}- indexed genome available ${white}\n"
else

	echo -e "${magenta}- indexing reference subgenome-1 ${white}\n"
	awk '{ sub("\r$",""); print}' $ref1 | awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' > ref.txt
	n=">${ref1%.f*}_"
	awk '{ sub("\r$",""); print}' ref.txt | awk -v n="$n" '{gsub(n,">"); print}' | awk -v n="$n" '{gsub(/>/,n); print}' > $ref1
	rm ref.txt
	$samtools faidx $ref1
	$java -jar $picard CreateSequenceDictionary REFERENCE=$ref1 OUTPUT=${ref1%.f*}.dict
	$ngm -r $ref1

  echo -e "${magenta}- indexing reference subgenome-2 ${white}\n"
  awk '{ sub("\r$",""); print}' $ref2 | awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' > ref.txt
  n=">${ref2%.f*}_"
  awk '{ sub("\r$",""); print}' ref.txt | awk -v n="$n" '{gsub(n,">"); print}' | awk -v n="$n" '{gsub(/>/,n); print}' > $ref2
  rm ref.txt
  $samtools faidx $ref2
  $java -jar $picard CreateSequenceDictionary REFERENCE=$ref2 OUTPUT=${ref2%.f*}.dict
  $ngm -r $ref2

  echo -e "${magenta}- indexing reference subgenome-3 ${white}\n"
  awk '{ sub("\r$",""); print}' $ref3 | awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' > ref.txt
  n=">${ref3%.f*}_"
  awk '{ sub("\r$",""); print}' ref.txt | awk -v n="$n" '{gsub(n,">"); print}' | awk -v n="$n" '{gsub(/>/,n); print}' > $ref3
  rm ref.txt
  $samtools faidx $ref3
  $java -jar $picard CreateSequenceDictionary REFERENCE=$ref3 OUTPUT=${ref3%.f*}.dict
  $ngm -r $ref3

fi

declare -a arr=("${ref1%.f*}.dict" "${ref1}" "${ref1}.fai" "${ref2%.f*}.dict" "${ref2}" "${ref2}.fai"  "${ref3%.f*}.dict" "${ref3}" "${ref3}.fai")
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
  if [[ "$(ls -A *.f* 2> /dev/null | wc -l)" -eq 0 ]] || [[ "$(ls -A ./preprocess/alignment/*sam* | wc -l)" -gt 0 ]]; then
    :> filename_reformatted.txt
    :> flushed_reads.txt
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
			fi
		fi

		cd ${projdir}/samples/pe
		if [ "$(ls -A ../pe 2> /dev/null)" ]; then
      if [ -z "$(ls -A ../se 2> /dev/null)" ]; then
				echo -e "${magenta}- only paired-end reads available in pe-folder ${white}\n"
				for i in *R1.f*; do (
					sleep $((RANDOM % 2))
          mv ${i%R1.f*}R2.f* ../ 2> /dev/null && rm ${i%R1.f*}R2.f* 2> /dev/null
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
			fi
		fi

		cd ${projdir}/samples/pe
		if [ "$(ls -A ../se 2> /dev/null)" ]; then
			if [ "$(ls -A ../pe 2> /dev/null)" ]; then
				for i in *R1.f*; do (
					sleep $((RANDOM % 2))
          mv ${i%R1.f*}R2.f* ../
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
  						awk 'NR%2==0' <(zcat $i) | awk 'NR%2==1' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=1000; i++){x=int(rand()*NR) + 1; print a[x];}}' > ${i}_length_distribution.txt
  					fi
  					if [[ "${fa_fq}" == ">" ]]; then
  						awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' <(zcat $i) | \
  						awk 'NR%2==0' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=1000; i++){x=int(rand()*NR) + 1; print a[x];}}' > ${i}_length_distribution.txt
  					fi
  				else
  					if [[ "${fa_fq}" == "@" ]]; then
  						awk 'NR%2==0' $i | awk 'NR%2==1' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=1000; i++){x=int(rand()*NR) + 1; print a[x];}}' > ${i}_length_distribution.txt
  					fi
  					if [[ "${fa_fq}" == ">" ]]; then
  						awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' $i | \
  						awk 'NR%2==0' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=1000; i++){x=int(rand()*NR) + 1; print a[x];}}' > ${i}_length_distribution.txt
  					fi
  				fi
        fi ) &
        if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
          wait
        fi
			done

      :> length_distribution.txt
      for lenfile in *_length_distribution.txt; do
        cat $lenfile >> length_distribution.txt
      done
  		wait
      rm *_length_distribution.txt
      awk '{print length($0)}' length_distribution.txt | awk '$1>=64' | sort -n > tmp.txt && :> length_distribution.txt &&
      mv tmp.txt length_distribution.txt &&
			export max_seqread_len=$(awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}' length_distribution.txt)
			rm length_distribution.txt

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
                cat "${i%.f*}"_RE1.fasta.gz "${i%.f*}"_RE2.fasta.gz "${i%.f*}"_RE1RE2.fasta.gz > "${i%.f*}".fasta.gz &&
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
                cat ${i%.f*}_RE1.fasta.gz ${i%.f*}_RE2.fasta.gz ${i%.f*}_RE1RE2.fasta.gz > ${i%.f*}.fasta.gz &&
                rm ${i%.f*}_RE1.fasta.gz ${i%.f*}_RE2.fasta.gz ${i%.f*}_RE1RE2.fasta.gz &&
                wait
              fi
              if [[ "${fa_fq}" == ">" ]]; then
                awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $i | awk 'NR%2==0' | awk '{gsub(/ATGCAT/,"ATGCAT\nATGCAT");}1' | \
                awk 'length >= 64 && length <= 600' | grep '^ATGCAT.*ATGCAT$\|^CATG.*CATG$' | awk '{print ">frag"NR"\n"$0}' | $gzip > ${i%.f*}.tmp1.gz &&
                awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $i | awk 'NR%2==0' | awk '{gsub(/CATG/,"CATG\nCATG");}1' | \
                awk 'length >= 64 && length <= 600' | grep '^CATG.*CATG$' | awk '{print ">frag"NR"\n"$0}' | $gzip > ${i%.f*}.tmp2.gz &&
                awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $i | awk 'NR%2==0' | awk '{gsub(/ATGCAT/,"ATGCAT\nATGCAT");gsub(/CATG/,"CATG\nCATG");}1' | \
                awk 'length >= 64 && length <= 600' | grep '^ATGCAT.*CATG$\|^CATG.*ATGCAT$' | awk '{print ">frag"NR"\n"$0}' | $gzip > ${i%.f*}.tmp3.gz &&
                cat "${i%.f*}".tmp1.gz "${i%.f*}".tmp2.gz "${i%.f*}".tmp3.gz > "${i%.f*}".tmp.gz &&
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
			time main &>> ${projdir}/log.out
		fi
	else
		echo -e "${magenta}- samples' PE fastq filenames requires formatting (i.e. needs to end in "_R1.fastq" or ".R1.fastq" and "_R2.fastq" or ".R2.fastq") ${white}\n"
		sleep 5 && exit 1
	fi
else
	if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
		time main &>> ${projdir}/log.out
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

	if [[ "$samples_list" == "samples_list_node_1.txt" ]] && [[ "$(ls -A ./samples/*.f* 2> /dev/null | wc -l)" -gt 0 ]] && [[ "$(ls -A ./preprocess/alignment/*sam* | wc -l)" -eq 0 ]]; then
		for i in samples_list_node_*.txt; do
			:> ${i%.txt}_hold.txt
			while read line; do
				sleep $((RANDOM % 2))
        ls -l ./samples/$line | awk '{print $5"\t"$9}' >> ${i%.txt}_hold.txt
			done < $i
			sort -nr -k1 ${i%.txt}_hold.txt | awk '{gsub(/.\/samples\//,""); print $2}' > $i
			rm "${i%.txt}"_hold.txt
		done
	fi
}
if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
	time main &>> ${projdir}/log.out
fi


######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${yellow}\n- GBSapp is Performing Read Alignments & Alignment Post-Processing\n${blue}##############################################################################${white}\n"

cd $projdir
if [[ $nodes -gt 1 ]] && [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
	rm -rf /tmp/"${samples_list%.txt}" 2> /dev/null
	mkdir -p /tmp/${samples_list%.txt}/refgenomes /tmp/${samples_list%.txt}/samples
  mkdir -p /tmp/${samples_list%.txt}/preprocess/tmp /tmp/${samples_list%.txt}/preprocess/alignment /tmp/${samples_list%.txt}/snpcall/tmp
	cp -r ${projdir}/refgenomes/* /tmp/${samples_list%.txt}/refgenomes/
	while IFS="" read -r i || [ -n "$i" ]; do
		cp ${i%.f*}_uniq_R1.fasta.gz /tmp/${samples_list%.txt}/samples/ 2> /dev/null
		cp ${projdir}/preprocess/alignment/${i%.f*}_redun.sam.gz /tmp/${samples_list%.txt}/preprocess/alignment/ 2> /dev/null
		cp ${projdir}/preprocess/${i%.f*}_*_precall.bam* /tmp/${samples_list%.txt}/preprocess/ 2> /dev/null
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
          cat ${i%.f*}_R1_uniq.txt.gz ${i%.f*}_R2_uniq.txt.gz | zcat | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{$1=$1};1' | \
          awk '{print ">seq"NR"_pe-"$1"\n"$2}' | gzip > ${i%.f*}_uniq.fasta.gz&&
          rm "${i%.f*}"_R1_uniq.txt.gz "${i%.f*}"_R2_uniq.txt.gz &&
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
    if [[ $nodes -eq 1 ]]; then cd ${projdir}/samples/ ; fi
    if [[ $nodes -gt 1 ]] && test -f ${projdir}/GBSapp_run_node_1.sh; then cd /tmp/${samples_list%.txt}/samples/ ; fi
    if test ! -f ${projdir}/precall_done.txt && test ! -f ${projdir}/alignment_done; then
      while IFS="" read -r alignfq || [ -n "$alignfq" ]; do
        sleep $((RANDOM % 2))
        if test ! -f ../preprocess/alignment/${alignfq%.f*}_redun.sam.gz; then
            $ngm -r ../refgenomes/$ref1 --qry ${alignfq%.f*}_uniq.fasta.gz -o ../preprocess/alignment/${alignfq%.f*}_redun.sam -t $ngmthreads --topn 6 --strata 6 --affine &&
            awk '/@HD/ || /@SQ/{print}' ../preprocess/alignment/${alignfq%.f*}_redun.sam 2> /dev/null > ../preprocess/alignment/${alignfq%.f*}_redun_head.sam
            grep -v '^@' ../preprocess/alignment/${alignfq%.f*}_redun.sam 2> /dev/null | awk -F"\t" 'BEGIN{FS=OFS="\t"} {$11="*"; print $0}' | cat ../preprocess/alignment/${alignfq%.f*}_redun_head.sam - | gzip  > ../preprocess/alignment/${alignfq%.f*}_redun_${ref1%.f*}.sam.gz &&
            rm ../preprocess/alignment/"${alignfq%.f*}"_redun.sam &&
            cp -rn ../preprocess/alignment/${alignfq%.f*}_redun_${ref1%.f*}.sam.gz ${projdir}/preprocess/alignment/ &&
            rm "${alignfq%.f*}"_uniq.fastq.gz ../preprocess/alignment/"${alignfq%.f*}"_redun_head.sam &&
            wait
            $ngm -r ../refgenomes/$ref2 --qry ${alignfq%.f*}_uniq.fasta.gz -o ../preprocess/alignment/${alignfq%.f*}_redun.sam -t $ngmthreads --topn 6 --strata 6 --affine &&
            awk '/@HD/ || /@SQ/{print}' ../preprocess/alignment/${alignfq%.f*}_redun.sam 2> /dev/null > ../preprocess/alignment/${alignfq%.f*}_redun_head.sam
            grep -v '^@' ../preprocess/alignment/${alignfq%.f*}_redun.sam 2> /dev/null | awk -F"\t" 'BEGIN{FS=OFS="\t"} {$11="*"; print $0}' | cat ../preprocess/alignment/${alignfq%.f*}_redun_head.sam - | gzip  > ../preprocess/alignment/${alignfq%.f*}_redun_${ref2%.f*}.sam.gz &&
            rm ../preprocess/alignment/"${alignfq%.f*}"_redun.sam &&
            cp -rn ../preprocess/alignment/${alignfq%.f*}_redun_${ref2%.f*}.sam.gz ${projdir}/preprocess/alignment/ &&
            rm "${alignfq%.f*}"_uniq.fastq.gz ../preprocess/alignment/"${alignfq%.f*}"_redun_head.sam &&
            wait
            $ngm -r ../refgenomes/$ref3 --qry ${alignfq%.f*}_uniq.fasta.gz -o ../preprocess/alignment/${alignfq%.f*}_redun.sam -t $ngmthreads --topn 6 --strata 6 --affine &&
            awk '/@HD/ || /@SQ/{print}' ../preprocess/alignment/${alignfq%.f*}_redun.sam 2> /dev/null > ../preprocess/alignment/${alignfq%.f*}_redun_head.sam
            grep -v '^@' ../preprocess/alignment/${alignfq%.f*}_redun.sam 2> /dev/null | awk -F"\t" 'BEGIN{FS=OFS="\t"} {$11="*"; print $0}' | cat ../preprocess/alignment/${alignfq%.f*}_redun_head.sam - | gzip  > ../preprocess/alignment/${alignfq%.f*}_redun_${ref3%.f*}.sam.gz &&
            rm ../preprocess/alignment/"${alignfq%.f*}"_redun.sam &&
            cp -rn ../preprocess/alignment/${alignfq%.f*}_redun_${ref3%.f*}.sam.gz ${projdir}/preprocess/alignment/ &&
            rm "${alignfq%.f*}"_uniq.fastq.gz ../preprocess/alignment/"${alignfq%.f*}"_redun_head.sam &&
            wait
        fi
      done < <( cat ${projdir}/${samples_list} )
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

    while IFS="" read -r i || [ -n "$i" ]; do (
      printf '\n###---'${i%.f*}'---###\n' > ${projdir}/alignment_summaries/${i%.f*}_summ_${ref1%.f*}.txt &&
      zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz | grep -v '^@PG' | tr ' ' '\t' | $samtools flagstat - >> ${projdir}/alignment_summaries/${i%.f*}_summ_${ref1%.f*}.txt &&
      printf '########################################################################################################\n\n' >> ${projdir}/alignment_summaries/${i%.f*}_summ_${ref1%.f*}.txt &&
      printf 'copy_number\tFrequency\tPercentage\n' > ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_Read_histogram_${ref1%.f*}.txt &&
      $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) | awk '{print $1}' | awk '{gsub(/_pe-/,"\t");gsub(/seq/,"");}1' | \
      awk '{while ($2-- > 0) print $1}' | awk '{!seen[$0]++}END{for (i in seen) print seen[i]}' | awk '{!seen[$0]++}END{for (i in seen) print i, seen[i]}' > ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_${ref1%.f*}.txt  &&
      awk 'NR==FNR{sum+= $2; next;} {printf("%s\t%s\t%3.3f%%\t%3.0f\n",$1,$2,100*$2/sum,100*$2/sum)}' ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_${ref1%.f*}.txt ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_${ref1%.f*}.txt | awk '$4 > 0' > ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot_${ref1%.f*}.txt &&
      unset IFS; printf "%s\t%s\t%s\t%*s\n" $(sed 's/$/ |/' ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot_${ref1%.f*}.txt) | tr ' ' '|' | sort -T ./tmp/ -k1,1 -n >> ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_Read_histogram_${ref1%.f*}.txt &&
      rm ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_${ref1%.f*}.txt ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot_${ref1%.f*}.txt
      printf '\n###---'${i%.f*}'---###\n' > ${projdir}/alignment_summaries/${i%.f*}_summ_${ref2%.f*}.txt &&
      zcat ./alignment/${i%.f*}_redun_${ref2%.f*}.sam.gz | grep -v '^@PG' | tr ' ' '\t' | $samtools flagstat - >> ${projdir}/alignment_summaries/${i%.f*}_summ_${ref2%.f*}.txt &&
      printf '########################################################################################################\n\n' >> ${projdir}/alignment_summaries/${i%.f*}_summ_${ref2%.f*}.txt &&
      printf 'copy_number\tFrequency\tPercentage\n' > ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_Read_histogram_${ref2%.f*}.txt &&
      $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun_${ref2%.f*}.sam.gz 2> /dev/null) | awk '{print $1}' | awk '{gsub(/_pe-/,"\t");gsub(/seq/,"");}1' | \
      awk '{while ($2-- > 0) print $1}' | awk '{!seen[$0]++}END{for (i in seen) print seen[i]}' | awk '{!seen[$0]++}END{for (i in seen) print i, seen[i]}' > ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_${ref2%.f*}.txt  &&
      awk 'NR==FNR{sum+= $2; next;} {printf("%s\t%s\t%3.3f%%\t%3.0f\n",$1,$2,100*$2/sum,100*$2/sum)}' ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_${ref2%.f*}.txt ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_${ref2%.f*}.txt | awk '$4 > 0' > ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot_${ref2%.f*}.txt &&
      unset IFS; printf "%s\t%s\t%s\t%*s\n" $(sed 's/$/ |/' ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot_${ref2%.f*}.txt) | tr ' ' '|' | sort -T ./tmp/ -k1,1 -n >> ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_Read_histogram_${ref2%.f*}.txt &&
      rm ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_${ref2%.f*}.txt ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot_${ref2%.f*}.txt
      printf '\n###---'${i%.f*}'---###\n' > ${projdir}/alignment_summaries/${i%.f*}_summ_${ref2%.f*}.txt &&
      zcat ./alignment/${i%.f*}_redun_${ref3%.f*}.sam.gz | grep -v '^@PG' | tr ' ' '\t' | $samtools flagstat - >> ${projdir}/alignment_summaries/${i%.f*}_summ_${ref3%.f*}.txt &&
      printf '########################################################################################################\n\n' >> ${projdir}/alignment_summaries/${i%.f*}_summ_${ref3%.f*}.txt &&
      printf 'copy_number\tFrequency\tPercentage\n' > ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_Read_histogram_${ref3%.f*}.txt &&
      $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun_${ref3%.f*}.sam.gz 2> /dev/null) | awk '{print $1}' | awk '{gsub(/_pe-/,"\t");gsub(/seq/,"");}1' | \
      awk '{while ($2-- > 0) print $1}' | awk '{!seen[$0]++}END{for (i in seen) print seen[i]}' | awk '{!seen[$0]++}END{for (i in seen) print i, seen[i]}'  > ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_${ref3%.f*}.txt  &&
      awk 'NR==FNR{sum+= $2; next;} {printf("%s\t%s\t%3.3f%%\t%3.0f\n",$1,$2,100*$2/sum,100*$2/sum)}' ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_${ref3%.f*}.txt ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_${ref3%.f*}.txt | awk '$4 > 0' > ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot_${ref3%.f*}.txt &&
      unset IFS; printf "%s\t%s\t%s\t%*s\n" $(sed 's/$/ |/' ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot_${ref3%.f*}.txt) | tr ' ' '|' | sort -T ./tmp/ -k1,1 -n >> ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_Read_histogram_${ref3%.f*}.txt &&
      rm ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_${ref3%.f*}.txt ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot_${ref3%.f*}.txt


      if test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai && test ! -f ${projdir}/preprocess/${i%.f*}_${ref2%.f*}_precall.bam.bai && test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_${ref2%.f*}_precall.bam.bai; then
        $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) | grep -v '^@' > ./alignment/${i%.f*}_redun.sam
        $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun_${ref2%.f*}.sam.gz 2> /dev/null) | grep -v '^@' >> ./alignment/${i%.f*}_redun.sam
        $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun_${ref3%.f*}.sam.gz 2> /dev/null) | grep -v '^@' >> ./alignment/${i%.f*}_redun.sam
        awk '$3 != "*"' ./alignment/${i%.f*}_redun.sam 2> /dev/null | awk '$6 != "*"' 2> /dev/null | awk '{print $3"\t"$3"\t"$0}' | \
        awk '{gsub(/_.*$/,"",$1); gsub(/_.*$/,"",$2)}1' > ${i%.f*}_Index0_subgenome.txt
        awk '{$2=$2"_"$3}1' ${i%.f*}_Index0_subgenome.txt | awk '!h[$2] { g[$2]=$0 } { h[$2]++ } END { for(k in g) print h[k], g[k] }' | \
        awk '!h[$4] { g[$4]=$0 } { h[$4]++ } END { for(k in g) print h[k], g[k] }' | awk '{print $1"\t"$5"\t"$3}' > ${i%.f*}_Index_subgenome.txt
        rm ./alignment/${i%.f*}_redun.sam

        awk '{if($1==3) print $2}' ${i%.f*}_Index_subgenome.txt | grep -Fw -f - <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) | cat <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null | grep '^@') - > ${i%.f*}_del_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam

        awk -v ref3=${ref3%.f*} '{if($1==2 && $3 != ref3) print $2}' ${i%.f*}_Index_subgenome.txt | grep -Fw -f - <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) | cat <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null | grep '^@') - > ${i%.f*}_del_${ref1%.f*}_${ref2%.f*}.sam
        awk -v ref2=${ref2%.f*} '{if($1==2 && $3 != ref2) print $2}' ${i%.f*}_Index_subgenome.txt | grep -Fw -f - <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) | cat <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null | grep '^@') - > ${i%.f*}_del_${ref1%.f*}_${ref3%.f*}.sam
        awk -v ref1=${ref1%.f*} '{if($1==2 && $3 != ref1) print $2}' ${i%.f*}_Index_subgenome.txt | grep -Fw -f - <(zcat ./alignment/${i%.f*}_redun_${ref2%.f*}.sam.gz 2> /dev/null) | cat <(zcat ./alignment/${i%.f*}_redun_${ref2%.f*}.sam.gz 2> /dev/null | grep '^@') - > ${i%.f*}_del_${ref2%.f*}_${ref3%.f*}.sam

        awk -v ref1=${ref1%.f*} '{if($1==1 && $3 == ref1) print $2}' ${i%.f*}_Index_subgenome.txt | grep -Fw -f - <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) | cat <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null | grep '^@') - > ${i%.f*}_del_${ref1%.f*}.sam
        awk -v ref2=${ref2%.f*} '{if($1==1 && $3 == ref2) print $2}' ${i%.f*}_Index_subgenome.txt | grep -Fw -f - <(zcat ./alignment/${i%.f*}_redun_${ref2%.f*}.sam.gz 2> /dev/null) | cat <(zcat ./alignment/${i%.f*}_redun_${ref2%.f*}.sam.gz 2> /dev/null | grep '^@') - > ${i%.f*}_del_${ref2%.f*}.sam
        awk -v ref3=${ref3%.f*} '{if($1==1 && $3 == ref3) print $2}' ${i%.f*}_Index_subgenome.txt | grep -Fw -f - <(zcat ./alignment/${i%.f*}_redun_${ref3%.f*}.sam.gz 2> /dev/null) | cat <(zcat ./alignment/${i%.f*}_redun_${ref3%.f*}.sam.gz 2> /dev/null | grep '^@') - > ${i%.f*}_del_${ref3%.f*}.sam
        rm ${i%.f*}_Index0_subgenome.txt ${i%.f*}_Index_subgenome.txt

        {
          if [[ "$ploidy" -eq 2 ]]; then
            if [[ -z "$downsample_2x" ]]; then
              export downsample=50
            else
              export downsample="$downsample_2x"
            fi
          fi
          if [[ "$ploidy" -eq 4 ]]; then
            if [[ -z "$downsample_4x" ]]; then
              export downsample=100
            else
              export downsample="$downsample_4x"
            fi
          fi
          if [[ "$ploidy" -eq 6 ]]; then
            if [[ -z "$downsample_6x" ]]; then
              export downsample=150
            else
              export downsample=$downsample_6x
            fi
          fi
          if [[ "$ploidy" -eq 8 ]]; then
            if [[ -z "$downsample_8x" ]]; then
              export downsample=200
            else
              export downsample=$downsample_8x
            fi
          fi

          if [[ "$paralogs" == false ]] && [[ "$uniquely_mapped" == true ]]; then
            awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_heading.sam &&
            $samtools view -F4 ${i%.f*}_del_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
            awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | awk '$1==1{print $0}' | \
            awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | \
            awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' | \
            awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_uniq.sam &&
            wait
          fi
          if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == true ]]; then
            awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_heading.sam &&
            $samtools view -F4 ${i%.f*}_del_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
            awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {print $0}}' 2> /dev/null | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_uniqeq.sam
            $samtools view -F4 ${i%.f*}_del_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
            awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 > 0) {print $0}}' 2> /dev/null | cat - ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_uniqeq.sam | \
            awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | \
            awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_uniq.sam &&
            rm ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_uniqeq.sam &&
            wait
          fi
          if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == false ]]; then
            awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_heading.sam &&
            $samtools view -F4 ${i%.f*}_del_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
            awk '!h[$1] { g[$1]=$20 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | \
            awk '$1==1{print $0}' | awk '{print $2}' > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_uniqsingle.sam &&
            $samtools view -F4 ${i%.f*}_del_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
            awk -F '\t' 'NR==FNR{a[$0];next} !($1 in a)' ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_uniqsingle.sam - | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' | \
            awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | \
            awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_uniq.sam &&
            rm ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_uniqsingle.sam &&
            wait
          fi
          wait

          awk '{while ($1-- > 0) print $0}' ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_uniq.sam | awk '{print "seq"NR"_"$0}' | tr -s ' ' | tr ' ' '\t' > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam &&
          rm ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_uniq.sam &&
          awk '{print $4}' ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam | awk '{printf "%d00\n", $1/100}' | paste - ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam | shuf | \
          awk -v max="$downsample" 'a[$1$4]++ < max' | awk '{$1=""}1' | awk '{$1=$1};1' | tr -s " " | tr " " "\t" | awk '{$11=$10}1' | \
          awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$11); print}' | cat ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_heading.sam - > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_downsample.sam &&
          :> ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam &&
          mv ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_downsample.sam ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam &&
          rm ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_heading.sam
          wait
        } & PIDexp_panref=$!
        wait $PIDexp_panref

        if [[ "${ploidy_ref1}" && "${ploidy_ref2}" ]]; then
          {
            if [[ "$ploidy" -eq 2 ]]; then
              if [[ -z "$downsample_2x" ]]; then
                export downsample=50
              else
                export downsample="$downsample_2x"
              fi
            fi
            if [[ "$ploidy" -eq 4 ]]; then
              if [[ -z "$downsample_4x" ]]; then
                export downsample=100
              else
                export downsample="$downsample_4x"
              fi
            fi
            if [[ "$ploidy" -eq 6 ]]; then
              if [[ -z "$downsample_6x" ]]; then
                export downsample=150
              else
                export downsample=$downsample_6x
              fi
            fi
            if [[ "$ploidy" -eq 8 ]]; then
              if [[ -z "$downsample_8x" ]]; then
                export downsample=200
              else
                export downsample=$downsample_8x
              fi
            fi

            if [[ "$paralogs" == false ]] && [[ "$uniquely_mapped" == true ]]; then
              awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_heading.sam &&
              $samtools view -F4 ${i%.f*}_del_${ref1%.f*}_${ref2%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
              awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | awk '$1==1{print $0}' | \
              awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | \
              awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' | \
              awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_uniq.sam &&
              wait
            fi
            if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == true ]]; then
              awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_heading.sam &&
              $samtools view -F4 ${i%.f*}_del_${ref1%.f*}_${ref2%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
              awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {print $0}}' 2> /dev/null | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_uniqeq.sam
              $samtools view -F4 ${i%.f*}_del_${ref1%.f*}_${ref2%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
              awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 > 0) {print $0}}' 2> /dev/null | cat - ${i%.f*}_${ref1%.f*}_${ref2%.f*}_uniqeq.sam | \
              awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | \
              awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_uniq.sam &&
              rm ${i%.f*}_${ref1%.f*}_${ref2%.f*}_uniqeq.sam &&
              wait
            fi
            if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == false ]]; then
              awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_heading.sam &&
              $samtools view -F4 ${i%.f*}_del_${ref1%.f*}_${ref2%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
              awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | \
              awk '$1==1{print $0}' | awk '{print $2}' > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_uniqsingle.sam &&
              $samtools view -F4 ${i%.f*}_del_${ref1%.f*}_${ref2%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
              awk -F '\t' 'NR==FNR{a[$0];next} !($1 in a)' ${i%.f*}_${ref1%.f*}_${ref2%.f*}_uniqsingle.sam - | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' | \
              awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | \
              awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_uniq.sam &&
              rm ${i%.f*}_${ref1%.f*}_${ref2%.f*}_uniqsingle.sam &&
              wait
            fi
            wait

            awk '{while ($1-- > 0) print $0}' ${i%.f*}_${ref1%.f*}_${ref2%.f*}_uniq.sam | awk '{print "seq"NR"_"$0}' | tr -s ' ' | tr ' ' '\t' > ${i%.f*}_${ref1%.f*}_${ref2%.f*}.sam &&
            rm ${i%.f*}_${ref1%.f*}_${ref2%.f*}_uniq.sam &&
            awk '{print $4}' ${i%.f*}_${ref1%.f*}_${ref2%.f*}.sam | awk '{printf "%d00\n", $1/100}' | paste - ${i%.f*}_${ref1%.f*}_${ref2%.f*}.sam | shuf | \
            awk -v max="$downsample" 'a[$1$4]++ < max' | awk '{$1=""}1' | awk '{$1=$1};1' | tr -s " " | tr " " "\t" | awk '{$11=$10}1' | \
            awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$11); print}' | cat ${i%.f*}_${ref1%.f*}_${ref2%.f*}_heading.sam - > ${i%.f*}_${ref1%.f*}_${ref2%.f*}_downsample.sam &&
            :> ${i%.f*}_${ref1%.f*}_${ref2%.f*}.sam &&
            mv ${i%.f*}_${ref1%.f*}_${ref2%.f*}_downsample.sam ${i%.f*}_${ref1%.f*}_${ref2%.f*}.sam &&
            rm ${i%.f*}_${ref1%.f*}_${ref2%.f*}_heading.sam
            wait
          } & PIDexp_ref1_ref2=$!
          wait $PIDexp_ref1_ref2
        fi

        if [[ "${ploidy_ref1}" && "${ploidy_ref3}" ]]; then
          {
            if [[ "$ploidy" -eq 2 ]]; then
              if [[ -z "$downsample_2x" ]]; then
                export downsample=50
              else
                export downsample="$downsample_2x"
              fi
            fi
            if [[ "$ploidy" -eq 4 ]]; then
              if [[ -z "$downsample_4x" ]]; then
                export downsample=100
              else
                export downsample="$downsample_4x"
              fi
            fi
            if [[ "$ploidy" -eq 6 ]]; then
              if [[ -z "$downsample_6x" ]]; then
                export downsample=150
              else
                export downsample=$downsample_6x
              fi
            fi
            if [[ "$ploidy" -eq 8 ]]; then
              if [[ -z "$downsample_8x" ]]; then
                export downsample=200
              else
                export downsample=$downsample_8x
              fi
            fi

            if [[ "$paralogs" == false ]] && [[ "$uniquely_mapped" == true ]]; then
              awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref1%.f*}_${ref3%.f*}_heading.sam &&
              $samtools view -F4 ${i%.f*}_del_${ref1%.f*}_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
              awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | awk '$1==1{print $0}' | \
              awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | \
              awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' | \
              awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref1%.f*}_${ref3%.f*}_uniq.sam &&
              wait
            fi
            if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == true ]]; then
              awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref1%.f*}_${ref3%.f*}_heading.sam &&
              $samtools view -F4 ${i%.f*}_del_${ref1%.f*}_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
              awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {print $0}}' 2> /dev/null | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' > ${i%.f*}_${ref1%.f*}_${ref3%.f*}_uniqeq.sam
              $samtools view -F4 ${i%.f*}_del_${ref1%.f*}_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
              awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 > 0) {print $0}}' 2> /dev/null | cat - ${i%.f*}_${ref1%.f*}_${ref3%.f*}_uniqeq.sam | \
              awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | \
              awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref1%.f*}_${ref3%.f*}_uniq.sam &&
              rm ${i%.f*}_${ref1%.f*}_${ref3%.f*}_uniqeq.sam &&
              wait
            fi
            if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == false ]]; then
              awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref1%.f*}_${ref3%.f*}_heading.sam &&
              $samtools view -F4 ${i%.f*}_del_${ref1%.f*}_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
              awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | \
              awk '$1==1{print $0}' | awk '{print $2}' > ${i%.f*}_${ref1%.f*}_${ref3%.f*}_uniqsingle.sam &&
              $samtools view -F4 ${i%.f*}_del_${ref1%.f*}_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
              awk -F '\t' 'NR==FNR{a[$0];next} !($1 in a)' ${i%.f*}_${ref1%.f*}_${ref3%.f*}_uniqsingle.sam - | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' | \
              awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | \
              awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref1%.f*}_${ref3%.f*}_uniq.sam &&
              rm ${i%.f*}_${ref1%.f*}_${ref3%.f*}_uniqsingle.sam &&
              wait
            fi
            wait

            awk '{while ($1-- > 0) print $0}' ${i%.f*}_${ref1%.f*}_${ref3%.f*}_uniq.sam | awk '{print "seq"NR"_"$0}' | tr -s ' ' | tr ' ' '\t' > ${i%.f*}_${ref1%.f*}_${ref3%.f*}.sam &&
            rm ${i%.f*}_${ref1%.f*}_${ref3%.f*}_uniq.sam &&
            awk '{print $4}' ${i%.f*}_${ref1%.f*}_${ref3%.f*}.sam | awk '{printf "%d00\n", $1/100}' | paste - ${i%.f*}_${ref1%.f*}_${ref3%.f*}.sam | shuf | \
            awk -v max="$downsample" 'a[$1$4]++ < max' | awk '{$1=""}1' | awk '{$1=$1};1' | tr -s " " | tr " " "\t" | awk '{$11=$10}1' | \
            awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$11); print}' | cat ${i%.f*}_${ref1%.f*}_${ref3%.f*}_heading.sam - > ${i%.f*}_${ref1%.f*}_${ref3%.f*}_downsample.sam &&
            :> ${i%.f*}_${ref1%.f*}_${ref3%.f*}.sam &&
            mv ${i%.f*}_${ref1%.f*}_${ref3%.f*}_downsample.sam ${i%.f*}_${ref1%.f*}_${ref3%.f*}.sam &&
            rm ${i%.f*}_${ref1%.f*}_${ref3%.f*}_heading.sam
            wait
          } & PIDexp_ref1_ref3=$!
          wait $PIDexp_ref1_ref3
        fi

        if [[ "${ploidy_ref2}" && "${ploidy_ref3}" ]]; then
          {
            if [[ "$ploidy" -eq 2 ]]; then
              if [[ -z "$downsample_2x" ]]; then
                export downsample=50
              else
                export downsample="$downsample_2x"
              fi
            fi
            if [[ "$ploidy" -eq 4 ]]; then
              if [[ -z "$downsample_4x" ]]; then
                export downsample=100
              else
                export downsample="$downsample_4x"
              fi
            fi
            if [[ "$ploidy" -eq 6 ]]; then
              if [[ -z "$downsample_6x" ]]; then
                export downsample=150
              else
                export downsample=$downsample_6x
              fi
            fi
            if [[ "$ploidy" -eq 8 ]]; then
              if [[ -z "$downsample_8x" ]]; then
                export downsample=200
              else
                export downsample=$downsample_8x
              fi
            fi

            if [[ "$paralogs" == false ]] && [[ "$uniquely_mapped" == true ]]; then
              awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref2%.f*}_${ref3%.f*}_heading.sam &&
              $samtools view -F4 ${i%.f*}_del_${ref2%.f*}_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
              awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | awk '$1==1{print $0}' | \
              awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | \
              awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' | \
              awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref2%.f*}_${ref3%.f*}_uniq.sam &&
              wait
            fi
            if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == true ]]; then
              awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref2%.f*}_${ref3%.f*}_heading.sam &&
              $samtools view -F4 ${i%.f*}_del_${ref2%.f*}_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
              awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {print $0}}' 2> /dev/null | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' > ${i%.f*}_${ref2%.f*}_${ref3%.f*}_uniqeq.sam
              $samtools view -F4 ${i%.f*}_del_${ref2%.f*}_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
              awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 > 0) {print $0}}' 2> /dev/null | cat - ${i%.f*}_${ref2%.f*}_${ref3%.f*}_uniqeq.sam | \
              awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | \
              awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref2%.f*}_${ref3%.f*}_uniq.sam &&
              rm ${i%.f*}_${ref2%.f*}_${ref3%.f*}_uniqeq.sam &&
              wait
            fi
            if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == false ]]; then
              awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref2%.f*}_${ref3%.f*}_heading.sam &&
              $samtools view -F4 ${i%.f*}_del_${ref2%.f*}_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
              awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | \
              awk '$1==1{print $0}' | awk '{print $2}' > ${i%.f*}_${ref2%.f*}_${ref3%.f*}_uniqsingle.sam &&
              $samtools view -F4 ${i%.f*}_del_${ref2%.f*}_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
              awk -F '\t' 'NR==FNR{a[$0];next} !($1 in a)' ${i%.f*}_${ref2%.f*}_${ref3%.f*}_uniqsingle.sam - | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' | \
              awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | \
              awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref2%.f*}_${ref3%.f*}_uniq.sam &&
              rm ${i%.f*}_${ref2%.f*}_${ref3%.f*}_uniqsingle.sam &&
              wait
            fi
            wait

            awk '{while ($1-- > 0) print $0}' ${i%.f*}_${ref2%.f*}_${ref3%.f*}_uniq.sam | awk '{print "seq"NR"_"$0}' | tr -s ' ' | tr ' ' '\t' > ${i%.f*}_${ref2%.f*}_${ref3%.f*}.sam &&
            rm ${i%.f*}_${ref2%.f*}_${ref3%.f*}_uniq.sam &&
            awk '{print $4}' ${i%.f*}_${ref2%.f*}_${ref3%.f*}.sam | awk '{printf "%d00\n", $1/100}' | paste - ${i%.f*}_${ref2%.f*}_${ref3%.f*}.sam | shuf | \
            awk -v max="$downsample" 'a[$1$4]++ < max' | awk '{$1=""}1' | awk '{$1=$1};1' | tr -s " " | tr " " "\t" | awk '{$11=$10}1' | \
            awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$11); print}' | cat ${i%.f*}_${ref2%.f*}_${ref3%.f*}_heading.sam - > ${i%.f*}_${ref2%.f*}_${ref3%.f*}_downsample.sam &&
            :> ${i%.f*}_${ref2%.f*}_${ref3%.f*}.sam &&
            mv ${i%.f*}_${ref2%.f*}_${ref3%.f*}_downsample.sam ${i%.f*}_${ref2%.f*}_${ref3%.f*}.sam &&
            rm ${i%.f*}_${ref2%.f*}_${ref3%.f*}_heading.sam
            wait
          } & PIDexp_ref2_ref3=$!
          wait $PIDexp_ref2_ref3
        fi

        if [[ "${ploidy_ref1}" ]]; then
          {
            if [[ "$ploidy_ref1" ]]; then
              if [[ "$ploidy_ref1" -eq 2 ]]; then
                if [[ -z "$downsample_2x" ]]; then
                  export downsample=50
                else
                  export downsample="$downsample_2x"
                fi
              fi
              if [[ "$ploidy_ref1" -eq 4 ]]; then
                if [[ -z "$downsample_4x" ]]; then
                  export downsample=100
                else
                  export downsample="$downsample_4x"
                fi
              fi
              if [[ "$ploidy_ref1" -eq 6 ]]; then
                if [[ -z "$downsample_6x" ]]; then
                  export downsample=150
                else
                  export downsample=$downsample_6x
                fi
              fi
              if [[ "$ploidy_ref1" -eq 8 ]]; then
                if [[ -z "$downsample_8x" ]]; then
                  export downsample=200
                else
                  export downsample=$downsample_8x
                fi
              fi
              if [[ "$paralogs" == false ]] && [[ "$uniquely_mapped" == true ]]; then
                awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref1%.f*}_heading.sam &&
                $samtools view -F4 ${i%.f*}_del_${ref1%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
                awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | awk '$1==1{print $0}' | \
                awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | \
                awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' | \
                awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref1%.f*}_uniq.sam &&
                wait
              fi
              if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == true ]]; then
                awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref1%.f*}_heading.sam &&
                $samtools view -F4 ${i%.f*}_del_${ref1%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
                awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {print $0}}' 2> /dev/null | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' > ${i%.f*}_${ref1%.f*}_uniqeq.sam
                $samtools view -F4 ${i%.f*}_del_${ref1%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
                awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 > 0) {print $0}}' 2> /dev/null | cat - ${i%.f*}_${ref1%.f*}_uniqeq.sam | \
                awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | \
                awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref1%.f*}_uniq.sam &&
                rm ${i%.f*}_${ref1%.f*}_uniqeq.sam &&
                wait
              fi
              if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == false ]]; then
                awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref1%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref1%.f*}_heading.sam &&
                $samtools view -F4 ${i%.f*}_del_${ref1%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
                awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | \
                awk '$1==1{print $0}' | awk '{print $2}' > ${i%.f*}_${ref1%.f*}_uniqsingle.sam &&
                $samtools view -F4 ${i%.f*}_del_${ref1%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
                awk -F '\t' 'NR==FNR{a[$0];next} !($1 in a)' ${i%.f*}_${ref1%.f*}_uniqsingle.sam - | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' | \
                awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | \
                awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref1%.f*}_uniq.sam &&
                rm ${i%.f*}_${ref1%.f*}_uniqsingle.sam &&
                wait
              fi
              wait

              awk '{while ($1-- > 0) print $0}' ${i%.f*}_${ref1%.f*}_uniq.sam | awk '{print "seq"NR"_"$0}' | tr -s ' ' | tr ' ' '\t' > ${i%.f*}_${ref1%.f*}.sam &&
              rm ${i%.f*}_${ref1%.f*}_uniq.sam &&
              awk '{print $4}' ${i%.f*}_${ref1%.f*}.sam | awk '{printf "%d00\n", $1/100}' | paste - ${i%.f*}_${ref1%.f*}.sam | shuf | \
              awk -v max="$downsample" 'a[$1$4]++ < max' | awk '{$1=""}1' | awk '{$1=$1};1' | tr -s " " | tr " " "\t" | awk '{$11=$10}1' | \
              awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$11); print}' | cat ${i%.f*}_${ref1%.f*}_heading.sam - > ${i%.f*}_${ref1%.f*}_downsample.sam &&
              :> ${i%.f*}_${ref1%.f*}.sam &&
              mv ${i%.f*}_${ref1%.f*}_downsample.sam ${i%.f*}_${ref1%.f*}.sam &&
              rm ${i%.f*}_${ref1%.f*}_heading.sam
              wait
            fi
          } & PIDexp_ref1=$!
          wait $PIDexp_ref1
        fi

        if [[ "${ploidy_ref2}" ]]; then
          {
            if [[ "$ploidy_ref2" ]]; then
              if [[ "$ploidy_ref2" -eq 2 ]]; then
                if [[ -z "$downsample_2x" ]]; then
                  export downsample=100
                else
                  export downsample="$downsample_2x"
                fi
              fi
              if [[ "$ploidy_ref2" -eq 4 ]]; then
                if [[ -z "$downsample_4x" ]]; then
                  export downsample=200
                else
                  export downsample="$downsample_4x"
                fi
              fi
              if [[ "$ploidy_ref2" -eq 6 ]]; then
                if [[ -z "$downsample_6x" ]]; then
                  export downsample=300
                else
                  export downsample=$downsample_6x
                fi
              fi
              if [[ "$ploidy_ref2" -eq 8 ]]; then
                if [[ -z "$downsample_8x" ]]; then
                  export downsample=400
                else
                  export downsample=$downsample_8x
                fi
              fi
              if [[ "$paralogs" == false ]] && [[ "$uniquely_mapped" == true ]]; then
                awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref2%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref2%.f*}_heading.sam &&
                $samtools view -F4 ${i%.f*}_del_${ref2%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
                awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | awk '$1==1{print $0}' | \
                awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | \
                awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' | \
                awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref2%.f*}_uniq.sam &&
                wait
              fi
              if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == true ]]; then
                awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref2%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref2%.f*}_heading.sam &&
                $samtools view -F4 ${i%.f*}_del_${ref2%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
                awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {print $0}}' 2> /dev/null | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' > ${i%.f*}_${ref2%.f*}_uniqeq.sam
                $samtools view -F4 ${i%.f*}_del_${ref2%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
                awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 > 0) {print $0}}' 2> /dev/null | cat - ${i%.f*}_${ref2%.f*}_uniqeq.sam | \
                awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | \
                awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref2%.f*}_uniq.sam &&
                rm ${i%.f*}_${ref2%.f*}_uniqeq.sam &&
                wait
              fi
              if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == false ]]; then
                awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref2%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref2%.f*}_heading.sam &&
                $samtools view -F4 ${i%.f*}_del_${ref2%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
                awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | \
                awk '$1==1{print $0}' | awk '{print $2}' > ${i%.f*}_${ref2%.f*}_uniqsingle.sam &&
                $samtools view -F4 ${i%.f*}_del_${ref2%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
                awk -F '\t' 'NR==FNR{a[$0];next} !($1 in a)' ${i%.f*}_${ref2%.f*}_uniqsingle.sam - | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' | \
                awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | \
                awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref2%.f*}_uniq.sam &&
                rm ${i%.f*}_${ref2%.f*}_uniqsingle.sam &&
                wait
              fi
              wait

              awk '{while ($1-- > 0) print $0}' ${i%.f*}_${ref2%.f*}_uniq.sam | awk '{print "seq"NR"_"$0}' | tr -s ' ' | tr ' ' '\t' > ${i%.f*}_${ref2%.f*}.sam &&
              rm ${i%.f*}_${ref2%.f*}_uniq.sam &&
              awk '{print $4}' ${i%.f*}_${ref2%.f*}.sam | awk '{printf "%d00\n", $1/100}' | paste - ${i%.f*}_${ref2%.f*}.sam | shuf | \
              awk -v max="$downsample" 'a[$1$4]++ < max' | awk '{$1=""}1' | awk '{$1=$1};1' | tr -s " " | tr " " "\t" | awk '{$11=$10}1' | \
              awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$11); print}' | cat ${i%.f*}_${ref2%.f*}_heading.sam - > ${i%.f*}_${ref2%.f*}_downsample.sam &&
              :> ${i%.f*}_${ref2%.f*}.sam &&
              mv ${i%.f*}_${ref2%.f*}_downsample.sam ${i%.f*}_${ref2%.f*}.sam &&
              rm ${i%.f*}_${ref2%.f*}_heading.sam
              wait
            fi
          } & PIDexp_ref2=$!
          wait $PIDexp_ref2
        fi

        if [[ "${ploidy_ref3}" ]]; then
          {
            if [[ "$ploidy_ref3" ]]; then
              if [[ "$ploidy_ref3" -eq 2 ]]; then
                if [[ -z "$downsample_2x" ]]; then
                  export downsample=100
                else
                  export downsample="$downsample_2x"
                fi
              fi
              if [[ "$ploidy_ref3" -eq 4 ]]; then
                if [[ -z "$downsample_4x" ]]; then
                  export downsample=200
                else
                  export downsample="$downsample_4x"
                fi
              fi
              if [[ "$ploidy_ref3" -eq 6 ]]; then
                if [[ -z "$downsample_6x" ]]; then
                  export downsample=300
                else
                  export downsample=$downsample_6x
                fi
              fi
              if [[ "$ploidy_ref3" -eq 8 ]]; then
                if [[ -z "$downsample_8x" ]]; then
                  export downsample=400
                else
                  export downsample=$downsample_8x
                fi
              fi
              if [[ "$paralogs" == false ]] && [[ "$uniquely_mapped" == true ]]; then
                awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref3%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref3%.f*}_heading.sam &&
                $samtools view -F4 ${i%.f*}_del_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
                awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | awk '$1==1{print $0}' | \
                awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | \
                awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' | \
                awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref3%.f*}_uniq.sam &&
                wait
              fi
              if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == true ]]; then
                awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref3%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref3%.f*}_heading.sam &&
                $samtools view -F4 ${i%.f*}_del_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
                awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {print $0}}' 2> /dev/null | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' > ${i%.f*}_${ref3%.f*}_uniqeq.sam
                $samtools view -F4 ${i%.f*}_del_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
                awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 > 0) {print $0}}' 2> /dev/null | cat - ${i%.f*}_${ref3%.f*}_uniqeq.sam | \
                awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | \
                awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref3%.f*}_uniq.sam &&
                rm ${i%.f*}_${ref3%.f*}_uniqeq.sam &&
                wait
              fi
              if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == false ]]; then
                awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun_${ref3%.f*}.sam.gz 2> /dev/null) > ${i%.f*}_${ref3%.f*}_heading.sam &&
                $samtools view -F4 ${i%.f*}_del_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
                awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | \
                awk '$1==1{print $0}' | awk '{print $2}' > ${i%.f*}_${ref3%.f*}_uniqsingle.sam &&
                $samtools view -F4 ${i%.f*}_del_${ref3%.f*}.sam | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
                awk -F '\t' 'NR==FNR{a[$0];next} !($1 in a)' ${i%.f*}_${ref3%.f*}_uniqsingle.sam - | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {$5=60}}1' | \
                awk '{gsub(/_pe-/,"_pe-\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | \
                awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_${ref3%.f*}_uniq.sam &&
                rm ${i%.f*}_${ref3%.f*}_uniqsingle.sam &&
                wait
              fi
              wait

              awk '{while ($1-- > 0) print $0}' ${i%.f*}_${ref3%.f*}_uniq.sam | awk '{print "seq"NR"_"$0}' | tr -s ' ' | tr ' ' '\t' > ${i%.f*}_${ref3%.f*}.sam &&
              rm ${i%.f*}_${ref3%.f*}_uniq.sam &&
              awk '{print $4}' ${i%.f*}_${ref3%.f*}.sam | awk '{printf "%d00\n", $1/100}' | paste - ${i%.f*}_${ref3%.f*}.sam | shuf | \
              awk -v max="$downsample" 'a[$1$4]++ < max' | awk '{$1=""}1' | awk '{$1=$1};1' | tr -s " " | tr " " "\t" | awk '{$11=$10}1' | \
              awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$11); print}' | cat ${i%.f*}_${ref3%.f*}_heading.sam - > ${i%.f*}_${ref3%.f*}_downsample.sam &&
              :> ${i%.f*}_${ref3%.f*}.sam &&
              mv ${i%.f*}_${ref3%.f*}_downsample.sam ${i%.f*}_${ref3%.f*}.sam &&
              rm ${i%.f*}_${ref3%.f*}_heading.sam
              wait
            fi
          } & PIDexp_ref3=$!
          wait $PIDexp_ref3
        fi

        j="${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam"
        $java $Xmxp -XX:ParallelGCThreads=$prepthreads -Djava.io.tmpdir=${projdir}/preprocess/tmp -jar $picard SortSam I=$j O=${j%.sam*}.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=${projdir}/preprocess/tmp && \
        $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard BuildBamIndex INPUT=${j%.sam*}.bam VALIDATION_STRINGENCY=LENIENT && \
        $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard AddOrReplaceReadGroups I=${j%.sam*}.bam O=${j%.sam*}_precall.bam RGLB=${i%.f*} RGPL=illumina RGPU=run RGSM=${i%.f*} VALIDATION_STRINGENCY=LENIENT && \
        $samtools index ${j%.sam*}_precall.bam &&
        rm $j ${j%.sam*}.bam ${j%.sam*}.bai &&
        wait

        if [[ "${ploidy_ref1}" && "${ploidy_ref2}" ]]; then
          j="${i%.f*}_${ref1%.f*}_${ref2%.f*}.sam"
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -Djava.io.tmpdir=${projdir}/preprocess/tmp -jar $picard SortSam I=$j O=${j%.sam*}.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=${projdir}/preprocess/tmp && \
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard BuildBamIndex INPUT=${j%.sam*}.bam VALIDATION_STRINGENCY=LENIENT && \
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard AddOrReplaceReadGroups I=${j%.sam*}.bam O=${j%.sam*}_precall.bam RGLB=${i%.f*} RGPL=illumina RGPU=run RGSM=${i%.f*} VALIDATION_STRINGENCY=LENIENT && \
          $samtools index ${j%.sam*}_precall.bam &&
          rm $j ${j%.sam*}.bam ${j%.sam*}.bai &&
          wait
        fi

        if [[ "${ploidy_ref1}" && "${ploidy_ref3}" ]]; then
          j="${i%.f*}_${ref1%.f*}_${ref3%.f*}.sam"
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -Djava.io.tmpdir=${projdir}/preprocess/tmp -jar $picard SortSam I=$j O=${j%.sam*}.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=${projdir}/preprocess/tmp && \
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard BuildBamIndex INPUT=${j%.sam*}.bam VALIDATION_STRINGENCY=LENIENT && \
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard AddOrReplaceReadGroups I=${j%.sam*}.bam O=${j%.sam*}_precall.bam RGLB=${i%.f*} RGPL=illumina RGPU=run RGSM=${i%.f*} VALIDATION_STRINGENCY=LENIENT && \
          $samtools index ${j%.sam*}_precall.bam &&
          rm $j ${j%.sam*}.bam ${j%.sam*}.bai &&
          wait
        fi

        if [[ "${ploidy_ref2}" && "${ploidy_ref3}" ]]; then
          j="${i%.f*}_${ref2%.f*}_${ref3%.f*}.sam"
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -Djava.io.tmpdir=${projdir}/preprocess/tmp -jar $picard SortSam I=$j O=${j%.sam*}.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=${projdir}/preprocess/tmp && \
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard BuildBamIndex INPUT=${j%.sam*}.bam VALIDATION_STRINGENCY=LENIENT && \
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard AddOrReplaceReadGroups I=${j%.sam*}.bam O=${j%.sam*}_precall.bam RGLB=${i%.f*} RGPL=illumina RGPU=run RGSM=${i%.f*} VALIDATION_STRINGENCY=LENIENT && \
          $samtools index ${j%.sam*}_precall.bam &&
          rm $j ${j%.sam*}.bam ${j%.sam*}.bai &&
          wait
        fi

        if [[ "${ploidy_ref1}" ]]; then
          j="${i%.f*}_${ref1%.f*}.sam"
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -Djava.io.tmpdir=${projdir}/preprocess/tmp -jar $picard SortSam I=$j O=${j%.sam*}.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=${projdir}/preprocess/tmp && \
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard BuildBamIndex INPUT=${j%.sam*}.bam VALIDATION_STRINGENCY=LENIENT && \
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard AddOrReplaceReadGroups I=${j%.sam*}.bam O=${j%.sam*}_precall.bam RGLB=${i%.f*} RGPL=illumina RGPU=run RGSM=${i%.f*} VALIDATION_STRINGENCY=LENIENT && \
          $samtools index ${j%.sam*}_precall.bam &&
          rm $j ${j%.sam*}.bam ${j%.sam*}.bai &&
          wait
        fi

        if [[ "${ploidy_ref2}" ]]; then
          j="${i%.f*}_${ref2%.f*}.sam"
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -Djava.io.tmpdir=${projdir}/preprocess/tmp -jar $picard SortSam I=$j O=${j%.sam*}.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=${projdir}/preprocess/tmp && \
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard BuildBamIndex INPUT=${j%.sam*}.bam VALIDATION_STRINGENCY=LENIENT && \
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard AddOrReplaceReadGroups I=${j%.sam*}.bam O=${j%.sam*}_precall.bam RGLB=${i%.f*} RGPL=illumina RGPU=run RGSM=${i%.f*} VALIDATION_STRINGENCY=LENIENT && \
          $samtools index ${j%.sam*}_precall.bam &&
          rm $j ${j%.sam*}.bam ${j%.sam*}.bai &&
          wait
        fi

        if [[ "${ploidy_ref3}" ]]; then
          j="${i%.f*}_${ref3%.f*}.sam"
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -Djava.io.tmpdir=${projdir}/preprocess/tmp -jar $picard SortSam I=$j O=${j%.sam*}.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=${projdir}/preprocess/tmp && \
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard BuildBamIndex INPUT=${j%.sam*}.bam VALIDATION_STRINGENCY=LENIENT && \
          $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard AddOrReplaceReadGroups I=${j%.sam*}.bam O=${j%.sam*}_precall.bam RGLB=${i%.f*} RGPL=illumina RGPU=run RGSM=${i%.f*} VALIDATION_STRINGENCY=LENIENT && \
          $samtools index ${j%.sam*}_precall.bam &&
          rm $j ${j%.sam*}.bam ${j%.sam*}.bai &&
          wait
        fi

        if [[ $nodes -gt 1 ]]; then cp /tmp/${samples_list%.txt}/preprocess/${j%.sam*}_precall.bam* ${projdir}/preprocess/; fi
        wait
      fi ) &
      if [[ $(jobs -r -p | wc -l) -ge $prepN ]]; then
        wait
      fi
    done < <(cat ${projdir}/${samples_list})
    wait && touch ${projdir}/precall_done_${samples_list}
    find . -maxdepth 1 -type f -name "*" | grep -v 'precall' | xargs rm 2> /dev/null
    wait


    cd ${projdir}/samples
    if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
      precall=$(ls ${projdir}/precall_done_samples_list_node_* | wc -l)
      while [[ "$precall" -lt $nodes ]]; do sleep 300; precall=$(ls ${projdir}/precall_done_samples_list_node_* | wc -l); done
      sleep $((RANDOM % 2))
      if [[ $precall == $nodes ]] && test ! -f ${projdir}/alignment_summaries/refgenome_paralogs.txt; then
        cd ${projdir}/alignment_summaries
        cat *_summ_${ref1%.f*}.txt > alignment_summaries_reads_${ref1%.f*}.txt &&
        rm *_summ_${ref1%.f*}.txt &&
        # Total number of reads per samples
        awk '/###---/ || /QC-passed/{print}' alignment_summaries_reads_${ref1%.f*}.txt | cut -d\+ -f1 | tr -d '\n' | \
        awk  'gsub(/---###/, "\t", $0)' | awk  'gsub(/###---/, "", $0)' | tr ' ' '\n' > total_reads_${ref1%.f*}.txt &&
        # Total number of mapped reads per samples
        cat alignment_summaries_reads_${ref1%.f*}.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /0_mapped/{print}' |\
        tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
        awk 'gsub("\\+0mapped", "\t", $0)' | cut -d\: -f1 > total_reads_mapped_${ref1%.f*}.txt &&
        # Total number of mapped paired reads per samples
        cat alignment_summaries_reads_${ref1%.f*}.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /properly_paired/{print}' |\
        tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
        awk 'gsub("\\+0properlypaired", "\t", $0)' | cut -d\: -f1 > total_reads_paired_${ref1%.f*}.txt &&
        echo -e "Samples\tTotal\tMapped\tPerc_Mapped\tPE_Mapped\t%_PE_Mapped" > summary_precall_${ref1%.f*}.txt &&
        awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' total_reads_mapped_${ref1%.f*}.txt  total_reads_${ref1%.f*}.txt  | \
        awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' total_reads_paired_${ref1%.f*}.txt - | awk '{gsub(/ /,"\t"); print $0}' | \
        cat summary_precall_${ref1%.f*}.txt - | awk '{print $1"\t"$2"\t"$3"\t"$4}' > Tabulated_Alignment_Read_Summaries_${ref1%.f*}.txt &&
        rm total_* summary_precall_${ref1%.f*}.txt &> /dev/null &&
        rm ${projdir}/samples/metrics.txt ${projdir}/preprocess/metrics.txt &> /dev/null &&
        wait

        cat *_summ_${ref2%.f*}.txt > alignment_summaries_reads_${ref2%.f*}.txt &&
        rm *_summ_${ref2%.f*}.txt &&
        # Total number of reads per samples
        awk '/###---/ || /QC-passed/{print}' alignment_summaries_reads_${ref2%.f*}.txt | cut -d\+ -f1 | tr -d '\n' | \
        awk  'gsub(/---###/, "\t", $0)' | awk  'gsub(/###---/, "", $0)' | tr ' ' '\n' > total_reads_${ref2%.f*}.txt &&
        # Total number of mapped reads per samples
        cat alignment_summaries_reads_${ref2%.f*}.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /0_mapped/{print}' |\
        tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
        awk 'gsub("\\+0mapped", "\t", $0)' | cut -d\: -f1 > total_reads_mapped_${ref2%.f*}.txt &&
        # Total number of mapped paired reads per samples
        cat alignment_summaries_reads_${ref2%.f*}.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /properly_paired/{print}' |\
        tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
        awk 'gsub("\\+0properlypaired", "\t", $0)' | cut -d\: -f1 > total_reads_paired_${ref2%.f*}.txt &&
        echo -e "Samples\tTotal\tMapped\tPerc_Mapped\tPE_Mapped\t%_PE_Mapped" > summary_precall_${ref2%.f*}.txt &&
        awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' total_reads_mapped_${ref2%.f*}.txt  total_reads_${ref2%.f*}.txt  | \
        awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' total_reads_paired_${ref2%.f*}.txt - | awk '{gsub(/ /,"\t"); print $0}' | \
        cat summary_precall_${ref2%.f*}.txt - | awk '{print $1"\t"$2"\t"$3"\t"$4}' > Tabulated_Alignment_Read_Summaries_${ref2%.f*}.txt &&
        rm total_* summary_precall_${ref2%.f*}.txt &> /dev/null &&
        rm ${projdir}/samples/metrics.txt ${projdir}/preprocess/metrics.txt &> /dev/null &&
        wait

        cat *_summ_${ref3%.f*}.txt > alignment_summaries_reads_${ref3%.f*}.txt &&
        rm *_summ_${ref3%.f*}.txt &&
        # Total number of reads per samples
        awk '/###---/ || /QC-passed/{print}' alignment_summaries_reads_${ref3%.f*}.txt | cut -d\+ -f1 | tr -d '\n' | \
        awk  'gsub(/---###/, "\t", $0)' | awk  'gsub(/###---/, "", $0)' | tr ' ' '\n' > total_reads_${ref3%.f*}.txt &&
        # Total number of mapped reads per samples
        cat alignment_summaries_reads_${ref3%.f*}.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /0_mapped/{print}' |\
        tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
        awk 'gsub("\\+0mapped", "\t", $0)' | cut -d\: -f1 > total_reads_mapped_${ref3%.f*}.txt &&
        # Total number of mapped paired reads per samples
        cat alignment_summaries_reads_${ref3%.f*}.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /properly_paired/{print}' |\
        tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
        awk 'gsub("\\+0properlypaired", "\t", $0)' | cut -d\: -f1 > total_reads_paired_${ref3%.f*}.txt &&
        echo -e "Samples\tTotal\tMapped\tPerc_Mapped\tPE_Mapped\t%_PE_Mapped" > summary_precall_${ref3%.f*}.txt &&
        awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' total_reads_mapped_${ref3%.f*}.txt  total_reads_${ref3%.f*}.txt  | \
        awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' total_reads_paired_${ref3%.f*}.txt - | awk '{gsub(/ /,"\t"); print $0}' | \
        cat summary_precall_${ref3%.f*}.txt - | awk '{print $1"\t"$2"\t"$3"\t"$4}' > Tabulated_Alignment_Read_Summaries_${ref3%.f*}.txt &&
        rm total_* summary_precall_${ref3%.f*}.txt &> /dev/null &&
        rm ${projdir}/samples/metrics.txt ${projdir}/preprocess/metrics.txt &> /dev/null &&
        wait

      fi
      cd ${projdir}/alignment_summaries
      printf "Sample\tGenome_Coverage(percentage)\n" > summary_genomecov.txt
      genome_size=$(awk '{print $3}' ../refgenomes/${ref1%.f*}.dict | awk '{gsub(/LN:/,"");}1' | awk '{s+=$1}END{print s}')
      for i in ../preprocess/*bam; do
        cov=$($bedtools genomecov -ibam $i -bga | awk -v pat=$genome_size '{s+=$4}END{print (s/pat)*100}')
        printf "${i%_*_*}\t$cov\n"  | awk '{gsub(/..\/preprocess\/processed\//,"");}1' >> summary_genomecov.txt
      done
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
  			echo -e "${magenta}- read alignments and alignment post-processinga already performed ${white}\n"
  		else
  			echo -e "${magenta}- performing read alignments and alignment post-processing ${white}\n"
  			time main &>> log.out
  		fi
  	fi
  fi
  if [ "$walkaway" == true ]; then
  	if [ "$alignments" == 1 ]; then
  		if test -f ${projdir}/alignment_done.txt; then
  			echo -e "${magenta}- read alignments and alignment post-processinga already performed ${white}\n"
  		else
  			echo -e "${magenta}- performing read alignments and alignment post-processing ${white}\n"
  			time main &>> log.out
  		fi
  	else
  		echo -e "${magenta}- skipping read alignments and alignment post-processing ${white}\n"
  	fi
  fi
