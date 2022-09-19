
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
if [ -z "$nodes" ]; then
 export nodes=1
fi
if [ -z "$multilocus" ]; then
	export multilocus=false
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
	export p2=$p1
fi
if [ -z "$softclip" ]; then
	export softclip=false
fi
if [ -z "$downsample" ]; then
	export downsample=0
fi
if [ -z "$joint_calling" ]; then
	export joint_calling=false
fi
if [ -z "$keep_gVCF" ]; then
	export keep_gVCF=false
fi
if [ -z "$variant_intervals" ]; then
	export variant_intervals=false
fi



main () {
	cd $projdir
	cd samples
	export nfiles=$(ls -1 -p | grep -v R2.f | grep -v / |  wc -l)
	export totalk=$(awk '/^MemTotal:/{print $2}' /proc/meminfo)
	export loopthreads=2
	if [[ "$threads" -gt 1 ]]; then
	  export N=$((threads/2))
	  export ram1=$(($totalk/$N))
	else
	  export N=1 && export loopthreads=threads
    export ram1=$totalk
	fi
	export ram1=$((ram1/1000000))
	export Xmx1=-Xmx${ram1}G
	export ram2=$(echo "$totalk*0.0000008" | bc)
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


echo -e "${blue}\n############################################################################## ${yellow}\n- Index Reference Genome \n${blue}##############################################################################${white}\n"
main () {
cd $projdir
cd refgenomes
for i in $(ls *.gz 2> /dev/null); do
	gunzip $i >/dev/null 2>&1
done

if ls ./*.ngm 1> /dev/null 2>&1; then
	:
else
	checknfastafiles=$(ls *.f* | grep -v .fai | grep -v .ngm | grep -v _original.fasta | wc -l)
	if [[ $checknfastafiles -gt 1 ]]; then
		echo -e "${magenta}- expecting only 1 fasta file for reference genome ${white}\n"
		echo -e "${magenta}- GBSapp will quit in 5 seconds ${white}\n"
		sleep 5 && exit 1
	fi
	if [ -z "$ref1" ]; then
		for ref in $(ls *.f* 2> /dev/null); do
			ref1=${ref%.fa*}.fasta
		done
	fi
	export ncontigscaffold=$(grep '>' $ref1 | wc -l)
	if [[ $ncontigscaffold -gt 3000 ]]; then
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
			Nstitch=$(printf "A%.0s" $(seq 100))
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
			echo ">""${filename%.txt}" >> $ref1
			cat "$filename" >> $ref1
		done
		rm Chr*
	fi


	export ncontigscaffold=$(grep '>' $ref1 | wc -l)
	if [[ $ncontigscaffold -gt 3000 ]]; then
		mkdir split
		awk '/^>/{close("file.fasta"f);f++}{print $0 > "./split/file"f}' "${ref1}" & PID=$!
		wait $PID
		cd split
		export checksplit=$(wc -c file* | head -n -1 | awk '($1 > 500000000 )' | wc -l)
		if [[ "$checksplit" -gt 0 ]]; then
			for i in $(ls file*); do
				name=$(awk 'NR==1{print $1}' "$i" | awk '{gsub(/>/,""); print}')
				awk 'NR>1{print $0}' $i > ${name}.fasta
				count=$(wc -c ${name}.fasta | awk '{print $1}')
				tr -d '\n' < "${name}".fasta | fold -w 100 |\
				split -d -l 5000000  & PID=$!
			  wait $PID
				for outfile in x*; do
					awk -v name="$name" -v outfile="$outfile" -i inplace 'BEGINFILE{print ">"name"_"outfile}{print}' $outfile
					mv $outfile ${name}_${outfile}.txt
				done
			done
			cd ../
			for i in $(ls *.f* 2> /dev/null); do
				mv $i ./old_"${i%.f*}_fasta.txt"
				cat ./split/*Chr*.txt > $i
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
	if [ -z "$ref1" ]; then
		for ref in $(ls *.f* 2> /dev/null); do
			ref1=${ref%%.f*}.fasta
		done
	fi
else
	echo -e "${magenta}- indexing single reference subgenome ${white}\n"
	if [ -z "$ref1" ]; then
		for ref in $(ls *.f* 2> /dev/null); do
			ref1=${ref%%.f*}.fasta
		done
	fi
	awk '{ sub("\r$",""); print}' $ref1 | awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' > ref.txt
	n=">${ref1%.f*}_"
	awk '{ sub("\r$",""); print}' ref.txt | awk -v n="$n" '{gsub(n,">"); print}' | awk -v n="$n" '{gsub(/>/,n); print}' > $ref1
	rm ref.txt
	$samtools faidx $ref1
	$java -jar $picard CreateSequenceDictionary REFERENCE= $ref1 OUTPUT=${ref1%.f*}.dict
	$ngm -r $ref1

fi

declare -a arr=("${ref1%.f*}.dict" "${ref1}" "${ref1}.fai")
for file in "${arr[@]}"; do
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
	if test ! -f filename_reformatted.txt; then
		if [ -d "se" ]; then
			:
		else
			mkdir se
		fi
		if [ -d "pe" ]; then
			:
		else
			mkdir pe
		fi

		cd ${projdir}/samples/se
		if [ -z "$(ls -A ../pe)" ]; then
			if [ -z "$(ls -A ../se)" ]; then
				cd ../
				for i in $(ls *.f* 2> /dev/null | grep -v R2.f); do
					if [[ "$i" == *.R1* ]]; then
						mv $i ${i/.R1/}
					elif [[ "$i" == *_R1* ]]; then
						mv $i ${i/_R1/}
					else
						:
					fi
				done
			fi
		fi

		cd ${projdir}/samples/se
		if [ -z "$(ls -A ../pe)" ]; then
			if [ "$(ls -A ../se)" ]; then
				echo -e "${magenta}- only single-end reads available in se-folder ${white}\n"
				for i in $(ls *.f* 2> /dev/null ); do
					if [[ "$i" == *.R1* ]]; then
						mv $i ../${i/.R1/}
					elif [[ "$i" == *_R1* ]]; then
						mv $i ../${i/_R1/}
					else
						mv $i ../$i
					fi
				done
			fi
		fi

		cd ${projdir}/samples/pe
		if [ -z "$(ls -A ../se)" ]; then
			if [ "$(ls -A ../pe)" ]; then
				echo -e "${magenta}- only paired-end reads available in pe-folder ${white}\n"
				for i in $(ls *R1.f* 2> /dev/null); do
					mv ${i%R1.f*}R2.f* ../
					if [[ "$i" == *.R1* ]]; then
						mv $i ../${i/.R1/}
					elif [[ "$i" == *_R1* ]]; then
						mv $i ../${i/_R1/}
					else
						echo -e "${magenta}- check paired-end filenames for proper filename format, i.e. .R1 or _R1 and .R2 or _R2  ${white}\n"
						echo -e "${magenta}- Do you want to continue running GBSapp? ${white}\n"
						read -p "- y(YES) or n(NO) " -n 1 -r
						if [[ ! $REPLY =~ ^[Yy]$ ]]; then
							printf '\n'
							exit 1
						fi
					fi
				done
			fi
		fi

		cd ${projdir}/samples/pe
		if [ "$(ls -A ../se)" ]; then
			if [ "$(ls -A ../pe)" ]; then
				for i in $(ls *R1.f* 2> /dev/null); do
					mv ${i%R1.f*}R2.f* ../
					if [[ "$i" == *.R1* ]]; then
						cat $i ../se/${i%.R1.f*}* > ../${i}
						mv ../${i} ../${i/.R1/}
						rm ../pe/$i ../se/${i%.R1.f*}*
					elif [[ "$i" == *_R1* ]]; then
						cat $i ../se/${i%_R1.f*}* > ../${i}
						mv ../${i} ../${i/_R1/}
						rm ../pe/$i ../se/${i%_R1.f*}*
					else
						echo -e "${magenta}- check paired-end filenames for proper filename format, i.e. .R1 or _R1 and .R2 or _R2 ${white}\n"
						echo -e "${magenta}- Do you want to continue running GBSapp? ${white}\n"
						read -p "- y(YES) or n(NO) " -n 1 -r
						if [[ ! $REPLY =~ ^[Yy]$ ]]; then
							printf '\n'
							exit 1
						fi
					fi
				done
			fi
		fi

		cd ..
		sampno=$(ls -1 | wc -l)
		if [[ "$sampno" == "0" ]]; then
			echo -e "${magenta}- \n- samples folder is empty, exiting pipeline ${white}\n"
			exit 1
		fi
		find . -type d -empty -delete
		echo filename_reformatted > filename_reformatted.txt
	fi

	if test ! -f flushed_reads.txt; then
		if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
			:> length_distribution.txt
			for i in $(ls -S *.f* | grep -v _uniq.fasta | grep -v _uniq_R1.fasta | grep -v _uniq_R2.fasta | grep -v _uniq.hold.fasta | grep -v _uniq_R1.hold.fasta | grep -v _uniq_R2.hold.fasta | grep -v fq.gz); do
				if [[ $(file $i 2> /dev/null) =~ gzip ]]; then
					fa_fq=$(zcat $projdir/samples/$i 2> /dev/null | head -n1 | cut -c1-1)
				else
					fa_fq=$(cat $projdir/samples/$i | head -n1 | cut -c1-1)
				fi

				if [[ $(file $i 2> /dev/null) =~ gzip ]]; then
					if [[ "${fa_fq}" == "@" ]]; then
						awk 'NR%2==0' <(zcat $i) | awk 'NR%2==1' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=1000; i++){x=int(rand()*NR) + 1; print a[x];}}' >> length_distribution.txt
					fi
					if [[ "${fa_fq}" == ">" ]]; then
						awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' <(zcat $i) | \
						awk 'NR%2==0' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=1000; i++){x=int(rand()*NR) + 1; print a[x];}}' >> length_distribution.txt
					fi
				else
					if [[ "${fa_fq}" == "@" ]]; then
						awk 'NR%2==0' $i | awk 'NR%2==1' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=1000; i++){x=int(rand()*NR) + 1; print a[x];}}' >> length_distribution.txt
					fi
					if [[ "${fa_fq}" == ">" ]]; then
						awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' $i | \
						awk 'NR%2==0' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=1000; i++){x=int(rand()*NR) + 1; print a[x];}}' >> length_distribution.txt
					fi
				fi
			done

			awk '{print length($0)}' length_distribution.txt | sort -n > tmp.txt; mv tmp.txt length_distribution.txt
			export max_seqread_len=$(awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}' length_distribution.txt)
			rm length_distribution.txt

			for i in $(ls -S *.f* | grep -v _uniq.fasta | grep -v _uniq_R1.fasta | grep -v _uniq_R2.fasta | grep -v _uniq.hold.fasta | grep -v _uniq_R1.hold.fasta | grep -v _uniq_R2.hold.fasta | grep -v fq.gz); do (
				if [[ $(file $i 2> /dev/null) =~ gzip ]]; then
					fa_fq=$(zcat $projdir/samples/$i 2> /dev/null | head -n1 | cut -c1-1)
				else
					fa_fq=$(cat $projdir/samples/$i | head -n1 | cut -c1-1)
				fi

				if [[ $(file $i 2> /dev/null) =~ gzip ]]; then
					if [[ "${fa_fq}" == "@" ]]; then
						awk 'NR%2==0' <(zcat $i) | awk 'NR%2==1' | awk -v max=$max_seqread_len '{print substr($0,1,max)}' | awk '{print "@"NR"\t"$1"\t"$1}' | \
						awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | awk '{print $1"\n"$2"\n+\n"$3}' | $gzip > ${i%.f*}_tmp.fa.gz  && mv ${i%.f*}_tmp.fa.gz $i
					fi
					if [[ "${fa_fq}" == ">" ]]; then
						grep -v '^>' <(zcat $i) | awk -v max=$max_seqread_len '{print substr($0,1,max)}' | awk '{print ">"NR"\n"$1}' | $gzip > ${i%.f*}_tmp.fa.gz && mv ${i%.f*}_tmp.fa.gz $i
					fi
				else
					if [[ "${fa_fq}" == "@" ]]; then
						awk 'NR%2==0' $i | awk 'NR%2==1' | awk -v max=$max_seqread_len '{print substr($0,1,max)}' | awk '{print "@"NR"\t"$1"\t"$1}' | \
						awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | awk '{print $1"\n"$2"\n+\n"$3}' > ${i%.f*}_tmp.fa.gz && mv ${i%.f*}_tmp.fa.gz $i
					fi
					if [[ "${fa_fq}" == ">" ]]; then
						grep -v '^>' $i | awk -v max=$max_seqread_len '{print substr($0,1,max)}' | awk '{print ">"NR"\n"$1}' > ${i%.f*}_tmp.fa.gz && mv ${i%.f*}_tmp.fa.gz $i
					fi
				fi ) &
				if [[ $(jobs -r -p | wc -l) -ge $prepN ]]; then
					wait
				fi
			done
		fi
		find . -type d -empty -delete
		printf "Improvement in flushed reads already implemented""\n" > flushed_reads.txt
	fi
	if [[ "$lib_type" =~ "WGS" ]] || [[ "$lib_type" =~ "wgs" ]]; then
		printf "Improvement in flushed reads not required for shotgun WGS data""\n" > flushed_reads.txt
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
	mkdir -p preprocess/tmp
  mkdir -p preprocess/alignment
	mkdir -p snpcall/tmp
	mkdir -p alignment_summaries
	mkdir -p ./alignment_summaries/copy_number

	if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
		for i in samples_list_node_*.txt; do
			:> ${i%.txt}_hold.txt
			while read line; do
				ls -l ./samples/$line | awk '{print $5"\t"$9}' >> ${i%.txt}_hold.txt
			done < $i
			sort -nr -k1 ${i%.txt}_hold.txt | awk '{gsub(/.\/samples\//,""); print $2}' > $i
			rm ${i%.txt}_hold.txt
		done
	fi
}
if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
	if [[ ! -f "${projdir}/organize_files_done.txt" ]]; then time main &>> ${projdir}/log.out; fi
fi



######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${yellow}\n- Performing Read Alignments & Alignment Post-Processing\n${blue}##############################################################################${white}\n"

cd $projdir
if [[ $nodes -gt 1 ]] && [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
	rm -rf /tmp/${samples_list%.txt} 2> /dev/null &&
	mkdir -p /tmp/${samples_list%.txt}/refgenomes /tmp/${samples_list%.txt}/samples
  mkdir -p /tmp/${samples_list%.txt}/preprocess/tmp /tmp/${samples_list%.txt}/preprocess/alignment /tmp/${samples_list%.txt}/snpcall/tmp
	cp -r ${projdir}/refgenomes/* /tmp/${samples_list%.txt}/refgenomes/
	for i in $(cat ${projdir}/${samples_list} ); do
		cp ${i%.f*}_uniq_R1.fasta.gz /tmp/${samples_list%.txt}/samples/ 2> /dev/null &&
		cp ${projdir}/preprocess/alignment/${i%.f*}_redun.sam.gz /tmp/${samples_list%.txt}/preprocess/alignment/ 2> /dev/null &&
		cp ${projdir}/preprocess/${i%.f*}_*_precall.bam* /tmp/${samples_list%.txt}/preprocess/ 2> /dev/null &&
		wait
	done
fi

main () {

	cd $projdir
	cd samples

	if [[ "$samples_list" == "samples_list_node_1.txt" ]] && test ! -f ${projdir}/organize_files_done.txt; then
		for i in $( cat ${projdir}/samples_list_node_* ); do
			if test ! -f ${projdir}/preprocess/${i%.f*}_redun.sam && test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai; then
				if test ! -f ${i%.f*}_uniq_R1.fasta.gz; then
					if [[ $(file $i | awk -F' ' '{print $2}') == gzip ]]; then
						zcat $i 2> /dev/null | awk 'NR%2==0' | awk 'NR%2' | gzip > ${i%.f*}_uniq.txt.gz 2> /dev/null &&
						wait
						if test -f ${i%.f*}_R2*; then
							zcat ${i%.f*}_R2* 2> /dev/null | awk 'NR%2==0' | awk 'NR%2' | gzip > ${i%.f*}_R2_uniq.txt.gz 2> /dev/null &&
							wait
						fi
						wait
						if test -f ${i%.f*}.R2*; then
							zcat ${i%.f*}.R2* 2> /dev/null | awk 'NR%2==0' | awk 'NR%2' | gzip > ${i%.f*}_R2_uniq.txt.gz 2> /dev/null &&
							wait
						fi
						wait
					else
						awk 'NR%2==0' $i | awk 'NR%2' | gzip > ${i%.f*}_uniq.txt.gz 2> /dev/null &&
						wait
						if test -f ${i%.f*}_R2*; then
							awk 'NR%2==0' ${i%.f*}_R2* | awk 'NR%2' | gzip > ${i%.f*}_R2_uniq.txt.gz 2> /dev/null &&
							wait
						fi
						wait
						if test -f ${i%.f*}.R2*; then
							awk 'NR%2==0' ${i%.f*}.R2* | awk 'NR%2' | gzip > ${i%.f*}_R2_uniq.txt.gz 2> /dev/null &&
							wait
						fi
						wait
					fi
					wait

					if test -f "${i%.f*}_R2_uniq.txt.gz"; then
						export LC_ALL=C; paste -d ~ <(zcat ${i%.f*}_uniq.txt.gz 2> /dev/null) <(zcat ${i%.f*}_R2_uniq.txt.gz 2> /dev/null) | expand -t $(( $(wc -L < $i ) + 2 )) | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{$1=$1};1' | \
						awk '{gsub(" /"," "); print}' | awk '{gsub("/\n","\n"); print}' | awk '{gsub("/"," "); print}' | awk '{gsub(" ","\t"); print}' | gzip > ${i%.f*}_rdrefseq.txt.gz 2> /dev/null &&
						repseqs=$(awk '{all[NR] = $0} END{print all[int(NR*0.999 - 0.5)]}' <(zcat ${i%.f*}_rdrefseq.txt.gz | awk '{print $1}' | sort -n))
						awk -v repseqs=$repseqs '$1<=repseqs{print $0}' <(zcat ${i%.f*}_rdrefseq.txt.gz) | gzip > ${i%.f*}_rdrefseq.tmp.gz; mv ${i%.f*}_rdrefseq.tmp.gz ${i%.f*}_rdrefseq.txt.gz
						wait
						awk 'NF==2 {print ">seq"NR"_se-"$1"\t"$2}' <(zcat ${i%.f*}_rdrefseq.txt.gz 2> /dev/null) | gzip > ${i%.f*}_rdrefseq_se.txt.gz 2> /dev/null &&
						wait
						awk 'NF==3 {print ">seq"NR"_pe-"$0}' <(zcat ${i%.f*}_rdrefseq.txt.gz 2> /dev/null) | awk '{print $1"\t"$3}' | gzip > ${i%.f*}_uniq_R2.fasta.gz 2> /dev/null &&
						wait
						awk 'NF==3 {print ">seq"NR"_pe-"$0}' <(zcat ${i%.f*}_rdrefseq.txt.gz 2> /dev/null) | awk '{print $1"\t"$2}' | cat - <(zcat ${i%.f*}_rdrefseq_se.txt.gz 2> /dev/null) | gzip > ${i%.f*}_uniq_R1.hold.fasta.gz 2> /dev/null &&
						wait
						rm ${i%.f*}*.txt* 2> /dev/null &&
						wait
						find . -size 0 -delete  2> /dev/null &&
						wait
						mv ${i%.f*}_uniq_R1.hold.fasta.gz ${i%.f*}_uniq_R1.fasta.gz  2> /dev/null &&
						wait
					else
						touch ${i%.f*}_R2_uniq.txt &&
						wait
						export LC_ALL=C; paste -d ~ <(zcat ${i%.f*}_uniq.txt.gz 2> /dev/null) ${i%.f*}_R2_uniq.txt | expand -t $(( $(wc -L < $i ) + 2 )) | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{$1=$1};1' | \
						awk '{gsub(" /"," "); print}' | awk '{gsub("/\n","\n"); print}' | awk '{gsub("/"," "); print}' | awk '{gsub(" ","\t"); print}' | gzip > ${i%.f*}_rdrefseq.txt.gz 2> /dev/null &&
						repseqs=$(awk '{all[NR] = $0} END{print all[int(NR*0.999 - 0.5)]}' <(zcat ${i%.f*}_rdrefseq.txt.gz | awk '{print $1}' | sort -n))
						awk -v repseqs=$repseqs '$1<=repseqs{print $0}' <(zcat ${i%.f*}_rdrefseq.txt.gz) | gzip > ${i%.f*}_rdrefseq.tmp.gz; mv ${i%.f*}_rdrefseq.tmp.gz ${i%.f*}_rdrefseq.txt.gz
						wait
						awk 'NF==2 {print ">seq"NR"_pe-"$0}' <(zcat ${i%.f*}_rdrefseq.txt.gz 2> /dev/null) | awk '{print $1"\t"$2}' | gzip > ${i%.f*}_uniq_R1.hold.fasta.gz 2> /dev/null &&
						wait
						rm ${i%.f*}*.txt* 2> /dev/null &&
						wait
						find . -size 0 -delete  2> /dev/null &&
						wait
						mv ${i%.f*}_uniq_R1.hold.fasta.gz ${i%.f*}_uniq_R1.fasta.gz  2> /dev/null &&
						wait
					fi
					wait
				fi
			fi
		done

		cd $projdir/samples
		find . -size 0 -delete  2> /dev/null &&
		touch ../report_fq_compress_index.txt
		for i in $( cat ${projdir}/samples_list_node_* ); do
			if [[ "$(zcat ${i%.f*}_uniq_R1.fasta.gz 2> /dev/null | head | wc -l)" -eq 0 ]] || [[ -z "${i%.f*}_uniq_R1.fasta.gz" ]]; then
				echo ${i%.f*}_uniq_R1.fasta.gz >> ../report_fq_compress_index.txt
			fi
		done
		if [[ -s ../report_fq_compress_index.txt ]]; then
			echo -e "${magenta}- $(wc -l ../report_fq_compress_index.txt | awk '{print $1}') fastq files not properly processed ${white}\n"
			END=10
			while [[ $END -gt 0 ]]; do
				for i in $( cat ${projdir}/${samples_list} ); do
					if [[ "$lib_type" =~ "RRS" || "$lib_type" =~ "rrs" ]] && test ! -f ${projdir}/compress_done.txt && test ! -f ${projdir}/organize_files_done.txt && test ! -f ${projdir}/preprocess/${i%.f*}_redun.sam && test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai; then
						if test ! -f ${i%.f*}_uniq_R1.fasta.gz; then
							if [[ $(file $i | awk -F' ' '{print $2}') == gzip ]]; then
								zcat $i 2> /dev/null | awk 'NR%2==0' | awk 'NR%2' | $gzip > ${i%.f*}_uniq.txt.gz 2> /dev/null &&
								wait
								if test -f ${i%.f*}_R2*; then
									zcat ${i%.f*}_R2* 2> /dev/null | awk 'NR%2==0' | awk 'NR%2' | $gzip > ${i%.f*}_R2_uniq.txt.gz 2> /dev/null &&
									wait
								fi
								wait
								if test -f ${i%.f*}.R2*; then
									zcat ${i%.f*}.R2* 2> /dev/null | awk 'NR%2==0' | awk 'NR%2' | $gzip > ${i%.f*}_R2_uniq.txt.gz 2> /dev/null &&
									wait
								fi
								wait
							else
								awk 'NR%2==0' $i | awk 'NR%2' | $gzip > ${i%.f*}_uniq.txt.gz 2> /dev/null &&
								wait
								if test -f ${i%.f*}_R2*; then
									awk 'NR%2==0' ${i%.f*}_R2* | awk 'NR%2' | $gzip > ${i%.f*}_R2_uniq.txt.gz 2> /dev/null &&
									wait
								fi
								wait
								if test -f ${i%.f*}.R2*; then
									awk 'NR%2==0' ${i%.f*}.R2* | awk 'NR%2' | $gzip > ${i%.f*}_R2_uniq.txt.gz 2> /dev/null &&
									wait
								fi
								wait
							fi
							wait

							if test -f "${i%.f*}_R2_uniq.txt.gz"; then
								export LC_ALL=C; paste -d ~ <(zcat ${i%.f*}_uniq.txt.gz 2> /dev/null) <(zcat ${i%.f*}_R2_uniq.txt.gz 2> /dev/null) | expand -t $(( $(wc -L < $i ) + 2 )) | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{$1=$1};1' | \
								awk '{gsub(" /"," "); print}' | awk '{gsub("/\n","\n"); print}' | awk '{gsub("/"," "); print}' | awk '{gsub(" ","\t"); print}' | $gzip > ${i%.f*}_rdrefseq.txt.gz 2> /dev/null &&
								repseqs=$(awk '{all[NR] = $0} END{print all[int(NR*0.999 - 0.5)]}' <(zcat ${i%.f*}_rdrefseq.txt.gz | awk '{print $1}' | sort -n))
								awk -v repseqs=$repseqs '$1<=repseqs{print $0}' <(zcat ${i%.f*}_rdrefseq.txt.gz) | gzip > ${i%.f*}_rdrefseq.tmp.gz; mv ${i%.f*}_rdrefseq.tmp.gz ${i%.f*}_rdrefseq.txt.gz
								wait
								awk 'NF==2 {print ">seq"NR"_se-"$1"\t"$2}' <(zcat ${i%.f*}_rdrefseq.txt.gz 2> /dev/null) | $gzip > ${i%.f*}_rdrefseq_se.txt.gz 2> /dev/null &&
								wait
								awk 'NF==3 {print ">seq"NR"_pe-"$0}' <(zcat ${i%.f*}_rdrefseq.txt.gz 2> /dev/null) | awk '{print $1"\t"$3}' | $gzip > ${i%.f*}_uniq_R2.fasta.gz 2> /dev/null &&
								wait
								awk 'NF==3 {print ">seq"NR"_pe-"$0}' <(zcat ${i%.f*}_rdrefseq.txt.gz 2> /dev/null) | awk '{print $1"\t"$2}' | cat - <(zcat ${i%.f*}_rdrefseq_se.txt.gz 2> /dev/null) | $gzip > ${i%.f*}_uniq_R1.hold.fasta.gz 2> /dev/null &&
								wait
								rm ${i%.f*}*.txt* 2> /dev/null &&
								wait
								find . -size 0 -delete  2> /dev/null &&
								wait
								mv ${i%.f*}_uniq_R1.hold.fasta.gz ${i%.f*}_uniq_R1.fasta.gz  2> /dev/null &&
								wait
							else
								touch ${i%.f*}_R2_uniq.txt &&
								wait
								export LC_ALL=C; paste -d ~ <(zcat ${i%.f*}_uniq.txt.gz 2> /dev/null) ${i%.f*}_R2_uniq.txt | expand -t $(( $(wc -L < $i ) + 2 )) | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{$1=$1};1' | \
								awk '{gsub(" /"," "); print}' | awk '{gsub("/\n","\n"); print}' | awk '{gsub("/"," "); print}' | awk '{gsub(" ","\t"); print}' | $gzip > ${i%.f*}_rdrefseq.txt.gz 2> /dev/null &&
								repseqs=$(awk '{all[NR] = $0} END{print all[int(NR*0.999 - 0.5)]}' <(zcat ${i%.f*}_rdrefseq.txt.gz | awk '{print $1}' | sort -n))
								awk -v repseqs=$repseqs '$1<=repseqs{print $0}' <(zcat ${i%.f*}_rdrefseq.txt.gz) | gzip > ${i%.f*}_rdrefseq.tmp.gz; mv ${i%.f*}_rdrefseq.tmp.gz ${i%.f*}_rdrefseq.txt.gz
								wait
								awk 'NF==2 {print ">seq"NR"_pe-"$0}' <(zcat ${i%.f*}_rdrefseq.txt.gz 2> /dev/null) | awk '{print $1"\t"$2}' | $gzip > ${i%.f*}_uniq_R1.hold.fasta.gz 2> /dev/null &&
								wait
								rm ${i%.f*}*.txt* 2> /dev/null &&
								wait
								find . -size 0 -delete  2> /dev/null &&
								wait
								mv ${i%.f*}_uniq_R1.hold.fasta.gz ${i%.f*}_uniq_R1.fasta.gz  2> /dev/null &&
								wait
							fi
							wait
						fi
					fi
				done
				END=$(($END-1))
				: > ../report_fq_compress_index.txt
				for j in $( cat ${projdir}/${samples_list} ); do
					if [[ "$(zcat ${j%.f*}_uniq_R1.fasta.gz 2> /dev/null | head -n1 | wc -l)" -eq 0 ]] || [[ -z "${j%.f*}_uniq_R1.fasta.gz" ]]; then
						echo ${j%.f*}_uniq_R1.fasta.gz >> ../report_fq_compress_index.txt
					fi
				done
			done
		fi
		find ../report_fq_compress_index.txt -size 0 -delete  2> /dev/null &&
		wait
	fi
	wait
	touch ${projdir}/organize_files_done.txt


  if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
    touch ${projdir}/compress_done.txt
  fi

  cd ${projdir}/samples

  if [[ $nodes -gt 1 ]]; then
    if [[ "$samples_list" != "samples_list_node_1.txt" ]]; then
      rm -rf /tmp/${samples_list%.txt} 2> /dev/null
    fi
    mkdir -p /tmp/${samples_list%.txt}/refgenomes /tmp/${samples_list%.txt}/samples
    mkdir -p /tmp/${samples_list%.txt}/preprocess/tmp /tmp/${samples_list%.txt}/preprocess/alignment /tmp/${samples_list%.txt}/snpcall/tmp
    :> touch ${projdir}/queue_move_${samples_list%.txt}.txt
    queue_move=$(ls ${projdir}/queue_move_samples_list_node_* | wc -l)
    while [[ "$queue_move" -gt 1 ]]; do
      rm ${projdir}/queue_move_${samples_list%.txt}.txt; sleep $[ ( $RANDOM % 120 )  + 30 ]s
      :> ${projdir}/queue_move_${samples_list%.txt}.txt
      queue_move=$(ls ${projdir}/queue_move_samples_list_node_* | wc -l)
    done
    cp -rn ${projdir}/refgenomes/* /tmp/${samples_list%.txt}/refgenomes/ &&
    for i in $(cat ${projdir}/${samples_list} ); do
      cp -rn ${i%.f*}_uniq_R*.fasta.gz /tmp/${samples_list%.txt}/samples/ 2> /dev/null &&
      cp -rn ${projdir}/preprocess/alignment/${i%.f*}_redun.sam.gz /tmp/${samples_list%.txt}/preprocess/alignment/ 2> /dev/null &&
      cp -rn ${projdir}/preprocess/${i%.f*}_*_precall.bam* /tmp/${samples_list%.txt}/preprocess/ 2> /dev/null &&
      wait
    done
    wait
    rm ${projdir}/queue_move_${samples_list%.txt}.txt
  fi

  cd ${projdir}
  if [[ $nodes -eq 1 ]]; then cd ${projdir}/samples/ ; fi
  if [[ $nodes -gt 1 ]] && test -f ${projdir}/GBSapp_run_node_1.sh; then cd /tmp/${samples_list%.txt}/samples/ ; fi
  if test ! -f ${projdir}/precall_done.txt && test ! -f ${projdir}/alignment_done; then
    for alignfq in $( cat ${projdir}/${samples_list} ); do
      if test ! -f ../preprocess/alignment/${alignfq%.f*}_redun.sam.gz; then
	        zcat ${alignfq%.f*}_uniq_*.fasta.gz | awk '{gsub(/>/,"@");}1'| awk '{print $1"\t"$2"\t"$2}' | \
					awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | awk '{print $1"\n"$2"\n+\n"$3}'| gzip > ${alignfq%.f*}_uniq.fastq.gz &&
	        $ngm -r ../refgenomes/$ref1 --qry ${alignfq%.f*}_uniq.fastq.gz -o ../preprocess/alignment/${alignfq%.f*}_redun.sam -t $threads --affine --min-identity 0 --topn 12 --strata 12 &&
          $gzip ../preprocess/alignment/${alignfq%.f*}_redun.sam &&
          cp -rn ../preprocess/alignment/${alignfq%.f*}_redun.sam.gz ${projdir}/preprocess/alignment/ &&
          rm ${alignfq%.f*}_uniq.fastq.gz &&
          wait
      fi
    done
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


  for i in $(cat ${projdir}/${samples_list} ); do (
    printf '\n###---'${i%.f*}'---###\n' > ${projdir}/alignment_summaries/${i%.f*}_summ.txt && \
    zcat ./alignment/${i%.f*}_redun.sam.gz | grep -v '^@PG' | tr ' ' '\t' | $samtools flagstat - >> ${projdir}/alignment_summaries/${i%.f*}_summ.txt && \
    printf '########################################################################################################\n\n' >> ${projdir}/alignment_summaries/${i%.f*}_summ.txt && \
    printf 'copy\tFrequency\tPercentage\n' > ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_Unique_Read_histogram.txt && \
    grep -v '^@' <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | awk -F' ' '{print $1}' | awk '{gsub(/_/,"\t");gsub(/\//,"\t");gsub(/pe-/,"");gsub(/se-/,""); print $2}' | \
    awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{print $2"\t"$1}' > ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number.txt  && \
    awk 'NR==FNR{sum+= $2; next;} {printf("%s\t%s\t%3.3f%%\t%3.0f\n",$1,$2,100*$2/sum,100*$2/sum)}' ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number.txt ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number.txt > ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot.txt && \
    unset IFS; printf "%s\t%s\t%s\t%*s\n" $(sed 's/$/ |/' ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot.txt) | tr ' ' '|' | sort -k2,2 -nr | awk '{gsub(/se-/,""); gsub(/pe-/,""); print}' >> ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_Unique_Read_histogram.txt &&
    rm ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number.txt ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot.txt


    if test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai; then
      # grep -v '^@' <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | awk '($3 != "\*")' 2> /dev/null  | awk '($6 != "\*")' 2> /dev/null  | awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | \
      # cat <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null | grep '^@') - > ${i%.f*}_del.sam
      # wait
      # echo "nloci~mapQ~CHROM~POS" > ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}.txt
      # grep -v '^@' ${i%.f*}_del.sam | awk -F '\t' '{print $1"\t"$3"\t"$4"\t"$5}' | awk '{gsub(/ /,"\t"); print}' | awk '{print $1"~"$5"~"$3"~"$4}' | grep -v '\*' >> ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}.txt &&
      # awk '!visited[$0]++' ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}.txt > ${projdir}/alignment_summaries/temp_${i%.f*}.txt &&
      # mv ${projdir}/alignment_summaries/temp_${i%.f*}.txt ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}.txt
      # rm ${i%.f*}_del.sam
      # wait

      if [[ "$multilocus" == true ]]; then
        awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) > ${i%.f*}_heading.sam &&
        grep -v '^@' <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | awk '($3 != "\*")' 2> /dev/null | awk '($6 != "\*")' 2> /dev/null  |  awk '{gsub(/_se-/,"_se-\t",$1); gsub(/_pe-/,"_pe-\t",$1)}1' | \
        awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 >= min) {print $0}}' > ${i%.f*}_uniq.sam
        wait
      fi
      if [[ "$multilocus" == false ]]; then
        awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) > ${i%.f*}_heading.sam &&
        grep -v '^@' <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | awk '($3 != "\*")' 2> /dev/null  | awk '($6 != "\*")' 2> /dev/null | awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | \
        tr " " "\t" | tr '\r' '\n' | awk '$1==1{print $0}' | awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | awk '{gsub(/_se-/,"_se-\t",$1); gsub(/_pe-/,"_pe-\t",$1)}1' | \
        awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" | awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 60) {print $0}}' > ${i%.f*}_uniq.sam
        wait
      fi

      for k in $(awk '{A[$1]++}END{for(i in A)print i}' ${i%.f*}_uniq.sam); do (
        awk -v n="${k}" '$1 == n{print $0}' ${i%.f*}_uniq.sam | awk -v n="$k" '{for(i=0;i<n;i++) print}' > ${i%.f*}_exp_${k}.sam ) &
        if [[ $(jobs -r -p | wc -l) -ge $prepthreads ]]; then
        wait
        fi
      done
      wait
      cat ${i%.f*}_exp_*.sam > ${i%.f*}_exp.sam && rm ${i%.f*}*_exp_*.sam &&
      awk '{print "seq"NR"_"$0}' ${i%.f*}_exp.sam | tr -s ' ' | \
      awk '!($3 ~ "\*")' 2> /dev/null | awk '!($6 ~ "\*")' 2> /dev/null | cat ${i%.f*}_heading.sam - | tr ' ' '\t' > ${i%.f*}_${ref1%.f*}.sam &&
      rm ${i%.f*}_exp.sam ${i%.f*}_uniq.sam ${i%.f*}_heading.sam


      j="${i%.f*}_${ref1%.f*}.sam"
      $java $Xmxp -XX:ParallelGCThreads=$prepthreads -Djava.io.tmpdir=`pwd`/tmp -jar $picard SortSam I=$j O=${j%.sam*}.bam  SORT_ORDER=coordinate  VALIDATION_STRINGENCY=LENIENT TMP_DIR=`pwd`/tmp && \
      $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard BuildBamIndex INPUT=${j%.sam*}.bam VALIDATION_STRINGENCY=LENIENT && \
      $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard AddOrReplaceReadGroups I=${j%.sam*}.bam O=${j%.sam*}_precall.bam RGLB=${i%.f*} RGPL=illumina RGPU=run RGSM=${i%.f*} VALIDATION_STRINGENCY=LENIENT && \
      $samtools index ${j%.sam*}_precall.bam &&
      rm $j ${j%.sam*}.bam ${j%.sam*}.bai &&
      if [[ $nodes -gt 1 ]]; then cp /tmp/${samples_list%.txt}/preprocess/${j%.sam*}_precall.bam* ${projdir}/preprocess/; fi
      wait
    fi ) &
    if [[ $(jobs -r -p | wc -l) -ge $prepN ]]; then
      wait
    fi
  done
  wait && touch ${projdir}/precall_done_${samples_list}
  wait
  ls * | grep -v precall | xargs rm 2> /dev/null
  wait


	cd ${projdir}/samples
	if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
		precall=$(ls ${projdir}/precall_done_samples_list_node_* | wc -l)
		while [[ "$precall" -lt $nodes ]]; do sleep 300; precall=$(ls ${projdir}/precall_done_samples_list_node_* | wc -l); done
		if [[ $precall == $nodes ]] && test ! -f ${projdir}/alignment_summaries/refgenome_paralogs.txt; then
			cd ${projdir}/alignment_summaries
			cat *_summ.txt > alignment_summaries_unique_reads.txt; rm -r *_summ.txt &&
			# Total number of reads per samples
			awk '/###---/ || /QC-passed/{print}' alignment_summaries_unique_reads.txt | cut -d\+ -f1 | tr -d '\n' | \
			awk  'gsub(/---###/, "\t", $0)' | awk  'gsub(/###---/, "", $0)' | tr ' ' '\n' > total_unique_reads.txt &&
			# Total number of mapped reads per samples
			cat alignment_summaries_unique_reads.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /0_mapped/{print}' |\
			tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
			awk 'gsub("\\+0mapped", "\t", $0)' | cut -d\: -f1 > total_unique_reads_mapped.txt &&
			# Total number of mapped paired reads per samples
			cat alignment_summaries_unique_reads.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /properly_paired/{print}' |\
			tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
			awk 'gsub("\\+0properlypaired", "\t", $0)' | cut -d\: -f1 > total_unique_reads_paired.txt &&
			echo -e "Samples\tTotal\tMapped\tPerc_Mapped\tPE_Mapped\t%_PE_Mapped" > summary_precall.txt &&
			awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' total_unique_reads_mapped.txt  total_unique_reads.txt  | \
			awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' total_unique_reads_paired.txt - | awk '{gsub(/ /,"\t"); print $0}' | \
			cat summary_precall.txt - > Tabulated_Alignment_Unique_Read_Summaries.txt &&
			awk '{gsub(/\t/,","); print $0}' Tabulated_Alignment_Unique_Read_Summaries.txt > Tabulated_Alignment_Unique_Read_Summaries.csv &&
			rm total_unique_* summary_precall.txt &> /dev/null &&
			rm ${projdir}/samples/metrics.txt ${projdir}/preprocess/metrics.txt &> /dev/null &&
      wait
			# cd $projdir/alignment_summaries
			# touch refgenome_paralogs.txt
			# for par in refgenome_paralogs_*.txt; do
			# 	cat refgenome_paralogs.txt $par | awk '!visited[$0]++' > temp_par.txt &&
			# 	mv temp_par.txt refgenome_paralogs.txt &&
			# 	wait
			# done
			# wait
			# rm refgenome_paralogs_*.txt
			# awk '{gsub(/~/,"\t"); print $0}' refgenome_paralogs.txt > temp.txt &&
			# awk '{print $3"\t"$4"\t"$1}' temp.txt | awk '!/CHROM|nloci|POS|mapQ/' | cat <(printf "CHROM\tPOS\tnloci\n") - > refgenome_paralogs.txt &&
			# rm temp.txt
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


touch ${projdir}/compress_done.txt


if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
  mkdir -p ${projdir}/preprocess/processed
	mv ${projdir}/preprocess/processed/*_precall.bam* ${projdir}/preprocess/ 2> /dev/null
	touch ${projdir}/call0
	sleep 5
fi
wait
while [[ ! -f "${projdir}/call0" ]]; do sleep 5; done
wait
if [[ $nodes -gt 1 ]]; then
  mkdir -p /tmp/${samples_list%.txt}/refgenomes /tmp/${samples_list%.txt}/samples
  mkdir -p /tmp/${samples_list%.txt}/preprocess/tmp /tmp/${samples_list%.txt}/preprocess/alignment /tmp/${samples_list%.txt}/snpcall/tmp
  cp -rn ${projdir}/refgenomes/* /tmp/${samples_list%.txt}/refgenomes/
  wait
  for i in $(cat ${projdir}/${samples_list} ); do
    cp -rn ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam* /tmp/${samples_list%.txt}/preprocess/
    wait
  done
fi


cd $projdir
echo -e "${magenta}- performing SNP calling ${white}\n"
cd $projdir
cd preprocess

if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
	if [[ "$joint_calling" == true ]]; then
		j=-I; input=""; k=""
		for i in $(ls *_precall.bam 2> /dev/null); do
			k="${j} ${i}"; input="${input} ${k}"
		done
		Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/${ref1%.*}.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -n | awk '{print $1}' | awk -v pat=${ref1%.f*} '$0 ~ pat')
		if [[ ! -z "$Exclude_Chromosome" ]]; then
			for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
				Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
			done
		fi
		if [[ -z "$Get_Chromosome" ]]; then
      if [[ -z "$interval_list" ]]; then
  			for selchr in $Get2_Chromosome; do (
  				if [[ -z "$(ls ${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf* 2> /dev/null)" ]]; then
  					$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/$ref1 -L ${selchr} ${input} -ploidy $ploidy -O ${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf.gz --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality $minmapq --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) &&
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
          if [[ -z "$(ls ${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf* 2> /dev/null)" ]]; then
            $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/$ref1 -L ${projdir}/variant_intervals_${selchr}.list ${input} -ploidy $ploidy -O ${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf.gz --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality $minmapq --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) &&
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
			grep -h '^#' ${pop}_${ploidy}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt
			cat ${pop}_${ploidy}x_*_raw.vcf | awk '!/^#/' > all.vcf
			cat vcf_header.txt all.vcf > ${pop}_${ploidy}x_raw.vcf
			rm vcf_header.txt all.vcf *.vcf.gz.tbi ${pop}_${ploidy}x_*_raw.vcf
		else
			if [[ -z "$(ls ${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf* 2> /dev/null)" ]]; then
				echo $Get_Chromosome | tr ',' '\n' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ${projdir}/refgenomes/${ref1%.fasta}.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2":1-"$3}' > ${projdir}/refgenomes/${ref1%.fasta}.list
				$GATK --java-options "$Xmx2 -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$threads" HaplotypeCaller -R ${projdir}/refgenomes/$ref1 -L ${projdir}/refgenomes/${ref1%.fasta}.list ${input} -ploidy $ploidy -O ${projdir}/snpcall/${pop}_${ploidy}x_raw.vcf.gz --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality $minmapq --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) &&
				cd ../snpcall
				gunzip ${pop}_${ploidy}x_raw.vcf.gz &&
				wait
			fi
		fi
	fi
fi

if [[ "$joint_calling" == false ]]; then

  if [[ $nodes -eq 1 ]]; then cd ${projdir}/preprocess/; fi
  if [[ $nodes -gt 1 ]] && test -f ${projdir}/GBSapp_run_node_1.sh; then cd /tmp/${samples_list%.txt}/preprocess/; fi

	for i in $(cat ${projdir}/${samples_list} ); do (
		if [[ -z "$(ls ${projdir}/snpcall/${pop}_${ploidy}x_raw.vcf* 2> /dev/null)" ]]; then
			if test ! -f "${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf"; then
					if [[ -z "$Get_Chromosome" ]]; then
            if [[ -z "$interval_list" ]]; then
  						$GATK --java-options "$Xmxg -Djava.io.tmpdir=../snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ../refgenomes/$ref1 -I ${i%.f*}_${ref1%.f*}_precall.bam -ploidy $ploidy -O ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality $minmapq --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) &&
              rm wait
            else
              $GATK --java-options "$Xmxg -Djava.io.tmpdir=../snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ../refgenomes/$ref1 -L ${projdir}/${interval_list} -I ${i%.f*}_${ref1%.f*}_precall.bam -ploidy $ploidy -O ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality $minmapq --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) &&
  						wait
            fi
					else
						echo $Get_Chromosome | tr ',' '\n' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ../refgenomes/${ref1%.fasta}.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2":1-"$3}' > ../refgenomes/${ref1%.fasta}.list
						$GATK --java-options "$Xmxg -Djava.io.tmpdir=../snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ../refgenomes/$ref1 -L ../refgenomes/${ref1%.fasta}.list -I ${i%.f*}_${ref1%.f*}_precall.bam -ploidy $ploidy -O ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality $minmapq --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) &&
						wait
					fi
          wait
					mv "${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz" "${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf.gz" &&
					rm "${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz.tbi" &&
          gunzip ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf.gz &&
          $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" IndexFeatureFile -I ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf &&
					wait
			fi
		fi
		mv ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam* ${projdir}/preprocess/processed/ 2> /dev/null &&
    rm ${i%.f*}_${ref1%.f*}_precall.bam* 2> /dev/null &&
    wait ) &
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
			wait
		fi
	done
	wait
  cd ${projdir}/preprocess/
	printf "variant calling completed on ${samples_list}" > ${projdir}/call1_${samples_list}

	if [[ "$samples_list" == "samples_list_node_1.txt" ]] && [[ -z "$(ls ${projdir}/snpcall/${pop}_${ploidy}x_raw.vcf* 2> /dev/null)" ]]; then
		call1=$(ls ${projdir}/call1_samples_list_node_* | wc -l)
		while [[ "$call1" -lt $nodes ]]; do sleep 300; call1=$(ls ${projdir}/call1_samples_list_node_* | wc -l); done
		if [[ $call1 == $nodes ]]; then
			cd ${projdir}/snpcall
      mkdir -p cohorts_1
      mv *.g.vcf* ./cohorts_1/

			for dir in cohorts*/; do
				cd $dir
				j=--variant; input=""; k=""
				for i in $(ls *.g.vcf 2> /dev/null); do
					k="${j} ${i}"; input="${input} ${k}"
				done
				if [[ -z "$Get_Chromosome" ]]; then
					Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/${ref1%.*}.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -n | awk '{print $1}' | awk -v pat=${ref1%.f*} '$0 ~ pat')
				else
					Get2_Chromosome=$(echo $Get_Chromosome | tr ',' '\n')
				fi
				if [[ ! -z "$Exclude_Chromosome" ]]; then
					for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
						Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
					done
				fi
				for selchr in $Get2_Chromosome; do
          if [[ -z "$interval_list" ]]; then
  					$GATK --java-options "$Xmx2 -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$threads" GenomicsDBImport ${input} -L ${selchr} --genomicsdb-workspace-path ${pop}_${ploidy}x_${selchr}_raw --genomicsdb-shared-posixfs-optimizations true --batch-size 50 --merge-input-intervals &&
  					wait
          else
            cat ${projdir}/${interval_list} | grep $selchr > ${projdir}/variant_intervals_${selchr}.list &&
            $GATK --java-options "$Xmx2 -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$threads" GenomicsDBImport ${input} -L ${projdir}/variant_intervals_${selchr}.list --genomicsdb-workspace-path ${pop}_${ploidy}x_${selchr}_raw --genomicsdb-shared-posixfs-optimizations true --batch-size 50 --merge-input-intervals &&
            rm ${projdir}/variant_intervals_${selchr}.list &&
            wait
          fi
				done
				for selchr in $Get2_Chromosome; do
					if test ! -f ${pop}_${ploidy}x_${selchr}_raw.vcf.gz; then
            if [[ -z "$interval_list" ]]; then
  						$GATK --java-options "$Xmx2 -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$threads" GenotypeGVCFs -R ${projdir}/refgenomes/$ref1 -L ${selchr} -V gendb://${pop}_${ploidy}x_${selchr}_raw -O ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz && \
  						rm -r ${pop}_${ploidy}x_${selchr}_raw && \
  						mv ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz ${pop}_${ploidy}x_${selchr}_raw.vcf.gz
  						mv ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz.tbi ${pop}_${ploidy}x_${selchr}_raw.vcf.gz.tbi &&
  						wait
            else
              cat ${projdir}/${interval_list} | grep $selchr > ${projdir}/variant_intervals_${selchr}.list &&
              $GATK --java-options "$Xmx2 -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$threads" GenotypeGVCFs -R ${projdir}/refgenomes/$ref1 -L ${projdir}/variant_intervals_${selchr}.list -V gendb://${pop}_${ploidy}x_${selchr}_raw -O ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz && \
              rm ${projdir}/variant_intervals_${selchr}.list &&
              rm -r ${pop}_${ploidy}x_${selchr}_raw &&
              mv ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz ${pop}_${ploidy}x_${selchr}_raw.vcf.gz
              mv ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz.tbi ${pop}_${ploidy}x_${selchr}_raw.vcf.gz.tbi &&
              wait
            fi
					fi
					if LC_ALL=C gzip -l ${pop}_${ploidy}x_${selchr}_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
						:
					else
						rm ../cohorts*/${pop}_${ploidy}x_*_raw.vcf.gz*
						rm ../cohorts*/${pop}_${ploidy}x_raw.vcf*
						rm ../${pop}_${ploidy}x_raw_cohorts*.vcf*
						echo -e "${magenta}- \n- SNP calling failed probably due to insufficient memory ${white}\n"
						echo -e "${magenta}- \n- Exiting pipeline in 5 seconds ${white}\n"
						sleep 5 && exit 1
					fi
          wait
				done
				wait
				for g in $(ls ${pop}_${ploidy}x_*_raw.vcf.gz); do
					gunzip $g
				done
				grep -h '^#' ${pop}_${ploidy}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt
				cat ${pop}_${ploidy}x_*_raw.vcf | awk '!/^#/' > all.vcf
				cat vcf_header.txt all.vcf > ${pop}_${ploidy}x_raw.vcf
				rm vcf_header.txt all.vcf
				rm ${pop}_${ploidy}x_*_raw.vcf ${pop}_${ploidy}x_*_raw.vcf.gz.tbi
				$bcftools view -I ${pop}_${ploidy}x_raw.vcf -O z -o ${pop}_${ploidy}x_raw.vcf.gz
				$bcftools index ${pop}_${ploidy}x_raw.vcf.gz

				$bcftools annotate -x FORMAT/PL ${pop}_${ploidy}x_raw.vcf.gz > ../${pop}_${ploidy}x_raw_${dir%/}.vcf
				cd ../
				$bcftools view -I ${pop}_${ploidy}x_raw_${dir%/}.vcf -O z -o ${pop}_${ploidy}x_raw_${dir%/}.vcf.gz
				$bcftools index ${pop}_${ploidy}x_raw_${dir%/}.vcf.gz
				wait
			done
			wait
			if [[ `ls -1 *cohorts*.vcf.gz 2>/dev/null | wc -l` -gt 1 ]]; then
				$bcftools merge *cohorts*.vcf.gz --force-samples -m all > ${pop}_${ploidy}x_raw.vcf
			else
				cp *cohorts*.vcf.gz ${pop}_${ploidy}x_raw.vcf.gz
				gunzip ${pop}_${ploidy}x_raw.vcf.gz
			fi

			if [[ "$keep_gVCF" == true ]]; then
				mkdir -p keep_gVCF
				mv ./cohorts*/*.g.vcf ./keep_gVCF
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

}
cd $projdir
# while [[ ! -f "$projdir/alignment_summaries/refgenome_paralogs.txt" ]]; do sleep 30; done
# sleep 10
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

if [[ -d "snpfilter" ]]; then
	number_snpfilter=$( ls -d snpfilter* | wc -l )
	mv snpfilter snpfilter_"${number_snpfilter}"
fi

mkdir snpfilter
cd snpfilter

large_numerous_chrom () {
  #retrieve the original SNP positions from unsplit chromosomes/pseudomolecules.
  cd $projdir
  cd snpcall
  if [[ "$checksplit" -gt 0 ]]; then
  	for ln in $(ls *_raw0.vcf 2> /dev/null); do
  		cp $ln ${ln%.vcf}_positions_split.vcf
  		awk '/^#CHROM/{close("file.vcf"f);f++}{print $0 > "file"f}' $ln
  		awk 'NR==1{print}' file1 | cat file - > file0.vcf
  		for j in $(seq 1 10); do
  			awk 'NR>1{print}' file1 | awk '{print $1,"\t",$2,"\t",$0}' | awk '{gsub(/_x0/,"\t"); print}' | \
  			awk -v splits=$j -F '\t' 'BEGIN{OFS="\t"} $2 ~ splits {$3=$3+500000000}1' | awk 'BEGIN{OFS="\t"} !($2="")' | awk 'BEGIN{OFS="\t"} !($3="")' | \
  			awk 'BEGIN{OFS="\t"} !($3="")' | awk 'BEGIN{OFS="\t"} !($3="")' | awk '{gsub(/\t\t/,"\t"); print }' | \
  			sort -k1,1 -k2n,2 | awk NF > file1.vcf
  		done
  		cat file0.vcf file1.vcf > $ln
  		rm file file1 file0.vcf file1.vcf
  	done

  	cd ${projdir}/alignment_summaries
  	lnr=refgenome_paralogs.txt
  	cp $lnr ${lnr%.txt}_positions_split.txt
  	awk '/^CHROM/{close("file.txt"f);f++}{print $0 > "file"f}' $lnr
  	awk 'NR==1{print}' file1 | cat file - > file0.txt
  	for j in $(seq 1 10); do
  		awk 'NR>1{print}' file1 | awk '{print $1,"\t",$2,"\t",$0}' | awk '{gsub(/_x0/,"\t"); print}' | \
  		awk -v splits=$j -F '\t' 'BEGIN{OFS="\t"} $2 ~ splits {$3=$3+500000000}1' | awk 'BEGIN{OFS="\t"} !($2="")' | awk 'BEGIN{OFS="\t"} !($3="")' | \
  		awk 'BEGIN{OFS="\t"} !($3="")' | awk 'BEGIN{OFS="\t"} !($3="")' | awk '{gsub(/\t\t/,"\t"); print }' | \
  		sort -k1,1 -k2n,2 | awk NF > file1.txt
  	done
  	cat file0.txt file1.txt > $lnr
  	rm file file1 file0.txt file1.txt
  	cd ${projdir}/snpcall
  else
  	:
  fi

  #retrieve the original SNP positions from conigs/scaffolds before concatenation
  cd $projdir
  cd snpcall
  export ncontigscaffold=$(grep '>' ${projdir}/refgenomes/${ref1%.fasta}_original.fasta &> /dev/null | wc -l)
  if [[ ! -f ./index_code.txt ]]; then
  	if [[ $ncontigscaffold -gt 3000 ]]; then
  		echo -e "${magenta}- retrieving SNP positions based on contigs/scaffold annotation ${white}\n"
  		for nc in $(ls *_raw0.vcf 2> /dev/null); do
  			if [[ "${nc}" =~ "vcf" ]]; then
  				awk '/^#CHROM/{close("file.vcf"f);f++}{print $0 > "file"f}' $nc
  				awk 'NR==1{print}' file1 | cat file - > file0.txt; rm file
  				awk 'NR>1{print $1,"\t",$2}' file1 > file_snp_stitch.txt
  				awk 'NR>1{print $0}' file1 > file1.txt; rm file1
  			fi
  			awk 'NR%10000==1{x="splitindex"++i;}{print > x}'  file_snp_stitch.txt & PID=$!
			  wait $PID
  			for s in $( ls splitindex* ); do
  				touch ${s}_index_code.txt
  				while IFS="" read -r p || [ -n "$p" ]; do
  					fchr=$(echo $p | awk '{print $1}'); fpos=$(echo $p | awk '{print $2}')
  					scf=$(awk -v fchr=$fchr '$1 == fchr' ${projdir}/refgenomes/contigscaffold_index.txt | \
  					awk -v t=$fpos 'BEGIN {var=t; highest=0}{ j = $NF;if ( j <= var && j > highest ) { highest=j} } END {print highest}' )
  					awk -v fchr=$fchr '$1 == fchr' ${projdir}/refgenomes/contigscaffold_index.txt | awk -v scf=$scf -v fchr=$fchr -v fpos=$fpos '$3 == scf {print fchr"\t"fpos"\t"$0}' >> ${s}_index_code.txt
  				done < $s
  			done
  			wait
  			cat *_index_code.txt > index_code.txt && rm splitindex*
  			awk 'BEGIN {FS=OFS="\t"}{ $5 = $2 - $5 } 1' index_code.txt > file2.txt
  			awk 'BEGIN {FS=OFS="\t"}!($3="")' file2.txt | awk '{gsub(/\t\t/,"\t"); print}' | awk '{gsub(/ /,""); print}' > index_code.txt
  			paste -d "\t" index_code.txt file1.txt | awk -F"\t" 'BEGIN {FS=OFS="\t"}{$1=$2=$5=$6=""; print $0}' | tr -s '\t' | sed -e 's/^[ \t]*//' > file2.txt
  			cat file0.txt file2.txt > $nc
  			rm file0.txt file1.txt file2.txt index_code.txt
  		done
  		wait
  		rm file_snp_stitch.txt


  		cd ${projdir}/alignment_summaries
  		awk 'NR==1{print}' refgenome_paralogs.txt > file0.txt
  		awk 'NR>1{print $1,"\t",$2}' refgenome_paralogs.txt > file_snp_stitch.txt
  		awk 'NR>1{print $0}' refgenome_paralogs.txt > file1.txt
  		awk 'NR%10000==1{x="splitindex"++i;}{print > x}'  file_snp_stitch.txt & PID=$!
		  wait $PID
  		for s in $(ls splitindex*); do
  			while IFS="" read -r p || [ -n "$p" ]; do
  			fchr=$(echo $p | awk '{print $1}'); fpos=$(echo $p | awk '{print $2}')
  			scf=$(awk -v fchr=$fchr '$1 == fchr' ${projdir}/refgenomes/contigscaffold_index.txt | \
  			awk -v t=$fpos 'BEGIN {var=t; highest=0}{ j = $NF;if ( j <= var && j > highest ) { highest=j} } END {print highest}' )
  			awk -v fchr=$fchr '$1 == fchr' ${projdir}/refgenomes/contigscaffold_index.txt | awk -v scf=$scf -v fchr=$fchr -v fpos=$fpos '$3 == scf {print fchr"\t"fpos"\t"$0}' >> ${s}_index_code.txt
  			done < $s
  		done
  		wait
  		cat *_index_code.txt > index_code.txt && rm splitindex*
  		awk 'BEGIN {FS=OFS="\t"}{ $5 = $2 - $5 } 1' index_code.txt > file2.txt
  		awk 'BEGIN {FS=OFS="\t"}!($3="")' file2.txt | awk '{gsub(/\t\t/,"\t"); print}' | awk '{gsub(/ /,""); print}' > index_code.txt
  		paste -d "\t" index_code.txt file1.txt | awk -F"\t" 'BEGIN {FS=OFS="\t"}{$1=$2=$5=$6=""; print $0}' | tr -s '\t' | sed -e 's/^[ \t]*//' > file2.txt
  		cat file0.txt file2.txt > refgenome_paralogs.txt
  		rm file0.txt file1.txt file2.txt file_snp_stitch.txt
  		cd ${projdir}/snpcall
  	fi
  fi
}

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
	genotype_missingness=0.2
else
	genotype_missingness=$( echo $genotype_missingness | awk '{gsub(/,/," ")}1' )
fi
if [ -z $sample_missingness ]; then
	sample_missingness=0.2
else
	sample_missingness=$( echo $sample_missingness | awk '{gsub(/,/," ")}1' )
fi
if [ -z $exclude_samples ]; then
	exclude_samples=NULL
fi
if [ -z $minRD_1x ]; then
	minRD_1x=2
fi
if [ -z $minRD_2x ]; then
	minRD_6x=6
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
if [ -z "$(ls -A *_DP_GT.txt 2>/dev/null)" ]; then
	if [ -z "$(ls -A *_x.vcf 2>/dev/null)" ]; then
		if [ -z "$(ls -A *_raw.vcf 2>/dev/null)" ]; then
			for g in $(ls *_raw.vcf.gz 2> /dev/null); do gunzip $g;	done
		fi
	fi
fi
wait

file1xG=$( if [ "$(ls -A *_1x_DP_GT.txt 2>/dev/null)" ]; then ls *_1x_DP_GT.txt | wc -l;  else echo 0; fi )
file1xV=$( if [ "$(ls -A *_1x_raw.vcf 2>/dev/null)" ]; then ls *_1x_raw.vcf | head -n 1 | wc -l; else echo 0; fi )
if [[ "${file1xG}" -lt 1 ]]; then
	if [[ "${file1xV}" -gt 0 ]]; then
    if test ! -f ${projdir}/vcf1x_trimmed.txt; then
  		for i in $(ls *_1x_raw.vcf 2> /dev/null); do
        samz=$(( $(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') / 100 * 80 ))
        awk -v pat="0:0,0:0" -v samz="$samz" 'gsub(pat,pat) < samz' $i | awk -v pat=".:0,0:0"  -v samz="$samz" 'gsub(pat,pat) < samz' > ${i}.tmp
        mv ${i}.tmp ${i}
  			$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$ref1 -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
  			wait
  			large_numerous_chrom &>> ${projdir}/log.out
  			wait
  			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf
  			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%100000==2{x=file"1x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
  		  wait $PID
  			rm ${i%.vcf}0.vcf* ${i%.vcf}trim.vcf*
  		done
      :> ${projdir}/vcf1x_trimmed.txt
    fi
    samz=$(( $(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') / 100 * 80 ))
    for ptrimvcf in $(ls *rawSPLIT*.vcf); do
      awk -v pat="0:0,0:0" -v samz="$samz" 'gsub(pat,pat) < samz' $ptrimvcf | awk -v pat=".:0,0:0"  -v samz="$samz" 'gsub(pat,pat) < samz' > ${ptrimvcf}.tmp
      mv ${ptrimvcf}.tmp ${ptrimvcf}
    done
		Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 1x "${GBSapp_dir}/tools/R" "1"
    rm ${projdir}/vcf1x_trimmed.txt 2> /dev/null
	fi
fi
wait
file2xG=$( if [ "$(ls -A *_2x_DP_GT.txt 2>/dev/null)" ]; then ls *_2x_DP_GT.txt | wc -l;  else echo 0; fi )
file2xV=$( if [ "$(ls -A *_2x_raw.vcf 2>/dev/null)" ]; then ls *_2x_raw.vcf | head -n 1 | wc -l; else echo 0; fi )
if [[ "${file2xG}" -lt 1 ]]; then
	if [[ "${file2xV}" -gt 0 ]]; then
    if test ! -f ${projdir}/vcf2x_trimmed.txt; then
  		for i in $(ls *_2x_raw.vcf 2> /dev/null); do
        samz=$(( $(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') / 100 * 80 ))
        awk -v pat="0/0:0,0:0" -v samz="$samz" 'gsub(pat,pat) < samz' $i | awk -v pat="./.:0,0:0"  -v samz="$samz" 'gsub(pat,pat) < samz' > ${i}.tmp
        mv ${i}.tmp ${i}
  			$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$ref1 -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
  			wait
  			large_numerous_chrom &>> ${projdir}/log.out
  			wait
  			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf
  			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%100000==2{x=file"2x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
  			wait $PID
  			rm ${i%.vcf}0.vcf* ${i%.vcf}trim.vcf*
  		done
      :> ${projdir}/vcf2x_trimmed.txt
    fi
    samz=$(( $(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') / 100 * 80 ))
    for ptrimvcf in $(ls *rawSPLIT*.vcf); do
      awk -v pat="0/0:0,0:0" -v samz="$samz" 'gsub(pat,pat) < samz' $ptrimvcf | awk -v pat="./.:0,0:0"  -v samz="$samz" 'gsub(pat,pat) < samz' > ${ptrimvcf}.tmp
      mv ${ptrimvcf}.tmp ${ptrimvcf}
    done
		Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 2x "${GBSapp_dir}/tools/R" "1"
    rm ${projdir}/vcf2x_trimmed.txt 2> /dev/null
	fi
fi
wait
file4xG=$( if [ "$(ls -A *_4x_DP_GT.txt 2>/dev/null)" ]; then ls *_4x_DP_GT.txt | wc -l;  else echo 0; fi )
file4xV=$( if [ "$(ls -A *_4x_raw.vcf 2>/dev/null)" ]; then ls *_4x_raw.vcf | head -n 1 | wc -l; else echo 0; fi )
if [[ "${file4xG}" -lt 1 ]]; then
	if [[ "${file4xV}" -gt 0 ]]; then
    if test ! -f ${projdir}/vcf4x_trimmed.txt; then
  		for i in $(ls *_4x_raw.vcf 2> /dev/null); do
        samz=$(( $(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') / 100 * 80 ))
        awk -v pat="0/0/0/0:0,0:0" -v samz="$samz" 'gsub(pat,pat) < samz' $i | awk -v pat="./././.:0,0:0"  -v samz="$samz" 'gsub(pat,pat) < samz' > ${i}.tmp
        mv ${i}.tmp ${i}
  			$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$ref1 -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
  			wait
  			large_numerous_chrom &>> ${projdir}/log.out
  			wait
  			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf
  			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%100000==2{x=file"4x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
  			wait $PID
  			wait
  			rm ${i%.vcf}0.vcf* ${i%.vcf}trim.vcf*
  		done
      :> ${projdir}/vcf4x_trimmed.txt
    fi
    samz=$(( $(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') / 100 * 80 ))
    for ptrimvcf in $(ls *rawSPLIT*.vcf); do
      awk -v pat="0/0/0/0:0,0:0" -v samz="$samz" 'gsub(pat,pat) < samz' $ptrimvcf | awk -v pat="./././.:0,0:0"  -v samz="$samz" 'gsub(pat,pat) < samz' > ${ptrimvcf}.tmp
      mv ${ptrimvcf}.tmp ${ptrimvcf}
    done
		Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 4x "${GBSapp_dir}/tools/R" "1"
    rm ${projdir}/vcf4x_trimmed.txt 2> /dev/null
	fi
fi
wait
file6xG=$( if [ "$(ls -A *_6x_DP_GT.txt 2>/dev/null)" ]; then ls *_6x_DP_GT.txt | wc -l;  else echo 0; fi )
file6xV=$( if [ "$(ls -A *_6x_raw.vcf 2>/dev/null)" ]; then ls *_6x_raw.vcf | head -n 1 | wc -l; else echo 0; fi )
if [[ "${file6xG}" -lt 1 ]]; then
	if [[ "${file6xV}" -gt 0 ]]; then
    if test ! -f ${projdir}/vcf6x_trimmed.txt; then
  		for i in $(ls *_6x_raw.vcf 2> /dev/null); do
        samz=$(( $(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') / 100 * 80 ))
        awk -v pat="0/0/0/0/0/0:0,0:0" -v samz="$samz" 'gsub(pat,pat) < samz' $i | awk -v pat="./././././.:0,0:0"  -v samz="$samz" 'gsub(pat,pat) < samz' > ${i}.tmp
        mv ${i}.tmp ${i}
  			$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$ref1 -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
  			wait
  			large_numerous_chrom &>> ${projdir}/log.out
  			wait
  			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf
  			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%100000==2{x=file"6x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
  			wait $PID
  			rm ${i%.vcf}0.vcf* ${i%.vcf}trim.vcf*
  		done
      :> ${projdir}/vcf6x_trimmed.txt
    fi
    samz=$(( $(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') / 100 * 80 ))
    for ptrimvcf in $(ls *rawSPLIT*.vcf); do
      awk -v pat="0/0/0/0/0/0:0,0:0" -v samz="$samz" 'gsub(pat,pat) < samz' $ptrimvcf | awk -v pat="./././././.:0,0:0"  -v samz="$samz" 'gsub(pat,pat) < samz' > ${ptrimvcf}.tmp
      mv ${ptrimvcf}.tmp ${ptrimvcf}
    done
		Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 6x "${GBSapp_dir}/tools/R" "1"
    rm ${projdir}/vcf6x_trimmed.txt 2> /dev/null
	fi
fi
wait
file8xG=$( if [ "$(ls -A *_8x_DP_GT.txt 2>/dev/null)" ]; then ls *_8x_DP_GT.txt | wc -l;  else echo 0; fi )
file8xV=$( if [ "$(ls -A *_8x_raw.vcf 2>/dev/null)" ]; then ls *_8x_raw.vcf | head -n 1 | wc -l; else echo 0; fi )
if [[ "${file8xG}" -lt 1 ]]; then
	if [[ "${file8xV}" -gt 0 ]]; then
    if test ! -f ${projdir}/vcf8x_trimmed.txt; then
  		for i in $(ls *_8x_raw.vcf 2> /dev/null); do
        samz=$(( $(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') / 100 * 80 ))
        awk -v pat="0/0/0/0/0/0/0/0:0,0:0" -v samz="$samz" 'gsub(pat,pat) < samz' $i | awk -v pat="./././././././.:0,0:0"  -v samz="$samz" 'gsub(pat,pat) < samz' > ${i}.tmp
        mv ${i}.tmp ${i}
  			$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$ref1 -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
  			wait
  			large_numerous_chrom &>> ${projdir}/log.out
  			wait
  			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf
  			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%100000==2{x=file"8x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
  			wait $PID ${i%.vcf}trim.vcf*
  			wait
  			rm ${i%.vcf}0.vcf*
  		done
      :> ${projdir}/vcf8x_trimmed.txt
    fi
    samz=$(( $(wc -l ${projdir}/samples_list_node_* | awk '{print $1}') / 100 * 80 ))
    for ptrimvcf in $(ls *rawSPLIT*.vcf); do
      awk -v pat="0/0/0/0/0/0/0/0:0,0:0" -v samz="$samz" 'gsub(pat,pat) < samz' $ptrimvcf | awk -v pat="./././././././.:0,0:0"  -v samz="$samz" 'gsub(pat,pat) < samz' > ${ptrimvcf}.tmp
      mv ${ptrimvcf}.tmp ${ptrimvcf}
    done
		Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 8x "${GBSapp_dir}/tools/R" "1"
    rm ${projdir}/vcf8x_trimmed.txt 2> /dev/null
	fi
fi
wait
filetest=*x.vcf*
if [ -z "$(ls -A *x.vcf* 2>/dev/null)" ]; then
	for v in $(ls *_DP_GT.txt 2> /dev/null); do
		vcfdose=${v%_DP*}; vcfdose=${vcfdose#*_}; out=$(ls *${vcfdose}_raw.vcf)
		for raw in $out; do
			grep '^#' $raw  > ${raw%_raw.vcf}.vcf
			awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$1,$2]){print b}}' $raw $v >> ${raw%_raw.vcf}.vcf
		done
	done
	wait
fi

if [ "$(ls -A *.vcf 2>/dev/null)" ]; then
	for v in $(ls *.vcf 2> /dev/null); do gzip $v; done
fi
wait


cd $projdir
cd samples
window1=$(ls -S | head -1 | xargs zcat -fq | awk '{ print length }' | sort -n | tail -1)
window=$((window1 + 20))
cd ${projdir}


for smiss in ${sample_missingness[@]}; do
for gmiss in ${genotype_missingness[@]}; do
if [[ -z "$p1" ]]; then
	if [ -d "${projdir}/snpfilter/1x" ]; then
		cd ${projdir}/snpfilter
		cp -r 1x 1x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./1x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_1x.R $pop $gmiss $smiss $minRD_1x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome_number
		wait
		rm ${pop}_1x_rawRD${minRD_1x}_DP_GT.txt ${pop}_1x_DP_GT.txt ${pop}_1x_rd${minRD_1x}.txt 2> /dev/null
		mkdir visualizations; mv *.tiff ./visualizations/

    if [[ "$variant_intervals" == true ]]; then
      mkdir -p variant_intervals
      cd variant_intervals
      while IFS="" read -r p || [ -n "$p" ]; do
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
        for vbreak in vbreak_*; do
          awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
          awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt
          wait
        done
        sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
        mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
        wait && rm vbreak_*
        wait
      done < <(awk 'NR>1{print $2}' ../${pop}_1x_rd${minRD_1x}_maf${maf}_dose.txt | sort | uniq)
      cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
      rm variant_intervals_*.txt && cd ../
      rm -rf variant_intervals
    fi

    # Extract sequence context of variants
    if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        awk '{print $1"\t"$2"\t"$3}' ${pop}_1x_rd${minRD_1x}_maf${maf}_dose.txt > snplist.txt &&
        awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
        awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
        for wint in {100,200}; do
          awk -v wint=$wint '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
          awk -v wint=$wint '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
          wait
        done
        awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
        mv snplist_haps.tmp snplist_haps.txt &&
        rm snplist_round.txt
      fi
      if [[ ! -z "$RE1" ]]; then
        mkdir -p seq_context_${RE1}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
          rm ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
        done
      fi
      if [[ ! -z "$RE2" ]]; then
        mkdir -p seq_context_${RE2}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
          rm ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
        done
      fi
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        mkdir consensus_seq_context
        for i in $(cat ${projdir}/samples_list_node_*); do
          cat ./seq_context_${RE1}/${i%.f*}* ./seq_context_${RE2}/${i%.f*}* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
          maxlen=$(awk '{print $2}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}')
          awk -v maxlen=$maxlen '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          mv ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta
          for aln in $(awk '{print $1}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | sort | uniq); do
            awk -v aln=$aln '$1 == aln {print}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
            bash $cons -sequence ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp -outseq ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp
            cat ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp | awk -v aln="${i%.f*}_${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' >> ./consensus_seq_context/${i%.f*}_seqcontext.cons
            rm ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          done
        done
        for nn in $(ls ./consensus_seq_context/*cons);do
          for run in {1..10}; do
            awk '{gsub(/n$|N$/,"");}1' $nn > ${nn}.tmp && mv ${nn}.tmp ${nn}
          done
        done
        seqid=$(cat ./consensus_seq_context/*cons | grep '>' | awk '{gsub(/>/,""); gsub(/_Chr/,"\t_Chr");}1' | awk '{print $2}' | sort | uniq)
        mkdir -p ./consensus_seq_context/seqid_combine_samples
        for splitseqid in $seqid; do
          filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}"
          cat ./consensus_seq_context/*cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
          awk -v pat=$splitseqid '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/${filename}
          minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{print length}' | sort -n | head -n 1)
          awk -v minlen=$minlen '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{gsub(/\t/,"\n");}1' > ./consensus_seq_context/seqid_combine_samples/${filename}.tmp
          mv ./consensus_seq_context/seqid_combine_samples/${filename}.tmp ./consensus_seq_context/seqid_combine_samples/${filename}
        done
        mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/
        mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*cons ./consensus_seq_context/sample_consensus_seqs
        mv seq_context_TGCAT ./consensus_seq_context/
        mv seq_context_CATG ./consensus_seq_context/
      fi
    fi
		find . -type f -empty -delete
	fi
	wait
	if [ -d "${projdir}/snpfilter/2x" ]; then
		cd ${projdir}/snpfilter
		cp -r 2x 2x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./2x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_2x.R $pop $gmiss $smiss $minRD_2x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome_number
		wait
		rm ${pop}_2x_rawRD${minRD_2x}_DP_GT.txt ${pop}_2x_DP_GT.txt ${pop}_2x_rd${minRD_2x}.txt 2> /dev/null
		mkdir visualizations; mv *.tiff ./visualizations/

    if [[ "$variant_intervals" == true ]]; then
      mkdir -p variant_intervals
      cd variant_intervals
      while IFS="" read -r p || [ -n "$p" ]; do
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
        for vbreak in vbreak_*; do
          awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
          awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt
          wait
        done
        sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
        mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
        wait && rm vbreak_*
        wait
      done < <(awk 'NR>1{print $2}' ../${pop}_2x_rd${minRD_2x}_maf${maf}_dose.txt | sort | uniq)
      cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
      rm variant_intervals_*.txt && cd ../
      rm -rf variant_intervals
    fi

    # Extract sequence context of variants
    if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        awk '{print $1"\t"$2"\t"$3}' ${pop}_2x_rd${minRD_2x}_maf${maf}_dose.txt > snplist.txt &&
        awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
        awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
        for wint in {100,200}; do
          awk -v wint=$wint '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
          awk -v wint=$wint '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
          wait
        done
        awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
        mv snplist_haps.tmp snplist_haps.txt &&
        rm snplist_round.txt
      fi
      if [[ ! -z "$RE1" ]]; then
        mkdir -p seq_context_${RE1}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
          rm ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
        done
      fi
      if [[ ! -z "$RE2" ]]; then
        mkdir -p seq_context_${RE2}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
          rm ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
        done
      fi
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        mkdir consensus_seq_context
        for i in $(cat ${projdir}/samples_list_node_*); do
          cat ./seq_context_${RE1}/${i%.f*}* ./seq_context_${RE2}/${i%.f*}* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
          maxlen=$(awk '{print $2}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}')
          awk -v maxlen=$maxlen '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          mv ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta
          for aln in $(awk '{print $1}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | sort | uniq); do
            awk -v aln=$aln '$1 == aln {print}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
            bash $cons -sequence ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp -outseq ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp
            cat ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp | awk -v aln="${i%.f*}_${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' >> ./consensus_seq_context/${i%.f*}_seqcontext.cons
            rm ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          done
        done
        for nn in $(ls ./consensus_seq_context/*cons);do
          for run in {1..10}; do
            awk '{gsub(/n$|N$/,"");}1' $nn > ${nn}.tmp && mv ${nn}.tmp ${nn}
          done
        done
        seqid=$(cat ./consensus_seq_context/*cons | grep '>' | awk '{gsub(/>/,""); gsub(/_Chr/,"\t_Chr");}1' | awk '{print $2}' | sort | uniq)
        mkdir -p ./consensus_seq_context/seqid_combine_samples
        for splitseqid in $seqid; do
          filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}"
          cat ./consensus_seq_context/*cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
          awk -v pat=$splitseqid '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/${filename}
          minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{print length}' | sort -n | head -n 1)
          awk -v minlen=$minlen '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{gsub(/\t/,"\n");}1' > ./consensus_seq_context/seqid_combine_samples/${filename}.tmp
          mv ./consensus_seq_context/seqid_combine_samples/${filename}.tmp ./consensus_seq_context/seqid_combine_samples/${filename}
        done
        mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/
        mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*cons ./consensus_seq_context/sample_consensus_seqs
        mv seq_context_TGCAT ./consensus_seq_context/
        mv seq_context_CATG ./consensus_seq_context/
      fi
    fi
		find . -type f -empty -delete
	fi
	wait
	if [ -d "${projdir}/snpfilter/4x" ]; then
		cd ${projdir}/snpfilter
		cp -r 4x 4x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./4x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_4x.R $pop $gmiss $smiss $minRD_4x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome_number
		wait
		rm ${pop}_4x_rawRD${minRD_4x}_DP_GT.txt ${pop}_4x_DP_GT.txt ${pop}_4x_rd${minRD_4x}.txt 2> /dev/null
		mkdir visualizations; mv *.tiff ./visualizations/

    if [[ "$variant_intervals" == true ]]; then
      mkdir -p variant_intervals
      cd variant_intervals
      while IFS="" read -r p || [ -n "$p" ]; do
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
          awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt
          wait
        done
        sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
        mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
        wait && rm vbreak_*
        wait
      done < <(awk 'NR>1{print $2}' ../${pop}_4x_rd${minRD_4x}_maf${maf}_dose.txt | sort | uniq)
      cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
      rm variant_intervals_*.txt && cd ../
      rm -rf variant_intervals
    fi

    # Extract sequence context of variants
    if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        awk '{print $1"\t"$2"\t"$3}' ${pop}_4x_rd${minRD_4x}_maf${maf}_dose.txt > snplist.txt &&
        awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
        awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
        for wint in {100,200}; do
          awk -v wint=$wint '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
          awk -v wint=$wint '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
          wait
        done
        awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
        mv snplist_haps.tmp snplist_haps.txt &&
        rm snplist_round.txt
      fi
      if [[ ! -z "$RE1" ]]; then
        mkdir -p seq_context_${RE1}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
          rm ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
        done
      fi
      if [[ ! -z "$RE2" ]]; then
        mkdir -p seq_context_${RE2}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
          rm ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
        done
      fi
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        mkdir consensus_seq_context
        for i in $(cat ${projdir}/samples_list_node_*); do
          cat ./seq_context_${RE1}/${i%.f*}* ./seq_context_${RE2}/${i%.f*}* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
          maxlen=$(awk '{print $2}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}')
          awk -v maxlen=$maxlen '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          mv ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta
          for aln in $(awk '{print $1}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | sort | uniq); do
            awk -v aln=$aln '$1 == aln {print}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
            bash $cons -sequence ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp -outseq ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp
            cat ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp | awk -v aln="${i%.f*}_${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' >> ./consensus_seq_context/${i%.f*}_seqcontext.cons
            rm ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          done
        done
        for nn in $(ls ./consensus_seq_context/*cons);do
          for run in {1..10}; do
            awk '{gsub(/n$|N$/,"");}1' $nn > ${nn}.tmp && mv ${nn}.tmp ${nn}
          done
        done
        seqid=$(cat ./consensus_seq_context/*cons | grep '>' | awk '{gsub(/>/,""); gsub(/_Chr/,"\t_Chr");}1' | awk '{print $2}' | sort | uniq)
        mkdir -p ./consensus_seq_context/seqid_combine_samples
        for splitseqid in $seqid; do
          filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}"
          cat ./consensus_seq_context/*cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
          awk -v pat=$splitseqid '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/${filename}
          minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{print length}' | sort -n | head -n 1)
          awk -v minlen=$minlen '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{gsub(/\t/,"\n");}1' > ./consensus_seq_context/seqid_combine_samples/${filename}.tmp
          mv ./consensus_seq_context/seqid_combine_samples/${filename}.tmp ./consensus_seq_context/seqid_combine_samples/${filename}
        done
        mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/
        mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*cons ./consensus_seq_context/sample_consensus_seqs
        mv seq_context_TGCAT ./consensus_seq_context/
        mv seq_context_CATG ./consensus_seq_context/
      fi
    fi
		find . -type f -empty -delete
	fi
	wait
	if [ -d "${projdir}/snpfilter/6x" ]; then
		cd ${projdir}/snpfilter
		cp -r 6x 6x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./6x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_6x.R $pop $gmiss $smiss $minRD_6x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome_number
		wait
		rm ${pop}_6x_rawRD${minRD_6x}_DP_GT.txt ${pop}_6x_DP_GT.txt ${pop}_6x_rd${minRD_6x}.txt 2> /dev/null
		mkdir visualizations; mv *.tiff ./visualizations/

    if [[ "$variant_intervals" == true ]]; then
      mkdir -p variant_intervals
      cd variant_intervals
      while IFS="" read -r p || [ -n "$p" ]; do
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
        for vbreak in vbreak_*; do
          awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
          awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt
          wait
        done
        sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
        mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
        wait && rm vbreak_*
        wait
      done < <(awk 'NR>1{print $2}' ../${pop}_6x_rd${minRD_6x}_maf${maf}_dose.txt | sort | uniq)
      cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
      rm variant_intervals_*.txt && cd ../
      rm -rf variant_intervals
    fi

    # Extract sequence context of variants
    if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        awk '{print $1"\t"$2"\t"$3}' ${pop}_6x_rd${minRD_6x}_maf${maf}_dose.txt > snplist.txt &&
        awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
        awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
        for wint in {100,200}; do
          awk -v wint=$wint '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
          awk -v wint=$wint '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
          wait
        done
        awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
        mv snplist_haps.tmp snplist_haps.txt &&
        rm snplist_round.txt
      fi
      if [[ ! -z "$RE1" ]]; then
        mkdir -p seq_context_${RE1}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' | \
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt
        done
      fi
      if [[ ! -z "$RE2" ]]; then
        mkdir -p seq_context_${RE2}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' | \
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt
        done
      fi
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        mkdir consensus_seq_context
        for i in $(cat ${projdir}/samples_list_node_*); do
          cat ./seq_context_${RE1}/${i%.f*}* ./seq_context_${RE2}/${i%.f*}* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
          maxlen=$(awk '{print $2}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}')
          awk -v maxlen=$maxlen '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          mv ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta
          for aln in $(awk '{print $1}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | sort | uniq); do
            awk -v aln=$aln '$1 == aln {print}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
            bash $cons -sequence ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp -outseq ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp
            cat ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp | awk -v aln="${i%.f*}_${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' >> ./consensus_seq_context/${i%.f*}_seqcontext.cons
            rm ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          done
        done
        for nn in $(ls ./consensus_seq_context/*cons);do
          for run in {1..10}; do
            awk '{gsub(/n$|N$/,"");}1' $nn > ${nn}.tmp && mv ${nn}.tmp ${nn}
          done
        done
        seqid=$(cat ./consensus_seq_context/*cons | grep '>' | awk '{gsub(/>/,""); gsub(/_Chr/,"\t_Chr");}1' | awk '{print $2}' | sort | uniq)
        mkdir -p ./consensus_seq_context/seqid_combine_samples
        for splitseqid in $seqid; do
          filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}"
          cat ./consensus_seq_context/*cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
          awk -v pat=$splitseqid '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/${filename}
          minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{print length}' | sort -n | head -n 1)
          awk -v minlen=$minlen '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{gsub(/\t/,"\n");}1' > ./consensus_seq_context/seqid_combine_samples/${filename}.tmp
          mv ./consensus_seq_context/seqid_combine_samples/${filename}.tmp ./consensus_seq_context/seqid_combine_samples/${filename}
        done
        mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/
        mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*cons ./consensus_seq_context/sample_consensus_seqs
        mv seq_context_TGCAT ./consensus_seq_context/
        mv seq_context_CATG ./consensus_seq_context/
      fi
    fi
		find . -type f -empty -delete
	fi
	wait
	if [ -d "${projdir}/snpfilter/8x" ]; then
		cd ${projdir}/snpfilter
		cp -r 8x 8x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./8x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_8x.R $pop $gmiss $smiss $minRD_8x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome_number
		wait
		rm ${pop}_8x_rawRD${minRD_8x}_DP_GT.txt ${pop}_8x_DP_GT.txt ${pop}_8x_rd${minRD_8x}.txt 2> /dev/null
		mkdir visualizations; mv *.tiff ./visualizations/

    if [[ "$variant_intervals" == true ]]; then
      mkdir -p variant_intervals
      cd variant_intervals
      while IFS="" read -r p || [ -n "$p" ]; do
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
        for vbreak in vbreak_*; do
          awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
          awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt
          wait
        done
        sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
        mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
        wait && rm vbreak_*
        wait
      done < <(awk 'NR>1{print $2}' ../${pop}_8x_rd${minRD_8x}_maf${maf}_dose.txt | sort | uniq)
      cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
      rm variant_intervals_*.txt && cd ../
      rm -rf variant_intervals
    fi

    # Extract sequence context of variants
    if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        awk '{print $1"\t"$2"\t"$3}' ${pop}_8x_rd${minRD_8x}_maf${maf}_dose.txt > snplist.txt &&
        awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
        awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
        for wint in {100,200}; do
          awk -v wint=$wint '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
          awk -v wint=$wint '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
          wait
        done
        awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
        mv snplist_haps.tmp snplist_haps.txt &&
        rm snplist_round.txt
      fi
      if [[ ! -z "$RE1" ]]; then
        mkdir -p seq_context_${RE1}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
          rm ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
        done
      fi
      if [[ ! -z "$RE2" ]]; then
        mkdir -p seq_context_${RE2}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
          rm ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
        done
      fi
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        mkdir consensus_seq_context
        for i in $(cat ${projdir}/samples_list_node_*); do
          cat ./seq_context_${RE1}/${i%.f*}* ./seq_context_${RE2}/${i%.f*}* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
          maxlen=$(awk '{print $2}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}')
          awk -v maxlen=$maxlen '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          mv ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta
          for aln in $(awk '{print $1}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | sort | uniq); do
            awk -v aln=$aln '$1 == aln {print}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
            bash $cons -sequence ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp -outseq ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp
            cat ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp | awk -v aln="${i%.f*}_${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' >> ./consensus_seq_context/${i%.f*}_seqcontext.cons
            rm ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          done
        done
        for nn in $(ls ./consensus_seq_context/*cons);do
          for run in {1..10}; do
            awk '{gsub(/n$|N$/,"");}1' $nn > ${nn}.tmp && mv ${nn}.tmp ${nn}
          done
        done
        seqid=$(cat ./consensus_seq_context/*cons | grep '>' | awk '{gsub(/>/,""); gsub(/_Chr/,"\t_Chr");}1' | awk '{print $2}' | sort | uniq)
        mkdir -p ./consensus_seq_context/seqid_combine_samples
        for splitseqid in $seqid; do
          filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}"
          cat ./consensus_seq_context/*cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
          awk -v pat=$splitseqid '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/${filename}
          minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{print length}' | sort -n | head -n 1)
          awk -v minlen=$minlen '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{gsub(/\t/,"\n");}1' > ./consensus_seq_context/seqid_combine_samples/${filename}.tmp
          mv ./consensus_seq_context/seqid_combine_samples/${filename}.tmp ./consensus_seq_context/seqid_combine_samples/${filename}
        done
        mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/
        mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*cons ./consensus_seq_context/sample_consensus_seqs
        mv seq_context_TGCAT ./consensus_seq_context/
        mv seq_context_CATG ./consensus_seq_context/
      fi
    fi
		find . -type f -empty -delete
	fi
	wait
fi

if [[ "$p1" ]]; then
	if [ -d "${projdir}/snpfilter/2x" ]; then
		cd ${projdir}/snpfilter
		cp -r 2x 2x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./2x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_2x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_2x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$pseg" "$haplome_number"
		wait
		rm ${pop}_2x_rawRD${minRD_2x}_DP_GT.txt ${pop}_2x_DP_GT.txt ${pop}_2x_rd${minRD_2x}.txt 2> /dev/null
		mkdir visualizations; mv *.tiff ./visualizations/

    if [[ "$variant_intervals" == true ]]; then
      mkdir -p variant_intervals
      cd variant_intervals
      while IFS="" read -r p || [ -n "$p" ]; do
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
        for vbreak in vbreak_*; do
          awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
          awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt
          wait
        done
        sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
        mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
        wait && rm vbreak_*
        wait
      done < <(awk 'NR>1{print $2}' ../${pop}_2x_rd${minRD_2x}_noSDdose.txt | sort | uniq)
      cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
      rm variant_intervals_*.txt && cd ../
      rm -rf variant_intervals
    fi

    # Extract sequence context of variants
    if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        awk '{print $1"\t"$2"\t"$3}' ${pop}_2x_rd${minRD_2x}_noSDdose.txt > snplist.txt &&
        awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
        awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
        for wint in {100,200}; do
          awk -v wint=$wint '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
          awk -v wint=$wint '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
          wait
        done
        awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
        mv snplist_haps.tmp snplist_haps.txt &&
        rm snplist_round.txt
      fi
      if [[ ! -z "$RE1" ]]; then
        mkdir -p seq_context_${RE1}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
          rm ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
        done
      fi
      if [[ ! -z "$RE2" ]]; then
        mkdir -p seq_context_${RE2}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
          rm ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
        done
      fi
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        mkdir consensus_seq_context
        for i in $(cat ${projdir}/samples_list_node_*); do
          cat ./seq_context_${RE1}/${i%.f*}* ./seq_context_${RE2}/${i%.f*}* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
          maxlen=$(awk '{print $2}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}')
          awk -v maxlen=$maxlen '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          mv ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta
          for aln in $(awk '{print $1}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | sort | uniq); do
            awk -v aln=$aln '$1 == aln {print}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
            bash $cons -sequence ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp -outseq ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp
            cat ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp | awk -v aln="${i%.f*}_${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' >> ./consensus_seq_context/${i%.f*}_seqcontext.cons
            rm ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          done
        done
        for nn in $(ls ./consensus_seq_context/*cons);do
          for run in {1..10}; do
            awk '{gsub(/n$|N$/,"");}1' $nn > ${nn}.tmp && mv ${nn}.tmp ${nn}
          done
        done
        seqid=$(cat ./consensus_seq_context/*cons | grep '>' | awk '{gsub(/>/,""); gsub(/_Chr/,"\t_Chr");}1' | awk '{print $2}' | sort | uniq)
        mkdir -p ./consensus_seq_context/seqid_combine_samples
        for splitseqid in $seqid; do
          filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}"
          cat ./consensus_seq_context/*cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
          awk -v pat=$splitseqid '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/${filename}
          minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{print length}' | sort -n | head -n 1)
          awk -v minlen=$minlen '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{gsub(/\t/,"\n");}1' > ./consensus_seq_context/seqid_combine_samples/${filename}.tmp
          mv ./consensus_seq_context/seqid_combine_samples/${filename}.tmp ./consensus_seq_context/seqid_combine_samples/${filename}
        done
        mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/
        mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*cons ./consensus_seq_context/sample_consensus_seqs
        mv seq_context_TGCAT ./consensus_seq_context/
        mv seq_context_CATG ./consensus_seq_context/
      fi
    fi
		find . -type f -empty -delete
	fi
	wait
	if [ -d "${projdir}/snpfilter/4x" ]; then
		cd ${projdir}/snpfilter
		cp -r 4x 4x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./4x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_4x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_4x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$pseg" "$haplome_number"
		wait
		rm ${pop}_4x_rawRD${minRD_4x}_DP_GT.txt ${pop}_4x_DP_GT.txt ${pop}_4x_rd${minRD_4x}.txt 2> /dev/null
		mkdir visualizations; mv *.tiff ./visualizations/

    if [[ "$variant_intervals" == true ]]; then
      mkdir -p variant_intervals
      cd variant_intervals
      while IFS="" read -r p || [ -n "$p" ]; do
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
        for vbreak in vbreak_*; do
          awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
          awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt
          wait
        done
        sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
        mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
        wait && rm vbreak_*
        wait
      done < <(awk 'NR>1{print $2}' ../${pop}_4x_rd${minRD_4x}_noSDdose.txt | sort | uniq)
      cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
      rm variant_intervals_*.txt && cd ../
      rm -rf variant_intervals
    fi

    # Extract sequence context of variants
    if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        awk '{print $1"\t"$2"\t"$3}' ${pop}_4x_rd${minRD_4x}_noSDdose.txt > snplist.txt &&
        awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
        awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
        for wint in {100,200}; do
          awk -v wint=$wint '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
          awk -v wint=$wint '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
          wait
        done
        awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
        mv snplist_haps.tmp snplist_haps.txt &&
        rm snplist_round.txt
      fi
      if [[ ! -z "$RE1" ]]; then
        mkdir -p seq_context_${RE1}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
          rm ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
        done
      fi
      if [[ ! -z "$RE2" ]]; then
        mkdir -p seq_context_${RE2}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
          rm ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
        done
      fi
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        mkdir consensus_seq_context
        for i in $(cat ${projdir}/samples_list_node_*); do
          cat ./seq_context_${RE1}/${i%.f*}* ./seq_context_${RE2}/${i%.f*}* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
          maxlen=$(awk '{print $2}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}')
          awk -v maxlen=$maxlen '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          mv ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta
          for aln in $(awk '{print $1}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | sort | uniq); do
            awk -v aln=$aln '$1 == aln {print}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
            bash $cons -sequence ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp -outseq ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp
            cat ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp | awk -v aln="${i%.f*}_${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' >> ./consensus_seq_context/${i%.f*}_seqcontext.cons
            rm ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          done
        done
        for nn in $(ls ./consensus_seq_context/*cons);do
          for run in {1..10}; do
            awk '{gsub(/n$|N$/,"");}1' $nn > ${nn}.tmp && mv ${nn}.tmp ${nn}
          done
        done
        seqid=$(cat ./consensus_seq_context/*cons | grep '>' | awk '{gsub(/>/,""); gsub(/_Chr/,"\t_Chr");}1' | awk '{print $2}' | sort | uniq)
        mkdir -p ./consensus_seq_context/seqid_combine_samples
        for splitseqid in $seqid; do
          filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}"
          cat ./consensus_seq_context/*cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
          awk -v pat=$splitseqid '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/${filename}
          minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{print length}' | sort -n | head -n 1)
          awk -v minlen=$minlen '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{gsub(/\t/,"\n");}1' > ./consensus_seq_context/seqid_combine_samples/${filename}.tmp
          mv ./consensus_seq_context/seqid_combine_samples/${filename}.tmp ./consensus_seq_context/seqid_combine_samples/${filename}
        done
        mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/
        mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*cons ./consensus_seq_context/sample_consensus_seqs
        mv seq_context_TGCAT ./consensus_seq_context/
        mv seq_context_CATG ./consensus_seq_context/
      fi
    fi
		find . -type f -empty -delete
	fi
	wait
	if [ -d "${projdir}/snpfilter/6x" ]; then
		cd ${projdir}/snpfilter
		cp -r 6x 6x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./6x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_6x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_6x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$pseg" "$haplome_number"
		wait
		rm ${pop}_6x_rawRD${minRD_6x}_DP_GT.txt ${pop}_6x_DP_GT.txt ${pop}_6x_rd${minRD_6x}.txt 2> /dev/null
		mkdir visualizations; mv *.tiff ./visualizations/

    if [[ "$variant_intervals" == true ]]; then
      mkdir -p variant_intervals
      cd variant_intervals
      while IFS="" read -r p || [ -n "$p" ]; do
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
          awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt
          wait
        done
        sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
        mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
        wait && rm vbreak_*
        wait
      done < <(awk 'NR>1{print $2}' ../${pop}_6x_rd${minRD_6x}_noSDdose.txt | sort | uniq)
      cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
      rm variant_intervals_*.txt && cd ../
      rm -rf variant_intervals
    fi

    # Extract sequence context of variants
    if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        awk '{print $1"\t"$2"\t"$3}' ${pop}_6x_rd${minRD_6x}_noSDdose.txt > snplist.txt &&
        awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
        awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
        for wint in {100,200}; do
          awk -v wint=$wint '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
          awk -v wint=$wint '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
          wait
        done
        awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
        mv snplist_haps.tmp snplist_haps.txt &&
        rm snplist_round.txt
      fi
      if [[ ! -z "$RE1" ]]; then
        mkdir -p seq_context_${RE1}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
          rm ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
        done
      fi
      if [[ ! -z "$RE2" ]]; then
        mkdir -p seq_context_${RE2}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
          rm ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
        done
      fi
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        mkdir consensus_seq_context
        for i in $(cat ${projdir}/samples_list_node_*); do
          cat ./seq_context_${RE1}/${i%.f*}* ./seq_context_${RE2}/${i%.f*}* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
          maxlen=$(awk '{print $2}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}')
          awk -v maxlen=$maxlen '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          mv ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta
          for aln in $(awk '{print $1}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | sort | uniq); do
            awk -v aln=$aln '$1 == aln {print}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
            bash $cons -sequence ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp -outseq ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp
            cat ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp | awk -v aln="${i%.f*}_${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' >> ./consensus_seq_context/${i%.f*}_seqcontext.cons
            rm ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          done
        done
        for nn in $(ls ./consensus_seq_context/*cons);do
          for run in {1..10}; do
            awk '{gsub(/n$|N$/,"");}1' $nn > ${nn}.tmp && mv ${nn}.tmp ${nn}
          done
        done
        seqid=$(cat ./consensus_seq_context/*cons | grep '>' | awk '{gsub(/>/,""); gsub(/_Chr/,"\t_Chr");}1' | awk '{print $2}' | sort | uniq)
        mkdir -p ./consensus_seq_context/seqid_combine_samples
        for splitseqid in $seqid; do
          filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}"
          cat ./consensus_seq_context/*cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
          awk -v pat=$splitseqid '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/${filename}
          minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{print length}' | sort -n | head -n 1)
          awk -v minlen=$minlen '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{gsub(/\t/,"\n");}1' > ./consensus_seq_context/seqid_combine_samples/${filename}.tmp
          mv ./consensus_seq_context/seqid_combine_samples/${filename}.tmp ./consensus_seq_context/seqid_combine_samples/${filename}
        done
        mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/
        mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*cons ./consensus_seq_context/sample_consensus_seqs
        mv seq_context_TGCAT ./consensus_seq_context/
        mv seq_context_CATG ./consensus_seq_context/
      fi
    fi
		find . -type f -empty -delete
	fi
	wait
	if [ -d "${projdir}/snpfilter/8x" ]; then
		cd ${projdir}/snpfilter
		cp -r 8x 8x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./8x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_8x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_8x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$pseg" "$haplome_number"
		wait
		rm ${pop}_8x_rawRD${minRD_8x}_DP_GT.txt ${pop}_8x_DP_GT.txt ${pop}_8x_rd${minRD_8x}.txt 2> /dev/null
		mkdir visualizations; mv *.tiff ./visualizations/

    if [[ "$variant_intervals" == true ]]; then
      mkdir -p variant_intervals
      cd variant_intervals
      while IFS="" read -r p || [ -n "$p" ]; do
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
        for vbreak in vbreak_*; do
          awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
          awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt
          wait
        done
        sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
        mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
        wait && rm vbreak_*
        wait
      done < <(awk 'NR>1{print $2}' ../${pop}_8x_rd${minRD_8x}_noSDdose.txt | sort | uniq)
      cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
      rm variant_intervals_*.txt && cd ../
      rm -rf variant_intervals
    fi

    # Extract sequence context of variants
    if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        awk '{print $1"\t"$2"\t"$3}' ${pop}_8x_rd${minRD_8x}_noSDdose.txt > snplist.txt &&
        awk '{print $2"\t"$3}' snplist.txt | awk '{$2=sprintf("%d00",$2/100)}1' > snplist_round.txt &&
        awk '{$1":"$2}' snplist_round.txt > snplist_haps.txt &&
        for wint in {100,200}; do
          awk -v wint=$wint '{print $1":"$2-wint}' snplist_round.txt >> snplist_haps.txt &&
          awk -v wint=$wint '{print $1":"$2+wint}' snplist_round.txt >> snplist_haps.txt &&
          wait
        done
        awk '!seen[$0] {print} {++seen[$0]}' snplist_haps.txt | grep -v '^CHROM' > snplist_haps.tmp &&
        mv snplist_haps.tmp snplist_haps.txt &&
        rm snplist_round.txt
      fi
      if [[ ! -z "$RE1" ]]; then
        mkdir -p seq_context_${RE1}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
          rm ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
        done
      fi
      if [[ ! -z "$RE2" ]]; then
        mkdir -p seq_context_${RE2}
        for i in $(cat ${projdir}/samples_list_node_*); do
          $samtools view ../../preprocess/processed/${i%.f*}_${ref1%.f*}_precall.bam | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
          awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' snplist_haps.txt - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
          rm ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
        done
      fi
      if [[ ! -z "$RE1" ]] || [[ ! -z "$RE2" ]]; then
        mkdir consensus_seq_context
        for i in $(cat ${projdir}/samples_list_node_*); do
          cat ./seq_context_${RE1}/${i%.f*}* ./seq_context_${RE2}/${i%.f*}* | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta &&
          maxlen=$(awk '{print $2}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print length}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.5 - 0.5)]}')
          awk -v maxlen=$maxlen '{print $1"\t"substr($2, 1, maxlen)}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '!seen[$0]++' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          mv ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta
          for aln in $(awk '{print $1}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | sort | uniq); do
            awk -v aln=$aln '$1 == aln {print}' ./consensus_seq_context/${i%.f*}_seqcontext.fasta | awk '{print ">seq-"NR"_"$1"\n"$2}' > ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
            bash $cons -sequence ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp -outseq ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp
            cat ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp | awk -v aln="${i%.f*}_${aln}_locus" '{gsub(/EMBOSS_001/,aln);}1' >> ./consensus_seq_context/${i%.f*}_seqcontext.cons
            rm ./consensus_seq_context/${i%.f*}_seqcontext.cons.tmp ./consensus_seq_context/${i%.f*}_seqcontext.fasta.tmp
          done
        done
        for nn in $(ls ./consensus_seq_context/*cons);do
          for run in {1..10}; do
            awk '{gsub(/n$|N$/,"");}1' $nn > ${nn}.tmp && mv ${nn}.tmp ${nn}
          done
        done
        seqid=$(cat ./consensus_seq_context/*cons | grep '>' | awk '{gsub(/>/,""); gsub(/_Chr/,"\t_Chr");}1' | awk '{print $2}' | sort | uniq)
        mkdir -p ./consensus_seq_context/seqid_combine_samples
        for splitseqid in $seqid; do
          filename=${splitseqid%_locus} && filename=${filename#_*}.fasta && filename="${filename/:/_}"
          cat ./consensus_seq_context/*cons | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
          awk -v pat=$splitseqid '$1 ~ pat' > ./consensus_seq_context/seqid_combine_samples/${filename}
          minlen=$(awk '{print $2}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{print length}' | sort -n | head -n 1)
          awk -v minlen=$minlen '{print $1"\t"substr($2,1,minlen)}' ./consensus_seq_context/seqid_combine_samples/${filename} | awk '{gsub(/\t/,"\n");}1' > ./consensus_seq_context/seqid_combine_samples/${filename}.tmp
          mv ./consensus_seq_context/seqid_combine_samples/${filename}.tmp ./consensus_seq_context/seqid_combine_samples/${filename}
        done
        mkdir -p ./consensus_seq_context/sequences && mv ./consensus_seq_context/*fasta ./consensus_seq_context/sequences/
        mkdir -p ./consensus_seq_context/sample_consensus_seqs && mv ./consensus_seq_context/*cons ./consensus_seq_context/sample_consensus_seqs
        mv seq_context_TGCAT ./consensus_seq_context/
        mv seq_context_CATG ./consensus_seq_context/
      fi
    fi
		find . -type f -empty -delete
	fi
fi
done
wait
done


cd ${projdir}/snpfilter/
find . -type f -empty -delete
find . -type d -empty -delete
for snpfilter_dir in $(ls -d */); do
	cd $snpfilter_dir
	smmiss_thresh=${snpfilter_dir#*smiss}
	smmiss_thresh=${smmiss_thresh%*/}
	smmiss_thresh=$((${smmiss_thresh#*.}*10))
	awk -v smisst=$smmiss_thresh '(NR>1) && ($2 <= smisst)' sample_missing_rate* > retained_samples.txt 2> /dev/null
	cd ../
done
wait
wc -l *gmiss*/*dose.txt | awk '{print $2"\t"$1-1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > gmiss_smiss_titration.txt
wc -l *gmiss*/eliminated* | awk '{print $2"\t"$1-1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > eliminated_samples.txt
wc -l *gmiss*/retained_samples.txt | awk '{print $2"\t"$1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > retained_samples.txt
echo -e "gmiss_smiss_thresholds\t#_of_retained_samples\t#_of_SNPs\t#_of_eliminated_samples\n----------------------\t-----------------------\t---------\t-----------------------" > summary_precall.txt
awk 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,"\t",a[$1],"\t",$2}' gmiss_smiss_titration.txt eliminated_samples.txt | \
awk 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,"\t",a[$1],"\t",$2,"\t",$3}' retained_samples.txt - |
cat summary_precall.txt - > gmiss_smiss.txt
rm gmiss_smiss_titration.txt eliminated_samples.txt retained_samples.txt summary_precall.txt
# wc -l *gmiss*/unique_mapped/*dose*.txt | awk '{print $2"\t"$1-1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > gmiss_smiss_titration.txt
# wc -l *gmiss*/eliminated* | awk '{print $2"\t"$1-1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > eliminated_samples.txt
# wc -l *gmiss*/retained_samples.txt | awk '{print $2"\t"$1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > retained_samples.txt
# echo -e "gmiss_smiss_thresholds\t#_of_retained_samples\t#_of_SNPs\t#_of_eliminated_samples\n----------------------\t-----------------------\t---------\t-----------------------" > summary_precall.txt
# awk 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,"\t",a[$1],"\t",$2}' gmiss_smiss_titration.txt eliminated_samples.txt | \
# awk 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,"\t",a[$1],"\t",$2,"\t",$3}' retained_samples.txt - |
# cat summary_precall.txt - > gmiss_smiss_unique_mapped.txt
# rm gmiss_smiss_titration.txt eliminated_samples.txt retained_samples.txt summary_precall.txt
wait
for snpfilter_dir in $(ls -d */); do
if [[ "$(ls ${snpfilter_dir}/*dose.txt 2> /dev/null | wc -l)" -eq 0 ]]; then
rm -rf ${snpfilter_dir}
fi
done
wait
ls ./*/*maf*.txt 2> /dev/null | grep -v 'maf0.txt' | grep -v 'dose' | grep -v 'binary' | xargs rm
ls ./*/*_plusSD.txt 2> /dev/null | xargs rm
ls ./*/*SD_1_G*G*.txt 2> /dev/null | xargs rm



cd "$projdir"/snpfilter
n="${ref1%.f*}_"
for snpfilter_dir in $(ls -d */); do
	if [ -d "$snpfilter_dir" ]; then
		cd $snpfilter_dir
		# mkdir -p unique_mapped

		for i in $(ls *dose.txt 2> /dev/null); do
			ARselect=${i%rd*}
			ARfile=$(ls ../../snpcall/${ARselect}*AR.txt 2> /dev/null)
      arr=$(cat ${projdir}/samples_list_node_* | awk '{gsub(/.f/,"\t.f");}1' | awk '{print $1}' | tr '\n' ',' | awk '{gsub(/,$/,"");}1')
      arr2=$(grep "CHROM" $i | awk '{$1=$2=$3=$4=$5=""}1' | tr -s ' ' | awk '{gsub(/ pvalue/,"");}1' | awk '{gsub(/\t/,",");gsub(/ /,",");gsub(/^,/,"");}1')
      darr=$(echo ${arr[@]},${arr2[@]} | tr ',' '\n' | sort | uniq -u | tr '\n' ',' | awk '{gsub(/,$/,"");}1')

			Rscript "${GBSapp_dir}"/scripts/R/heterozygote_vs_allele_ratio.R "$i" "$ARfile" "${ploidy}x" "1" "$darr" "${GBSapp_dir}/tools/R"
			awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$2,$3]){print b}}' $ARfile $i | awk '{gsub(/NA/,"na"); print $1"_"$2"\t"$0}' | \
			awk -v n="$n" '{gsub(n,""); gsub(/CHROM_POS/,"SNP");}1' > ${i%.txt}_AR_metric.txt

			vcfdose=${i%_rd*}; vcfdose=${vcfdose#*_}
			zcat ../../snpcall/*${vcfdose}.vcf.gz | grep '^#' > ${i%.txt}.vcf
			awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$2,$3]){print b}}' <(zcat ../../snpcall/*${vcfdose}.vcf.gz) $i >> ${i%.txt}.vcf
			$bcftools view -s "$arr" ${i%.txt}.vcf > tmp.vcf && mv tmp.vcf ${i%.txt}.vcf

			grep -v '^##' ${i%.txt}.vcf | awk '{gsub(/#CHROM/,"CHROM");}1' > ${i%.txt}_tmp.vcf
			Rscript "${GBSapp_dir}"/scripts/R/recode_vcf.R "${i%.txt}_tmp.vcf" "$i" "${i%.txt}_AR_metric.txt" "${ploidy}x" "$darr" "${GBSapp_dir}/tools/R"
			grep '^##' ${i%.txt}.vcf | cat - <(awk '{gsub(/CHROM/,"#CHROM");}1' dose_temp.vcf) > ${i%.txt}.vcf
			mv AR_temp.txt ${i%.txt}_AR_metric.txt
			rm ${i%.txt}_tmp.vcf dose_temp.vcf
			gzip ${i%.txt}.vcf

			grep -v 'CHROM' ${i%.txt}_AR_metric.txt | awk -F'\t' '{$1=$2=$3=$4=$5=""}1' | awk -F'\t' '{gsub(/na/,"0");}1' | \
			awk '{gsub(/\t/," ");}1' | awk '{gsub(/-/,"");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/0,/,",");}1' | awk '{gsub(/,0$/,",");}1' | \
			awk -F',' -v OFS=',' -v OFMT='%0.3g' '{s=0; numFields=0; for(i=2; i<=NF;i++){if(length($i)){s+=$i; numFields++}} print (numFields ? s/numFields : 0)}' | \
			cat <(printf "Allele_ratio_mean\n") - | paste <(awk '{print $1"\t"$2"\t"$3}' ${i%.txt}_AR_metric.txt) - > ${i%.txt}_AR_mean.txt

			Rscript "${GBSapp_dir}"/scripts/R/hapmap_format.R "$i" "${GBSapp_dir}/tools/R"
      mv outfile.hmp.txt "${i%.txt}.hmp.txt"
		done
		wait

		for i in $(ls *dose* 2> /dev/null); do
			awk -v n="$n" '{gsub(n,""); print $0}' $i > ${i%.txt}_hold.txt
			mv ${i%.txt}_hold.txt $i
		done
    wait


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
		# 	$bcftools view -s "$arr" ${i%.txt}.vcf > tmp.vcf && mv tmp.vcf ${i%.txt}.vcf
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
		# 	$bcftools view -s "$arr" ${i%.txt}.vcf > tmp.vcf && mv tmp.vcf ${i%.txt}.vcf
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

		cd "$projdir"/snpfilter
	fi
done


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
for snpfilter_dir in $(ls -d */); do
	cd "$snpfilter_dir"
	dose=${snpfilter_dir%_gmiss*} && \
	dose=${dose%_*} && \
	awk 'BEGIN{OFS="\t"}{print $2,$3}' *dose_RefPrefix.txt > CHROM_POS.txt && \
	cat "$projdir"/snpcall/*"${dose}"*.vcf | grep -Fwf CHROM_POS.txt - | awk '{gsub(/#/,""); print $0}' > snp_allele_depth.txt && \
	cd "$projdir"/snpfilter
done
wait
for snpfilter_dir in $(ls -d */); do (
	cd "$snpfilter_dir"
	dose=${snpfilter_dir%_gmiss*} && \
	dose=${dose%_*} && \
	poptype=${snpfilter_dir%_*_*} && \
	poptype=${poptype//*_} && \
	getrd=$( echo *dose.txt | awk  'BEGIN{OFS=FS="_rd"};{print $2}' | awk  'BEGIN{OFS=FS="_"};{print $1}' ) && \
	if [[ $poptype == biparental ]]; then
		cd ../"$snpfilter_dir"
		Rscript "${GBSapp_dir}"/scripts/R/SNPaccuracy_biparental_ReadDepth.R $p1 $p2 $getrd "${GBSapp_dir}"/tools/R
		wait
	fi
	if [[ $poptype == diversity ]]; then
		cd ../"$snpfilter_dir"
		Rscript "${GBSapp_dir}"/scripts/R/SNPaccuracy_diversity_ReadDepth.R $getrd "${GBSapp_dir}"/tools/R
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
for snpfilter_dir in $(ls -d */); do
	cd "$snpfilter_dir"  && \
	cd genotype_accuracy  && \
	i=0  && \
	for f in $(ls *); do
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
for snpfilter_dir in $(ls -d */); do
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
	for sample in $(ls -S *_precall.bam); do
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
	for sample in $(ls -S ./subsamples/*.fasta.gz); do
		$bwa mem -t $threads snpseq_context.fasta $sample | $samtools view -F 4 > ${sample%.fasta.gz}_align.txt
		awk '{print $1"\t"$3}' ${sample%.fasta.gz}_align.txt | awk 'BEGIN{OFS="\t"} {gsub(/^.*-/,"",$1); print}' > ${sample%.fasta.gz}_alignh.txt
		$samtools bam2fq ${sample%.fasta.gz}_align.txt | awk 'NR%2==0' | awk 'NR%2==1' > ${sample%.fasta.gz}_aligns.txt
		paste -d '\t' ${sample%.fasta.gz}_alignh.txt ${sample%.fasta.gz}_aligns.txt > ${sample%.fasta.gz}_align.txt
		rm $sample ${sample%.fasta.gz}_alignh.txt ${sample%.fasta.gz}_aligns.txt
	done
	rm snpseq_context.sh snpseq_context*
	rm -rf subref

	for sample in $(ls -S ./subsamples/*_align.txt); do
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
					compress=$(echo $line | tr '\t' '\n' | sort | uniq | awk 'NF > 0' | tr '\n' '\t')
					printf '..%s..' "$compress\n" >> consensus_haplos2.txt
				done
				cat consensus_haplos2.txt > ./paralog_haplo_filter/consensus_haplos.txt
				rm consensus_haplos1.txt && rm consensus_haplos2.txt
			done
		fi
	fi
	cat ./paralog_haplo_filter/consensus_haplos.txt | while read -r line; do
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
			for i in $(ls -S ./paralog_haplo_filter/*_hapmissingSNP.txt); do
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
			for i in $(ls -S ./paralog_haplo_filter/*_hapmissingSNP.txt); do
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
			for i in $(ls -S ./paralog_haplo_filter/*_hapmissingSNP.txt); do
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
			for i in $(ls -S ./paralog_haplo_filter/*_hapmissingSNP.txt); do
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
	find ../ -size 0 -delete >/dev/null 2>&1
	touch Analysis_Complete
	rm *node*.txt
	mv ${projdir}/GBSapp_run_node_1.sh ${projdir}/GBSapp_run_node_1_done.sh 2> /dev/null
	mv ${projdir}/GBSapp_run_node.sh ${projdir}/GBSapp_run_node_done.sh 2> /dev/null
else
	touch Analysis_Complete_${samples_list}
fi
wait
echo -e "${magenta}- Run Complete. ${white}\n"
