
if [ -z "$threads" ]; then
	threads=$(nproc --all)
	if [[ "$threads" -ge 4 ]]; then
		threads=$((threads-2))
	fi
fi
if [ -z "$nodes" ]; then
 nodes=1
fi
if [[ -z $maxindel ]]; then
	maxindel=100
fi
if [[ -z $PEdist ]]; then
	PEdist=250
fi
if [ -z "$maxHaplotype" ]; then
	maxHaplotype=128
fi
if [ -z "$haplome_number" ]; then
	haplome_number=1
fi
if [ -z "$p2" ]; then
	p2=$p1
fi
if [ -z "$mhap_freq" ]; then
	mhap_freq=1
fi
if [ -z "$softclip" ]; then
	softclip=true
fi
if [ -z "$ncohorts" ]; then
	ncohorts=1
fi
if [ -z "$exit_before_jointcall" ]; then
	exit_before_jointcall=false
fi
if [ -z "$keep_gVCF" ]; then
	keep_gVCF=false
fi


echo -e "${blue}\n############################################################################## ${yellow}\n- Index Reference Genome \n${blue}##############################################################################${white}\n"
main () {
cd $projdir
cd refgenomes
for i in *.gz; do
	$gunzip $i >/dev/null 2>&1
done

if [[ -d "ref" ]]; then
	:
else
	checknfastafiles=$(ls *.f* | wc -l)
	if [[ $checknfastafiles -gt 1 ]]; then
		echo -e "${magenta}- expecting only 1 fasta file for reference genome ${white}\n"
		echo -e "${magenta}- GBSapp will quit in 5 seconds ${white}\n"
		sleep 5 && exit 1
	fi
	if [ -z "$ref1" ]; then
		for ref in $(ls *.f*); do
			ref1=${ref%.fa*}.fasta
		done
	fi
	export ncontigscaffold=$(grep '>' $ref1 | wc -l)
	if [[ $ncontigscaffold -gt 300 ]]; then
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
	if [[ $ncontigscaffold -gt 300 ]]; then
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
			for i in *.f*; do
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
if [[ -d "ref" ]]; then
	echo -e "${magenta}- indexed genome available ${white}\n"
	if [ -z "$ref1" ]; then
		for ref in $(ls *.f*); do
			ref1=${ref%%.f*}.fasta
		done
	fi
else
	echo -e "${magenta}- indexing single reference subgenome ${white}\n"
	if [ -z "$ref1" ]; then
		for ref in $(ls *.f*); do
			ref1=${ref%%.f*}.fasta
		done
	fi
	awk '{ sub("\r$",""); print}' $ref1 | awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' > ref.txt
	n=">${ref1%.f*}_"
	awk '{ sub("\r$",""); print}' ref.txt | awk -v n="$n" '{gsub(n,">"); print}' | awk -v n="$n" '{gsub(/>/,n); print}' > $ref1
	rm ref.txt
	$samtools faidx $ref1
	$java -jar $picard CreateSequenceDictionary REFERENCE= $ref1 OUTPUT=${ref1%.f*}.dict
	$bbmap ref=$ref1 k=4 threads=$threads
	$java -ea $Xmx2 -cp ${GBSapp_dir}/tools/bbmap/current/ align2.BBMap ref=$ref1 k=4 threads=$threads

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
			for i in $(ls *.f* | grep -v R2.f); do
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
			for i in $(ls *.f* ); do
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
			for i in $(ls *R1.f*); do
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
			for i in $(ls *R1.f*); do
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
}
cd $projdir
cd samples
if [ -d "pe" ]; then
	fqpass=$(find ./pe -maxdepth 1 -name '*_R1.f*' -o -name '*_R2.f*' -o -name '*.R1.f*' -o -name '*.R2.f*' | wc -l)
	fqfail=$(ls ./pe/* | wc -l)
	fqfail=$((fqfail-fqpass))
	if [[ "$fqfail" -lt 1 ]]; then
		if [[ "$samples_list" == "samples_list_node_1.txt" ]] && test ! -f ${projdir}/alignment_summaries/total_read_count.txt; then
			time main &>> ${projdir}/log.out
		fi
	else
		echo -e "${magenta}- samples' PE fastq filenames requires formatting (i.e. needs to end in "_R1.fastq" or ".R1.fastq" and "_R2.fastq" or ".R2.fastq") ${white}\n"
		sleep 5 && exit 1
	fi
else
	if [[ "$samples_list" == "samples_list_node_1.txt" ]] && test ! -f ${projdir}/alignment_summaries/total_read_count.txt; then
		time main &>> ${projdir}/log.out
	fi
fi


main () {
	cd ${projdir}
	mkdir -p preprocess
	mkdir -p snpcall
	mkdir -p alignment_summaries
	mkdir -p ./alignment_summaries/copy_number
	touch ${projdir}/alignment_summaries/total_read_count.txt
	printf 'sample\tnumber_of_reads\n' > ${projdir}/alignment_summaries/total_read_count.txt
	mkdir -p ${projdir}/alignment_summaries/background_mutation_test
	printf 'sample\trare_hap\tpop_hap\ttotal_hap\tmutation_load\n' > ${projdir}/alignment_summaries/pop_mutation_load.txt
}
if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
	time main &>> ${projdir}/log.out
fi


main () {
	cd $projdir
	cd samples
	nfiles=$(ls -1 -p | grep -v R2.f | grep -v / |  wc -l)
	totalk=$(awk '/^MemTotal:/{print $2}' /proc/meminfo)
	loopthreads=2
	if [[ "$threads" -gt 1 ]]; then
	  N=$((threads/2))
	  ram1=$(($totalk/$N))
	else
	  N=1 && loopthreads=threads
	fi
	ram1=$((ram1/1000000))
	Xmx1=-Xmx${ram1}G
	ram2=$(echo "$totalk*0.00000095" | bc)
	ram2=${ram2%.*}
	Xmx2=-Xmx${ram2}G
	if [[ "$nfiles" -lt "$N" ]]; then
	  N=$nfiles && loopthreads=$threads
	fi

	if [[ "$threads" -le 6 ]]; then
		prepthreads=threads
		Xmxp=$Xmx2
		prepN=1
	else
		prepthreads=6
		prepN=$(( threads / prepthreads ))
		ramprep=$(( ram2 / prepN ))
		Xmxp=-Xmx${ramprep}G
	fi

	if [[ "$threads" -le 4 ]]; then
		gthreads=threads
		Xmxg=$Xmx2
		gN=1
	else
		gthreads=4
		gN=$(( threads / gthreads ))
		ramg=$(( ram2 / gN ))
		Xmxg=-Xmx${ramg}G
	fi
}
cd $projdir
main &>> log.out

######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${yellow}\n- Performing Read Alignments & Alignment Post-Processing\n${blue}##############################################################################${white}\n"
main () {
	cd $projdir
	cd samples

	for i in $( cat ${projdir}/${samples_list} ); do (
		if test ! -f ${projdir}/compress_done.txt  && test ! -f ${projdir}/preprocess/${i%.f*}_redun.sam && test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai; then
			while test ! -f ${i%.f*}_uniq_R1.fasta; do
				sleep $[ ( $RANDOM % 30 )  + 10 ]s
				if [[ $(file $i | awk -F' ' '{print $2}') == gzip ]]; then
					$gunzip -c $i | awk 'NR%2==0' | awk 'NR%2' > ${i%.f*}_uniq.txt 2> /dev/null
					wait
					if test -f ${i%.f*}_R2*; then
						$gunzip -c ${i%.f*}_R2* | awk 'NR%2==0' | awk 'NR%2' > ${i%.f*}_R2_uniq.txt 2> /dev/null
						wait
					fi
					if test -f ${i%.f*}.R2*; then
						$gunzip -c ${i%.f*}.R2* | awk 'NR%2==0' | awk 'NR%2' > ${i%.f*}_R2_uniq.txt 2> /dev/null
						wait
					fi
				else
					awk 'NR%2==0' $i | awk 'NR%2' > ${i%.f*}_uniq.txt 2> /dev/null
					wait
					if test -f ${i%.f*}_R2*; then
						awk 'NR%2==0' ${i%.f*}_R2* | awk 'NR%2' > ${i%.f*}_R2_uniq.txt 2> /dev/null
						wait
					fi
					if test -f ${i%.f*}.R2*; then
						awk 'NR%2==0' ${i%.f*}.R2* | awk 'NR%2' > ${i%.f*}_R2_uniq.txt 2> /dev/null
						wait
					fi
				fi
				if test -f "${i%.f*}_R2_uniq.txt"; then
				:
				else
				touch "${i%.f*}_R2_uniq.txt"
				fi
				wait

				cat ${i%.f*}_uniq.txt | printf "${i%.f*}""\t""$(wc -l)""\n" > ${projdir}/alignment_summaries/${i%.f*}_total_read_count.txt  2> /dev/null &&
				export LC_ALL=C; paste -d ~ ${i%.f*}_uniq.txt ${i%.f*}_R2_uniq.txt | expand -t $(( $(wc -L < $i ) + 2 )) | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{$1=$1};1' | \
				awk '{gsub(" /"," "); print}' | awk '{gsub("/\n","\n"); print}' | awk '{gsub("/"," "); print}' | awk '{gsub(" ","\t"); print}' > ${i%.f*}_rdrefseq.txt 2> /dev/null &&
				awk 'NF==2 {print ">seq"NR"_se-"$1"\t"$2}' ${i%.f*}_rdrefseq.txt > ${i%.f*}_rdrefseq_se.txt 2> /dev/null &&
				awk 'NF==3 {print ">seq"NR"_pe-"$0}' ${i%.f*}_rdrefseq.txt | awk '{print $1"\t"$3}' > ${i%.f*}_uniq_R2.fasta 2> /dev/null &&
				awk 'NF==3 {print ">seq"NR"_pe-"$0}' ${i%.f*}_rdrefseq.txt | awk '{print $1"\t"$2}' | cat - ${i%.f*}_rdrefseq_se.txt > ${i%.f*}_uniq_R1.hold.fasta 2> /dev/null &&
				cat ${i%.f*}_uniq_R1.hold.fasta > ${projdir}/alignment_summaries/background_mutation_test/${i%.f*}_pop_haps.fasta 2> /dev/null &&
				rm ${i%.f*}*.txt 2> /dev/null
				find . -size 0 -delete  2> /dev/null
				mv ${i%.f*}_uniq_R1.hold.fasta ${i%.f*}_uniq_R1.fasta  2> /dev/null
			done
		fi ) &
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
			wait
		fi
	done
	wait

	if [[ "$(wc -l ${projdir}/alignment_summaries/total_read_count.txt | awk '{print $1}')" -le 1 ]]; then
		find ${projdir}/alignment_summaries/*_total_read_count.txt | xargs cat > ${projdir}/alignment_summaries/total_read_count.hold.txt &&
		cat ${projdir}/alignment_summaries/total_read_count.hold.txt >> ${projdir}/alignment_summaries/total_read_count.txt &&
		rm ${projdir}/alignment_summaries/*_total_read_count.txt ${projdir}/alignment_summaries/total_read_count.hold.txt
	fi
	wait

	if [[ ! -f "${projdir}/align1_${samples_list}" ]]; then touch "${projdir}/align1_${samples_list}"; fi
	if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
		align=$(ls ${projdir}/align1_samples_list_node_* | wc -l)
		while [[ "$align" -lt $nodes ]]; do sleep 300; align=$(ls ${projdir}/align1_samples_list_node_* | wc -l); done
		if [[ $align == $nodes ]] && test ! -f ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqPass.txt; then
			find ${projdir}/alignment_summaries/background_mutation_test/*_pop_haps.fasta | xargs cat > ${projdir}/alignment_summaries/background_mutation_test/pop_haps.txt &&
			rm ${projdir}/alignment_summaries/background_mutation_test/*_pop_haps.fasta ${projdir}/align1_${samples_list} &&
			awk -F "\t" 'BEGIN { OFS=FS }; { print $1, substr($2, 1, 64); }' ${projdir}/alignment_summaries/background_mutation_test/pop_haps.txt  | \
			awk '{a[$2]++} END{for(s in a){print a[s]" "s}}' | awk -F'\t' '{gsub(/ /,"\t"); print}' > ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freq.txt &&
			awk -F "\t" 'BEGIN { OFS=FS }; { print $1, substr($2, 1, 64); }' ${projdir}/alignment_summaries/background_mutation_test/pop_haps.txt | \
			awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$2]=$0;next} ($2) in a{print $0, a[$2]}' ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freq.txt - | \
			awk '{print $1"\t"$2"\t"$3}' > ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqall.txt &&
			rm ${projdir}/alignment_summaries/background_mutation_test/pop_haps.txt ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freq.txt &&
			awk -v phap=$mhap_freq '($3 == phap)' ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqall.txt > ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqFail.txt &&
			awk -v phap=$mhap_freq '($3 > phap)' ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqall.txt > ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqPass.txt &&
			rm ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqall.txt
		fi
	fi
	wait

	while [[ ! -f "${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqPass.txt" ]]; do sleep 300; done
	sleep 5
	cd ${projdir}/samples

	for i in $(cat ${projdir}/${samples_list} ); do (
		if test ! -f ${projdir}/hapfilter_done.txt && test ! -f ${projdir}/preprocess/${i%.f*}_redun.sam && test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai; then
			while test ! -f ${projdir}/alignment_summaries/${i%.f*}_pop_mutation_load.txt; do
				sleep $[ ( $RANDOM % 30 )  + 10 ]s
				export nempty=$( wc -l ${i%.f*}_uniq_R2.fasta &> /dev/null | awk '{print $1}' )
				cd ${projdir}/refgenomes/
				p_hap=$(awk '{print $1}' ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqPass.txt | awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' - ${projdir}/samples/${i%.f*}_uniq_R1.fasta | wc -l)
				r_hap=$(awk '{print $1}' ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqFail.txt | awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' - ${projdir}/samples/${i%.f*}_uniq_R1.fasta | wc -l)

				if [[ "$nempty" -gt 0 ]]; then
					awk '{print $1}' ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqPass.txt | \
					awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' - ${projdir}/samples/${i%.f*}_uniq_R1.fasta | \
					awk 'NF{NF-=1};1' > ${projdir}/samples/${i%.f*}_uniq_R1_singleton.txt &&
					grep '_se-' ${projdir}/samples/${i%.f*}_uniq_R1_singleton.txt | awk '{gsub(/>/,"@"); print}' | awk '{print $1"\t"$2"\t"$2}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | awk '{print $1"\n"$2"\n+\n"$3}' | $gzip > ${projdir}/samples/${i%.f*}_uniq_singleton.fq.gz &&
					grep '_pe-' ${projdir}/samples/${i%.f*}_uniq_R1_singleton.txt | awk '{gsub(/>/,"@"); print}' | awk '{print $1"\t"$2"\t"$2}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | awk '{print $1"/1\n"$2"\n+\n"$3}' | $gzip > ${projdir}/samples/${i%.f*}_uniq_R1.fq.gz &&
					rm ${projdir}/samples/${i%.f*}_uniq_R1_singleton.txt &&
					awk '{print $1}' ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqPass.txt | \
					awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' - ${projdir}/samples/${i%.f*}_uniq_R2.fasta | awk 'NF{NF-=1};1' | awk '{gsub(/>/,"@"); print}' | \
					awk '{print $1"\t"$2"\t"$2}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | awk '{print $1"/2\n"$2"\n+\n"$3}' | $gzip  > ${projdir}/samples/${i%.f*}_uniq_R2.fq.gz
					wait
				else
					awk '{print $1}' ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqPass.txt | \
					awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' - ${projdir}/samples/${i%.f*}_uniq_R1.fasta | awk 'NF{NF-=1};1' | awk '{gsub(/>/,"@"); print}' | \
					awk '{print $1"\t"$2"\t"$2}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | awk '{print $1"\n"$2"\n+\n"$3}' | $gzip  > ${projdir}/samples/${i%.f*}_uniq_R1.fq.gz
					wait
				fi

				thap=$(awk -v var1=$r_hap -v var2=$p_hap 'BEGIN { print  ( var1 + var2 ) }' )
				mload=$(awk -v var1=$r_hap -v var2=$thap 'BEGIN { print  ( var1 / var2 ) }' )
				printf "${i%.f*}\t${r_hap}\t${p_hap}\t${thap}\t${mload}\n" > ${projdir}/alignment_summaries/${i%.f*}_pop_mutation_load.hold.txt &&
				mv ${projdir}/alignment_summaries/${i%.f*}_pop_mutation_load.hold.txt ${projdir}/alignment_summaries/${i%.f*}_pop_mutation_load.txt &&
				cd ${projdir}/samples
			done
		fi ) &
		if [[ $(jobs -r -p | wc -l) -ge $((gN*2)) ]]; then
			wait
		fi
	done
	wait && touch ${projdir}/compress_done.txt


	for i in $(cat ${projdir}/${samples_list} ); do (
		cd ${projdir}/refgenomes
		if test ! -f ${projdir}/alignment_done.txt && test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai; then
			while test ! -f ${projdir}/preprocess/${i%.f*}_redun.sam; do
				if [[ "$nempty" -gt 0 ]]; then
					$java -ea $Xmxg -cp ${GBSapp_dir}/tools/bbmap/current/ align2.BBMap fast=t threads=$gthreads averagepairdist=$PEdist deterministic=t maxindel=$maxindel local=t keepnames=t maxsites=12 saa=f secondary=t ambiguous=all ref=$ref1 in1=${projdir}/samples/${i%.f*}_uniq_R1.fq.gz in2=${projdir}/samples/${i%.f*}_uniq_R2.fq.gz out=${projdir}/preprocess/${i%.f*}_redun_R1R2.sam &&
					$java -ea $Xmxg -cp ${GBSapp_dir}/tools/bbmap/current/ align2.BBMap fast=t threads=$gthreads averagepairdist=$PEdist deterministic=t maxindel=$maxindel local=t keepnames=t maxsites=12 saa=f secondary=t ambiguous=all ref=$ref1 in1=${projdir}/samples/${i%.f*}_uniq_singleton.fq.gz out=${projdir}/preprocess/${i%.f*}_redun_singleton.sam &&
					grep -v '^@' ${projdir}/preprocess/${i%.f*}_redun_singleton.sam | cat ${projdir}/preprocess/${i%.f*}_redun_R1R2.sam - > ${projdir}/preprocess/${i%.f*}_redun.hold.sam &&
					rm ${projdir}/preprocess/${i%.f*}_redun_singleton.sam ${projdir}/preprocess/${i%.f*}_redun_R1R2.sam
					wait
				else
					$java -ea $Xmxg -cp ${GBSapp_dir}/tools/bbmap/current/ align2.BBMap fast=t threads=$gthreads maxindel=$maxindel local=t keepnames=t maxsites=12 saa=f secondary=t ambiguous=all ref=$ref1 in1=${projdir}/samples/${i%.f*}_uniq_R1.fq.gz out=${projdir}/preprocess/${i%.f*}_redun.hold.sam
					wait
				fi
				rm ${projdir}/samples/${i%.f*}_uniq_*.fq.gz && mv ${projdir}/preprocess/${i%.f*}_redun.hold.sam ${projdir}/preprocess/${i%.f*}_redun.sam
				wait
				if [[ "$nempty" -gt 0 ]]; then
					rm ${projdir}/samples/${i%.f*}_uniq_R1.fasta ${projdir}/samples/${i%.f*}_uniq_R2.fasta
					wait
				else
					rm ${projdir}/samples/${i%.f*}_uniq_R1.fasta
					wait
				fi
			done
			wait
		fi
		cd ${projdir}/samples ) &
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
			wait
		fi
	done
	wait && touch ${projdir}/hapfilter_done.txt


	if [[ ! -f "${projdir}/align2_${samples_list}" ]]; then touch "${projdir}/align2_${samples_list}"; fi
	if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
		align=$(ls ${projdir}/align2_samples_list_node_* | wc -l)
		while [[ "$align" -lt $nodes ]]; do sleep 300; align=$(ls ${projdir}/align2_samples_list_node_* | wc -l); done
		if [[ $align == $nodes &&  "$(wc -l ${projdir}/alignment_summaries/pop_mutation_load.txt | awk '{print $1}')" -eq 1 ]]; then
			find ${projdir}/alignment_summaries/*_pop_mutation_load.txt | xargs cat > ${projdir}/alignment_summaries/pop_mutation_load.hold.txt &&
			cat ${projdir}/alignment_summaries/pop_mutation_load.hold.txt >> ${projdir}/alignment_summaries/pop_mutation_load.txt &&
			rm ${projdir}/alignment_summaries/pop_mutation_load.hold.txt &&
			rm ${projdir}/alignment_summaries/*_pop_mutation_load.txt ${projdir}/align2_${samples_list}
			wait
		fi
	fi

	align2_end=$( wc -l ${projdir}/alignment_summaries/pop_mutation_load.txt | awk '{print $1}')
	while [[ "$align2_end" -lt 2 ]]; do sleep 300; done
	sleep 60
	cd ${projdir}/samples


	for i in $(cat ${projdir}/${samples_list} ); do (
		if test ! -f ${projdir}/precall_done.txt; then
			printf '\n###---'${i%.f*}'---###\n' > ${projdir}/alignment_summaries/${i%.f*}_summ.txt && \
			$samtools flagstat ${projdir}/preprocess/${i%.f*}_redun.sam >> ${projdir}/alignment_summaries/${i%.f*}_summ.txt && \
			printf '########################################################################################################\n\n' >> ${projdir}/alignment_summaries/${i%.f*}_summ.txt &&
			printf 'copy\tFrequency\tPercentage\n' > ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_Unique_Read_histogram.txt && \
			grep -v '^@' ${projdir}/preprocess/${i%.f*}_redun.sam | awk -F' ' '{print $1}' | awk '{gsub(/_/,"\t"); print $2}' | \
			sort | uniq -c | awk '{print $2"\t"$1}' > ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number.txt  && \
			awk 'NR==FNR{sum+= $2; next;} {printf("%s\t%s\t%3.3f%%\t%3.0f\n",$1,$2,100*$2/sum,100*$2/sum)}' ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number.txt ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number.txt > ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot.txt && \
			unset IFS; printf "%s\t%s\t%s\t%*s\n" $(sed 's/$/ |/' ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot.txt) | tr ' ' '|' | sort -k2,2 -nr | awk '{gsub(/se-/,""); gsub(/pe-/,""); print}' >> ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_Unique_Read_histogram.txt && rm ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number.txt ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot.txt


			grep -v '^@' ${projdir}/preprocess/${i%.f*}_redun.sam | awk '($3 != "\*")' | awk '($6 != "\*")' | awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | \
			cat <(grep '^@' ${projdir}/preprocess/${i%.f*}_redun.sam) - > ${projdir}/preprocess/${i%.f*}_del.sam
			wait

			echo "nloci~mapQ~CHROM~POS" > ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}.txt
			grep -v '^@' ${projdir}/preprocess/${i%.f*}_del.sam | awk -F'\t' '{print $1"\t"$3"\t"$4"\t"$5}' | awk '{gsub(/ /,"\t"); print}' | awk '{print $1"~"$5"~"$3"~"$4}' | grep -v '*' >> ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}.txt &&
			awk '!visited[$0]++' ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}.txt > ${projdir}/alignment_summaries/temp_${i%.f*}.txt &&
			mv ${projdir}/alignment_summaries/temp_${i%.f*}.txt ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}.txt
			rm ${projdir}/preprocess/${i%.f*}_del.sam
			wait


			awk '/@HD/ || /@SQ/{print}' ${projdir}/preprocess/${i%.f*}_redun.sam > ${projdir}/preprocess/${i%.f*}_heading.sam
			grep -v '^@' ${projdir}/preprocess/${i%.f*}_redun.sam | awk '($3 != "\*")' | awk '{gsub(/_se-/,"_se-\t",$1); gsub(/_pe-/,"_pe-\t",$1)}1' | \
			awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' > ${projdir}/preprocess/${i%.f*}_uniq.sam &&
			for j in $(LC_ALL=C; sort -n -k1,1 ${projdir}/preprocess/${i%.f*}_uniq.sam | awk '{print $1}' | uniq); do
				awk -v n="^${j}" '$0~n{print $0}' ${projdir}/preprocess/${i%.f*}_uniq.sam | awk -v n="$j" '{for(i=0;i<n;i++) print}' >> ${projdir}/preprocess/${i%.f*}_exp.sam &&
				wait
			done
			awk '{print "seq"NR"_"$0}' ${projdir}/preprocess/${i%.f*}_exp.sam | tr -s ' ' | awk '{gsub(/256/,"0",$2);gsub(/272/,"16",$2); print}' | \
			awk '!($3 ~ "\*")' | awk '!($6 ~ "\*")' | cat ${projdir}/preprocess/${i%.f*}_heading.sam - | tr ' ' '\t' > ${projdir}/preprocess/${i%.f*}_${ref1%.f*}.sam &&
			rm ${projdir}/preprocess/${i%.f*}_exp.sam ${projdir}/preprocess/${i%.f*}_uniq.sam ${projdir}/preprocess/${i%.f*}_heading.sam

			j="${i%.f*}_${ref1%.f*}.sam"
			cd ${projdir}/preprocess
			$java $Xmx2 -XX:ParallelGCThreads=$gthreads -jar $picard SortSam I=$j O=${j%.sam*}.bam  SORT_ORDER=coordinate && \
			$java $Xmx2 -XX:ParallelGCThreads=$gthreads -jar $picard BuildBamIndex INPUT=${j%.sam*}.bam && \
			$java $Xmx2 -XX:ParallelGCThreads=$gthreads -jar $picard AddOrReplaceReadGroups I=${j%.sam*}.bam O=${j%.sam*}_precall.bam RGLB=${i%.f*} RGPL=illumina VALIDATION_STRINGENCY=LENIENT RGPU=run RGSM=${i%.f*} && \
			$samtools index ${j%.sam*}_precall.bam
			ls ${i%.f*}_* | grep -v precall | xargs rm

			cd ${projdir}/samples
		fi ) &
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
			wait
		fi
	done

	wait && touch ${projdir}/alignment_done.txt && touch ${projdir}/precall_done.txt

	if [[ ! -f "${projdir}/align3_${samples_list}" ]]; then touch "${projdir}/align3_${samples_list}"; fi
	if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
		align=$(ls ${projdir}/align3_samples_list_node_* | wc -l)
		while [[ "$align" -lt $nodes ]]; do sleep 300; align=$(ls ${projdir}/align3_samples_list_node_* | wc -l); done
		if [[ $align == $nodes ]] && test ! -f ${projdir}/alignment_summaries/refgenome_paralogs.txt; then
			rm ${projdir}/align3_${samples_list}
			cd ${projdir}/alignment_summaries
			cat *_summ.txt > alignment_summaries_unique_reads.txt; rm -r *_summ.txt
			# Total number of reads per samples
			awk '/###---/ || /QC-passed/{print}' alignment_summaries_unique_reads.txt | cut -d\+ -f1 | tr -d '\n' | \
			awk  'gsub(/---###/, "\t", $0)' | awk  'gsub(/###---/, "", $0)' | tr ' ' '\n' > total_unique_reads.txt
			# Total number of mapped reads per samples
			cat alignment_summaries_unique_reads.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /0_mapped/{print}' |\
			tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
			awk 'gsub("\\+0mapped", "\t", $0)' | cut -d\: -f1 > total_unique_reads_mapped.txt
			# Total number of mapped paired reads per samples
			cat alignment_summaries_unique_reads.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /properly_paired/{print}' |\
			tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
			awk 'gsub("\\+0properlypaired", "\t", $0)' | cut -d\: -f1 > total_unique_reads_paired.txt
			echo -e "Samples\tTotal\tMapped\tPerc_Mapped\tPE_Mapped\t%_PE_Mapped" > summary_precall.txt
			awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' total_unique_reads_mapped.txt  total_unique_reads.txt  | \
			awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' total_unique_reads_paired.txt - | awk '{gsub(/ /,"\t"); print $0}' | \
			cat summary_precall.txt - > Tabulated_Alignment_Unique_Read_Summaries.txt
			awk '{gsub(/\t/,","); print $0}' Tabulated_Alignment_Unique_Read_Summaries.txt > Tabulated_Alignment_Unique_Read_Summaries.csv
			rm total_unique_* summary_precall.txt &> /dev/null
			rm ${projdir}/samples/metrics.txt ${projdir}/preprocess/metrics.txt &> /dev/null

			cd $projdir/alignment_summaries

			touch refgenome_paralogs.txt
			for par in refgenome_paralogs_*.txt; do (
				cat refgenome_paralogs.txt $par | awk '!visited[$0]++' > temp_par.txt
				mv temp_par.txt refgenome_paralogs.txt ) &
			 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
				 wait
			 fi
			done
			wait
			rm refgenome_paralogs_*.txt
			awk '{gsub(/~/,"\t"); print $0}' refgenome_paralogs.txt | awk 'BEGIN{OFS="\t"; }; {if($2==0) $1 = "multilocus"; else $1 = $1; }; 1' | \
			awk 'BEGIN{OFS="\t"; };{print $3,$4,$1}' | awk '$3>max[$1,$2]{max[$1,$2]=$3; row[$1,$2]=$0} END{for (i in row) print row[i]}' | \
			awk '{print $0,substr($2, 1, length($2)-2)}' | awk '$3>max[$1,$4]{max[$1,$4]=$3; row[$1,$4]=$0} END{for (i in row) print row[i]}' > temp.txt
			awk '{print $1"\t"$2"\t"$3}' temp.txt | awk '!/CHROM/' | cat <(printf "CHROM\tPOS\tnloci\n") - > refgenome_paralogs.txt
			rm temp.txt
		fi
	fi
}
cd $projdir
while [[ ! -f "${projdir}/alignment_summaries/pop_mutation_load.txt" ]]; do sleep 300; done
sleep 60
if [ "$walkaway" == false ]; then
	echo -e "${magenta}- Do you want to perform read alignments and alignment post-processing? ${white}\n"
	read -p "- y(YES) or n(NO) " -t 36000 -n 1 -r
	if [[ ! $REPLY =~ ^[Yy]$ ]]; then
		printf '\n'
		echo -e "${magenta}- skipping read alignments and alignment post-processing ${white}\n"
	else
		printf '\n'
		echo -e "${magenta}- performing read alignments and alignment post-processing ${white}\n"
		time main &>> log.out
	fi
fi
if [ "$walkaway" == true ]; then
	if [ "$alignments" == 1 ]; then
		echo -e "${magenta}- performing read alignments and alignment post-processing ${white}\n"
		time main &>> log.out
	else
		echo -e "${magenta}- skipping read alignments and alignment post-processing ${white}\n"
	fi
fi


######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${yellow}\n- Performing Variant Calling with GATK HaplotypeCaller\n${blue}##############################################################################${white}\n"
main () {
cd $projdir
echo -e "${magenta}- performing SNP calling ${white}\n"
cd $projdir
cd preprocess
for i in $(cat ${projdir}/${samples_list} ); do (
	if test ! -f "${pop}_${ploidy}x_raw.vcf.gz" && test ! -f ${projdir}/GVCF_done.txt; then
		while test ! -f "${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf.gz"; do
				$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  HaplotypeCaller -R ${projdir}/refgenomes/$ref1 -I "${i%.f*}_${ref1%.f*}_precall.bam" -ploidy $ploidy -O "${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz" -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start 0 --minimum-mapping-quality 0 --max-num-haplotypes-in-population "$((ploidy * maxHaplotype))" && \
				mv "${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz" "${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf.gz"
				mv "${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz.tbi" "${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf.gz.tbi"
		done
	fi ) &
	if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		wait
	fi
done
wait

if [[ ! -f "${projdir}/call1_${samples_list}" ]]; then touch "${projdir}/call1_${samples_list}"; fi
if [[ "$samples_list" == "samples_list_node_1.txt" ]] && test ! -f ${projdir}/snpcall/${pop}_${ploidy}x_raw.vcf*; then
	align=$(ls ${projdir}/call1_samples_list_node_* | wc -l)
	while [[ "$align" -lt $nodes ]]; do sleep 300; align=$(ls ${projdir}/call1_samples_list_node_* | wc -l); done
	if [[ $align == $nodes ]]; then
		rm ${projdir}/call1_${samples_list}
		cd ${projdir}/snpcall
		if [[ "$ncohorts" -eq 1 ]]; then
			cz=$(ls *.g.vcf.gz | wc -l)
		fi
		if [[ "$ncohorts" -gt 1 ]]; then
			cz=$(ls *.g.vcf.gz | wc -l)
			cz=$(( (cz / ncohorts ) + (cz % ncohorts > 0)))
		fi

		i=0
		for f in `find . -maxdepth 1 -iname '*.g.vcf.gz' -type f | shuf`; do
			d=cohorts_$(printf %02d $((i/cz+1)))
			mkdir -p $d
			mv "$f" $d; mv "${f}.tbi" $d
			let i++
		done

		for dir in cohorts*/; do
			cd $dir
			j=--variant; input=""; k=""
			for i in $(ls *.g.vcf.gz); do
				k="${j} ${i}"; input="${input} ${k}"
			done
			Get_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/${ref1%.*}.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' | awk -v pat=${ref1%.f*} '$0 ~ pat')
			for selchr in $Get_Chromosome; do (
				while test ! -f ${pop}_${ploidy}x_"${selchr}"_raw.vcf.gz; do
					$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  GenomicsDBImport $input -L $selchr --genomicsdb-workspace-path ${pop}_${ploidy}x_"${selchr}"_raw
					$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  GenotypeGVCFs -R ${projdir}/refgenomes/$ref1 -L $selchr -V gendb://${pop}_${ploidy}x_"${selchr}"_raw -O ${pop}_${ploidy}x_"${selchr}"_raw.hold.vcf.gz && \
					rm -r ${pop}_${ploidy}x_"${selchr}"_raw && \
					mv ${pop}_${ploidy}x_"${selchr}"_raw.hold.vcf.gz ${pop}_${ploidy}x_"${selchr}"_raw.vcf.gz
					mv ${pop}_${ploidy}x_"${selchr}"_raw.hold.vcf.gz.tbi ${pop}_${ploidy}x_"${selchr}"_raw.vcf.gz.tbi
				done
				if LC_ALL=C gzip -l ${pop}_${ploidy}x_"${selchr}"_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
					:
				else
					rm ../cohorts*/${pop}_${ploidy}x_*_raw.vcf.gz*
					rm ../cohorts*/${pop}_${ploidy}x_raw.vcf*
					rm ../${pop}_${ploidy}x_raw_cohorts*.vcf*
					echo -e "${magenta}- \n- SNP calling failed probably due to insufficient memory ${white}\n"
					echo -e "${magenta}- \n- Exiting pipeline in 5 seconds ${white}\n"
					sleep 5 && exit 1
				fi
				)&
				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
				fi
			done
			wait
			for g in $(ls ${pop}_${ploidy}x_*_raw.vcf.gz); do
				$gunzip $g
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
			$gunzip ${pop}_${ploidy}x_raw.vcf.gz
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
wait && touch ${projdir}/GVCF_done.txt
}
cd $projdir
while [[ ! -f "$projdir/alignment_summaries/refgenome_paralogs.txt" ]]; do sleep 300; done
sleep 60
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
  	for ln in $(ls *_raw0.vcf); do
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
  	if [[ $ncontigscaffold -gt 300 ]]; then
  		echo -e "${magenta}- retrieving SNP positions based on contigs/scaffold annotation ${white}\n"
  		for nc in $(ls *_raw0.vcf); do
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
if [ "$paleopolyploid" == "true" ]; then
	mkdir 2x
	mkdir 4x
	mkdir 6x
	mkdir 8x
fi
cd ${projdir}/snpcall
if [ -z "$(ls -A *_DP_GT.txt 2>/dev/null)" ]; then
	if [ -z "$(ls -A *_x.vcf 2>/dev/null)" ]; then
		if [ -z "$(ls -A *_raw.vcf 2>/dev/null)" ]; then
			for g in *_raw.vcf.gz; do $gunzip $g;	done
		fi
	fi
fi
wait

file1xG=$( if [ "$(ls -A *_1x_DP_GT.txt 2>/dev/null)" ]; then ls *_1x_DP_GT.txt | wc -l;  else echo 0; fi )
file1xV=$( if [ "$(ls -A *_1x_raw.vcf 2>/dev/null)" ]; then ls *_1x_raw.vcf | head -n 1 | wc -l; else echo 0; fi )
if [[ "${file1xG}" -lt 1 ]]; then
	if [[ "${file1xV}" -gt 0 ]]; then
		for i in *_1x_raw.vcf; do
			$GATK LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$ref1 -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
			wait
			large_numerous_chrom &>> ${projdir}/log.out
			wait
			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf
			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%100000==2{x=file"1x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
		  wait $PID
			rm ${i%.vcf}0.vcf* ${i%.vcf}trim.vcf*
		done
		Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 2x
	fi
fi
wait
file2xG=$( if [ "$(ls -A *_2x_DP_GT.txt 2>/dev/null)" ]; then ls *_2x_DP_GT.txt | wc -l;  else echo 0; fi )
file2xV=$( if [ "$(ls -A *_2x_raw.vcf 2>/dev/null)" ]; then ls *_2x_raw.vcf | head -n 1 | wc -l; else echo 0; fi )
if [[ "${file2xG}" -lt 1 ]]; then
	if [[ "${file2xV}" -gt 0 ]]; then
		for i in *_2x_raw.vcf; do
			$GATK LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$ref1 -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
			wait
			large_numerous_chrom &>> ${projdir}/log.out
			wait
			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf
			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%100000==2{x=file"2x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
			wait $PID
			rm ${i%.vcf}0.vcf* ${i%.vcf}trim.vcf*
		done
		Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 2x
	fi
fi
wait
file4xG=$( if [ "$(ls -A *_4x_DP_GT.txt 2>/dev/null)" ]; then ls *_4x_DP_GT.txt | wc -l;  else echo 0; fi )
file4xV=$( if [ "$(ls -A *_4x_raw.vcf 2>/dev/null)" ]; then ls *_4x_raw.vcf | head -n 1 | wc -l; else echo 0; fi )
if [[ "${file4xG}" -lt 1 ]]; then
	if [[ "${file4xV}" -gt 0 ]]; then
		for i in *_4x_raw.vcf; do
			$GATK LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$ref1 -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
			wait
			large_numerous_chrom &>> ${projdir}/log.out
			wait
			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf
			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%100000==2{x=file"4x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
			wait $PID
			wait
			rm ${i%.vcf}0.vcf* ${i%.vcf}trim.vcf*
		done
		Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 4x
	fi
fi
wait
file6xG=$( if [ "$(ls -A *_6x_DP_GT.txt 2>/dev/null)" ]; then ls *_6x_DP_GT.txt | wc -l;  else echo 0; fi )
file6xV=$( if [ "$(ls -A *_6x_raw.vcf 2>/dev/null)" ]; then ls *_6x_raw.vcf | head -n 1 | wc -l; else echo 0; fi )
if [[ "${file6xG}" -lt 1 ]]; then
	if [[ "${file6xV}" -gt 0 ]]; then
		for i in *_6x_raw.vcf; do
			$GATK LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$ref1 -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
			wait
			large_numerous_chrom &>> ${projdir}/log.out
			wait
			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf
			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%100000==2{x=file"6x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
			wait $PID
			rm ${i%.vcf}0.vcf* ${i%.vcf}trim.vcf*
		done
		Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 6x
	fi
fi
wait
file8xG=$( if [ "$(ls -A *_8x_DP_GT.txt 2>/dev/null)" ]; then ls *_8x_DP_GT.txt | wc -l;  else echo 0; fi )
file8xV=$( if [ "$(ls -A *_8x_raw.vcf 2>/dev/null)" ]; then ls *_8x_raw.vcf | head -n 1 | wc -l; else echo 0; fi )
if [[ "${file8xG}" -lt 1 ]]; then
	if [[ "${file8xV}" -gt 0 ]]; then
		for i in *_8x_raw.vcf; do
			$GATK LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$ref1 -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
			wait
			large_numerous_chrom &>> ${projdir}/log.out
			wait
			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf
			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%100000==2{x=file"8x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
			wait $PID ${i%.vcf}trim.vcf*
			wait
			rm ${i%.vcf}0.vcf*
		done
		Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 8x
	fi
fi
wait
filetest=*x.vcf*
if [ -z "$(ls -A *x.vcf* 2>/dev/null)" ]; then
	for v in *_DP_GT.txt; do
		vcfdose=${v%_DP*}; vcfdose=${vcfdose#*_}; out=$(ls *${vcfdose}_raw.vcf)
		for raw in $out; do
			grep '^#' $raw  > ${raw%_raw.vcf}.vcf
			awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$1,$2]){print b}}' $raw $v >> ${raw%_raw.vcf}.vcf
		done
	done
	wait
fi

if [ "$(ls -A *.vcf 2>/dev/null)" ]; then
	for v in *.vcf; do $gzip $v; done
fi
wait


cd $projdir
cd samples
window1=$(ls -S | head -1 | xargs $zcat -fq | awk '{ print length }' | sort -n | tail -1)
window=$((window1 + 20))
cd ${projdir}

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
if [ -z $snpformats ]; then
	snpformats=false
fi

for smiss in ${sample_missingness[@]}; do
for gmiss in ${genotype_missingness[@]}; do
if [[ -z "$p1" ]]; then
	if [ -d "${projdir}/snpfilter/1x" ]; then
		cd ${projdir}/snpfilter
		cp -r 1x 1x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./1x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_1x.R $pop $gmiss $smiss $minRD_1x $exclude_samples "${GBSapp_dir}/tools/R" $maf $snpformats $haplome_number
		wait
		rm ${pop}_1x_rawRD${minRD_1x}_DP_GT.txt ${pop}_1x_DP_GT.txt ${pop}_1x_rd${minRD_1x}.txt
		awk 'NR>1{print $1,"\t",$2,"\t",$3}' *dose.txt | sort -V -k2,2 -k3,3 >  snplist_rd.txt
		chrid=$(LC_ALL=C; sort -n -k2,2 -S 50% snplist_rd.txt | awk '{print $2}' | uniq)
		for i in $chrid; do
			awk -v n="$i" '$0~n{print $0}' snplist_rd.txt | awk 'NR==0{old = $3; next} {print $1,"\t",$2,"\t",$3,"\t",$3 - old; old = $3}' |\
			awk -v f=$window1 '{print $1,"\t",$2,"\t",($4>f?$3:"")}' | awk 'NF==2{print $1,"\t",$2,"\t",p "\t" $3; next} {p=$3} 1' |\
			awk -v s=$window '{print $1"\t"$2":"$3-s"-"$3+s}' >> snplist.txt
		done
		sort -V -u -k2,2 snplist.txt > snplist_nonredun.txt; rm snplist_rd.txt
		find . -type f -empty -delete
	fi
	wait
	if [ -d "${projdir}/snpfilter/2x" ]; then
		cd ${projdir}/snpfilter
		cp -r 2x 2x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./2x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_2x.R $pop $gmiss $smiss $minRD_2x $exclude_samples "${GBSapp_dir}/tools/R" $maf $snpformats $haplome_number
		wait
		rm ${pop}_2x_rawRD${minRD_2x}_DP_GT.txt ${pop}_2x_DP_GT.txt ${pop}_2x_rd${minRD_2x}.txt
		awk 'NR>1{print $1,"\t",$2,"\t",$3}' *dose.txt | sort -V -k2,2 -k3,3 >  snplist_rd.txt
		chrid=$(LC_ALL=C; sort -n -k2,2 -S 50% snplist_rd.txt | awk '{print $2}' | uniq)
		for i in $chrid; do
			awk -v n="$i" '$0~n{print $0}' snplist_rd.txt | awk 'NR==0{old = $3; next} {print $1,"\t",$2,"\t",$3,"\t",$3 - old; old = $3}' |\
			awk -v f=$window1 '{print $1,"\t",$2,"\t",($4>f?$3:"")}' | awk 'NF==2{print $1,"\t",$2,"\t",p "\t" $3; next} {p=$3} 1' |\
			awk -v s=$window '{print $1"\t"$2":"$3-s"-"$3+s}' >> snplist.txt
		done
		sort -V -u -k2,2 snplist.txt > snplist_nonredun.txt; rm snplist_rd.txt
		find . -type f -empty -delete
	fi
	wait
	if [ -d "${projdir}/snpfilter/4x" ]; then
		cd ${projdir}/snpfilter
		cp -r 4x 4x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./4x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_4x.R $pop $gmiss $smiss $minRD_4x $exclude_samples "${GBSapp_dir}/tools/R" $maf $snpformats $haplome_number
		wait
		rm ${pop}_4x_rawRD${minRD_4x}_DP_GT.txt ${pop}_4x_DP_GT.txt ${pop}_4x_rd${minRD_4x}.txt
		awk 'NR>1{print $1,"\t",$2,"\t",$3}' *dose.txt | sort -V -k2,2 -k3,3 >  snplist_rd.txt
		chrid=$(LC_ALL=C; sort -n -k2,2 -S 50% snplist_rd.txt | awk '{print $2}' | uniq)
		for i in $chrid; do
			awk -v n="$i" '$0~n{print $0}' snplist_rd.txt | awk 'NR==0{old = $3; next} {print $1,"\t",$2,"\t",$3,"\t",$3 - old; old = $3}' |\
			awk -v f=$window1 '{print $1,"\t",$2,"\t",($4>f?$3:"")}' | awk 'NF==2{print $1,"\t",$2,"\t",p "\t" $3; next} {p=$3} 1' |\
			awk -v s=$window '{print $1"\t"$2":"$3-s"-"$3+s}' >> snplist.txt
		done
		sort -V -u -k2,2 snplist.txt > snplist_nonredun.txt; rm snplist_rd.txt
		find . -type f -empty -delete
	fi
	wait
	if [ -d "${projdir}/snpfilter/6x" ]; then
		cd ${projdir}/snpfilter
		cp -r 6x 6x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./6x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_6x.R $pop $gmiss $smiss $minRD_6x $exclude_samples "${GBSapp_dir}/tools/R" $maf $snpformats $haplome_number
		wait
		rm ${pop}_6x_rawRD${minRD_6x}_DP_GT.txt ${pop}_6x_DP_GT.txt ${pop}_6x_rd${minRD_6x}.txt
		awk 'NR>1{print $1,"\t",$2,"\t",$3}' *dose.txt | sort -V -k2,2 -k3,3 >  snplist_rd.txt
		chrid=$(LC_ALL=C; sort -n -k2,2 -S 50% snplist_rd.txt | awk '{print $2}' | uniq)
		for i in $chrid; do
			awk -v n="$i" '$0~n{print $0}' snplist_rd.txt | awk 'NR==0{old = $3; next} {print $1,"\t",$2,"\t",$3,"\t",$3 - old; old = $3}' |\
			awk -v f=$window1 '{print $1,"\t",$2,"\t",($4>f?$3:"")}' | awk 'NF==2{print $1,"\t",$2,"\t",p "\t" $3; next} {p=$3} 1' |\
			awk -v s=$window '{print $1"\t"$2":"$3-s"-"$3+s}' >> snplist.txt
		done
		sort -V -u -k2,2 snplist.txt > snplist_nonredun.txt; rm snplist_rd.txt
		find . -type f -empty -delete
	fi
	wait
	if [ -d "${projdir}/snpfilter/8x" ]; then
		cd ${projdir}/snpfilter
		cp -r 8x 8x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./8x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_8x.R $pop $gmiss $smiss $minRD_8x $exclude_samples "${GBSapp_dir}/tools/R" $maf $snpformats $haplome_number
		wait
		rm ${pop}_8x_rawRD${minRD_8x}_DP_GT.txt ${pop}_8x_DP_GT.txt ${pop}_8x_rd${minRD_8x}.txt
		awk 'NR>1{print $1,"\t",$2,"\t",$3}' *dose.txt | sort -V -k2,2 -k3,3 >  snplist_rd.txt
		chrid=$(LC_ALL=C; sort -n -k2,2 -S 50% snplist_rd.txt | awk '{print $2}' | uniq)
		for i in $chrid; do
			awk -v n="$i" '$0~n{print $0}' snplist_rd.txt | awk 'NR==0{old = $3; next} {print $1,"\t",$2,"\t",$3,"\t",$3 - old; old = $3}' |\
			awk -v f=$window1 '{print $1,"\t",$2,"\t",($4>f?$3:"")}' | awk 'NF==2{print $1,"\t",$2,"\t",p "\t" $3; next} {p=$3} 1' |\
			awk -v s=$window '{print $1"\t"$2":"$3-s"-"$3+s}' >> snplist.txt
		done
		sort -V -u -k2,2 snplist.txt > snplist_nonredun.txt; rm snplist_rd.txt
		find . -type f -empty -delete
	fi
	wait
fi

if [[ "$p1" ]]; then
	if [ -d "${projdir}/snpfilter/2x" ]; then
		cd ${projdir}/snpfilter
		cp -r 2x 2x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./2x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_2x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_2x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$snpformats" "$pseg" "$haplome_number"
		wait
		rm ${pop}_2x_rawRD${minRD_2x}_DP_GT.txt ${pop}_2x_DP_GT.txt ${pop}_2x_rd${minRD_2x}.txt
		awk 'NR>1{print $1,"\t",$2,"\t",$3}' *dose.txt | sort -V -k2,2 -k3,3 >  snplist_rd.txt
		chrid=$(LC_ALL=C; sort -n -k2,2 -S 50% snplist_rd.txt | awk '{print $2}' | uniq)
		for i in $chrid; do
			awk -v n="$i" '$0~n{print $0}' snplist_rd.txt | awk 'NR==0{old = $3; next} {print $1,"\t",$2,"\t",$3,"\t",$3 - old; old = $3}' |\
			awk -v f=$window1 '{print $1,"\t",$2,"\t",($4>f?$3:"")}' | awk 'NF==2{print $1,"\t",$2,"\t",p "\t" $3; next} {p=$3} 1' |\
			awk -v s=$window '{print $1"\t"$2":"$3-s"-"$3+s}' >> snplist.txt
		done
		sort -V -u -k2,2 snplist.txt > snplist_nonredun.txt; rm snplist_rd.txt
		find . -type f -empty -delete
	fi
	wait
	if [ -d "${projdir}/snpfilter/4x" ]; then
		cd ${projdir}/snpfilter
		cp -r 4x 4x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./4x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_4x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_4x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$snpformats" "$pseg" "$haplome_number"
		wait
		rm ${pop}_4x_rawRD${minRD_4x}_DP_GT.txt ${pop}_4x_DP_GT.txt ${pop}_4x_rd${minRD_4x}.txt
		awk 'NR>1{print $1,"\t",$2,"\t",$3}' *dose.txt | sort -V -k2,2 -k3,3 >  snplist_rd.txt
		chrid=$(LC_ALL=C; sort -n -k2,2 -S 50% snplist_rd.txt | awk '{print $2}' | uniq)
		for i in $chrid; do
			awk -v n="$i" '$0~n{print $0}' snplist_rd.txt | awk 'NR==0{old = $3; next} {print $1,"\t",$2,"\t",$3,"\t",$3 - old; old = $3}' |\
			awk -v f=$window1 '{print $1,"\t",$2,"\t",($4>f?$3:"")}' | awk 'NF==2{print $1,"\t",$2,"\t",p "\t" $3; next} {p=$3} 1' |\
			awk -v s=$window '{print $1"\t"$2":"$3-s"-"$3+s}' >> snplist.txt
		done
		sort -V -u -k2,2 snplist.txt > snplist_nonredun.txt; rm snplist_rd.txt
		find . -type f -empty -delete
	fi
	wait
	if [ -d "${projdir}/snpfilter/6x" ]; then
		cd ${projdir}/snpfilter
		cp -r 6x 6x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./6x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_6x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_6x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$snpformats" "$pseg" "$haplome_number"
		wait
		rm ${pop}_6x_rawRD${minRD_6x}_DP_GT.txt ${pop}_6x_DP_GT.txt ${pop}_6x_rd${minRD_6x}.txt
		awk 'NR>1{print $1,"\t",$2,"\t",$3}' *dose.txt | sort -V -k2,2 -k3,3 >  snplist_rd.txt
		chrid=$(LC_ALL=C; sort -n -k2,2 -S 50% snplist_rd.txt | awk '{print $2}' | uniq)
		for i in $chrid; do
			awk -v n="$i" '$0~n{print $0}' snplist_rd.txt | awk 'NR==0{old = $3; next} {print $1,"\t",$2,"\t",$3,"\t",$3 - old; old = $3}' |\
			awk -v f=$window1 '{print $1,"\t",$2,"\t",($4>f?$3:"")}' | awk 'NF==2{print $1,"\t",$2,"\t",p "\t" $3; next} {p=$3} 1' |\
			awk -v s=$window '{print $1"\t"$2":"$3-s"-"$3+s}' >> snplist.txt
		done
		sort -V -u -k2,2 snplist.txt > snplist_nonredun.txt; rm snplist_rd.txt
		find . -type f -empty -delete
	fi
	wait
	if [ -d "${projdir}/snpfilter/8x" ]; then
		cd ${projdir}/snpfilter
		cp -r 8x 8x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./8x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_8x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_8x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$snpformats" "$pseg" "$haplome_number"
		wait
		rm ${pop}_8x_rawRD${minRD_8x}_DP_GT.txt ${pop}_8x_DP_GT.txt ${pop}_8x_rd${minRD_8x}.txt
		awk 'NR>1{print $1,"\t",$2,"\t",$3}' *dose.txt | sort -V -k2,2 -k3,3 >  snplist_rd.txt
		chrid=$(LC_ALL=C; sort -n -k2,2 -S 50% snplist_rd.txt | awk '{print $2}' | uniq)
		for i in $chrid; do
			awk -v n="$i" '$0~n{print $0}' snplist_rd.txt | awk 'NR==0{old = $3; next} {print $1,"\t",$2,"\t",$3,"\t",$3 - old; old = $3}' |\
			awk -v f=$window1 '{print $1,"\t",$2,"\t",($4>f?$3:"")}' | awk 'NF==2{print $1,"\t",$2,"\t",p "\t" $3; next} {p=$3} 1' |\
			awk -v s=$window '{print $1"\t"$2":"$3-s"-"$3+s}' >> snplist.txt
		done
		sort -V -u -k2,2 snplist.txt > snplist_nonredun.txt; rm snplist_rd.txt
		find . -type f -empty -delete
	fi
fi
done
wait
done


cd ${projdir}/snpfilter/
find . -type f -empty -delete
find . -type d -empty -delete
wc -l *gmiss*/*dose.txt | awk '{print $2"\t"$1-1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > gmiss_smiss_titration.txt
wc -l *gmiss*/eliminated* | awk '{print $2"\t"$1-1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > eliminated_samples.txt
echo -e "gmiss_smiss_thresholds\t#_of_SNPs\t#_of_eliminated_samples\n----------------------\t---------\t-----------------------" > summary_precall.txt
awk 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,"\t",a[$1],"\t",$2}' gmiss_smiss_titration.txt eliminated_samples.txt | cat summary_precall.txt - > gmiss_smiss_titration_eliminated_samples.txt
rm gmiss_smiss_titration.txt eliminated_samples.txt summary_precall.txt


cd "$projdir"/snpfilter
n="${ref1%.f*}_"
for snpfilter_dir in $(ls -d */); do
	if [ -d "$snpfilter_dir" ]; then
		cd $snpfilter_dir
		for v in *dose.txt; do
			vcfdose=${v%_rd*}; vcfdose=${vcfdose#*_}
			$zcat ../../snpcall/*${vcfdose}.vcf.gz | grep '^#' > ${v%.txt}.vcf
			awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$2,$3]){print b}}' <($gzip -dc ../../snpcall/*${vcfdose}.vcf.gz) $v >> ${v%.txt}.vcf
			$gzip ${v%.txt}.vcf
		done
		for i in *dose.txt *binary.txt; do
			awk -v n="$n" '{gsub(n,""); print $0}' $i > ${i%.txt}_hold.txt
			mv ${i%.txt}_hold.txt $i
		done
		if [[ "${snpformats}" == "true" ]]; then
		for i in *nucleotide.txt *nucleotidedeg.txt; do
			awk -v n="$n" '{gsub(n,""); print $0}' $i > ${i%.txt}_hold.txt
			mv ${i%.txt}_hold.txt $i
		done
		fi
		wait
		cd ../
	fi
done
}
cd $projdir
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
		awk -F '\t' 'BEGIN{OFS="\t"}{print $1}' $i | cat - ./paralog_haplo_filter/SNP_files.txt | sort | uniq -u | awk '{print $0"\tNA\tNA"}' | cat $i - | sort -u | sort -k1,1 >> ${i%_haplotypes.txt}_hapmissingSNP.txt
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


######################################################################################################################################################
cd ${projdir}
rm compress_done.txt hapfilter_done.txt alignment_done.txt precall_done.txt GVCF_done.txt
touch Analysis_Complete
wait
echo -e "${magenta}- Run Complete. ${white}\n"
