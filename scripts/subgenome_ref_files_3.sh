
if [ -z "$threads" ]; then
	threads=$(nproc --all)
	if [[ "$threads" -ge 4 ]]; then
		threads=$((threads-2))
	fi
fi
if [ -z "$lib_type" ]; then
 lib_type=RRS
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
if [ -z "$multilocus" ]; then
	multilocus=true
fi
if [ "$multilocus" == "true" ]; then
	multilocus=0
fi
if [ "$multilocus" == "false" ]; then
	multilocus=null
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
	mhap_freq=0
fi
if [ -z "$softclip" ]; then
	softclip=false
fi
if [ -z "$downsample" ]; then
	downsample=0
fi
if [ -z "$joint_calling" ]; then
	joint_calling=false
fi
if [ -z "$keep_gVCF" ]; then
	keep_gVCF=false
fi



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
	checknfastafiles=$(ls *.f* | grep -v .fai | grep -v .ngm | grep -v panref | grep -v _original.fasta | wc -l)
	if [[ $checknfastafiles -gt 3 ]]; then
		echo -e "${magenta}- expecting only 2 fasta files for reference genome ${white}\n"
		echo -e "${magenta}- GBSapp will quit in 5 seconds ${white}\n"
		sleep 5 && exit 1
	fi

	reffilenames=($ref1 $ref2)
	for refg in "${reffilenames[@]}"; do
		export ncontigscaffold=$(grep '>' $refg | wc -l)
		if [[ $ncontigscaffold -gt 300 ]]; then
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
				Nstitch=$(printf "A%.0s" $(seq 100))
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
	if [[ $ncontigscaffold -gt 300 ]]; then
		mkdir split
		for reffile in $(ls *.f*); do
			awk -v reffile=${reffile%.f*} '/^>/{close("file.fasta"f);f++}{print $0 > "./split/file_"reffile"_"f}' "$reffile" & PID=$!
		  wait $PID
		done
		cd split
		export checksplit=$( wc -c file_* | head -n -1 | awk '($1 > 500000000 )' | wc -l )
		if [[ "$checksplit" -gt 0 ]]; then
			for i in $(ls file*); do
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
file=panref.fasta
if test -f "$file"; then
	echo -e "${magenta}- combined reference assemblies of subgenomes already exists ${white}\n"
else
	echo -e "${magenta}- combining reference assembly of subgenomes ${white}\n"
	touch panref.txt
	for i in $(ls *.f*); do
		n=">${i%.f*}_"
		awk '{ sub("\r$",""); print}' $i | awk -v n="$n" '{gsub(n,""); print}' | awk -v n="$n" '{gsub(/>/,n); print}' >> panref.txt
	done
	awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' panref.txt > panref.fasta
	rm panref.txt
fi

if ls ./*.ngm 1> /dev/null 2>&1; then
	echo -e "${magenta}- indexed genome already exist ${white}\n"
else
	echo -e "${magenta}- indexing reference genome ${white}\n"
	$samtools faidx panref.fasta
	$java -jar $picard CreateSequenceDictionary REFERENCE= panref.fasta OUTPUT=panref.dict
	$ngm -r panref.fasta
fi

declare -a arr=("panref.dict" "panref.fasta" "panref.fasta.fai")
for file in "${arr[@]}"; do
if test -f $file; then
	:
else
	echo -e "${cyan}- reference genome not indexed: missing $file ${white}\n"
	echo -e "${cyan}- make sure reference genome is indexed ${white}\n"
	echo -e "${cyan}- GBSapp will quit in 5 seconds ${white}\n"
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
		if [[ "$fqpass" -gt 0 ]]; then
			if [[ "$samples_list" == "samples_list_node_1.txt" ]] && test ! -f ${projdir}/alignment_summaries/total_read_count.txt; then
				time main &>> ${projdir}/log.out
			fi
		fi
	else
		echo -e "${magenta}- samples' PE fastq filenames requires formatting (i.e. needs to end in "_R1.fastq" and "_R2.fastq") ${white}\n"
		sleep 5 && exit 1
	fi
else
	if [[ "$samples_list" == "samples_list_node_1.txt" ]] && test ! -f ${projdir}/alignment_summaries/total_read_count.txt; then
		if [[ ! -f "${projdir}/organize_files_done.txt" ]]; then time main &>> ${projdir}/log.out; fi
	fi
fi


main () {
	cd ${projdir}
	mkdir -p preprocess
	mkdir -p snpcall
	mkdir -p alignment_summaries
	mkdir -p ./alignment_summaries/copy_number
	if [[ ! -f ${projdir}/alignment_summaries/total_read_count.txt ]]; then
		touch ${projdir}/alignment_summaries/total_read_count.txt
	fi
	printf 'sample\tnumber_of_reads\n' > ${projdir}/alignment_summaries/total_read_count.txt
	if [[ "$mhap_freq" -gt 0 ]]; then
		mkdir -p ${projdir}/alignment_summaries/background_mutation_test
		if [[ ! -f ${projdir}/alignment_summaries/pop_mutation_load.txt ]]; then
			printf 'sample\trare_hap\tpop_hap\ttotal_hap\tmutation_load\n' > ${projdir}/alignment_summaries/pop_mutation_load.txt
		fi
	fi
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
	ram2=$(echo "$totalk*0.0000009" | bc)
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
		gthread=threads
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
echo -e "${blue}\n############################################################################## ${yellow}\n- GBSapp is Performing Read Alignments & Alignment Post-Processing\n${blue}##############################################################################${white}\n"

main () {
	cd $projdir
	cd samples

	for i in $( cat ${projdir}/${samples_list} ); do
		if [[ "$lib_type" == "RRS" ]] && test ! -f ${projdir}/compress_done.txt && test ! -f ${projdir}/organize_files_done.txt && test ! -f ${projdir}/preprocess/${i%.f*}_redun.sam && test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai; then
			if test ! -f ${i%.f*}_uniq_R1.fasta.gz; then
				if [[ $(file $i | awk -F' ' '{print $2}') == gzip ]]; then
					zcat $i | awk 'NR%2==0' | awk 'NR%2' | $gzip > ${i%.f*}_uniq.txt.gz 2> /dev/null &&
					wait
					if test -f ${i%.f*}_R2*; then
						zcat ${i%.f*}_R2* | awk 'NR%2==0' | awk 'NR%2' | $gzip > ${i%.f*}_R2_uniq.txt.gz 2> /dev/null &&
						wait
					fi
					wait
					if test -f ${i%.f*}.R2*; then
						zcat ${i%.f*}.R2* | awk 'NR%2==0' | awk 'NR%2' | $gzip > ${i%.f*}_R2_uniq.txt.gz 2> /dev/null &&
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
					zcat ${i%.f*}_uniq.txt.gz | printf "${i%.f*}""\t""$(wc -l)""\n" > ${projdir}/alignment_summaries/${i%.f*}_total_read_count.txt  2> /dev/null &&
					wait
					export LC_ALL=C; paste -d ~ <(zcat ${i%.f*}_uniq.txt.gz) <(zcat ${i%.f*}_R2_uniq.txt.gz 2> /dev/null) | expand -t $(( $(wc -L < $i ) + 2 )) | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{$1=$1};1' | \
					awk '{gsub(" /"," "); print}' | awk '{gsub("/\n","\n"); print}' | awk '{gsub("/"," "); print}' | awk '{gsub(" ","\t"); print}' | $gzip > ${i%.f*}_rdrefseq.txt.gz 2> /dev/null &&
					wait
					awk 'NF==2 {print ">seq"NR"_se-"$1"\t"$2}' <(zcat ${i%.f*}_rdrefseq.txt.gz 2> /dev/null) | $gzip > ${i%.f*}_rdrefseq_se.txt.gz 2> /dev/null &&
					wait
					awk 'NF==3 {print ">seq"NR"_pe-"$0}' <(zcat ${i%.f*}_rdrefseq.txt.gz 2> /dev/null) | awk '{print $1"\t"$3}' | $gzip > ${i%.f*}_uniq_R2.fasta.gz 2> /dev/null &&
					wait
					awk 'NF==3 {print ">seq"NR"_pe-"$0}' <(zcat ${i%.f*}_rdrefseq.txt.gz 2> /dev/null) | awk '{print $1"\t"$2}' | $gzip | cat - <(zcat ${i%.f*}_rdrefseq_se.txt.gz) > ${i%.f*}_uniq_R1.hold.fasta.gz 2> /dev/null &&
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
					zcat ${i%.f*}_uniq.txt.gz | printf "${i%.f*}""\t""$(wc -l)""\n" > ${projdir}/alignment_summaries/${i%.f*}_total_read_count.txt  2> /dev/null &&
					wait
					export LC_ALL=C; paste -d ~ <(zcat ${i%.f*}_uniq.txt.gz) ${i%.f*}_R2_uniq.txt | expand -t $(( $(wc -L < $i ) + 2 )) | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{$1=$1};1' | \
					awk '{gsub(" /"," "); print}' | awk '{gsub("/\n","\n"); print}' | awk '{gsub("/"," "); print}' | awk '{gsub(" ","\t"); print}' | $gzip > ${i%.f*}_rdrefseq.txt.gz 2> /dev/null &&
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


	cd $projdir/samples
	find . -size 0 -delete  2> /dev/null &&
	touch ../report_fq_compress_index.txt
	for i in $( cat ${projdir}/${samples_list} ); do
		if [[ "$(zcat ${i%.f*}_uniq_R1.fasta.gz 2> /dev/null | head -n1 | wc -l)" -eq 0 ]] || [[ -z "${i%.f*}_uniq_R1.fasta.gz" ]]; then
			echo ${i%.f*}_uniq_R1.fasta.gz >> ../report_fq_compress_index.txt
		fi
	done
	if [[ $(wc -l ../report_fq_compress_index.txt | awk '{print $1}' 2> /dev/null) -eq 0 ]]; then
		echo -e "${magenta}- $(wc -l ../report_fq_compress_index.txt | awk '{print $1}') fastq files not properly processed ${white}\n"
		END=10
		while [[ $END -gt 0 ]]; do
			for i in $( cat ${projdir}/${samples_list} ); do
				if [[ "$lib_type" == "RRS" ]] && test ! -f ${projdir}/compress_done.txt && test ! -f ${projdir}/organize_files_done.txt && test ! -f ${projdir}/preprocess/${i%.f*}_redun.sam && test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai; then
					if test ! -f ${i%.f*}_uniq_R1.fasta.gz; then
						if [[ $(file $i | awk -F' ' '{print $2}') == gzip ]]; then
							zcat $i | awk 'NR%2==0' | awk 'NR%2' | $gzip > ${i%.f*}_uniq.txt.gz 2> /dev/null &&
							wait
							if test -f ${i%.f*}_R2*; then
								zcat ${i%.f*}_R2* | awk 'NR%2==0' | awk 'NR%2' | $gzip > ${i%.f*}_R2_uniq.txt.gz 2> /dev/null &&
								wait
							fi
							wait
							if test -f ${i%.f*}.R2*; then
								zcat ${i%.f*}.R2* | awk 'NR%2==0' | awk 'NR%2' | $gzip > ${i%.f*}_R2_uniq.txt.gz 2> /dev/null &&
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
							zcat ${i%.f*}_uniq.txt.gz | printf "${i%.f*}""\t""$(wc -l)""\n" > ${projdir}/alignment_summaries/${i%.f*}_total_read_count.txt  2> /dev/null &&
							wait
							export LC_ALL=C; paste -d ~ <(zcat ${i%.f*}_uniq.txt.gz) <(zcat ${i%.f*}_R2_uniq.txt.gz 2> /dev/null) | expand -t $(( $(wc -L < $i ) + 2 )) | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{$1=$1};1' | \
							awk '{gsub(" /"," "); print}' | awk '{gsub("/\n","\n"); print}' | awk '{gsub("/"," "); print}' | awk '{gsub(" ","\t"); print}' | $gzip > ${i%.f*}_rdrefseq.txt.gz 2> /dev/null &&
							wait
							awk 'NF==2 {print ">seq"NR"_se-"$1"\t"$2}' <(zcat ${i%.f*}_rdrefseq.txt.gz 2> /dev/null) | $gzip > ${i%.f*}_rdrefseq_se.txt.gz 2> /dev/null &&
							wait
							awk 'NF==3 {print ">seq"NR"_pe-"$0}' <(zcat ${i%.f*}_rdrefseq.txt.gz 2> /dev/null) | awk '{print $1"\t"$3}' | $gzip > ${i%.f*}_uniq_R2.fasta.gz 2> /dev/null &&
							wait
							awk 'NF==3 {print ">seq"NR"_pe-"$0}' <(zcat ${i%.f*}_rdrefseq.txt.gz 2> /dev/null) | awk '{print $1"\t"$2}' | $gzip | cat - <(zcat ${i%.f*}_rdrefseq_se.txt.gz) > ${i%.f*}_uniq_R1.hold.fasta.gz 2> /dev/null &&
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
							zcat ${i%.f*}_uniq.txt.gz | printf "${i%.f*}""\t""$(wc -l)""\n" > ${projdir}/alignment_summaries/${i%.f*}_total_read_count.txt  2> /dev/null &&
							wait
							export LC_ALL=C; paste -d ~ <(zcat ${i%.f*}_uniq.txt.gz) ${i%.f*}_R2_uniq.txt | expand -t $(( $(wc -L < $i ) + 2 )) | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{$1=$1};1' | \
							awk '{gsub(" /"," "); print}' | awk '{gsub("/\n","\n"); print}' | awk '{gsub("/"," "); print}' | awk '{gsub(" ","\t"); print}' | $gzip > ${i%.f*}_rdrefseq.txt.gz 2> /dev/null &&
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
	touch ${projdir}/organize_files_done.txt



	if [[ "$lib_type" == "RRS" ]] && [[ "$(wc -l ${projdir}/alignment_summaries/total_read_count.txt | awk '{print $1}')" -le 1 ]]; then
		cd ${projdir}/alignment_summaries/
		find -type f -name "*_total_read_count.txt" | xargs cat > total_read_count.hold.txt &&
		cat total_read_count.txt total_read_count.hold.txt > total_read_count.hold2.txt &&
		mv total_read_count.hold2.txt total_read_count.txt &&
		rm ${projdir}/alignment_summaries/*_total_read_count.txt ${projdir}/alignment_summaries/total_read_count.hold.txt &&
		wait
	fi
	wait


	if [[ "$lib_type" == "RRS" ]] && [[ ! -f "${projdir}/align1_${samples_list}" ]] && test ! -f ${projdir}/compress_done.txt; then touch "${projdir}/align1_${samples_list}"; fi
	if [[ "$samples_list" == "samples_list_node_1.txt" ]] && [[ "$mhap_freq" -gt 0 ]] && [[ ! -f "${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqFail.txt" ]]; then
		align=$(ls ${projdir}/align1_samples_list_node_* | wc -l)
		while [[ "$align" -lt $nodes ]]; do
			sleep 30; align=$(ls ${projdir}/align1_samples_list_node_* | wc -l)
		done
		for i in $( cat ${projdir}/${samples_list} ); do
			zcat ${projdir}/samples/${i%.f*}_uniq_R1.fasta.gz | awk -F "\t" 'BEGIN { OFS=FS }; { print $1, substr($2, 1, 64); }' | awk '{a[$2]++} END{for(s in a){print a[s]" "s}}' | \
			awk -F '\t' '{gsub(/ /,"\t"); print}' | awk -F '\t' 'BEGIN{OFS="\t"} NR==FNR{a[$2]=$0;next} ($2) in a{print $0, a[$2]}' - <(zcat ${projdir}/samples/${i%.f*}_uniq_R1.fasta.gz) | \
			awk '$3==1{print $1"\t"$2}' | $gzip > ${projdir}/alignment_summaries/background_mutation_test/${i%.f*}_pop_haps.fasta.gz 2> /dev/null &&
			wait
		done

		if [[ $align == $nodes ]] && test ! -f ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqFail.txt; then
			cd ${projdir}/alignment_summaries/background_mutation_test/
			touch ${projdir}/alignment_summaries/background_mutation_test/pop_haps.fasta
			find -type f -name "*_pop_haps.fasta.gz" | xargs cat > pop_haps.txt.gz &&
			rm *_pop_haps.fasta.gz ${projdir}/align1_${samples_list} &&
			zcat pop_haps.txt.gz | awk '{a[$2]++} END{for(s in a){print a[s]" "s}}' | awk -v phap=$mhap_freq '($1 <= phap)' | awk -F '\t' '{gsub(/ /,"\t"); print}' | \
			awk -F '\t' 'BEGIN{OFS="\t"} NR==FNR{a[$2]=$0;next} ($2) in a{print $0, a[$2]}' - <(zcat pop_haps.txt.gz) | \
			awk '{print $1"\t"$2"\t"$3}'  > pop_haps_freqFail.txt &&
			cd ../
		fi
	fi
	wait

	cd ${projdir}/samples

	for i in $(cat ${projdir}/${samples_list} ); do
		if [[ "$lib_type" == "RRS" ]] && test ! -f ${projdir}/hapfilter_done.txt && test ! -f ${projdir}/compress_done.txt && test ! -f "${projdir}/preprocess/${i%.f*}_redun.sam.gz" && test ! -f "${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai"; then
				sleep $[ ( $RANDOM % 30 )  + 10 ]s
				export nempty=$( ls ${projdir}/samples/${i%.f*}_uniq_R2.fasta.gz 2> /dev/null | wc -l | awk '{print $1}' )
				 if [[ "$mhap_freq" -gt 0 ]]; then
					 if [[ "$nempty" -gt 0 ]]; then
					  awk '{print $2}' ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqFail.txt | \
					  grep -v -x -f - <(zcat ${projdir}/samples/${i%.f*}_uniq_R1.fasta.gz) | $gzip > ${projdir}/samples/${i%.f*}_uniq_R1_singleton.txt.gz &&
						grep '_se-' <(zcat ${projdir}/samples/${i%.f*}_uniq_R1_singleton.txt.gz) | awk '{gsub(/>/,"@"); print}' | awk '{print $1"\t"$2"\t"$2}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | awk '{print $1"\n"$2"\n+\n"$3}' | $gzip > ${projdir}/samples/${i%.f*}_uniq_singleton.fq.gz &&
						grep '_pe-' <(zcat ${projdir}/samples/${i%.f*}_uniq_R1_singleton.txt.gz) | awk '{gsub(/>/,"@"); print}' | awk '{print $1"\t"$2"\t"$2}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | awk '{print $1"/1\n"$2"\n+\n"$3}' | $gzip > ${projdir}/samples/${i%.f*}_uniq_R1.fq.gz &&
						rm ${projdir}/samples/${i%.f*}_uniq_R1_singleton.txt.gz &&
						awk '{print $2}' ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqFail.txt | \
						grep -v -x -f - <(zcat ${projdir}/samples/${i%.f*}_uniq_R2.fasta.gz) | awk '{gsub(/>/,"@"); print}' | \
						awk '{print $1"\t"$2"\t"$2}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | awk '{print $1"/2\n"$2"\n+\n"$3}' | $gzip > ${projdir}/samples/${i%.f*}_uniq_R2.fq.gz &&
						rm ${projdir}/samples/${i%.f*}_uniq_R1.fasta.gz ${projdir}/samples/${i%.f*}_uniq_R2.fasta.gz
						wait
					else
						awk '{print $2}' ${projdir}/alignment_summaries/background_mutation_test/pop_haps_freqFail.txt | \
						grep -v -x -f - <(zcat ${projdir}/samples/${i%.f*}_uniq_R1.fasta.gz) | awk '{gsub(/>/,"@"); print}' | \
						awk '{print $1"\t"$2"\t"$2}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | awk '{print $1"\n"$2"\n+\n"$3}' | $gzip > ${projdir}/samples/${i%.f*}_uniq_R1.fq.gz &&
						rm ${projdir}/samples/${i%.f*}_uniq_R1.fasta.gz
						wait
					fi
				 else
					 if [[ "$nempty" -gt 0 ]]; then
						grep '_se-' <(zcat ${projdir}/samples/${i%.f*}_uniq_R1.fasta.gz) | awk '{gsub(/>/,"@"); print}' | awk '{print $1"\t"$2"\t"$2}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | awk '{print $1"\n"$2"\n+\n"$3}' | $gzip > ${projdir}/samples/${i%.f*}_uniq_singleton.fq.gz &&
						grep '_pe-' <(zcat ${projdir}/samples/${i%.f*}_uniq_R1.fasta.gz) | awk '{gsub(/>/,"@"); print}' | awk '{print $1"\t"$2"\t"$2}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | awk '{print $1"/1\n"$2"\n+\n"$3}' | $gzip > ${projdir}/samples/${i%.f*}_uniq_R1.fq.gz &&
						awk '{gsub(/>/,"@"); print}' <(zcat ${projdir}/samples/${i%.f*}_uniq_R2.fasta.gz) | awk '{print $1"\t"$2"\t"$2}' | \
						awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | awk '{print $1"/2\n"$2"\n+\n"$3}' | $gzip > ${projdir}/samples/${i%.f*}_uniq_R2.fq.gz &&
						rm ${projdir}/samples/${i%.f*}_uniq_R1.fasta.gz ${projdir}/samples/${i%.f*}_uniq_R2.fasta.gz
						wait
					else
						awk '{gsub(/>/,"@"); print}' <(zcat ${projdir}/samples/${i%.f*}_uniq_R1.fasta.gz) | awk '{print $1"\t"$2"\t"$2}' | \
						awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$3); print}' | awk '{print $1"\n"$2"\n+\n"$3}' | gzip  > ${projdir}/samples/${i%.f*}_uniq_R1.fq.gz &&
						rm ${projdir}/samples/${i%.f*}_uniq_R1.fasta.gz
						wait
					fi
				 fi
		fi
	done
	wait && touch ${projdir}/compress_done.txt


	for i in $(cat ${projdir}/${samples_list} ); do
		cd ${projdir}/refgenomes
		if [[ "$lib_type" == "RRS" ]] && test ! -f ${projdir}/precall_done.txt && test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai; then
			export nempty=$( ls ${projdir}/samples/${i%.f*}_uniq_R2.fq.gz 2> /dev/null | wc -l | awk '{print $1}' )
			if test ! -f ${projdir}/preprocess/${i%.f*}_redun.sam.gz; then
				if [[ "$nempty" -gt 0 ]]; then
					$ngm -r panref.fasta --qry ${projdir}/samples/${i%.f*}_uniq_R1.fq.gz -o ${projdir}/preprocess/${i%.f*}_redun_R1.sam -t $threads --min-identity 0 --topn 12 -strata 12 &&
					$ngm -r panref.fasta --qry ${projdir}/samples/${i%.f*}_uniq_R2.fq.gz -o ${projdir}/preprocess/${i%.f*}_redun_R2.sam -t $threads --min-identity 0 --topn 12 --strata 12 &&
					if [[ -s ${projdir}/preprocess/${i%.f*}_redun_singleton.sam ]]; then $ngm -r $ref1 --qry ${projdir}/samples/${i%.f*}_uniq_singleton.fq.gz -o ${projdir}/preprocess/${i%.f*}_redun_singleton.sam -t $threads --min-identity 0 --topn 12 --strata 12; fi &&
					mv ${projdir}/preprocess/${i%.f*}_redun_R1.sam ${projdir}/preprocess/${i%.f*}_redun.hold.sam &&
					grep -v '^@' ${projdir}/preprocess/${i%.f*}_redun_R2.sam >> ${projdir}/preprocess/${i%.f*}_redun.hold.sam &&
					if [[ -s ${projdir}/preprocess/${i%.f*}_redun_singleton.sam ]]; then grep -v '^@' ${projdir}/preprocess/${i%.f*}_redun_singleton.sam >> ${projdir}/preprocess/${i%.f*}_redun.hold.sam; fi &&
					$gzip ${projdir}/preprocess/${i%.f*}_redun.hold.sam &&
					mv ${projdir}/preprocess/${i%.f*}_redun.hold.sam.gz ${projdir}/preprocess/${i%.f*}_redun.sam.gz 2> /dev/null &&
					wait
				else
					$ngm -r panref.fasta -q ${projdir}/samples/${i%.f*}_uniq_R1.fq.gz -o ${projdir}/preprocess/${i%.f*}_redun.hold.sam -t $threads --min-identity 0 --topn 12 --strata 12 &&
					$gzip ${projdir}/preprocess/${i%.f*}_redun.hold.sam &&
					mv ${projdir}/preprocess/${i%.f*}_redun.hold.sam.gz ${projdir}/preprocess/${i%.f*}_redun.sam.gz 2> /dev/null &&
					rm ${projdir}/samples/${i%.f*}_uniq_*.fq.gz &&
					wait
				fi
				rm ${projdir}/preprocess/${i%.f*}_redun_singleton.sam ${projdir}/preprocess/${i%.f*}_redun_R*.sam 2> /dev/null &&
				wait
			fi
			wait
		fi
		cd ${projdir}/samples
	done

	for i in $(cat ${projdir}/${samples_list} ); do
		cd ${projdir}/refgenomes
		if [[ "$lib_type" == "WGS" ]] && test ! -f ${projdir}/precall_done.txt && test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai; then
			export nempty=$( ls ${projdir}/samples/${i%.f*}_R2.f*.gz 2> /dev/null | wc -l | awk '{print $1}' )
			if test ! -f ${projdir}/preprocess/${i%.f*}_redun.sam.gz; then
				if [[ "$nempty" -gt 0 ]]; then
					$ngm -r panref.fasta --qry ${projdir}/samples/$i -o ${projdir}/preprocess/${i%.f*}_R1.sam -t $threads --min-identity 0 --topn 12 --strata 12 &&
					$ngm -r panref.fasta --qry ${projdir}/samples/${i%.f*}_R2.fastq.gz -o ${projdir}/preprocess/${i%.f*}_R2.sam -t $threads --min-identity 0 --topn 12 --strata 12 &&
					mv ${projdir}/preprocess/${i%.f*}_redun_R1.sam ${projdir}/preprocess/${i%.f*}_redun.hold.sam &&
					grep -v '^@' ${projdir}/preprocess/${i%.f*}_redun_R2.sam >> ${projdir}/preprocess/${i%.f*}_redun.hold.sam &&
					$gzip ${projdir}/preprocess/${i%.f*}_redun.hold.sam &&
					rm ${projdir}/preprocess/${i%.f*}_redun_R*.sam &&
					wait
				else
					$ngm -r panref.fasta --qry ${projdir}/samples/$i -o ${projdir}/preprocess/${i%.f*}_redun.hold.sam -t $threads --min-identity 0 --topn 12 --strata 12 &&
					$gzip ${projdir}/preprocess/${i%.f*}_redun.hold.sam &&
					wait
				fi
				mv ${projdir}/preprocess/${i%.f*}_redun.hold.sam.gz ${projdir}/preprocess/${i%.f*}_redun.sam.gz 2> /dev/null &&
				wait
			fi
			wait
		fi
		cd ${projdir}/samples
	done
	wait && touch ${projdir}/hapfilter_done.txt


	for i in $( cat ${projdir}/${samples_list} ); do
		if test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai && test ! -f ${projdir}/preprocess/${i%.f*}_${ref2%.f*}_precall.bam.bai && test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_${ref2%.f*}_precall.bam.bai; then
			if test ! -f ${projdir}/alignment_done.txt; then
				printf '\n###---'${i%.f*}'---###\n' > ${projdir}/alignment_summaries/${i%.f*}_summ.txt && \
				$samtools flagstat <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz) >> ${projdir}/alignment_summaries/${i%.f*}_summ.txt && \
				printf '########################################################################################################\n\n' >> ${projdir}/alignment_summaries/${i%.f*}_summ.txt && \
				printf 'copy\tFrequency\tPercentage\n' > ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_Unique_Read_histogram.txt && \
				grep -v '^@' <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz) | awk -F' ' '{print $1}' | awk '{gsub(/_/,"\t");gsub(/\//,"\t");gsub(/pe-/,"");gsub(/se-/,""); print $2}' | \
				sort | uniq -c | awk '{print $2"\t"$1}' > ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number.txt  && \
				awk 'NR==FNR{sum+= $2; next;} {printf("%s\t%s\t%3.3f%%\t%3.0f\n",$1,$2,100*$2/sum,100*$2/sum)}' ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number.txt ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number.txt > ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot.txt && \
				unset IFS; printf "%s\t%s\t%s\t%*s\n" $(sed 's/$/ |/' ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot.txt) | tr ' ' '|' | sort -k2,2 -nr | awk '{gsub(/se-/,""); gsub(/pe-/,""); print}' >> ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number_Unique_Read_histogram.txt && rm ${projdir}/alignment_summaries/copy_number/${i%.f*}_copy_number.txt ${projdir}/alignment_summaries/copy_number/${i%.f*}_plot.txt

				grep -v '^@' <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz) | awk '($3 != "\*")' | awk '($6 != "\*")' | awk '{print $3"\t"$3"\t"$0}' | \
				awk '{gsub(/_.*$/,"",$1); gsub(/_.*$/,"",$2)}1' > ${projdir}/preprocess/${i%.f*}_Index0_subgenome.txt
				awk '{$2=$2"_"$3}1' ${projdir}/preprocess/${i%.f*}_Index0_subgenome.txt | awk '!h[$2] { g[$2]=$0 } { h[$2]++ } END { for(k in g) print h[k], g[k] }' | \
				awk '!h[$4] { g[$4]=$0 } { h[$4]++ } END { for(k in g) print h[k], g[k] }' | awk '{print $1"\t"$5"\t"$3}' > ${projdir}/preprocess/${i%.f*}_Index_subgenome.txt

        awk '{if($1==3) print $2}' ${projdir}/preprocess/${i%.f*}_Index_subgenome.txt | grep -Fw -f - <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz) | cat <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz | grep '^@') - > ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam

        awk '{if($1==2) print $0}' ${projdir}/preprocess/${i%.f*}_Index_subgenome.txt | awk -v ref1=${ref1%.f*} -v ref2=${ref2%.f*} '{if($3 == ref1 || $3 == ref2) print $2}' | grep -Fw -f - <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz) > ${projdir}/preprocess/${i%.f*}_del0_${ref1%.f*}_${ref2%.f*}.sam
				awk -v ref3=${ref3%.f*} -F '\t' '$3 ~ ref3{print $1}' ${projdir}/preprocess/${i%.f*}_del0_${ref1%.f*}_${ref2%.f*}.sam | \
				grep -vFwf - ${projdir}/preprocess/${i%.f*}_del0_${ref1%.f*}_${ref2%.f*}.sam | \
				cat <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz | grep '^@') - > ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}_${ref2%.f*}.sam
				rm ${projdir}/preprocess/${i%.f*}_del0_${ref1%.f*}_${ref2%.f*}.sam
        awk '{if($1==2) print $0}' ${projdir}/preprocess/${i%.f*}_Index_subgenome.txt | awk -v ref1=${ref1%.f*} -v ref2=${ref3%.f*} '{if($3 == ref1 || $3 == ref3) print $2}' | grep -Fw -f - <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz) > ${projdir}/preprocess/${i%.f*}_del0_${ref1%.f*}_${ref3%.f*}.sam
				awk -v ref2=${ref2%.f*} -F '\t' '$3 ~ ref2{print $1}' ${projdir}/preprocess/${i%.f*}_del0_${ref1%.f*}_${ref3%.f*}.sam | \
				grep -vFwf - ${projdir}/preprocess/${i%.f*}_del0_${ref1%.f*}_${ref3%.f*}.sam | \
				cat <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz | grep '^@') - > ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}_${ref3%.f*}.sam
				rm ${projdir}/preprocess/${i%.f*}_del0_${ref1%.f*}_${ref3%.f*}.sam
        awk '{if($1==2) print $0}' ${projdir}/preprocess/${i%.f*}_Index_subgenome.txt | awk -v ref1=${ref2%.f*} -v ref2=${ref3%.f*} '{if($3 == ref2 || $3 == ref3) print $2}' | grep -Fw -f - <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz) > ${projdir}/preprocess/${i%.f*}_del0_${ref2%.f*}_${ref3%.f*}.sam
				awk -v ref1=${ref1%.f*} -F '\t' '$3 ~ ref1{print $1}' ${projdir}/preprocess/${i%.f*}_del0_${ref2%.f*}_${ref3%.f*}.sam | \
				grep -vFwf - ${projdir}/preprocess/${i%.f*}_del0_${ref2%.f*}_${ref3%.f*}.sam | \
				cat <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz | grep '^@') - > ${projdir}/preprocess/${i%.f*}_del_${ref2%.f*}_${ref3%.f*}.sam
				rm ${projdir}/preprocess/${i%.f*}_del0_${ref2%.f*}_${ref3%.f*}.sam

        awk -v ref1=${ref1%.f*} '{if($1==1 && $3 == ref1) print $2}' ${projdir}/preprocess/${i%.f*}_Index_subgenome.txt | grep -Fw -f - <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz) | cat <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz | grep '^@') - > ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}.sam
        awk -v ref2=${ref2%.f*} '{if($1==1 && $3 == ref2) print $2}' ${projdir}/preprocess/${i%.f*}_Index_subgenome.txt | grep -Fw -f - <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz) | cat <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz | grep '^@') - > ${projdir}/preprocess/${i%.f*}_del_${ref2%.f*}.sam
        awk -v ref3=${ref3%.f*} '{if($1==1 && $3 == ref3) print $2}' ${projdir}/preprocess/${i%.f*}_Index_subgenome.txt | grep -Fw -f - <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz) | cat <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz | grep '^@') - > ${projdir}/preprocess/${i%.f*}_del_${ref3%.f*}.sam
				rm ${projdir}/preprocess/${i%.f*}_Index0_subgenome.txt ${projdir}/preprocess/${i%.f*}_Index_subgenome.txt

        echo "nloci~POS~mapQ~CHROM~CHROM" > ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.txt
				grep -v '^@' ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam | awk '($3 != "\*")' | awk '{$2=$2"_"$3}1' | awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | \
				cat <(grep '^@' ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam) - > ${projdir}/alignment_summaries/temp_${i%.f*}.txt
				awk -F '\t' '{print $1"\t"$3"\t"$4"\t"$5}' ${projdir}/alignment_summaries/temp_${i%.f*}.txt | awk '{gsub(/ /,"\t"); print}' | awk '{print $1"~"$5"~"$3"~"$4}' | grep -v '\*' | grep -v '^@' >> ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.txt
				awk '!visited[$0]++' ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.txt > ${projdir}/alignment_summaries/temp_${i%.f*}.txt
				mv ${projdir}/alignment_summaries/temp_${i%.f*}.txt ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.txt

        echo "nloci~POS~mapQ~CHROM~CHROM" > ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}_${ref2%.f*}.txt
				grep -v '^@' ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}_${ref2%.f*}.sam | awk '($3 != "\*")' | awk '{$2=$2"_"$3}1' | awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | \
				cat <(grep '^@' ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}_${ref2%.f*}.sam) - > ${projdir}/alignment_summaries/temp_${i%.f*}.txt
				awk -F '\t' '{print $1"\t"$3"\t"$4"\t"$5}' ${projdir}/alignment_summaries/temp_${i%.f*}.txt | awk '{gsub(/ /,"\t"); print}' | awk '{print $1"~"$5"~"$3"~"$4}' | grep -v '\*' | grep -v '^@' >> ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}_${ref2%.f*}.txt
				awk '!visited[$0]++' ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}_${ref2%.f*}.txt > ${projdir}/alignment_summaries/temp_${i%.f*}.txt
				mv ${projdir}/alignment_summaries/temp_${i%.f*}.txt ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}_${ref2%.f*}.txt
        echo "nloci~POS~mapQ~CHROM~CHROM" > ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}_${ref3%.f*}.txt
        grep -v '^@' ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}_${ref3%.f*}.sam | awk '($3 != "\*")' | awk '{$2=$2"_"$3}1' | awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | \
        cat <(grep '^@' ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}_${ref3%.f*}.sam) - > ${projdir}/alignment_summaries/temp_${i%.f*}.txt
        awk -F '\t' '{print $1"\t"$3"\t"$4"\t"$5}' ${projdir}/alignment_summaries/temp_${i%.f*}.txt | awk '{gsub(/ /,"\t"); print}' | awk '{print $1"~"$5"~"$3"~"$4}' | grep -v '\*' | grep -v '^@' >> ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}_${ref3%.f*}.txt
        awk '!visited[$0]++' ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}_${ref3%.f*}.txt > ${projdir}/alignment_summaries/temp_${i%.f*}.txt
        mv ${projdir}/alignment_summaries/temp_${i%.f*}.txt ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}_${ref3%.f*}.txt
        echo "nloci~POS~mapQ~CHROM~CHROM" > ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref2%.f*}_${ref3%.f*}.txt
        grep -v '^@' ${projdir}/preprocess/${i%.f*}_del_${ref2%.f*}_${ref3%.f*}.sam | awk '($3 != "\*")' | awk '{$2=$2"_"$3}1' | awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | \
        cat <(grep '^@' ${projdir}/preprocess/${i%.f*}_del_${ref2%.f*}_${ref3%.f*}.sam) - > ${projdir}/alignment_summaries/temp_${i%.f*}.txt
        awk -F '\t' '{print $1"\t"$3"\t"$4"\t"$5}' ${projdir}/alignment_summaries/temp_${i%.f*}.txt | awk '{gsub(/ /,"\t"); print}' | awk '{print $1"~"$5"~"$3"~"$4}' | grep -v '\*' | grep -v '^@' >> ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref2%.f*}_${ref3%.f*}.txt
        awk '!visited[$0]++' ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref2%.f*}_${ref3%.f*}.txt > ${projdir}/alignment_summaries/temp_${i%.f*}.txt
        mv ${projdir}/alignment_summaries/temp_${i%.f*}.txt ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref2%.f*}_${ref3%.f*}.txt

				echo "nloci~POS~mapQ~CHROM~CHROM" > ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}.txt
				grep -v '^@' ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}.sam | awk '($3 != "\*")' | awk '{$2=$2"_"$3}1' | awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | \
				cat <(grep '^@' ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}.sam) - > ${projdir}/alignment_summaries/temp_${i%.f*}.txt
				awk -F '\t' '{print $1"\t"$3"\t"$4"\t"$5}' ${projdir}/alignment_summaries/temp_${i%.f*}.txt | awk '{gsub(/ /,"\t"); print}' | awk '{print $1"~"$5"~"$3"~"$4}' | grep -v '\*' | grep -v '^@'  >> ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}.txt
				awk '!visited[$0]++' ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}.txt > ${projdir}/alignment_summaries/temp_${i%.f*}.txt
				mv ${projdir}/alignment_summaries/temp_${i%.f*}.txt ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}.txt
				echo "nloci~POS~mapQ~CHROM~CHROM" > ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref2%.f*}.txt
				grep -v '^@' ${projdir}/preprocess/${i%.f*}_del_${ref2%.f*}.sam | awk '($3 != "\*")' | awk '{$2=$2"_"$3}1' | awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | \
				cat <(grep '^@' ${projdir}/preprocess/${i%.f*}_del_${ref2%.f*}.sam) - > ${projdir}/alignment_summaries/temp_${i%.f*}.txt
				awk -F '\t' '{print $1"\t"$3"\t"$4"\t"$5}' ${projdir}/alignment_summaries/temp_${i%.f*}.txt | awk '{gsub(/ /,"\t"); print}' | awk '{print $1"~"$5"~"$3"~"$4}' | grep -v '\*' | grep -v '^@'  >> ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref2%.f*}.txt
				awk '!visited[$0]++' ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref2%.f*}.txt > ${projdir}/alignment_summaries/temp_${i%.f*}.txt
				mv ${projdir}/alignment_summaries/temp_${i%.f*}.txt ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref2%.f*}.txt
				wait
        echo "nloci~POS~mapQ~CHROM~CHROM" > ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref3%.f*}.txt
				grep -v '^@' ${projdir}/preprocess/${i%.f*}_del_${ref3%.f*}.sam | awk '($3 != "\*")' | awk '{$2=$2"_"$3}1' | awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | \
				cat <(grep '^@' ${projdir}/preprocess/${i%.f*}_del_${ref3%.f*}.sam) - > ${projdir}/alignment_summaries/temp_${i%.f*}.txt
				awk -F '\t' '{print $1"\t"$3"\t"$4"\t"$5}' ${projdir}/alignment_summaries/temp_${i%.f*}.txt | awk '{gsub(/ /,"\t"); print}' | awk '{print $1"~"$5"~"$3"~"$4}' | grep -v '\*' | grep -v '^@'  >> ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref3%.f*}.txt
				awk '!visited[$0]++' ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}.txt > ${projdir}/alignment_summaries/temp_${i%.f*}.txt
				mv ${projdir}/alignment_summaries/temp_${i%.f*}.txt ${projdir}/alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}.txt
        wait

				{
	        awk '/@HD/ || /@SQ/{print}' <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz) > ${projdir}/preprocess/${i%.f*}_heading_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam
					grep -v '^@' ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam | awk '($3 != "\*")' | awk '{gsub(/_se-/,"_se-\t",$1); gsub(/_pe-/,"_pe-\t",$1)}1' | \
					awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | awk -v multilocus=$multilocus -F '\t' 'BEGIN{OFS="\t"} {if ($5==multilocus) {$5=$5+40}}1' > ${projdir}/preprocess/${i%.f*}_uniq_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam
					for j in $(LC_ALL=C; sort -n -k1,1 ${projdir}/preprocess/${i%.f*}_uniq_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam | awk '{print $1}' | uniq); do
					  awk -v n="^${j}" '$0~n{print $0}' ${projdir}/preprocess/${i%.f*}_uniq_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam | awk -v n="$j" '{for(i=0;i<n;i++) print}' >> ${projdir}/preprocess/${i%.f*}_exp_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam
					done; wait
					awk '{print "seq"NR"_"$0}' ${projdir}/preprocess/${i%.f*}_exp_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam | tr -s ' ' | awk '{gsub(/256/,"0",$2); gsub(/272/,"16",$2); print}' | \
					awk '($3 != "\*")' | awk '($6 != "\*")'  | cat ${projdir}/preprocess/${i%.f*}_heading_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam - | tr ' ' '\t' > ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam
					rm ${projdir}/preprocess/${i%.f*}_exp_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam ${projdir}/preprocess/${i%.f*}_uniq_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam ${projdir}/preprocess/${i%.f*}_heading_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam
				} & PIDexp1=$!

				{
	        awk '/@HD/ || /@SQ/{print}' <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz) > ${projdir}/preprocess/${i%.f*}_heading_${ref1%.f*}_${ref2%.f*}.sam
	        grep -v '^@' ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}_${ref2%.f*}.sam | awk '($3 != "\*")' | awk '{gsub(/_se-/,"_se-\t",$1); gsub(/_pe-/,"_pe-\t",$1)}1' | \
	        awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | awk -v multilocus=$multilocus -F '\t' 'BEGIN{OFS="\t"} {if ($5==multilocus) {$5=$5+40}}1' > ${projdir}/preprocess/${i%.f*}_uniq_${ref1%.f*}_${ref2%.f*}.sam
	        for j in $(LC_ALL=C; sort -n -k1,1 ${projdir}/preprocess/${i%.f*}_uniq_${ref1%.f*}_${ref2%.f*}.sam | awk '{print $1}' | uniq); do
	        	awk -v n="^${j}" '$0~n{print $0}' ${projdir}/preprocess/${i%.f*}_uniq_${ref1%.f*}_${ref2%.f*}.sam | awk -v n="$j" '{for(i=0;i<n;i++) print}' >> ${projdir}/preprocess/${i%.f*}_exp_${ref1%.f*}_${ref2%.f*}.sam
	        done; wait
	        awk '{print "seq"NR"_"$0}' ${projdir}/preprocess/${i%.f*}_exp_${ref1%.f*}_${ref2%.f*}.sam | tr -s ' ' | awk '{gsub(/256/,"0",$2); gsub(/272/,"16",$2); print}' | \
	        awk '($3 != "\*")' | awk '($6 != "\*")'  | cat ${projdir}/preprocess/${i%.f*}_heading_${ref1%.f*}_${ref2%.f*}.sam - | tr ' ' '\t' > ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_${ref2%.f*}.sam
	        rm ${projdir}/preprocess/${i%.f*}_exp_${ref1%.f*}_${ref2%.f*}.sam ${projdir}/preprocess/${i%.f*}_uniq_${ref1%.f*}_${ref2%.f*}.sam ${projdir}/preprocess/${i%.f*}_heading_${ref1%.f*}_${ref2%.f*}.sam ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}_${ref2%.f*}.sam
				} & PIDexp2=$!

				{
	        awk '/@HD/ || /@SQ/{print}' <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz) > ${projdir}/preprocess/${i%.f*}_heading_${ref1%.f*}_${ref3%.f*}.sam
	        grep -v '^@' ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}_${ref3%.f*}.sam | awk '($3 != "\*")' | awk '{gsub(/_se-/,"_se-\t",$1); gsub(/_pe-/,"_pe-\t",$1)}1' | \
	        awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | awk -v multilocus=$multilocus -F '\t' 'BEGIN{OFS="\t"} {if ($5==multilocus) {$5=$5+40}}1' > ${projdir}/preprocess/${i%.f*}_uniq_${ref1%.f*}_${ref3%.f*}.sam
	        for j in $(LC_ALL=C; sort -n -k1,1 ${projdir}/preprocess/${i%.f*}_uniq_${ref1%.f*}_${ref3%.f*}.sam | awk '{print $1}' | uniq); do
	        	awk -v n="^${j}" '$0~n{print $0}' ${projdir}/preprocess/${i%.f*}_uniq_${ref1%.f*}_${ref3%.f*}.sam | awk -v n="$j" '{for(i=0;i<n;i++) print}' >> ${projdir}/preprocess/${i%.f*}_exp_${ref1%.f*}_${ref3%.f*}.sam
	        done; wait
	        awk '{print "seq"NR"_"$0}' ${projdir}/preprocess/${i%.f*}_exp_${ref1%.f*}_${ref3%.f*}.sam | tr -s ' ' | awk '{gsub(/256/,"0",$2); gsub(/272/,"16",$2); print}' | \
	        awk '($3 != "\*")' | awk '($6 != "\*")'  | cat ${projdir}/preprocess/${i%.f*}_heading_${ref1%.f*}_${ref3%.f*}.sam - | tr ' ' '\t' > ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_${ref3%.f*}.sam
	        rm ${projdir}/preprocess/${i%.f*}_exp_${ref1%.f*}_${ref3%.f*}.sam ${projdir}/preprocess/${i%.f*}_uniq_${ref1%.f*}_${ref3%.f*}.sam ${projdir}/preprocess/${i%.f*}_heading_${ref1%.f*}_${ref3%.f*}.sam ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}_${ref3%.f*}.sam
				} & PIDexp3=$!

				{
	        awk '/@HD/ || /@SQ/{print}' <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz) > ${projdir}/preprocess/${i%.f*}_heading_${ref2%.f*}_${ref3%.f*}.sam
	        grep -v '^@' ${projdir}/preprocess/${i%.f*}_del_${ref2%.f*}_${ref3%.f*}.sam | awk '($3 != "\*")' | awk '{gsub(/_se-/,"_se-\t",$1); gsub(/_pe-/,"_pe-\t",$1)}1' | \
	        awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | awk -v multilocus=$multilocus -F '\t' 'BEGIN{OFS="\t"} {if ($5==multilocus) {$5=$5+40}}1' > ${projdir}/preprocess/${i%.f*}_uniq_${ref2%.f*}_${ref3%.f*}.sam
	        for j in $(LC_ALL=C; sort -n -k1,1 ${projdir}/preprocess/${i%.f*}_uniq_${ref2%.f*}_${ref3%.f*}.sam | awk '{print $1}' | uniq); do
	        	awk -v n="^${j}" '$0~n{print $0}' ${projdir}/preprocess/${i%.f*}_uniq_${ref2%.f*}_${ref3%.f*}.sam | awk -v n="$j" '{for(i=0;i<n;i++) print}' >> ${projdir}/preprocess/${i%.f*}_exp_${ref2%.f*}_${ref3%.f*}.sam
	        done; wait
	        awk '{print "seq"NR"_"$0}' ${projdir}/preprocess/${i%.f*}_exp_${ref2%.f*}_${ref3%.f*}.sam | tr -s ' ' | awk '{gsub(/256/,"0",$2); gsub(/272/,"16",$2); print}' | \
	        awk '($3 != "\*")' | awk '($6 != "\*")'  | cat ${projdir}/preprocess/${i%.f*}_heading_${ref2%.f*}_${ref3%.f*}.sam - | tr ' ' '\t' > ${projdir}/preprocess/${i%.f*}_${ref2%.f*}_${ref3%.f*}.sam
	        rm ${projdir}/preprocess/${i%.f*}_exp_${ref2%.f*}_${ref3%.f*}.sam ${projdir}/preprocess/${i%.f*}_uniq_${ref2%.f*}_${ref3%.f*}.sam ${projdir}/preprocess/${i%.f*}_heading_${ref2%.f*}_${ref3%.f*}.sam ${projdir}/preprocess/${i%.f*}_del_${ref2%.f*}_${ref3%.f*}.sam
				} & PIDexp4=$!

				wait $PIDexp1
				wait $PIDexp2
				wait $PIDexp3
				wait $PIDexp4
				wait

				{
					awk '/@HD/ || /@SQ/{print}' <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz) > ${projdir}/preprocess/${i%.f*}_heading_${ref1%.f*}.sam
					grep -v '^@' ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}.sam | awk '($3 != "\*")' | awk '{gsub(/_se-/,"_se-\t",$1); gsub(/_pe-/,"_pe-\t",$1)}1' | \
					awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | awk -v multilocus=$multilocus -F '\t' 'BEGIN{OFS="\t"} {if ($5==multilocus) {$5=$5+40}}1' > ${projdir}/preprocess/${i%.f*}_uniq_${ref1%.f*}.sam
					for j in $(LC_ALL=C; sort -n -k1,1 ${projdir}/preprocess/${i%.f*}_uniq_${ref1%.f*}.sam | awk '{print $1}' | uniq); do
					  awk -v n="^${j}" '$0~n{print $0}' ${projdir}/preprocess/${i%.f*}_uniq_${ref1%.f*}.sam | awk -v n="$j" '{for(i=0;i<n;i++) print}' >> ${projdir}/preprocess/${i%.f*}_exp_${ref1%.f*}.sam
					done; wait
					awk '{print "seq"NR"_"$0}' ${projdir}/preprocess/${i%.f*}_exp_${ref1%.f*}.sam | tr -s ' ' | awk '{gsub(/256/,"0",$2); gsub(/272/,"16",$2); print}' | \
					awk '($3 != "\*")' | awk '($6 != "\*")'  | cat ${projdir}/preprocess/${i%.f*}_heading_${ref1%.f*}.sam - | tr ' ' '\t' > ${projdir}/preprocess/${i%.f*}_${ref1%.f*}.sam
					rm ${projdir}/preprocess/${i%.f*}_exp_${ref1%.f*}.sam ${projdir}/preprocess/${i%.f*}_uniq_${ref1%.f*}.sam ${projdir}/preprocess/${i%.f*}_heading_${ref1%.f*}.sam ${projdir}/preprocess/${i%.f*}_del_${ref1%.f*}.sam
				} & PIDexp5=$!

				{
					awk '/@HD/ || /@SQ/{print}' <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz) > ${projdir}/preprocess/${i%.f*}_heading_${ref2%.f*}.sam
					grep -v '^@' ${projdir}/preprocess/${i%.f*}_del_${ref2%.f*}.sam | awk '($3 != "\*")' | awk '{gsub(/_se-/,"_se-\t",$1); gsub(/_pe-/,"_pe-\t",$1)}1' | \
					awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | awk -v multilocus=$multilocus -F '\t' 'BEGIN{OFS="\t"} {if ($5==multilocus) {$5=$5+40}}1' > ${projdir}/preprocess/${i%.f*}_uniq_${ref2%.f*}.sam
					for j in $(LC_ALL=C; sort -n -k1,1 ${projdir}/preprocess/${i%.f*}_uniq_${ref2%.f*}.sam | awk '{print $1}' | uniq); do
					  awk -v n="^${j}" '$0~n{print $0}' ${projdir}/preprocess/${i%.f*}_uniq_${ref2%.f*}.sam | awk -v n="$j" '{for(i=0;i<n;i++) print}' >> ${projdir}/preprocess/${i%.f*}_exp_${ref2%.f*}.sam
					done; wait
					awk '{print "seq"NR"_"$0}' ${projdir}/preprocess/${i%.f*}_exp_${ref2%.f*}.sam | tr -s ' ' | awk '{gsub(/256/,"0",$2); gsub(/272/,"16",$2); print}' | \
					awk '($3 != "\*")' | awk '($6 != "\*")'  | cat ${projdir}/preprocess/${i%.f*}_heading_${ref2%.f*}.sam - | tr ' ' '\t' > ${projdir}/preprocess/${i%.f*}_${ref2%.f*}.sam
					rm ${projdir}/preprocess/${i%.f*}_exp_${ref2%.f*}.sam ${projdir}/preprocess/${i%.f*}_uniq_${ref2%.f*}.sam ${projdir}/preprocess/${i%.f*}_heading_${ref2%.f*}.sam ${projdir}/preprocess/${i%.f*}_del_${ref2%.f*}.sam
				} & PIDexp6=$!

				{
	        awk '/@HD/ || /@SQ/{print}' <(zcat ${projdir}/preprocess/${i%.f*}_redun.sam.gz) > ${projdir}/preprocess/${i%.f*}_heading_${ref3%.f*}.sam
					grep -v '^@' ${projdir}/preprocess/${i%.f*}_del_${ref3%.f*}.sam | awk '($3 != "\*")' | awk '{gsub(/_se-/,"_se-\t",$1); gsub(/_pe-/,"_pe-\t",$1)}1' | \
					awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | awk -v multilocus=$multilocus -F '\t' 'BEGIN{OFS="\t"} {if ($5==multilocus) {$5=$5+40}}1' > ${projdir}/preprocess/${i%.f*}_uniq_${ref3%.f*}.sam
					for j in $(LC_ALL=C; sort -n -k1,1 ${projdir}/preprocess/${i%.f*}_uniq_${ref3%.f*}.sam | awk '{print $1}' | uniq); do
					  awk -v n="^${j}" '$0~n{print $0}' ${projdir}/preprocess/${i%.f*}_uniq_${ref3%.f*}.sam | awk -v n="$j" '{for(i=0;i<n;i++) print}' >> ${projdir}/preprocess/${i%.f*}_exp_${ref3%.f*}.sam
					done; wait
					awk '{print "seq"NR"_"$0}' ${projdir}/preprocess/${i%.f*}_exp_${ref3%.f*}.sam | tr -s ' ' | awk '{gsub(/256/,"0",$2); gsub(/272/,"16",$2); print}' | \
					awk '($3 != "\*")' | awk '($6 != "\*")'  | cat ${projdir}/preprocess/${i%.f*}_heading_${ref3%.f*}.sam - | tr ' ' '\t' > ${projdir}/preprocess/${i%.f*}_${ref3%.f*}.sam
					rm ${projdir}/preprocess/${i%.f*}_exp_${ref3%.f*}.sam ${projdir}/preprocess/${i%.f*}_uniq_${ref3%.f*}.sam ${projdir}/preprocess/${i%.f*}_heading_${ref3%.f*}.sam ${projdir}/preprocess/${i%.f*}_del_${ref3%.f*}.sam
				} & PIDexp7=$!

				wait $PIDexp5
				wait $PIDexp6
				wait $PIDexp7
				wait
				rm ${projdir}/preprocess/${i%.f*}_redun.sam.gz

				declare -a arr=("${i%.f*}_${ref1%.f*}.sam" "${i%.f*}_${ref2%.f*}.sam" "${i%.f*}_${ref3%.f*}.sam" \
        "${i%.f*}_${ref1%.f*}_${ref2%.f*}.sam" "${i%.f*}_${ref1%.f*}_${ref3%.f*}.sam" "${i%.f*}_${ref2%.f*}_${ref3%.f*}.sam" \
        "${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam" )  && \
			  cd ${projdir}/preprocess
			  for j in "${arr[@]}"; do
			        $java $Xmx2 -XX:ParallelGCThreads=$gthreads -jar $picard SortSam I=$j O=${j%.sam*}.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT && \
			        $java $Xmx2 -XX:ParallelGCThreads=$gthreads -jar $picard BuildBamIndex INPUT=${j%.sam*}.bam VALIDATION_STRINGENCY=LENIENT && \
			        $java $Xmx2 -XX:ParallelGCThreads=$gthreads -jar $picard AddOrReplaceReadGroups I=${j%.sam*}.bam O=${j%.sam*}_precall.bam RGLB=${i%.f*} RGPL=illumina RGPU=run RGSM=${i%.f*} VALIDATION_STRINGENCY=LENIENT && \
			        $samtools index ${j%.sam*}_precall.bam &&
							rm $j ${j%.sam*}.bam ${j%.sam*}.bai &&
							wait
			  done
				wait
			fi
		fi

	done
	wait && touch ${projdir}/precall_done.txt
	ls ${projdir}/preprocess/* | grep -v precall | xargs rm &&
	cd ${projdir}/samples

	if [[ ! -f "${projdir}/align3_${samples_list}" ]]; then touch "${projdir}/align3_${samples_list}"; fi
	if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
		align=$(ls ${projdir}/align3_samples_list_node_* | wc -l)
		while [[ "$align" -lt $nodes ]]; do sleep 30; align=$(ls ${projdir}/align3_samples_list_node_* | wc -l); done
		if [[ $align == $nodes ]] && test ! -f ${projdir}/alignment_summaries/refgenome_paralogs.txt; then
			rm ${projdir}/align3_${samples_list}
			cd ${projdir}/alignment_summaries
			cat ${projdir}/alignment_summaries/*_summ.txt > ${projdir}/alignment_summaries/alignment_summaries_unique_reads.txt; rm -r ${projdir}/alignment_summaries/*_summ.txt &&
			# Total number of reads per samples
			awk '/###---/ || /QC-passed/{print}' ${projdir}/alignment_summaries/alignment_summaries_unique_reads.txt | cut -d\+ -f1 | tr -d '\n' | \
			awk  'gsub(/---###/, "\t", $0)' | awk  'gsub(/###---/, "", $0)' | tr ' ' '\n' > ${projdir}/alignment_summaries/total_unique_reads.txt &&
			# Total number of mapped reads per samples
			cat ${projdir}/alignment_summaries/alignment_summaries_unique_reads.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /0_mapped/{print}' |\
			tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
			awk 'gsub("\\+0mapped", "\t", $0)' | tr ":" "\t" | cut -d\: -f1 | awk 'gsub(/ /, "\t")' > ${projdir}/alignment_summaries/total_unique_reads_mapped.txt &&
			# Total number of mapped reads per samples
			cat ${projdir}/alignment_summaries/alignment_summaries_unique_reads.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /properly_paired/{print}' |\
			tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
			awk 'gsub("\\+0properlypaired", "\t", $0)' | tr ":" "\t" | cut -d\: -f1 | awk 'gsub(/ /, "\t")' > ${projdir}/alignment_summaries/total_unique_reads_paired.txt &&
			echo -e "Samples\tTotal\tMapped\tPerc_Mapped\tPE_Mapped\t%_PE_Mapped" > ${projdir}/alignment_summaries/summary_precall.txt &&
			awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' ${projdir}/alignment_summaries/total_unique_reads_mapped.txt  ${projdir}/alignment_summaries/total_unique_reads.txt  | \
			awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' ${projdir}/alignment_summaries/total_unique_reads_paired.txt - | awk '{gsub(/ /,"\t"); print $0}' | \
			cat ${projdir}/alignment_summaries/summary_precall.txt - > ${projdir}/alignment_summaries/Tabulated_Alignment_Unique_Read_Summaries.txt &&
			awk '{gsub(/\t/,","); print $0}' ${projdir}/alignment_summaries/Tabulated_Alignment_Unique_Read_Summaries.txt > ${projdir}/alignment_summaries/Tabulated_Alignment_Unique_Read_Summaries.csv &&
			rm ${projdir}/alignment_summaries/total_unique_* ${projdir}/alignment_summaries/summary_precall.txt &> /dev/null &&
			rm ${projdir}/samples/metrics.txt ${projdir}/preprocess/metrics.txt &> /dev/null &&

			cd $projdir/alignment_summaries


      touch refgenome_paralogs_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.txt
      for par in refgenome_paralogs_*_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.txt; do
      	cat refgenome_paralogs_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.txt $par | awk '!visited[$0]++' > temp_par.txt &&
      	mv temp_par.txt refgenome_paralogs_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.txt &&
      	wait
      done
      wait
      rm refgenome_paralogs_*_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.txt
      awk '{gsub(/~/,"\t"); print $0}' refgenome_paralogs_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.txt > temp.txt &&
      awk '{print $4"\t"$2"\t"$1}' temp.txt | awk '!/CHROM/' | cat <(printf "CHROM\tPOS\tnloci\n") - > refgenome_paralogs_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.txt &&

      touch refgenome_paralogs_${ref1%.f*}_${ref2%.f*}.txt
      for par in refgenome_paralogs_*_${ref1%.f*}_${ref2%.f*}.txt; do
      	cat refgenome_paralogs_${ref1%.f*}_${ref2%.f*}.txt $par | awk '!visited[$0]++' > temp_par.txt &&
      	mv temp_par.txt refgenome_paralogs_${ref1%.f*}_${ref2%.f*}.txt &&
      	wait
      done
      wait
      rm refgenome_paralogs_*_${ref1%.f*}_${ref2%.f*}.txt
      awk '{gsub(/~/,"\t"); print $0}' refgenome_paralogs_${ref1%.f*}_${ref2%.f*}.txt > temp.txt &&
      awk '{print $4"\t"$2"\t"$1}' temp.txt | awk '!/CHROM/' | cat <(printf "CHROM\tPOS\tnloci\n") - > refgenome_paralogs_${ref1%.f*}_${ref2%.f*}.txt &&

      touch refgenome_paralogs_${ref1%.f*}_${ref3%.f*}.txt
      for par in refgenome_paralogs_*_${ref1%.f*}_${ref3%.f*}.txt; do
      	cat refgenome_paralogs_${ref1%.f*}_${ref3%.f*}.txt $par | awk '!visited[$0]++' > temp_par.txt &&
      	mv temp_par.txt refgenome_paralogs_${ref1%.f*}_${ref3%.f*}.txt &&
      	wait
      done
      wait
      rm refgenome_paralogs_*_${ref1%.f*}_${ref3%.f*}.txt
      awk '{gsub(/~/,"\t"); print $0}' refgenome_paralogs_${ref1%.f*}_${ref3%.f*}.txt > temp.txt &&
      awk '{print $4"\t"$2"\t"$1}' temp.txt | awk '!/CHROM/' | cat <(printf "CHROM\tPOS\tnloci\n") - > refgenome_paralogs_${ref1%.f*}_${ref3%.f*}.txt &&

      touch refgenome_paralogs_${ref2%.f*}_${ref3%.f*}.txt
      for par in refgenome_paralogs_*_${ref2%.f*}_${ref3%.f*}.txt; do
      	cat refgenome_paralogs_${ref2%.f*}_${ref3%.f*}.txt $par | awk '!visited[$0]++' > temp_par.txt &&
      	mv temp_par.txt refgenome_paralogs_${ref2%.f*}_${ref3%.f*}.txt &&
      	wait
      done
      wait
      rm refgenome_paralogs_*_${ref2%.f*}_${ref3%.f*}.txt
      awk '{gsub(/~/,"\t"); print $0}' refgenome_paralogs_${ref2%.f*}_${ref3%.f*}.txt > temp.txt &&
      awk '{print $4"\t"$2"\t"$1}' temp.txt | awk '!/CHROM/' | cat <(printf "CHROM\tPOS\tnloci\n") - > refgenome_paralogs_${ref2%.f*}_${ref3%.f*}.txt &&

			touch refgenome_paralogs_${ref1%.f*}.txt
			for par in refgenome_paralogs_*_${ref1%.f*}.txt; do
				cat refgenome_paralogs_${ref1%.f*}.txt $par | awk '!visited[$0]++' > temp_par.txt &&
				mv temp_par.txt refgenome_paralogs_${ref1%.f*}.txt &&
				wait
			done
			wait
			rm refgenome_paralogs_*_${ref1%.f*}.txt
			awk '{gsub(/~/,"\t"); print $0}' refgenome_paralogs_${ref1%.f*}.txt > temp.txt &&
			awk '{print $4"\t"$2"\t"$1}' temp.txt | awk '!/CHROM/' | cat <(printf "CHROM\tPOS\tnloci\n") - > refgenome_paralogs_${ref1%.f*}.txt &&

			touch refgenome_paralogs_${ref2%.f*}.txt
			for par in $(ls refgenome_paralogs_*_${ref2%.f*}.txt | grep -v ${ref1%.f*}); do
				cat refgenome_paralogs_${ref2%.f*}.txt $par | awk '!visited[$0]++' > temp_par.txt &&
				mv temp_par.txt refgenome_paralogs_${ref2%.f*}.txt &&
				wait
			done
			wait
			rm refgenome_paralogs_*_${ref2%.f*}.txt
			awk '{gsub(/~/,"\t"); print $0}' refgenome_paralogs_${ref2%.f*}.txt > temp.txt &&
			awk '{print $4"\t"$2"\t"$1}' temp.txt | awk '!/CHROM/' | cat <(printf "CHROM\tPOS\tnloci\n") - > refgenome_paralogs_${ref2%.f*}.txt &&

      touch refgenome_paralogs_${ref3%.f*}.txt
			for par in refgenome_paralogs_*_${ref3%.f*}.txt; do
				cat refgenome_paralogs_${ref3%.f*}.txt $par | awk '!visited[$0]++' > temp_par.txt &&
				mv temp_par.txt refgenome_paralogs_${ref3%.f*}.txt &&
				wait
			done
			wait
			rm refgenome_paralogs_*_${ref3%.f*}.txt
			awk '{gsub(/~/,"\t"); print $0}' refgenome_paralogs_${ref3%.f*}.txt > temp.txt &&
			awk '{print $4"\t"$2"\t"$1}' temp.txt | awk '!/CHROM/' | cat <(printf "CHROM\tPOS\tnloci\n") - > refgenome_paralogs_${ref3%.f*}.txt &&

			awk 'FNR==1 && NR!=1 { while (/^CHROM/) getline; }1 {print}' refgenome_paralogs_*.txt > refgenome_paralogs.txt &&
			rm refgenome_paralogs_* temp.txt
		fi
	fi

	wait && touch ${projdir}/alignment_done.txt
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



######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${yellow}\n- Performing Variant Calling with GATK HaplotypeCaller\n${blue}##############################################################################${white}\n"
main () {
	cd $projdir
	cd preprocess
	mkdir -p processed

	if [[ ! -f "${projdir}/call0_${samples_list}" ]]; then touch "${projdir}/call0_${samples_list}"; fi
	if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
		align=$(ls ${projdir}/call0_samples_list_node_* | wc -l)
		while [[ "$align" -lt $nodes ]]; do sleep 30; align=$(ls ${projdir}/call0_samples_list_node_* | wc -l); done
		if [[ $align == $nodes ]]; then
			rm ${projdir}/call0_${samples_list}
			if [[ -z "$(ls -A ./processed/ &> /dev/null)" ]]; then
				:
			else
				mv $projdir/preprocess/processed/*_precall.bam* $projdir/preprocess/
			fi
		fi
	fi


	if [[ "$joint_calling" == true ]]; then

    echo -e "${magenta}- performing SNP calling across entire genome (subgenome 1, 2, and 3) ${white}\n"

    j=-I; input=""; k=""
    for i in $(ls *_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_precall.bam); do
      k="${j} ${i}"; input="${input} ${k}"
    done
    Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' )
    if [[ ! -z "$Exclude_Chromosome" ]]; then
      for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
        Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
      done
    fi

    if [[ -z "$Get_Chromosome" ]]; then
      for selchr in $Get2_Chromosome; do (
        if test ! -f ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_${selchr}_raw.vcf; then
          $GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${selchr} ${input} -ploidy $ploidy -O ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_${selchr}_raw.vcf.gz --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) &&
          gunzip ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_${selchr}_raw.vcf.gz &&
          wait
        fi
        ) &
        if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
          wait
        fi
      done
      wait
      cd ../snpcall
      grep -h '^#' ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt
      cat ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_*_raw.vcf | awk '!/^#/' > all.vcf
      cat vcf_header.txt all.vcf > ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw.vcf
      rm vcf_header.txt all.vcf *.vcf.gz.tbi ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_*_raw.vcf
    else
      echo $Get_Chromosome | tr ',' '\n' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2}' > ${projdir}/refgenomes/panref.list
      $GATK --java-options "$Xmx2 -XX:+UseParallelGC -XX:ParallelGCThreads=$threads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${projdir}/refgenomes/panref.list ${input}-ploidy $ploidy -O ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw.vcf.gz --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) &&
      cd ../snpcall
      gunzip ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw.vcf &&
      wait
    fi

    cd ${projdir}/preprocess
    mv *_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_precall* ./processed/

    ######################
    echo -e "${magenta}- performing SNP calling across entire genome (subgenome 1 and 2) ${white}\n"

    j=-I; input=""; k=""
    for i in $(ls *_${ref1%.f*}_${ref2%.f*}_precall.bam); do
    	k="${j} ${i}"; input="${input} ${k}"
    done
    calcploidy=$(($ploidy_ref1 + $ploidy_ref2))
    Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' | awk -v pat1=${ref1%.f*} -v pat2=${ref2%.f*} 'pat1 || pat2')
    if [[ ! -z "$Exclude_Chromosome" ]]; then
    	for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
    		Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
    	done
    fi

    if [[ -z "$Get_Chromosome" ]]; then
    	for selchr in $Get2_Chromosome; do (
    		if test ! -f ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_${selchr}_raw.vcf; then
    			$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${selchr} ${input} -ploidy ${calcploidy} -O ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((calcploidy * maxHaplotype)) &&
    			gunzip ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz &&
    			wait
    		fi
    		) &
    		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
    			wait
    		fi
    	done
    	wait
    	cd ../snpcall
    	grep -h '^#' ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt
    	cat ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_*_raw.vcf | awk '!/^#/' > all.vcf
    	cat vcf_header.txt all.vcf > ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf
    	rm vcf_header.txt all.vcf *.vcf.gz.tbi ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_*_raw.vcf
    else
    	echo $Get_Chromosome | tr ',' '\n' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2}' > ${projdir}/refgenomes/panref.list
    	$GATK --java-options "$Xmx2 -XX:+UseParallelGC -XX:ParallelGCThreads=$threads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${projdir}/refgenomes/panref.list ${input}-ploidy $calcploidy -O ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf.gz --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((calcploidy * maxHaplotype)) &&
    	cd ../snpcall
    	gunzip ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf &&
    	wait
    fi

    cd ${projdir}/preprocess
    mv *_${ref1%.f*}_${ref2%.f*}_precall* ./processed/

    ######################
    echo -e "${magenta}- performing SNP calling across entire genome (subgenome 1 and 3) ${white}\n"

    j=-I; input=""; k=""
    for i in $(ls *_${ref1%.f*}_${ref3%.f*}_precall.bam); do
    	k="${j} ${i}"; input="${input} ${k}"
    done
    calcploidy=$(($ploidy_ref1 + $ploidy_ref3))
    Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' | awk -v pat1=${ref1%.f*} -v pat2=${ref3%.f*} 'pat1 || pat2')
    if [[ ! -z "$Exclude_Chromosome" ]]; then
    	for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
    		Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
    	done
    fi

    if [[ -z "$Get_Chromosome" ]]; then
    	for selchr in $Get2_Chromosome; do (
    		if test ! -f ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw.vcf; then
    			$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${selchr} ${input} -ploidy ${calcploidy} -O ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((calcploidy * maxHaplotype)) &&
    			gunzip ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz &&
    			wait
    		fi
    		) &
    		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
    			wait
    		fi
    	done
    	wait
    	cd ../snpcall
    	grep -h '^#' ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt
    	cat ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf | awk '!/^#/' > all.vcf
    	cat vcf_header.txt all.vcf > ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf
    	rm vcf_header.txt all.vcf *.vcf.gz.tbi ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf
    else
    	echo $Get_Chromosome | tr ',' '\n' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2}' > ${projdir}/refgenomes/panref.list
    	$GATK --java-options "$Xmx2 -XX:+UseParallelGC -XX:ParallelGCThreads=$threads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${projdir}/refgenomes/panref.list ${input}-ploidy $calcploidy -O ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf.gz --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((calcploidy * maxHaplotype)) &&
    	cd ../snpcall
    	gunzip ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf &&
    	wait
    fi

    cd ${projdir}/preprocess
    mv *_${ref1%.f*}_${ref3%.f*}_precall* ./processed/

    ######################
    echo -e "${magenta}- performing SNP calling across entire genome (subgenome 2 and 3) ${white}\n"

    j=-I; input=""; k=""
    for i in $(ls *_${ref1%.f*}_${ref4%.f*}_precall.bam); do
    	k="${j} ${i}"; input="${input} ${k}"
    done
    calcploidy=$(($ploidy_ref1 + $ploidy_ref4))
    Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' | awk -v pat1=${ref1%.f*} -v pat2=${ref4%.f*} 'pat1 || pat2')
    if [[ ! -z "$Exclude_Chromosome" ]]; then
    	for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
    		Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
    	done
    fi

    if [[ -z "$Get_Chromosome" ]]; then
    	for selchr in $Get2_Chromosome; do (
    		if test ! -f ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref4%.f*}_${calcploidy}x_${selchr}_raw.vcf; then
    			$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${selchr} ${input} -ploidy ${calcploidy} -O ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref4%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((calcploidy * maxHaplotype)) &&
    			gunzip ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref4%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz &&
    			wait
    		fi
    		) &
    		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
    			wait
    		fi
    	done
    	wait
    	cd ../snpcall
    	grep -h '^#' ${pop}_${ref1%.f*}_${ref4%.f*}_${calcploidy}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt
    	cat ${pop}_${ref1%.f*}_${ref4%.f*}_${calcploidy}x_*_raw.vcf | awk '!/^#/' > all.vcf
    	cat vcf_header.txt all.vcf > ${pop}_${ref1%.f*}_${ref4%.f*}_${calcploidy}x_raw.vcf
    	rm vcf_header.txt all.vcf *.vcf.gz.tbi ${pop}_${ref1%.f*}_${ref4%.f*}_${calcploidy}x_*_raw.vcf
    else
    	echo $Get_Chromosome | tr ',' '\n' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2}' > ${projdir}/refgenomes/panref.list
    	$GATK --java-options "$Xmx2 -XX:+UseParallelGC -XX:ParallelGCThreads=$threads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${projdir}/refgenomes/panref.list ${input}-ploidy $calcploidy -O ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref4%.f*}_${calcploidy}x_raw.vcf.gz --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((calcploidy * maxHaplotype)) &&
    	cd ../snpcall
    	gunzip ${pop}_${ref1%.f*}_${ref4%.f*}_${calcploidy}x_raw.vcf &&
    	wait
    fi

    cd ${projdir}/preprocess
    mv *_${ref1%.f*}_${ref4%.f*}_precall* ./processed/

	  ######################

	  echo -e "${magenta}- performing SNP calling on subgenome-1 ${white}\n"
	  j=-I; input=""; k=""
	  for i in $(ls *_${ref1%.f*}_precall.bam); do
	    k="${j} ${i}"; input="${input} ${k}"
	  done
	  Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}'| awk -v pat=${ref1%.f*} '$0 ~ pat' )
		if [[ ! -z "$Exclude_Chromosome" ]]; then
			for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
				Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
			done
		fi

	  if [[ -z "$Get_Chromosome" ]]; then
	    for selchr in $Get2_Chromosome; do (
				if test ! -f ${projdir}/snpcall/${pop}_${ref1%.f*}_${ploidy_ref1}x_${selchr}_raw.vcf; then
	      	$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${selchr} ${input} -ploidy $ploidy_ref1 -O ${projdir}/snpcall/${pop}_${ref1%.f*}_${ploidy_ref1}x_${selchr}_raw.vcf.gz --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((ploidy_ref1 * maxHaplotype)) &&
	      	gunzip ${projdir}/snpcall/${pop}_${ref1%.f*}_${ploidy_ref1}x_${selchr}_raw.vcf.gz &&
					wait
				fi
	      ) &
	      if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
	        wait
	      fi
	    done
			wait
	    cd ../snpcall
	    grep -h '^#' ${pop}_${ref1%.f*}_${ploidy_ref1}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt
	    cat ${pop}_${ref1%.f*}_${ploidy_ref1}x_*_raw.vcf | awk '!/^#/' > all.vcf
	    cat vcf_header.txt all.vcf > ${pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf
	    rm vcf_header.txt all.vcf *.vcf.gz.tbi ${pop}_${ref1%.f*}_${ploidy_ref1}x_*_raw.vcf
	  else
	    echo $Get_Chromosome | tr ',' '\n' | awk -v pat=${ref1%.f*} '$0 ~ pat' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2}' > ${projdir}/refgenomes/panref.list
	    $GATK --java-options "$Xmx2 -XX:+UseParallelGC -XX:ParallelGCThreads=$threads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${projdir}/refgenomes/panref.list ${input}-ploidy $ploidy -O ${projdir}/snpcall/${pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf.gz --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((ploidy_ref1 * maxHaplotype)) &&
	    cd ../snpcall
	    gunzip ${pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf &&
	    wait
	  fi

	  cd ${projdir}/preprocess
	  mv *_${ref1%.f*}_precall* ./processed/

	  ######################

	  echo -e "${magenta}- performing SNP calling on subgenome-2 ${white}\n"
	  j=-I; input=""; k=""
	  for i in $(ls *_${ref2%.f*}_precall.bam); do
	    k="${j} ${i}"; input="${input} ${k}"
	  done
	  Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}'| awk -v pat=${ref2%.f*} '$0 ~ pat' )
		if [[ ! -z "$Exclude_Chromosome" ]]; then
			for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
				Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
			done
		fi

	  if [[ -z "$Get_Chromosome" ]]; then
	    for selchr in $Get2_Chromosome; do (
				if test ! -f ${projdir}/snpcall/${pop}_${ref2%.f*}_${ploidy_ref2}x_${selchr}_raw.vcf; then
	      	$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${selchr} ${input} -ploidy $ploidy_ref2 -O ${projdir}/snpcall/${pop}_${ref2%.f*}_${ploidy_ref2}x_${selchr}_raw.vcf.gz --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((ploidy_ref2 * maxHaplotype)) &&
	      	gunzip ${projdir}/snpcall/${pop}_${ref2%.f*}_${ploidy_ref2}x_${selchr}_raw.vcf.gz &&
					wait
				fi
	      ) &
	      if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
	        wait
	      fi
	    done
			wait
	    cd ../snpcall
	    grep -h '^#' ${pop}_${ref2%.f*}_${ploidy_ref2}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt
	    cat ${pop}_${ref2%.f*}_${ploidy_ref2}x_*_raw.vcf | awk '!/^#/' > all.vcf
	    cat vcf_header.txt all.vcf > ${pop}_${ref2%.f*}_${ploidy_ref2}x_raw.vcf
	    rm vcf_header.txt all.vcf *.vcf.gz.tbi ${pop}_${ref2%.f*}_${ploidy_ref2}x_*_raw.vcf
	  else
	    echo $Get_Chromosome | tr ',' '\n' | awk -v pat=${ref2%.f*} '$0 ~ pat' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2}' > ${projdir}/refgenomes/panref.list
	    $GATK --java-options "$Xmx2 -XX:+UseParallelGC -XX:ParallelGCThreads=$threads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${projdir}/refgenomes/panref.list ${input}-ploidy $ploidy -O ${projdir}/snpcall/${pop}_${ref2%.f*}_${ploidy_ref2}x_raw.vcf.gz --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((ploidy_ref2 * maxHaplotype)) &&
	    cd ../snpcall
	    gunzip ${pop}_${ref2%.f*}_${ploidy_ref2}x_raw.vcf &&
	    wait
	  fi

	  cd ${projdir}/preprocess
	  mv *_${ref2%.f*}_precall* ./processed/

    ######################
	  echo -e "${magenta}- performing SNP calling on subgenome-3 ${white}\n"

	  j=-I; input=""; k=""
	  for i in $(ls *_${ref3%.f*}_precall.bam); do
	    k="${j} ${i}"; input="${input} ${k}"
	  done
	  Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}'| awk -v pat=${ref3%.f*} '$0 ~ pat' )
		if [[ ! -z "$Exclude_Chromosome" ]]; then
			for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
				Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
			done
		fi

	  if [[ -z "$Get_Chromosome" ]]; then
	    for selchr in $Get2_Chromosome; do (
				if test ! -f ${projdir}/snpcall/${pop}_${ref3%.f*}_${ploidy_ref3}x_${selchr}_raw.vcf; then
	      	$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${selchr} ${input} -ploidy $ploidy_ref3 -O ${projdir}/snpcall/${pop}_${ref3%.f*}_${ploidy_ref3}x_${selchr}_raw.vcf.gz --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((ploidy_ref3 * maxHaplotype)) &&
	      	gunzip ${projdir}/snpcall/${pop}_${ref3%.f*}_${ploidy_ref3}x_${selchr}_raw.vcf.gz &&
					wait
				fi
	      ) &
	      if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
	        wait
	      fi
	    done
			wait
	    cd ../snpcall
	    grep -h '^#' ${pop}_${ref3%.f*}_${ploidy_ref3}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt
	    cat ${pop}_${ref3%.f*}_${ploidy_ref3}x_*_raw.vcf | awk '!/^#/' > all.vcf
	    cat vcf_header.txt all.vcf > ${pop}_${ref3%.f*}_${ploidy_ref3}x_raw.vcf
	    rm vcf_header.txt all.vcf *.vcf.gz.tbi ${pop}_${ref3%.f*}_${ploidy_ref3}x_*_raw.vcf
	  else
	    echo $Get_Chromosome | tr ',' '\n' | awk -v pat=${ref3%.f*} '$0 ~ pat' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2}' > ${projdir}/refgenomes/panref.list
	    $GATK --java-options "$Xmx2 -XX:+UseParallelGC -XX:ParallelGCThreads=$threads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${projdir}/refgenomes/panref.list ${input}-ploidy $ploidy -O ${projdir}/snpcall/${pop}_${ref3%.f*}_${ploidy_ref3}x_raw.vcf.gz --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((ploidy_ref3 * maxHaplotype)) &&
	    cd ../snpcall
	    gunzip ${pop}_${ref3%.f*}_${ploidy_ref3}x_raw.vcf &&
	    wait
	  fi

	  cd ${projdir}/preprocess
	  mv *_${ref3%.f*}_precall* ./processed/

	fi


	if [[ "$joint_calling" == false ]]; then

    echo -e "${magenta}- performing SNP calling on subgenome-1, 2 and 3 ${white}\n"

    for i in $(cat ${projdir}/${samples_list} ); do (
    	if test ! -f "${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw.vcf" && test ! -f "${projdir}/GVCF_1_2_3_done.txt"; then
    		if test ! -f "${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.g.vcf.gz"; then
    			if [[ -z "$Get_Chromosome" ]]; then
    				$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -I ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_precall.bam -ploidy $ploidy -O ${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.hold.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) &&
    				wait
    			else
    				echo $Get_Chromosome | tr ',' '\n' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2}' > ${projdir}/refgenomes/panref.list
    				$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${projdir}/refgenomes/panref.list -I ${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_precall.bam -ploidy $ploidy -O ${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.hold.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) &&
    				wait
    			fi
    			mv "${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.hold.g.vcf.gz" "${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.g.vcf.gz" && \
    			mv "${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.hold.g.vcf.gz.tbi" "${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.g.vcf.gz.tbi" &&
    			wait
    		fi
    	fi ) &
    	if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
    		wait
    	fi
    done
    wait

    if [[ ! -f "${projdir}/call123_${samples_list}" ]]; then touch "${projdir}/call123_${samples_list}"; fi
    if [[ "$samples_list" == "samples_list_node_1.txt" ]] && test ! -f ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw.vcf*; then
    	align=$(ls ${projdir}/call123_samples_list_node_* | wc -l)
    	while [[ "$align" -lt $nodes ]]; do sleep 30; align=$(ls ${projdir}/call123_samples_list_node_* | wc -l); done
    	if [[ $align == $nodes ]]; then
    		rm ${projdir}/call123_${samples_list}
    		cd ${projdir}/snpcall
    		cz=$(ls *.g.vcf.gz | wc -l)
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
    			for i in $(ls *_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.g.vcf.gz); do
    				k="${j} ${i}"; input="${input} ${k}"
    			done
    			if [[ -z "$Get_Chromosome" ]]; then
    				Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' )
    			else
    				Get2_Chromosome=$(echo $Get_Chromosome | tr ',' '\n')
    			fi
    			if [[ ! -z "$Exclude_Chromosome" ]]; then
    				for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
    					Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
    				done
    			fi
    			for selchr in $Get2_Chromosome; do
    				$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  GenomicsDBImport ${input} -L ${selchr} --genomicsdb-workspace-path ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_${selchr}_raw
    			done
    			for selchr in $Get2_Chromosome; do (
    				if test ! -f "${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_${selchr}_raw.vcf.gz"; then
    					$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  GenotypeGVCFs -R ${projdir}/refgenomes/panref.fasta -L ${selchr} -V gendb://${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_${selchr}_raw -O ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_${selchr}_raw.hold.vcf.gz && \
    					rm -r ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_${selchr}_raw && \
    					mv "${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_${selchr}_raw.hold.vcf.gz" "${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_${selchr}_raw.vcf.gz"
    					mv "${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_${selchr}_raw.hold.vcf.gz.tbi" "${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_${selchr}_raw.vcf.gz.tbi" &&
    					wait
    				fi
    				if LC_ALL=C gzip -l ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_${selchr}_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
    					:
    				else
    					rm ../cohorts*/${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_*_raw.vcf.gz*
    					rm ../cohorts*/${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw.vcf*
    					rm ../${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw_cohorts*.vcf*
    					echo -e "${magenta}- \n- SNP calling failed probably due to insufficient memory ${white}\n"
    					echo -e "${magenta}- \n- Exiting pipeline in 5 seconds ${white}\n"
    					sleep 5 && exit 1
    				fi ) &
    				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
    					wait
    				fi
    			done
    			wait
    			for g in $(ls ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_*_raw.vcf.gz); do (
    				gunzip $g ) &
    				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
    					wait
    				fi
    			done
    			grep -h '^#' ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt
    			cat ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_*_raw.vcf | awk '!/^#/' > all.vcf
    			cat vcf_header.txt all.vcf > ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw.vcf
    			rm vcf_header.txt all.vcf
    			rm ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_*_raw.vcf ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_*_raw.vcf.gz.tbi
    			$bcftools view -I ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw.vcf -O z -o ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw.vcf.gz
    			$bcftools index ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw.vcf.gz

    			$bcftools annotate -x FORMAT/PL ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw.vcf.gz > ../${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw_${dir%/}.vcf
    			cd ../
    			$bcftools view -I ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw_${dir%/}.vcf -O z -o ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw_${dir%/}.vcf.gz
    			$bcftools index ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw_${dir%/}.vcf.gz
    			wait
    		done
    		wait
    		if [[ `ls -1 *cohorts*.vcf.gz 2>/dev/null | wc -l` -gt 1 ]]; then
    			$bcftools merge *cohorts*.vcf.gz --force-samples -m all > ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw.vcf
    		else
    			cp *cohorts*.vcf.gz ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw.vcf.gz
    			gunzip ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw.vcf.gz
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
    		mv *_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_precall* ./processed/
    	fi
    fi
    wait && touch ${projdir}/GVCF_1_2_3_done.txt

    ######################
    echo -e "${magenta}- performing SNP calling on subgenome-1 and 2 ${white}\n"

    calcploidy=$(($ploidy_ref1 + $ploidy_ref2))
    for i in $(cat ${projdir}/${samples_list} ); do (
    	if test ! -f "${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf" && test ! -f "${projdir}/GVCF_1_2_done.txt"; then
    		if test ! -f "${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref2%.f*}.g.vcf.gz"; then
    			if [[ -z "$Get_Chromosome" ]]; then
    				$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -I ${i%.f*}_${ref1%.f*}_${ref2%.f*}_precall.bam -ploidy $calcploidy -O ${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref2%.f*}.hold.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((calcploidy * maxHaplotype)) &&
    				wait
    			else
    				echo $Get_Chromosome | tr ',' '\n' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2}' > ${projdir}/refgenomes/panref.list
    				$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${projdir}/refgenomes/panref.list -I ${i%.f*}_${ref1%.f*}_${ref2%.f*}_precall.bam -ploidy $calcploidy -O ${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref2%.f*}.hold.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((calcploidy * maxHaplotype)) &&
    				wait
    			fi
    			mv "${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref2%.f*}.hold.g.vcf.gz" "${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref2%.f*}.g.vcf.gz" && \
    			mv "${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref2%.f*}.hold.g.vcf.gz.tbi" "${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref2%.f*}.g.vcf.gz.tbi" &&
    			wait
    		fi
    	fi ) &
    	if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
    		wait
    	fi
    done
    wait

    if [[ ! -f "${projdir}/call12_${samples_list}" ]]; then touch "${projdir}/call12_${samples_list}"; fi
    if [[ "$samples_list" == "samples_list_node_1.txt" ]] && test ! -f ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf*; then
    	align=$(ls ${projdir}/call12_samples_list_node_* | wc -l)
    	while [[ "$align" -lt $nodes ]]; do sleep 30; align=$(ls ${projdir}/call12_samples_list_node_* | wc -l); done
    	if [[ $align == $nodes ]]; then
    		rm ${projdir}/call12_${samples_list}
    		cd ${projdir}/snpcall
    		cz=$(ls *.g.vcf.gz | wc -l)
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
    			for i in $(ls *_${ref1%.f*}_${ref2%.f*}.g.vcf.gz); do
    				k="${j} ${i}"; input="${input} ${k}"
    			done
    			if [[ -z "$Get_Chromosome" ]]; then
    				Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' )
    			else
    				Get2_Chromosome=$(echo $Get_Chromosome | tr ',' '\n')
    			fi
    			if [[ ! -z "$Exclude_Chromosome" ]]; then
    				for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
    					Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
    				done
    			fi
    			for selchr in $Get2_Chromosome; do
    				$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  GenomicsDBImport ${input} -L ${selchr} --genomicsdb-workspace-path ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_${selchr}_raw
    			done
    			for selchr in $Get2_Chromosome; do (
    				if test ! -f "${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz"; then
    					$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  GenotypeGVCFs -R ${projdir}/refgenomes/panref.fasta -L ${selchr} -V gendb://${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_${selchr}_raw -O ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_${selchr}_raw.hold.vcf.gz && \
    					rm -r ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_${selchr}_raw && \
    					mv "${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_${selchr}_raw.hold.vcf.gz" "${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz"
    					mv "${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_${selchr}_raw.hold.vcf.gz.tbi" "${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz.tbi" &&
    					wait
    				fi
    				if LC_ALL=C gzip -l ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
    					:
    				else
    					rm ../cohorts*/${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_*_raw.vcf.gz*
    					rm ../cohorts*/${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf*
    					rm ../${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw_cohorts*.vcf*
    					echo -e "${magenta}- \n- SNP calling failed probably due to insufficient memory ${white}\n"
    					echo -e "${magenta}- \n- Exiting pipeline in 5 seconds ${white}\n"
    					sleep 5 && exit 1
    				fi ) &
    				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
    					wait
    				fi
    			done
    			wait
    			for g in $(ls ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_*_raw.vcf.gz); do (
    				$gunzip $g ) &
    				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
    					wait
    				fi
    			done
    			grep -h '^#' ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt
    			cat ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_*_raw.vcf | awk '!/^#/' > all.vcf
    			cat vcf_header.txt all.vcf > ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf
    			rm vcf_header.txt all.vcf
    			rm ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_*_raw.vcf ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_*_raw.vcf.gz.tbi
    			$bcftools view -I ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf -O z -o ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf.gz
    			$bcftools index ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf.gz

    			$bcftools annotate -x FORMAT/PL ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf.gz > ../${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw_${dir%/}.vcf
    			cd ../
    			$bcftools view -I ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw_${dir%/}.vcf -O z -o ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw_${dir%/}.vcf.gz
    			$bcftools index ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw_${dir%/}.vcf.gz
    			wait
    		done
    		wait
    		if [[ `ls -1 *cohorts*.vcf.gz 2>/dev/null | wc -l` -gt 1 ]]; then
    			$bcftools merge *cohorts*.vcf.gz --force-samples -m all > ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf
    		else
    			cp *cohorts*.vcf.gz ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf.gz
    			$gunzip ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf.gz
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
    		mv *_${ref1%.f*}_${ref2%.f*}_precall* ./processed/
    	fi
    fi
    wait && touch ${projdir}/GVCF_1_2_done.txt

    ######################
    echo -e "${magenta}- performing SNP calling on subgenome-1 and 3 ${white}\n"

    calcploidy=$(($ploidy_ref1 + $ploidy_ref3))
    for i in $(cat ${projdir}/${samples_list} ); do (
    	if test ! -f "${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf" && test ! -f "${projdir}/GVCF_1_3_done.txt"; then
    		if test ! -f "${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref3%.f*}.g.vcf.gz"; then
    			if [[ -z "$Get_Chromosome" ]]; then
    				$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -I ${i%.f*}_${ref1%.f*}_${ref3%.f*}_precall.bam -ploidy $calcploidy -O ${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref3%.f*}.hold.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((calcploidy * maxHaplotype)) &&
    				wait
    			else
    				echo $Get_Chromosome | tr ',' '\n' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2}' > ${projdir}/refgenomes/panref.list
    				$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R "${projdir}/refgenomes/panref.fasta" -L ${projdir}/refgenomes/panref.list -I ${i%.f*}_${ref1%.f*}_${ref3%.f*}_precall.bam -ploidy $calcploidy -O ${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref3%.f*}.hold.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((calcploidy * maxHaplotype)) &&
    				wait
    			fi
    			mv "${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref3%.f*}.hold.g.vcf.gz" "${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref3%.f*}.g.vcf.gz" && \
    			mv "${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref3%.f*}.hold.g.vcf.gz.tbi" "${projdir}/snpcall/${i%.f*}_${ref1%.f*}_${ref3%.f*}.g.vcf.gz.tbi" &&
    			wait
    		fi
    	fi ) &
    	if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
    		wait
    	fi
    done
    wait

    if [[ ! -f "${projdir}/call13_${samples_list}" ]]; then touch "${projdir}/call13_${samples_list}"; fi
    if [[ "$samples_list" == "samples_list_node_1.txt" ]] && test ! -f ${projdir}/snpcall/${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf*; then
    	align=$(ls ${projdir}/call13_samples_list_node_* | wc -l)
    	while [[ "$align" -lt $nodes ]]; do sleep 30; align=$(ls ${projdir}/call13_samples_list_node_* | wc -l); done
    	if [[ $align == $nodes ]]; then
    		rm ${projdir}/call13_${samples_list}
    		cd ${projdir}/snpcall
    		cz=$(ls *.g.vcf.gz | wc -l)
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
    			for i in $(ls *_${ref1%.f*}_${ref3%.f*}.g.vcf.gz); do
    				k="${j} ${i}"; input="${input} ${k}"
    			done
    			if [[ -z "$Get_Chromosome" ]]; then
    				Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' )
    			else
    				Get2_Chromosome=$(echo $Get_Chromosome | tr ',' '\n')
    			fi
    			if [[ ! -z "$Exclude_Chromosome" ]]; then
    				for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
    					Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
    				done
    			fi
    			for selchr in $Get2_Chromosome; do
    				$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  GenomicsDBImport ${input} -L ${selchr} --genomicsdb-workspace-path ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw
    			done
    			for selchr in $Get2_Chromosome; do (
    				if test ! -f "${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz"; then
    					$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  GenotypeGVCFs -R ${projdir}/refgenomes/panref.fasta -L ${selchr} -V gendb://${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw -O ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw.hold.vcf.gz && \
    					rm -r ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw && \
    					mv "${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw.hold.vcf.gz" "${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz"
    					mv "${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw.hold.vcf.gz.tbi" "${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz.tbi" &&
    					wait
    				fi
    				if LC_ALL=C gzip -l ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
    					:
    				else
    					rm ../cohorts*/${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf.gz*
    					rm ../cohorts*/${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf*
    					rm ../${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw_cohorts*.vcf*
    					echo -e "${magenta}- \n- SNP calling failed probably due to insufficient memory ${white}\n"
    					echo -e "${magenta}- \n- Exiting pipeline in 5 seconds ${white}\n"
    					sleep 5 && exit 1
    				fi ) &
    				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
    					wait
    				fi
    			done
    			wait
    			for g in $(ls ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf.gz); do (
    				$gunzip $g ) &
    				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
    					wait
    				fi
    			done
    			grep -h '^#' ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt
    			cat ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf | awk '!/^#/' > all.vcf
    			cat vcf_header.txt all.vcf > ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf
    			rm vcf_header.txt all.vcf
    			rm ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf.gz.tbi
    			$bcftools view -I ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf -O z -o ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf.gz
    			$bcftools index ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf.gz

    			$bcftools annotate -x FORMAT/PL ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf.gz > ../${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw_${dir%/}.vcf
    			cd ../
    			$bcftools view -I ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw_${dir%/}.vcf -O z -o ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw_${dir%/}.vcf.gz
    			$bcftools index ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw_${dir%/}.vcf.gz
    			wait
    		done
    		wait
    		if [[ `ls -1 *cohorts*.vcf.gz 2>/dev/null | wc -l` -gt 1 ]]; then
    			$bcftools merge *cohorts*.vcf.gz --force-samples -m all > ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf
    		else
    			cp *cohorts*.vcf.gz ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf.gz
    			$gunzip ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf.gz
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
    		mv *_${ref1%.f*}_${ref3%.f*}_precall* ./processed/
    	fi
    fi
    wait && touch ${projdir}/GVCF_1_3_done.txt

    ######################
    echo -e "${magenta}- performing SNP calling on subgenome-2 and 3 ${white}\n"

    calcploidy=$(($ploidy_ref2 + $ploidy_ref3))
    for i in $(cat ${projdir}/${samples_list} ); do (
    	if test ! -f "${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf" && test ! -f "${projdir}/GVCF_2_3_done.txt"; then
    		if test ! -f "${projdir}/snpcall/${i%.f*}_${ref2%.f*}_${ref3%.f*}.g.vcf.gz"; then
    			if [[ -z "$Get_Chromosome" ]]; then
    				$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -I ${i%.f*}_${ref2%.f*}_${ref3%.f*}_precall.bam -ploidy $calcploidy -O ${projdir}/snpcall/${i%.f*}_${ref2%.f*}_${ref3%.f*}.hold.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((calcploidy * maxHaplotype)) &&
    				wait
    			else
    				echo $Get_Chromosome | tr ',' '\n' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2}' > ${projdir}/refgenomes/panref.list
    				$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R "${projdir}/refgenomes/panref.fasta" -L ${projdir}/refgenomes/panref.list -I ${i%.f*}_${ref2%.f*}_${ref3%.f*}_precall.bam -ploidy $calcploidy -O ${projdir}/snpcall/${i%.f*}_${ref2%.f*}_${ref3%.f*}.hold.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((calcploidy * maxHaplotype)) &&
    				wait
    			fi
    			mv "${projdir}/snpcall/${i%.f*}_${ref2%.f*}_${ref3%.f*}.hold.g.vcf.gz" "${projdir}/snpcall/${i%.f*}_${ref2%.f*}_${ref3%.f*}.g.vcf.gz" && \
    			mv "${projdir}/snpcall/${i%.f*}_${ref2%.f*}_${ref3%.f*}.hold.g.vcf.gz.tbi" "${projdir}/snpcall/${i%.f*}_${ref2%.f*}_${ref3%.f*}.g.vcf.gz.tbi" &&
    			wait
    		fi
    	fi ) &
    	if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
    		wait
    	fi
    done
    wait

    if [[ ! -f "${projdir}/call23_${samples_list}" ]]; then touch "${projdir}/call23_${samples_list}"; fi
    if [[ "$samples_list" == "samples_list_node_1.txt" ]] && test ! -f ${projdir}/snpcall/${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf*; then
    	align=$(ls ${projdir}/call23_samples_list_node_* | wc -l)
    	while [[ "$align" -lt $nodes ]]; do sleep 30; align=$(ls ${projdir}/call23_samples_list_node_* | wc -l); done
    	if [[ $align == $nodes ]]; then
    		rm ${projdir}/call23_${samples_list}
    		cd ${projdir}/snpcall
    		cz=$(ls *.g.vcf.gz | wc -l)
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
    			for i in $(ls *_${ref2%.f*}_${ref3%.f*}.g.vcf.gz); do
    				k="${j} ${i}"; input="${input} ${k}"
    			done
    			if [[ -z "$Get_Chromosome" ]]; then
    				Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' )
    			else
    				Get2_Chromosome=$(echo $Get_Chromosome | tr ',' '\n')
    			fi
    			if [[ ! -z "$Exclude_Chromosome" ]]; then
    				for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
    					Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
    				done
    			fi
    			for selchr in $Get2_Chromosome; do
    				$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  GenomicsDBImport ${input} -L ${selchr} --genomicsdb-workspace-path ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw
    			done
    			for selchr in $Get2_Chromosome; do (
    				if test ! -f "${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz"; then
    					$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  GenotypeGVCFs -R ${projdir}/refgenomes/panref.fasta -L ${selchr} -V gendb://${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw -O ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw.hold.vcf.gz && \
    					rm -r ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw && \
    					mv "${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw.hold.vcf.gz" "${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz"
    					mv "${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw.hold.vcf.gz.tbi" "${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz.tbi" &&
    					wait
    				fi
    				if LC_ALL=C gzip -l ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_${selchr}_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
    					:
    				else
    					rm ../cohorts*/${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf.gz*
    					rm ../cohorts*/${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf*
    					rm ../${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw_cohorts*.vcf*
    					echo -e "${magenta}- \n- SNP calling failed probably due to insufficient memory ${white}\n"
    					echo -e "${magenta}- \n- Exiting pipeline in 5 seconds ${white}\n"
    					sleep 5 && exit 1
    				fi ) &
    				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
    					wait
    				fi
    			done
    			wait
    			for g in $(ls ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf.gz); do (
    				$gunzip $g ) &
    				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
    					wait
    				fi
    			done
    			grep -h '^#' ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt
    			cat ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf | awk '!/^#/' > all.vcf
    			cat vcf_header.txt all.vcf > ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf
    			rm vcf_header.txt all.vcf
    			rm ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf.gz.tbi
    			$bcftools view -I ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf -O z -o ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf.gz
    			$bcftools index ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf.gz

    			$bcftools annotate -x FORMAT/PL ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf.gz > ../${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw_${dir%/}.vcf
    			cd ../
    			$bcftools view -I ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw_${dir%/}.vcf -O z -o ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw_${dir%/}.vcf.gz
    			$bcftools index ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw_${dir%/}.vcf.gz
    			wait
    		done
    		wait
    		if [[ `ls -1 *cohorts*.vcf.gz 2>/dev/null | wc -l` -gt 1 ]]; then
    			$bcftools merge *cohorts*.vcf.gz --force-samples -m all > ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf
    		else
    			cp *cohorts*.vcf.gz ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf.gz
    			$gunzip ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf.gz
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
    		mv *_${ref2%.f*}_${ref3%.f*}_precall* ./processed/
    	fi
    fi
    wait && touch ${projdir}/GVCF_2_3_done.txt

    ######################

			echo -e "${magenta}- performing SNP calling on subgenome-1 ${white}\n"
			for i in $(cat ${projdir}/${samples_list} ); do (
				if test ! -f "${pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf" && test ! -f "${projdir}/GVCF_1_done.txt"; then
					if test ! -f ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf.gz; then
						if [[ -z "$Get_Chromosome" ]]; then
							$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -I ${i%.f*}_${ref1%.f*}_precall.bam -ploidy $ploidy_ref1 -O ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((ploidy_ref1 * maxHaplotype)) &&
							wait
						else
						  echo $Get_Chromosome | tr ',' '\n' | awk -v pat=${ref1%.f*} '$0 ~ pat' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2}' > ${projdir}/refgenomes/panref.list
							$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${projdir}/refgenomes/panref.list -I ${i%.f*}_${ref1%.f*}_precall.bam -ploidy $ploidy_ref1 -O ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((ploidy_ref1 * maxHaplotype)) &&
							wait
						fi
						mv "${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz" "${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf.gz" && \
						mv "${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz.tbi" "${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf.gz.tbi" &&
						wait
					fi
				fi ) &
				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
					wait
				fi
			done
			wait


			if [[ ! -f "${projdir}/call1_${samples_list}" ]]; then touch "${projdir}/call1_${samples_list}"; fi
			if [[ "$samples_list" == "samples_list_node_1.txt" ]] && test ! -f $${projdir}/snpcall/{pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf*; then
			  align=$(ls ${projdir}/call1_samples_list_node_* | wc -l)
			  while [[ "$align" -lt $nodes ]]; do sleep 30; align=$(ls ${projdir}/call1_samples_list_node_* | wc -l); done
			  if [[ $align == $nodes ]]; then
			    rm ${projdir}/call1_${samples_list}
					cd ${projdir}/snpcall
					cz=$(ls *.g.vcf.gz | wc -l)
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
						for i in $(ls *_${ref1%.f*}.g.vcf.gz); do
							k="${j} ${i}"; input="${input} ${k}"
						done
						if [[ -z "$Get_Chromosome" ]]; then
							Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}'| awk -v pat=${ref1%.f*} '$0 ~ pat' )
						else
							Get2_Chromosome=$(echo $Get_Chromosome | tr ',' '\n')
						fi
						if [[ ! -z "$Exclude_Chromosome" ]]; then
							for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
								Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
							done
						fi
						for selchr in $Get2_Chromosome; do
							$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  GenomicsDBImport ${input} -L ${selchr} --genomicsdb-workspace-path ${pop}_${ref1%.f*}_${ploidy_ref1}x_${selchr}_raw
						done
						for selchr in $Get2_Chromosome; do (
							if test ! -f "${pop}_${ref1%.f*}_${ploidy_ref1}x_${selchr}_raw.vcf.gz"; then
								$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  GenomicsDBImport ${input} -L ${selchr} --genomicsdb-workspace-path ${pop}_${ref1%.f*}_${ploidy_ref1}x_${selchr}_raw
								$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  GenotypeGVCFs -R ${projdir}/refgenomes/panref.fasta -L ${selchr} -V gendb://${pop}_${ref1%.f*}_${ploidy_ref1}x_${selchr}_raw -O ${pop}_${ref1%.f*}_${ploidy_ref1}x_${selchr}_raw.hold.vcf.gz && \
								rm -r ${pop}_${ref1%.f*}_${ploidy_ref1}x_${selchr}_raw && \
								mv "${pop}_${ref1%.f*}_${ploidy_ref1}x_${selchr}_raw.hold.vcf.gz" "${pop}_${ref1%.f*}_${ploidy_ref1}x_${selchr}_raw.vcf.gz"
								mv "${pop}_${ref1%.f*}_${ploidy_ref1}x_${selchr}_raw.hold.vcf.gz.tbi" "${pop}_${ref1%.f*}_${ploidy_ref1}x_${selchr}_raw.vcf.gz.tbi" &&
								wait
							fi
							if LC_ALL=C gzip -l ${pop}_${ref1%.f*}_${ploidy_ref1}x_${selchr}_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
								:
							else
								rm ../cohorts*/${pop}_${ref1%.f*}_${ploidy_ref1}x_*_raw.vcf.gz*
								rm ../cohorts*/${pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf*
								rm ../${pop}_${ref1%.f*}_${ploidy_ref1}x_raw_cohorts*.vcf
								echo -e "${magenta}- \n- SNP calling failed probably due to insufficient memory ${white}\n"
								echo -e "${magenta}- \n- Exiting pipeline in 5 seconds ${white}\n"
								sleep 5 && exit 1
							fi ) &
							if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
								wait
							fi
						done
						wait
						for g in $(ls ${pop}_${ref1%.f*}_${ploidy_ref1}x_*_raw.vcf.gz); do (
							gunzip $g ) &
							if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
								wait
							fi
						done
						grep -h '^#' ${pop}_${ref1%.f*}_${ploidy_ref1}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt
						cat ${pop}_${ref1%.f*}_${ploidy_ref1}x_*_raw.vcf | awk '!/^#/' > all.vcf
						cat vcf_header.txt all.vcf > ${pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf
						rm vcf_header.txt all.vcf
						rm ${pop}_${ref1%.f*}_${ploidy_ref1}x_*_raw.vcf ${pop}_${ref1%.f*}_${ploidy_ref1}x_*_raw.vcf.gz.tbi
						$bcftools view -I ${pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf -O z -o ${pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf.gz
						$bcftools index ${pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf.gz

						$bcftools annotate -x FORMAT/PL ${pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf.gz > ../${pop}_${ref1%.f*}_${ploidy_ref1}x_raw_${dir%/}.vcf
						cd ../
						$bcftools view -I ${pop}_${ref1%.f*}_${ploidy_ref1}x_raw_${dir%/}.vcf -O z -o ${pop}_${ref1%.f*}_${ploidy_ref1}x_raw_${dir%/}.vcf.gz
						$bcftools index ${pop}_${ref1%.f*}_${ploidy_ref1}x_raw_${dir%/}.vcf.gz
						wait
					done
					wait
					if [[ `ls -1 *cohorts*.vcf.gz 2>/dev/null | wc -l` -gt 1 ]]; then
						$bcftools erge *cohorts*.vcf.gz --force-samples -m all > ${pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf
					else
						cp *cohorts*.vcf.gz ${pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf.gz
						gunzip ${pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf.gz
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
					mv *_${ref1%.f*}_precall* ./processed/
				fi
			fi
			wait && touch ${projdir}/GVCF_1_done.txt

		######################

			echo -e "${magenta}- performing SNP calling on subgenome-2 ${white}\n"
			for i in $(cat ${projdir}/${samples_list} ); do (
				if test ! -f "${pop}_${ref2%.f*}_${ploidy_ref2}x_raw.vcf" && test ! -f "${projdir}/GVCF_2_done.txt"; then
					if test ! -f "${projdir}/snpcall/${i%.f*}_${ref2%.f*}.g.vcf.gz"; then
						if [[ -z "$Get_Chromosome" ]]; then
							$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -I ${i%.f*}_${ref2%.f*}_precall.bam -ploidy $ploidy_ref2 -O ${projdir}/snpcall/${i%.f*}_${ref2%.f*}.hold.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((ploidy_ref2 * maxHaplotype)) &&
							wait
						else
						  echo $Get_Chromosome | tr ',' '\n' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2}' > ${projdir}/refgenomes/panref.list
							$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${projdir}/refgenomes/panref.list -I ${i%.f*}_${ref2%.f*}_precall.bam -ploidy $ploidy_ref2 -O ${projdir}/snpcall/${i%.f*}_${ref2%.f*}.hold.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((ploidy_ref2 * maxHaplotype)) &&
							wait
						fi
						mv "${projdir}/snpcall/${i%.f*}_${ref2%.f*}.hold.g.vcf.gz" "${projdir}/snpcall/${i%.f*}_${ref2%.f*}.g.vcf.gz" && \
						mv "${projdir}/snpcall/${i%.f*}_${ref2%.f*}.hold.g.vcf.gz.tbi" "${projdir}/snpcall/${i%.f*}_${ref2%.f*}.g.vcf.gz.tbi" &&
						wait
					fi
				fi ) &
				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
					wait
				fi
			done
			wait


			if [[ ! -f "${projdir}/call2_${samples_list}" ]]; then touch "${projdir}/call2_${samples_list}"; fi
			if [[ "$samples_list" == "samples_list_node_1.txt" ]] && test ! -f ${projdir}/snpcall/${pop}_${ref2%.f*}_${ploidy_ref2}x_raw.vcf*; then
			  align=$(ls ${projdir}/call2_samples_list_node_* | wc -l)
			  while [[ "$align" -lt $nodes ]]; do sleep 30; align=$(ls ${projdir}/call2_samples_list_node_* | wc -l); done
			  if [[ $align == $nodes ]]; then
			    rm ${projdir}/call2_${samples_list}
					cd ${projdir}/snpcall
					cz=$(ls *.g.vcf.gz | wc -l)
					i=0
					for f in `find . -maxdepth 1 -iname '*.g.vcf.gz' -type f | shuf`; do
						d=cohorts_$(printf %02d $((i/cz+1)))
						mkdir -p $d
						mv "$f" $d; mv "${f}.tbi" $d
						let i++
					done


					for dir in *cohorts*/; do
						cd $dir
						j=--variant; input=""; k=""
						for i in $(ls *_${ref2%.f*}.g.vcf.gz); do
							k="${j} ${i}"; input="${input} ${k}"
						done
						if [[ -z "$Get_Chromosome" ]]; then
							Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' | awk -v pat=${ref2%.f*} '$0 ~ pat')
						else
							Get2_Chromosome=$(echo $Get_Chromosome | tr ',' '\n')
						fi
						if [[ ! -z "$Exclude_Chromosome" ]]; then
							for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
								Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
							done
						fi
						for selchr in $Get2_Chromosome; do
							$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  GenomicsDBImport ${input} -L ${selchr} --genomicsdb-workspace-path ${pop}_${ref2%.f*}_${ploidy_ref2}x_${selchr}_raw
						done
						for selchr in $Get2_Chromosome; do (
							if test ! -f "${pop}_${ref2%.f*}_${ploidy_ref2}x_${selchr}_raw.vcf.gz"; then
								$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  GenotypeGVCFs -R ${projdir}/refgenomes/panref.fasta -L ${selchr} -V gendb://${pop}_${ref2%.f*}_${ploidy_ref2}x_${selchr}_raw -O ${pop}_${ref2%.f*}_${ploidy_ref2}x_${selchr}_raw.hold.vcf.gz && \
								rm -r ${pop}_${ref2%.f*}_${ploidy_ref2}x_${selchr}_raw && \
								mv "${pop}_${ref2%.f*}_${ploidy_ref2}x_${selchr}_raw.hold.vcf.gz" "${pop}_${ref2%.f*}_${ploidy_ref2}x_${selchr}_raw.vcf.gz" &&
								wait
							fi
							if LC_ALL=C gzip -l ${pop}_${ref2%.f*}_${ploidy_ref2}x_${selchr}_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
								:
							else
								rm ../cohorts*/${pop}_${ref2%.f*}_${ploidy_ref2}x_*_raw.vcf.gz*
								rm ../cohorts*/${pop}_${ref2%.f*}_${ploidy_ref2}x_raw.vcf*
								rm ../${pop}_${ref2%.f*}_${ploidy_ref2}x_raw_cohorts*.vcf*
								echo -e "${magenta}- \n- SNP calling failed probably due to insufficient memory ${white}\n"
								echo -e "${magenta}- \n- Exiting pipeline in 5 seconds ${white}\n"
								sleep 5 && exit 1
							fi ) &
							if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
								wait
							fi
						done
						wait
						for g in $(ls ${pop}_${ref2%.f*}_${ploidy_ref2}x_*_raw.vcf.gz); do (
							gunzip $g ) &
							if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
								wait
							fi
						done
						grep -h '^#' ${pop}_${ref2%.f*}_${ploidy_ref2}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt
						cat ${pop}_${ref2%.f*}_${ploidy_ref2}x_*_raw.vcf | awk '!/^#/' > all.vcf
						cat vcf_header.txt all.vcf > ${pop}_${ref2%.f*}_${ploidy_ref2}x_raw.vcf
						rm vcf_header.txt all.vcf
						rm ${pop}_${ref2%.f*}_${ploidy_ref2}x_*_raw.vcf ${pop}_${ref2%.f*}_${ploidy_ref2}x_*_raw.vcf.gz.tbi
						$bcftools view -I ${pop}_${ref2%.f*}_${ploidy_ref2}x_raw.vcf -O z -o ${pop}_${ref2%.f*}_${ploidy_ref2}x_raw.vcf.gz
						$bcftools index ${pop}_${ref2%.f*}_${ploidy_ref2}x_raw.vcf.gz

						$bcftools annotate -x FORMAT/PL ${pop}_${ref2%.f*}_${ploidy_ref2}x_raw.vcf.gz > ../${pop}_${ref2%.f*}_${ploidy_ref2}x_raw_${dir%/}.vcf
						cd ../
						$bcftools view -I ${pop}_${ref2%.f*}_${ploidy_ref2}x_raw_${dir%/}.vcf -O z -o ${pop}_${ref2%.f*}_${ploidy_ref2}x_raw_${dir%/}.vcf.gz
						$bcftools index ${pop}_${ref2%.f*}_${ploidy_ref2}x_raw_${dir%/}.vcf.gz
						wait
					done
					wait
					if [[ `ls -1 *cohorts*.vcf.gz 2>/dev/null | wc -l` -gt 1 ]]; then
						$bcftools merge *cohorts*.vcf.gz --force-samples -m all > ${pop}_${ref2%.f*}_${ploidy_ref2}x_raw.vcf
					else
						cp *cohorts*.vcf.gz ${pop}_${ref2%.f*}_${ploidy_ref2}x_raw.vcf.gz
						gunzip ${pop}_${ref2%.f*}_${ploidy_ref2}x_raw.vcf.gz
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
					mv *_${ref2%.f*}_precall* ./processed/
				fi
			fi
		wait && touch ${projdir}/GVCF_2_done.txt

    ######################

			echo -e "${magenta}- performing SNP calling on subgenome-3 ${white}\n"
			for i in $(cat ${projdir}/${samples_list} ); do (
				if test ! -f "${pop}_${ref3%.f*}_${ploidy_ref3}x_raw.vcf" && test ! -f "${projdir}/GVCF_3_done.txt"; then
					if test ! -f "${projdir}/snpcall/${i%.f*}_${ref3%.f*}.g.vcf.gz"; then
						if [[ -z "$Get_Chromosome" ]]; then
							$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -I ${i%.f*}_${ref3%.f*}_precall.bam -ploidy $ploidy_ref3 -O ${projdir}/snpcall/${i%.f*}_${ref3%.f*}.hold.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population $((ploidy_ref3 * maxHaplotype)) &&
							wait
						else
						  echo $Get_Chromosome | tr ',' '\n' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2}' > ${projdir}/refgenomes/panref.list
							$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ${projdir}/refgenomes/panref.fasta -L ${projdir}/refgenomes/panref.list -I ${i%.f*}_${ref3%.f*}_precall.bam -ploidy $ploidy_ref3 -O ${projdir}/snpcall/${i%.f*}_${ref3%.f*}.hold.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases $softclip --max-reads-per-alignment-start $downsample --minimum-mapping-quality 0 --max-num-haplotypes-in-population "$((ploidy_ref3 * maxHaplotype))" &&
							wait
						fi
						mv "${projdir}/snpcall/${i%.f*}_${ref3%.f*}.hold.g.vcf.gz" "${projdir}/snpcall/${i%.f*}_${ref3%.f*}.g.vcf.gz" && \
						mv "${projdir}/snpcall/${i%.f*}_${ref3%.f*}.hold.g.vcf.gz.tbi" "${projdir}/snpcall/${i%.f*}_${ref3%.f*}.g.vcf.gz.tbi" &&
						wait
					fi
				fi ) &
				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
					wait
				fi
			done
			wait


			if [[ ! -f "${projdir}/call3_${samples_list}" ]]; then touch "${projdir}/call3_${samples_list}"; fi
			if [[ "$samples_list" == "samples_list_node_1.txt" ]] && test ! -f ${projdir}/snpcall/${pop}_${ref3%.f*}_${ploidy_ref3}x_raw.vcf*; then
			  align=$(ls ${projdir}/call3_samples_list_node_* | wc -l)
			  while [[ "$align" -lt $nodes ]]; do sleep 30; align=$(ls ${projdir}/call2_samples_list_node_* | wc -l); done
			  if [[ $align == $nodes ]]; then
			    rm ${projdir}/call3_${samples_list}
					cd ${projdir}/snpcall
					cz=$(ls *.g.vcf.gz | wc -l)
					i=0
					for f in `find . -maxdepth 1 -iname '*.g.vcf.gz' -type f | shuf`; do
						d=cohorts_$(printf %02d $((i/cz+1)))
						mkdir -p $d
						mv "$f" $d; mv "${f}.tbi" $d
						let i++
					done


					for dir in *cohorts*/; do
						cd $dir
						j=--variant; input=""; k=""
						for i in $(ls *_${ref3%.f*}.g.vcf.gz); do
							k="${j} ${i}"; input="${input} ${k}"
						done
						if [[ -z "$Get_Chromosome" ]]; then
							Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' | awk -v pat=${ref3%.f*} '$0 ~ pat')
						else
							Get2_Chromosome=$(echo $Get_Chromosome | tr ',' '\n')
						fi
						if [[ ! -z "$Exclude_Chromosome" ]]; then
							for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
								Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
							done
						fi
						for selchr in $Get2_Chromosome; do
							$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  GenomicsDBImport ${input} -L ${selchr} --genomicsdb-workspace-path ${pop}_${ref3%.f*}_${ploidy_ref3}x_${selchr}_raw
						done
						for selchr in $Get2_Chromosome; do (
							if test ! -f "${pop}_${ref3%.f*}_${ploidy_ref3}x_${selchr}_raw.vcf.gz"; then
								$GATK --java-options "$Xmxg -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads"  GenotypeGVCFs -R ${projdir}/refgenomes/panref.fasta -L ${selchr} -V gendb://${pop}_${ref3%.f*}_${ploidy_ref3}x_${selchr}_raw -O ${pop}_${ref3%.f*}_${ploidy_ref3}x_${selchr}_raw.hold.vcf.gz && \
								rm -r ${pop}_${ref3%.f*}_${ploidy_ref3}x_${selchr}_raw && \
								mv "${pop}_${ref3%.f*}_${ploidy_ref3}x_${selchr}_raw.hold.vcf.gz" "${pop}_${ref3%.f*}_${ploidy_ref3}x_${selchr}_raw.vcf.gz" &&
								wait
							fi
							if LC_ALL=C gzip -l ${pop}_${ref3%.f*}_${ploidy_ref3}x_${selchr}_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
								:
							else
								rm ../cohorts*/${pop}_${ref3%.f*}_${ploidy_ref3}x_*_raw.vcf.gz*
								rm ../cohorts*/${pop}_${ref3%.f*}_${ploidy_ref3}x_raw.vcf*
								rm ../${pop}_${ref3%.f*}_${ploidy_ref3}x_raw_cohorts*.vcf*
								echo -e "${magenta}- \n- SNP calling failed probably due to insufficient memory ${white}\n"
								echo -e "${magenta}- \n- Exiting pipeline in 5 seconds ${white}\n"
								sleep 5 && exit 1
							fi ) &
							if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
								wait
							fi
						done
						wait
						for g in $(ls ${pop}_${ref3%.f*}_${ploidy_ref3}x_*_raw.vcf.gz); do (
							gunzip $g ) &
							if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
								wait
							fi
						done
						grep -h '^#' ${pop}_${ref3%.f*}_${ploidy_ref3}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt
						cat ${pop}_${ref3%.f*}_${ploidy_ref3}x_*_raw.vcf | awk '!/^#/' > all.vcf
						cat vcf_header.txt all.vcf > ${pop}_${ref3%.f*}_${ploidy_ref3}x_raw.vcf
						rm vcf_header.txt all.vcf
						rm ${pop}_${ref3%.f*}_${ploidy_ref3}x_*_raw.vcf ${pop}_${ref3%.f*}_${ploidy_ref3}x_*_raw.vcf.gz.tbi
						$bcftools view -I ${pop}_${ref3%.f*}_${ploidy_ref3}x_raw.vcf -O z -o ${pop}_${ref3%.f*}_${ploidy_ref3}x_raw.vcf.gz
						$bcftools index ${pop}_${ref3%.f*}_${ploidy_ref3}x_raw.vcf.gz

						$bcftools annotate -x FORMAT/PL ${pop}_${ref3%.f*}_${ploidy_ref3}x_raw.vcf.gz > ../${pop}_${ref3%.f*}_${ploidy_ref3}x_raw_${dir%/}.vcf
						cd ../
						$bcftools view -I ${pop}_${ref3%.f*}_${ploidy_ref3}x_raw_${dir%/}.vcf -O z -o ${pop}_${ref3%.f*}_${ploidy_ref3}x_raw_${dir%/}.vcf.gz
						$bcftools index ${pop}_${ref3%.f*}_${ploidy_ref3}x_raw_${dir%/}.vcf.gz
						wait
					done
					wait
					if [[ `ls -1 *cohorts*.vcf.gz 2>/dev/null | wc -l` -gt 1 ]]; then
						$bcftools merge *cohorts*.vcf.gz --force-samples -m all > ${pop}_${ref3%.f*}_${ploidy_ref3}x_raw.vcf
					else
						cp *cohorts*.vcf.gz ${pop}_${ref3%.f*}_${ploidy_ref3}x_raw.vcf.gz
						gunzip ${pop}_${ref3%.f*}_${ploidy_ref3}x_raw.vcf.gz
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
					mv *_${ref3%.f*}_precall* ./processed/
				fi
			fi
		wait && touch ${projdir}/GVCF_3_done.txt

		######################

		if [[ ! -f "${projdir}/callfinal_${samples_list}" ]]; then touch "${projdir}/callfinal_${samples_list}"; fi
		if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
		  align=$(ls ${projdir}/callfinal_samples_list_node_* | wc -l)
		  while [[ "$align" -lt $nodes ]]; do sleep 30; align=$(ls ${projdir}/callfinal_samples_list_node_* | wc -l); done
		  if [[ $align == $nodes ]]; then
		    rm ${projdir}/callfinal_${samples_list}
				echo -e "${magenta}- keeping *realign.bam & *realign.bam files in ./preprocess/processed/ ${white}\n"
				mv ${projdir}/preprocess/processed/* ${projdir}/preprocess/
				rmdir ${projdir}/preprocess/processed
			fi
		fi

	fi

}
cd $projdir
while [[ ! -f "$projdir/alignment_summaries/refgenome_paralogs.txt" ]]; do sleep 30; done
sleep 10
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
echo -e "${blue}\n############################################################################## ${yellow}\n- GBSapp is Performing SNP Filtering\n${blue}##############################################################################${white}\n"
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

mkdir 2x
mkdir 4x
mkdir 6x
mkdir 8x
cd ${projdir}/snpcall
if [ -z "$(ls -A *_DP_GT.txt 2>/dev/null)" ]; then
	if [ -z "$(ls -A *_x.vcf 2>/dev/null)" ]; then
		if [ -z "$(ls -A *_raw.vcf 2>/dev/null)" ]; then
			for g in *_raw.vcf.gz; do gunzip $g;	done
		fi
	fi
fi
wait

file2xG=$( if [ "$(ls -A *_2x_DP_GT.txt 2>/dev/null)" ]; then ls *_2x_DP_GT.txt | wc -l;  else echo 0; fi )
file2xV=$( if [ "$(ls -A *_2x_raw.vcf 2>/dev/null)" ]; then ls *_2x_raw.vcf | head -n 1 | wc -l; else echo 0; fi )
if [[ "${file2xG}" -lt 1 ]]; then
	if [[ "${file2xV}" -gt 0 ]]; then
		for i in *_2x_raw.vcf; do
			$GATK LeftAlignAndTrimVariants -R ${projdir}/refgenomes/panref.fasta -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
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
			$GATK LeftAlignAndTrimVariants -R ${projdir}/refgenomes/panref.fasta -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
			wait
			large_numerous_chrom &>> ${projdir}/log.out
			wait
			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf
			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%100000==2{x=file"4x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
			wait $PID
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
			$GATK LeftAlignAndTrimVariants -R ${projdir}/refgenomes/panref.fasta -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
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
			$GATK LeftAlignAndTrimVariants -R ${projdir}/refgenomes/panref.fasta -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
			wait
			large_numerous_chrom &>> ${projdir}/log.out
			wait
			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf
			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%100000==2{x=file"8x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
			wait $PID
			rm ${i%.vcf}0.vcf* ${i%.vcf}trim.vcf*
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
	for v in *.vcf; do gzip $v; done
fi
wait



cd $projdir
cd samples
window1=$(ls -S | head -1 | xargs zcat -fq | awk '{ print length }' | sort -n | tail -1)
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


for smiss in ${sample_missingness[@]}; do
for gmiss in ${genotype_missingness[@]}; do
if [[ -z "$p1" ]]; then
	if [ -d "${projdir}/snpfilter/2x" ]; then
		cd ${projdir}/snpfilter
		cp -r 2x 2x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./2x_diversity_gmiss"${gmiss}"_smiss"${smiss}"

		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_2x.R $pop $gmiss $smiss $minRD_2x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome_number
		wait
		rm ${pop}_2x_rawRD${minRD_2x}_DP_GT.txt ${pop}_2x_DP_GT.txt ${pop}_2x_rd${minRD_2x}.txt
		mkdir visualizations; mv *.tiff ./visualizations/
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
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_4x.R $pop $gmiss $smiss $minRD_4x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome_number
		wait
		rm ${pop}_4x_rawRD${minRD_4x}_DP_GT.txt ${pop}_4x_DP_GT.txt ${pop}_4x_rd${minRD_4x}.txt
		mkdir visualizations; mv *.tiff ./visualizations/
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
		haplome=$((haplome_number * 2))
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_6x.R $pop $gmiss $smiss $minRD_6x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome
		wait
		rm ${pop}_6x_rawRD${minRD_6x}_DP_GT.txt ${pop}_6x_DP_GT.txt ${pop}_6x_rd${minRD_6x}.txt
		mkdir visualizations; mv *.tiff ./visualizations/
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
		haplome=$((haplome_number * 4))
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_8x.R $pop $gmiss $smiss $minRD_8x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome
		wait
		rm ${pop}_8x_rawRD${minRD_8x}_DP_GT.txt ${pop}_8x_DP_GT.txt ${pop}_8x_rd${minRD_8x}.txt
		mkdir visualizations; mv *.tiff ./visualizations/
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
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_2x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_2x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$pseg" "$haplome_number"
		wait
		rm ${pop}_2x_rawRD${minRD_2x}_DP_GT.txt ${pop}_2x_DP_GT.txt ${pop}_2x_rd${minRD_2x}.txt
		mkdir visualizations; mv *.tiff ./visualizations/
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
		haplome=$((haplome_number * 2))
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_4x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_4x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$pseg" "$haplome"
		wait
		rm ${pop}_4x_rawRD${minRD_4x}_DP_GT.txt ${pop}_4x_DP_GT.txt ${pop}_4x_rd${minRD_4x}.txt
		mkdir visualizations; mv *.tiff ./visualizations/
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
		haplome=$((haplome_number * 3))
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_6x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_6x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$pseg" "$haplome"
		wait
		rm ${pop}_6x_rawRD${minRD_6x}_DP_GT.txt ${pop}_6x_DP_GT.txt ${pop}_6x_rd${minRD_6x}.txt
		mkdir visualizations; mv *.tiff ./visualizations/
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
		haplome=$((haplome_number * 4))
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_8x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_8x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$pseg" "$haplome"
		wait
		rm ${pop}_8x_rawRD${minRD_8x}_DP_GT.txt ${pop}_8x_DP_GT.txt ${pop}_8x_rd${minRD_8x}.txt
		mkdir visualizations; mv *.tiff ./visualizations/
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
done
done


cd ${projdir}/snpfilter/
find . -type f -empty -delete
find . -type d -empty -delete
wc -l *gmiss*/*dose.txt | awk '{print $2"\t"$1-1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > gmiss_smiss_titration.txt
wc -l *gmiss*/eliminated* | awk '{print $2"\t"$1-1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > eliminated_samples.txt
echo -e "gmiss_smiss_thresholds\t#_of_SNPs\t#_of_eliminated_samples\n----------------------\t---------\t-----------------------" > summary_precall.txt
awk 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,"\t",a[$1],"\t",$2}' gmiss_smiss_titration.txt eliminated_samples.txt | cat summary_precall.txt - > gmiss_smiss.txt
rm gmiss_smiss_titration.txt eliminated_samples.txt summary_precall.txt

wc -l *gmiss*/unique_mapped/*dose*.txt | awk '{print $2"\t"$1-1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > gmiss_smiss_titration.txt
wc -l *gmiss*/eliminated* | awk '{print $2"\t"$1-1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > eliminated_samples.txt
echo -e "gmiss_smiss_thresholds\t#_of_SNPs\t#_of_eliminated_samples\n----------------------\t---------\t-----------------------" > summary_precall.txt
awk 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,"\t",a[$1],"\t",$2}' gmiss_smiss_titration.txt eliminated_samples.txt | cat summary_precall.txt - > gmiss_smiss_unique_mapped.txt
rm gmiss_smiss_titration.txt eliminated_samples.txt summary_precall.txt

ls ./*/*maf*.txt | grep -v 'maf0.txt' | grep -v 'dose' | xargs rm
ls ./*/*_plusSD.txt | xargs rm 2> /dev/null
ls ./*/*SD_1_G*G*.txt | xargs rm 2> /dev/null
ls ./*/*binary*.txt | xargs rm 2> /dev/null

cd "$projdir"/snpfilter
for snpfilter_dir in $(ls -d */); do
	if [ -d "$snpfilter_dir" ]; then
		cd $snpfilter_dir
		for v in *dose.txt; do
			vcfdose=${v%_rd*}; vcfdose=${vcfdose#*_}
			zcat ../../snpcall/*${vcfdose}.vcf.gz | grep '^#' > ${v%.txt}.vcf
			awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$2,$3]){print b}}' <(gzip -dc ../../snpcall/*${vcfdose}.vcf.gz) $v >> ${v%.txt}.vcf
			gzip ${v%.txt}.vcf
		done
		wait
		for i in *dose*; do
			Rscript "${GBSapp_dir}"/scripts/R/hapmap_format.R "$i" "${i%.txt}.hmp.txt"
		done
		wait

		cd unique_mapped
		for v in *dose*; do
			vcfdose=${v%_rd*}; vcfdose=${vcfdose#*_}
			zcat ../../../snpcall/*${vcfdose}.vcf.gz | grep '^#' > ${v%.txt}.vcf
			awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$2,$3]){print b}}' <(gzip -dc ../../../snpcall/*${vcfdose}.vcf.gz) $v >> ${v%.txt}.vcf
			gzip ${v%.txt}.vcf
		done
		wait
		for i in *dose*; do
			Rscript "${GBSapp_dir}"/scripts/R/hapmap_format.R "$i" "${i%.txt}.hmp.txt"
		done
		wait

		cd ../../
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
	awk 'BEGIN{OFS="\t"}{print $2,$3}' *dose.txt > CHROM_POS.txt && \
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
echo -e "${magenta}- Under Development ${white}\n"
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
		time main &>> log.out
	fi
fi
if [ "$walkaway" == true ]; then
	if [ "$paralogfiltering_haplotyping" == 1 ]; then
		echo -e "${magenta}- performing sequence-based haplotyping and filtering ${white}\n"
		time main &>> log.out
	else
		echo -e "${magenta}- skipping sequence-based haplotyping and filtering ${white}\n"
	fi
fi


#####################################################################################################################################################
cd ${projdir}
rm compress_done.txt hapfilter_done.txt alignment_done.txt precall_done.txt GVCF*_done.txt samples_list_node_*.txt >/dev/null 2>&1
touch Analysis_Complete
wait
echo -e "${magenta}- Run Complete. ${white}\n"
