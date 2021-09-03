if R --version; then
	:
else
	module add R
fi

if [ -z "$threads" ]; then
	threads=$(nproc --all)
	if [[ "$threads" -ge 4 ]]; then
		threads=$((threads-2))
	fi
fi
echo threads: "$threads"
if [ -z "$haplome_number" ]; then
	haplome_number=1
else
	haplome_number=$((haplome_number))
fi
echo haplome_number: "$haplome_number"
if [ -z "$multilocus" ]; then
	paralogs=10
fi
if [ "$multilocus" == true ]; then
	paralogs=10
fi
if [ "$multilocus" == false ]; then
	paralogs=1
fi
if [ -z "$maxHaplotype" ]; then
	maxHaplotype=6
fi

if [ -z "$p2" ]; then
	p2=$p1
fi
if [ -z "$paleopolyploid" ]; then
	paleopolyploid=false
fi

if [[ "$biased_downsample" -gt 0 ]]; then
	if [[ "$ubbiased_downsample" -gt 0 ]]; then
		echo -e "${magenta}- warning: GBSapp will perform unbiased downsampling, and then biased downsampling ${white}\n"
		sleep 5 && exit 1
	fi
fi
if [ -z "$biased_downsample" ]; then
	biased_downsample=0
	if [ -z "$unbiased_downsample" ]; then
		unbiased_downsample=25
		echo -e "${magenta}- GBSapp will perform unbiased downsampling ${white}\n"
	fi
else
	if [[ "$biased_downsample" -eq 0 ]]; then
		if [[ "$unbiased_downsample" -eq 0 ]]; then
			echo -e "${magenta}- downsampling will not be performed, i.e. using all reads ${white}\n"
		fi
		if [ -z "$unbiased_downsample" ]; then
			unbiased_downsample=25
			echo -e "${magenta}- GBSapp will perform unbiased downsampling ${white}\n"
		fi
	fi
	if [[ "$biased_downsample" -gt 0 ]]; then
		if [[ "$unbiased_downsample" -eq 0 ]]; then
			echo -e "${magenta}- biased downsampling be performed ${white}\n"
		fi
		if [ -z "$unbiased_downsample" ]; then
			unbiased_downsample=0
			echo -e "${magenta}- unbiased downsampling will not be performed ${white}\n"
		fi
		if [ "$unbiased_downsample" -gt 0 ]; then
						echo -e "${magenta}- warning: GBSapp will perform unbiased downsampling, and then biased downsampling ${white}\n"
		fi
	fi
fi
if [ -z "$softclip" ]; then
	softclip=true
fi



cd $projdir
echo -e "${blue}\n############################################################################## ${yellow}\n- Index Reference Genome \n${blue}##############################################################################${white}\n"
main () {
cd $projdir
cd refgenomes
for i in *.gz; do
	gunzip $i >/dev/null 2>&1
done

file=panref.dict
if [[ -f "$file" ]]; then
	:
else
	checknfastafiles=$(ls *.f* | wc -l)
	if [[ $checknfastafiles -gt 3 ]]; then
		echo -e "${magenta}- expecting only 3 fasta files for reference genome ${white}\n"
		echo -e "${magenta}- GBSapp will quit in 5 seconds ${white}\n"
		sleep 5 && exit 1
	fi
	reffilenames=($ref1 $ref2 $ref3)
	for refg in "${reffilenames[@]}"; do
		export ncontigscaffold=$(grep '>' $refg | wc -l)
		if [[ $ncontigscaffold -gt 200 ]]; then
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


	mkdir split
	for reffile in $(ls *.f*); do
		awk -v reffile=${reffile%.f*} '/^>/{close("file.fasta"f);f++}{print $0 > "./split/file_"reffile"_"f}' "$reffile"
	done
	cd split
	checksplit=$( wc -c file_* | head -n -1 | awk '($1 > 500000000 )' | wc -l )
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
		for i in *.f*; do (
			mv $i ./old_"${i%.f*}_fasta.txt"
			cat ./split/${i%.f*}_Chr* > $i ) &
		done
		wait
		rm -r split
	else
		cd ${projdir}/refgenomes
		rm -r split
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

file=panref.dict
if test -f "$file"; then
	echo -e "${magenta}- indexed genome already exist ${white}\n"
else
	echo -e "${magenta}- indexing reference genome ${white}\n"
	$bwa index -a bwtsw panref.fasta
	$samtools faidx panref.fasta
	$java $Xmxg -XX:ParallelGCThreads=$gthreads -jar $picard CreateSequenceDictionary REFERENCE= panref.fasta OUTPUT=panref.dict
fi

declare -a arr=("panref*.dict" "panref*.amb" "panref*.bwt" "panref*.pac" "panref.fasta" "panref*.ann" "panref*.fai" "panref*.sa")
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
time main 2> log.out


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
			time main &>> ../log.out
		fi
	else
		echo -e "${magenta}- samples' PE fastq filenames requires formatting (i.e. needs to end in "_R1.fastq" and "_R2.fastq") ${white}\n"
		sleep 5 && exit 1
	fi
else
	time main &>> ../log.out
fi



main () {
	cd $projdir
	mkdir preprocess
	mkdir snpcall
	mkdir alignment_summaries
	mkdir ./alignment_summaries/copy_number
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
echo -e "${blue}\n############################################################################## ${yellow}\n- GBSapp is Performing Read Alignments & Alignment Post-Processing\n${blue}##############################################################################${white}\n"
main () {
cd $projdir

cd samples
for i in $(ls -S *.f* | grep -v R2.f); do (
  if gzip -t $i; then
	gunzip -c $i | awk 'NR%2==0' | awk 'NR%2' > ${i%.f*}_uniq.txt &
	if test -f ${i%.f*}_R2*; then gunzip -c ${i%.f*}_R2* | awk 'NR%2==0' | awk 'NR%2' > ${i%.f*}_R2_uniq.txt &
	fi
	if test -f ${i%.f*}.R2*; then gunzip -c ${i%.f*}.R2* | awk 'NR%2==0' | awk 'NR%2' > ${i%.f*}_R2_uniq.txt &
	fi
  else
	awk 'NR%2==0' $i | awk 'NR%2' > ${i%.f*}_uniq.txt &
	if test -f ${i%.f*}_R2*; then awk 'NR%2==0' ${i%.f*}_R2* | awk 'NR%2' > ${i%.f*}_R2_uniq.txt &
	fi
	if test -f ${i%.f*}.R2*; then awk 'NR%2==0' ${i%.f*}.R2* | awk 'NR%2' > ${i%.f*}_R2_uniq.txt &
	fi
  fi
  wait
  if test -f "${i%.f*}_R2_uniq.txt"; then
	:
  else
	touch "${i%.f*}_R2_uniq.txt"
  fi
  export LC_ALL=C; paste -d ~ ${i%.f*}_uniq.txt ${i%.f*}_R2_uniq.txt | expand -t $(( $(wc -L < $i ) + 2 )) | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{$1=$1};1' | \
  awk '{gsub(" /"," "); print}' | awk '{gsub("/\n","\n"); print}' | awk '{gsub("/"," "); print}' | awk '{gsub(" ","\t"); print}' > ${i%.f*}_rdrefseq.txt
  awk 'NF==2 {print ">seq"NR"_se-"$1"\n"$2}' ${i%.f*}_rdrefseq.txt > ${i%.f*}_rdrefseq_se.txt
  awk 'NF==3 {print ">seq"NR"_pe-"$0}' ${i%.f*}_rdrefseq.txt | awk '{print $1"\n"$3}' > ${i%.f*}_uniq_R2.fasta
  awk 'NF==3 {print ">seq"NR"_pe-"$0}' ${i%.f*}_rdrefseq.txt | awk '{print $1"\n"$2}' | cat - ${i%.f*}_rdrefseq_se.txt > ${i%.f*}_uniq_R1.fasta
  nempty=$( wc -l ${i%.f*}_uniq_R2.fasta | awk '{print $1}' )
  if [[ "$nempty" -gt 0 ]]; then
	$bwa mem -t $prepthreads ../refgenomes/panref.fasta ${i%.f*}_uniq_R1.fasta ${i%.f*}_uniq_R2.fasta > ../preprocess/${i%.f*}_del.sam
  else
	$bwa mem -t $prepthreads ../refgenomes/panref.fasta ${i%.f*}_uniq_R1.fasta > ../preprocess/${i%.f*}_del.sam
  fi
  rm ${i%.f*}_*_uniq* ${i%.f*}_uniq* ${i%.f*}_rdrefseq*
  printf '\n###---'${i%.f*}'---###\n' > ../alignment_summaries/${i%.f*}_summ.txt && \
  $samtools flagstat ../preprocess/${i%.f*}_del.sam >> ../alignment_summaries/${i%.f*}_summ.txt && \
  printf '########################################################################################################\n\n' >> ../alignment_summaries/${i%.f*}_summ.txt && \
  $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard SortSam I=../preprocess/${i%.f*}_del.sam O=../preprocess/${i%.f*}_del.bam SORT_ORDER=coordinate && \
  $samtools view -h ../preprocess/${i%.f*}_del.bam | awk '{if ($1=="@HD" || $1=="@SQ" || $1=="@PG" || $7=="*" || $7=="=") { print }}' > ../preprocess/${i%.f*}_del.sam  && \
  awk -F ${ref1%.f*}_ '{print NF-1, NR}' ../preprocess/${i%.f*}_del.sam > ../preprocess/${i%.f*}_${ref1%.f*}_count.txt  && \
  awk -F ${ref2%.f*}_ '{print NF-1, NR}' ../preprocess/${i%.f*}_del.sam > ../preprocess/${i%.f*}_${ref2%.f*}_count.txt  && \
  awk -F ${ref3%.f*}_ '{print NF-1, NR}' ../preprocess/${i%.f*}_del.sam > ../preprocess/${i%.f*}_${ref3%.f*}_count.txt  && \
  paste <(awk '{print $1}' ../preprocess/${i%.f*}_${ref1%.f*}_count.txt ) <(awk '{print $1}' ../preprocess/${i%.f*}_${ref2%.f*}_count.txt ) <(awk '{print $1}' ../preprocess/${i%.f*}_${ref3%.f*}_count.txt ) > ../preprocess/${i%.f*}_pancount.txt && \
  paste -d "\t" ../preprocess/${i%.f*}_pancount.txt ../preprocess/${i%.f*}_del.sam > ../preprocess/${i%.f*}_2del.sam  && \

  printf 'copy\tFrequency\tPercentage\n' > ../alignment_summaries/copy_number/${i%.f*}_copy_number_Unique_Read_histogram.txt && \
  grep -vwE "(@HD|@SQ|@PG)" ../preprocess/${i%.f*}_2del.sam | awk '{print $1"\t"$2"\t"$3}' | awk 'BEGIN {FS=OFS="\t"} {sum=0; n=0; for(i=1;i<=NF;i++) {sum+=$i; ++n} print sum}' |\
  sort | uniq -c | awk '{print $2"\t"$1}' > ../alignment_summaries/copy_number/${i%.f*}_copy_number.txt  && \
  awk 'NR==FNR{sum+= $2; next;} {printf("%s\t%s\t%3.3f%%\t%3.0f\n",$1,$2,100*$2/sum,100*$2/sum)}' ../alignment_summaries/copy_number/${i%.f*}_copy_number.txt ../alignment_summaries/copy_number/${i%.f*}_copy_number.txt > ../alignment_summaries/copy_number/${i%.f*}_plot.txt && \
  unset IFS; printf "%s\t%s\t%s\t%*s\n" $(sed 's/$/ |/' ../alignment_summaries/copy_number/${i%.f*}_plot.txt) | tr ' ' '|' | sort -k1 -n >> ../alignment_summaries/copy_number/${i%.f*}_copy_number_Unique_Read_histogram.txt && rm ../alignment_summaries/copy_number/${i%.f*}_copy_number.txt ../alignment_summaries/copy_number/${i%.f*}_plot.txt && \

  awk '{if ($4=="@HD" || $4=="@SQ" || $4=="@PG"); print $0}' ../preprocess/${i%.f*}_2del.sam > ../preprocess/${i%.f*}_3del.sam  && \
  grep -vwE "(@HD|@SQ|@PG)" ../preprocess/${i%.f*}_2del.sam | awk -v paralogs=$paralogs '{ if ($1 >= 1 && $1 <= paralogs && $2 == 0 && $3 == 0); print $0}' >> ../preprocess/${i%.f*}_3del.sam   && \
  awk '!($1=$2=$3="")' ../preprocess/${i%.f*}_3del.sam | awk '{$1=$1};1' > ../preprocess/${i%.f*}_${ref1%.f*}.sam  && \
	awk '/@HD/ || /@SQ/{print}' ../preprocess/${i%.f*}_${ref1%.f*}.sam | awk '{gsub(/ /,"\t"); print}' > ../preprocess/${i%.f*}_heading.sam
  awk '!/^@/ { print }' ../preprocess/${i%.f*}_${ref1%.f*}.sam | awk '{gsub(/^.*-/,"",$1); print $0}' | awk '{gsub(/^.*_/,"",$2); print $0}' > ../preprocess/${i%.f*}_uniq.sam

	if [ "$multilocus" == "true" ]; then
		echo "nloci~mapQ~CHROM~POS" > ../alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}.txt
		grep -vwE "(@HD|@SQ|@PG)" ../preprocess/${i%.f*}_2del.sam | awk '{print $1"~"$8"~"$6"~"$7}' >> ../alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}.txt
		awk '!visited[$0]++' ../alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}.txt > ../alignment_summaries/temp_${i%.f*}.txt
		mv ../alignment_summaries/temp_${i%.f*}.txt ../alignment_summaries/refgenome_paralogs_${i%.f*}_${ref1%.f*}.txt
	fi

	unds=$(( unbiased_downsample * ploidy_ref1 ))
	awk -F' ' 'BEGIN{FS=OFS=" "} {if (a[$3,$4]<$1) {a[$3,$4]=$1; data[$3,$4]=$0}} END{for (i in a) print data[i]}' ../preprocess/${i%.f*}_uniq.sam | \
	awk '{print $1,$3,$4}' | grep -v '*' > ../preprocess/${i%.f*}_chrpos.txt
	awk '{print $1,$0}' ../preprocess/${i%.f*}_uniq.sam |
	awk -F" " 'BEGIN{FS=OFS=" "} {if (NR==FNR) {a[$2FS$3]=$1; next} if ($4FS$5 in a) {$1=a[$4FS$5]} print}' ../preprocess/${i%.f*}_chrpos.txt - | \
	awk '!( $4 == "*" )' | awk '{print $1,$0}' | awk -v unds=$unds 'BEGIN{FS=OFS=" "} { if($2 > unds) {printf("%.3f", $1=$1/unds)} else {printf($1=$3);} $1="" }1' | \
	awk -v unds=$unds 'BEGIN{FS=OFS=" "} { if($2 > unds) {printf("%.f", $1=$3/$1)} else {printf($1=$3);} $1="" }1' | \
	awk '!( $1 < "1" )' |  awk '!($2=$3="")' > ../preprocess/${i%.f*}_uniqsamp.sam

  for j in $(LC_ALL=C; sort -n -k1,1 ../preprocess/${i%.f*}_uniqsamp.sam | grep -v @ | awk '{print $1}' | uniq); do
  	awk -v n="^${j}" '$0~n{print $0}' ../preprocess/${i%.f*}_uniqsamp.sam | awk -v n="$j" '{for(i=0;i<n;i++) print}' >> ../preprocess/${i%.f*}_exp.sam
  done; wait
  awk '{print "seq"NR"_"$0}' ../preprocess/${i%.f*}_exp.sam | awk '{gsub(/ /,"\t"); print}' | awk '{$11 = $10; print}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$11); print}' |\
  awk -v paralogs=$paralogs 'BEGIN { OFS = "\t"; }; { if ($5 == 0 || paralogs > 1) $5 = 60; else $5 = $5; }; 1' | awk '!( $6 == "*" )' | cat ../preprocess/${i%.f*}_heading.sam - > ../preprocess/${i%.f*}_${ref1%.f*}.sam
	rm ../preprocess/${i%.f*}_exp.sam ../preprocess/${i%.f*}_chrpos.txt ../preprocess/${i%.f*}_uniqsamp.sam ../preprocess/${i%.f*}_uniq.sam ../preprocess/${i%.f*}_heading.sam

  awk '{if ($4=="@HD" || $4=="@SQ" || $4=="@PG"); print $0}' ../preprocess/${i%.f*}_2del.sam > ../preprocess/${i%.f*}_3del.sam  && \
  grep -vwE "(@HD|@SQ|@PG)" ../preprocess/${i%.f*}_2del.sam | awk -v paralogs=$paralogs '{ if ($1 == 0 && $2 >= 1 && $2 <= paralogs && $3 == 0); print $0}' >> ../preprocess/${i%.f*}_3del.sam   && \
  awk '!($1=$2=$3="")' ../preprocess/${i%.f*}_3del.sam | awk '{$1=$1};1' > ../preprocess/${i%.f*}_${ref2%.f*}.sam  && \
	awk '/@HD/ || /@SQ/{print}' ../preprocess/${i%.f*}_${ref2%.f*}.sam | awk '{gsub(/ /,"\t"); print}' > ../preprocess/${i%.f*}_heading.sam
  awk '!/^@/ { print }' ../preprocess/${i%.f*}_${ref2%.f*}.sam | awk '{gsub(/^.*-/,"",$1); print $0}' | awk '{gsub(/^.*_/,"",$2); print $0}' > ../preprocess/${i%.f*}_uniq.sam

	if [ "$multilocus" == "true" ]; then
		echo "nloci~mapQ~CHROM~POS" > ../alignment_summaries/refgenome_paralogs_${i%.f*}_${ref2%.f*}.txt
		grep -vwE "(@HD|@SQ|@PG)" ../preprocess/${i%.f*}_2del.sam | awk '{print $2"~"$8"~"$6"~"$7}' >> ../alignment_summaries/refgenome_paralogs_${i%.f*}_${ref2%.f*}.txt
		awk '!visited[$0]++' ../alignment_summaries/refgenome_paralogs_${i%.f*}_${ref2%.f*}.txt > ../alignment_summaries/temp_${i%.f*}.txt
		mv ../alignment_summaries/temp_${i%.f*}.txt ../alignment_summaries/refgenome_paralogs_${i%.f*}_${ref2%.f*}.txt
	fi

	unds=$(( unbiased_downsample * ploidy_ref2 ))
	awk -F' ' 'BEGIN{FS=OFS=" "} {if (a[$3,$4]<$1) {a[$3,$4]=$1; data[$3,$4]=$0}} END{for (i in a) print data[i]}' ../preprocess/${i%.f*}_uniq.sam | \
	awk '{print $1,$3,$4}' | grep -v '*' > ../preprocess/${i%.f*}_chrpos.txt
	awk '{print $1,$0}' ../preprocess/${i%.f*}_uniq.sam |
	awk -F" " 'BEGIN{FS=OFS=" "} {if (NR==FNR) {a[$2FS$3]=$1; next} if ($4FS$5 in a) {$1=a[$4FS$5]} print}' ../preprocess/${i%.f*}_chrpos.txt - | \
	awk '!( $4 == "*" )' | awk '{print $1,$0}' | awk -v unds=$unds 'BEGIN{FS=OFS=" "} { if($2 > unds) {printf("%.3f", $1=$1/unds)} else {printf($1=$3);} $1="" }1' | \
	awk -v unds=$unds 'BEGIN{FS=OFS=" "} { if($2 > unds) {printf("%.f", $1=$3/$1)} else {printf($1=$3);} $1="" }1' | \
	awk '!( $1 < "1" )' |  awk '!($2=$3="")' > ../preprocess/${i%.f*}_uniqsamp.sam

  for j in $(LC_ALL=C; sort -n -k1,1 ../preprocess/${i%.f*}_uniqsamp.sam | grep -v @ | awk '{print $1}' | uniq); do
  	awk -v n="^${j}" '$0~n{print $0}' ../preprocess/${i%.f*}_uniqsamp.sam | awk -v n="$j" '{for(i=0;i<n;i++) print}' >> ../preprocess/${i%.f*}_exp.sam
  done; wait
  awk '{print "seq"NR"_"$0}' ../preprocess/${i%.f*}_exp.sam | awk '{gsub(/ /,"\t"); print}' | awk '{$11 = $10; print}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$11); print}' |\
  awk -v paralogs=$paralogs 'BEGIN { OFS = "\t"; }; { if ($5 == 0 || paralogs > 1) $5 = 60; else $5 = $5; }; 1' | awk '!( $6 == "*" )' | cat ../preprocess/${i%.f*}_heading.sam - > ../preprocess/${i%.f*}_${ref2%.f*}.sam
	rm ../preprocess/${i%.f*}_exp.sam ../preprocess/${i%.f*}_chrpos.txt ../preprocess/${i%.f*}_uniqsamp.sam ../preprocess/${i%.f*}_uniq.sam ../preprocess/${i%.f*}_heading.sam

  awk '{if ($4=="@HD" || $4=="@SQ" || $4=="@PG"); print $0}' ../preprocess/${i%.f*}_2del.sam > ../preprocess/${i%.f*}_3del.sam  && \
  grep -vwE "(@HD|@SQ|@PG)" ../preprocess/${i%.f*}_2del.sam | awk -v paralogs=$paralogs '{ if ($1 == 0 && $2 == 0 && $3 >= 1 && $3 <= paralogs); print $0}' >> ../preprocess/${i%.f*}_3del.sam   && \
  awk '!($1=$2=$3="")' ../preprocess/${i%.f*}_3del.sam | awk '{$1=$1};1' > ../preprocess/${i%.f*}_${ref3%.f*}.sam  && \
	awk '/@HD/ || /@SQ/{print}' ../preprocess/${i%.f*}_${ref3%.f*}.sam | awk '{gsub(/ /,"\t"); print}' > ../preprocess/${i%.f*}_heading.sam
  awk '!/^@/ { print }' ../preprocess/${i%.f*}_${ref3%.f*}.sam | awk '{gsub(/^.*-/,"",$1); print $0}' | awk '{gsub(/^.*_/,"",$2); print $0}' > ../preprocess/${i%.f*}_uniq.sam

	if [ "$multilocus" == "true" ]; then
		echo "nloci~mapQ~CHROM~POS" > ../alignment_summaries/refgenome_paralogs_${i%.f*}_${ref3%.f*}.txt
		grep -vwE "(@HD|@SQ|@PG)" ../preprocess/${i%.f*}_2del.sam | awk '{print $3"~"$8"~"$6"~"$7}' >> ../alignment_summaries/refgenome_paralogs_${i%.f*}_${ref3%.f*}.txt
		awk '!visited[$0]++' ../alignment_summaries/refgenome_paralogs_${i%.f*}_${ref3%.f*}.txt > ../alignment_summaries/temp_${i%.f*}.txt
		mv ../alignment_summaries/temp_${i%.f*}.txt ../alignment_summaries/refgenome_paralogs_${i%.f*}_${ref3%.f*}.txt
	fi

	unds=$(( unbiased_downsample * ploidy_ref3 ))
	awk -F' ' 'BEGIN{FS=OFS=" "} {if (a[$3,$4]<$1) {a[$3,$4]=$1; data[$3,$4]=$0}} END{for (i in a) print data[i]}' ../preprocess/${i%.f*}_uniq.sam | \
	awk '{print $1,$3,$4}' | grep -v '*' > ../preprocess/${i%.f*}_chrpos.txt
	awk '{print $1,$0}' ../preprocess/${i%.f*}_uniq.sam |
	awk -F" " 'BEGIN{FS=OFS=" "} {if (NR==FNR) {a[$2FS$3]=$1; next} if ($4FS$5 in a) {$1=a[$4FS$5]} print}' ../preprocess/${i%.f*}_chrpos.txt - | \
	awk '!( $4 == "*" )' | awk '{print $1,$0}' | awk -v unds=$unds 'BEGIN{FS=OFS=" "} { if($2 > unds) {printf("%.3f", $1=$1/unds)} else {printf($1=$3);} $1="" }1' | \
	awk -v unds=$unds 'BEGIN{FS=OFS=" "} { if($2 > unds) {printf("%.f", $1=$3/$1)} else {printf($1=$3);} $1="" }1' | \
	awk '!( $1 < "1" )' |  awk '!($2=$3="")' > ../preprocess/${i%.f*}_uniqsamp.sam

  for j in $(LC_ALL=C; sort -n -k1,1 ../preprocess/${i%.f*}_uniqsamp.sam | grep -v @ | awk '{print $1}' | uniq); do
  	awk -v n="^${j}" '$0~n{print $0}' ../preprocess/${i%.f*}_uniqsamp.sam| awk -v n="$j" '{for(i=0;i<n;i++) print}' >> ../preprocess/${i%.f*}_exp.sam
  done; wait
  awk '{print "seq"NR"_"$0}' ../preprocess/${i%.f*}_exp.sam | awk '{gsub(/ /,"\t"); print}' | awk '{$11 = $10; print}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$11); print}' |\
  awk -v paralogs=$paralogs 'BEGIN { OFS = "\t"; }; { if ($5 == 0 || paralogs > 1) $5 = 60; else $5 = $5; }; 1' | awk '!( $6 == "*" )' | cat ../preprocess/${i%.f*}_heading.sam - > ../preprocess/${i%.f*}_${ref3%.f*}.sam
	rm ../preprocess/${i%.f*}_exp.sam ../preprocess/${i%.f*}_chrpos.txt ../preprocess/${i%.f*}_uniqsamp.sam ../preprocess/${i%.f*}_uniq.sam ../preprocess/${i%.f*}_heading.sam

  awk '{if ($4=="@HD" || $4=="@SQ" || $4=="@PG"); print $0}' ../preprocess/${i%.f*}_2del.sam > ../preprocess/${i%.f*}_3del.sam  && \
  grep -vwE "(@HD|@SQ|@PG)" ../preprocess/${i%.f*}_2del.sam | awk -v paralogs=$paralogs '{ if ($1 >= 1 && $1 <= paralogs && $2 >= 1 && $2 <= paralogs && $3 == 0); print $0}' >> ../preprocess/${i%.f*}_3del.sam   && \
  awk '!($1=$2=$3="")' ../preprocess/${i%.f*}_3del.sam | awk '{$1=$1};1' > ../preprocess/${i%.f*}_${ref1%.f*}_${ref2%.f*}.sam  && \
	awk '/@HD/ || /@SQ/{print}' ../preprocess/${i%.f*}_${ref1%.f*}_${ref2%.f*}.sam | awk '{gsub(/ /,"\t"); print}' > ../preprocess/${i%.f*}_heading.sam
  awk '!/^@/ { print }' ../preprocess/${i%.f*}_${ref1%.f*}_${ref2%.f*}.sam | awk '{gsub(/^.*-/,"",$1); print $0}' | awk '{gsub(/^.*_/,"",$2); print $0}' > ../preprocess/${i%.f*}_uniq.sam

	calcploidy=$(( ploidy_ref1 + ploidy_ref2 ))
	unds=$(( unbiased_downsample * calcploidy))
	awk -F' ' 'BEGIN{FS=OFS=" "} {if (a[$3,$4]<$1) {a[$3,$4]=$1; data[$3,$4]=$0}} END{for (i in a) print data[i]}' ../preprocess/${i%.f*}_uniq.sam | \
	awk '{print $1,$3,$4}' | grep -v '*' > ../preprocess/${i%.f*}_chrpos.txt
	awk '{print $1,$0}' ../preprocess/${i%.f*}_uniq.sam |
	awk -F" " 'BEGIN{FS=OFS=" "} {if (NR==FNR) {a[$2FS$3]=$1; next} if ($4FS$5 in a) {$1=a[$4FS$5]} print}' ../preprocess/${i%.f*}_chrpos.txt - | \
	awk '!( $4 == "*" )' | awk '{print $1,$0}' | awk -v unds=$unds 'BEGIN{FS=OFS=" "} { if($2 > unds) {printf("%.3f", $1=$1/unds)} else {printf($1=$3);} $1="" }1' | \
	awk -v unds=$unds 'BEGIN{FS=OFS=" "} { if($2 > unds) {printf("%.f", $1=$3/$1)} else {printf($1=$3);} $1="" }1' | \
	awk '!( $1 < "1" )' |  awk '!($2=$3="")' > ../preprocess/${i%.f*}_uniqsamp.sam

  for j in $(LC_ALL=C; sort -n -k1,1 ../preprocess/${i%.f*}_uniqsamp.sam | grep -v @ | awk '{print $1}' | uniq); do
  	awk -v n="^${j}" '$0~n{print $0}' ../preprocess/${i%.f*}_uniqsamp.sam | awk -v n="$j" '{for(i=0;i<n;i++) print}' >> ../preprocess/${i%.f*}_exp.sam
  done; wait
  awk '{print "seq"NR"_"$0}' ../preprocess/${i%.f*}_exp.sam | awk '{gsub(/ /,"\t"); print}' | awk '{$11 = $10; print}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$11); print}' |\
  awk -v paralogs=$paralogs 'BEGIN { OFS = "\t"; }; { if ($5 == 0 || paralogs > 1) $5 = 60; else $5 = $5; }; 1' | awk '!( $6 == "*" )' | cat ../preprocess/${i%.f*}_heading.sam - > ../preprocess/${i%.f*}_${ref1%.f*}_${ref2%.f*}.sam
	rm ../preprocess/${i%.f*}_exp.sam ../preprocess/${i%.f*}_chrpos.txt ../preprocess/${i%.f*}_uniqsamp.sam ../preprocess/${i%.f*}_uniq.sam ../preprocess/${i%.f*}_heading.sam

  awk '{if ($4=="@HD" || $4=="@SQ" || $4=="@PG"); print $0}' ../preprocess/${i%.f*}_2del.sam > ../preprocess/${i%.f*}_3del.sam  && \
  grep -vwE "(@HD|@SQ|@PG)" ../preprocess/${i%.f*}_2del.sam | awk -v paralogs=$paralogs '{ if ($1 >= 1 && $1 <= paralogs && $2 == 0 && $3 >= 1 && $3 <= paralogs); print $0}' >> ../preprocess/${i%.f*}_3del.sam   && \
  awk '!($1=$2=$3="")' ../preprocess/${i%.f*}_3del.sam | awk '{$1=$1};1' > ../preprocess/${i%.f*}_${ref1%.f*}_${ref3%.f*}.sam  && \
	awk '/@HD/ || /@SQ/{print}' ../preprocess/${i%.f*}_${ref1%.f*}_${ref3%.f*}.sam | awk '{gsub(/ /,"\t"); print}' > ../preprocess/${i%.f*}_heading.sam
  awk '!/^@/ { print }' ../preprocess/${i%.f*}_${ref1%.f*}_${ref3%.f*}.sam | awk '{gsub(/^.*-/,"",$1); print $0}' | awk '{gsub(/^.*_/,"",$2); print $0}' > ../preprocess/${i%.f*}_uniq.sam

	calcploidy=$(( ploidy_ref1 + ploidy_ref3 ))
	unds=$(( unbiased_downsample * calcploidy))
	awk -F' ' 'BEGIN{FS=OFS=" "} {if (a[$3,$4]<$1) {a[$3,$4]=$1; data[$3,$4]=$0}} END{for (i in a) print data[i]}' ../preprocess/${i%.f*}_uniq.sam | \
	awk '{print $1,$3,$4}' | grep -v '*' > ../preprocess/${i%.f*}_chrpos.txt
	awk '{print $1,$0}' ../preprocess/${i%.f*}_uniq.sam |
	awk -F" " 'BEGIN{FS=OFS=" "} {if (NR==FNR) {a[$2FS$3]=$1; next} if ($4FS$5 in a) {$1=a[$4FS$5]} print}' ../preprocess/${i%.f*}_chrpos.txt - | \
	awk '!( $4 == "*" )' | awk '{print $1,$0}' | awk -v unds=$unds 'BEGIN{FS=OFS=" "} { if($2 > unds) {printf("%.3f", $1=$1/unds)} else {printf($1=$3);} $1="" }1' | \
	awk -v unds=$unds 'BEGIN{FS=OFS=" "} { if($2 > unds) {printf("%.f", $1=$3/$1)} else {printf($1=$3);} $1="" }1' | \
	awk '!( $1 < "1" )' |  awk '!($2=$3="")' > ../preprocess/${i%.f*}_uniqsamp.sam

  for j in $(LC_ALL=C; sort -n -k1,1 ../preprocess/${i%.f*}_uniqsamp.sam | grep -v @ | awk '{print $1}' | uniq); do
  	awk -v n="^${j}" '$0~n{print $0}' ../preprocess/${i%.f*}_uniqsamp.sam | awk -v n="$j" '{for(i=0;i<n;i++) print}' >> ../preprocess/${i%.f*}_exp.sam
  done; wait
  awk '{print "seq"NR"_"$0}' ../preprocess/${i%.f*}_exp.sam | awk '{gsub(/ /,"\t"); print}' | awk '{$11 = $10; print}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$11); print}' |\
  awk -v paralogs=$paralogs 'BEGIN { OFS = "\t"; }; { if ($5 == 0 || paralogs > 1) $5 = 60; else $5 = $5; }; 1' | awk '!( $6 == "*" )' | cat ../preprocess/${i%.f*}_heading.sam - > ../preprocess/${i%.f*}_${ref1%.f*}_${ref3%.f*}.sam
	rm ../preprocess/${i%.f*}_exp.sam ../preprocess/${i%.f*}_chrpos.txt ../preprocess/${i%.f*}_uniqsamp.sam ../preprocess/${i%.f*}_uniq.sam ../preprocess/${i%.f*}_heading.sam

  awk '{if ($4=="@HD" || $4=="@SQ" || $4=="@PG"); print $0}' ../preprocess/${i%.f*}_2del.sam > ../preprocess/${i%.f*}_3del.sam  && \
  grep -vwE "(@HD|@SQ|@PG)" ../preprocess/${i%.f*}_2del.sam | awk -v paralogs=$paralogs '{ if ($1 == 0 && $2 >= 1 && $2 <= paralogs && $3 >= 1 && $3 <= paralogs); print $0}' >> ../preprocess/${i%.f*}_3del.sam   && \
  awk '!($1=$2=$3="")' ../preprocess/${i%.f*}_3del.sam | awk '{$1=$1};1' > ../preprocess/${i%.f*}_${ref2%.f*}_${ref3%.f*}.sam  && \
	awk '/@HD/ || /@SQ/{print}' ../preprocess/${i%.f*}_${ref2%.f*}_${ref3%.f*}.sam | awk '{gsub(/ /,"\t"); print}' > ../preprocess/${i%.f*}_heading.sam
  awk '!/^@/ { print }' ../preprocess/${i%.f*}_${ref2%.f*}_${ref3%.f*}.sam | awk '{gsub(/^.*-/,"",$1); print $0}' | awk '{gsub(/^.*_/,"",$2); print $0}' > ../preprocess/${i%.f*}_uniq.sam

	calcploidy=$(( ploidy_ref2 + ploidy_ref3 ))
	awk -F' ' 'BEGIN{FS=OFS=" "} {if (a[$3,$4]<$1) {a[$3,$4]=$1; data[$3,$4]=$0}} END{for (i in a) print data[i]}' ../preprocess/${i%.f*}_uniq.sam | \
	awk '{print $1,$3,$4}' | grep -v '*' > ../preprocess/${i%.f*}_chrpos.txt
	awk '{print $1,$0}' ../preprocess/${i%.f*}_uniq.sam |
	awk -F" " 'BEGIN{FS=OFS=" "} {if (NR==FNR) {a[$2FS$3]=$1; next} if ($4FS$5 in a) {$1=a[$4FS$5]} print}' ../preprocess/${i%.f*}_chrpos.txt - | \
	awk '!( $4 == "*" )' | awk '{print $1,$0}' | awk -v unds=$unds 'BEGIN{FS=OFS=" "} { if($2 > unds) {printf("%.3f", $1=$1/unds)} else {printf($1=$3);} $1="" }1' | \
	awk -v unds=$unds 'BEGIN{FS=OFS=" "} { if($2 > unds) {printf("%.f", $1=$3/$1)} else {printf($1=$3);} $1="" }1' | \
	awk '!( $1 < "1" )' |  awk '!($2=$3="")' > ../preprocess/${i%.f*}_uniqsamp.sam

  for j in $(LC_ALL=C; sort -n -k1,1 ../preprocess/${i%.f*}_uniqsamp.sam| grep -v @ | awk '{print $1}' | uniq); do
  	awk -v n="^${j}" '$0~n{print $0}' ../preprocess/${i%.f*}_uniqsamp.sam | awk -v n="$j" '{for(i=0;i<n;i++) print}' >> ../preprocess/${i%.f*}_exp.sam
  done; wait
  awk '{print "seq"NR"_"$0}' ../preprocess/${i%.f*}_exp.sam | awk '{gsub(/ /,"\t"); print}' | awk '{$11 = $10; print}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$11); print}' |\
  awk -v paralogs=$paralogs 'BEGIN { OFS = "\t"; }; { if ($5 == 0 || paralogs > 1) $5 = 60; else $5 = $5; }; 1' | awk '!( $6 == "*" )' | cat ../preprocess/${i%.f*}_heading.sam - > ../preprocess/${i%.f*}_${ref2%.f*}_${ref3%.f*}.sam
	rm ../preprocess/${i%.f*}_exp.sam ../preprocess/${i%.f*}_chrpos.txt ../preprocess/${i%.f*}_uniqsamp.sam ../preprocess/${i%.f*}_uniq.sam ../preprocess/${i%.f*}_heading.sam


  awk '{if ($4=="@HD" || $4=="@SQ" || $4=="@PG"); print $0}' ../preprocess/${i%.f*}_2del.sam > ../preprocess/${i%.f*}_3del.sam  && \
  grep -vwE "(@HD|@SQ|@PG)" ../preprocess/${i%.f*}_2del.sam | awk -v paralogs=$paralogs '{ if ($1 >= 1 && $1 <= paralogs && $2 >= 1 && $2 <= paralogs && $3 >= 1 && $3 <= paralogs); print $0}' >> ../preprocess/${i%.f*}_3del.sam   && \
  awk '!($1=$2=$3="")' ../preprocess/${i%.f*}_3del.sam | awk '{$1=$1};1' > ../preprocess/${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam  && \
	awk '/@HD/ || /@SQ/{print}' ../preprocess/${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam | awk '{gsub(/ /,"\t"); print}' > ../preprocess/${i%.f*}_heading.sam
  awk '!/^@/ { print }' ../preprocess/${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam | awk '{gsub(/^.*-/,"",$1); print $0}' | awk '{gsub(/^.*_/,"",$2); print $0}' > ../preprocess/${i%.f*}_uniq.sam

	calcploidy=$(( ploidy_ref1 + ploidy_ref2 + ploidy_ref3 ))
	unds=$(( unbiased_downsample * calcploidy))
	awk -F' ' 'BEGIN{FS=OFS=" "} {if (a[$3,$4]<$1) {a[$3,$4]=$1; data[$3,$4]=$0}} END{for (i in a) print data[i]}' ../preprocess/${i%.f*}_uniq.sam | \
	awk '{print $1,$3,$4}' | grep -v '*' > ../preprocess/${i%.f*}_chrpos.txt
	awk '{print $1,$0}' ../preprocess/${i%.f*}_uniq.sam |
	awk -F" " 'BEGIN{FS=OFS=" "} {if (NR==FNR) {a[$2FS$3]=$1; next} if ($4FS$5 in a) {$1=a[$4FS$5]} print}' ../preprocess/${i%.f*}_chrpos.txt - | \
	awk '!( $4 == "*" )' | awk '{print $1,$0}' | awk -v unds=$unds 'BEGIN{FS=OFS=" "} { if($2 > unds) {printf("%.3f", $1=$1/unds)} else {printf($1=$3);} $1="" }1' | \
	awk -v unds=$unds 'BEGIN{FS=OFS=" "} { if($2 > unds) {printf("%.f", $1=$3/$1)} else {printf($1=$3);} $1="" }1' | \
	awk '!( $1 < "1" )' |  awk '!($2=$3="")' > ../preprocess/${i%.f*}_uniqsamp.sam

  for j in $(LC_ALL=C; sort -n -k1,1 ../preprocess/${i%.f*}_uniqsamp.sam | grep -v @ | awk '{print $1}' | uniq); do
  	awk -v n="^${j}" '$0~n{print $0}' ../preprocess/${i%.f*}_uniqsamp.sam | awk -v n="$j" '{for(i=0;i<n;i++) print}' >> ../preprocess/${i%.f*}_exp.sam
  done; wait
  awk '{print "seq"NR"_"$0}' ../preprocess/${i%.f*}_exp.sam | awk '{gsub(/ /,"\t"); print}' | awk '{$11 = $10; print}' | awk 'BEGIN{OFS="\t"}{gsub(/A|a|C|c|G|g|T|t|N|n/,"I",$11); print}' |\
  awk -v paralogs=$paralogs 'BEGIN { OFS = "\t"; }; { if ($5 == 0 || paralogs > 1) $5 = 60; else $5 = $5; }; 1' | awk '!( $6 == "*" )' | cat ../preprocess/${i%.f*}_heading.sam - > ../preprocess/${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam
	rm ../preprocess/${i%.f*}_exp.sam ../preprocess/${i%.f*}_chrpos.txt ../preprocess/${i%.f*}_uniqsamp.sam ../preprocess/${i%.f*}_uniq.sam ../preprocess/${i%.f*}_heading.sam


  declare -a arr=("${i%.f*}_${ref1%.f*}.sam" "${i%.f*}_${ref2%.f*}.sam" "${i%.f*}_${ref3%.f*}.sam" \
  "${i%.f*}_${ref1%.f*}_${ref2%.f*}.sam" "${i%.f*}_${ref1%.f*}_${ref3%.f*}.sam" "${i%.f*}_${ref2%.f*}_${ref3%.f*}.sam" \
  "${i%.f*}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.sam") && \
  cd $projdir
  cd preprocess
  for j in "${arr[@]}"; do
        $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard SortSam I=$j O=${j%.sam*}.bam SORT_ORDER=coordinate && \
        $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard BuildBamIndex INPUT=${j%.sam*}.bam && \
        $java $Xmxp -XX:ParallelGCThreads=$prepthreads -jar $picard AddOrReplaceReadGroups I=${j%.sam*}.bam O=${j%.sam*}_precall.bam RGLB=${i%.f*} RGPL=illumina VALIDATION_STRINGENCY=LENIENT RGPU=run RGSM=${i%.f*} && \
        $samtools index ${j%.sam*}_precall.bam
  done
  ls ${i%.f*}_* | grep -v precall | xargs rm && \
  cd ../samples
	wait )&
	if [[ $(jobs -r -p | wc -l) -ge $prepN ]]; then
		wait
	fi
done
wait
cat ../alignment_summaries/*_summ.txt > ../alignment_summaries/alignment_summaries_unique_reads.txt; rm -r ../alignment_summaries/*_summ.txt
# Total number of reads per samples
awk '/###---/ || /QC-passed/{print}' ../alignment_summaries/alignment_summaries_unique_reads.txt | cut -d\+ -f1 | tr -d '\n' | \
awk  'gsub(/---###/, "\t", $0)' | awk  'gsub(/###---/, "", $0)' | tr ' ' '\n' > ../alignment_summaries/total_unique_reads.txt
# Total number of mapped reads per samples
cat ../alignment_summaries/alignment_summaries_unique_reads.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /0_mapped/{print}' |\
tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
awk 'gsub("\\+0mapped", "\t", $0)' | tr ":" "\t" | cut -d\: -f1 | awk 'gsub(/ /, "\t")' > ../alignment_summaries/total_unique_reads_mapped.txt
# Total number of mapped reads per samples
cat ../alignment_summaries/alignment_summaries_unique_reads.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /properly_paired/{print}' |\
tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
awk 'gsub("\\+0properlypaired", "\t", $0)' | tr ":" "\t" | cut -d\: -f1 | awk 'gsub(/ /, "\t")' > ../alignment_summaries/total_unique_reads_paired.txt
echo -e "Samples\tTotal\tMapped\tPerc_Mapped\tPE_Mapped\t%_PE_Mapped" > ../alignment_summaries/summary_precall.txt
awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' ../alignment_summaries/total_unique_reads_mapped.txt  ../alignment_summaries/total_unique_reads.txt  | \
awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' ../alignment_summaries/total_unique_reads_paired.txt - | awk '{gsub(/ /,"\t"); print $0}' | \
cat ../alignment_summaries/summary_precall.txt - > ../alignment_summaries/Tabulated_Alignment_Unique_Read_Summaries.txt
awk '{gsub(/\t/,","); print $0}' ../alignment_summaries/Tabulated_Alignment_Unique_Read_Summaries.txt > ../alignment_summaries/Tabulated_Alignment_Unique_Read_Summaries.csv
rm ../alignment_summaries/total* ../alignment_summaries/summary_precall.txt
rm ../samples/metrics.txt ../preprocess/metrics.txt

cd $projdir/alignment_summaries

touch refgenome_paralogs_${ref1%.f*}.txt
for par in refgenome_paralogs_*_${ref1%.f*}.txt;then
	cat refgenome_paralogs_${ref1%.f*}.txt $par | awk '!visited[$0]++' > temp_par.txt
	mv temp_par.txt refgenome_paralogs_${ref1%.f*}.txt
fi
rm refgenome_paralogs_*_${ref1%.f*}.txt
awk '{gsub(/~/,"\t"); print $0}' refgenome_paralogs_${ref1%.f*}.txt | awk 'BEGIN{OFS="\t"; }; {if($2==0) $1 = "multilocus"; else $1 = $1; }; 1' | \
awk 'BEGIN{OFS="\t"; };{print $3,$4,$1}' | awk '$3>max[$1,$2]{max[$1,$2]=$3; row[$1,$2]=$0} END{for (i in row) print row[i]}' | \
awk '{print $0,substr($2, 1, length($2)-2)}' | awk '$3>max[$1,$4]{max[$1,$4]=$3; row[$1,$4]=$0} END{for (i in row) print row[i]}' > temp.txt
awk '{print $1"\t"$2"\t"$3}' temp.txt | awk '!/CHROM/' | cat <(printf "CHROM\tPOS\tnloci\n") - > refgenome_paralogs_${ref1%.f*}.txt

touch refgenome_paralogs_${ref2%.f*}.txt
for par in refgenome_paralogs_*_${ref2%.f*}.txt;then
	cat refgenome_paralogs_${ref2%.f*}.txt $par | awk '!visited[$0]++' > temp_par.txt
	mv temp_par.txt refgenome_paralogs_${ref2%.f*}.txt
fi
rm refgenome_paralogs_*_${ref2%.f*}.txt
awk '{gsub(/~/,"\t"); print $0}' refgenome_paralogs_${ref2%.f*}.txt | awk 'BEGIN{OFS="\t"; }; {if($2==0) $1 = "multilocus"; else $1 = $1; }; 1' | \
awk 'BEGIN{OFS="\t"; };{print $3,$4,$1}' | awk '$3>max[$1,$2]{max[$1,$2]=$3; row[$1,$2]=$0} END{for (i in row) print row[i]}' | \
awk '{print $0,substr($2, 1, length($2)-2)}' | awk '$3>max[$1,$4]{max[$1,$4]=$3; row[$1,$4]=$0} END{for (i in row) print row[i]}' > temp.txt
awk '{print $1"\t"$2"\t"$3}' temp.txt | awk '!/CHROM/' | cat <(printf "CHROM\tPOS\tnloci\n") - > refgenome_paralogs_${ref2%.f*}.txt

touch refgenome_paralogs_${ref3%.f*}.txt
for par in refgenome_paralogs_*_${ref3%.f*}.txt;then
	cat refgenome_paralogs_${ref3%.f*}.txt $par | awk '!visited[$0]++' > temp_par.txt
	mv temp_par.txt refgenome_paralogs_${ref3%.f*}.txt
fi
rm refgenome_paralogs_*_${ref3%.f*}.txt
awk '{gsub(/~/,"\t"); print $0}' refgenome_paralogs_${ref3%.f*}.txt | awk 'BEGIN{OFS="\t"; }; {if($2==0) $1 = "multilocus"; else $1 = $1; }; 1' | \
awk 'BEGIN{OFS="\t"; };{print $3,$4,$1}' | awk '$3>max[$1,$2]{max[$1,$2]=$3; row[$1,$2]=$0} END{for (i in row) print row[i]}' | \
awk '{print $0,substr($2, 1, length($2)-2)}' | awk '$3>max[$1,$4]{max[$1,$4]=$3; row[$1,$4]=$0} END{for (i in row) print row[i]}' > temp.txt
awk '{print $1"\t"$2"\t"$3}' temp.txt | awk '!/CHROM/' | cat <(printf "CHROM\tPOS\tnloci\n") - > refgenome_paralogs_${ref3%.f*}.txt

awk 'FNR==1 && NR!=1 { while (/^CHROM/) getline; }1 {print}' refgenome_paralogs_*.txt > refgenome_paralogs.txt
rm refgenome_paralogs_*


}
cd $projdir
dirpreprocess=./preprocess
if [ "$walkaway" == false ]; then
	echo -e "${magenta}- Do you want to perform read alignments and alignment post-processing? ${white}\n"
	read -p "- y(YES) or n(NO) " -t 36000 -n 1 -r
	if [[ ! $REPLY =~ ^[Yy]$ ]]; then
		printf '\n'
		echo -e "${magenta}- skipping read alignments and alignment post-processing ${white}\n"
	else
		printf '\n'
		echo -e "${magenta}- performing read alignments and alignment post-processing ${white}\n"
		dirpreprocess=./preprocess
		if [ "$(ls -A $dirpreprocess)" ]; then
			echo -e "${magenta}- \n- preprocess folder should be empty, exiting pipeline ${white}\n"
			sleep 10 && exit 1
		fi
		time main &>> log.out
	fi
fi
if [ "$walkaway" == true ]; then
	if [ "$alignments" == 1 ]; then
		echo -e "${magenta}- performing read alignments and alignment post-processing ${white}\n"
		dirpreprocess=./preprocess
		if [ "$(ls -A $dirpreprocess)" ]; then
			echo -e "${magenta}- \n- preprocess folder should be empty, exiting pipeline ${white}\n"
			sleep 10 && exit 1
		fi
		time main &>> log.out
	else
		echo -e "${magenta}- skipping read alignments and alignment post-processing ${white}\n"
	fi
fi




######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${yellow}\n- Performing Variant Calling with GATK HaplotypeCaller\n${blue}##############################################################################${white}\n"
main () {
	cd $projdir
	ncohorts=1
	cd preprocess
	mkdir processed
	if [[ -z "$(ls -A ./processed/)" ]]; then
		:
	else
		mv ./processed/*_precall.bam* ./
	fi

######################

echo -e "${magenta}- performing SNP calling across entire genome (subgenome-1, 2, and 3) ${white}\n"
for i in $(ls -S *_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_precall.bam); do (
	$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads"  HaplotypeCaller -R ../refgenomes/panref.fasta -I $i -ploidy $ploidy -O ../snpcall/${i%_precall.bam}.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases true --max-reads-per-alignment-start "$(( biased_downsample * ploidy ))" --max-num-haplotypes-in-population "$((ploidy * maxHaplotype))"
	) &
	if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		wait
	fi
done
wait
cd ../snpcall
if [[ "$ncohorts" -eq  1 ]]; then
	cz=$(ls *.g.vcf.gz | wc -l)
fi
if [[ "$ncohorts" -gt  1 ]]; then
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
	for i in $(ls *_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}.g.vcf.gz); do
		k="${j} ${i}"; input="${input} ${k}"
	done
	Get_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ../../refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' )
	for selchr in $Get_Chromosome; do (
		$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" GenomicsDBImport $input -L $selchr --genomicsdb-workspace-path ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_"${selchr}"_raw
		$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" GenotypeGVCFs -R ../../refgenomes/panref.fasta -L $selchr -V gendb://${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_"${selchr}"_raw -O ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_"${selchr}"_raw.vcf.gz
		rm -r ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_"${selchr}"_raw
		if LC_ALL=C gzip -l ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_"${selchr}"_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
			:
		else
			rm ../cohorts*/${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_*_raw.vcf.gz*
			rm ../cohorts*/${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_raw.vcf*
			rm ../${pop}_${ref1%.f*}_$if [[ "$ncohorts" -gt 0 ]]; then{ref2%.f*}_${ref3%.f*}_${ploidy}x_raw_cohorts*.vcf*
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
	for g in $(ls ${pop}_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_${ploidy}x_*_raw.vcf.gz); do (
	  gunzip $g )&
	  if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
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

rm -r cohorts*
rm *cohorts*
cd ${projdir}/preprocess
mv *_${ref1%.f*}_${ref2%.f*}_${ref3%.f*}_precall* ./processed/

######################

echo -e "${magenta}- performing SNP calling on subgenomes 1 and 2 ${white}\n"
calcploidy=$(( ploidy_ref1 + ploidy_ref2 ))
awk 'NR>1{print $2,"\t",$3}' ../refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}'| awk -v pat1=${ref1%.f*} -v pat2=${ref2%.f*} 'pat1 || pat2' > ../refgenomes/intervals.list
for i in $(ls -S *_${ref1%.f*}_${ref2%.f*}_precall.bam); do (
	$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ../refgenomes/panref.fasta -L ../refgenomes/intervals.list -I $i -ploidy ${calcploidy} -O ../snpcall/${i%_precall.bam}.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases true --max-reads-per-alignment-start "$(( biased_downsample * calcploidy ))" --max-num-haplotypes-in-population "$((calcploidy * maxHaplotype))"
	) &
	if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		wait
	fi
done
wait
rm ../refgenomes/intervals.list
cd ../snpcall
if [[ "$ncohorts" -eq  1 ]]; then
	cz=$(ls *.g.vcf.gz | wc -l)
fi
if [[ "$ncohorts" -gt  1 ]]; then
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
	for i in $(ls *_${ref1%.f*}_${ref2%.f*}.g.vcf.gz); do
		k="${j} ${i}"; input="${input} ${k}"
	done
	Get_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ../../refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' | awk -v pat1=${ref1%.f*} -v pat2=${ref2%.f*} 'pat1 || pat2')
	for selchr in $Get_Chromosome; do (
		$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" GenomicsDBImport $input -L $selchr --genomicsdb-workspace-path ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_"${selchr}"_raw
		$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" GenotypeGVCFs -R ../../refgenomes/panref.fasta -L $selchr -V gendb://${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_"${selchr}"_raw -O ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_"${selchr}"_raw.vcf.gz
		rm -r ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_"${selchr}"_raw
		if LC_ALL=C gzip -l ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_"${selchr}"_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
			:
		else
			rm ../cohorts*/${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_*_raw.vcf.gz*
			rm ../cohorts*/${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf*
			rm ../${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw_cohorts*.vcf*
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
	for g in $(ls ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_*_raw.vcf.gz); do (
	  gunzip $g )&
	  if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
	  wait
	  fi
	done
	grep -h '^#' ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt
	cat ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_*_raw.vcf | awk '!/^#/' > all.vcf
	cat vcf_header.txt all.vcf > ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf
	rm vcf_header.txt all.vcf
	rm ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_*_raw.vcf ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_*_raw.vcf.gz.tbi
	$bcftools view -I ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf -O z -o ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf.gz
	$bcftools index ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf

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
	gunzip ${pop}_${ref1%.f*}_${ref2%.f*}_${calcploidy}x_raw.vcf.gz
fi

rm -r cohorts*
rm *cohorts*
cd ${projdir}/preprocess
mv *_${ref1%.f*}_${ref2%.f*}_precall* ./processed/

######################

echo -e "${magenta}- performing SNP calling on subgenomes 1 and 3 ${white}\n"
calcploidy=$(( ploidy_ref1 + ploidy_ref3 ))
awk 'NR>1{print $2,"\t",$3}' ../refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}'| awk -v pat1=${ref1%.f*} -v pat2=${ref3%.f*} 'pat1 || pat2' > ../refgenomes/intervals.list
for i in $(ls -S *_${ref1%.f*}_${ref3%.f*}_precall.bam); do (
	$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ../refgenomes/panref.fasta -L ../refgenomes/intervals.list -I $i -ploidy ${calcploidy} -O ../snpcall/${i%_precall.bam}.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases true --max-reads-per-alignment-start "$(( biased_downsample * calcploidy ))" --max-num-haplotypes-in-population "$((calcploidy * maxHaplotype))"
	) &
	if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		wait
	fi
done
wait
rm ../refgenomes/intervals.list
cd ../snpcall
if [[ "$ncohorts" -eq  1 ]]; then
	cz=$(ls *.g.vcf.gz | wc -l)
fi
if [[ "$ncohorts" -gt  1 ]]; then
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
	for i in $(ls *_${ref1%.f*}_${ref3%.f*}.g.vcf.gz); do
		k="${j} ${i}"; input="${input} ${k}"
	done
	Get_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ../../refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' | awk -v pat1=${ref1%.f*} -v pat2=${ref3%.f*} 'pat1 || pat2')
	for selchr in $Get_Chromosome; do (
		$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" GenomicsDBImport $input -L $selchr --genomicsdb-workspace-path ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_"${selchr}"_raw
		$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" GenotypeGVCFs -R ../../refgenomes/panref.fasta -L $selchr -V gendb://${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_"${selchr}"_raw -O ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_"${selchr}"_raw.vcf.gz
		rm -r ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_*_raw
		if LC_ALL=C gzip -l ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_"${selchr}"_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
			:
		else
			rm ../cohorts*/${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf.gz*
			rm ../cohorts*/${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf*
			rm ../${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw_cohorts*.vcf*
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
	for g in $(ls ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf.gz); do (
	  gunzip $g )&
	  if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
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
	gunzip ${pop}_${ref1%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf.gz
fi

rm -r cohorts*
rm *cohorts*
cd ${projdir}/preprocess
mv *_${ref1%.f*}_${ref3%.f*}_precall* ./processed/

######################

echo -e "${magenta}- performing SNP calling on subgenomes 2 and 3 ${white}\n"
calcploidy=$(( ploidy_ref2 + ploidy_ref3 ))
awk 'NR>1{print $2,"\t",$3}' ../refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}'| awk -v pat1=${ref2%.f*} -v pat2=${ref3%.f*} 'pat1 || pat2' > ../refgenomes/intervals.list
for i in $(ls -S *_${ref2%.f*}_${ref3%.f*}_precall.bam); do (
	$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ../refgenomes/panref.fasta -L ../refgenomes/intervals.list -I $i -ploidy ${calcploidy} -O ../snpcall/${i%_precall.bam}.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases true --max-reads-per-alignment-start "$(( biased_downsample * calcploidy ))" --max-num-haplotypes-in-population "$((calcploidy * maxHaplotype))"
	) &
	if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		wait
	fi
done
wait
rm  ../refgenomes/intervals.list
cd ../snpcall
if [[ "$ncohorts" -eq  1 ]]; then
	cz=$(ls *.g.vcf.gz | wc -l)
fi
if [[ "$ncohorts" -gt  1 ]]; then
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
	for i in $(ls *_${ref2%.f*}_${ref3%.f*}.g.vcf.gz); do
		k="${j} ${i}"; input="${input} ${k}"
	done
	Get_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ../../refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' | awk -v pat1=${ref2%.f*} -v pat2=${ref3%.f*} 'pat1 || pat2')
	for selchr in $Get_Chromosome; do (
		$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" GenomicsDBImport $input -L $selchr --genomicsdb-workspace-path ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_"${selchr}"_raw
		$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" GenotypeGVCFs -R ../../refgenomes/panref.fasta -L $selchr -V gendb://${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_"${selchr}"_raw -O ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_"${selchr}"_raw.vcf.gz
		rm -r ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_*_raw
		if LC_ALL=C gzip -l ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_"${selchr}"_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
			:
		else
			rm ../cohorts*/${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf.gz*
			rm ../cohorts*/${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf*
			rm ../${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw_cohorts*.vcf*
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
	for g in $(ls ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_*_raw.vcf.gz); do (
	  gunzip $g )&
	  if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
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
	gunzip ${pop}_${ref2%.f*}_${ref3%.f*}_${calcploidy}x_raw.vcf.gz
fi

rm -r cohorts*
rm *cohorts*
cd ${projdir}/preprocess
mv *_${ref2%.f*}_${ref3%.f*}_precall* ./processed/

######################

echo -e "${magenta}- performing SNP calling on subgenome-1 ${white}\n"
awk 'NR>1{print $2,"\t",$3}' ../refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}'| awk -v pat=${ref1%.f*} '$0 ~ pat' > ../refgenomes/intervals.list
for i in $(ls -S *_${ref1%.f*}_precall.bam); do (
	$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ../refgenomes/panref.fasta -L ../refgenomes/intervals.list -I $i -ploidy ${ploidy_ref1} -O ../snpcall/${i%_precall.bam}.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases true --max-reads-per-alignment-start "$(( biased_downsample * ploidy_ref1 ))" --max-num-haplotypes-in-population "$((ploidy_ref1 * maxHaplotype))"
	) &
	if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		wait
	fi
done
wait
rm ../refgenomes/intervals.list
cd ../snpcall
if [[ "$ncohorts" -eq  1 ]]; then
	cz=$(ls *.g.vcf.gz | wc -l)
fi
if [[ "$ncohorts" -gt  1 ]]; then
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
	for i in $(ls *_${ref1%.f*}.g.vcf.gz); do
		k="${j} ${i}"; input="${input} ${k}"
	done
	Get_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ../../refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' | awk -v pat=${ref1%.f*} '$0 ~ pat')
	for selchr in $Get_Chromosome; do (
		$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" GenomicsDBImport $input -L $selchr --genomicsdb-workspace-path ${pop}_${ref1%.f*}_${ploidy_ref1}x_"${selchr}"_raw
		$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" GenotypeGVCFs -R ../../refgenomes/panref.fasta -L $selchr -V gendb://${pop}_${ref1%.f*}_${ploidy_ref1}x_"${selchr}"_raw -O ${pop}_${ref1%.f*}_${ploidy_ref1}x_"${selchr}"_raw.vcf.gz
		rm -r ${pop}_${ref1%.f*}_${ploidy_ref1}x_*_raw
		if LC_ALL=C gzip -l ${pop}_${ref1%.f*}_${ploidy_ref1}x_"${selchr}"_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
			:
		else
			rm ../cohorts*/${pop}_${ref1%.f*}_${ploidy_ref1}x_*_raw.vcf.gz*
			rm ../cohorts*/${pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf*
			rm ../${pop}_${ref1%.f*}_${ploidy_ref1}x_raw_cohorts*.vcf*
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
	for g in $(ls ${pop}_${ref1%.f*}_${ploidy_ref1}x_*_raw.vcf.gz); do (
	  gunzip $g )&
	  if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
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
	$bcftools merge *cohorts*.vcf.gz --force-samples -m all > ${pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf
else
	cp *cohorts*.vcf.gz ${pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf.gz
	gunzip ${pop}_${ref1%.f*}_${ploidy_ref1}x_raw.vcf.gz
fi

rm -r cohorts*
rm *cohorts*
cd ${projdir}/preprocess
mv *_${ref1%.f*}_precall* ./processed/

######################

echo -e "${magenta}- performing SNP calling on subgenome-2 ${white}\n"
awk 'NR>1{print $2,"\t",$3}' ../refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}'| awk -v pat=${ref2%.f*} '$0 ~ pat' > ../refgenomes/intervals.list
for i in $(ls -S *_${ref2%.f*}_precall.bam); do (
	$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ../refgenomes/panref.fasta -L ../refgenomes/intervals.list -I $i -ploidy ${ploidy_ref2} -O ../snpcall/${i%_precall.bam}.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases true --max-reads-per-alignment-start "$(( biased_downsample * ploidy_ref2 ))" --max-num-haplotypes-in-population "$((ploidy_ref2 * maxHaplotype))"
	) &
	if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		wait
	fi
done
wait
rm ../refgenomes/intervals.list
cd ../snpcall
if [[  "$ncohorts" -eq  1 ]]; then
	cz=$(ls *.g.vcf.gz | wc -l)
fi
if [[ "$ncohorts" -gt  1 ]]; then
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
	for i in $(ls *_${ref2%.f*}.g.vcf.gz); do
		k="${j} ${i}"; input="${input} ${k}"
	done
	Get_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ../../refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' | awk -v pat=${ref2%.f*} '$0 ~ pat')
	for selchr in $Get_Chromosome; do (
		$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" GenomicsDBImport $input -L $selchr --genomicsdb-workspace-path ${pop}_${ref2%.f*}_${ploidy_ref2}x_"${selchr}"_raw
		$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" GenotypeGVCFs -R ../../refgenomes/panref.fasta -L $selchr -V gendb://${pop}_${ref2%.f*}_${ploidy_ref2}x_"${selchr}"_raw -O ${pop}_${ref2%.f*}_${ploidy_ref2}x_"${selchr}"_raw.vcf.gz
		rm -r ${pop}_${ref2%.f*}_${ploidy_ref2}x_*_raw
		if LC_ALL=C gzip -l ${pop}_${ref2%.f*}_${ploidy_ref2}x_"${selchr}"_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
			:
		else
			rm ../cohorts*/${pop}_${ref2%.f*}_${ploidy_ref2}x_*_raw.vcf.gz*
			rm ../cohorts*/${pop}_${ref2%.f*}_${ploidy_ref2}x_raw.vcf*
			rm ../${pop}_${ref2%.f*}_${ploidy_ref2}x_raw_cohorts*.vcf*
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
	for g in $(ls ${pop}_${ref2%.f*}_${ploidy_ref2}x_*_raw.vcf.gz); do (
	  gunzip $g )&
	  if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
	  wait

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

rm -r cohorts*
rm *cohorts*
cd ${projdir}/preprocess
mv *_${ref2%.f*}_precall* ./processed/


######################

echo -e "${magenta}- performing SNP calling on subgenome-3 ${white}\n"
awk 'NR>1{print $2,"\t",$3}' ../refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}'| awk -v pat=${ref3%.f*} '$0 ~ pat' > ../refgenomes/intervals.list
for i in $(ls -S *_${ref3%.f*}_precall.bam); do (
	$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" HaplotypeCaller -R ../refgenomes/panref.fasta -L ../refgenomes/intervals.list -I $i -ploidy ${ploidy_ref3} -O ../snpcall/${i%_precall.bam}.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases true --max-reads-per-alignment-start "$(( biased_downsample * ploidy_ref3 ))" --max-num-haplotypes-in-population "$((ploidy_ref3 * maxHaplotype))"
	) &
	if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		wait
	fi
done
wait
rm ../refgenomes/intervals.list
cd ../snpcall
if [[ "$ncohorts" -eq  1 ]]; then
	cz=$(ls *.g.vcf.gz | wc -l)
fi
if [[ "$ncohorts" -gt  1 ]]; then
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
	for i in $(ls *_${ref3%.f*}.g.vcf.gz); do
		k="${j} ${i}"; input="${input} ${k}"
	done
	Get_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ../../refgenomes/panref.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}' | awk -v pat=${ref3%.f*} '$0 ~ pat')
	for selchr in $Get_Chromosome; do (
		$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" GenomicsDBImport $input -L $selchr --genomicsdb-workspace-path ${pop}_${ref3%.f*}_${ploidy_ref3}x_"${selchr}"_raw
		$GATK --java-options "$Xmxg -XX:ParallelGCThreads=$gthreads" GenotypeGVCFs -R ../../refgenomes/panref.fasta -L $selchr -V gendb://${pop}_${ref3%.f*}_${ploidy_ref3}x_"${selchr}"_raw -O ${pop}_${ref3%.f*}_${ploidy_ref3}x_"${selchr}"_raw.vcf.gz
		rm -r ${pop}_${ref3%.f*}_${ploidy_ref3}x_*_raw
		if LC_ALL=C gzip -l ${pop}_${ref3%.f*}_${ploidy_ref3}x_"${selchr}"_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
			:
		else
			rm ../cohorts*/${pop}_${ref3%.f*}_${ploidy_ref3}x_*_raw.vcf.gz*
			rm ../cohorts*/${pop}_${ref3%.f*}_${ploidy_ref3}x_raw.vcf*
			rm ../${pop}_${ref3%.f*}_${ploidy_ref3}x_raw_cohorts*.vcf*
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
	for g in $(ls ${pop}_${ref3%.f*}_${ploidy_ref3}x_*_raw.vcf.gz); do (
	  gunzip $g )&
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

rm -r cohorts*
rm *cohorts*
cd ${projdir}/preprocess
mv *_${ref3%.f*}_precall* ./processed/


########################
echo -e "${magenta}- keeping *realign.bam & *realign.bam files in ./preprocess/processed/ ${white}\n"
mv ../preprocess/processed/* ../preprocess/
rmdir ../preprocess/processed


#retrieve the original SNP positions from unsplit chromosomes/pseudomolecules.
cd $projdir
cd snpcall
if [[ "$checksplit" -gt 0 ]]; then
	for i in $(ls *.vcf); do
		cp $i ${i%.vcf}_positions_split.vcf
		awk '/^#CHROM/{close("file.vcf"f);f++}{print $0 > "file"f}' $i
		awk 'NR==1{print}' file1 | cat file - > file0.vcf
		for j in $(seq 1 10); do
			awk 'NR>1{print}' file1 | awk '{print $1,"\t",$2,"\t",$0}' | awk '{gsub(/_x0/,"\t"); print}' | \
			awk -v splits=$j -F '\t' 'BEGIN{OFS="\t"} $2 ~ splits {$3=$3+500000000}1' | awk 'BEGIN{OFS="\t"} !($2="")' | awk 'BEGIN{OFS="\t"} !($3="")' | \
			awk 'BEGIN{OFS="\t"} !($3="")' | awk 'BEGIN{OFS="\t"} !($3="")' | awk '{gsub(/\t\t/,"\t"); print }' | \
			sort -k1,1 -k2n,2 | awk NF > file1.vcf
		done
		cat file0.vcf file1.vcf > $i
		rm file file1 file0.vcf file1.vcf
	done

		cd ../alignment_summaries
		i=refgenome_paralogs.txt
		cp $i ${i%.txt}_positions_split.txt
		awk '/^CHROM/{close("file.txt"f);f++}{print $0 > "file"f}' $i
		awk 'NR==1{print}' file1 | cat file - > file0.txt
		for j in $(seq 1 10); do
			awk 'NR>1{print}' file1 | awk '{print $1,"\t",$2,"\t",$0}' | awk '{gsub(/_x0/,"\t"); print}' | \
			awk -v splits=$j -F '\t' 'BEGIN{OFS="\t"} $2 ~ splits {$3=$3+500000000}1' | awk 'BEGIN{OFS="\t"} !($2="")' | awk 'BEGIN{OFS="\t"} !($3="")' | \
			awk 'BEGIN{OFS="\t"} !($3="")' | awk 'BEGIN{OFS="\t"} !($3="")' | awk '{gsub(/\t\t/,"\t"); print }' | \
			sort -k1,1 -k2n,2 | awk NF > file1.txt
		done
		cat file0.txt file1.txt > $i
		rm file file1 file0.txt file1.txt
		cd ../snpcall

else
	:
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
		dirsnpcall=./snpcall/*
		if [ ${#dirsnpcall[@]} -gt 1 ]; then
			echo -e "${magenta}- \n- snpcall folder should be empty or contain only 1 directory, exiting pipeline ${white}\n"
			sleep 10 && exit 1
		fi
		time main &>> log.out
	fi
fi
if [ "$walkaway" == true ]; then
	if [ "$snp_calling" == 1 ]; then
		echo -e "${magenta}- performing SNP/variant calling ${white}\n"
		dirsnpcall=./snpcall/*
		if [ ${#dirsnpcall[@]} -gt 1 ]; then
			echo -e "${magenta}- \n- snpcall folder should be empty or contain only 1 directory, exiting pipeline ${white}\n"
			sleep 10 && exit 1
		fi
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

mkdir 2x
mkdir 4x
mkdir 6x
mkdir 8x
cd ../snpcall
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
			$GATK LeftAlignAndTrimVariants -R ../refgenomes/panref.fasta -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
			wait
			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' | awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%100000==2{x=file"SPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' -
			wait
			rm ${i%.vcf}0.vcf*
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
			$GATK LeftAlignAndTrimVariants -R ../refgenomes/panref.fasta -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
			wait
			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' | awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%100000==2{x=file"SPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' -
			wait
			rm ${i%.vcf}0.vcf*
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
			$GATK LeftAlignAndTrimVariants -R ../refgenomes/panref.fasta -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
			wait
			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' | awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%100000==2{x=file"SPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' -
			wait
			rm ${i%.vcf}0.vcf*
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
			$GATK LeftAlignAndTrimVariants -R ../refgenomes/panref.fasta -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles --keep-original-ac
			wait
			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' | awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%100000==2{x=file"SPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' -
			wait
			rm ${i%.vcf}0.vcf*
		done
		Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 8x
	fi
fi
wait
filetest=*x.vcf*
if [ -z "$(ls -A *x.vcf* 2>/dev/null)" ]; then
	for v in *_DP_GT.txt; do (
		vcfdose=${v%_DP*}; vcfdose=${vcfdose#*_}; out=$(ls *${vcfdose}_raw.vcf)
		for raw in $out; do
			grep '^#' $raw  > ${raw%_raw.vcf}.vcf
			awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$1,$2]){print b}}' $raw $v >> ${raw%_raw.vcf}.vcf
		done )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait
fi

#retrieve the original SNP positions from conigs/scaffolds before concatenation
cd $projdir
cd snpcall

reffilenames=($ref1 $ref2 $ref3)
for refg in "${reffilenames[@]}"; do
	export ncontigscaffold=$(grep '>' ../refgenomes/${refg%.fasta}_original.fasta | wc -l)
	if [[ ! -f ./index_code.txt ]]; then
		if [[ $ncontigscaffold -gt 200 ]]; then
			mkdir snps_stitched_ref
			gzip *_raw.vcf*; mv *_raw.vcf* ./snps_stitched_ref/
			echo -e "${magenta}- retrieving SNP positions based on contigs/scaffold annotation ${white}\n".
			for i in $(ls *.vcf *_DP_GT.txt); do
				if [[ "${i##*.}" = "vcf" ]]; then
					awk '/^#CHROM/{close("file.vcf"f);f++}{print $0 > "file"f}' $i
					awk 'NR==1{print}' file1 | cat file - > file0.txt; rm file
					awk 'NR>1{print $1,"\t",$2}' file1 > file_snp_stitch.txt
					awk 'NR>1{print $0}' file1 > file1.txt; rm file1
				fi
				if [[ "${i##*.}" = "txt" ]]; then
					awk 'NR==1{print}' $i > file0.txt
					awk 'NR>1{print $1,"\t",$2}' $i > file_snp_stitch.txt
					awk 'NR>1{print $0}' $i > file1.txt
				fi
				> index_code.txt
				while IFS="" read -r p || [ -n "$p" ]; do
					fchr=$(echo $p | awk '{print $1}'); fpos=$(echo $p | awk '{print $2}')
					scf=$(awk -v fchr=$fchr '$1 == fchr' ../refgenomes/contigscaffold_index.txt | \
					awk -v t=$fpos 'BEGIN {var=t; highest=0}{ j = $NF;if ( j <= var && j > highest ) { highest=j} } END {print highest}' )
					awk -v fchr=$fchr '$1 == fchr' ../refgenomes/contigscaffold_index.txt | awk -v scf=$scf -v fchr=$fchr -v fpos=$fpos '$3 == scf {print fchr"\t"fpos"\t"$0}' >> index_code.txt
				done < file_snp_stitch.txt
				awk 'BEGIN {FS=OFS="\t"}{ $5 = $2 - $5 } 1' index_code.txt > file2.txt
				awk 'BEGIN {FS=OFS="\t"}!($3="")' file2.txt | awk '{gsub(/\t\t/,"\t"); print}' | awk '{gsub(/ /,""); print}' > index_code.txt
				paste -d "\t" index_code.txt file1.txt | awk -F"\t" 'BEGIN {FS=OFS="\t"}{$1=$2=$5=$6=""; print $0}' | tr -s '\t' | sed -e 's/^[ \t]*//' > file2.txt
				cat file0.txt file2.txt > $i
				rm file0.txt file1.txt file2.txt file_snp_stitch.txt
			done

			cd ../alignment_summaries
			awk 'NR==1{print}' refgenome_paralogs.txt > file0.txt
			awk 'NR>1{print $1,"\t",$2}' refgenome_paralogs.txt > file_snp_stitch.txt
			awk 'NR>1{print $0}' refgenome_paralogs.txt > file1.txt
			> index_code.txt
			while IFS="" read -r p || [ -n "$p" ]; do
			fchr=$(echo $p | awk '{print $1}'); fpos=$(echo $p | awk '{print $2}')
			scf=$(awk -v fchr=$fchr '$1 == fchr' ../refgenomes/contigscaffold_index.txt | \
			awk -v t=$fpos 'BEGIN {var=t; highest=0}{ j = $NF;if ( j <= var && j > highest ) { highest=j} } END {print highest}' )
			awk -v fchr=$fchr '$1 == fchr' ../refgenomes/contigscaffold_index.txt | awk -v scf=$scf -v fchr=$fchr -v fpos=$fpos '$3 == scf {print fchr"\t"fpos"\t"$0}' >> index_code.txt
			done < file_snp_stitch.txt
			awk 'BEGIN {FS=OFS="\t"}{ $5 = $2 - $5 } 1' index_code.txt > file2.txt
			awk 'BEGIN {FS=OFS="\t"}!($3="")' file2.txt | awk '{gsub(/\t\t/,"\t"); print}' | awk '{gsub(/ /,""); print}' > index_code.txt
			paste -d "\t" index_code.txt file1.txt | awk -F"\t" 'BEGIN {FS=OFS="\t"}{$1=$2=$5=$6=""; print $0}' | tr -s '\t' | sed -e 's/^[ \t]*//' > file2.txt
			cat file0.txt file2.txt > refgenome_paralogs.txt
			rm file0.txt file1.txt file2.txt file_snp_stitch.txt index_code.txt
			cd ../snpcall
		fi
	fi
done

if [ "$(ls -A *.vcf 2>/dev/null)" ]; then
	for v in *.vcf; do gzip $v; done
fi
wait


cd $projdir
cd samples
window1=$(ls -S | head -1 | xargs zcat -fq | awk '{ print length }' | sort -n | tail -1)
window=$((window1 + 20))
cd ../

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
	if [ -d "${projdir}/snpfilter/2x" ]; then
		cd ${projdir}/snpfilter
		cp -r 2x 2x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./2x_diversity_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_2x.R $pop $gmiss $smiss $minRD_2x $exclude_samples "${GBSapp_dir}/tools/R" $maf $snpformats
		wait
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
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_4x.R $pop $gmiss $smiss $minRD_4x $exclude_samples "${GBSapp_dir}/tools/R" $maf $snpformats
		wait
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
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_6x.R $pop $gmiss $smiss $minRD_6x $exclude_samples "${GBSapp_dir}/tools/R" $maf $snpformats
		wait
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
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_8x.R $pop $gmiss $smiss $minRD_8x $exclude_samples "${GBSapp_dir}/tools/R" $maf $snpformats
		wait
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
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_2x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_2x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$snpformats" "$pseg"
		wait
		awk 'NR>1{print $1,"\t",$2,"\t",$3}' *dose.txt | sort -V -k2,2 -k3,3 >  snplist_rd.txt
		chrid=$(LC_ALL=C; sort -n -k2,2 -S 50% snplist_rd.txt | awk '{print $2}' | uniq)
		for i in $chrid; do
			awk -v n="$i" '$0~n{print $0}' snplist_rd.txt | awk 'NR==0{old = $3; next} {print $1,"\t",$2,"\t",$3,"\t",$3 - old; old = $3}' |\
			awk -v f=$window1 '{print $1,"\t",$2,"\t",($4>f?$3:"")}' | awk 'NF==2{print $1,"\t",$2,"\t",p "\t" $3; next} {p=$3} 1' |\
			awk -v s=$window '{print $1"\t"$2":"$3-s"-"$3+s}' >> snplist.txt
		done
		sort -V -u -k2,2 snplist.txt > snplist_nonredun.txt; rm snplist_rd.txt
		find . -type f -empty -delete
		wait
	fi
	if [ -d "${projdir}/snpfilter/4x" ]; then
		cd ${projdir}/snpfilter
		cp -r 4x 4x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./4x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_4x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_4x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$snpformats" "$pseg"
		wait
		awk 'NR>1{print $1,"\t",$2,"\t",$3}' *dose.txt | sort -V -k2,2 -k3,3 >  snplist_rd.txt
		chrid=$(LC_ALL=C; sort -n -k2,2 -S 50% snplist_rd.txt | awk '{print $2}' | uniq)
		for i in $chrid; do
			awk -v n="$i" '$0~n{print $0}' snplist_rd.txt | awk 'NR==0{old = $3; next} {print $1,"\t",$2,"\t",$3,"\t",$3 - old; old = $3}' |\
			awk -v f=$window1 '{print $1,"\t",$2,"\t",($4>f?$3:"")}' | awk 'NF==2{print $1,"\t",$2,"\t",p "\t" $3; next} {p=$3} 1' |\
			awk -v s=$window '{print $1"\t"$2":"$3-s"-"$3+s}' >> snplist.txt
		done
		sort -V -u -k2,2 snplist.txt > snplist_nonredun.txt; rm snplist_rd.txt
		find . -type f -empty -delete
		wait
	fi
	if [ -d "${projdir}/snpfilter/6x" ]; then
		cd ${projdir}/snpfilter
		cp -r 6x 6x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./6x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_6x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_6x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$snpformats" "$pseg"
		wait
		awk 'NR>1{print $1,"\t",$2,"\t",$3}' *dose.txt | sort -V -k2,2 -k3,3 >  snplist_rd.txt
		chrid=$(LC_ALL=C; sort -n -k2,2 -S 50% snplist_rd.txt | awk '{print $2}' | uniq)
		for i in $chrid; do
			awk -v n="$i" '$0~n{print $0}' snplist_rd.txt | awk 'NR==0{old = $3; next} {print $1,"\t",$2,"\t",$3,"\t",$3 - old; old = $3}' |\
			awk -v f=$window1 '{print $1,"\t",$2,"\t",($4>f?$3:"")}' | awk 'NF==2{print $1,"\t",$2,"\t",p "\t" $3; next} {p=$3} 1' |\
			awk -v s=$window '{print $1"\t"$2":"$3-s"-"$3+s}' >> snplist.txt
		done
		sort -V -u -k2,2 snplist.txt > snplist_nonredun.txt; rm snplist_rd.txt
		find . -type f -empty -delete
		wait
	fi
	if [ -d "${projdir}/snpfilter/8x" ]; then
		cd ${projdir}/snpfilter
		cp -r 8x 8x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		cd ./8x_biparental_gmiss"${gmiss}"_smiss"${smiss}"
		Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_8x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_8x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$snpformats" "$pseg"
		wait
		awk 'NR>1{print $1,"\t",$2,"\t",$3}' *dose.txt | sort -V -k2,2 -k3,3 >  snplist_rd.txt
		chrid=$(LC_ALL=C; sort -n -k2,2 -S 50% snplist_rd.txt | awk '{print $2}' | uniq)
		for i in $chrid; do
			awk -v n="$i" '$0~n{print $0}' snplist_rd.txt | awk 'NR==0{old = $3; next} {print $1,"\t",$2,"\t",$3,"\t",$3 - old; old = $3}' |\
			awk -v f=$window1 '{print $1,"\t",$2,"\t",($4>f?$3:"")}' | awk 'NF==2{print $1,"\t",$2,"\t",p "\t" $3; next} {p=$3} 1' |\
			awk -v s=$window '{print $1"\t"$2":"$3-s"-"$3+s}' >> snplist.txt
		done
		sort -V -u -k2,2 snplist.txt > snplist_nonredun.txt; rm snplist_rd.txt
		find . -type f -empty -delete
		wait
	fi
fi
done
done


cd ${projdir}/snpfilter/
find . -type f -empty -delete
find . -type d -empty -delete
wc -l *gmiss*/*dose.txt | awk '{print $2"\t"$1-1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > gmiss_smiss_titration.txt
wc -l *gmiss*/eliminated* | awk '{print $2"\t"$1-1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > eliminated_samples.txt
echo -e "gmiss_smiss_thresholds\t#_of_SNPs\t#_of_eliminated_samples\n----------------------\t---------\t-----------------------" > summary_header.txt
awk 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,"\t",a[$1],"\t",$2}' gmiss_smiss_titration.txt eliminated_samples.txt | cat summary_header.txt - > gmiss_smiss_titration_eliminated_samples.txt
rm gmiss_smiss_titration.txt eliminated_samples.txt summary_header.txt


echo -e "${magenta}- Anchoring SNP/variants to a reference subgenome: ${ref1%.fa*} ${white}\n"
cd $projdir
cd refgenomes
file=${ref1%.f*}.dict
if test -f "$file"; then
	echo -e "${magenta}- indexed subgenome already exist ${white}\n"
else
	$bwa index -a bwtsw ${ref1}
	$samtools faidx ${ref1}
	$java $Xmxg -XX:ParallelGCThreads=$gthreads -jar $picard CreateSequenceDictionary REFERENCE= ${ref1} OUTPUT=${ref1%.f*}.dict
fi

cd $projdir
cd snpfilter
if [ "$ploidy" == 6 ]; then
	for snpfilter_dir in $(ls -d 6x*/); do
		if [ -d "$snpfilter_dir" ]; then
			cd $snpfilter_dir
			touch  refpos.fasta
			awk 'NR >= 2 {print "$samtools faidx ${projdir}/refgenomes/panref.fasta "$2":"$3"-"($3 + 1000)" >> refpos.fasta;"}' *dose.txt > snpseq_context.sh
			bash snpseq_context.sh
			wait
			$bwa mem -t $loopthreads ${projdir}/refgenomes/${ref1} refpos.fasta > snpseq_context.sam
			grep -v '^@' snpseq_context.sam | awk -F'\t' '{print $1"\t"$3"_"$4"\t"$3"\t"$4}' | awk '{gsub(/-/,"\t")}1' | awk '{gsub(/:/,"_")}1' | awk '{gsub(/__/,"_")}1' > SNP_ID_ref1.txt
			echo -e "SNP\textend\tSNP\tCHROM\tPOS\n$(cat SNP_ID_ref1.txt)" > SNP_ID_ref1.txt
			for i in *dose.txt *binary.txt; do
				awk -F"\t" 'FNR==NR{a[$1]=$3"\t"$4"\t"$5;next} ($1 in a) {print $1,a[$1],$0}' OFS="\t" SNP_ID_ref1.txt  $i | \
				awk '{$1=$5=$6=$7=""; print}' OFS="\t" | tr -s '\t' | awk '{$1=$1;print}' OFS="\t" > ${i%.txt*}_${ref1%.fa*}_anchored.txt
			done
			if [[ "${snpformats}" == "true" ]]; then
			for i in *nucleotide.txt; do
				awk -F"\t" 'FNR==NR{a[$1]=$3"\t"$4"\t"$5;next} ($1 in a) {print $1,a[$1],$0}' OFS="\t" SNP_ID_ref1.txt  $i | \
				awk '{$1=$5=$6=$7=""; print}' OFS="\t" | tr -s '\t' | awk '{$1=$1;print}' OFS="\t" > ${i%.txt*}_${ref1%.fa*}_anchored.txt
			done
			fi
		rm refpos* snpseq_context.* SNP_ID_ref1.txt
		fi
		wait
		cd ../
	done
fi
if [ "$ploidy" == 8 ]; then
	for snpfilter_dir in $(ls -d 8x*/); do
		if [ -d "$snpfilter_dir" ]; then
			cd $snpfilter_dir
			touch  refpos.fasta
			awk 'NR >= 2 {print "$samtools faidx ${projdir}/refgenomes/panref.fasta "$2":"$3"-"($3 + 1000)" >> refpos.fasta;"}' *dose.txt > snpseq_context.sh
			bash snpseq_context.sh
			wait
			$bwa mem -t $loopthreads ${projdir}/refgenomes/${ref1} refpos.fasta > snpseq_context.sam
			grep -v '^@' snpseq_context.sam | awk -F'\t' '{print $1"\t"$3"_"$4"\t"$3"\t"$4}' | awk '{gsub(/-/,"\t")}1' | awk '{gsub(/:/,"_")}1' | awk '{gsub(/__/,"_")}1' > SNP_ID_ref1.txt
			echo -e "SNP\textend\tSNP\tCHROM\tPOS\n$(cat SNP_ID_ref1.txt)" > SNP_ID_ref1.txt
			for i in *dose.txt *binary.txt; do
				awk -F"\t" 'FNR==NR{a[$1]=$3"\t"$4"\t"$5;next} ($1 in a) {print $1,a[$1],$0}' OFS="\t" SNP_ID_ref1.txt  $i | \
				awk '{$1=$5=$6=$7=""; print}' OFS="\t" | tr -s '\t' | awk '{$1=$1;print}' OFS="\t" > ${i%.txt*}_${ref1%.fa*}_anchored.txt
			done
			if [[ "${snpformats}" == "true" ]]; then
			for i in *nucleotide.txt; do
				awk -F"\t" 'FNR==NR{a[$1]=$3"\t"$4"\t"$5;next} ($1 in a) {print $1,a[$1],$0}' OFS="\t" SNP_ID_ref1.txt  $i | \
				awk '{$1=$5=$6=$7=""; print}' OFS="\t" | tr -s '\t' | awk '{$1=$1;print}' OFS="\t" > ${i%.txt*}_${ref1%.fa*}_anchored.txt
			done
			fi
		rm refpos* snpseq_context.* SNP_ID_ref1.txt
		fi
		wait
		cd ../
	done
fi
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


######################################################################################################################################################
wait
echo -e "${magenta}- Run Complete. ${white}\n"
