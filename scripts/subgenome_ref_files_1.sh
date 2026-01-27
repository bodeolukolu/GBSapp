
if [ -z "${slurm_module:-}" ]; then
 export slurm_module=true
fi
if [ -z "${threads:-}" ]; then
	export threads=$(nproc --all)
	if [[ "$threads" -ge 4 ]]; then
		export threads=$((threads-2))
	fi
fi
export aligner=minimap2

if [ -z "${RNA:-}" ]; then
 export RNA=false
fi
if [ -z "${variant_caller:-}" ]; then
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
fi
if [ -z "${lib_type:-}" ]; then
 export lib_type=RRS
fi
if [ -z "${subsample_WGS_in_silico_qRRS:-}" ]; then
 export subsample_WGS_in_silico_qRRS=false
fi
if [ -z "${nodes:-}" ]; then
 export nodes=1
fi
if [ -z "${biallelic:-}" ]; then
	export biallelic=false
fi
if [ -z "${paralogs:-}" ]; then
	export paralogs=false
fi
if [ -z "${max_pseudoMol:-}" ]; then
	export max_pseudoMol=5000
fi
if [ -z "${uniquely_mapped:-}" ]; then
	export uniquely_mapped=true
fi
if [ -z "${minmapq:-}" ]; then
	export minmapq=20
fi
if [ "$RNA" == "true" ]; then
	export minmapq=1
fi
if [ -z "${maxHaplotype:-}" ]; then
	export maxHaplotype=128
fi
if [ -z "${haplome_number:-}" ]; then
	export haplome_number=1
fi
if [ -z "${p2:-}" ]; then
  if   [ "${p1:-}" ]; then
	 export p2=$p1
 fi
fi
if [ -z "${use_softclip:-}" ]; then
	export use_softclip=false
fi
if [ "$use_softclip" == "false" ]; then
	export dont_use_softclip=true
fi
if [ "$use_softclip" == "true" ]; then
	export dont_use_softclip=false
fi

if [ -z "${joint_calling:-}" ]; then
	export joint_calling=false
fi
if [ -z "${keep_gVCF:-}" ]; then
	export keep_gVCF=false
fi
if [ -z "${filter_ExcHet:-}" ]; then
  filter_ExcHet=false
fi
if [ -z "${filtered_vcf:-}" ]; then
  filtered_vcf=true
fi
if [ -z "${genomecov_est:-}" ]; then
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

  if [[ "$threads" -le 8 ]]; then
    export alnthreads=$threads
    export alnN=1
  else
    export alnthreads=8
    export alnN=$(( threads / alnthreads ))
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
  shopt -s nullglob
  for i in *.gz; do
      sleep $((RANDOM % 2))
      gunzip "$i" >/dev/null 2>&1 || true
  done
  shopt -u nullglob

  cd $projdir
  cd refgenomes
  if [[ -f "${ref1%.f*}_unstitched.fasta" ]]; then
    mv "${ref1%.f*}_unstitched.fasta" ../
  fi
  if [[ -f "unsplit_${ref1%.f*}_fasta.txt" ]]; then
    mv "unsplit_${ref1%.f*}_fasta.txt" ../
  fi
  originalREF=$(ls *_original.fasta 2>/dev/null)
  if [[ -f "$originalREF" ]]; then
    mv "${ref1%.f*}_original.fasta" ../$ref1
    if [[ -d pangenomes ]]; then mv pangenomes ../; fi
    rm -rf ./*
    mv "../$ref1" ./
    mv "../pangenomes" ./
  fi
  if [[ -f "../${ref1%.f*}_unstitched.fasta"  ]]; then
    mv "../${ref1%.f*}_unstitched.fasta" ./
  fi
  if [[ -f "../unsplit_${ref1%.f*}_fasta.txt" ]]; then
    mv "../unsplit_${ref1%.f*}_fasta.txt" ./
  fi

  if ls ./*.dict 1> /dev/null 2>&1; then
  	:
  else
  	checknfastafiles=$(ls *.f* | grep -v .fai | grep -v .ngm | grep -v _original.fasta | wc -l)
  	if [[ $checknfastafiles -gt 1 ]]; then
  		echo -e "${magenta}- expecting only 1 fasta file for reference genome ${white}\n"
  		echo -e "${magenta}- GBSapp will quit in 5 seconds ${white}\n"
  		sleep 5 && exit 1
  	fi
  	if [ -z "${ref1:-}" ]; then
  		for ref in *.f*; do
  			sleep $((RANDOM % 2))
        ref1=${ref%.fa*}.fasta
  		done
  	fi

    # Filter contigs >= 1000 bp (robust)
    if awk '/^>/{if(NR>1 && len<1000) exit 1; len=0; next}{len+=length($0)}END{if(len<1000) exit 1}' "$ref1"; then
      echo "All contigs/scaffolds/chromosomes sequences >= 1000 bp"
    else
      awk 'BEGIN{RS=">"; ORS=""}
      NR>1{
        n = split($0, a, "\n")
        seq = ""
        for (i = 2; i <= n; i++) seq = seq a[i]
        if (length(seq) >= 1000) print ">" $0
      }' "$ref1" > "${ref1}.tmp" &&
      mv "${ref1}.tmp" "$ref1"
    fi
    ncontigscaffold=$(grep -c '^>' "$ref1")
    if [[ $ncontigscaffold -gt $max_pseudoMol ]]; then
      nfakechr=$(( threads > 1 ? threads/2 : 1 ))
      # Convert FASTA to single-line-per-contig (tab-separated)
      awk '/^>/{if(N>0) printf("\n"); ++N; printf("%s\t",$0); next}
           {printf("%s",$0)}
           END{printf("\n")}' "$ref1" > panref0.txt
      # Shuffle contigs (deterministic order not required here)
      awk '{ a[NR]=$0 }
      END{
        srand()
        for(i=NR;i>0;i--){
          j=int(rand()*i)+1
          tmp=a[i]; a[i]=a[j]; a[j]=tmp
        }
        for(i=1;i<=NR;i++) print a[i]
      }' panref0.txt > panref0.fasta
      # Split into pseudochromosomes
      flength=$(wc -l panref0.fasta | awk '{print $1}')
      nsplit=$(( (flength + nfakechr - 1) / nfakechr ))
      split -a 2 -d -l "$nsplit" panref0.fasta Chr
      rm -f panref0.txt panref0.fasta
      for i in Chr*; do
        sleep $((RANDOM % 2))
        # Reconstruct FASTA
        awk '{print $1"\n"$2}' "$i" > "${i}.fasta"
        # Wrap sequence at 100 bp
        awk '/^>/{print; next}
             {gsub(/\s+/,""); print}' "${i}.fasta" | fold -w 100 > "${i}.txt"
        # Generate contig → pseudochromosome index
        awk '/^>/{if (l!="") print l; print; l=0; next}
        {l+=length($0)}
        END{print l}' "${i}.fasta" |
        awk '{gsub(/>/,""); print}' |
        awk '{ ORS = NR % 2 ? "\t" : "\n" } 1' |
        awk '{print $1"\t"$2+100}' |
        awk -F "\t" '{sum+=$2; print $1"\t"(sum-$2)}' |
        awk -v chrid="$i" -v refgenome="${ref1%.f*}" \
            '{print refgenome"_"chrid"\t"$1"\t"$2}' >> contigscaffold_index.txt
      done

      # Preserve filtered-but-unstitched reference
      cp "$ref1" "${ref1%.f*}_unstitched.fasta"
      # Build stitched pseudochromosome reference
      : > "$ref1"
      for filename in Chr*.txt; do
        sleep $((RANDOM % 2))
        echo ">${filename%.txt}" >> "$ref1"
        cat "$filename" >> "$ref1"
      done
      # Cleanup (safe)
      rm -f Chr?? Chr??.txt Chr??.fasta
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
      mv $ref1 "./unsplit_${ref1%.f*}_fasta.txt"
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
  if [[ "$aligner" == "minimap2" ]]; then
    if [[ -n "$(compgen -G "./*.mmi")" ]]; then
    	echo -e "${magenta}- indexed genome available ${white}\n"
    	if [ -z "${ref1:-}" ]; then
    		for ref in *.f*; do
    			sleep $((RANDOM % 2))
          ref1=${ref%%.f*}.fasta
    		done
    	fi
    else
    	echo -e "${magenta}- indexing single reference subgenome ${white}\n"
    	if [ -z "${ref1:-}" ]; then
    		for ref in *.f*; do
    			sleep $((RANDOM % 2))
          ref1=${ref%%.f*}.fasta
    		done
    	fi
    	awk '{ sub("\r$",""); print}' $ref1 | awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' | \
      awk '/>/{gsub(/a$/,"A");gsub(/b$/,"B");gsub(/c$/,"C");gsub(/d$/,"D");gsub(/e$/,"E");gsub(/f$/,"F");gsub(/g$/,"G");gsub(/h$/,"H");}1' > ref.txt &&
    	n=">${ref1%.f*}_" &&
    	awk '{ sub("\r$",""); print}' ref.txt | awk -v n="$n" '{gsub(n,">"); print}' | awk -v n="$n" '{gsub(/>/,n); print}' > $ref1 &&
    	rm -f ref.txt

      # processing pangenomes
      if [[ -d pangenomes ]]; then
        cp "$ref1" "${ref1%.f*}_original.fasta" &&
        $samtools faidx "${ref1%.f*}_original.fasta" &&
        $java -jar $picard CreateSequenceDictionary REFERENCE="${ref1%.f*}_original.fasta" OUTPUT="${ref1%.f*}_original.dict" &&
        $minimap2 -d "${ref1%.f*}_original.mmi" "${ref1%.f*}_original.fasta" &&
        wait

        cd pangenomes
        shopt -s nullglob
        for panref in *.fa.gz *.fna.gz *.fasta.gz; do
          gunzip -f "$panref"
        done
        for panref in *.fa *.fna *.fasta; do
          mv "$panref" "${panref%.f*}.fasta"
        done

        # remove small contigs (<1kb) and limit pseudomolecules to 5,000 sequences
        shopt -s nullglob
        for panref in *.fasta; do (
          ## 1. Filter sequences ≥1000 bp (sequence-only length)
          if awk '/^>/{if(NR>1 && len<1000) exit 1; len=0; next}{len+=length($0)}END{if(len<1000) exit 1}' "$panref"; then
            echo "All contigs/scaffolds/chromosomes sequences >= 1000 bp"
          else
            awk 'BEGIN{RS=">"; ORS=""}
            NR>1{
              n = split($0, a, "\n")
              seq = ""
              for (i = 2; i <= n; i++) seq = seq a[i]
              if (length(seq) >= 1000) print ">" $0
            }' "$panref" > "${panref}.tmp" &&
            mv "${panref}.tmp" "$panref"
          fi
          ncontigscaffold=$(grep '>' "$panref" | wc -l)
          if [[ "$ncontigscaffold" -gt "$max_pseudoMol" ]]; then
            mkdir original_pangenomes
            nfakechr=$((threads/2))
            awk 'BEGIN{RS=">"; FS="\n"}
            NR>1{
              header=$1
              seq=""
              for(i=2;i<=NF;i++) seq=seq $i
              print ">"header"\t"seq
            }' "$panref" > "${panref%.f*}_1.fasta"
            awk '{ a[NR]=$0 }
            END{
              srand()
              for(i=NR;i>0;i--){
                j=int(rand()*i)+1
                tmp=a[i]; a[i]=a[j]; a[j]=tmp
              }
              for(i=1;i<=NR;i++) print a[i]
            }' "${panref%.f*}_1.fasta" > "${panref%.f*}_2.fasta"
            flength=$(wc -l "${panref%.f*}_2.fasta" | awk '{print $1}')
            nsplit=$(( (flength + nfakechr - 1) / nfakechr ))
            ## 5. Split into pseudochromosomes
            split -a 2 -d -l "$nsplit" "${panref%.f*}_2.fasta" Chr
            rm -f "${panref%.f*}_1.fasta" "${panref%.f*}_2.fasta"
            ## index file per panref (FIXED: no race condition)
            index_tmp="contigscaffold_index_${panref%.f*}.txt"
            : > "$index_tmp"

            for i in Chr*; do
              sleep $((RANDOM % 2))
              ## 6. Restore FASTA formatting
              awk '{print $1"\n"$2}' "$i" > "${i}.fasta"
              grep -v '^>' "${i}.fasta" | tr -d '\n' | fold -w 100 > "${i}.txt"
              ## 7. Build coordinate index (APPEND, not overwrite)
              awk '/^>/{if (l!="") print l; print; l=0; next}
                   {l+=length($0)}
                   END{print l}' "${i}.fasta" | awk '{gsub(/>/,""); print $0}' | awk '{ ORS = NR % 2 ? "\t" : "\n" } 1' | \
              awk '{print $1"\t"$2+100}' | awk -F "\t" '{sum+=$2; print $1"\t"sum-$2}' | \
              awk -v chrid="$i" -v refgenome="${panref%.f*}" '{print refgenome"_"chrid"\t"$1"\t"$2}' >> "$index_tmp"
            done
            cp "$panref" original_pangenomes/
            ## 9. Replace FASTA with stitched pseudochromosomes
            : > "$panref"
            for filename in Chr*.txt; do
              sleep $((RANDOM % 2))
              echo ">${filename%.txt}" >> "$panref"
              cat "$filename" >> "$panref"
            done
            rm -f Chr* "${index_tmp}"
          fi
        ) &
        while (( $(jobs -rp | wc -l) >= $gN )); do sleep 2; done
        done
        wait
        ## 10. Merge all index files SAFELY
        if [[ "$ncontigscaffold" -gt "$max_pseudoMol" ]]; then
          cat contigscaffold_index_*.txt > contigscaffold_index.txt
          rm -f contigscaffold_index_*.txt
        fi

        # add pangenome as prefix to pangenomes fasta header
        for panref in *.fasta; do (
          awk '{ sub("\r$",""); print}' "$panref" | awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' | \
          awk '/>/{gsub(/a$/,"A");gsub(/b$/,"B");gsub(/c$/,"C");gsub(/d$/,"D");gsub(/e$/,"E");gsub(/f$/,"F");gsub(/g$/,"G");gsub(/h$/,"H");}1' > "${panref%%.fasta}.txt" &&
          n=">pangenome_${panref%.fasta}_" &&
          awk '{ sub("\r$",""); print}' "${panref%%.fasta}.txt" | awk -v n="$n" '{gsub(n,">"); print}' | awk -v n="$n" '{gsub(/>/,n); print}' > "${panref}.tmp" &&
          mv "${panref}.tmp" "$panref" &&
          rm -f "${panref%%.fasta}.txt"
          ) &
          while (( $(jobs -rp | wc -l) >= $gN )); do sleep 2; done
        done
        wait

        cd ../
        $samtools faidx $ref1 &&
        $java -jar $picard CreateSequenceDictionary REFERENCE=$ref1 OUTPUT=${ref1%.f*}.dict &&
        $minimap2 -d ${ref1%.f*}.mmi $ref1 &&

        cd pangenomes
        : > panref.raw.fa
        # 1. Extract divergent regions
        for panref in *.fasta; do (
          $minimap2 -x asm5 -t "$gthreads" "../$ref1" "$panref" > "${panref%.fasta}.paf" &&
          awk '{print $6,$8,$9}' OFS='\t' "${panref%.fasta}.paf" > "${panref%.fasta}.aligned.bed" &&
          $samtools faidx "$panref" &&
          awk -v OFS='\t' '{print $1,0,$2}' "${panref}.fai" > "${panref%.fasta}.full.bed" &&
          $bedtools subtract -a "${panref%.fasta}.full.bed" -b "${panref%.fasta}.aligned.bed" \
          | $bedtools getfasta -fi "$panref" -bed - -fo "${panref%.fasta}.div.fa"
        ) &
        while (( $(jobs -rp | wc -l) >= $gN )); do sleep 1; done
        done
        wait
        cat *.div.fa > panref.raw.fa
        rm -f *.paf *.bed *.fai *.div.fa
        # 2. Self-alignment (for near-identical collapse)
        $minimap2 -x asm20 -c --cs -t "$threads" panref.raw.fa panref.raw.fa > panref.self.paf
        # 3. Identify redundant sequences (≥95% coverage, not self)
        awk '$1 != $6 {
          qlen=$2; tlen=$7
          aln=$10
          cov=aln/(qlen<tlen?qlen:tlen)
          if(cov >= 0.95){
            print $1
          }
        }' panref.self.paf | sort -u > redundant.ids
        # 4. Keep only non-redundant (longest representative survives)
        awk -v dropfile="redundant.ids" 'BEGIN {
          RS=">"; FS="\n"; ORS=""
          while ((getline < dropfile) > 0) {
            drop[$1] = 1
          }
          close(dropfile)
        }
        NR > 1 {
          id = $1
          if (!(id in drop)) {
            print ">" $0
          }
        }' panref.raw.fa | \
        awk '/^>/ { sub(/:.*/, "", $0); print; next } { print }' > panref.fasta
        rm -f panref.self.paf redundant.ids panref.raw.fa

        $minimap2 -x asm5 -t $threads "../${ref1%.f*}_original.fasta" panref.fasta > ../pangenome2primary.paf
        mv panref.fasta ../
        cat "../${ref1}" ../panref.fasta > ref_combined.fasta &&
        mv ref_combined.fasta "../${ref1}" && wait
      fi

      cd $projdir
      cd refgenomes
      if [[ -d pangenomes ]]; then
        $samtools faidx panref.fasta &&
        $java -jar $picard CreateSequenceDictionary REFERENCE=panref.fasta OUTPUT=panref.dict &&
        $minimap2 -d panref.mmi panref.fasta
        rm -f ${ref1%.f*}.dict ${ref1%.f*}.fai ${ref1%.f*}.mmi
      fi

      $samtools faidx $ref1 &&
      $java -jar $picard CreateSequenceDictionary REFERENCE=$ref1 OUTPUT=${ref1%.f*}.dict &&
      $minimap2 -d ${ref1%.f*}.mmi $ref1 &&
      $genmap index -F $ref1 -I ./genmap_out >/dev/null 2>&1
      wait
      $genmap map -K 100 -E 2 -I ./genmap_out -O ./genmap_out/ -t -w -bg --threads $threads >/dev/null 2>&1
      wait
      if [[ "$paralogs" == "false" ]]; then
        threshold=1
      else
        threshold=8
      fi
      wait
      ref_base=$(basename "$ref1" .fasta) &&
      python3 "$wig2bed" "${projdir}/refgenomes/genmap_out/${ref_base}.genmap.wig" "$threshold" "${projdir}/refgenomes/genmap_out/lowmap_merged.bed" > /dev/null 2>&1
      wait
      $bedtools maskfasta -fi "${projdir}/refgenomes/$ref1" -bed "${projdir}/refgenomes/genmap_out/lowmap_merged.bed" -fo "${projdir}/refgenomes/${ref1%.f*}.hardmasked.fasta" >/dev/null 2>&1
      wait
      mv "${ref1%.f*}.hardmasked.fasta" "$ref1" &&
      rm -f ${ref1%.f*}.dict ${ref1%.f*}.fai ${ref1%.f*}.mmi
      $samtools faidx "$ref1" &&
      $java -jar $picard CreateSequenceDictionary REFERENCE=$ref1 OUTPUT=${ref1%.f*}.dict &&
      $minimap2 -d "${ref1%.f*}.mmi" "$ref1"
      wait
    fi
  fi
  wait
  if [[ "$RNA" == "true" ]]; then
    mkdir -p star_index
    if [[ -z "$(compgen -G "./star_index/*")" ]]; then
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
if [[ "$alignments" == 1 ]] && [[ "$snp_calling" == 1 ]]; then
  if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
    if [[ ! -f compress_done.txt && ! -f organize_files_done.txt ]]; then
      time main 2>> ${projdir}/log.out
    fi
  fi
fi

set -Eeuo pipefail
trap 'echo "FAILED at line $LINENO"' ERR
echo -e "${blue}\n############################################################################## ${yellow}\n- Organizing sample fastq files \n${blue}##############################################################################${white}\n"
main () {
	cd $projdir
	cd samples
  if [[ -z "$(compgen -G "*.f*")" ]] || [[ -n "$(compgen -G "./preprocess/alignment/*sam*")" ]]; then
    if [[ ! -d "pe" ]] && [[ ! -d "se" ]]; then
      :> filename_reformatted.txt
      :> flushed_reads.txt
    fi
  fi

  cd "${projdir}/samples"

  cd "${projdir}/samples" || exit 1

  if [[ ! -f filename_reformatted.txt ]]; then

    normalize_ext() {
      local f="$1" core gz=""
      [[ "$f" == *.gz ]] && gz=".gz" && core="${f%.gz}" || core="$f"
      case "$core" in
        *.fa|*.fna|*.fasta) core="${core%.*}.fasta" ;;
        *.fq|*.fastq)       core="${core%.*}.fastq" ;;
      esac
      echo "${core}${gz}"
    }

    process_dir() {
      local dir="$1"
      shopt -s nullglob
      for f in "$dir"/*.*; do (
        [[ -f "$f" ]] || continue
        base=$(basename "$f")
        normfile=$(normalize_ext "$base")
        if [[ "$base" != "$normfile" ]]; then
          mv "$f" "$dir/$normfile"
        else
          normfile="$base"
        fi
        [[ "$dir/$normfile" != *.gz ]] && gzip -f "$dir/$normfile"
      ) &
      while (( $(jobs -rp | wc -l) >= $gN )); do sleep 2; done
      done
      shopt -u nullglob
      wait
    }

    merge_se_orphans() {
      local dir="$1"
      local enforce_r1r2="$2"   # true / false
      shopt -s nullglob
      for f in "$dir"/*.R1.*.gz "$dir"/*_R1.*.gz; do (
        [[ -f "$f" ]] || continue
        name=$(basename "$f")
        base="${name%%.R1*}"
        base="${base%%_R1*}"
        case "$name" in
          *.fastq.gz) out_ext="fastq.gz" ;;
          *.fasta.gz) out_ext="fasta.gz" ;;
          *) exit 0 ;;
        esac
        R2=$(ls "$dir/${base}.R2."*.gz "$dir/${base}_R2."*.gz 2>/dev/null | head -n1 || true)
        if [[ -n "${R2:-}" ]]; then
          # Normal orphan merge
          zcat "$f" "$R2" | gzip > "$dir/${base}_R1R2.$out_ext"
          rm -f "$f" "$R2"
        elif [[ "$enforce_r1r2" == true ]]; then
          mv "$f" "$dir/${base}_R1R2.$out_ext"
        fi
      ) &
      while (( $(jobs -rp | wc -l) >= $gN )); do sleep 2; done
      done
      shopt -u nullglob
      wait
    }

    mkdir -p se pe
    process_dir se
    process_dir pe

    shopt -s nullglob
    for dir in "$projdir"/samples/se "$projdir"/samples/pe "$projdir"/samples/; do
        [[ -d "$dir" ]] || continue
        echo "Processing directory: $dir"
        ############################
        # PASS 1: renaming R1 file with _R1 for Paired-end (R2-driven,)
        ############################
        for r2 in "$dir"/*R2.*.gz "$dir"/*_R2.*.gz; do
            [[ -f "$r2" ]] || continue
            filename=$(basename "$r2")
            base="${filename%%.R2*}"
            base="${base%%_R2*}"
            r1=$(ls "$dir/${base}.R1."*.gz "$dir/${base}_R1."*.gz 2>/dev/null | head -n1 || true)
            if [[ -z "${r1:-}" ]]; then
                echo "Warning: found R2 file $filename without matching R1!"
            else
                r1_base=$(basename "$r1")
                if [[ "$r1_base" != *".R1."* && "$r1_base" != *"_R1."* ]]; then
                    ext="${r1_base#${base}.}"
                    new_r1="$dir/${base}_R1.${ext}"
                    echo "Renaming $r1 → $new_r1"
                    mv "$r1" "$new_r1"
                fi
            fi
        done
        ############################
        # PASS 2: renaming R1 file with _R1 for Single-end R1-only
        ############################
        for r1 in "$dir"/*.gz; do
            [[ -f "$r1" ]] || continue
            r1_base=$(basename "$r1")
            # Skip files already containing R1/R2
            [[ "$r1_base" == *".R1."* || "$r1_base" == *"_R1."* ]] && continue
            [[ "$r1_base" == *".R2."* || "$r1_base" == *"_R2."* ]] && continue
            base="${r1_base%%.*}"
            ext="${r1_base#${base}.}"
            new_r1="$dir/${base}_R1.${ext}"
            echo "Renaming single-end R1 $r1 → $new_r1"
            mv "$r1" "$new_r1"
        done
    done
    shopt -s nullglob

    has_pe=false
    shopt -s nullglob
    if compgen -G "pe/*R1*.gz" > /dev/null; then
      has_pe=true
    fi
    shopt -u nullglob
    merge_se_orphans se "$has_pe"

    for d in se pe; do
      shopt -s nullglob
      files=( "$d"/* )
      shopt -u nullglob
      (( ${#files[@]} > 0 )) && mv -n "${files[@]}" .
    done

    rm -rf se pe
    touch filename_reformatted.txt
  fi

  # in silico RRS sub-sampling of WGS
  if [[ ! -f in_silico_RRS.txt ]]; then
    if [[ "$lib_type" =~ [Ww][Gg][Ss] ]]; then
      if [[ "$subsample_WGS_in_silico_qRRS" == false ]]; then
        printf "Improvement in flushed reads not required for shotgun WGS data\n"
        touch in_silico_RRS.txt
      else
        shopt -s nullglob
        files=(*.f*.gz)
        for i in "${files[@]}"; do
          (
            # Skip already processed or unwanted files
            [[ "$i" == *_uniq.fasta.gz ]] && continue
            [[ "$i" == *_R1_uniq.fasta.gz ]] && continue
            [[ "$i" == *_R2_uniq.fasta.gz ]] && continue
            [[ "$i" == *_uniq.hold.fasta.gz ]] && continue
            [[ "$i" == *_R1_uniq.hold.fasta.gz ]] && continue
            [[ "$i" == *_R2_uniq.hold.fasta.gz ]] && continue
            [[ "$i" == *_tmp.fasta ]] && continue

            echo "Processing $i for in silico RRS"

            # Determine reader command
            if file "$i" 2>/dev/null | grep -q gzip; then
              reader="zcat $i"
            else
              reader="cat $i"
            fi

            tmpfile="${i%.f*}_tmp.fasta"

            $reader | awk -v MODE="$MODE" '
            BEGIN {
              MIN=64; MAX=600
              H="AAGCTT"; M="CCGG"
              if(MODE=="medium"){
                keep["H,H"]=keep["H,M"]=keep["M,H"]=keep["M,M"]=1
                keep["H,N"]=keep["N,H"]=keep["M,N"]=keep["N,M"]=1
              } else if(MODE=="low"){
                keep["H,H"]=keep["H,M"]=keep["M,H"]=1
                keep["H,N"]=keep["N,H"]=1
              }
            }
            function emit(seq,n,i,o,tmp,p,cut,t,l,r,len,frag){
              n=1; cut[1]=0; t[1]="N"; tmp=seq; o=0
              while(match(tmp,/(AAGCTT|CCGG)/)){
                p=o+RSTART+RLENGTH-1
                cut[++n]=p
                t[n]=(substr(tmp,RSTART,RLENGTH)==H?"H":"M")
                o+=RSTART
                tmp=substr(tmp,RSTART+1)
              }
              cut[++n]=length(seq); t[n]="N"
              for(i=1;i<n;i++){
                l=t[i]; r=t[i+1]
                len=cut[i+1]-cut[i]
                if(len>=MIN && len<=MAX && keep[l","r]){
                  frag=substr(seq,cut[i]+1,len)
                  # Trim homopolymers at start/end
                  gsub(/^(A{10,}|C{10,}|G{10,}|T{10,})/,"",frag)
                  gsub(/(A{10,}|C{10,}|G{10,}|T{10,})$/,"",frag)
                  if(length(frag)>=MIN){printf(">frag_%d_%s_%s\n%s\n",NR,l,r,frag)}
                }
              }
            }
            /^@/ {next}      # skip FASTQ header lines
            /^>/ {next}      # skip FASTA headers
            {emit($0)}' > "$tmpfile" && mv "$tmpfile" "${i%.f*}.fasta"

          ) &

          while (( $(jobs -rp | wc -l) >= gN )); do sleep 2; done
        done
        wait
        echo "In silico RRS-based resample of WGS implemented" > in_silico_RRS.txt
      fi
    fi
  fi
  wait

  # flush at read end, compress and index fasta file
  shopt -s nullglob
  if [[ ! -f flushed_reads.txt ]]; then
    # Helpers
    # Detect format from first character
    seq_format() {
      zcat "$1" 2>/dev/null | head -n1 | cut -c1 || echo ""
    }
    # Convert FASTQ or FASTA (multi-line allowed) → single-line sequences
    to_singleline_seqs() {
      local file="$1"
      local fmt="$2"

      if [[ "$fmt" == "@" ]]; then
        # FASTQ → sequences only
        zcat "$file" | awk 'NR % 4 == 2'
      elif [[ "$fmt" == ">" ]]; then
        # FASTA → flatten sequences
        zcat "$file" | awk '
          /^>/ {if (seq) print seq; seq=""; next}
          {seq=seq $0}
          END {if (seq) print seq}
        '
      else
        return 1
      fi
    }
    # Parameters
    MIN_LEN=64
    SAMPLE_PERCENT=1
    : "${gN:=4}"

    # Step 1: Convert → single-line, trim homopolymers
    for i in *.f*.gz; do
      [[ "$i" == *_uniq.fasta.gz ]] && continue
      (
        fmt=$(seq_format "$i")
        [[ "$fmt" == "@" || "$fmt" == ">" ]] || exit 0
        echo "Processing $i ($fmt)"
        to_singleline_seqs "$i" "$fmt" |
        awk '{
          seq=$0
          gsub(/^(A{10,}|C{10,}|G{10,}|T{10,})/,"",seq)
          gsub(/(A{10,}|C{10,}|G{10,}|T{10,})$/,"",seq)
          print seq
        }' | gzip > "tmp_flat_$i"
      ) &
      while (( $(jobs -rp | wc -l) >= gN )); do sleep 1; done
    done
    wait

    # Step 2: Length distribution sampling
    for i in *.f*.gz; do
      [[ "$i" == *_uniq.fasta.gz ]] && continue
      (
        tmp="tmp_flat_$i"
        total_reads=$(zcat "$tmp" | wc -l)
        (( total_reads > 0 )) || exit 0
        sample_size=$(( total_reads * SAMPLE_PERCENT / 100 ))
        (( sample_size < 1 )) && sample_size=1
        zcat "$tmp" | shuf -n "$sample_size" > "${i}_length_distribution.txt"
      ) &
      while (( $(jobs -rp | wc -l) >= gN )); do sleep 1; done
    done
    wait

    # Step 3: Compute max read length (95th percentile)
    awk '{print length($0)}' *_length_distribution.txt | awk -v minlen="$MIN_LEN" '$1 >= minlen' | sort -n > length_distribution.txt
    max_seqread_len=$(awk '{a[NR]=$1} END{print a[int(NR*0.95)]}' length_distribution.txt)
    if [[ -z "${max_seqread_len:-}" ]]; then
      echo "ERROR: no reads passed MIN_LEN=$MIN_LEN" >&2
      exit 1
    fi
    rm -f *_length_distribution.txt length_distribution.txt

    # Step 4: Filter, trim, deduplicate → FASTA with counts
    for i in *.f*.gz; do
      [[ "$i" == tmp_flat_* ]] && continue
      [[ "$i" == *_uniq.fasta.gz ]] && continue
      (
        echo "Filtering, trimming, counting duplicates for $i"
        zcat "tmp_flat_$i" |
        awk -v max="$max_seqread_len" -v minlen="$MIN_LEN" '
          length($0) >= minlen {
            seq = substr($0, 1, max)
            count[seq]++
          }
          END {
            n = 0
            for (s in count) {
              n++
              printf(">seq%d_%d\n%s\n", n, count[s], s)
            }
          }' | gzip > "${i%.f*}_uniq.fasta.gz"
      ) &
      while (( $(jobs -rp | wc -l) >= $gN )); do sleep 1; done
    done
    wait

    # Cleanup & bookkeeping
    rm -f tmp_flat_*.gz
    {
      echo "Converted FASTA/FASTQ → single-line sequences"
      echo "Trimmed homopolymers"
      echo "MIN_LEN = $MIN_LEN"
      echo "95th percentile max length = $max_seqread_len"
    } > "${projdir}/organize_files_done.txt"
    echo "Improvement in flushed reads already implemented" > flushed_reads.txt
  fi
  shopt -u nullglob
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
  shopt -s nullglob
  files=(./samples/*.f*)
  if [[ "$samples_list" == "samples_list_node_1.txt" && ${#files[@]} -gt 0 ]]; then
    for i in samples_list_node_*.txt; do
      :> ${i%.txt}_hold.txt
      while IFS="" read -r line; do
        ls -l ./samples/${line}_R1_uniq.fasta.gz | awk '{print $5"\t"$9}' >> ${i%.txt}_hold.txt
      done < <(grep -v '_tmp.fa' $i | grep -v _R2.f | grep -v _uniq.fasta | grep -v _uniq_R1.fasta | grep -v _uniq_R2.fasta | grep -v _uniq.hold.fasta | grep -v _uniq_R1.hold.fasta | grep -v _uniq_R2.hold.fasta | grep -v fq.gz |  awk '{ sub(/_R1.*/, "", $0); sub(/\.f.*/, "", $0)}1' )
      sort -nr -k1 ${i%.txt}_hold.txt | awk '{gsub(/.\/samples\//,""); print $2}' | awk 'NR>1{print prev} {prev=$0} END{printf "%s", prev}' | awk '{ sub(/_R1.*/, "", $0);}1' > $i
      rm -f "${i%.txt}"_hold.txt
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
if [[ "$(cat ${projdir}/samples_list_node_* | wc -l)" -lt 50 ]] && [[ -z ${maf:-} ]]; then
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
      rm -f "${projdir}"/queue_move_${samples_list%.txt}.txt; sleep $[ ( $RANDOM % 120 )  + 30 ]s
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
    rm -f "${projdir}"/queue_move_${samples_list%.txt}.txt
    wait
  fi

  # Perform read alignments
  cd ${projdir}
  if [[ $nodes -eq 1 ]]; then cd ${projdir}/samples ; fi
  if [[ $nodes -gt 1 ]] && test -f ${projdir}/GBSapp_run_node_1.sh; then cd /tmp/${samples_list%.txt}/samples ; fi
  if [[ "$RNA" == "false" ]]; then
    if [[ "$aligner" == "minimap2" ]]; then
      if test ! -f "${projdir}/precall_done.txt" && test ! -f "${projdir}/alignment_done.txt"; then
        # Automated alignment with structure detection and separate mapping
        while IFS="" read -r alignfq || [ -n "$alignfq" ]; do (
            sleep $((RANDOM % 2))
            sample_base="${alignfq%.f*}"
            r1_file="${sample_base}_R1_uniq.fasta.gz"
            r2_file="${sample_base}_R2_uniq.fasta.gz"
            r1r2_file="${sample_base}_R1R2_uniq.fasta.gz"
            tmp_sv=$(mktemp "${projdir}/samples/tmp.XXXXXX")
            sv_out="${projdir}/sv_mode.txt"

            #### Step 0: Structure detection (sample 1% reads, min 1000)
            detect_structure_mode="low"  # default
            if [[ -f "$r1_file" ]]; then
                total_reads=$(zcat "$r1_file" | awk '/^>/ {c++} END {print c}')
                sample_size=$(( total_reads / 100 ))
                if (( sample_size < 1000 )); then sample_size=1000; fi

                tmp_sample_r1="${sample_base}_structure_sample_R1.fasta"
                zcat "$r1_file" | awk -v n="$sample_size" '/^>/ {if(seq!="" && count<n){print header"\n"seq; count++} header=$0; seq=""; next} {seq=seq $0} END{if(seq!="" && count<n){print header"\n"seq}}' > "$tmp_sample_r1"
                if [[ -f "$r2_file" && -s "$r2_file" ]]; then
                  tmp_sample_r2="${sample_base}_structure_sample_R2.fasta"
                  zcat "$r2_file" | awk -v n="$sample_size" '/^>/ {if(seq!="" && count<n){print header"\n"seq; count++} header=$0; seq=""; next} {seq=seq $0} END{if(seq!="" && count<n){print header"\n"seq}}' > "$tmp_sample_r2"
                fi
                tmp_aln="${sample_base}_structure_check.sam"
                # Use PE if R2 exists, else SE
                if [[ -f "${tmp_sample_r2:-}" ]]; then
                    echo "Structure detection: PE mode for $sample_base"
                    $minimap2 -t "$alnthreads" -ax pe "../refgenomes/${ref1%.f*}.mmi" "$tmp_sample_r1" "$tmp_sample_r2" > "$tmp_aln" 2>/dev/null
                else
                    echo "Structure detection: SE mode for $sample_base"
                    $minimap2 -t "$alnthreads" -ax sr "../refgenomes/${ref1%.f*}.mmi" "$tmp_sample_r1" > "$tmp_aln" 2>/dev/null
                fi
                discordant_pct=$(grep -v '^@' "$tmp_aln" | awk '{if($2>=256) d++} END{if(NR>0){print (d/NR)*100}else{print 0}}')
                if (( $(echo "$discordant_pct > 10" | bc -l) )); then
                    detect_structure_mode="high"
                fi
                rm -f "${tmp_sample_r1:-}" "${tmp_sample_r2:-}" "${tmp_aln:-}"
                echo "Structure mode for $sample_base: $detect_structure_mode"
            fi
            printf "${alignfq%.f*}:\tdetect_structure_mde = ${detect_structure_mode} \t discordant_pct = $discordant_pct \n" > "$tmp_sv"
            cat "$tmp_sv" >> "$sv_out" && rm -f "$tmp_sv"

            if [[ "$detect_structure_mode" == "high" ]]; then
              #### Step 1: Align PE or SE reads
              if [[ -f "$r1_file" && -f "$r2_file" && -s "$r2_file" ]]; then
                  echo "Aligning PE reads: $r1_file + $r2_file (structure: $detect_structure_mode)"
                  $minimap2 -t $alnthreads -ax pe --secondary=yes -k15 -w5 -Y -f 0.0005 -N 8 -n 2 -m 25 "../refgenomes/${ref1%.f*}.mmi" \
                  "$r1_file" "$r2_file" > "${sample_base}.sam"
              else
                  echo "Aligning SE reads: $r1_file (structure: $detect_structure_mode)"
                  $minimap2 -t $alnthreads -ax sr --secondary=yes -k15 -w5 "../refgenomes/${ref1%.f*}.mmi" \
                  "$r1_file" > "${sample_base}.sam"
              fi
              #### Step 2: Align stitched R1R2 reads separately
              if [[ -f "$r1r2_file" && -s "$r1r2_file" ]]; then
                  echo "Aligning stitched reads: $r1r2_file"
                  $minimap2 -t $alnthreads -ax sr --secondary=yes -k15 -w5 "../refgenomes/${ref1%.f*}.mmi" \
                  "$r1r2_file" > "${sample_base}_R1R2.sam"
                  # Merge SAMs
                  cat "${sample_base}.sam" <(grep -v '^@' "${sample_base}_R1R2.sam") > "${sample_base}_all.sam"
                  rm -f "${sample_base}_R1R2.sam"
              else
                  mv "${sample_base}.sam" "${sample_base}_all.sam"
              fi
            fi

            if [[ "$detect_structure_mode" == "low" ]]; then
              #### Step 1: Align PE or SE reads
              if [[ -f "$r1_file" && -f "$r2_file" && -s "$r2_file" ]]; then
                  echo "Aligning PE reads: $r1_file + $r2_file (structure: $detect_structure_mode)"
                  $minimap2 -t $alnthreads -ax pe --secondary=no --sam-hit-only -Y -f 0.0005 -N 8 -n 2 -m 25 "../refgenomes/${ref1%.f*}.mmi" \
                  "$r1_file" "$r2_file" > "${sample_base}.sam"
              else
                  echo "Aligning SE reads: $r1_file (structure: $detect_structure_mode)"
                  $minimap2 -t $alnthreads -ax sr --secondary=no --sam-hit-only "../refgenomes/${ref1%.f*}.mmi" \
                  "$r1_file" > "${sample_base}.sam"
              fi
              #### Step 2: Align stitched R1R2 reads separately
              if [[ -f "$r1r2_file" && -s "$r1r2_file" ]]; then
                  echo "Aligning stitched reads: $r1r2_file"
                  $minimap2 -t $alnthreads -ax sr --secondary=no --sam-hit-only "../refgenomes/${ref1%.f*}.mmi" \
                  "$r1r2_file" > "${sample_base}_R1R2.sam"
                  # Merge SAMs
                  cat "${sample_base}.sam" <(grep -v '^@' "${sample_base}_R1R2.sam") > "${sample_base}_all.sam"
                  rm -f "${sample_base}_R1R2.sam"
              else
                  mv "${sample_base}.sam" "${sample_base}_all.sam"
              fi
            fi


            #### Step 3: Sort and compress
            samtools sort -O SAM "${sample_base}_all.sam" | gzip > "../preprocess/alignment/${sample_base}_redun.sam.gz"
            cp -rn "../preprocess/alignment/${sample_base}_redun.sam.gz" "${projdir}/preprocess/alignment/"
            rm -f "${sample_base}_all.sam") &
            while (( $(jobs -rp | wc -l) >= $alnN )); do sleep 2; done
        done < <(cat "${projdir}/${samples_list}")
        wait
      fi
    fi
  fi
  wait

  cd ${projdir}
  if [[ "$RNA" == "true" ]]; then
    mkdir -p ../preprocess/staralign
    if test ! -f ${projdir}/precall_done.txt && test ! -f ${projdir}/alignment_done.txt; then
      cd samples
      while IFS="" read -r alignfq || [ -n "$alignfq" ]; do (
          sleep $((RANDOM % 2))
          sample_base="${alignfq%.f*}"
          r1_file="${sample_base}_R1_uniq.fasta.gz"
          r2_file="${sample_base}_R2_uniq.fasta.gz"
          r1r2_file="${sample_base}_R1R2_uniq.fasta.gz"
          align_dir="../preprocess/alignment"
          star_dir="../preprocess/staralign"
          out_prefix="${star_dir}/${sample_base}_redun_"
          # Skip if already processed
          if [[ -f "$align_dir/${sample_base}_redun.sam.gz" ]]; then
              continue
          fi

          # Function to convert FASTA -> FASTQ
          fasta_to_fastq() {
              local fasta="$1"
              local fastq_out="$2"
              awk 'BEGIN {RS=">"; ORS=""} NR>1 {
                  header=substr($0,1,index($0,"\n")-1);
                  seq=substr($0,index($0,"\n")+1);
                  gsub(/\n/,"",seq);
                  if(length(seq) >= 64) print ">" header "\n" seq "\n";
              }' <(zcat "$fasta") | \
              awk 'BEGIN {OFS="\n"} /^>/ {header=substr($0,2); next} {print "@" header,$0,"+",
              gensub(/./,"I","g",$0)}' | gzip > "$fastq_out"
          }

          # Step 0: Determine input reads to align
          if [[ -f "$r1_file" && -f "$r2_file" && -s "$r2_file" ]]; then
              echo "PE reads detected for $sample_base: $r1_file + $r2_file"
              fasta_to_fastq "$r1_file" "${sample_base}_R1_uniq.fastq.gz"
              fasta_to_fastq "$r2_file" "${sample_base}_R2_uniq.fastq.gz"
              fastq_in=("${sample_base}_R1_uniq.fastq.gz" "${sample_base}_R2_uniq.fastq.gz")
              star_args="--readFilesIn ${fastq_in[*]}"
          elif [[ -f "$r1_file" ]]; then
              echo "SE reads detected for $sample_base: $r1_file"
              fasta_to_fastq "$r1_file" "${sample_base}_R1_uniq.fastq.gz"
              star_args="--readFilesIn ${sample_base}_R1_uniq.fastq.gz"
          elif [[ -f "$r1r2_file" && -s "$r1r2_file" ]]; then
              echo "Stitched reads detected for $sample_base: $r1r2_file (high-SV mode)"
              fasta_to_fastq "$r1r2_file" "${sample_base}_R1R2_uniq.fastq.gz"
              star_args="--readFilesIn ${sample_base}_R1R2_uniq.fastq.gz"
          else
              echo "No valid reads found for $sample_base, skipping..."
              continue
          fi

          # Step 1: Run STAR
          $star --runThreadN "$alnthreads" \
              --genomeDir "${projdir}/refgenomes/star_index" \
              $star_args \
              --readFilesCommand zcat \
              --outFileNamePrefix "$out_prefix" \
              --outSAMtype SAM \
              --outFilterMultimapNmax 8 \
              --twopassMode Basic
          wait

          # Step 2: Capture header
          awk '/@HD/ || /@SQ/ {print}' "${out_prefix}Aligned.out.sam" 2>/dev/null \
              > "$align_dir/${sample_base}_redun_head.sam"
          # Step 3: Extract aligned reads only
          grep -v '^@' "${out_prefix}Aligned.out.sam" 2>/dev/null | \
          awk 'BEGIN{FS=OFS="\t"} !($10=="*" && $6 !~ /^\*$/){print}' | \
          cat "$align_dir/${sample_base}_redun_head.sam" - > "$align_dir/${sample_base}_redun.sam"
          # Step 4: Convert to BAM and split N cigars
          $samtools view -S -b "$align_dir/${sample_base}_redun.sam" > "$align_dir/${sample_base}_redun.bam"
          $GATK SplitNCigarReads \
              -R "${projdir}/refgenomes/$ref1" \
              -I "$align_dir/${sample_base}_redun.bam" \
              -O "$align_dir/${sample_base}_redunsplit.bam" --verbosity ERROR
          # Step 5: Compress final SAM
          $samtools view -h "$align_dir/${sample_base}_redunsplit.bam" | gzip > "$align_dir/${sample_base}_redun.sam.gz"
          # Step 6: Cleanup
          rm -f "$align_dir/${sample_base}_redun.sam" \
                "$align_dir/${sample_base}_redun.bam" \
                "$align_dir/${sample_base}_redunsplit.bam" \
                "$align_dir/${sample_base}_redunsplit.bai" \
                "$align_dir/${sample_base}_redun_head.sam"
          rm -f "${sample_base}_R1.fastq.gz" "${sample_base}_R2.fastq.gz" "${sample_base}_R1R2.fastq.gz"
          rm -rf "$star_dir/*"
          cp -rn "$align_dir/${sample_base}_redun.sam.gz" "${projdir}/preprocess/alignment/"
          ) &
          while (( $(jobs -rp | wc -l) >= $alnN )); do sleep 2; done
      done < <(cat "${projdir}/${samples_list}")
    fi
  fi
  wait

  cat ${projdir}sv_mode.txt > ${projdir}/alignment_done_${samples_list}

  cd ${projdir}/preprocess
  if [[ "$samples_list" == "samples_list_node_1.txt" ]] && test ! -f ${projdir}/alignment_done.txt; then
    shopt -s nullglob
    while files=("${projdir}"/alignment_done_samples_list_node_*); [[ ${#files[@]} -lt "$nodes" ]]; do
      sleep 300
    done
    cat ${projdir}/sv_mode.txt > ${projdir}/alignment_done.txt
  fi
  rm -f ${projdir}/samples/tmp.* ${projdir}/sv_mode.txt

  while [[ ! -f ${projdir}/alignment_done.txt ]]; do
    sleep 300
  done
  wait
  if [[ $nodes -eq 1 ]]; then cd ${projdir}/preprocess/; fi
  if [[ $nodes -gt 1 ]] && test -f ${projdir}/GBSapp_run_node_1.sh; then cd /tmp/${samples_list%.txt}/preprocess/; fi

  if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
    touch ${projdir}/compress_done.txt
  fi
  if [[ "$ploidy" -eq 1 ]]; then
    if [[ -z "${downsample_1x:-}" ]]; then
      downsample=50
    else
      downsample="$downsample_1x"
    fi
  fi
  if [[ "$ploidy" -eq 2 ]]; then
    if [[ -z "${downsample_2x:-}" ]]; then
      downsample=100
    else
      downsample="$downsample_2x"
    fi
  fi
  if [[ "$ploidy" -eq 4 ]]; then
    if [[ -z "${downsample_4x:-}" ]]; then
      downsample=200
    else
      downsample="$downsample_4x"
    fi
  fi
  if [[ "$ploidy" -eq 6 ]]; then
    if [[ -z "${downsample_6x:-}" ]]; then
      downsample=300
    else
      downsample=$downsample_6x
    fi
  fi
  if [[ "$ploidy" -eq 8 ]]; then
    if [[ -z "${downsample_8x:-}" ]]; then
      downsample=400
    else
      downsample=$downsample_8x
    fi
  fi

  while IFS="" read -r i || [ -n "$i" ]; do (
    printf '\n###---'${i%.f*}'---###\n' > ${projdir}/alignment_summaries/${i%.f*}_summ.txt &&
    zcat ./alignment/${i%.f*}_redun.sam.gz | grep -v '^@' | awk '{gsub(/_/,"\t",$1);}1' | tr ' ' '\t' > ${i%.f*}_full.sam &&
    split -l 10000 ${i%.f*}_full.sam ${i%.f*}_chunk_ && rm -f ${i%.f*}_full.sam
    find . -name '${i%.f*}_chunk_*' -print0 | xargs -0 -P "$gthreads" -I{} bash -c 'awk '\''{for(i=0;i<=$2-1;i++) print $0}'\'' "$1" > "$1.out"' _ {} &&
    cat ${i%.f*}_chunk_* | awk '!($2="")1' | awk '{$1=$1"_"NR}1' | awk '{gsub(/ /,"\t");}1' > ${i%.f*}_full.sam &&
    cat <(zcat ./alignment/${i%.f*}_redun.sam.gz | grep '^@') ${i%.f*}_full.sam | gzip > ${i%.f*}_full.sam.gz &&
    rm -f ${i%.f*}_chunk_* ${i%.f*}_full.sam
    $samtools flagstat ${i%.f*}_full.sam.gz >> ${projdir}/alignment_summaries/${i%.f*}_summ.txt &&
    printf '########################################################################################################\n\n' >> ${projdir}/alignment_summaries/${i%.f*}_summ.txt &&

    if test ! -f ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam.bai; then
      if [[ "$paralogs" == false ]] && [[ "$uniquely_mapped" == true ]]; then
        awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) > ${i%.f*}_heading.sam &&
        $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
        awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0 || $5 >= min) {print $0}}' | awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | \
        tr " " "\t" | tr '\r' '\n' | awk '$1==1{print $0}' | awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | awk '{gsub(/_/,"_\t",$1);}1' | awk '{print $2"\t"$0}' | \
        awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" > ${i%.f*}_uniq.sam
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
        awk '{gsub(/_/,"_\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" > ${i%.f*}_uniq.sam &&
        rm -f ${i%.f*}_uniqeq.sam
      fi
      if [[ "$paralogs" == true ]] && [[ "$uniquely_mapped" == false ]]; then
        awk '/@HD/ || /@SQ/{print}' <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) > ${i%.f*}_heading.sam &&
        $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
        awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0 || $5 >= min) {print $0}}' | awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | \
        tr " " "\t" | tr '\r' '\n' | awk '$1==1{print $0}' | awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | awk '{gsub(/_/,"_\t",$1);}1' | awk '{print $2"\t"$0}' | \
        awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" > ${i%.f*}_uniq.sampart &&
        $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
        awk -v min=$minmapq -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0 || $5 >= min) {print $0}}' | awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | \
        tr " " "\t" | tr '\r' '\n' | awk '$1==1 || $1<=6{print $0}' | awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | \
        awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 == 0) {print $0}}' 2> /dev/null > ${i%.f*}_uniqeq.sam
        $samtools view -F4 <(zcat ./alignment/${i%.f*}_redun.sam.gz 2> /dev/null) | grep -v '^@' | awk '$3 != "*"' 2> /dev/null | awk '$6 != "*"' 2> /dev/null | \
        awk '!h[$1] { g[$1]=$0 } { h[$1]++ } END { for(k in g) print h[k], g[k] }' | tr " " "\t" | tr '\r' '\n' | awk '$1==1 || $1<=6{print $0}' | awk '{$1=""}1' | awk '$1=$1' | tr " " "\t" | \
        awk -F '\t' 'BEGIN{OFS="\t"} {if ($5 > 0) {print $0}}' 2> /dev/null | cat - ${i%.f*}_uniqeq.sam | \
        awk '{gsub(/_/,"_\t",$1);}1' | awk '{print $2"\t"$0}' | awk '{$2=$3=""}1' | tr -s " " | tr " " "\t" > ${i%.f*}_uniqAll.sam &&
        rm -f ${i%.f*}_uniqeq.sam
        awk 'NR==FNR{a[$0]=1;next}!a[$0]' ${i%.f*}_uniqpart.sam ${i%.f*}_uniqAll.sam > ${i%.f*}_uniq.sam &&
        ${i%.f*}_uniqpart.sam ${i%.f*}_uniqAll.sam
      fi

      awk '{while ($1-- > 0) print $0}' ${i%.f*}_uniq.sam | \
      awk 'BEGIN{OFS="\t"}{split($0,a,"\t"); base=a[1]; copy[base]++; a[1]=base"_copy"copy[base]; print a[1],$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' > ${i%.f*}_${ref1%.f*}.sam &&
      rm -f ${i%.f*}_uniq.sam
      BIN_SIZE=50   # bp bin for locus-aware downsampling
      # NOTE: Synthetic base qualities (Q=40); BQSR must NOT be run downstream
      awk -v b=$BIN_SIZE '{printf "%d\n", int($4/b)}' ${i%.f*}_${ref1%.f*}.sam | paste - ${i%.f*}_${ref1%.f*}.sam | \
      shuf | awk -v max="$downsample" 'a[$1$4]++ < max' | awk '{$1=""}1' | awk '{$1=$1};1' | tr -s " " | tr " " "\t" | awk '{$11=$10}1' | \
      awk -v Q="I" 'BEGIN{OFS="\t"}{gsub(/[ACGTNacgtn]/,Q,$11); print}' | cat ${i%.f*}_heading.sam - > ${i%.f*}_${ref1%.f*}_downsample.sam &&
      :> ${i%.f*}_${ref1%.f*}.sam &&
      mv ${i%.f*}_${ref1%.f*}_downsample.sam ${i%.f*}_${ref1%.f*}.sam &&
      rm -f ${i%.f*}_heading.sam

      j="${i%.f*}_${ref1%.f*}.sam"
      CIGAR_FILTER='I[0-9]+I|D[0-9]+D|D[0-9]+I|I[0-9]+D'
      cat ${j} | grep -Ev "$CIGAR_FILTER"  > cleaned_${j};  mv cleaned_${j} ${j} &&
      $java $Xmx2 -XX:ParallelGCThreads=$gthreads -Djava.io.tmpdir=${projdir}/preprocess/tmp -jar $picard SortSam I=$j O=${j%.sam*}.bam  SORT_ORDER=coordinate  VALIDATION_STRINGENCY=LENIENT TMP_DIR=${projdir}/preprocess/tmp >/dev/null 2>&1 && \
      $java $Xmx2 -XX:ParallelGCThreads=$gthreads -jar $picard BuildBamIndex INPUT=${j%.sam*}.bam VALIDATION_STRINGENCY=LENIENT >/dev/null 2>&1 && \
      $java $Xmx2 -XX:ParallelGCThreads=$gthreads -jar $picard AddOrReplaceReadGroups I=${j%.sam*}.bam O=${j%.sam*}_precall.bam RGLB=${i%.f*} RGPL=illumina RGPU=run RGSM=${i%.f*} VALIDATION_STRINGENCY=LENIENT >/dev/null 2>&1 && \
      $samtools index ${j%.sam*}_precall.bam &&
      rm -f $j ${j%.sam*}.bam ${j%.sam*}.bai
      if [[ $nodes -gt 1 ]]; then cp /tmp/${samples_list%.txt}/preprocess/${j%.sam*}_precall.bam* ${projdir}/preprocess/; fi
    fi
    ) &
    while (( $(jobs -rp | wc -l) >= $gN )); do sleep 2; done
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
      rm -f *_summ.txt
			# Total number of reads per samples
      awk '/###---/ || /QC-passed/{print}' alignment_summaries_reads.txt | cut -d\+ -f1 | tr -d '\n' | \
      awk  'gsub(/---###/, "\t", $0)' | awk  'gsub(/###---/, "", $0)' | tr ' ' '\n' > total_reads.txt &&
      # Total number of mapped reads per samples
      cat alignment_summaries_reads.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /0_primary_mapped/{print}' |\
      tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
      awk 'gsub("\\+0primarymapped", "\t", $0)' | cut -d\: -f1 > total_reads_mapped.txt &&
      # Total number of mapped paired reads per samples
      cat alignment_summaries_reads.txt | tr ' ' '_' | tr '(' '_' |  tr ')' '_' | awk '/###---/ || /properly_paired/{print}' |\
      tr -d '\n' | awk 'gsub(/###---/, "\n", $0)' | awk 'gsub(/---###/, "\t", $0)' | awk 'gsub(/_/, "", $2)' | \
      awk 'gsub("\\+0properlypaired", "\t", $0)' | cut -d\: -f1 > total_reads_paired.txt &&
      echo -e "Samples\tTotal\tMapped\tPerc_Mapped\tPE_Mapped\t%_PE_Mapped" > summary_precall.txt &&
      awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' total_reads_mapped.txt  total_reads.txt  | \
      awk 'FNR==NR{a[$1]=$2 FS $3;next} ($1 in a) {print $0,a[$1]}' total_reads_paired.txt - | awk '{gsub(/ /,"\t"); print $0}' | \
      cat summary_precall.txt - | awk '{print $1"\t"$2"\t"$3"\t"$4}' | \
      awk 'BEGIN{OFS="\t"} NR==1 {print; next} {printf "%s\t%d\t%d\t%.2f%%\n", $1, $2, $3, ($3/$2)*100}' > Tabulated_Alignment_Read_Summaries.txt &&
			rm -f total_* summary_precall.txt
			rm -f ${projdir}/samples/metrics.txt ${projdir}/preprocess/metrics.txt
      wait
		fi
    if [[ "$genomecov_est" == true ]] && [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
      cd ${projdir}/alignment_summaries
      printf "Sample\tGenome_Coverage(percentage)\n" > summary_genomecov.txt
      genome_size=$(awk '{print $3}' ../refgenomes/${ref1%.f*}.dict | awk '{gsub(/LN:/,"");}1' | awk '{s+=$1}END{print s}')
      for i in ../preprocess/alignment/*_redun.sam.gz; do
          while (( $(jobs -rp | wc -l) >= $gN )); do sleep 2; done
          (
              tmp_out=$(mktemp "${projdir}/alignment_summaries/tmp.XXXXXX")
              $samtools view -bS <(zcat "$i" 2> /dev/null) | $samtools sort - > "${i%*_redun.sam.gz}.bam"
              cov=$($bedtools genomecov -ibam "${i%*_redun.sam.gz}.bam" -bga | awk '{print ($4>1)?1:$4}' \
              | awk -v pat=$genome_size '{s+=$1}END{print (s/pat)*100}')
              printf "%s\t%s\n" "${i%*_redun.sam.gz}" "$cov" | awk '{gsub(/..\/preprocess\/alignment\//,"");}1' > "$tmp_out"
              cat "$tmp_out" >> summary_genomecov.txt
              rm -f "$tmp_out" "${i%*_redun.sam.gz}.bam" ) &
      done
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
  if [[ ! -s ${interval_list:-} ]]; then export interval_list=""; fi

  if [[ "$samples_list" == "samples_list_node_1.txt" ]]; then
  	if [[ "$joint_calling" == true ]] && [[ "$variant_caller" == "gatk" ]]; then
  		j=-I; input=""; k=""
  		for i in *_precall.bam; do
  			k="${j} ${i}"; input="${input} ${k}"
  		done
  		Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/${ref1%.*}.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}')
  		if [[ -n "${Exclude_Chromosome:-}" ]]; then
  			for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
  				Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
  			done
  		fi
  		if [[ -z "${Get_Chromosome:-}" ]]; then
        if [[ -z "${interval_list:-}" ]]; then
    			for selchr in $Get2_Chromosome; do (
            if [[ -z "$(compgen -G "${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf*")" ]]; then
    					$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller \
              -R ${projdir}/refgenomes/$ref1 -L ${selchr} ${input} -ploidy $ploidy \
              -O ${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf.gz --max-reads-per-alignment-start 0 \
              --minimum-mapping-quality $minmapq --dont-use-soft-clipped-bases $dont_use_softclip \
              --disable-bam-index-caching true --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) --verbosity ERROR &&
    					gunzip ${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf.gz
    				fi
    				) &
    				while (( $(jobs -rp | wc -l) >= $gN )); do sleep 2; done
    			done
        else
          for selchr in $Get2_Chromosome; do (
            cat ${projdir}/${interval_list} | grep $selchr > ${projdir}/variant_intervals_${selchr}.list &&
            if [[ -z "$(compgen -G "${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf*")" ]]; then
              $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller \
              -R ${projdir}/refgenomes/$ref1 -L ${projdir}/variant_intervals_${selchr}.list ${input} -ploidy $ploidy \
              -O ${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf.gz --max-reads-per-alignment-start 0 \
              --minimum-mapping-quality $minmapq --dont-use-soft-clipped-bases $dont_use_softclip --disable-bam-index-caching true \
              --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) --verbosity ERROR&&
              gunzip ${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf.gz &&
              rm -f ${projdir}/variant_intervals_${selchr}.list
            fi
            ) &
            while (( $(jobs -rp | wc -l) >= $gN )); do sleep 2; done
          done
        fi
  			wait
  			cd ../snpcall
  			grep -h '^#' ${pop}_${ploidy}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt &&
  			cat ${pop}_${ploidy}x_*_raw.vcf | awk '!/^#/' > all.vcf &&
  			cat vcf_header.txt all.vcf > ${pop}_${ref1%.f*}_${ploidy}x_raw.vcf
  			rm -f vcf_header.txt all.vcf *.vcf.gz.tb* ${pop}_${ploidy}x_*_raw.vcf
        wait
  		else
        if [[ -z "$(compgen -G "${projdir}/snpcall/${pop}_${ploidy}x_${selchr}_raw.vcf*")" ]]; then
  				echo $Get_Chromosome | tr ',' '\n' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ${projdir}/refgenomes/${ref1%.fasta}.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2":1-"$3}' > ${projdir}/refgenomes/${ref1%.fasta}.list &&
  				$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller \
          -R ${projdir}/refgenomes/$ref1 -L ${projdir}/refgenomes/${ref1%.fasta}.list ${input} -ploidy $ploidy \
          -O ${projdir}/snpcall/${pop}_${ploidy}x_raw.vcf.gz --max-reads-per-alignment-start 0 --minimum-mapping-quality $minmapq \
          --dont-use-soft-clipped-bases $dont_use_softclip --disable-bam-index-caching true \
          --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) --verbosity ERROR &&
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
          if [[ -z "$(compgen -G "${projdir}/snpcall/${pop}_${ref1%.f*}_${ploidy}x_raw.vcf*")" ]]; then
      			if test ! -f "${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf"; then
      					if [[ -z "${Get_Chromosome:-}" ]]; then
                  if [[ -z "${interval_list:-}" ]]; then
        						$GATK --java-options "$Xmxg -Djava.io.tmpdir=../snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller \
                    -R ../refgenomes/$ref1 -I ${i%.f*}_${ref1%.f*}_precall.bam -ploidy $ploidy -O ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz \
                    -ERC GVCF --max-reads-per-alignment-start 0 --minimum-mapping-quality $minmapq --dont-use-soft-clipped-bases $dont_use_softclip \
                    --disable-bam-index-caching true --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) --verbosity ERROR
                  else
                    $GATK --java-options "$Xmxg -Djava.io.tmpdir=../snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller \
                    -R ../refgenomes/$ref1 -L ${projdir}/${interval_list} -I ${i%.f*}_${ref1%.f*}_precall.bam -ploidy $ploidy \
                    -O ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz -ERC GVCF --max-reads-per-alignment-start 0 \
                    --minimum-mapping-quality $min--dont-use-soft-clipped-bases --disable-bam-index-caching $dont_use_softclip \
                    --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) --verbosity ERROR
                  fi
      					else
      						echo $Get_Chromosome | tr ',' '\n' | awk '{print "SN:"$1}' | awk 'NR==FNR{a[$1];next}$2 in a{print $0}' - ../refgenomes/${ref1%.fasta}.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $2":1-"$3}' > ../refgenomes/${ref1%.fasta}.list
      						$GATK --java-options "$Xmxg -Djava.io.tmpdir=../snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" HaplotypeCaller \
                  -R ../refgenomes/$ref1 -L ../refgenomes/${ref1%.fasta}.list -I ${i%.f*}_${ref1%.f*}_precall.bam -ploidy $ploidy \
                  -O ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz -ERC GVCF --max-reads-per-alignment-start 0 \
                  --minimum-mapping-quality $minmapq --dont-use-soft-clipped-bases --disable-bam-index-caching $dont_use_softclip \
                  --max-num-haplotypes-in-population $((ploidy * maxHaplotype)) --verbosity ERROR
      					fi
      					mv ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf.gz &&
                rm -f ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.hold.g.vcf.gz.tb*
                gunzip ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf.gz &&
                $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" IndexFeatureFile \
                -I ${projdir}/snpcall/${i%.f*}_${ref1%.f*}.g.vcf --verbosity ERROR
      			fi
      		fi
      		mv ${projdir}/preprocess/${i%.f*}_${ref1%.f*}_precall.bam* ${projdir}/preprocess/processed/ 2> /dev/null
          rm -f ${i%.f*}_${ref1%.f*}_precall.bam*
          ) &
    		while (( $(jobs -rp | wc -l) >= $gN )); do sleep 2; done
    	done < <(cat ${projdir}/${samples_list})
    fi
  	wait
    cd ${projdir}/preprocess/
  	printf "variant calling completed on ${samples_list}" > ${projdir}/call1_${samples_list}

    shopt -s nullglob
    files=("${projdir}/snpcall/${pop}_${ref1%.f*}_${ploidy}x_raw.vcf"*)
    if [[ "$samples_list" == "samples_list_node_1.txt" && ${#files[@]} -eq 0 ]]; then
      shopt -s nullglob
      call1=0
      while [[ "$call1" -lt "$nodes" ]]; do
          files=("${projdir}/call1_samples_list_node_"*)
          call1=${#files[@]}
          sleep 300
      done
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
  				if [[ -z "${Get_Chromosome:-}" ]]; then
  					Get2_Chromosome=$(awk 'NR>1{print $2,"\t",$3}' ${projdir}/refgenomes/${ref1%.*}.dict | awk '{gsub(/SN:/,"");gsub(/LN:/,""); print $0}' | sort -k2,2 -nr | awk '{print $1}')
  				else
  					Get2_Chromosome=$(echo $Get_Chromosome | tr ',' '\n')
  				fi
  				if [[ -n "${Exclude_Chromosome:-}" ]]; then
  					for i in $(echo "$Exclude_Chromosome" | tr ',' '\n'); do
  						Get2_Chromosome=$(echo $Get2_Chromosome | awk -v i=$i '{gsub(i,"");}1')
  					done
  				fi

  				for selchr in $Get2_Chromosome; do (
            if [[ ! -f "${pop}_${ploidy}x_${selchr}_raw.vcf.gz.tbi" ]]; then
              rm -rf ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz
              rm -rf ${pop}_${ploidy}x_${selchr}_raw
              export TILEDB_DISABLE_FILE_LOCKING=1
              if [[ -z "${interval_list:-}" ]]; then
                $GATK --java-options "$Xmx1 -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$loopthreads" GenomicsDBImport \
                ${input} -L ${selchr} --genomicsdb-workspace-path ${pop}_${ploidy}x_${selchr}_raw --genomicsdb-shared-posixfs-optimizations true \
                --batch-size 50 --merge-input-intervals --verbosity ERROR
              else
                cat ${projdir}/${interval_list} | grep $selchr > ${projdir}/variant_intervals_${selchr}.list &&
                $GATK --java-options "$Xmx1 -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$loopthreads" GenomicsDBImport \
                ${input} -L ${projdir}/variant_intervals_${selchr}.list --genomicsdb-workspace-path ${pop}_${ploidy}x_${selchr}_raw \
                --genomicsdb-shared-posixfs-optimizations true --batch-size 50 --merge-input-intervals --verbosity ERROR &&
                rm -f ${projdir}/variant_intervals_${selchr}.list
              fi
            fi
            ) &
        		while (( $(jobs -rp | wc -l) >= $N )); do sleep 2; done
  				done
          wait


  				for selchr in $Get2_Chromosome; do (
            export TILEDB_DISABLE_FILE_LOCKING=1
  					if [[ ! -f "${pop}_${ploidy}x_${selchr}_raw.vcf.gz.tbi" ]]; then
              if [[ -z "${interval_list:-}" ]]; then
    						$GATK --java-options "$Xmx1 -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$loopthreads" GenotypeGVCFs \
                -R ${projdir}/refgenomes/$ref1 -L ${selchr} -V gendb://${pop}_${ploidy}x_${selchr}_raw \
                -O ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz --verbosity ERROR &&
    						rm -rf ${pop}_${ploidy}x_${selchr}_raw
    						mv ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz ${pop}_${ploidy}x_${selchr}_raw.vcf.gz &&
    						mv ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz.tbi ${pop}_${ploidy}x_${selchr}_raw.vcf.gz.tbi
              else
                cat ${projdir}/${interval_list} | grep $selchr > ${projdir}/variant_intervals_${selchr}.list &&
                $GATK --java-options "$Xmx1 -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$loopthreads" \
                GenotypeGVCFs -R ${projdir}/refgenomes/$ref1 -L ${projdir}/variant_intervals_${selchr}.list \
                -V gendb://${pop}_${ploidy}x_${selchr}_raw -O ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz --verbosity ERROR &&
                rm -f ${projdir}/variant_intervals_${selchr}.list
                rm -rf ${pop}_${ploidy}x_${selchr}_raw
                mv ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz ${pop}_${ploidy}x_${selchr}_raw.vcf.gz &&
                mv ${pop}_${ploidy}x_${selchr}_raw.hold.vcf.gz.tbi ${pop}_${ploidy}x_${selchr}_raw.vcf.gz.tbi
              fi
            fi
  					if LC_ALL=C gzip -l ${pop}_${ploidy}x_${selchr}_raw.vcf.gz | awk 'NR==2 {exit($2!=0)}'; then
  						:
  					else
  						rm -f ../cohorts*/${pop}_${ploidy}x_*_raw.vcf.gz* 2> /dev/null
  						rm -f ../cohorts*/${pop}_${ploidy}x_raw.vcf* 2> /dev/null
  						rm -f ../${pop}_${ploidy}x_raw_cohorts*.vcf* 2> /dev/null
  						echo -e "${magenta}- \n- SNP calling failed probably due to insufficient memory ${white}\n"
  						echo -e "${magenta}- \n- Exiting pipeline in 5 seconds ${white}\n"
  						sleep 5 && exit 1
  					fi
            ) &
        		while (( $(jobs -rp | wc -l) >= $N )); do sleep 2; done
  				done
  				wait
  				for g in ${pop}_${ploidy}x_*_raw.vcf.gz; do
  					gunzip $g
  				done
  				grep -h '^#' ${pop}_${ploidy}x_*_raw.vcf | awk '!visited[$0]++' | awk '!/^##GATKCommandLine/' > vcf_header.txt &&
  				cat ${pop}_${ploidy}x_*_raw.vcf | awk '!/^#/' > all.vcf &&
  				cat vcf_header.txt all.vcf > ${pop}_${ploidy}x_raw.vcf &&
  				rm -f vcf_header.txt all.vcf
  				rm -f ${pop}_${ploidy}x_*_raw.vcf ${pop}_${ploidy}x_*_raw.vcf.gz.tb* 2> /dev/null
  				$bcftools view -I ${pop}_${ploidy}x_raw.vcf -O z -o ${pop}_${ploidy}x_raw.vcf.gz &&
  				$bcftools index ${pop}_${ploidy}x_raw.vcf.gz &&

  				$bcftools annotate -x FORMAT/PL ${pop}_${ploidy}x_raw.vcf.gz > ../${pop}_${ploidy}x_raw_${dir%/}.vcf &&
  				cd ../
  				$bcftools view -I ${pop}_${ploidy}x_raw_${dir%/}.vcf -O z -o ${pop}_${ploidy}x_raw_${dir%/}.vcf.gz &&
  				$bcftools index ${pop}_${ploidy}x_raw_${dir%/}.vcf.gz
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
  				rm -rf *cohorts*
  			else
  				rm -rf cohorts*
  				rm -rf *cohorts*
  			fi
  			cd ${projdir}/preprocess
  		fi
  	fi
  	wait

  	###########
    # perform lift over for pangenome-aware variant calling
    set -euo pipefail
    cd $projdir
    if [[ -d "./refgenomes/pangenomes" ]] && find "./refgenomes/pangenomes" -type f -size +0c -print -quit | grep -q .; then
      if [[ "$(zcat "./snpcall/${pop}_${ref1%.f*}_${ploidy}x_raw.vcf"* | awk '$1 ~ /^pangenome_/ {print 1}' | head | wc -l)" -gt 0 ]] ; then
        if [[ ! -f "${projdir}/projection_done.txt" ]]; then
        cd snpcall
        primary_prefix="${ref1%.f*}_Chr"
        primary_ref="../refgenomes/${ref1%.f*}_original.fasta"
        secondary_ref="../refgenomes/panref.fasta"
        paf_file="../refgenomes/pangenome2primary.paf"
        chain_out="../refgenomes/pangenome_to_primary.chain"
        vcf_file="${pop}_${ref1%.f*}_${ploidy}x_raw.vcf"
        lifted_vcf="${vcf_file%.vcf}_lifted.vcf"
        failed_vcf="${vcf_file%.vcf}_failed.vcf"
        # compress and index
        if [[ -f "${vcf_file}.gz" ]]; then
          gunzip "${vcf_file}.gz"
        fi
        $bcftools view -Oz -o "${vcf_file}.gz" "$vcf_file" &&
        $bcftools index "${vcf_file}.gz"
        rm -f "$vcf_file"
        wait
        sort -k6,6 -k8,8n "$paf_file" | awk '!seen[$6,$8,$9]++' > "${paf_file%.paf}.clean.paf" &&
        python3 ${GBSapp_dir}/tools/paf2chain.py "${paf_file%.paf}.clean.paf" "$secondary_ref" "$primary_ref" "$chain_out" &&

        # extract TF contig names from the primary reference
        awk '{print $1 "\t1\t" $2}' "${primary_ref}.fai" | grep "^${primary_prefix}" > ${ref1%.f*}.regions &&
        # subset variants strictly to primary reference contigs &&
        $bcftools view -R ${ref1%.f*}.regions -Oz -o ${ref1%.f*}_only.tmp.vcf.gz "${vcf_file}.gz" &&
        $bcftools index ${ref1%.f*}_only.tmp.vcf.gz &&
        # reheader safely with full TF FASTA
        $bcftools reheader -f "${primary_ref}.fai" -o ${ref1%.f*}_only.vcf.gz ${ref1%.f*}_only.tmp.vcf.gz &&
        $bcftools index ${ref1%.f*}_only.vcf.gz &&

        lifted_pangenome_vcf="${vcf_file%.vcf}_lifted.vcf.gz"
        ${GBSapp_dir}/tools/transanno/target/release/transanno liftvcf \
        --original-assembly "$primary_ref" --new-assembly "$secondary_ref" \
        --chain "$chain_out" --vcf ${ref1%.f*}_only.vcf.gz --output "$lifted_pangenome_vcf" --fail failed.vcf.gz &&
        $bcftools sort "$lifted_pangenome_vcf" -Oz -o "${lifted_pangenome_vcf%.vcf.gz}.sorted.vcf.gz" &&
        $bcftools norm -f "$secondary_ref" -m-any "${lifted_pangenome_vcf%.vcf.gz}.sorted.vcf.gz" \
        -Oz -o "${lifted_pangenome_vcf%.vcf.gz}.sorted.norm.vcf.gz" &&
        $bcftools view -e 'GT="./."' "${lifted_pangenome_vcf%.vcf.gz}.sorted.norm.vcf.gz" -Oz -o "$lifted_pangenome_vcf" &&
        $bcftools index "$lifted_pangenome_vcf" &&

        # Filter variants with max 80% missingness
        $bcftools view -i 'F_MISSING < 0.8' "${vcf_file}.gz" -Oz -o combined.filtmiss20perc.vcf.gz &&
        $bcftools index combined.filtmiss20perc.vcf.gz &&
        $bcftools view -i 'F_MISSING < 0.8' "$lifted_pangenome_vcf" -Oz -o pangenome.filtmiss20perc.vcf.gz &&
        $bcftools index pangenome.filtmiss20perc.vcf.gz &&
        # Pangenome SNPs
        $bcftools view -r $($bcftools query -f '%CHROM\n' pangenome.filtmiss20perc.vcf.gz | sort -u | grep '^pangenome_' | paste -sd ',') \
        pangenome.filtmiss20perc.vcf.gz -Oz -o pangenome.vcf.gz &&
        bcftools index pangenome.vcf.gz &&
        # Primary reference SNPs
        $bcftools view -r $($bcftools query -f '%CHROM\n' combined.filtmiss20perc.vcf.gz | sort -u | grep "^${primary_prefix}" | paste -sd ',') \
        combined.filtmiss20perc.vcf.gz -Oz -o primary.vcf.gz &&
        $bcftools index primary.vcf.gz &&
        # Normalize primary
        $bcftools norm -f "$primary_ref" -m-any primary.vcf.gz -Oz -o primary.norm.vcf.gz &&
        $bcftools index primary.norm.vcf.gz &&
        # Normalize pangenome
        $bcftools norm -f "$secondary_ref" -m-any pangenome.vcf.gz -Oz -o pangenome.norm.vcf.gz &&
        $bcftools index pangenome.norm.vcf.gz &&
        # Merge AFTER normalization
        final_vcf="${pop}_${ref1%.f*}_${ploidy}x_raw.vcf.gz"
        $bcftools concat -a -Oz -o "$final_vcf" primary.norm.vcf.gz pangenome.norm.vcf.gz &&
        $bcftools index  "$final_vcf"
        printf "Pangenome projection with structural absence completed\n" > "${projdir}/projection_done.txt"
        rm -f "${ref1%.f*}_only.tmp.vcf.gz"* "${ref1%.f*}_only.vcf.gz"* merged_final.vcf.gz* \
        pangenome.filtmiss20perc.vcf.gz* combined.filtmiss20perc.vcf.gz* failed.vcf.gz \
        primary.vcf.gz* pangenome.vcf.gz* primary.norm.vcf.gz* pangenome.norm.vcf.gz* \
        "${vcf_file%.vcf}_lifted"* "${ref1%.f*}.regions"
        wait
        fi
      fi
    fi
    ###########

  	if [[ $nodes -gt 1 ]] && test -f ${projdir}/GBSapp_run_node_1.sh; then
  		touch ${projdir}/queue_move_${samples_list%.txt}.txt
  		queue_move=$(ls ${projdir}/queue_move_samples_list_node_* | wc -l)
  		while [[ "$queue_move" -gt 1 ]]; do
  			sleep $((RANDOM % 2))
        rm -f ${projdir}/queue_move_${samples_list%.txt}.txt; sleep $[ ( $RANDOM % 120 )  + 30 ]s
  			:> ${projdir}/queue_move_${samples_list%.txt}.txt
  			queue_move=$(ls ${projdir}/queue_move_samples_list_node_* | wc -l)
  		done
  		cd /tmp/${samples_list%.txt}/preprocess/
  		mv * ${projdir}/preprocess/ && cd ${projdir}
  		rm -rf /tmp/${samples_list%.txt} 2> /dev/null
  		rm -f ${projdir}/queue_move_${samples_list%.txt}.txt
  	fi
  fi

  if [[ "$variant_caller" == "bcftools" ]]; then
    cd ${projdir}/preprocess
    ls *.bam > bam_list.txt
    $bcftools mpileup -Ou -f ../refgenomes/$ref1 -b bam_list.txt | \
    $bcftools call -mv -Ov -o ${projdir}/snpcall/${pop}_${ref1%.f*}_${ploidy}x_raw.vcf &&
    rm -f bam_list.txt
    cd ${projdir}/snpcall/
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

  for number_snpfilter in snpfilter snpfilter_biallelic; do
    [[ -d $number_snpfilter ]] && mv "$number_snpfilter" "${number_snpfilter}_$(ls -d snpfilter* 2>/dev/null | wc -l)"
  done


  mkdir snpfilter
  cd snpfilter

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

  if [ -z "${genotype_missingness:-}" ]; then
  	genotype_missingness=1
  else
  	genotype_missingness=$( echo $genotype_missingness | awk '{gsub(/,/," ")}1' )
  fi
  if [ -z "${sample_missingness:-}" ]; then
  	sample_missingness=1
  else
  	sample_missingness=$( echo $sample_missingness | awk '{gsub(/,/," ")}1' )
  fi
  if [ -z "${exclude_samples:-}" ]; then
  	exclude_samples=NULL
  fi
  if [ -z "${select_samples:-}" ]; then
  	select_samples=NULL
  fi
  if [[ -z "${select_samples:-}" ]]; then export select_samples=NULL; fi
  if [[ -f "${projdir}/$select_samples" ]]; then
    awk 'BEGIN {FS = "\t";OFS = ""}{print $0,".fasta.gz"}' "${projdir}/$select_samples" | \
    awk '{gsub(/.fasta.gz.fasta.gz/,".fasta.gz");}1' > "${projdir}"/fetch_samples_seq.txt
  fi
  if [[ ! -f "${projdir}/$select_samples" ]]; then cat "${projdir}"/samples_list_node_* > "${projdir}"/fetch_samples_seq.txt; fi
  if [[ ! -f "${projdir}/$select_samples" ]] && [[ "$exclude_samples" ]]; then
    echo $exclude_samples | tr ',' '\n' | awk 'BEGIN {FS = "\t";OFS = ""}{print $0,".fasta.gz"}' > "${projdir}"/fetch_excluded.txt
    grep -vFf "${projdir}"/fetch_excluded.txt "${projdir}"/fetch_samples_seq.txt > "${projdir}"/fetch_samples_seq0.txt
    mv "${projdir}"/fetch_samples_seq0.txt "${projdir}"/fetch_samples_seq.txt
  fi

  if [ -z "${minRD_1x:-}" ]; then
  	minRD_1x=2
  fi
  if [ -z "${minRD_2x:-}" ]; then
  	minRD_2x=6
  fi
  if [ -z "${minRD_4x:-}" ]; then
  	minRD_4x=25
  fi
  if [ -z "${minRD_6x:-}" ]; then
  	minRD_6x=45
  fi
  if [ -z "${minRD_8x:-}" ]; then
  	minRD_8x=100
  fi
  if [ -z "${pseg:-}" ]; then
  	pseg=0.001
  fi
  if [ -z "${maf:-}" ]; then
  	maf=0.02
  fi

  cd ${projdir}/snpcall
  for g in *_raw.vcf.gz; do gunzip $g;	done
  rm -f *_raw.vcf.gz.csi || true
  wait

  files1xG=(*_1x_DP_GT.txt)
  files1xV=(*_1x_raw.vcf)
  file1xG=${#files1xG[@]}
  file1xV=${#files1xV[@]}
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
          $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" \
          IndexFeatureFile -I $i  --verbosity ERROR&&
    			$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" \
          LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$refg -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  \
          --dont-trim-alleles --keep-original-ac --verbosity ERROR &&
    			wait
    			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf &&
    			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%10000==2{x=file"1x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
    		  wait $PID
    			rm -f "${i%.vcf}"0.vcf* "${i%.vcf}"trim.vcf*
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
  		if ! Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 1x "${GBSapp_dir}/tools/R" "1" "$filter_ExcHet"; then
        echo "Warning: Rscript failed for $pop — continuing pipeline"
      fi
      wait
      rm -f "${projdir}"/vcf1x_trimmed.txt
      wait
  	fi
  fi
  wait
  files2xG=(*_2x_DP_GT.txt)
  files2xV=(*_2x_raw.vcf)
  file2xG=${#files2xG[@]}
  file2xV=${#files2xV[@]}
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
          $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" \
          IndexFeatureFile -I $i --verbosity ERROR &&
    			$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" \
          LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$refg -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles \
          --keep-original-ac --verbosity ERROR &&
    			wait
    			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf &&
    			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%10000==2{x=file"2x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
    			wait $PID
    			rm -f "${i%.vcf}"0.vcf* "${i%.vcf}"trim.vcf*
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
  		if ! Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 2x "${GBSapp_dir}/tools/R" "1" "$filter_ExcHet"; then
        echo "Warning: Rscript failed for $pop — continuing pipeline"
      fi
      wait
      rm -f ${projdir}/vcf2x_trimmed.txt
  	fi
  fi
  wait
  files4xG=(*_4x_DP_GT.txt)
  files4xV=(*_4x_raw.vcf)
  file4xG=${#files4xG[@]}
  file4xV=${#files4xV[@]}
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
          $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" \
          IndexFeatureFile -I $i --verbosity ERROR &&
    			$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" \
          LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$refg -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles \
          --keep-original-ac --verbosity ERROR &&
    			wait
    			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf &&
    			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%10000==2{x=file"4x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
    			wait $PID
    			rm -f "${i%.vcf}"0.vcf* "${i%.vcf}"trim.vcf*
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
  		if ! Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 4x "${GBSapp_dir}/tools/R" "1" "$filter_ExcHet"; then
        echo "Warning: Rscript failed for $pop — continuing pipeline"
      fi
      wait
      rm -f ${projdir}/vcf4x_trimmed.txt
      wait
  	fi
  fi
  wait
  files6xG=(*_6x_DP_GT.txt)
  files6xV=(*_6x_raw.vcf)
  file6xG=${#files6xG[@]}
  file6xV=${#files6xV[@]}
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
          $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" \
          IndexFeatureFile -I $i --verbosity ERROR &&
    			$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" \
          LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$refg -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles \
          --keep-original-ac --verbosity ERROR &&
          wait
    			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf &&
    			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%10000==2{x=file"6x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
    			wait $PID
          rm -f "${i%.vcf}"0.vcf* "${i%.vcf}"trim.vcf*
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
  		if ! Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 6x "${GBSapp_dir}/tools/R" "1" "$filter_ExcHet"; then
        echo "Warning: Rscript failed for $pop — continuing pipeline"
      fi
      wait
      rm -f ${projdir}/vcf6x_trimmed.txt
      wait
  	fi
  fi
  wait
  files8xG=(*_8x_DP_GT.txt)
  files8xV=(*_8x_raw.vcf)
  file8xG=${#files8xG[@]}
  file8xV=${#files8xV[@]}
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
          $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" \
          IndexFeatureFile -I $i --verbosity ERROR &&
    			$GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" \
          LeftAlignAndTrimVariants -R ${projdir}/refgenomes/$refg -V $i -O ${i%.vcf}0.vcf --split-multi-allelics  --dont-trim-alleles \
          --keep-original-ac --verbosity ERROR &&
          wait
    			awk '!/^##/' ${i%.vcf}0.vcf | awk '{gsub(/^#/,""); print $0}' > ${i%.vcf}trim.vcf &&
    			awk -v file=${i%.vcf} 'BEGIN{getline f;}NR%10000==2{x=file"8x_rawSPLIT"++i".vcf";a[i]=x;print f>x;}{print > x}' ${i%.vcf}trim.vcf & PID=$!
          wait $PID
          rm -f "${i%.vcf}"0.vcf* "${i%.vcf}"trim.vcf*
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
  		if ! Rscript "${GBSapp_dir}"/scripts/R/VCF_2_DP_GT.R "${pop}" 8x "${GBSapp_dir}/tools/R" "1" "$filter_ExcHet"; then
        echo "Warning: Rscript failed for $pop — continuing pipeline"
      fi
      wait
      rm -f ${projdir}/vcf8x_trimmed.txt
      wait
  	fi
  fi
  wait


  shopt -s nullglob
  files=( *x.vcf* )
  shopt -u nullglob
  if (( ${#files[@]} == 0 )); then
    for v in *_DP_GT.txt; do
      vcfdose=${v%_DP*}
      vcfdose=${vcfdose#*_}
      shopt -s nullglob
      raw_vcfs=( *"${vcfdose}"_raw.vcf )
      shopt -u nullglob
      if (( ${#raw_vcfs[@]} == 0 )); then
        echo "WARNING: no raw VCF for dose ${vcfdose}" >&2
        continue
      fi
      for raw in "${raw_vcfs[@]}"; do
        grep '^#' "$raw" > "${raw%_raw.vcf}.vcf" &&
        awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$1,$2]){print b}}' \
            "$raw" <(awk '{print $1"\t"$2}' "$v" | awk '!seen[$0]++') \
            >> "${raw%_raw.vcf}.vcf" &&
        $bcftools sort "${raw%_raw.vcf}.vcf" | gzip > "${raw%_raw.vcf}.vcf.gz" &&
        rm -f "${raw%_raw.vcf}.vcf"
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
      if [[ -z "${p1:-}" ]]; then
      	if [ -d "${projdir}/snpfilter/1x" ]; then
      		cd ${projdir}/snpfilter &&
      		cp -r 1x 1x_diversity_gmiss"${gmiss}"_smiss"${smiss}" &&
      		cd ./1x_diversity_gmiss"${gmiss}"_smiss"${smiss}" &&
      		if ! Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_1x.R $pop $gmiss $smiss $minRD_1x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome_number $biallelic $select_samples; then
            echo "Warning: Rscript failed for $pop — continuing pipeline"
          fi
          wait
      		rm -f "${pop}"_1x_rawRD"${minRD_1x}"_DP_GT.txt "${pop}"_1x_DP_GT.txt "${pop}"_1x_rd"${minRD_1x}".txt
      		mkdir -p visualizations && mv ./*.tiff ./visualizations/ 2>/dev/null || true
          wait

          if [[ "${variant_intervals:-false}" == "true" ]]; then
            mkdir -p variant_intervals &&
            cd variant_intervals &&
            wait
            while IFS="" read -r p || [ -n "$p" ]; do
              sleep $((RANDOM % 2))
              cat ../${pop}_1x_rd${minRD_1x}_maf${maf}_dose.txt | grep "$p" | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
              paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
              awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
              rm -f variant_intervals_${p}.txt
              :> variant_intervals_${p}.txt &&
              grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
              awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
              chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep "$p" | awk '{gsub(/LN:/,""); print $3}') &&
              awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
              rm -f vbreak_break.txt variant_intervals_${p}.tmp
              wait
              shopt -s nullglob
              for vbreak in vbreak_*; do
                awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
                awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt &&
                wait
              done
              shopt -s nullglob
              sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
              mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
              rm -f vbreak_*
              wait
            done < <(awk 'NR>1{print $2}' ../${pop}_1x_rd${minRD_1x}_maf${maf}_dose.txt | sort | uniq || true)
            cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
            rm -f variant_intervals_*.txt
            cd ../
            rm -rf variant_intervals
          fi

          # Extract sequence context of variants
          shopt -s nullglob
          bam_files=( "${projdir}/preprocess/processed/"*_precall.bam )
          if (( ${#bam_files[@]} )); then
            mv "${bam_files[@]}" "${projdir}/preprocess/"
          fi
          shopt -u nullglob

          if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
              rm -f snplist_round.txt
            fi
            if [[ -n "${RE1:-}" ]]; then
              mkdir -p seq_context_"${RE1}" &&
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
                rm -f ./seq_context_"${RE1}"/${i%.f*}_seqcontext.tmp 2> /dev/null
              done < <(cat ${projdir}/fetch_samples_seq.txt)
            fi
            if [[ -n "${RE2:-}" ]]; then
              mkdir -p seq_context_"${RE2}" &&
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
                rm -f ./seq_context_"${RE2}"/${i%.f*}_seqcontext.tmp
              done < <(cat ${projdir}/fetch_samples_seq.txt)
            fi
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
                rm -f ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf
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
              mv seq_context_CCGG ./consensus_seq_context/ 2> /dev/null &&
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
      		if ! Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_2x.R $pop $gmiss $smiss $minRD_2x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome_number $biallelic $select_samples; then
            echo "Warning: Rscript failed for $pop — continuing pipeline"
          fi
          wait
      		rm -f "${pop}"_2x_rawRD"${minRD_2x}"_DP_GT.txt "${pop}"_2x_DP_GT.txt "${pop}"_2x_rd"${minRD_2x}".txt
      		mkdir -p visualizations && mv ./*.tiff ./visualizations/ 2>/dev/null || true
          wait

          if [[ "${variant_intervals:-false}" == "true" ]]; then
            mkdir -p variant_intervals &&
            cd variant_intervals &&
            wait
            while IFS="" read -r p || [ -n "$p" ]; do
              sleep $((RANDOM % 2))
              cat ../${pop}_2x_rd${minRD_2x}_maf${maf}_dose.txt | grep "$p" | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
              paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
              awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
              rm -f variant_intervals_${p}.txt
              :> variant_intervals_${p}.txt &&
              grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
              awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
              chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep "$p" | awk '{gsub(/LN:/,""); print $3}') &&
              awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
              rm -f vbreak_break.txt variant_intervals_${p}.tmp
              shopt -s nullglob
              for vbreak in vbreak_*; do
                awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
                awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt &&
                wait
              done
              shopt -s nullglob
              sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
              mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
              rm -f vbreak_*
            done < <(awk 'NR>1{print $2}' ../${pop}_2x_rd${minRD_2x}_maf${maf}_dose.txt | sort | uniq || true)
            cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
            rm -f variant_intervals_*.txt
            cd ../
            rm -rf variant_intervals
          fi

          # Extract sequence context of variants
          shopt -s nullglob
          bam_files=( "${projdir}/preprocess/processed/"*_precall.bam )
          if (( ${#bam_files[@]} )); then
            mv "${bam_files[@]}" "${projdir}/preprocess/"
          fi
          shopt -u nullglob

          if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
              rm -f snplist_round.txt
              wait
            fi
            if [[ -n "${RE1:-}" ]]; then
              mkdir -p seq_context_"${RE1}" &&
              wait
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
                rm -f ./seq_context_"${RE1}"/${i%.f*}_seqcontext.tmp
              done < <(awk 'NR==1{print $0}' ${pop}_2x_rd${minRD_2x}_maf${maf}_dose.txt | awk '{gsub(/SNP\tCHROM\tPOS\tREF\tALT\t/,""); gsub(/\t/,".fastq.gz\n");}1')
            fi
            if [[ -n "${RE2:-}" ]]; then
              mkdir -p seq_context_"${RE2}" &&
              wait
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
                rm -f ./seq_context_"${RE2}"/${i%.f*}_seqcontext.tmp
              done < <(awk 'NR==1{print $0}' ${pop}_2x_rd${minRD_2x}_maf${maf}_dose.txt | awk '{gsub(/SNP\tCHROM\tPOS\tREF\tALT\t/,""); gsub(/\t/,".fastq.gz\n");}1')
            fi
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
                rm -f ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf
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
              mv seq_context_CCGG ./consensus_seq_context/ 2> /dev/null &&
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
      		if ! Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_4x.R $pop $gmiss $smiss $minRD_4x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome_number $biallelic $select_samples; then
            echo "Warning: Rscript failed for $pop — continuing pipeline"
          fi
          wait
      		rm -f ${pop}_4x_rawRD${minRD_4x}_DP_GT.txt ${pop}_4x_DP_GT.txt ${pop}_4x_rd${minRD_4x}.txt
      		mkdir -p visualizations && mv ./*.tiff ./visualizations/ 2>/dev/null || true
          wait

          if [[ "${variant_intervals:-false}" == "true" ]]; then
            mkdir -p variant_intervals &&
            cd variant_intervals &&
            wait
            while IFS="" read -r p || [ -n "$p" ]; do
              sleep $((RANDOM % 2))
              cat ../${pop}_4x_rd${minRD_4x}_maf${maf}_dose.txt | grep "$p" | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
              paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
              awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
              rm -f variant_intervals_${p}.txt
              :> variant_intervals_${p}.txt &&
              grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
              awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
              chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep "$p" | awk '{gsub(/LN:/,""); print $3}') &&
              awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
              rm -f vbreak_break.txt variant_intervals_${p}.tmp
              shopt -s nullglob
              for vbreak in vbreak_*; do
                awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
                awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt &&
                wait
              done
              shopt -s nullglob
              sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
              mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
              rm -f vbreak_*
              wait
            done < <(awk 'NR>1{print $2}' ../${pop}_4x_rd${minRD_4x}_maf${maf}_dose.txt | sort | uniq || true)
            cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
            rm -f variant_intervals_*.txt
            cd ../
            rm -rf variant_intervals
          fi

          # Extract sequence context of variants
          shopt -s nullglob
          bam_files=( "${projdir}/preprocess/processed/"*_precall.bam )
          if (( ${#bam_files[@]} )); then
            mv "${bam_files[@]}" "${projdir}/preprocess/"
          fi
          shopt -u nullglob

          if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
              rm -f snplist_round.txt
            fi
            if [[ -n "${RE1:-}" ]]; then
              mkdir -p seq_context_${RE1} &&
              wait
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
                rm -f ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
              done < <(cat ${projdir}/fetch_samples_seq.txt)
            fi
            if [[ -n "${RE2:-}" ]]; then
              mkdir -p seq_context_${RE2} &&
              wait
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
                rm -f ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
              done < <(cat ${projdir}/fetch_samples_seq.txt)
            fi
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
                rm -f ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf
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
              mv seq_context_CCGG ./consensus_seq_context/ 2> /dev/null &&
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
      		if ! Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_6x.R $pop $gmiss $smiss $minRD_6x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome_number $biallelic $select_samples; then
            echo "Warning: Rscript failed for $pop — continuing pipeline"
          fi
          wait
      		rm -f ${pop}_6x_rawRD${minRD_6x}_DP_GT.txt ${pop}_6x_DP_GT.txt ${pop}_6x_rd${minRD_6x}.txt
      		mkdir -p visualizations && mv ./*.tiff ./visualizations/ 2>/dev/null || true
          wait

          if [[ "${variant_intervals:-false}" == "true" ]]; then
            mkdir -p variant_intervals
            cd variant_intervals
            while IFS="" read -r p || [ -n "$p" ]; do
              sleep $((RANDOM % 2))
              cat ../${pop}_6x_rd${minRD_6x}_maf${maf}_dose.txt | grep "$p" | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
              paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
              awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
              rm -f variant_intervals_${p}.txt
              :> variant_intervals_${p}.txt &&
              grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
              awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
              chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep "$p" | awk '{gsub(/LN:/,""); print $3}') &&
              awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
              rm -f vbreak_break.txt variant_intervals_${p}.tmp
              shopt -s nullglob
              for vbreak in vbreak_*; do
                awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
                awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt
                wait
              done
              shopt -s nullglob
              sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
              mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
              rm -f vbreak_*
            done < <(awk 'NR>1{print $2}' ../${pop}_6x_rd${minRD_6x}_maf${maf}_dose.txt | sort | uniq || true)
            cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
            rm -f variant_intervals_*.txt
            cd ../
            rm -rf variant_intervals
          fi

          # Extract sequence context of variants
          shopt -s nullglob
          bam_files=( "${projdir}/preprocess/processed/"*_precall.bam )
          if (( ${#bam_files[@]} )); then
            mv "${bam_files[@]}" "${projdir}/preprocess/"
          fi
          shopt -u nullglob

          if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
              rm -f snplist_round.txt
            fi
            if [[ -n "${RE1:-}" ]]; then
              mkdir -p seq_context_${RE1} &&
              wait
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' | \
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
                wait
              done < <(cat ${projdir}/fetch_samples_seq.txt)
            fi
            if [[ -n "${RE2:-}" ]]; then
              mkdir -p seq_context_${RE2}
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' | \
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
                wait
              done < <(cat ${projdir}/fetch_samples_seq.txt)
            fi
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
                rm -f ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf
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
              mv seq_context_CCGG ./consensus_seq_context/ 2> /dev/null &&
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
      		if ! Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_diversity_8x.R $pop $gmiss $smiss $minRD_8x $exclude_samples "${GBSapp_dir}/tools/R" $maf $haplome_number $biallelic $select_samples; then
            echo "Warning: Rscript failed for $pop — continuing pipeline"
          fi
          wait
      		rm -f ${pop}_8x_rawRD${minRD_8x}_DP_GT.txt ${pop}_8x_DP_GT.txt ${pop}_8x_rd${minRD_8x}.txt
      		mkdir -p visualizations && mv ./*.tiff ./visualizations/ 2>/dev/null || true
          wait

          if [[ "${variant_intervals:-false}" == "true" ]]; then
            mkdir -p variant_intervals &&
            cd variant_intervals &&
            wait
            while IFS="" read -r p || [ -n "$p" ]; do
              sleep $((RANDOM % 2))
              cat ../${pop}_8x_rd${minRD_8x}_maf${maf}_dose.txt | grep "$p" | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
              paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
              awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
              rm -f variant_intervals_${p}.txt
              :> variant_intervals_${p}.txt &&
              grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
              awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
              chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep "$p" | awk '{gsub(/LN:/,""); print $3}') &&
              awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
              rm -f vbreak_break.txt variant_intervals_${p}.tmp
              shopt -s nullglob
              for vbreak in vbreak_*; do
                awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
                awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt &&
                wait
              done
              shopt -s nullglob
              sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
              mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
              rm -f vbreak_*
            done < <(awk 'NR>1{print $2}' ../${pop}_8x_rd${minRD_8x}_maf${maf}_dose.txt | sort | uniq || true)
            cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
            rm -f variant_intervals_*.txt
            cd ../
            rm -rf variant_intervals
          fi

          # Extract sequence context of variants
          shopt -s nullglob
          bam_files=( "${projdir}/preprocess/processed/"*_precall.bam )
          if (( ${#bam_files[@]} )); then
            mv "${bam_files[@]}" "${projdir}/preprocess/"
          fi
          shopt -u nullglob

          if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
              rm -f snplist_round.txt
            fi
            if [[ -n "${RE1:-}" ]]; then
              mkdir -p seq_context_${RE1} &&
              wait
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
                rm -f ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
              done < <(cat ${projdir}/fetch_samples_seq.txt)
            fi
            if [[ -n "${RE2:-}" ]]; then
              mkdir -p seq_context_${RE2} &&
              wait
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
                rm -f ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
              done < <(cat ${projdir}/fetch_samples_seq.txt)
            fi
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
                rm -f ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf
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
              mv seq_context_CCGG ./consensus_seq_context/ 2> /dev/null &&
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
      		if ! Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_2x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_2x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$pseg" "$haplome_number" $biallelic $select_samples; then
            echo "Warning: Rscript failed for $pop — continuing pipeline"
          fi
          wait
      		rm -f "${pop}"_2x_rawRD"${minRD_2x}"_DP_GT.txt "${pop}"_2x_DP_GT.txt "${pop}"_2x_rd"${minRD_2x}".txt
      		mkdir -p visualizations && mv ./*.tiff ./visualizations/ 2>/dev/null || true
          wait

          if [[ "${variant_intervals:-false}" == "true" ]]; then
            mkdir -p variant_intervals &&
            cd variant_intervals &&
            wait
            while IFS="" read -r p || [ -n "$p" ]; do
              sleep $((RANDOM % 2))
              cat ../${pop}_2x_rd${minRD_2x}_noSDdose.txt | grep "$p" | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
              paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
              awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
              rm -f variant_intervals_${p}.txt
              :> variant_intervals_${p}.txt &&
              grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
              awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
              chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep "$p" | awk '{gsub(/LN:/,""); print $3}') &&
              awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
              rm -f vbreak_break.txt variant_intervals_${p}.tmp
              shopt -s nullglob
              for vbreak in vbreak_*; do
                awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
                awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt &&
                wait
              done
              shopt -s nullglob
              sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
              mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
              rm vbreak_*
            done < <(awk 'NR>1{print $2}' ../${pop}_2x_rd${minRD_2x}_noSDdose.txt | sort | uniq)
            cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
            rm -f variant_intervals_*.txt
            cd ../
            rm -rf variant_intervals
          fi

          # Extract sequence context of variants
          shopt -s nullglob
          bam_files=( "${projdir}/preprocess/processed/"*_precall.bam )
          if (( ${#bam_files[@]} )); then
            mv "${bam_files[@]}" "${projdir}/preprocess/"
          fi
          shopt -u nullglob

          if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
              rm -f snplist_round.txt
            fi
            if [[ -n "${RE1:-}" ]]; then
              mkdir -p seq_context_"${RE1}" &&
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
                rm -f ./seq_context_"${RE1}"/${i%.f*}_seqcontext.tmp
              done < <(cat ${projdir}/fetch_samples_seq.txt)
            fi
            if [[ -n "${RE2:-}" ]]; then
              mkdir -p seq_context_"${RE2}" &&
              wait
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
                rm -f ./seq_context_"${RE2}"/${i%.f*}_seqcontext.tmp
              done < <(cat ${projdir}/fetch_samples_seq.txt)
            fi
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
                rm -f ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf
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
              mv seq_context_CCGG ./consensus_seq_context/ 2> /dev/null &&
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
      		if ! Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_4x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_4x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$pseg" "$haplome_number" $biallelic $select_samples; then
            echo "Warning: Rscript failed for $pop — continuing pipeline"
          fi
          wait
      		rm -f ${pop}_4x_rawRD${minRD_4x}_DP_GT.txt ${pop}_4x_DP_GT.txt ${pop}_4x_rd${minRD_4x}.txt
      		mkdir -p visualizations && mv ./*.tiff ./visualizations/ 2>/dev/null || true
          wait

          if [[ "${variant_intervals:-false}" == "true" ]]; then
            mkdir -p variant_intervals &&
            cd variant_intervals &&
            wait
            while IFS="" read -r p || [ -n "$p" ]; do
              sleep $((RANDOM % 2))
              cat ../${pop}_4x_rd${minRD_4x}_noSDdose.txt | grep "$p" | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
              paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
              awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
              rm -f variant_intervals_${p}.txt
              :> variant_intervals_${p}.txt &&
              grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
              awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
              chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep "$p" | awk '{gsub(/LN:/,""); print $3}') &&
              awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
              rm -f vbreak_break.txt variant_intervals_${p}.tmp
              shopt -s nullglob
              for vbreak in vbreak_*; do
                awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
                awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt &&
                wait
              done
              shopt -s nullglob
              sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
              mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
              rm -f vbreak_*
            done < <(awk 'NR>1{print $2}' ../${pop}_4x_rd${minRD_4x}_noSDdose.txt | sort | uniq)
            cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
            rm -f variant_intervals_*.txt
            cd ../
            rm -rf variant_intervals
          fi

          # Extract sequence context of variants
          shopt -s nullglob
          bam_files=( "${projdir}/preprocess/processed/"*_precall.bam )
          if (( ${#bam_files[@]} )); then
            mv "${bam_files[@]}" "${projdir}/preprocess/"
          fi
          shopt -u nullglob

          if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
              rm -f snplist_round.txt
            fi
            if [[ -n "${RE1:-}" ]]; then
              mkdir -p seq_context_${RE1} &&
              wait
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
                rm -f ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
              done < <(cat ${projdir}/fetch_samples_seq.txt)
            fi
            if [[ -n "${RE2:-}" ]]; then
              mkdir -p seq_context_${RE2} &&
              wait
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
                rm -f ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
              done < <(cat ${projdir}/fetch_samples_seq.txt)
            fi
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
                rm -f ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf
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
              mv seq_context_CCGG ./consensus_seq_context/ 2> /dev/null &&
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
      		if ! Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_6x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_6x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$pseg" "$haplome_number" $biallelic $select_samples; then
            echo "Warning: Rscript failed for $pop — continuing pipeline"
          fi
          wait
      		rm -f ${pop}_6x_rawRD${minRD_6x}_DP_GT.txt ${pop}_6x_DP_GT.txt ${pop}_6x_rd${minRD_6x}.txt
      		mkdir -p visualizations && mv ./*.tiff ./visualizations/ 2>/dev/null || true
          wait

          if [[ "${variant_intervals:-false}" == "true" ]]; then
            mkdir -p variant_intervals &&
            cd variant_intervals &&
            wait
            while IFS="" read -r p || [ -n "$p" ]; do
              sleep $((RANDOM % 2))
              grep $p | ../${pop}_6x_rd${minRD_6x}_noSDdose.txt | grep "$p" | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
              paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
              awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
              rm -f variant_intervals_${p}.txt
              :> variant_intervals_${p}.txt &&
              grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
              awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
              chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep "$p" | awk '{gsub(/LN:/,""); print $3}') &&
              awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
              rm -f vbreak_break.txt variant_intervals_${p}.tmp
              shopt -s nullglob
              for vbreak in vbreak_*; do
                awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
                awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt &&
                wait
              done
              shopt -s nullglob
              sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
              mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
              rm -f vbreak_*
            done < <(awk 'NR>1{print $2}' ../${pop}_6x_rd${minRD_6x}_noSDdose.txt | sort | uniq)
            cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
            rm -f variant_intervals_*.txt
            cd ../
            rm -rf variant_intervals
          fi

          # Extract sequence context of variants
          shopt -s nullglob
          bam_files=( "${projdir}/preprocess/processed/"*_precall.bam )
          if (( ${#bam_files[@]} )); then
            mv "${bam_files[@]}" "${projdir}/preprocess/"
          fi
          shopt -u nullglob

          if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
              rm -f snplist_round.txt
            fi
            if [[ -n "${RE1:-}" ]]; then
              mkdir -p seq_context_${RE1}
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
                rm -f ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
              done < <(cat ${projdir}/fetch_samples_seq.txt)
            fi
            if [[ -n "${RE2:-}" ]]; then
              mkdir -p seq_context_${RE2} &&
              wait
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
                rm -f ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
              done < <(cat ${projdir}/fetch_samples_seq.txt)
            fi
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
                rm -f ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf
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
              mv seq_context_CCGG ./consensus_seq_context/ 2> /dev/null &&
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
      		if ! Rscript "${GBSapp_dir}"/scripts/R/GBSapp_Filter_8x.R "$pop" "$p1" "$p2" "$gmiss" "$smiss" "$minRD_8x" "$exclude_samples" "${GBSapp_dir}/tools/R" "$pseg" "$haplome_number" $biallelic $select_samples; then
            echo "Warning: Rscript failed for $pop — continuing pipeline"
          fi
          wait
      		rm -f ${pop}_8x_rawRD${minRD_8x}_DP_GT.txt ${pop}_8x_DP_GT.txt ${pop}_8x_rd${minRD_8x}.txt
      		mkdir -p visualizations && mv ./*.tiff ./visualizations/ 2>/dev/null || true
          wait

          if [[ "${variant_intervals:-false}" == "true" ]]; then
            mkdir -p variant_intervals &&
            cd variant_intervals &&
            wait
            while IFS="" read -r p || [ -n "$p" ]; do
              sleep $((RANDOM % 2))
              cat ../${pop}_8x_rd${minRD_8x}_noSDdose.txt | grep "$p" | awk 'NR>1{print $2"\t"$3}' | sort -n -k2,2 | awk 'NF' > variant_intervals_${p}.txt &&
              paste variant_intervals_${p}.txt <(awk 'NR>1{print $2}' variant_intervals_${p}.txt) | awk '$3==""{$3=0}1' | tr ' ' '\t' | awk '{print $0"\t"$3-$2}' | \
              awk '{gsub(/-/,"");}1' | awk -F"\t" '$4>2000{$4="break"}1' | tr ' ' '\t' > variant_intervals_${p}.tmp &&
              rm -f variant_intervals_${p}.txt
              :> variant_intervals_${p}.txt &&
              grep 'break' variant_intervals_${p}.tmp > vbreak_break.txt &&
              awk '/break/{x="vbreak_"++i;next}{print > x;}' <( printf "X\tx\tx\tbreak\n" | cat - variant_intervals_${p}.tmp) &&
              chrend=$(cat ${projdir}/refgenomes/${ref1%.fasta}.dict | grep "$p" | awk '{gsub(/LN:/,""); print $3}') &&
              awk '{print $2-500"\t"$2+500}' vbreak_break.txt | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | tr ' ' '\t' >> variant_intervals_${p}.txt &&
              rm -f vbreak_break.txt variant_intervals_${p}.tmp
              shopt -s nullglob
              for vbreak in vbreak_*; do
                awk 'NR==1{print $2;} END {print $2;}' $vbreak | tr '\n' ' ' | awk -v chrom="$p" 'BEGIN{FS=OFS=chrom"\t"}{print value OFS $0}' | \
                awk '{print $1"\t"$2-500"\t"$2+500}' >> variant_intervals_${p}.txt &&
                wait
              done
              shopt -s nullglob
              sort -n -k2,2 variant_intervals_${p}.txt | awk '$2<1{$2=1}1' | awk -v chrend=$chrend '$2>chrend{$2=chrend}1' > variant_intervals_${p}.tmp &&
              mv variant_intervals_${p}.tmp variant_intervals_${p}.txt &&
              rm vbreak_* &&
              wait
            done < <(awk 'NR>1{print $2}' ../${pop}_8x_rd${minRD_8x}_noSDdose.txt | sort | uniq)
            cat variant_intervals_*.txt | awk '{print $1":"$2"-"$3}' > ../variant_intervals.list &&
            rm -f variant_intervals_*.txt
            cd ../
            rm -rf variant_intervals
          fi

          # Extract sequence context of variants
          shopt -s nullglob
          bam_files=( "${projdir}/preprocess/processed/"*_precall.bam )
          if (( ${#bam_files[@]} )); then
            mv "${bam_files[@]}" "${projdir}/preprocess/"
          fi
          shopt -u nullglob

          if [[ "$lib_type" =~ "RRS" ]] || [[ "$lib_type" =~ "rrs" ]]; then
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
              rm -f snplist_round.txt
            fi
            if [[ -n "${RE1:-}" ]]; then
              mkdir -p seq_context_${RE1} &&
              wait
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE1}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp &&
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE1}/${i%.f*}_seqcontext.txt &&
                rm -f ./seq_context_${RE1}/${i%.f*}_seqcontext.tmp
              done < <(cat ${projdir}/fetch_samples_seq.txt)
            fi
            if [[ -n "${RE2:-}" ]]; then
              mkdir -p seq_context_${RE2} &&
              wait
              while IFS="" read -r i || [ -n "$i" ]; do
                sleep $((RANDOM % 2))
                $samtools view ../../preprocess/${i%.f*}_${ref1%.f*}_precall.bam 2> /dev/null | awk '{print $3"\t"$4"\t"$10}' | awk -v seq="^${RE2}" '$3 ~ seq' | awk '{$2=sprintf("%d00",$2/100)}1' > ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp &&
                awk -v pat="${ref1%.f*}_" '{gsub(pat,""); print $1":"$2"\t"$3}' ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp | awk 'NR==FNR {a[$1]++; next} $1 in a' <(awk -v pat="${ref1%.f*}_" '{gsub(pat,"");}1' snplist_haps.txt) - > ./seq_context_${RE2}/${i%.f*}_seqcontext.txt &&
                rm -f ./seq_context_${RE2}/${i%.f*}_seqcontext.tmp
              done < <(cat ${projdir}/fetch_samples_seq.txt)
            fi
            if [[ -n "${RE1:-}" ]] || [[ -n "${RE2:-}" ]]; then
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
                rm -f ./consensus_seq_context/"${i%.f*}"_seqcontext.cons.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.fasta.tmp ./consensus_seq_context/"${i%.f*}"_seqcontext.msf
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
              mv seq_context_CCGG ./consensus_seq_context/ 2> /dev/null &&
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
  rm -f ${projdir}/fetch_samples_seq.txt
  find . -type f -empty -delete
  find . -type d -empty -delete
  for snpfilter_dir in */; do
  	cd $snpfilter_dir &&
    mkdir -p visualizations &&
    mkdir -p visualizations
    mv ./*.tiff ./visualizations/ 2>/dev/null || true
    smmiss_thresh=${snpfilter_dir#*smiss} &&
  	smmiss_thresh=${smmiss_thresh%*/} &&
  	smmiss_thresh=$(echo "$smmiss_thresh * 100" | bc 2>/dev/null) || smmiss_t
    if compgen -G "sample_missing_rate*" > /dev/null; then
        awk -v smisst=$smmiss_thresh '(NR>1) && ($2 <= smisst)' sample_missing_rate* > retained_samples.txt
    else
        echo "" > retained_samples.txt
    fi
  	cd ../
    wait
  done
  wait
  shopt -s nullglob
  files=(*gmiss*/*dose.txt)
  if (( ${#files[@]} )); then
      wc -l "${files[@]}" | awk '{print $2"\t"$1-1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > gmiss_smiss_titration.txt
  fi
  files=(*gmiss*/eliminated*)
  if (( ${#files[@]} )); then
      wc -l "${files[@]}" | awk '{print $2"\t"$1-1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > eliminated_samples.txt
  fi
  files=(*gmiss*/retained_samples.txt)
  if (( ${#files[@]} )); then
      wc -l "${files[@]}" | awk '{print $2"\t"$1}' | grep -v "total" | awk '{sub(/\/.*$/,"",$1); print $1"\t"$2}' > retained_samples.txt
  fi
  echo -e "gmiss_smiss_thresholds\t#_of_retained_samples\t#_of_SNPs\t#_of_eliminated_samples\n----------------------\t-----------------------\t---------\t-----------------------" > summary_precall.txt &&
  awk 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,"\t",a[$1],"\t",$2}' gmiss_smiss_titration.txt eliminated_samples.txt | \
  awk 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,"\t",a[$1],"\t",$2,"\t",$3}' retained_samples.txt - | \
  cat summary_precall.txt - > gmiss_smiss.txt &&
  rm -f gmiss_smiss_titration.txt eliminated_samples.txt retained_samples.txt summary_precall.txt
  if [[ -z "$(compgen -G "${snpfilter_dir}/*dose.txt")" ]]; then
    rm -rf ${snpfilter_dir}
  fi
  wait
  shopt -s nullglob
  for f in ./*/*maf*.txt; do
      [[ "$f" == *maf0.txt ]] && continue
      [[ "$f" == *dose* ]] && continue
      [[ "$f" == *binary* ]] && continue
      rm -f "$f"
  done
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

      			if ! Rscript "${GBSapp_dir}"/scripts/R/heterozygote_vs_allele_ratio.R "$i" "$ARfile" "${ploidydir}x" "1" "darr2.txt" "${GBSapp_dir}/tools/R"; then
              echo "Warning: Rscript failed for $pop — continuing pipeline"
            fi
            wait
      			awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$2,$3]){print b}}' $ARfile $i | awk '{gsub(/NA/,"na"); print $1"_"$2"\t"$0}' | \
      			awk -v pat1="${n}_Chr" -v pat2="${n}_chr" '{gsub(pat1,"Chr");gsub(pat2,"Chr");gsub(/CHROM_POS/,"SNP");}1' > ${i%.txt}_AR_metric.txt &&
            wait


            shopt -s nullglob
            if [[ ! -f "$projdir/split_done.txt" ]]; then
                split_files=( "$projdir"/snpcall/*x.vcf* )
                if (( ${#split_files[@]} == 0 )); then
                    echo "WARNING: No VCFs to process for dose ${vcfdose}" >&2
                else
                    for split in "${split_files[@]}"; do
                        # Extract reference genome name
                        refg="${split##*/}"          # filename only
                        refg="${refg#*_}"            # remove prefix up to first underscore
                        refg="${refg%%_*}.fasta"     # remove suffix after next underscore, append .fasta
                        ungzipped="${split%.gz}"
                        gunzip -c "$split" > "$ungzipped"
                        $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" \
                            IndexFeatureFile -I "$ungzipped" --verbosity ERROR
                        # Left-align and trim variants
                        $GATK --java-options "$Xmxg -Djava.io.tmpdir=${projdir}/snpcall/tmp -XX:+UseParallelGC -XX:ParallelGCThreads=$gthreads" \
                            LeftAlignAndTrimVariants -R "${projdir}/refgenomes/$refg" -V "$ungzipped" -O "${split%.vcf}_split.vcf" \
                            --split-multi-allelics --dont-trim-alleles --keep-original-ac --verbosity ERROR
                        $bcftools view -I "${split%.vcf}_split.vcf" -O z -o "${split%.vcf}_split.vcf.gz"
                        $bcftools index "${split%.vcf}_split.vcf.gz"
                        rm -f "${split%.vcf}_split.vcf"
                        [[ "$split" != "$ungzipped" ]] && rm -f "$ungzipped"
                    done
                fi
                :> "$projdir/split_done.txt"
            fi
            shopt -u nullglob


            # Merge logic
            split_vcfs=( "$projdir"/snpcall/*_split.vcf.gz )
            if (( ${#split_vcfs[@]} > 1 )); then
                $bcftools merge "${split_vcfs[@]}" --force-samples -m all -O z -o "${i%rd*}split.vcf.gz"
            elif (( ${#split_vcfs[@]} == 1 )); then
                cp "${split_vcfs[0]}" "./${i%rd*}split.vcf.gz"
            else
                echo "WARNING: No *_split.vcf.gz files found in $projdir/snpcall/, creating empty placeholder" >&2
                touch "./${i%rd*}split.vcf.gz"
            fi
            wait || true


            for i in *dose*; do
              awk -v pat1="${n}_Chr" -v pat2="${n}_chr" '{gsub(pat1,"Chr");gsub(pat2,"Chr"); print $0}' $i > ${i%.txt}_hold.txt &&
              mv ${i%.txt}_hold.txt $i &&
              wait
            done
            wait
            zcat ${i%rd*}split.vcf.gz | grep '^#' | awk -v pat1="${n}_Chr" -v pat2="${n}_chr" '{gsub(pat1,"Chr");gsub(pat2,"Chr");gsub(/chr/,"Chr");}1' > ${i%.txt}_header.vcf &&
            awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$2,$3]){print b}}' <(zcat ${i%rd*}split.vcf.gz | grep -v '^#' | awk -v pat1="${n}_Chr" -v pat2="${n}_chr" '{gsub(pat1,"Chr");gsub(pat2,"Chr");gsub(/chr/,"Chr");}1') <(awk '!seen[$0] {print} {++seen[$0]}' $i) | \
            sort -Vk1,1 -Vk2,2 | cat ${i%.txt}_header.vcf - > ${i%dose.txt}.vcf &&
            rm -f ${i%.txt}_header.vcf
            rm -f ${i%rd*}split.vcf.gz
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
            if ! Rscript "${GBSapp_dir}"/scripts/R/hapmap_format.R "$i" "${GBSapp_dir}/tools/R"; then
              echo "Warning: Rscript failed for $pop — continuing pipeline"
            fi
            wait
            mv outfile.hmp.txt "${i%dose.txt}.hmp.txt" 2> /dev/null &&
            wait
          fi

          mv ${i%.txt}_AR_metric.txt ${i%dose.txt}_AR_metric.txt 2> /dev/null &&
          mv ${i%.txt}_AR_mean.txt ${i%dose.txt}_AR_mean.txt 2> /dev/null &&
          find . -type f -empty -delete
          rm -f AR*.txt
          wait
    		done
    		wait
    		rm -f *tmp.vcf *header.vcf darr.txt darr2.txt
    	fi
      wait
    done
    wait
  } & PID_genvcf=$!
  [[ -n "$PID_genvcf" ]] && wait "$PID_genvcf"

  rm -f ${projdir}/snpcall/*_split.vcf*
  rm -f $projdir/split_done.txt
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
	rm -f CHROM_POS.txt snp_allele_depth.txt
	cd "$projdir"/snpfilter ) &
	while (( $(jobs -rp | wc -l) >= $N )); do sleep 2; done
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
	if [ -d subref ]; then rm -rf subref; fi
	if [ -d subsamples ]; then rm -rf subsamples; fi
	if [ -d paralog_haplo_filter ]; then rm -rf paralog_haplo_filter; fi
	if test -f "$file"; then rm -rf haplo.txt; fi
	mkdir subref
	mkdir subsamples
	mkdir paralog_haplo_filter
	cd subref
	export ref1=${ref1}
	awk '{print "$samtools faidx ../../../refgenomes/${ref1} " $2 " >> refpos.fasta"}' ../snplist_nonredun.txt > snpseq_context.sh
	bash snpseq_context.sh
	wait
	rm -f *.sh
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
		rm -f ${outfile}_precall.fasta.gz ${outfile}_precall.fasta.gz ${outfile}.sam ${outfile}_uniq* ${outfile}.fasta ${outfile}_FR.fasta
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
		rm -f $sample ${sample%.fasta.gz}_alignh.txt ${sample%.fasta.gz}_aligns.txt
	done
	rm -f snpseq_context.sh snpseq_context*
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
		rm -f ${sample%_align.txt}_snp.txt ${sample%_align.txt}_haplotypes1.txt ${sample}
	done

	rmdir subsamples
	cat ./paralog_haplo_filter/*_haplotypes.txt | awk -F '\t' 'BEGIN{OFS="\t"}{print $1}' | sort -u > ./paralog_haplo_filter/SNP_files.txt
	for i in ./paralog_haplo_filter/*_haplotypes.txt; do
		happrecall=${i#*/*/}; happrecall=${happrecall%_haplotypes.txt}
		printf "SNP\t${happrecall}_haplosRD\t${happrecall}_haplos\n" > ${i%_haplotypes.txt}_hapmissingSNP.txt
		awk -F '\t' 'BEGIN{OFS="\t"}{print $1}' $i | cat - ./paralog_haplo_filter/SNP_files.txt | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{print $0"\tNA\tNA"}' | cat $i - | sort -u | sort -k1,1 >> ${i%_haplotypes.txt}_hapmissingSNP.txt
	done
	rm -f ./paralog_haplo_filter/SNP_files.txt

	touch ./paralog_haplo_filter/consensus_haplos.txt
	if [ -z "${p1:-}" ]; then
		if [ -z "${p2:-}" ]; then
			for i in ./paralog_haplo_filter/*_hapmissingSNP.txt; do
				awk '{print $3}' $i  | awk '{gsub("/","\t"); print}' | paste -d '\t' ./paralog_haplo_filter/consensus_haplos.txt - > consensus_haplos1.txt
				cat consensus_haplos1.txt | while read -r line; do
					sleep $((RANDOM % 2))
          compress=$(echo $line | tr '\t' '\n' | sort | uniq | awk 'NF > 0' | tr '\n' '\t')
					printf '..%s..' "$compress\n" >> consensus_haplos2.txt
				done
				cat consensus_haplos2.txt > ./paralog_haplo_filter/consensus_haplos.txt
				rm -f consensus_haplos1.txt && rm consensus_haplos2.txt
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
				rm -f consensus_haplos1.txt && rm -f consensus_haplos2.txt
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
	rm -f ./paralog_haplo_filter/consensus_haplo2.txt

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

	rm -f ./paralog_haplo_filter/consensus_haplo_line.txt
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
  rm -f call*
  ls ./snpfilter/*/*_plusSD.txt 2> /dev/null | xargs rm 2> /dev/null &&
  ls ./snpfilter/*/*SD_1_G*G*.txt 2> /dev/null | xargs rm 2> /dev/null &&
  ls *node*.txt | grep -v 'samples_list' | xargs rm 2> /dev/null &&
	rm -f steps.txt
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
