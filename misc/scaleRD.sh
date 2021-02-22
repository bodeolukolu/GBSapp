
nfiles=6
cores=6
mkdir scaleRD
cd samples

for i in $(ls -S *.f* | grep -v R2.f | head -n $nfiles); do (
	if gzip -t $i; then
		printf 'Unique_RD\tFrequency\tTotal_Reads\tPercentage_Cumulative\n---------\t---------\t-----------\t---------------------\n' > ../scaleRD/${i%.f*}_uniqRD_histogram.txt
		gunzip -c $i | awk 'NR%2==0' | awk 'NR%2' | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' |\
		awk '{print $1}' | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{print $2"\t"$1}' | sort -k1,1 -n |\
		awk '{print $1,"\t",$2,"\t"$1*$2}' > ../scaleRD/${i%.f*}_uniqRD_hold1.txt
		awk 'NR==FNR{sum+= $3; next;} {printf("%s\t%s\t%3.0f\t%3.2f%%\t%3.0f\n",$1,$2,$3,100*$3/sum,100*$3/sum)}' ../scaleRD/${i%.f*}_uniqRD_hold1.txt ../scaleRD/${i%.f*}_uniqRD_hold1.txt |\
		awk '{total += $4; $4 = total}1' | awk '{total += $5; $5 = total}1' > ../scaleRD/${i%.f*}_uniqRD_hold2.txt
		unset IFS; printf "%s\t%s\t%s\t%s\t%*s\n" $(sed 's/$/ |/' ../scaleRD/${i%.f*}_uniqRD_hold2.txt) | tr ' ' '|' >> ../scaleRD/${i%.f*}_uniqRD_histogram.txt
		rm ../scaleRD/${i%.f*}_uniqRD_hold*
	else
		printf 'Unique_RD\tFrequency\tTotal_Reads\tPercentage_Cumulative\n---------\t---------\t-----------\t---------------------\n' > ../scaleRD/${i%.f*}_uniqRD_histogram.txt
		awk 'NR%2==0' $i | awk 'NR%2' | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' |\
		awk '{print $1}' | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk '{print $2"\t"$1}' | sort -k1,1 -n |\
		awk '{print $1,"\t",$2,"\t"$1*$2}' > ../scaleRD/${i%.f*}_uniqRD_hold1.txt
		awk 'NR==FNR{sum+= $3; next;} {printf("%s\t%s\t%3.0f\t%3.2f%%\t%3.0f\n",$1,$2,$3,100*$3/sum,100*$3/sum)}' ../scaleRD/${i%.f*}_uniqRD_hold1.txt ../scaleRD/${i%.f*}_uniqRD_hold1.txt  |\
		awk '{total += $4; $4 = total}1' | awk '{total += $5; $5 = total}1' > ../scaleRD/${i%.f*}_uniqRD_hold2.txt
		unset IFS; printf "%s\t%s\t%s\t%s\t%*s\n" $(sed 's/$/ |/' ../scaleRD/${i%.f*}_uniqRD_hold2.txt) | tr ' ' '|' >> ../scaleRD/${i%.f*}_uniqRD_histogram.txt
		rm ../scaleRD/${i%.f*}_uniqRD_hold* 
	fi ) &
	if [[ $(jobs -r -p | wc -l) -ge $cores ]]; then
		wait
	fi
done
