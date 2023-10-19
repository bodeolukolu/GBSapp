

java -jar snpEff.jar databases | grep 'Ipomoea'
java -jar snpEff.jar download -v Ipomoea_triloba

# Build database



# format vcf file
awk '{gsub(/TF_Chr0/,"");gsub(/TF_Chr/,"");}1' <(zcat mgs_TF_TL_6x.vcf.gz) | grep -v '^##' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' > mgs_6x.vcf 


java -Xmx16g -jar snpEff.jar -c snpEff.config -v Ipomoea_triloba mgs_6x.vcf > test.ann.vcf
