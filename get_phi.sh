#!/bin/bash

phi_lim=$1

echo "PHI_LIM=$phi_lim"

#cat $(find fasta/ | grep "stopout") | awk '{print NR":"$5}'

#cat $(find fasta/ | grep "stopout") | awk -v phi=$phi_lim '{if($5>"$phi"){print NR}}' > line_numbers.txt

#while IFS= read -r line; 
#do 
#awk -v line=$line '{if(NR==$line){print $1}}' Stop_codons; 
#echo "$line"
#done < line_numbers.txt

cat $(find fasta/ | grep "stopout") | awk '{if($5>5)print NR}' > line_numbers.txt

