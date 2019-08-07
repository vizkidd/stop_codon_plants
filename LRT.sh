#!/bin/bash

cat $(find fasta/ | grep "stopout")  > diff.txt

while IFS= read -r line;
do
l=$(awk '{print $6}')
lrt=$(bc <<< 'scale=10; 2*'$l)
if []
done < diff.txt


