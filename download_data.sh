#!/bin/bash
#SBATCH -J "STOP_CDN"
#SBATCH --output=/data/vkarthik/ortho/logs/slurm.%j.out
#SBATCH -N 1
#SBATCH -n 8

download_data()
{
#Header is trimmed
#cat $1
gene_acc=$( cat $1 ) #$(awk -F "\t" '{print $2}' $1)) #-F "\t"
count=$(cat $1 | wc -l)
i=0
echo "Number of accessions: $count"
ext="${1#*.}"
echo "Type: $ext"
filename=$3 #"${1%.$ext}"
echo "Output: $filename.fasta & $2"
#echo $filename
if [ ! -f fasta/$filename/$filename.clean.fasta ]; then
	for gene in ${gene_acc[@]}; do
	i=$[i+1]
	#echo $i
	#echo "$gene"
	og_id=$(echo "$gene" | awk '{split($0,a,":"); print a[1]}')
	if [ $(grep -w "$og_id" names.txt | wc -l) -ne 0 ]; then
		gid=$(python orthodb_api.py $gene)
		if [ "$gid" != "-1" ]; then
			sp=$(grep -w "$og_id" species.tab | awk '{print $3"_"$4}')
			printf "$gid\t$sp"
			python cds_from_ncbi.py $gid "$og_id:$gene" >> fasta/$filename/$filename.fasta
			echo "-> $(grep "$og_id:$gene" fasta/$filename/$filename.fasta | wc -l)"
			sleep .75
		else
			echo "$gene::$gid"
			printf "$filename\t$gene\n" >> $2 ##Couldn't download this sequence so store it in provided output path
			#echo "$i - $gene couldnt be found"
		fi
		#printf .
	fi
	done #< "$1" #| pv -pt -i0.2 -s $count -w 80 > /dev/null
fi
}

##ENTRYPOINT

if [ $# -eq 0 ]
  then
    echo "Give homolog groups file"
    echo "eg, ./get_data.sh groups.tsv [edit shell file and give proper delimiters]"
    exit
fi

mkdir fasta

##download_data <gene_accesseions file> <info_output>
groups=$(awk -F "\t" '{print $1}' $1 | uniq ) ##Assume 1st column has group information, use sort if you want
#file_list=$(ls *.$ext)
printf "group\tgenes.not.downloaded\n" > download_manually.txt
#echo $groups
for group in ${groups[@]}; do
filename="$group"
echo $filename
mkdir fasta/$filename
accs=$(awk -F "\t" '{if($1=="'$filename'") {print $2}}' $1 | sort | uniq )
printf "%s\n""${accs[@]}"  > fasta/$filename.txt #"%s\n"
if [ -s "fasta/$filename/DONE" ]
then ##Cluster already processed
echo "$filename Already Done"
else ##Cluster not processed so process it
download_data fasta/$filename.txt fasta/$filename.unfinished.txt $filename
cat fasta/$filename.unfinished.txt >> download_manually.txt
#cat fasta/$filename/$filename.fasta
if [ -f fasta/$filename/$filename.fasta ]; then
fastaclean fasta/$filename/$filename.fasta &> fasta/$filename/$filename.filter.fasta
#fastavalidcds fasta/$filename/$filename.filter.fasta > fasta/$filename/$filename.clean.fasta

##Remove duplicate entries
sed -e '/^>/s/$/@/' -e 's/^>/#/' fasta/$filename/$filename.filter.fasta | tr -d '\n' | tr "#" "\n" | tr "@" "\t" | sort -u -k1,1 | sed -e 's/^/>/' -e 's/\t/\n/' > fasta/temp.fa
cat fasta/temp.fa | sed '1d' > fasta/$filename/$filename.nodup.fasta
rm fasta/temp.fa

python remove_empty_sequences.py fasta/$filename/$filename.nodup.fasta > fasta/$filename/$filename.noempty.fasta

seqkit sort fasta/$filename/$filename.noempty.fasta > fasta/$filename/$filename.clean.fasta
rm fasta/$filename/$filename.noempty.fasta
printf "\nSaved $(fastalength fasta/$filename/$filename.clean.fasta | wc -l) sequences\n"

python select_sequences.py fasta/$filename/$filename.clean.fasta names.txt.dict fasta/$filename/$filename.clean.fasta

fi

printf "\nSelected $(fastalength fasta/$filename/$filename.clean.fasta | wc -l) sequences\n"

gene_cnts=$(fastalength fasta/$filename/$filename.clean.fasta | wc -l)
printf "$filename\t$gene_cnts\n" >> counts.txt
sbatch --output=fasta/$filename/progress.txt -J $filename -N 1 -n 8 process_data.sh $filename $2 #$2 is the min number of NAs allowed in the cluster
echo "------------------------------------------"
#CLEANUP
rm fasta/$filename.unfinished.txt fasta/$filename.txt fasta/$filename/$filename.filter.fasta #fasta/$filename.clean.fasta
fi
done

