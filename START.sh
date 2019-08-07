#!/bin/bash
#SBATCH -J "STP_ORT"
#SBATCH --output=/data/vkarthik/ortho/logs/slurm.%j.out
#SBATCH -N 1
#SBATCH -n 8

extract_clusters()
{
grep -F -w --file=$1 OG2genes_orgs_sorted.tab | awk '{print $1"\t"$2}' > cluster_families.txt

echo "Extracted Cluster Families:"
echo $(awk '{print $1}' cluster_families.txt | sort | uniq | wc -l)

awk '{print $1}' cluster_families.txt | sort | uniq > clusters.txt

org_cnt=$(wc -l $1 | awk '{print $1}')
lower_lim=$((org_cnt * 2)) ##minimum number of transcripts you want in a cluster
##For organism specific count info use OG2genes_orgs_sorted.tab & for gene specific count info use OG2genes.tab
python get_cluster_counts.py OG2genes_orgs_sorted.tab $lower_lim> expected_counts_all.txt 

grep -F -w --file=clusters.txt expected_counts_all.txt | sort -k2 -n > expected_counts_chosen.txt

##Ideal clusters have one or more than one sequence for each organism
##so we choose clusters which are more than count of organisms and less than a scaled upper limit
upper_lim=$((org_cnt * $2))
#lower_lim=$((org_cnt * 2))
gawk -v org_cnt=$lower_lim -v upper_lim=$upper_lim '{if($2 > lower_lim && $2 < upper_lim){print}}' expected_counts_chosen.txt > chosen_clusters.txt

grep "$3" chosen_clusters.txt > chosen_specific_clusters.txt
}

get_top_orgs()
{
rm org_count.txt
while IFS= read -r line
do
cnt=$(grep -w -F --file=chosen_specific_clusters.txt OG2genes_orgs_sorted.tab | grep -w "$line" | wc -l)
printf "\n$line\t$cnt" >> org_count.txt
#printf .
done < org_list.txt
}

##ENTRYPOINT

##PARAMS- 1. NCBI Tax ID, (#3193 - Embryophyta,#33090 - Viridiplantae) must be taken from levels.tab file
##$2 is the scale for upper limit (kinda like folds, so if it's 10 then it'll be like 10*number of organisms)
##$3 i a highre level tax ID of clusters which you want to analyze ,#33090 - Viridiplantae (something from the first column of level.tab). This is useful when you want to analyze only a specific taxa level
##$4=0 if you want to include all clusters
##$4 is min number of orgs you want in a cluster

rm org_list.txt

python extract_org_ids.py levels2species.tab OGs.tab OG2genes.tab $1 > org_list.txt

###To reduce time, we convert gene based groups to organism based groups.
###We do this by splitting the second column conaining the gene IDs and checking it with org IDs
cat OG2genes.tab | awk '{split($2,a,":"); print $1 "\t" a[1]}' > OG2genes_orgs.tab

cat OG2genes_orgs.tab | sort | uniq > OG2genes_orgs_sorted.tab
#rm OG2genes_orgs.tab

extract_clusters org_list.txt $2 $3 

get_top_orgs

##CHOOSE top contributing organisms
sort -k2 -r -n org_count.txt | awk '{print $1}'| head -n $4 > top_orgs.txt

##RECHOOSE CLUSTERS AGAIN BASED ON THESE TOP ORGS
extract_clusters top_orgs.txt $2 $3

mkdir fasta
rm groups.tsv

##OLD CODE
##Sort in descending order and pick N big clusters
##sort -r -nk2 expected_counts_chosen.txt | head -n $4 | awk '{print $1}' > chosen_clusters.txt

grep -F -w --file=top_orgs.txt species.tab | awk '{print $2"\t"$3"_"$4}' > names.txt

python names2dict.py names.txt > Stop_codons

#proc_ids=()

org_cnt=$(cat top_orgs.txt | wc -l)
##WE NEED 75% OF THE CLUSTER TO BE WITHOUT NAs
min_na_cnt=$(( $org_cnt*1/2 ))
##min_org_cnt=${min_org_cnt_float%.*}
echo "Min NAs: $min_na_cnt"

awk '{print $1}' chosen_specific_clusters.txt > final_clusters.txt

input="final_clusters.txt"
while IFS= read -r line
do
 mkdir fasta/$line
 #cat species.tab | grep -E "(^|\s)$line($|\s)" | awk '{print $3}'
 #cat cluster_names.txt | grep -E "(^|\s)$line($|\s)"  > fasta/$line/name.txt
 cat OG2genes.tab | grep -E "(^|\s)$line($|\s)" | awk '{print $1"\t"$2}' > fasta/$line/group.tsv
 #uniq_org_cnt=$(awk '{split($2,a,":"); print a[1]}' fasta/$line/group.tsv | sort | uniq | wc -l)
 awk '{split($2,a,":"); print a[1]}' fasta/$line/group.tsv | sort | uniq > tmp.txt
 uniq_org_cnt=$(grep -F -w --file=top_orgs.txt tmp.txt | wc -l)
 min_org_cnt=$(( $org_cnt*1/2 ))
 if (( uniq_org_cnt > min_org_cnt ));then 
 echo "Unique org count: $uniq_org_cnt"
 cat fasta/$line/group.tsv >> groups.tsv
 ./download_data.sh fasta/$line/group.tsv $min_na_cnt #&> clusters/$line/progress.txt & 
 fi
 #proc_ids+=("$!")
done < "$input"

#for id in ${proc_ids[@]}
#do
#  wait $id
#done
