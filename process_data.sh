#!/bin/bash

filename=$1
#fastaindex -f fasta/$filename/$filename.clean.fasta --index fasta/$filename/$filename.clean.idx
fastacomposition fasta/$filename/$filename.clean.fasta ##Will show higher GC composition because it's CDS

mkdir fasta/$filename/out

##MAKE it into one single line (linearize)
#awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < fasta/$filename.clean.fasta > fasta/$filename.final.fasta

python2 trim_fasta_header.py fasta/$filename/$filename.clean.fasta 50 > fasta/$filename/$filename.trimd.fasta

#rm fasta/$filename/$filename.clean.fasta

#java -jar /data/vkarthik/tools/macse2.jar -prog translateNT2AA -keep_final_stop_ON true -out_NT fasta/$filename/$filename.trimd_NT.fasta -out_AA fasta/$filename/$filename.trimd_AA.fasta -seq fasta/$filename/$filename.trimd.fasta

#rm fasta/$filename/$filename.trimd.fasta

#~/anaconda3/lib/R/bin/Rscript align_seqs.R fasta/$filename/$filename.trimd_AA.fasta fasta/$filename/$filename.AA.aln #fasta/$filename/$filename.newick

~/anaconda3/lib/R/bin/Rscript align_seqs_init.R fasta/$filename/$filename.trimd.fasta fasta/$filename/$filename.AA.aln

##python2 /data/vkarthik/tools/RevTrans-1.4/revtrans.py fasta/$filename/$filename.clean.fasta fasta/$filename/$filename.AA.aln  -match trans > fasta/$filename/out/$filename.NT.aln
##nohup ~/anaconda3/lib/R/bin/Rscript stopcodon.R fasta/$filename/out/$filename.NT.aln fasta/$filename/out/$filename.newick &> fasta/$filename/out/$filename.log &
##wait $!
##exit

#python2 /data/vkarthik/tools/RevTrans-1.4/revtrans.py fasta/$filename/$filename.fasta fasta/$filename/$filename.AA.aln  -match trans > fasta/$filename/$filename.NT.al$

#CLUSTERING based on stop codons
python2 stop_codon_locations.py fasta/$filename/$filename.AA.aln AA fasta/$filename/$filename.dict

#ALSO removes clusters with only one sequence
python extract_cluster_sequences.py fasta/$filename/$filename.trimd.fasta fasta/$filename/$filename.dict fasta/$filename/out 3 > fasta/$filename/filelist.txt

#rm fasta/$filename/$filename.AA.aln fasta/$filename/$filename.newick fasta/$filename/$filename.fasta

#echo $clust_keys > fasta/$filename/filelist.txt

while IFS= read -r key
do
  #java -jar /data/vkarthik/tools/macse2.jar -prog translateNT2AA -keep_final_stop_ON true -out_AA fasta/$filename/out/$key.AA.fasta -seq fasta/$filename/out/$key.fasta
  ~/anaconda3/lib/R/bin/Rscript align_seqs.R fasta/$filename/out/$key.fasta fasta/$filename/$key.NT.aln
  #python2 /data/vkarthik/tools/RevTrans-1.4/revtrans.py fasta/$filename/out/$key.fasta fasta/$filename/out/$key.AA.aln  -match trans > fasta/$filename/$key.NT.aln
  sed 's/ /_/g' fasta/$filename/$key.NT.aln > fasta/$filename/out/$key.NT.aln
  rm fasta/$filename/$key.NT.aln
  nohup ~/anaconda3/lib/R/bin/Rscript stopcodon.R fasta/$filename/out/$key.NT.aln fasta/$filename/out/$key.newick "$filename.$key" &> fasta/$filename/out/$key.log &
  wait $!
  treescale=1
  if [ -s "fasta/$filename/out/$key.NT.aln.stopout" ]
  then
    cat fasta/$filename/out/$key.NT.aln.stopout >> extMG_params
    treescale=$(awk '{print $4}' fasta/$filename/out/$key.NT.aln.stopout)
  fi
  ##if .stopout exists
   #if [ -f fasta/$filename/out/$filename_$key.NT.aln.stopout ]; then
   # compseq $key.fasta -word 1 -outfile $key.freq
   # cat freq.txt | tail -n 6 | awk '{print $3}' | head -n 4 > fasta/$filename/out/$key.tmp
   # ar=()
   # while read -r fr; 
   # do 
   # ar+=($fr) 
   # done < fasta/$filename/out/$key.tmp
   # echo ${ar[@]} > fasta/$filename/out/$key.freq
   #fi

  python cut_stops.py fasta/$filename/out/$key.NT.aln "$filename.$key" fasta/$filename/out/$key.convert names.txt.dict $2 > fasta/$filename/out/$key.stops
  if [ -s "fasta/$filename/out/$key.stops" ]
  then
  ##Sub-cluster passed min NA count check & is FILE is NOT EMPTY
  perl fasta2relaxedPhylip.pl -f fasta/$filename/out/$key.convert -o fasta/$filename/out/$key.phy
  codonphyml -i fasta/$filename/out/$key.phy -m MG --fmodel F1X4 -t e -f empirical -w DM0 -q -o lr -u fasta/$filename/out/$key.newick
  
  if [ -s "fasta/$filename/out/$key.phy_codonphyml_stats.txt" ] && [ -s "fasta/$filename/out/$key.phy_codonphyml_tree.txt" ]
  then 
  omega=$(grep "Nonsynonymous/synonymous" fasta/$filename/out/$key.phy_codonphyml_stats.txt | awk '{print $4}')
  #omega=$(cat fasta/$filename/out/$key.phy_codonphyml_stats.txt | tail -n +25 | head -n 1 | awk '{print $4}')
  kappa=$(grep "Transition/transversion" fasta/$filename/out/$key.phy_codonphyml_stats.txt | awk '{print $4}')
  #kappa=$(cat fasta/$filename/out/$key.phy_codonphyml_stats.txt | tail -n +24 | head -n 1 | awk '{print $4}')
  #treescale=$( grep "MLTreeSize" $key.phy_codonphyml_stats.txt | awk '{print $3}') tree size is not the same as tree scale
  pi_G=$(grep "Nucleotides frequencies" fasta/$filename/out/$key.phy_codonphyml_stats.txt | awk '{split($6,a,"=");print a[2]}')
  #pi_G=$(cat fasta/$filename/out/$key.phy_codonphyml_stats.txt | tail -n +26 | head -n 1 | awk '{split($6,a,"=");print a[2]}')
  pi_A=$(grep "Nucleotides frequencies" fasta/$filename/out/$key.phy_codonphyml_stats.txt | awk '{split($5,a,"=");print a[2]}')
  #pi_A=$(cat fasta/$filename/out/$key.phy_codonphyml_stats.txt | tail -n +26 | head -n 1 | awk '{split($5,a,"=");print a[2]}')
  pi_C=$(grep "Nucleotides frequencies" fasta/$filename/out/$key.phy_codonphyml_stats.txt | awk '{split($4,a,"=");print a[2]}')
  #pi_C=$(cat fasta/$filename/out/$key.phy_codonphyml_stats.txt | tail -n +26 | head -n 1 | awk '{split($4,a,"=");print a[2]}')
  pi_T=$(grep "Nucleotides frequencies" fasta/$filename/out/$key.phy_codonphyml_stats.txt | awk '{split($3,a,"=");print a[2]}')
  python parameter_formatter.py "$filename.$key" $kappa $omega $pi_A $pi_C $pi_G $pi_T $treescale >> Model_parameters.mgf1x4
  cat fasta/$filename/out/$key.stops >> Stop_codons
  printf "\n$filename.$key\t$(cat fasta/$filename/out/$key.phy_codonphyml_tree.txt)" >> Trees
  else
  echo "codonphyml error"
  echo $filename >> unfinished_clusters.txt
  fi
  else
  ##min NA count check failed
  echo "min NA count check failed"
  echo $filename >> unfinished_clusters.txt
  fi
  #cat $filename_$key.newick >> Trees
  rm fasta/$filename/out/$key.AA.aln fasta/$filename/out/$key.AA.fasta fasta/$filename/out/$key.phy fasta/$filename/out/$key.fasta fasta/$filename/out/$key.convert
  echo "DONE" > fasta/$filename/DONE
done < "fasta/$filename/filelist.txt"

rm fasta/$filename/$filename.AA.aln fasta/$filename/$filename.nodup.fasta fasta/$filename/$filename.trimd.fasta fasta/$filename/$filename.fasta fasta/$filename/$filename.trimd_AA.fasta
