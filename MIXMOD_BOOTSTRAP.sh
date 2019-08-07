cp Trees Trees.orig
awk '{print $2}' Trees.orig > Trees ##Trees will have only trees, Trees.orig will have both trees and cluster names

##UNCOMMENT IF you want to sort, but it'll be unnecessary
##SAVE Trees as trees, Stop_codons as stops, Model_parameters.mgf1x4 as params.mgf1x4 so that we can process them for the mixture model
#cp Stop_codons stops
#cp Trees trees
#cp Model_parameters.mgf1x4 params.mgf1x4

#head -n 1 stops > Stop_codons
#sed '1d' stops | sort -k1 >> Stop_codons
#cat trees | sort -k1 | awk '{print $2}' > Trees
#cat params.mgf1x4 | sort -k1 > Model_parameters.mgf1x4

#OPTIONAL- SET  treescale values to 1
#awk -F" " -vOFS=" " '{$8=1} {print}' Model_parameters.mgf1x4

##RUN mixture model
Rscript mixture_model.rscript 

mkdir UAG
mkdir UAA
mkdir UGA

Rscript stop_specific_boot_proc.R

cp which_stop.txt mixture_model.rscript UAG/
cp which_stop.txt mixture_model.rscript UAA/
cp which_stop.txt mixture_model.rscript UGA/

Rscript UAG/mixture_model.rscript
Rscript UGA/mixture_model.rscript
Rscript UAA/mixture_model.rscript

Rscript bootstrap.r

Rscript plot.r

