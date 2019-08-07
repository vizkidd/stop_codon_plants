#!/bin/bash

rm Trees Stop_codons Model_parameters.mgf1x4 counts.txt unfinished_clusters.txt logs/* Rplots.pdf unfinished_clusters.txt groups.tsv extMG_params
scancel $(squeue | grep "vkarthik" |awk '{print $1}' | sed '1d')

