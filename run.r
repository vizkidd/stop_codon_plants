#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

model = 'MG'
codon_frequencies_model = 'f1x4'

#####Requires ape and expm
require(ape)
require(expm)

tree_file = args[1]
seq_file = args[2]
out_file = paste(seq_file,".stopout",sep="")

source('setup_generatormatrix.r')
source('felsenstein.r')
source('like_func.r')
source('output.r')
