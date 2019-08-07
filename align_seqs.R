if(length(args)==0){
  print("Give seq_file,aa alignment outpt, tree_file output names as args")
  quit(status=-1)
}

suppressMessages(library(DECIPHER))
suppressMessages(library(Biostrings))
suppressMessages(library(phylotools))

args = commandArgs(trailingOnly=TRUE)

header_trim_len<- 50

seq_file = args[1]
dna_aln_file=args[2]
#tree_file = args[3]

#clean.fasta.name(infile = seq_file,outfile = seq_file)
data.dna<-readDNAStringSet(seq_file,format = "fasta")

#data.aa@ranges@NAMES=substr(toupper(data.aa@ranges@NAMES),1,header_trim_len)
#data.aa<-translate(data.nt,genetic.code = GENETIC_CODE,if.fuzzy.codon = "solve")
#data.nt<- CorrectFrameshifts(data.nt,myAAStringSet = data.aa,type = "both")
#data.aa<-translate(data.nt,genetic.code = GENETIC_CODE,if.fuzzy.codon = "solve")
#data.aa.hec<-PredictHEC(data.aa, type="probabilities")
#data.aa.aln<- AlignSeqs(data.aa,iterations = length(data.aa),refinements = length(data.aa),alphabet = AA_STANDARD,useStructures = TRUE,structures = data.aa.hec)
#data.aa.aln<-RemoveGaps(data.aa.aln,removeGaps = "common")

data.dna.aln<-AlignTranslation(data.dna,iterations=0,refinements=0)
writeXStringSet(data.dna.aln,filepath = dna_aln_file,format = "fasta")



