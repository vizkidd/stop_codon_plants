#!/usr/bin/env Rscript


##############
#Given a Tree with branch lengths and a parameter file (including kappa, omega and nucleotide equilibrium frequencies), 
#simulate coding sequences from the stop-extended codon model, based on MG94 X HKY85 X F1x4
#Usage: Rscript sim_stop.rscript treefile parameterfile
#Format of the parameter file: textfile containing the following 8 values on a single line: gene_name kappa omega pi_A pi_C pi_G pi_T alignment_length_in_codons
#(where the pis represent equilibrium nucleotide frequencies)
#
#Outputs a file called <gene_name>.sim.fastas containing a simulated coding sequence alignment, including
#stop codon
##############

#Requires ape and expm
require(ape)
require(expm)

treescale = 1

args = commandArgs(trailingOnly=TRUE)
tree_file = args[1]
tree = read.tree(tree_file)
numtaxa = length(tree$tip.label)
param_file = args[2]
stats = read.table(param_file,row.names=1)
gene = rownames(stats)
kappa = stats[1,1]
omega = stats[1,2]
phi = stats[1,3]
ff = unlist(stats[1,4:7]) # equilibrium nucleotide frequencies
cds_len = stats[1,8]
cdAlign = array(rep(0,numtaxa*cds_len),c(numtaxa,cds_len))
dimnames(cdAlign)[[1]]=tree$tip.label
out_file = paste(gene,".sim.mgf1x4.fastas",sep="")

amino_acids = c("F", "L", "S", "Y", "C", "W", "P", "H", "Q", "R", "I", "M", "T", "K", "N", "V", "A", "D", "E", "G","X")
codons = c("ttt", "ttc", "tta", "ttg", "tct", "tcc", "tca", "tcg", "tat", "tac", "tgt", "tgc", "tgg", "ctt", "ctc", "cta", "ctg", "cct", "ccc", "cca", "ccg", "cat", "cac","caa", "cag", "cgt", "cgc", "cga", "cgg", "att", "atc", "ata", "atg", "act", "acc", "aca", "acg", "aat", "aac","aaa", "aag", "agt", "agc", "aga", "agg", "gtt", "gtc","gta", "gtg", "gct", "gcc", "gca", "gcg", "gat", "gac","gaa", "gag", "ggt", "ggc", "gga", "ggg", "tag", "tga","taa")
aa = c("F", "F", "L", "L", "S", "S", "S", "S", "Y", "Y","C", "C", "W", "L", "L", "L", "L", "P", "P", "P", "P","H", "H", "Q", "Q", "R", "R", "R", "R", "I", "I", "I","M", "T", "T", "T", "T", "N", "N", "K", "K", "S", "S","R", "R", "V", "V", "V", "V", "A", "A", "A", "A", "D","D", "E", "E", "G", "G", "G", "G", "X", "X", "X")
codon_numbers = list(ttt = 1, ttc = 2, tta = 3, ttg = 4,tct = 5, tcc = 6, tca = 7, tcg = 8, tat = 9, tac = 10,tgt = 11, tgc = 12, tgg = 13, ctt = 14, ctc = 15, cta = 16,ctg = 17, cct = 18, ccc = 19, cca = 20, ccg = 21, cat = 22,cac = 23, caa = 24, cag = 25, cgt = 26, cgc = 27, cga = 28,cgg = 29, att = 30, atc = 31, ata = 32, atg = 33, act = 34,acc = 35, aca = 36, acg = 37, aat = 38, aac = 39, aaa = 40,aag = 41, agt = 42, agc = 43, aga = 44, agg = 45, gtt = 46,gtc = 47, gta = 48, gtg = 49, gct = 50, gcc = 51, gca = 52,gcg = 53, gat = 54, gac = 55, gaa = 56, gag = 57, ggt = 58,ggc = 59, gga = 60, ggg = 61, tag = 62, tga = 63, taa = 64)
nucleotides = list(a = 1, c = 2, g = 3, t = 4)
purine = c(1, 0, 1, 0)

######Setup########
setup = function(omega,phi,kappa,treefile) {
tree = read.tree(treefile)
tree=reorder(tree,"postorder")
edges = tree$edge
#set up PI
PI = array(rep(0, 64*64), dim = c(64, 64))
PI1x4 = array(rep(0, 64*64), dim = c(64, 64))

	for (i in 1:64) {
	PI1x4[i, i] = ff[nucleotides[[substr(codons[i], 1, 1)]]] * ff[nucleotides[[substr(codons[i],2, 2)]]]* ff[nucleotides[[substr(codons[i],3, 3)]]]
	}

PI = PI1x4/sum(PI1x4)
Mg_mult = array(c(rep(0, 4096)), dim = c(64, 64))

###############Set up the generator matrix############
R = array(c(rep(0, 4096)), dim = c(64, 64))
	for (i in 1:64) {
		for (j in 1:64) {
		diffs = 0
			for (k in 1:3) {
			nuc1 = nucleotides[[substr(codons[i], k, k)]]
			nuc2 = nucleotides[[substr(codons[j], k, k)]]
				if (nuc1 != nuc2) {
				diffs = diffs + 1
				Mg_mult[i,j] = ff[nuc2]
					if (purine[nuc1] == purine[nuc2]) {
					R[i, j] = kappa
					}
					else {
					R[i,j] = 1
					}
					if (aa[i] != aa[j]) {
					R[i,j] = R[i,j] * omega 
					}
					if(i > 61 & j <= 61) { # stop to non stop
					R[i,j] = 0
					}
					if(i <= 61 & j > 61) { # non stop to stop
					R[i,j] = 0
					}
					if(i > 61 & j > 61 & i != j) {
					R[i,j] = phi * R[i,j]
					}
				}
			}
			if (diffs > 1) {
			R[i, j] = 0
			}
		}
	}

R = R * Mg_mult
scale_fac = sum(PI%*%R)
R = (1/scale_fac)*R

	for(i in 1:64) {
	R[i,i] = -sum(R[i,-i])
	}

return(list(edges=edges,tree=tree,R=R,PI=PI))
}

#####Felsenstein pruning algorithm
felsen_sim=function(anc_codon, pos, tree, node, TPM) {
edges = tree$edge
s = which(edges[,1] == node)
	for(k in s) {
		if(edges[k,2] <= length(tree$tip.label)){
			if(!is.null(dim(cdAlign))) {
				if(!is.na(cdAlign[tree$tip.label[edges[k,2]],pos])) {
				P = t(t(TPM[k,anc_codon,]))
	  			desc_codon = sample(1:64,1,prob=P)
	  			cdAlign[tree$tip.label[edges[k,2]],pos] <<- desc_codon
				}
			}
			else {
				if(!is.na(cdAlign[tree$tip.label[edges[k,2]]])) {
				P = t(t(TPM[k,anc_codon,]))
	  			desc_codon = sample(1:64,1,prob=P)
				}
			}
		}
		else {
		P = t(t(TPM[k,anc_codon,]))
	  	desc_codon = sample(1:64,1,prob=P)
		felsen_sim(desc_codon,pos,tree,edges[k,2], TPM)
		}
	}
}


########End of setup########

out = setup(omega,phi,kappa,tree_file)
edges = out$edges
tree = out$tree
R = out$R
PI = out$PI
TPM = array(rep(0,64*64*nrow(edges)),dim=c(nrow(edges),64,64))

for(i in 1:nrow(edges)){
TPM[i,,] = expm(treescale*tree$edge.length[i]*R)
}

pi = diag(PI)
pi1 = pi[1:61]/sum(pi[1:61])
pi2 = pi[62:64]/sum(pi[62:64])

for(i in 1:(ncol(cdAlign)-1)) {
root_codon = sample(1:61,1,prob=pi1)
felsen_sim(root_codon,i,tree,edges[nrow(edges),1],TPM)
}

root_codon = sample(1:3,1,prob=pi2)
root_codon = root_codon + 61
felsen_sim(root_codon,i+1,tree,edges[nrow(edges),1],TPM)
pst = paste(codons[unlist(cdAlign[1,])],sep="")
speciesnames = rownames(cdAlign)

write("",out_file)
for(i in 1:nrow(cdAlign)) {
write(paste(">",speciesnames[i],sep=""),out_file,append=T)
write(paste(codons[unlist(cdAlign[i,])],collapse=""),out_file,append=T)
}
