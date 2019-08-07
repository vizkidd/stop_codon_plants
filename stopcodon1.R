#!/usr/bin/env Rscript

## Usage: Rscript stopcodon.rscript <treefile.ph> <seqfile.fasta>
## The tree file should include branch lengths. Note that relative branch lengths are treated as fixed.
## The sequence file should contain a codon-aware alignment in fasta format and must include the 
## stop codon as the last position of the alignment. Sequences that include a gap or a codon other than
## a stop codon (i.e. sequences for which the stop codon is not positionally homologous with the last 
## position in the alignment) are excluded.
## Output: a file called stopcodon.rscript.out, which contains the following: The maximum log likelihood;
## ML estimate of kappa; ML estimate of omega; ML estimate of the treescaling parameter; ML estimate of phi
## delta_lnL (difference in log likelihood compared to a model with phi=1; convergence of optimizer (0 = success,
## 1 = failure)

##Kappa is transition/transversion rate

args = commandArgs(trailingOnly=TRUE)

model = 'MG'
codon_frequencies_model = 'f1x4'
header_trim_len<- 50

#####Requires ape and expm
suppressMessages(require(ape))
suppressMessages(require(expm))
suppressMessages(library(DECIPHER))
suppressMessages(library(Biostrings))
suppressMessages(library(phylotools))
suppressMessages(library(lmtest))

######Setup########
setup = function(omega,phi,kappa,treefile, seqfile) {
  amino_acids = c("F", "L", "S", "Y", "C", "W", "P", "H", "Q", "R", "I", "M", "T", "K", "N", "V", "A", "D", "E", "G","X")
  codons = c("ttt", "ttc", "tta", "ttg", "tct", "tcc", "tca", "tcg", "tat", "tac", "tgt", "tgc", "tgg", "ctt", "ctc", "cta", "ctg", "cct", "ccc", "cca", "ccg", "cat", "cac","caa", "cag", "cgt", "cgc", "cga", "cgg", "att", "atc", "ata", "atg", "act", "acc", "aca", "acg", "aat", "aac","aaa", "aag", "agt", "agc", "aga", "agg", "gtt", "gtc","gta", "gtg", "gct", "gcc", "gca", "gcg", "gat", "gac","gaa", "gag", "ggt", "ggc", "gga", "ggg", "tag", "tga","taa")
  aa = c("F", "F", "L", "L", "S", "S", "S", "S", "Y", "Y","C", "C", "W", "L", "L", "L", "L", "P", "P", "P", "P","H", "H", "Q", "Q", "R", "R", "R", "R", "I", "I", "I","M", "T", "T", "T", "T", "N", "N", "K", "K", "S", "S","R", "R", "V", "V", "V", "V", "A", "A", "A", "A", "D","D", "E", "E", "G", "G", "G", "G", "X", "X", "X")
  codon_numbers = list(ttt = 1, ttc = 2, tta = 3, ttg = 4,tct = 5, tcc = 6, tca = 7, tcg = 8, tat = 9, tac = 10,tgt = 11, tgc = 12, tgg = 13, ctt = 14, ctc = 15, cta = 16,ctg = 17, cct = 18, ccc = 19, cca = 20, ccg = 21, cat = 22,cac = 23, caa = 24, cag = 25, cgt = 26, cgc = 27, cga = 28,cgg = 29, att = 30, atc = 31, ata = 32, atg = 33, act = 34,acc = 35, aca = 36, acg = 37, aat = 38, aac = 39, aaa = 40,aag = 41, agt = 42, agc = 43, aga = 44, agg = 45, gtt = 46,gtc = 47, gta = 48, gtg = 49, gct = 50, gcc = 51, gca = 52,gcg = 53, gat = 54, gac = 55, gaa = 56, gag = 57, ggt = 58,ggc = 59, gga = 60, ggg = 61, tag = 62, tga = 63, taa = 64)
  nucleotides = list(a = 1, c = 2, g = 3, t = 4)
  purine = c(1, 0, 1, 0)
  
  tree = read.tree(tree_file)
  tree=reorder(tree,"postorder")
  #tree=reorder(tree_file,"postorder")
  edges = tree$edge
  
  #seqs = read.dna(seq_file,fo="fasta")
  chars = seqfile #as.character(seqs)
  if(class(chars)=="list"){
    alignment_len<-length(chars[[1]])
    data.clean<-matrix(,nrow=length(chars),ncol = alignment_len)
    for (i in 1:length(chars)) {
      if(length(chars[[i]])==alignment_len){
        temp<-t(matrix(unlist(chars[[i]])))
        #rownames(temp)<-names(chars[i])
        data.clean[i,]<-temp[1,]
      }
    }
    rownames(data.clean)<-names(chars)
    chars<-data.clean
  }
  if(!all(sort(tree$tip.label)==sort(rownames(chars)))){
    ##Trimming the headers
    tree$tip.label=substr(toupper(tree$tip.label),1,header_trim_len)
    rownames(chars)=substr(toupper(rownames(chars)),1,header_trim_len)
  }
  chars = chars[tree$tip.label,]
  
  #  if(ncol(chars)%%3==1)
  #  {
  ##Shift frames to left/right because they are not in frame (adding a column of gaps at 1st column)
  ##CHEATING
  #    chars<-cbind(as.matrix(rep("-",nrow(chars)),nrow=nrow(chars)),as.matrix(rep("-",nrow(chars)),nrow=nrow(chars)),chars)
  #  }
  
  #set up PI
  f1 = c(0, 0, 0, 0) 
  f2 = c(0, 0, 0, 0)
  f3 = c(0, 0, 0, 0)
  ff = c(0,0,0,0)
  PI = array(rep(0, 64*64), dim = c(64, 64))
  PI3x4 = array(rep(0, 64*64), dim = c(64, 64))
  PI1x4 = array(rep(0, 64*64), dim = c(64, 64))
  s_c=0
  #convert to codons - produces cdAlign which is a matrix representing the coding sequence alignment. Each
  #row is a sequence, each column is a numeric representation (1-64) of the codon
  cdAlign = array(rep(0,nrow(chars)*ncol(chars)/3),c(nrow(chars),ncol(chars)/3))
  for(i in 1:nrow(chars)) {
    s = paste(chars[i,],collapse="")
    for(j in 1:(ncol(chars)/3)) {
      ss = substr(s,3*(j-1)+1,3*j)
      ss=tolower(ss)
      #if(ss!="+++" && ss!="---")
      #{print(ss)}
      if(!is.null(codon_numbers[[ss]]) && ss!="+++" && ss!="---") {
        if(ss=="taa"||ss=="tag"||ss=="tga"){
          ##counting stop codons ##print(ss)
          s_c=s_c+1
        }
        cdAlign[i,j] = codon_numbers[[ss]]
        f1[nucleotides[[substr(ss,1,1)]]] = f1[nucleotides[[substr(ss,1,1)]]] + 1
        f2[nucleotides[[substr(ss,2,2)]]] = f2[nucleotides[[substr(ss,2,2)]]] + 1
        f3[nucleotides[[substr(ss,3,3)]]] = f3[nucleotides[[substr(ss,3,3)]]] + 1
        ff[nucleotides[[substr(ss,1,1)]]] = ff[nucleotides[[substr(ss,1,1)]]] + 1
        ff[nucleotides[[substr(ss,2,2)]]] = ff[nucleotides[[substr(ss,2,2)]]] + 1
        ff[nucleotides[[substr(ss,3,3)]]] = ff[nucleotides[[substr(ss,3,3)]]] + 1
      }
      else {
        cdAlign[i,j] = NA
      }
    }
  }
  #print(c("SC:",s_c,"nSeqs",nrow(chars)))
  dimnames(cdAlign)[[1]]=tree$tip.label
  ##cdAlign cleanup empty/masked rows
  mskd_cols<-c()
  for(i in 1:ncol(cdAlign)){
    if(all(is.na(cdAlign[,i]))){
      mskd_cols<-c(mskd_cols,i)
    }
  }
  #mskd_cols=unique(mskd_cols)
  ##IF mskd_cols is not empty
  if(!is.null(mskd_cols)){
    cdAlign<-cdAlign[,-mskd_cols]
  }
  
  #Frame based proportions of nucleotides(on each position of the frame)?
  f1 = f1/sum(f1)
  f2 = f2/sum(f2)
  f3 = f3/sum(f3)
  ff = ff/sum(ff)
  ##Construct PI based on model (1x4,3x4) (1x4=same rates based on global proption, 3x4=induvidual rates)
  for (i in 1:64) {
    PI3x4[i, i] = f1[nucleotides[[substr(codons[i], 1, 1)]]] * f2[nucleotides[[substr(codons[i],2, 2)]]]* f3[nucleotides[[substr(codons[i],3, 3)]]]
    PI1x4[i, i] = ff[nucleotides[[substr(codons[i], 1, 1)]]] * ff[nucleotides[[substr(codons[i],2, 2)]]]* ff[nucleotides[[substr(codons[i],3, 3)]]]
  }
  
  # Toggle f3x4 f1x4
  #PI = PI3x4/sum(PI3x4)
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
          if(codon_frequencies_model == 'f1x4') {
            Mg_mult[i,j] = ff[nuc2]
          } else if(codon_frequencies_model == 'f3x4') {
            if(k==1) {Mg_mult[i,j] = f1[nuc2]}
            if(k==2) {Mg_mult[i,j] = f2[nuc2]}
            if(k==3) {Mg_mult[i,j] = f3[nuc2]}
          } else {print("Codon frequencies model undefined")
            return(NA)
          }
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
  
  #model defines instantaneous rate matrix (or TPMs?) [more info: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4379412/]
  ##GY- target frequencies in TPM, MG- target frequencies in PI
  ##GY- target frequencies are suseptible to change, MG- frequencies are constant
  ##MG-style models; the first, MG1, employs four parameters for nucleotide frequencies (one per nucleotide; [Muse and Gaut 1994](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4379412/#msv003-B43))
  ##
  ##dN/dS-based model matrices account for nucleotide mutational bias by incorporating either target codon (Goldman and Yang 1994) or target nucleotide (Muse and Gaut 1994) frequencies; these frameworks are known, respectively, as Goldman–Yang (GY)-style and Muse–Gaut (MG)-style models (Kosakovsky Pond et al. 2010). For example, the instantaneous rate matrix element giving the substitution probability from codon AAA to AAG would contain the target codon frequency PAAG in GY-style models but the target nucleotide frequency πG in MG-style models
  ##GY-style matrices may be expressed in the framework of the GTR model, in which the instantaneous matrix Q can be decomposed into a 61 × 61 symmetric substitution rate matrix and a 61-dimensional vector containing the equilibrium codon frequencies. The latter corresponds to the stationary distribution of the Markov chain. In contrast, MG-style rate matrices are written in terms of nucleotide frequencies rather than codon frequencies. Therefore, whether these models fit into the GTR framework is unclear a priori. We now describe how the MG-style matrix can be rewritten in terms of a symmetric matrix and a vector of equilibrium codon frequencies, thus demonstrating that these matrices also fit into the GTR(General time reversible) framework.
  
  ##NOTE:ASK CATHAL ABOUT THIS PART OF CODE
  if(model == 'MG') {
    R = R * Mg_mult
  } else if(model == 'GY') {
    R = R %*% PI
  } else {print("Model undefined")
    return(NA)
  }
  
  #Scale the generator matrix
  scale_fac = sum(PI%*%R)
  ##Normalize it so that it doesn't exceed 1
  R = (1/scale_fac)*R
  
  ##Re-estimate generator matrix
  ##Sum each codon row without the column corresponding to the same codon an negate it
  ##Probability of not getting that codons is the reason for negating the sum? is it shorthand way of writing , p(codon) = -(codon occuring - codon not occuring)
  for(i in 1:64) {
    R[i,i] = -sum(R[i,-i])
  }
  
  return(list(edges=edges,tree=tree,R=R,cdAlign=cdAlign,PI=PI))
}

#####Felsenstein pruning algorithm
felsen=function(pos, tree, node, cdAlign, TPM) {
  edges = tree$edge
  prod = t(t(rep(1,64)))
  ##We see which edges are connected to the node
  s = which(edges[,1] == node)
  for(k in s) {
    if(edges[k,2] <= length(tree$tip.label)){
      if(!is.null(dim(cdAlign))) {
        if(!is.na(cdAlign[tree$tip.label[edges[k,2]],pos])) {
          prod = prod * t(t(TPM[k,,cdAlign[tree$tip.label[edges[k,2]],pos]]))
        }
      }
      else {
        if(!is.na(cdAlign[tree$tip.label[edges[k,2]]])) {
          prod = prod * t(t(TPM[k,,cdAlign[tree$tip.label[edges[k,2]]]]))
        }
      }
    }
    else {
      prod = prod * TPM[k,,]%*%felsen(pos,tree,edges[k,2],cdAlign, TPM)
    }
  }
  return(prod)
}


########End of setup########

#Takes parameters: kappa, omega, treescale, phi
lik_fun = function(pars,treefile,seqfile,phifixed=0) {
  kappa = pars[1]
  omega = pars[2]
  treescale = pars[3]
  #phi = pars[4]
  if(phifixed==1) {
    phi = 1
  } else {
    phi = pars[4]
  }
  
  if(phi < 0 | omega < 0 | kappa < 0 | treescale < 0) {
    return(NA)
  }
  plot_ll<<-c()
  out = setup(omega,phi,kappa,treefile,seqfile)
  edges = out$edges
  tree = out$tree
  ##cdAlign returning all NAs
  cdAlign = out$cdAlign
  R = out$R
  PI = out$PI
  ##N sequences have 2N edges. we get the probability of each edge as a TPM, ie, the probability of edge transitioning and transforming into another sequence.
  TPM = array(rep(0,64*64*nrow(edges)),dim=c(nrow(edges),64,64))
  ##R is crashing here
  for(i in 1:nrow(edges)){
    ##TPM is the probability of transition of each codon in the column of all sequences mutating to another codon. there are 64 slices, one slice for each codon. (why do we take edge though?)
    TPM[i,,] = expm(treescale*tree$edge.length[i]*R)
  }
  likelihood = 0
  for(i in 1:ncol(cdAlign)) {
    ##NAs produced qhen i=4 and trying to add cdAlign & TPM doesn't work(non conformable?), pruning works
    ##NODE=edges[nrow(edges),1],(why do we take only the last row of edges, which corresponds to the longest edge), guess we prune the longest branch
    likelihood = likelihood + log(sum(diag(PI)*felsen(i, tree, edges[nrow(edges),1], cdAlign,TPM)))
    #print(likelihood)
    plot_ll<<-c(plot_ll,likelihood)
    #if(is.infinite(likelihood)==TRUE){
    # print(i)
    # break ##Remove this
    #}
  }
  
  #count=count+1
  print(c(pars,likelihood))
  return(likelihood)
}

##ENTRYPOINT
#MAKE SURE TO CHANGE '!'/'?' to N before proceeding
#tree_file = "cathal.ph" #args[1]
#seq_file = "cathal.fasta" #args[2] "../fasta/temp/3228.tree.newick"  "../fasta/temp/3228.NT.fasta"
tree_file = args[2]
seq_file = args[1]
#seq_file = "../fasta/PHG-01950.clean_NT.fasta" #args[2]
#aa_aln="../fasta/PHG-01950.clean_AA1.fasta"
#tree_file = "temp.newick" #args[1]
#seq_file = "temp.fasta" #args[2]
out_file = paste(seq_file,".stopout",sep="")

##continued from ALIGN_SEQS.R

##To reverse peptide alignment, use 
##revtrans fasta/PHG-01950.clean.fasta fasta/PHG.AA.aln.fasta  -match trans > fasta/PHG-01950.clean_NT.fasta
#python trim_fasta_header.py fasta/PHG-01950.clean_NT_NT.fasta 10 > fasta/PHG-01950.clean_NT_NT1.fasta 

##Current ENTRYPOINT
plot_ll<<-c()
clean.fasta.name(infile = seq_file,outfile = seq_file)
data.aln<-readDNAMultipleAlignment(seq_file,format = "fasta")
d <- DistanceMatrix(myXStringSet = data.aln@unmasked,correction = "Jukes-Cantor", verbose = FALSE)
MODELS

##MAKE SURE to change '!'/'?' to Ns before proceeding
#ML approach only works for DNA seqs
#clu <- IdClusters(d, method="ML",  showPlot = TRUE,model="TN93+G4", myXStringSet=data.aln@unmasked)
#print(clu)
dend<- IdClusters(d, method="NJ", type="dendrogram", showPlot = TRUE)
plot(dend)
str(dend)

#Dont mask because it masks stops
#dna.masked<-MaskAlignment(data.aln@unmasked)
#masks <- lapply(width(colmask(dna.masked)), rep, x="+")
#masks <- unlist(lapply(masks, paste, collapse=""))
#masked_dna <- replaceAt(dna.masked@unmasked, at=IRanges(colmask(dna.masked)), value=masks)
#BrowseSeqs(masked_dna)

#masked_dna<-as.matrix(masked_dna)
masked_dna<-as.matrix(as(data.aln,"DNAStringSet"))
WriteDendrogram(dend,file=tree_file,quoteLabels = FALSE)

##phi=1 is neutral selection
##phi<1 is purifying selection
m0.out = optim(c(2,0.2,1),lik_fun,treefile=tree_file,seqfile=masked_dna,phifixed=1,control=list(fnscale=-1))
print(m0.out)
plot(plot_ll)
#m1.out = optim(c(2,0.2,1,0),lik_fun,treefile=tree_file,seqfile=seq_file,control=list(fnscale=-1))
m1.out = optim(c(m0.out$par,0.5),lik_fun,treefile=tree_file,seqfile=masked_dna,phifixed=0,control=list(fnscale=-1))
print(m1.out)
plot(plot_ll)

m2.out = optim(c(m1.out$par),lik_fun,treefile=tree_file,seqfile=masked_dna,control=list(fnscale=-1))
print(m2.out)

#m3.out = optim(c(m0.out$par,1.5),lik_fun,treefile=tree_file,seqfile=masked_dna,control=list(fnscale=-1))
#print(m3.out)

lrtest(m2.out,m1.out)

diff = m2.out$value - m0.out$value
write(c(m2.out$value,m2.out$par,diff,m2.out$convergence),out_file,ncol=10)

#m4.out = optim(c(m2.out$par),lik_fun,treefile=tree_file,seqfile=masked_dna,control=list(fnscale=-1))
#print(m4.out)

