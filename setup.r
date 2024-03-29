#!/usr/bin/env Rscript
setup_all = function() {
  require(ape)
  require(expm)
  
  #writeout = "fast_discrete.out"
  
  amino_acids = c("F", "L", "S", "Y", "C", "W", "P", "H", "Q", "R", "I", "M", "T", "K", "N", "V", "A", "D", "E", "G","X")
  codons = c("ttt", "ttc", "tta", "ttg", "tct", "tcc", "tca", "tcg", "tat", "tac", "tgt", "tgc", "tgg", "ctt", "ctc", "cta", "ctg", "cct", "ccc", "cca", "ccg", "cat", "cac","caa", "cag", "cgt", "cgc", "cga", "cgg", "att", "atc", "ata", "atg", "act", "acc", "aca", "acg", "aat", "aac","aaa", "aag", "agt", "agc", "aga", "agg", "gtt", "gtc","gta", "gtg", "gct", "gcc", "gca", "gcg", "gat", "gac","gaa", "gag", "ggt", "ggc", "gga", "ggg", "tag", "tga","taa")
  aa = c("F", "F", "L", "L", "S", "S", "S", "S", "Y", "Y","C", "C", "W", "L", "L", "L", "L", "P", "P", "P", "P","H", "H", "Q", "Q", "R", "R", "R", "R", "I", "I", "I","M", "T", "T", "T", "T", "N", "N", "K", "K", "S", "S","R", "R", "V", "V", "V", "V", "A", "A", "A", "A", "D","D", "E", "E", "G", "G", "G", "G", "X", "X", "X")
  codon_numbers = list(ttt = 1, ttc = 2, tta = 3, ttg = 4,tct = 5, tcc = 6, tca = 7, tcg = 8, tat = 9, tac = 10,tgt = 11, tgc = 12, tgg = 13, ctt = 14, ctc = 15, cta = 16,ctg = 17, cct = 18, ccc = 19, cca = 20, ccg = 21, cat = 22,cac = 23, caa = 24, cag = 25, cgt = 26, cgc = 27, cga = 28,cgg = 29, att = 30, atc = 31, ata = 32, atg = 33, act = 34,acc = 35, aca = 36, acg = 37, aat = 38, aac = 39, aaa = 40,aag = 41, agt = 42, agc = 43, aga = 44, agg = 45, gtt = 46,gtc = 47, gta = 48, gtg = 49, gct = 50, gcc = 51, gca = 52,gcg = 53, gat = 54, gac = 55, gaa = 56, gag = 57, ggt = 58,ggc = 59, gga = 60, ggg = 61, tag = 62, tga = 63, taa = 64)
  nucleotides = list(a = 1, c = 2, g = 3, t = 4)
  purine = c(1, 0, 1, 0)
  
  PI = array(rep(0,64*64),c(64,64))
  PIlist = list()
  TPM_n_list = c()
  CDAlign = read.table(CDAlign_file,h=T)
  print("read stops")
  
  #######Read pre-processed output from codonphyml (kappa, omega and equilibrium frequencies)
  stats = read.table(prep_stats_file)
  print("read model params")
  kappa = stats[,2]
  omega = stats[,3]
  treescale = stats[,8]
  
  
  ###############Set up the template generator matrix############
  RS = list() # A list with a separate generator matrix for each orthologue family
  
  ###Basic stop-extended R matrix (with omega,kappa = 1)
  R = array(c(rep(0, 4096)), dim = c(64, 64))
  for (i in 1:64) {
    for (j in 1:64) {
      diffs = 0
      for (k in 1:3) {
        nuc1 = nucleotides[[substr(codons[i], k, k)]]
        nuc2 = nucleotides[[substr(codons[j], k, k)]]
        if (nuc1 != nuc2) {
          diffs = diffs + 1
          R[i,j] = 1
          if(i > 61 & j <= 61) { # stop to non stop
            R[i,j] = 0
          }
          if(i <= 61 & j > 61) { # non stop to stop
            R[i,j] = 0
          }
        }
      }
      if (diffs > 1) {
        R[i, j] = 0
      }
    }
  }
  #####################
  
  ################Identify the positions in the matrix to be multipled by gene-specific factors (equilibrium nucleotide/codon frequencies,kappa,omega)
  f1rows = c()
  f2rows = c()
  f3rows = c()
  f1cols = c()
  f2cols = c()
  f3cols = c()
  f1nuc = c()
  f2nuc = c()
  f3nuc = c()
  kapparows = c()
  kappacols = c()
  omegarows = c()
  omegacols = c()
  for (i in 1:64) {
    for (j in 1:64) {
      diffs = 0
      for (k in 1:3) {
        nuc1 = nucleotides[[substr(codons[i], k, k)]]
        nuc2 = nucleotides[[substr(codons[j], k, k)]]
        if (nuc1 != nuc2 & R[i,j] != 0) {
          if (purine[nuc1] == purine[nuc2]) {
            kapparows = c(kapparows,i)
            kappacols = c(kappacols,j)
          }
          if (aa[i] != aa[j]) {
            omegarows = c(omegarows,i)
            omegacols = c(omegacols,j)
          }
          if(k==1) {
            f1rows = c(f1rows,i)
            f1cols = c(f1cols,j)
            f1nuc = c(f1nuc,nuc2)
          }
          if(k==2) {
            f2rows = c(f2rows,i)
            f2cols = c(f2cols,j)
            f2nuc = c(f2nuc,nuc2)
          }
          if(k==3) {
            f3rows = c(f3rows,i)
            f3cols = c(f3cols,j)
            f3nuc = c(f3nuc,nuc2)
          }
        }
      }
    }
  }
  
  ####################################################################
  
  #Trees = read.tree("prep_trees.temp")
  #Trees = read.tree("rm.test.tree")
  Trees = read.tree(prep_trees_file) 
  Trees = reorder(Trees,"postorder")
  
  
  ##############Calculate a matrix of codon equilibrium frequencies (rows are genes, columns are codons)
  if(codon_frequencies_model == 'f3x4') {
    PI3x4 = array(rep(0,nrow(stats)*64),c(nrow(stats),64))
    f1 = stats[,4:7]
    f2 = stats[,8:11]
    f3 = stats[,12:15]
    for (i in 1:64) {
      PI3x4[, i] = f1[,nucleotides[[substr(codons[i], 1, 1)]]] * f2[,nucleotides[[substr(codons[i],2, 2)]]]* f3[,nucleotides[[substr(codons[i],3, 3)]]]
    }
  } else if(codon_frequencies_model == 'f1x4') {
    PI1x4 = array(rep(0,nrow(stats)*64),c(nrow(stats),64))
    ff = stats[,4:7]
    f1 = ff
    f2 = ff
    f3 = ff
    for (i in 1:64) {
      PI1x4[, i] = ff[,nucleotides[[substr(codons[i], 1, 1)]]] * ff[,nucleotides[[substr(codons[i],2, 2)]]]* ff[,nucleotides[[substr(codons[i],3, 3)]]]
    }
  } else {
    print("Codon frequencies model undefined")
    return(NA)
  }
  
  
  ############### The following multiplier matrix, Rmult, is required for MG model only. Rmult_ko is for kappa and omega (required for MG and GY)
  for(ii in 1:nrow(CDAlign)) {
    Rmult = array(c(rep(1, 4096)), dim = c(64, 64))
    Rmult_ko = array(c(rep(1, 4096)), dim = c(64, 64))
    for(i in 1:length(f1rows)) {
      Rmult[f1rows[i],f1cols[i]] = f1[ii,f1nuc[i]]
    }
    for(i in 1:length(f2rows)) {
      Rmult[f2rows[i],f2cols[i]] = f2[ii,f2nuc[i]]
    }
    for(i in 1:length(f3rows)) {
      Rmult[f3rows[i],f3cols[i]] = f3[ii,f3nuc[i]]
    }
    for(i in 1:length(kapparows)) {
      Rmult_ko[kapparows[i],kappacols[i]] = kappa[ii] * Rmult_ko[kapparows[i],kappacols[i]]
    }
    for(i in 1:length(omegarows)) {
      Rmult_ko[omegarows[i],omegacols[i]] = omega[ii] * Rmult_ko[omegarows[i],omegacols[i]]
    }
    
    if(codon_frequencies_model == 'f3x4') {
      diag(PI) = PI3x4[ii,]
    } else if(codon_frequencies_model == 'f1x4') {
      diag(PI) = PI1x4[ii,]
    }
    
    
    #RR is just a temporary intermediate between the R template and RS
    RR = Rmult_ko * R
    if(model == 'MG') {
      RR = Rmult * RR
    } else if(model == 'GY') {
      RR = RR%*%PI
    }
    
    
    #Calculate the scale factor on the coding part of the rate matrix only (to make consistent with tree length estimates from codonphyml)
    scale_fac = sum(PI[1:61,1:61]%*%RR[1:61,1:61])
    
    R_p = RR[62:64,62:64]
    for(i in 1:3) {
      R_p[i,i] = -sum(R_p[i,-i])
    }
    
    R_p = (1/scale_fac)*R_p
    
    RS[[ii]] = R_p
    PIlist[[ii]] = diag(PI)[62:64]
    
    ###### Calculate the TPM matrix for phi = 1
    tree = Trees[[ii]]
    edges = tree$edge
    tree$edge.length = tree$edge.length * treescale[ii]
    TPM_n = array(rep(0,3*3*nrow(edges)),dim=c(nrow(edges),3,3))
    for(i in 1:nrow(edges)){ ###### Probably inefficiency here 
      TPM_n[i,,] = expm(tree$edge.length[i]*RS[[ii]])
    }
    TPM_n_list[[ii]] = TPM_n
  }
  return(list(RS=RS,CDAlign = CDAlign,Trees=Trees,PIlist=PIlist,TPM_n_list=TPM_n_list))
  
}
############### 
