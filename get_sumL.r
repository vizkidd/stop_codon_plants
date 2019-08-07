get_sumL = function(setup.out,discrete_phis) {
  require(expm)
  CDAlign = setup.out$CDAlign
  Trees = setup.out$Trees
  TPM = setup.out$TPM
  PIlist = setup.out$PIlist
  n=length(PIlist)
  m=length(discrete_phis)
  sumL=array(c(rep(0, n*m)), dim = c(n, m))
  for(ii in 1:n) {
    tree = Trees[[ii]]
    cdAlign = CDAlign[ii,]#Nov 10, 2018
    pi = PIlist[[ii]]
    edges = tree$edge
    
    for(jj in 1:m) {
      phi_low=discrete_phis[jj]
      TPM_p = array(rep(0,3*3*nrow(edges)),dim=c(nrow(edges),3,3))#Oct 7, 2018
      for(i in 1:nrow(edges)){#Oct 7, 2018
        TPM_p[i,,] = expm(tree$edge.length[i]*setup.out$RS[[ii]]*phi_low)
      }
      sumL[ii,jj]=sum(pi*felsen(tree, edges[nrow(edges),1], cdAlign,TPM_p))#Nov 11, 2018 -- for new felsen 
    }
    print(ii)
  }
  return(sumL)
}
