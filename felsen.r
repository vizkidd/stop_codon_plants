
felsen=function(tree, node, cdAlign, TPM) {
  edges = tree$edge
  s = which(edges[,1] == node)
  prod = t(t(rep(1,3)))
  for(k in s) {
      if(edges[k,2] <= length(tree$tip.label)){
          if(!is.na(cdAlign[,tree$tip.label[edges[k,2]]])) {
          prod = prod * t(t(TPM[k,,cdAlign[,tree$tip.label[edges[k,2]]]]))
          }
      }
      else {
        prod = prod * TPM[k,,]%*%felsen(tree,edges[k,2],cdAlign, TPM)
      }
  }
  return(prod)
}
