#!/usr/bin/env Rscript


lik_fun = function(pars) {
phi_low = pars[1]
w1 = pars[2]
 	if(phi_low < 0 | phi_low >= 1 | w1 < 0 | w1 > 1) {
	return(NA)
	}
LL = 0
	for(ii in 1:nrow(CDAlign)) {
	tree = Trees[[ii]]
	tree$edge.length = tree$edge.length * treescale[ii]
  	cdAlign = CDAlign[ii,]
	pi = PIlist[[ii]]
  	edges = tree$edge
  	TPM_p = array(rep(0,3*3*nrow(edges)),dim=c(nrow(edges),3,3))
	TPM_n = TPM_n_list[[ii]]
  		for(i in 1:nrow(edges)){ ###### Probably inefficiency here 
  		TPM_p[i,,] = expm(tree$edge.length[i]*RS[[ii]]*phi_low)
		}
	likelihood = w1 * sum(pi*felsen(tree, edges[nrow(edges),1], cdAlign,TPM_p))
	likelihood = likelihood + (1-w1) * sum(pi*felsen(tree, edges[nrow(edges),1], cdAlign,TPM_n))
  	LL = LL + log(likelihood)
	}
print(c(pars,LL))
return(LL)
}

