#Takes parameters: omega, phi
lik_fun_spline = function(pars,splinef, L2) {
  phi_low = pars[1]
  w1 = pars[2]
  if(phi_low < 0 | phi_low >= 1 | w1 < 0 | w1 > 1) {
    return(NA)
  }
  
  L1=matrix(0,length(splinef),1)
  for(i in 1:length(splinef)){
    f=splinef[[i]]
    L1[i]=f(phi_low)  
  }
  
  L=w1*L1+(1-w1)*L2
  LL=sum(log(L+.Machine$double.xmin))
  #print(c(pars,LL))
  return(LL)
}


#L3=linear_matrix_interpolation(sumL,discrete_phis, phi_low)
#L=(1-w1)*L2+ w1*L3
#LL=sum(log(L+.Machine$double.xmin))

