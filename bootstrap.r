#!/usr/bin/env Rscript

#######Script to perform a bootstrap over alignments to assess uncertainty in the estimates of p and phi.

#install the following packages if you have not.
#install.packages("ape")
#install.packages("expm")

source('setup.r')
source('felsen.r')
source('lik_fun.r')
source('lik_fun_spline.r')
source('get_sumL.r')

require(expm)
require(ape)

#configuration
CDAlign_file="Stop_codons"
prep_trees_file="Trees"
writeout = "mixture_model"
print(writeout)

#model specification
#codon_frequencies_model = 'f3x4'
codon_frequencies_model = 'f1x4'
#model = 'GY'
model = 'MG'

prep_stats_file="Model_parameters.mgf1x4"


step=0.05
discrete_phis=seq(0,1,step)


setup.out = setup_all()
save(setup.out, file='setup.out.RData')

sumL=get_sumL(setup.out,discrete_phis)
save(sumL, file=paste('sumLstep', step,'.RData',sep=""))


n=dim(sumL)[1]
print ("Dims:")
print (n)
splinef=list()#vector("list", n) #list()
for(i in 1:n){
  fun_c=splinefun(discrete_phis, sumL[i,], method = "fmm", ties = mean)
  splinef[[i]]<-fun_c
  print(i)
}
save(splinef, file=paste('splinefstep', step,'.RData',sep="")) 

print("Wrote Spline")

set.seed(1)


sumL_all=sumL
splinef_all=splinef
CDAlign = setup.out$CDAlign
olist=row.names(CDAlign) 


om1=read.table("which_stop.txt",as.is=T)
print("read stops of a single org")
olist1 = om1$V1
codon = om1$V2
IND=charmatch(olist1, olist)


stopcodonlist=c("all","UGA", "UAG","UAA")

bootstrap_sample=1000
repeat_iter=1#

for(si in 1:length(stopcodonlist)){

  
  file=paste(writeout, "_se_bs",bootstrap_sample,stopcodonlist[si],".RData",sep="")
  
  if(si==1){
    SUML=sumL_all
    splinef=splinef_all
  } else {
    ind=which(codon %in% stopcodonlist[si])
    mapind=IND[ind]#Nov 15, 2018
    ind1= which(mapind!= 'NA')#Nov 15, 2018
    mapind=mapind[ind1]
    SUML=sumL_all[mapind,]#Nov 14, 2018
    splinef=splinef_all[mapind]#Nov 14, 2018
    
  }
  
  ws=rep(0,bootstrap_sample)
  phis=rep(0,bootstrap_sample)
  
  sn=dim(SUML)[1]
  
  for(bi in 1:bootstrap_sample){
    
    
    LL=rep(0,repeat_iter)
    Out<-matrix(list(), 1, repeat_iter)
    
    ind=sample(1:sn,sn,replace=TRUE)
    sumL=SUML[ind,]
    
    for(i in 1:repeat_iter){
      
      if(i==1) {
        #pars_init=c(.35,.7)
        pars_init=c(.5,.01)
      }
      if(i>1) {
        pars_init=runif(2)
      }
      
      m0.out = optim(pars_init,fn = lik_fun_spline,splinef=splinef[ind],L2=sumL[,ncol(sumL)], control=list(fnscale=-1)) 
      Out[[i]]=m0.out
      LL[i]=m0.out$value
    }
    maxi=which.max(LL)
    out=Out[[maxi]]
    
    phis[bi]=out$par[1]
    ws[bi]=out$par[2]
    
    print(c(stopcodonlist[si], out$par,out$value)) 
    
  }
  
  save(phis, ws, file=file)
  
  se=sd(ws)
}



