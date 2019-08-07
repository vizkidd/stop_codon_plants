bootstraps = 100

weights<-data.frame("UGA"=c(1:1000),"UAA"=c(1:1000),"UAG"=c(1:1000),"All"=c(1:1000))
phi_est<-data.frame("UGA"=c(1:1000),"UAA"=c(1:1000),"UAG"=c(1:1000),"All"=c(1:1000))
load("mixture_model_se_bs1000UGA.RData")
weights[,1]<-ws
phi_est[,1]<-phis
load("mixture_model_se_bs100UAA.RData")
weights[,2]<-ws
phi_est[,2]<-phis
load("mixture_model_se_bs100UAG.RData")
weights[,3]<-ws
phi_est[,3]<-phis
load("mixture_model_se_bs100all.RData")
weights[,4]<-ws
phi_est[,4]<-phis
write.table(weights,"Bootstraps.weights")
write.table(phi_est,"Bootstraps.phis")

##WEIGHTS
Ws = read.table("Bootstraps.weights")
Phis = read.table("Bootstraps.phis")

##UNCOMMENT if you have low proportion of any stop or simply to experiment
##Weight not estimable for phi close to one (which can occur for STOPS, due to low proportion under selection)
##Phi not estimable for weight close to zero

#Ws$UAG[Phis$UAG > 0.99] = NA  
#Phis$UAG[Ws$UAG < 0.01] = NA  

#Ws$UGA[Phis$UGA > 0.99] = NA  ##Weight not estimable for phi close to one (which can occur for UAG, due to low proportion under selection)
#Phis$UGA[Ws$UGA < 0.01] = NA

#Ws$UAA[Phis$UAA > 0.99] = NA  ##Weight not estimable for phi close to one (which can occur for UAG, due to low proportion under selection)
#Phis$UAA[Ws$UAA < 0.01] = NA

proportion = c(Ws$UGA,Ws$UAG,Ws$UAA,Ws$All)
stop_codon = c(rep('UGA',bootstraps),rep('UAG',bootstraps),rep('UAA',bootstraps),rep('All',bootstraps))
df = data.frame(proportion=proportion,stop_codon = stop_codon)
colnames(df) = c("proportion","stop codon")

##ADD xlims, ylims if necessary
ggplot(df, aes(proportion, fill = `stop codon`)) + theme_bw(base_size = 20) + geom_histogram( aes(y = ..density..), position = 'identity', binwidth=0.005) + scale_fill_manual(values=c("black", "orange", "gold","chocolate4"))

##PHI
all.phis = c(Phis$UGA,Phis$UAG,Phis$UAA,Phis$All)
df.phi = data.frame(phi=all.phis,stop_codon = stop_codon)
colnames(df.phi) = c("phi","stop codon")

##ADD xlims, ylims if necessary
ggplot(df.phi, aes(phi, fill = `stop codon`)) + theme_bw(base_size = 20) + geom_histogram( aes(y = ..density..), position = 'identity', binwidth=0.005) + scale_fill_manual(values=c("black", "orange", "gold","chocolate4"))
