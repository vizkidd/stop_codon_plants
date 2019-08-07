
load('mixture_model_se_bs1000UGA.RData')
bsw = ws
bsphi = phis

load('mixture_model_se_bs1000UAA.RData')
bsw = cbind(bsw,ws)
bsphi = cbind(bsphi,phis)

load('mixture_model_se_bs1000UAG.RData')
bsw = cbind(bsw,ws)
bsphi = cbind(bsphi,phis)

load('mixture_model_se_bs1000all.RData')
bsw = cbind(bsw,ws)
bsphi = cbind(bsphi,phis)

colnames(bsw) = c("UGA","UAA","UAG","All")
colnames(bsphi) = c("UGA","UAA","UAG","All")

write.table(bsw,"Bootstraps.weights")
write.table(bsphi,"Bootstraps.phis")

